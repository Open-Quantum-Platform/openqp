"""TRAH is refused for MOM and rstctmo calculations.

TRAH minimises over orbital rotations and never re-diagonalises, so neither
occupation rule can act: from a non-aufbau stationary point it descended to
a lower state (H2O/6-31G UHF core hole -> -75.985 Eh) while reporting
convergence.  The input checker and scf_driver refuse the combination.
Runs OpenQP in a subprocess.
"""

import json
import os
import subprocess
import sys
import textwrap

import pytest

np = pytest.importorskip("numpy")

# Distorted water, no symmetry; 13 basis functions in 6-31G.
SYSTEM = """\
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.633194329  -0.483194329  -0.554469223"""

DECK = """\
[input]
system=
{system}
charge=0
runtype=energy
basis=6-31g
method=hf

[guess]
type=huckel
swapmo={swap}
save_mol=True

[scf]
type={scf}
multiplicity={mult}
maxit=100
conv=1.0e-8
{extra}
"""

# scf type -> (multiplicity, swapmo, n_alpha, n_beta)
STATES = {
    "uhf": (1, "1,6", 5, 5),    # alpha O 1s -> LUMO core excitation
    "rhf": (1, "5,6", 5, 5),    # HOMO^2 -> LUMO^2 double excitation
    "rohf": (3, "1,5", 6, 4),   # O 1s singly occupied core-hole triplet
}


def _have_native_oqp():
    try:
        import oqp  # noqa: F401
    except Exception:
        return False
    return True


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(),
    reason="native OpenQP library not importable",
)


def _env(threads):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    return env


def _run_deck(tmp_path, name, deck, threads=1):
    inp = tmp_path / f"{name}.inp"
    inp.write_text(deck)
    proc = subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", str(inp)],
        capture_output=True, text=True, cwd=str(tmp_path), env=_env(threads),
        timeout=900,
    )
    log_path = inp.with_suffix(".log")
    log = log_path.read_text() if log_path.exists() else ""
    return proc, log


def _energy(tmp_path, name, deck, threads=1):
    proc, log = _run_deck(tmp_path, name, deck, threads)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert "SCF convergence achieved" in log, log[-3000:]
    return json.loads((tmp_path / f"{name}.json").read_text())["energy"]


def _run_script(tmp_path, body, threads=1):
    script = tmp_path / "probe.py"
    script.write_text(textwrap.dedent(body))
    return subprocess.run(
        [sys.executable, str(script)],
        capture_output=True, text=True, cwd=str(tmp_path), env=_env(threads),
        timeout=900,
    )


@pytest.mark.parametrize("rule", ["mom=True", "rstctmo=True"])
def test_input_with_trah_and_occupation_rule_is_refused(tmp_path, rule):
    deck = DECK.format(system=SYSTEM, swap="1,6", scf="uhf", mult=1,
                       extra=f"{rule}\nconverger_type=trah")
    proc, log = _run_deck(tmp_path, "refused", deck)
    assert proc.returncode != 0
    assert "SCF convergence achieved" not in log
    assert "TRAH cannot keep" in proc.stdout + proc.stderr + log


DRIVER_SCRIPT = """
from oqp.openqp import OPENQP

cfg = {{"input.runtype": "energy", "input.method": "hf", "input.basis": "6-31g",
        "input.charge": 0, "input.system": {system!r}, "scf.type": "uhf",
        "scf.multiplicity": 1, "guess.type": "huckel", "guess.swapmo": "1,6",
        "scf.{rule}": True}}
op = OPENQP(cfg, True)
op.sp._prep_guess()
op.sp.swapmo()
op.mol.data.set_scf_converger_type("trah")   # as an escalation stage would
op.sp.scf()
print("CONVERGED", op.mol.mol_energy.SCF_converged)
"""


@pytest.mark.parametrize("rule", ["mom", "rstctmo"])
def test_native_driver_refuses_trah_with_occupation_rule(tmp_path, rule):
    proc = _run_script(tmp_path, DRIVER_SCRIPT.format(system="\n" + SYSTEM, rule=rule))
    assert proc.returncode != 0
    assert "CONVERGED" not in proc.stdout
    logs = "".join(p.read_text() for p in tmp_path.glob("*.log"))
    assert "TRAH cannot keep" in proc.stdout + proc.stderr + logs
