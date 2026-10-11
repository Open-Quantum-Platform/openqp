"""The first MOM step must compare with the guess orbitals.

handle_mom skipped its first call and stored the aufbau fill of iteration 1
as the reference, so even an active MOM (mom_switch above the first DIIS
error) followed the refilled state and lost a non-aufbau guess.  H2O/6-31G
UHF, swapmo=1,6, mom_switch=100: -75.985 Eh (ground state) instead of the
core-excited -56.162 Eh that rstctmo reaches.  Runs OpenQP in a subprocess.
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


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("scf", ["uhf", "rhf"])
def test_mom_from_swapmo_guess_reaches_restricted_mo_state(tmp_path, scf, threads):
    mult, swap, _, _ = STATES[scf]
    deck = DECK.format(system=SYSTEM, swap=swap, scf=scf, mult=mult, extra="{extra}")
    e_rstct = _energy(tmp_path, "rstct", deck.format(extra="rstctmo=True"), threads)
    e_mom = _energy(tmp_path, "mom", deck.format(extra="mom=True\nmom_switch=100.0"),
                    threads)
    assert e_mom == pytest.approx(e_rstct, abs=1e-6)
