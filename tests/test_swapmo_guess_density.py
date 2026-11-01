"""swapmo on a supplied guess must swap the guess density with the orbitals.

A supplied guess (JSON, resident orbitals) brings its own OQP::DM_A/DM_B,
and the SCF may accept orbitals and density together on its first build.
swapmo swapped only the orbital rows, so the density stayed that of the
unswapped occupation: a JSON ground state with swapmo=1,6 converged on the
first build at the ground-state energy with core-hole orbitals.

A computed guess (Huckel here) keeps its density: the swap only names the
occupation that rstctmo/MOM then follow, and the delta-CHP core-hole examples
(examples/XAS/delta-CHP-MRSF) converge from that start.  Runs OpenQP in a
subprocess.
"""

import json
import os
import subprocess
import sys
import textwrap

import pytest

np = pytest.importorskip("numpy")

SYSTEM = """
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.633194329  -0.483194329  -0.554469223"""

# scf type -> (multiplicity, swapmo, n_alpha, n_beta)
STATES = {
    "rhf": (1, "5,6", 5, 5),
    "rohf": (3, "1,5", 6, 4),
    "uhf": (1, "1,6", 5, 5),
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

SCRIPT = """
import json
import numpy as np
from oqp.openqp import OPENQP

cfg = {{"input.runtype": "energy", "input.method": "hf", "input.basis": "6-31g",
        "input.charge": 0, "input.system": {system!r}, "scf.type": {scf!r},
        "scf.multiplicity": {mult}}}
cfg.update({guess!r})
op = OPENQP(cfg, True)
op.sp._prep_guess()
before = {{k: np.asarray(op.mol.data["OQP::" + k], dtype=float).copy()
          for k in ("VEC_MO_A", "DM_A", "DM_B")}}
op.sp.swapmo()
d = op.mol.data
out = {{k: np.asarray(d["OQP::" + k], dtype=float).tolist()
        for k in ("VEC_MO_A", "VEC_MO_B", "DM_A", "DM_B")}}
out.update({{"before_" + k: v.tolist() for k, v in before.items()}})
print("JSON" + json.dumps(out))
"""


DECK = """\
[input]
system={system}
charge=0
runtype=energy
basis=6-31g
method=hf

[guess]
type=huckel
save_mol=True

[scf]
type={scf}
multiplicity={mult}
conv=1.0e-8
"""


def _env():
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    return env


def _supplied_json(tmp_path, scf):
    """Converged SCF of the unswapped state, saved for a type=json guess."""
    inp = tmp_path / "ground.inp"
    inp.write_text(DECK.format(system=SYSTEM, scf=scf, mult=STATES[scf][0]))
    proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(inp)], capture_output=True,
                          text=True, cwd=str(tmp_path), env=_env(), timeout=900)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    return tmp_path / "ground.json"


def _script(tmp_path, scf, guess):
    mult = STATES[scf][0]
    body = SCRIPT.format(system=SYSTEM, scf=scf, mult=mult, guess=guess)
    path = tmp_path / "probe.py"
    path.write_text(textwrap.dedent(body))
    proc = subprocess.run([sys.executable, str(path)], capture_output=True, text=True,
                          cwd=str(tmp_path), env=_env(), timeout=900)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    line = [l for l in proc.stdout.splitlines() if l.startswith("JSON")][-1]
    return json.loads(line[4:])


def _unpack(packed, n):
    m = np.zeros((n, n))
    m[np.tril_indices(n)] = np.asarray(packed, dtype=float).ravel()
    return m + np.tril(m, -1).T


def _expected(scf, va, vb, na, nb):
    if scf == "rhf":
        return 2.0 * va[:na].T @ va[:na], None
    vb = va if scf == "rohf" else vb
    return va[:na].T @ va[:na], vb[:nb].T @ vb[:nb]


@pytest.mark.parametrize("scf", ["rhf", "rohf", "uhf"])
def test_swapmo_on_a_supplied_guess_swaps_the_density(tmp_path, scf):
    _, swap, na, nb = STATES[scf]
    js = _supplied_json(tmp_path, scf)
    d = _script(tmp_path, scf, {"guess.type": "json", "guess.file": str(js),
                                "guess.swapmo": swap})
    va, vb = np.asarray(d["VEC_MO_A"]), np.asarray(d["VEC_MO_B"])
    n = va.shape[1]
    exp_a, exp_b = _expected(scf, va, vb, na, nb)
    assert np.abs(_unpack(d["DM_A"], n) - exp_a).max() < 1e-12
    if exp_b is not None:
        assert np.abs(_unpack(d["DM_B"], n) - exp_b).max() < 1e-12
    if scf in ("rhf", "rohf"):
        assert np.abs(vb - va).max() == 0.0     # beta follows the swapped alpha set
    # The swap moved occupied weight, so the density really changed (ROHF
    # swapmo=1,5 stays inside the alpha set and changes only DM_B).
    changed = max(np.abs(_unpack(d[k], n) - _unpack(d["before_" + k], n)).max()
                  for k in ("DM_A", "DM_B"))
    assert changed > 1e-2


@pytest.mark.parametrize("scf", ["rhf", "rohf", "uhf"])
def test_swapmo_on_a_computed_guess_keeps_the_guess_density(tmp_path, scf):
    _, swap, _, _ = STATES[scf]
    d = _script(tmp_path, scf, {"guess.type": "huckel", "guess.swapmo": swap})
    assert np.asarray(d["DM_A"]) == pytest.approx(np.asarray(d["before_DM_A"]), abs=0.0)
    assert np.asarray(d["DM_B"]) == pytest.approx(np.asarray(d["before_DM_B"]), abs=0.0)
    # ... while the orbital rows were swapped.
    assert np.abs(np.asarray(d["VEC_MO_A"]) - np.asarray(d["before_VEC_MO_A"])).max() > 1e-2
