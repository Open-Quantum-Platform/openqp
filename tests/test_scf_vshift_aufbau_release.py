"""A level shift must not decide which state the SCF converges to.

The level shift raises the virtuals of the orbitals in hand, so it keeps the
occupation the iterations started from.  The zero-VSHIFT build then tested
convergence (FDS - SDF, small at any stationary point) before any unshifted
Fock had chosen the occupation again, and accepted a non-aufbau state:

    UHF triplet H2O/6-31G PBE, vshift 0.3:  -75.9351089 Eh (alpha b2 occupied
    above an empty a1) instead of -76.0143957 Eh; DIIS and SOSCF alike.
    ROHF triplet, vshift 1.0:               -75.9343919 instead of -76.0139994.
    RHF restarted from a HOMO^2 -> LUMO^2 solution, vshift 0.3: stayed at
    -74.9932296 instead of the ground state -75.9852335.

Each shifted run must reach the energy of the unshifted SCF from a normal
guess, and its stored density must be that of the lowest orbitals of the
final Fock.  Runs OpenQP in a subprocess.
"""

import json
import os
import subprocess
import sys

import pytest

np = pytest.importorskip("numpy")

DECK = """\
[input]
system=
 8   0.000000000   0.000000000  -0.041061554
 1  -0.533194329   0.533194329  -0.614469223
 1   {h2}
charge=0
method=hf
basis=6-31g
runtype=energy
{func}

[guess]
{guess}
save_mol=true

[scf]
type={scf}
multiplicity={mult}
maxit=200
conv=1.0e-7
{extra}

[dftgrid]
rad_npts=96
ang_npts=302
"""
H2_C2V = "0.533194329  -0.533194329  -0.614469223"
H2_C1 = "0.633194329  -0.483194329  -0.554469223"
RELEASED = "occupation released"


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


def _run(tmp_path, name, threads=1, **fmt):
    inp = tmp_path / f"{name}.inp"
    inp.write_text(DECK.format(**fmt))
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    proc = subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", str(inp)],
        capture_output=True, text=True, cwd=str(tmp_path), env=env, timeout=1800,
    )
    log_path = inp.with_suffix(".log")
    log = log_path.read_text() if log_path.exists() else ""
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert "SCF convergence achieved" in log, log[-3000:]
    data = json.loads((tmp_path / f"{name}.json").read_text())
    return data, log


def _unpack(packed, n):
    m = np.zeros((n, n))
    m[np.tril_indices(n)] = np.asarray(packed, dtype=float).ravel()
    return m + np.tril(m, -1).T


def _assert_aufbau_uhf(data, nalpha, nbeta):
    """UHF publishes the aufbau orbitals of the final Fock: they must hold DM."""
    smat = _unpack(data["OQP::SM"], 13)
    for spin, nocc in (("A", nalpha), ("B", nbeta)):
        vec = np.asarray(data[f"OQP::VEC_MO_{spin}"], dtype=float)[:nocc]
        dens = _unpack(data[f"OQP::DM_{spin}"], 13)
        captured = np.trace(vec @ smat @ dens @ smat @ vec.T)
        assert captured == pytest.approx(nocc, abs=1e-4), (spin, captured)


@pytest.mark.parametrize("converger,threads", [("diis", 1), ("diis", 4), ("soscf", 1)])
def test_uhf_triplet_level_shift_reaches_the_unshifted_state(tmp_path, converger, threads):
    common = dict(h2=H2_C2V, func="functional=pbe", guess="type=huckel", scf="uhf", mult=3)
    conv = "" if converger == "diis" else f"converger_type={converger}\n"
    plain, _ = _run(tmp_path, "plain", threads, extra=conv, **common)
    shifted, log = _run(tmp_path, "shifted", threads, extra=conv + "vshift=0.3", **common)
    assert RELEASED in log
    assert shifted["energy"] == pytest.approx(plain["energy"], abs=1e-6)
    _assert_aufbau_uhf(shifted, 6, 4)


@pytest.mark.parametrize("converger", ["diis", "soscf"])
def test_rohf_triplet_level_shift_reaches_the_unshifted_state(tmp_path, converger):
    common = dict(h2=H2_C2V, func="functional=pbe", guess="type=huckel", scf="rohf", mult=3)
    conv = "" if converger == "diis" else f"converger_type={converger}\n"
    plain, _ = _run(tmp_path, "plain", extra=conv, **common)
    shifted, log = _run(tmp_path, "shifted", extra=conv + "vshift=1.0", **common)
    assert RELEASED in log
    assert shifted["energy"] == pytest.approx(plain["energy"], abs=1e-6)


@pytest.mark.parametrize("converger,threads", [("diis", 1), ("diis", 4), ("soscf", 1)])
def test_rhf_level_shift_releases_a_supplied_excited_occupation(tmp_path, converger, threads):
    common = dict(h2=H2_C1, func="", scf="rhf", mult=1)
    ground, _ = _run(tmp_path, "ground", threads, guess="type=huckel", extra="", **common)
    excited, _ = _run(tmp_path, "excited", threads,
                      guess="type=huckel\nswapmo=5,6", extra="rstctmo=True", **common)
    assert excited["energy"] > ground["energy"] + 0.5        # HOMO^2 -> LUMO^2
    conv = "" if converger == "diis" else f"converger_type={converger}\n"
    shifted, log = _run(tmp_path, "shifted", threads,
                        guess=f"type=json\nfile={tmp_path / 'excited.json'}",
                        extra=conv + "vshift=0.3", **common)
    assert RELEASED in log
    assert shifted["energy"] == pytest.approx(ground["energy"], abs=1e-6)
