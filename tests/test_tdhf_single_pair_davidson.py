"""TDA/RPA on a one-pair response space: SCF orbital energies, Davidson, exit status.

H2/STO-3G at R = 1.4 bohr has one occupied and one virtual orbital, so the
TDHF response space is a single sigma_g -> sigma_u pair (nocc*nvir = 1).
The subspace is then the whole space from the first iteration and the Ritz
value of the 1x1 reduced matrix is the exact answer.

Before the fix, ``rpaeig`` took the TDA Ritz value of a 1x1 reduced matrix
from its *output* eigenvector array instead of from A, so every TDA run with
one trial vector reported E = 0.  On H2 the subspace cannot grow
(nvec == mxvec), the Davidson stopped "not converged", and pyoqp then called a
bare ``exit()`` -- status 0 -- so the failure passed as a successful run.

Even RPA, which converged, was wrong (22.41 eV instead of 25.30 eV): the guess
density of H2/STO-3G is exact by symmetry, the DIIS loop converged on its first
Fock build before diagonalising any SCF Fock, and the published orbital
energies were the Huckel guess eigenvalues (-1.0721, +0.0721 Eh; Hcore gave
-1.2528, -0.4756 Eh and 11.98 eV).  The excitation energy depended on the
guess.  The tests below run both guesses.

Reference values are closed-form, not taken from another program.  With the
analytic s-Gaussian integrals of BSE STO-3G hydrogen (Szabo & Ostlund,
Appendix A; their Sec. 3.5.2 gives eps_u - eps_g = 1.2485, J_gu = 0.6636,
K_gu = 0.1813 Eh) and the single pair g -> u:

    singlet  A = de - J + 2K,  B =  K
    TDA  omega = A
    RPA  omega = sqrt((A - B)(A + B))

Evaluated to 1e-10 from the analytic integrals: eps = (-0.57820302,
0.67026784), J_gu = 0.663564, K_gu = 0.181258, E_SCF = -1.1167143238 Eh.
"""

import json
import os
import re
import subprocess
import sys

import pytest

HARTREE_EV = 27.211386245988  # CODATA 2018, as source/physical_constants.F90

# Singlet omega (Eh) for H2/STO-3G, R = 0.740848 Angstrom = 1.4 bohr.
# Closed-shell triplet RPA/TDA is refused by the input checker on this base.
REFERENCE = {
    "tda": 0.9474226807,
    "rpa": 0.9299222055,
}
EPS_REFERENCE = (-0.57820302, 0.67026784)

H2 = """\
[input]
system=
   1   0.000000000   0.000000000   0.000000000
   1   0.000000000   0.000000000   0.740848000
charge=0
runtype=energy
basis=sto-3g
method=tdhf

[guess]
type={guess}

[scf]
multiplicity=1
type=rhf
conv=1.0e-10

[tdhf]
type={td}
multiplicity=1
nstate={nstate}
conv=1.0e-10
"""

# Distorted water (no symmetry, so no root is alone in a symmetry block and
# the unit start vector cannot be an exact eigenvector) in 6-31G:
# nocc*nvir = 5*8 = 40.  maxit=1 with an unreachable tolerance cannot
# converge, even after the driver's two automatic restarts (maxit 1 -> 2 -> 4).
H2O_UNCONVERGED = """\
[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.633194329  -0.483194329  -0.554469223
charge=0
runtype=energy
basis=6-31g
method=tdhf

[guess]
type=huckel

[scf]
multiplicity=1
type=rhf

[tdhf]
type=tda
nstate=1
maxit=1
conv=1.0e-30
"""

# Triplet ROHF reference for MRSF.  Both STO-3G orbitals are singly occupied,
# so the density is fixed by symmetry and every guess converges on the first
# Fock build -- the case where the guess orbital energies used to survive.
H2_ROHF_MRSF = """\
[input]
system=
   1   0.000000000   0.000000000   0.000000000
   1   0.000000000   0.000000000   0.740848000
charge=0
runtype=energy
basis=sto-3g
method=tdhf

[guess]
type={guess}

[scf]
multiplicity=3
type=rohf
conv=1.0e-10

[tdhf]
type=mrsf
multiplicity=1
nstate=1
conv=1.0e-10
"""

# One-electron UHF reference (H2+ doublet): the density is fixed by symmetry
# and the alpha virtual is empty, so a level shift on it is visible.
H2P_UHF = """\
[input]
system=
   1   0.000000000   0.000000000   0.000000000
   1   0.000000000   0.000000000   0.740848000
charge=1
runtype=energy
basis=sto-3g
method=hf

[guess]
type={guess}

[scf]
multiplicity=2
type=uhf
conv=1.0e-10
"""

SUMMARY_ROOT = re.compile(r"^\s*0 -> 1\s+(-?\d+\.\d+)", re.MULTILINE)
ALPHA_EPS = re.compile(
    r"Alpha Orbitals -+\s*\n\s*1\s+2\s*\n\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)")
# Davidson progress lines; the last occurrence of each root is its final value.
DAVIDSON_STATE = re.compile(r"State\s+(\d+)\s+E =\s+(-?\d+\.\d+) eV")


def _have_native_oqp():
    try:
        import oqp  # noqa: F401
    except Exception:
        return False
    import oqp
    return hasattr(oqp, "tdhf_energy")


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(),
    reason="native OpenQP library with tdhf_energy not importable",
)


SCF_ITER = re.compile(
    r"^\s+\d+\s+-?\d+\.\d+E[+-]\d+\s+-?\d+\.\d+E[+-]\d+\s+\d+\s+(\d+\.\d+E[+-]\d+)",
    re.MULTILINE)


def _diis_errors(log):
    return [float(e) for e in SCF_ITER.findall(log)]


def _run(tmp_path, name, deck, threads=1):
    inp = tmp_path / f"{name}.inp"
    inp.write_text(deck)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    proc = subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", str(inp)],
        capture_output=True, text=True, cwd=str(tmp_path), env=env, timeout=600,
    )
    log_path = inp.with_suffix(".log")
    log = log_path.read_text() if log_path.exists() else ""
    return proc, log


@pytest.mark.parametrize("guess", ["huckel", "hcore"])
def test_one_iteration_scf_publishes_fock_eigenvalues(tmp_path, guess):
    proc, log = _run(tmp_path, f"h2_scf_{guess}",
                     H2.format(td="rpa", guess=guess, nstate=1))
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    # The case this guards: the guess density is exact, so the first Fock
    # build already meets the tolerance.  A model-Hamiltonian guess is not
    # accepted there; the one extra iteration diagonalises the SCF Fock.
    assert _diis_errors(log)[0] < 1e-10, log[-3000:]
    assert re.search(r"after\s+2 iterations", log), log[-3000:]
    eps = ALPHA_EPS.findall(log)
    assert eps, log[-3000:]
    assert [float(e) for e in eps[-1]] == pytest.approx(EPS_REFERENCE, abs=1e-7)


def test_one_iteration_scf_refreshes_supplied_orbital_energies(tmp_path):
    """Supplied orbitals keep coefficients and phases; their energies are refreshed.

    The supplied JSON carries the converged H2 orbitals with deliberately wrong
    orbital energies and the virtual orbital's sign flipped.  The SCF converges
    on its first Fock build, so before the fix the wrong energies were published
    and fed the response; the orbitals themselves must not be rotated or
    re-phased (restart, MD and NAMD phase tracking rely on that).
    """
    first = H2.format(td="rpa", guess="huckel", nstate=1).replace(
        "type=huckel\n", "type=huckel\nsave_mol=True\n")
    proc, _ = _run(tmp_path, "h2_first", first)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    data = json.loads((tmp_path / "h2_first.json").read_text())
    vec = [list(row) for row in data["OQP::VEC_MO_A"]]
    vec[1] = [-c for c in vec[1]]
    data["OQP::VEC_MO_A"] = vec
    data["OQP::E_MO_A"] = [9.0, 11.0]
    (tmp_path / "supplied.json").write_text(json.dumps(data))

    second = H2.format(td="rpa", guess="json", nstate=1).replace(
        "type=json\n", "type=json\nfile=supplied.json\nsave_mol=True\n")
    proc, log = _run(tmp_path, "h2_supplied", second)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert re.search(r"after\s+1 iterations", log), log[-3000:]
    out = json.loads((tmp_path / "h2_supplied.json").read_text())
    assert out["OQP::E_MO_A"] == pytest.approx(EPS_REFERENCE, abs=1e-7)
    flat = [c for row in out["OQP::VEC_MO_A"] for c in row]
    assert flat == pytest.approx([c for row in vec for c in row], abs=1e-12)
    omega = float(SUMMARY_ROOT.findall(log)[-1]) / HARTREE_EV
    assert omega == pytest.approx(REFERENCE["rpa"], abs=1e-6)


def _occupied_density(vec, nocc, factor):
    nbf = len(vec[0])
    return [factor * sum(vec[k][i] * vec[k][j] for k in range(nocc))
            for i in range(nbf) for j in range(i + 1)]


# scf type -> (multiplicity, [(density key, orbital key, occupied, factor)])
LOOSE_REFS = {
    "rhf": (1, [("OQP::DM_A", "OQP::VEC_MO_A", 5, 2.0)]),
    "uhf": (1, [("OQP::DM_A", "OQP::VEC_MO_A", 5, 1.0),
                ("OQP::DM_B", "OQP::VEC_MO_B", 5, 1.0)]),
    "rohf": (3, [("OQP::DM_A", "OQP::VEC_MO_A", 6, 1.0),
                 ("OQP::DM_B", "OQP::VEC_MO_A", 4, 1.0)]),
}


@pytest.mark.parametrize("converger", ["diis", "soscf", "trah"])
@pytest.mark.parametrize("ref", ["rhf", "uhf", "rohf"])
def test_loose_convergence_keeps_orbitals_and_density_together(tmp_path, ref, converger):
    """A guess that only meets a loose tolerance must not split the SCF state.

    Water/STO-3G with a Huckel guess and conv=1.0: the first Fock build (DIIS
    error 0.73 for RHF) passes the tolerance although the guess is not
    self-consistent, and the run exits on an early iteration.  Whatever is
    published, the stored density must be the one built from the stored
    occupied orbitals.  Diagonalising the Fock at the exit replaced the
    orbitals while keeping the density: DIIS on iteration 1, the SOSCF/TRAH
    exit handler and the UHF publication after the loop on any iteration.
    """
    mult, checks = LOOSE_REFS[ref]
    extra = "" if converger == "diis" else f"converger_type={converger}\n"
    deck = H2O_UNCONVERGED.replace("basis=6-31g", "basis=sto-3g").replace(
        "type=huckel\n", "type=huckel\nsave_mol=True\n").replace(
        "[scf]\nmultiplicity=1\ntype=rhf\n",
        f"[scf]\nmultiplicity={mult}\ntype={ref}\nconv=1.0\n" + extra)
    deck = deck.replace("method=tdhf", "method=hf")
    deck = deck[:deck.index("[tdhf]")]
    proc, log = _run(tmp_path, "h2o_loose", deck)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    if ref == "rhf" and converger != "trah":
        assert _diis_errors(log)[0] < 1.0, log[-3000:]
    data = json.loads((tmp_path / "h2o_loose.json").read_text())
    for dkey, vkey, nocc, factor in checks:
        packed = _occupied_density(data[vkey], nocc, factor)     # rows are MOs
        assert data[dkey] == pytest.approx(packed, abs=1e-10), (dkey, vkey)


@pytest.mark.parametrize("route", ["cold", "supplied"])
@pytest.mark.parametrize("ref", ["rhf", "rohf", "uhf"])
def test_level_shift_never_reaches_published_orbital_energies(tmp_path, ref, route):
    """V-DIIS adds its level shift to the first Fock before switching it off.

    The first build already meets the tolerance, so the published orbital
    energies come either from the extra iteration a cold guess takes or from
    the refresh of supplied orbitals.  Both must use the unshifted Fock: the
    same orbital energies as a run without any shift (a shifted virtual would
    read 0.7703 instead of 0.6703 Eh for RHF H2).
    """
    deck = {"rhf": H2.format(td="rpa", guess="huckel", nstate=1),
            "rohf": H2_ROHF_MRSF.format(guess="huckel"),
            "uhf": H2P_UHF.format(guess="huckel")}[ref]
    plain = deck.replace("type=huckel\n", "type=huckel\nsave_mol=True\n")
    proc, log = _run(tmp_path, f"{ref}_plain", plain)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    eps_plain = [float(e) for e in ALPHA_EPS.findall(log)[-1]]

    shifted = deck.replace(f"type={ref}\nconv=1.0e-10\n",
                           f"type={ref}\nconv=1.0e-10\ndiis_type=vdiis\nvshift=0.1\n")
    assert shifted != deck
    if route == "supplied":
        shifted = shifted.replace(
            "type=huckel\n", f"type=json\nfile={ref}_plain.json\n")
    proc, log = _run(tmp_path, f"{ref}_{route}_vdiis", shifted)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert _diis_errors(log)[0] < 1e-10, log[-3000:]
    expected_iters = 1 if route == "supplied" else 2
    assert re.search(rf"after\s+{expected_iters} iterations", log), log[-3000:]
    eps = [float(e) for e in ALPHA_EPS.findall(log)[-1]]
    assert eps == pytest.approx(eps_plain, abs=1e-8)


def test_cold_guess_is_not_accepted_without_an_scf_iteration(tmp_path):
    """With [scf] maxit=1 there is no iteration left to diagonalise an SCF Fock.

    The exact H2 guess density meets the tolerance on the first build, but its
    orbitals are the guess Hamiltonian's.  That SCF must not be reported as
    converged; if the run recovers (the TRAH escalation), whatever it finally
    publishes must be the Fock eigenvalues, never the guess eigenvalues.
    """
    deck = H2.format(td="rpa", guess="huckel", nstate=1).replace(
        "type=rhf\nconv=1.0e-10\n", "type=rhf\nconv=1.0e-10\nmaxit=1\n")
    proc, log = _run(tmp_path, "h2_maxit1", deck)
    assert "SCF did not converge" in log, log[-3000:]
    eps = ALPHA_EPS.findall(log)
    assert eps, log[-3000:]
    assert [float(e) for e in eps[-1]] != pytest.approx([-1.0720752606, 0.0720752606], abs=1e-6)
    if proc.returncode == 0:
        assert [float(e) for e in eps[-1]] == pytest.approx(EPS_REFERENCE, abs=1e-7)


def test_one_iteration_rohf_is_guess_independent(tmp_path):
    results = {}
    for guess in ("huckel", "hcore"):
        proc, log = _run(tmp_path, f"h2_rohf_{guess}",
                         H2_ROHF_MRSF.format(guess=guess))
        assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
        assert _diis_errors(log)[0] < 1e-10, log[-3000:]
        assert re.search(r"ROHF energy is\s+\S+\s+after\s+2 iterations", log), \
            log[-3000:]
        eps = [float(e) for e in ALPHA_EPS.findall(log)[-1]]
        assert "energies converged" in log, log[-3000:]
        states = {int(n): float(e) for n, e in DAVIDSON_STATE.findall(log)}
        assert states, log[-3000:]
        results[guess] = (eps, states)
    (eps_h, st_h), (eps_c, st_c) = results["huckel"], results["hcore"]
    # The ROHF orbitals of this reference are the same whatever the guess, so
    # their energies, and every MRSF state built on them, must be too.
    assert eps_h == pytest.approx(eps_c, abs=1e-8)
    assert sorted(st_h) == sorted(st_c)
    for n in st_h:
        assert st_h[n] == pytest.approx(st_c[n], abs=1e-6)  # eV, printed to 1e-6


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("guess", ["huckel", "hcore"])
@pytest.mark.parametrize("td", sorted(REFERENCE))
def test_h2_single_pair_matches_closed_form(tmp_path, td, guess, threads):
    proc, log = _run(tmp_path, f"h2_{td}_{guess}",
                     H2.format(td=td, guess=guess, nstate=1), threads)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    # The Davidson itself must report convergence, not just the process.
    assert "TD-DFT energies converged" in log
    assert "NOT CONVERGED" not in log
    assert "auto-restart" not in log
    roots = SUMMARY_ROOT.findall(log)
    assert roots, log[-2000:]
    omega = float(roots[-1]) / HARTREE_EV
    # The summary prints 6 decimals in eV (2e-8 Eh); the reference is exact.
    assert omega == pytest.approx(REFERENCE[td], abs=1e-6)


def test_nstate_above_response_dimension_is_refused(tmp_path):
    proc, log = _run(tmp_path, "h2_nstate2",
                     H2.format(td="tda", guess="huckel", nstate=2))
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr + log
    assert "exceeds the response-space dimension" in out


def test_unconverged_davidson_exits_nonzero(tmp_path):
    proc, log = _run(tmp_path, "h2o_unconverged", H2O_UNCONVERGED)
    # The Davidson must really have failed, otherwise this tests nothing ...
    assert "NOT CONVERGED" in log, log[-2000:]
    assert "TD-DFT energies converged" not in log
    # ... and that failure must not leave the process looking successful.
    assert proc.returncode != 0
