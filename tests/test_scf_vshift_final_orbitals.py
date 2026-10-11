"""Published SCF orbital energies must not carry the level shift.

The DIIS loop tests convergence before it diagonalises the new Fock, so the
orbitals it publishes are the eigenpairs of the previous iteration's (DIIS)
Fock.  With a level shift the loop converged once with the shift, printed
"Performing a last SCF with zero VSHIFT.", built one more Fock without it and
exited on that build: the orbitals in hand still came from the shifted Fock,
so every virtual energy was raised by vshift (RHF) and, for ROHF, the open
shells by vshift/2 and the virtuals by vshift.  After v-DIIS switches the
shift off, shifted Focks stay in the DIIS subspace.  Every consumer of the
orbital energies (TDHF/MRSF A-matrix diagonals, MP2, Koopmans, the printed
orbitals) inherited the shift.  The converged orbitals are kept (they built
the converged density and energy); only their energies are refreshed against
the final unshifted Fock.

H2/STO-3G, R = 1.4 bohr, vshift = 0.1, before the fix:

    RHF   -0.5782030162  0.7702678589   (correct  0.6702678589)
    ROHF  -0.2925650159  0.4960804870   (unshifted -0.3425650159 0.4460804870)

The RHF reference is closed form (Szabo & Ostlund, Appendix A integrals; see
tests/test_tdhf_single_pair_davidson.py).  Every other case compares a shifted
run with the same run at vshift = 0: a level shift may change the path to the
converged density, never the converged orbital energies.  The H2O cases
converge in many iterations, so the DIIS subspace holds several shifted
Focks.  Every shifted run must also store the density built from the stored
occupied orbitals, i.e. publish one consistent SCF state.  The two v-DIIS
cases guard that route; they also pass without this fix on a base that
already handles the iteration-1 v-DIIS exit (fock_noshift).
"""

import json
import os
import re
import subprocess
import sys

import pytest

H2_RHF_EPS = (-0.5782030162, 0.6702678589)

H2 = """\
[input]
system=
   1   0.000000000   0.000000000   0.000000000
   1   0.000000000   0.000000000   0.740848000
charge=0
runtype=energy
basis=sto-3g
method=hf

[guess]
type={guess}
save_mol=True

[scf]
multiplicity={mult}
type={scf}
conv=1.0e-10
{extra}
"""

# Distorted water, no symmetry; 13 basis functions in 6-31G.
H2O = """\
[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.633194329  -0.483194329  -0.554469223
charge={charge}
runtype=energy
basis=6-31g
method=hf

[guess]
type=huckel
save_mol=True

[scf]
multiplicity={mult}
type={scf}
conv=1.0e-10
{extra}
"""

LAST_SCF = "Performing a last SCF with zero VSHIFT."
NUM = re.compile(r"^-?\d+\.\d+$")


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
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert "SCF convergence achieved" in log, log[-3000:]
    data = json.loads((tmp_path / f"{name}.json").read_text())
    return log, data


def _alpha_energies(log):
    """All alpha orbital energies of the last printed orbital table."""
    head = log.rfind("Alpha Orbitals")
    assert head >= 0, log[-3000:]
    eps = []
    for line in log[head:].splitlines()[1:]:
        tok = line.split()
        if not tok:
            continue
        if all(t.isdigit() for t in tok):
            continue
        if all(NUM.match(t) for t in tok):
            eps.extend(float(t) for t in tok)
            continue
        break
    assert eps, log[head:head + 2000]
    return eps


def _assert_density_from_orbitals(data, nocc_a, nocc_b, rhf):
    """Stored densities equal those built from the stored occupied orbitals."""
    vec = data["OQP::VEC_MO_A"]                      # rows are MOs
    nbf = len(vec[0])

    def packed(nocc, factor):
        return [factor * sum(vec[k][i] * vec[k][j] for k in range(nocc))
                for i in range(nbf) for j in range(i + 1)]

    if rhf:
        assert data["OQP::DM_A"] == pytest.approx(packed(nocc_a, 2.0), abs=1e-10)
    else:
        assert data["OQP::DM_A"] == pytest.approx(packed(nocc_a, 1.0), abs=1e-10)
        assert data["OQP::DM_B"] == pytest.approx(packed(nocc_b, 1.0), abs=1e-10)


def _pair(tmp_path, name, template, shift_line, threads=1, **fmt):
    """Run a deck with and without a level shift."""
    shifted = _run(tmp_path, f"{name}_shift",
                   template.format(extra=shift_line, **fmt), threads)
    plain = _run(tmp_path, f"{name}_plain",
                 template.format(extra="", **fmt), threads)
    return shifted, plain


@pytest.mark.parametrize("guess", ["huckel", "hcore"])
def test_h2_rhf_vshift_matches_closed_form(tmp_path, guess):
    (log_s, data_s), (log_p, _) = _pair(tmp_path, f"h2_rhf_{guess}", H2,
                                        "vshift=0.1", guess=guess,
                                        mult=1, scf="rhf")
    # The shifted run must really have gone through the zero-shift step.
    assert LAST_SCF in log_s
    assert _alpha_energies(log_s) == pytest.approx(H2_RHF_EPS, abs=1e-8)
    assert _alpha_energies(log_p) == pytest.approx(H2_RHF_EPS, abs=1e-8)
    _assert_density_from_orbitals(data_s, 1, 1, rhf=True)


def test_h2_rhf_vdiis_matches_closed_form(tmp_path):
    # v-DIIS sets vshift = 0.1 and zeroes it only after the first build, whose
    # DIIS error is already below the switch: no "last SCF" step is printed.
    log, data = _run(tmp_path, "h2_vdiis",
                     H2.format(guess="huckel", mult=1, scf="rhf",
                               extra="diis_type=vdiis"))
    assert LAST_SCF not in log
    assert _alpha_energies(log) == pytest.approx(H2_RHF_EPS, abs=1e-8)
    _assert_density_from_orbitals(data, 1, 1, rhf=True)


def test_h2_rohf_triplet_vshift_matches_unshifted(tmp_path):
    (log_s, data_s), (log_p, _) = _pair(tmp_path, "h2_rohf", H2, "vshift=0.1",
                                        guess="huckel", mult=3, scf="rohf")
    assert LAST_SCF in log_s
    eps_s, eps_p = _alpha_energies(log_s), _alpha_energies(log_p)
    assert len(eps_s) == len(eps_p) == 2
    assert eps_s == pytest.approx(eps_p, abs=1e-8)
    _assert_density_from_orbitals(data_s, 2, 0, rhf=False)


@pytest.mark.parametrize("threads", [1, 4])
def test_h2o_rhf_vshift_matches_unshifted(tmp_path, threads):
    (log_s, data_s), (log_p, _) = _pair(tmp_path, f"h2o_rhf_t{threads}", H2O,
                                        "vshift=0.3", threads,
                                        charge=0, mult=1, scf="rhf")
    assert LAST_SCF in log_s
    eps_s, eps_p = _alpha_energies(log_s), _alpha_energies(log_p)
    assert len(eps_s) == len(eps_p) == 13
    assert eps_s == pytest.approx(eps_p, abs=1e-6)
    _assert_density_from_orbitals(data_s, 5, 5, rhf=True)


@pytest.mark.parametrize("threads", [1, 4])
def test_h2o_cation_rohf_vshift_matches_unshifted(tmp_path, threads):
    (log_s, data_s), (log_p, _) = _pair(tmp_path, f"h2op_rohf_t{threads}", H2O,
                                        "vshift=0.3", threads,
                                        charge=1, mult=2, scf="rohf")
    assert LAST_SCF in log_s
    eps_s, eps_p = _alpha_energies(log_s), _alpha_energies(log_p)
    assert len(eps_s) == len(eps_p) == 13
    assert eps_s == pytest.approx(eps_p, abs=1e-6)
    _assert_density_from_orbitals(data_s, 5, 4, rhf=False)


def test_h2o_rhf_vdiis_matches_cdiis(tmp_path):
    # v-DIIS switches the shift off inside the loop after several shifted
    # builds, which stay in the DIIS subspace.
    (log_s, data_s), (log_p, _) = _pair(tmp_path, "h2o_vdiis", H2O,
                                        "diis_type=vdiis",
                                        charge=0, mult=1, scf="rhf")
    assert LAST_SCF not in log_s
    eps_s, eps_p = _alpha_energies(log_s), _alpha_energies(log_p)
    assert len(eps_s) == len(eps_p) == 13
    assert eps_s == pytest.approx(eps_p, abs=1e-6)
    _assert_density_from_orbitals(data_s, 5, 5, rhf=True)
