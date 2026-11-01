"""Native TRAH with the augmented Roothaan-Hall (ARH) model Hessian, and ARH
itself (``[scf] converger_type=arh``).

``[scf] trh_arh=n`` replaces the exact response Fock in the TRAH Hessian-vector
product by the ARH model built from the last n (density, Fock) pairs (Host et
al., J. Chem. Phys. 129, 124106 (2008); open-shell: Zhang and Yang,
arXiv:2606.03709).  The orbital gradient is unchanged, so the converged energy
must be the exact-Hessian one; what changes is the number of Fock builds.

Each case runs the same deck with the exact Hessian and with ARH and checks:

* both converge to the same energy (the model affects the path only);
* the ARH run needs no more response-Fock builds, and at most the two Fock
  builds of one rejected model step more in total;
* the Davidson micro-solve converges: the model is the symmetric secant form
  G (X^T G)^+ G^T.  The one-sided model F'P of the paper is not symmetric, and
  with it every macro-iteration ran to the micro-iteration limit (mean 37-40
  per macro-iteration measured on these and larger systems).  Its symmetric
  part (F'P + P F')/2 is indefinite with an XC kernel, and ARH alone then
  stalled on every DFT system tried (test_arh_converger_converges_dft_without_fallback);
* for Hartree-Fock, F is linear in D, so every stored pair satisfies
  F_i - F_ref = F'(D_ref).(D_i - D_ref) to screening accuracy.  A density
  convention mismatch between the stored densities and the trial densities of
  calc_h_op (RHF total vs spin density, a missing factor 2) makes this O(1).

RHF, UHF and ROHF Hartree-Fock and ROKS cover the three parameter layouts of
calc_h_op and the XC-kernel response.
"""

import os
import re
import subprocess
import sys

import pytest

H2O = """
O   0.000000000   0.000000000   0.000000000
H   0.000000000  -0.757000000   0.587000000
H   0.000000000   0.757000000   0.587000000"""

CH2 = """
C   0.000000000   0.000000000   0.105000000
H   0.000000000   0.990000000  -0.315000000
H   0.000000000  -0.990000000  -0.315000000"""

GRID = "dftgrid(rad_npts=50,ang_npts=194)"

CASES = {
    "rhf-hf": dict(geom=H2O, route="rhf/6-31g*", grid=""),
    "uhf-hf": dict(geom=H2O, route="uhf/6-31g* charge=1 mult=2", grid=""),
    "rohf-hf": dict(geom=H2O, route="rohf/6-31g* charge=1 mult=2", grid=""),
    "roks-bhhlyp": dict(geom=CH2, route="roks/bhhlyp/6-31g* mult=3", grid=GRID),
}


def _deck(case, converger="trah", arh=None, verbose=1, extra=""):
    """Canonical .oqp input for one case; arh=None leaves trh_arh at its default."""
    c = CASES[case]
    scf = f"converger_type={converger},conv=1e-07,maxit=200,verbose={verbose}"
    if arh is not None:
        scf += f",trh_arh={arh}"
    if extra:
        scf += "," + extra
    return f"{c['route']} guess(type=huckel) scf({scf}) {c['grid']}".rstrip(), c["geom"]


def _have_native_oqp():
    try:
        import oqp  # noqa: F401
    except Exception:
        return False
    return True


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(), reason="native OpenQP library not importable")

ENERGY = re.compile(r"Final \w+ energy is\s+(-?\d+\.\d+)")
BUILDS = re.compile(r"Native TRAH: Fock builds =\s*(\d+)\s+response Fock builds =\s*(\d+)")
MICRO = re.compile(r"\s(\d+)\s+(?:acc|rej)\s*$", re.M)
DIAG = re.compile(r"Native TRAH: diagnostic response Fock builds =\s*(\d+)")
LINEARITY = re.compile(r"ARH model:\s+\d+ pairs, linearity residual of the newest =\s*(\S+)")


def _run(tmp_path, name, deck):
    route, geom = deck
    xyz = tmp_path / f"{name}.xyz"
    atoms = geom.strip().splitlines()
    xyz.write_text(f"{len(atoms)}\n{name}\n" + "\n".join(atoms) + "\n")
    inp = tmp_path / f"{name}.oqp"
    inp.write_text(f'{route}\ngeom="{xyz.name}"\n')
    env = dict(os.environ, OMP_NUM_THREADS="2")
    proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", str(inp)],
                          capture_output=True, text=True, cwd=str(tmp_path),
                          env=env, timeout=900)
    log = inp.with_suffix(".log").read_text() if inp.with_suffix(".log").exists() else ""
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert "SCF convergence achieved" in log, log[-3000:]
    return log


@pytest.mark.parametrize("case", sorted(CASES))
def test_arh_hessian_same_energy_bounded_fock_builds(tmp_path, case):
    c = CASES[case]
    exact = _run(tmp_path, "exact", _deck(case, arh=0))
    arh = _run(tmp_path, "arh", _deck(case, arh=20, verbose=3))

    assert "Hessian: exact response" in exact
    assert "Hessian: augmented Roothaan-Hall model, memory =  20" in arh

    e_exact = float(ENERGY.findall(exact)[-1])
    e_arh = float(ENERGY.findall(arh)[-1])
    assert e_arh == pytest.approx(e_exact, abs=1e-8)

    # The model never costs response builds, and a model step that predicts
    # the energy poorly costs at most its trial energy and one gradient
    # refresh before the exact Hessian takes over.  These systems converge in
    # 4-6 macro-iterations, too few for the model to pay off everywhere: on
    # CH2 ROKS the only model step is the first one, made before any (D, F)
    # pair is stored, it is rejected and the rest follows the exact path
    # (14+39 exact vs 16+39 builds).  The speedup is a property measured by
    # the verification gate on larger open-shell systems, not asserted here.
    #
    # RHF Hartree-Fock water is exempt.  The secant model G (X^T G)^+ G^T is
    # the smallest positive-semidefinite completion of the stored span, so it
    # overestimates the curvature outside it (rho 1.3-1.5 on every water step)
    # and its linear convergence saves no exact build on this 8-step problem
    # (19+34 vs 8+24 builds).  RHF defaults to the exact Hessian for that reason.
    fock_x, resp_x = map(int, BUILDS.findall(exact)[-1])
    fock_a, resp_a = map(int, BUILDS.findall(arh)[-1])
    # The ARH run is at verbose=3 for the linearity residual below; each
    # residual costs one exact response build, reported apart from the
    # algorithm's counts so the two runs compare at equal footing.
    assert not DIAG.findall(exact)
    n_lin = len(LINEARITY.findall(arh))
    assert [int(n) for n in DIAG.findall(arh)] == [n_lin], (n_lin, DIAG.findall(arh))
    if case != "rhf-hf":
        assert resp_a <= resp_x, (fock_x, resp_x, fock_a, resp_a)
        assert fock_a + resp_a <= fock_x + resp_x + 2, (fock_x, resp_x, fock_a, resp_a)

    micro = [int(m) for m in MICRO.findall(arh)]
    assert micro and sum(micro) / len(micro) < 25, micro

    # Hartree-Fock: stored pairs obey G_i = F' X_i.  Agreement is ~1e-14 between
    # pairs from the same Fock-build mode; a pair straddling the switch from
    # incremental to full Fock builds carries their drift (up to 1.5e-5
    # measured).  A density-convention mismatch gives O(1).
    residuals = [float(r) for r in LINEARITY.findall(arh)]
    assert residuals, arh[-3000:]
    if not c["grid"]:
        assert min(residuals) < 1e-10, residuals
        assert max(residuals) < 1e-3, residuals


# Default (no trh_arh key): ARH for open shells (ROHF/ROKS, UHF/UKS); RHF
# keeps the exact Hessian, where the model saves no build.  Explicit
# trh_arh=0 / n>0 overrides this for any reference.
DEFAULT_HESSIAN = {
    "rhf-hf": "Hessian: exact response",
    "uhf-hf": "Hessian: augmented Roothaan-Hall model, memory =  20",
    "rohf-hf": "Hessian: augmented Roothaan-Hall model, memory =  20",
    "roks-bhhlyp": "Hessian: augmented Roothaan-Hall model, memory =  20",
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_arh_default_by_reference(tmp_path, case):
    deck = _deck(case)
    assert "trh_arh" not in deck[0]
    log = _run(tmp_path, "default", deck)
    assert DEFAULT_HESSIAN[case] in log, log[-3000:]


# [scf] converger_type=arh: ARH as a method of its own.  The model is used in
# every macro-iteration, the tight tail included, with one trust radius; a
# poorly predicted step shrinks that radius instead of being retried with the
# exact Hessian.  Like DIIS it does no stability check of its own (the model
# is positive semidefinite and cannot see a saddle), so no exact response
# product is formed at all; [scf] stability=true adds the TRAH safeguard pass.
def test_arh_converger_same_energy_no_exact_steps(tmp_path):
    for case in sorted(CASES):
        exact = _run(tmp_path, f"{case}-exact", _deck(case, arh=0))
        arh = _run(tmp_path, f"{case}-arh", _deck(case, converger="arh"))

        assert "Native ARH (augmented Roothaan-Hall trust region)" in arh, case
        assert ("Hessian: augmented Roothaan-Hall model in every step, memory =  20"
                " (no stability check; see [scf] stability)") in arh, case
        assert "model poor" not in arh, case

        e_exact = float(ENERGY.findall(exact)[-1])
        e_arh = float(ENERGY.findall(arh)[-1])
        assert e_arh == pytest.approx(e_exact, abs=1e-8), case

        _, resp_a = map(int, BUILDS.findall(arh)[-1])
        assert resp_a == 0, (case, resp_a)
        assert "unstable (Hess eig" not in arh, case
        assert "Verifying SCF stability" not in arh, case


# [scf] stability=true with converger_type=arh runs the post-SCF TRAH
# stability pass (exact Hessian), as it does after DIIS.  At this stable
# minimum the pass keeps the ARH solution.
def test_arh_converger_stability_opt_in_runs_trah_pass(tmp_path):
    log = _run(tmp_path, "uhf-hf-arh-stab", _deck("uhf-hf", converger="arh", extra="stability=true"))
    assert "Native ARH (augmented Roothaan-Hall trust region)" in log
    assert "PyOQP: Verifying SCF stability (TRAH)" in log, log[-3000:]
    # the safeguard pass is the exact TRAH, not ARH again
    assert "Native TRAH (trust-region Newton, augmented-Hessian Davidson)" in log
    assert "SCF point was unstable" not in log
    plain = _run(tmp_path, "uhf-hf-arh", _deck("uhf-hf", converger="arh"))
    assert float(ENERGY.findall(log)[-1]) == pytest.approx(
        float(ENERGY.findall(plain)[-1]), abs=1e-8)


SETTER_SCRIPT = """
from oqp.openqp import OPENQP

cfg = {"input.runtype": "energy", "input.method": "hf", "input.basis": "sto-3g",
       "input.charge": 0, "input.system": "\\n 1 0 0 0\\n 1 0 0 0.74",
       "scf.type": "rhf", "scf.multiplicity": 1}
op = OPENQP(cfg, True)
control = op.mol.data._data.control
seen = []
for name in ("arh", "trah", "arh", "diis", "arh", "soscf"):
    op.mol.data.set_scf_converger_type(name)
    seen.append((name, int(control.converger_type), bool(control.trh_arh_only)))
print(seen)
"""


def test_converger_name_sets_and_clears_arh_only(tmp_path):
    # Escalation and the stability safeguard switch the converger by name, so
    # leaving ARH must clear the flag the native driver reads.
    script = tmp_path / "setter.py"
    script.write_text(SETTER_SCRIPT)
    proc = subprocess.run([sys.executable, str(script)], capture_output=True,
                          text=True, cwd=str(tmp_path), timeout=300)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    seen = eval(proc.stdout.strip().splitlines()[-1])
    assert seen == [("arh", 2, True), ("trah", 2, False), ("arh", 2, True),
                    ("diis", 0, False), ("arh", 2, True), ("soscf", 1, False)]


# The symmetric part (F'P + P F')/2 of the one-sided ARH model keeps the
# coupling of the stored span to the rest of the space but not the block that
# bounds it; with an XC kernel this made the model Hessian indefinite and
# converger_type=arh collapsed its trust radius on benzene B3LYP (|g| 1.1e-4),
# FeCl4- ROKS, NO2 ROKS and Fe(H2O)6 ROKS, each finished only by the SOSCF
# fallback.  The escalation is disabled here so a stall fails the test.
DFT_ROUTE = ("rks/b3lypv1r/6-31g* guess(type=huckel) "
             "scf(converger_type=arh,escalation=diis,conv=1e-07,maxit=200) "
             "dftgrid(rad_npts=75,ang_npts=302)")

BENZENE = """
C   1.397000   0.000000   0.000000
C   0.698500   1.209800   0.000000
C  -0.698500   1.209800   0.000000
C  -1.397000   0.000000   0.000000
C  -0.698500  -1.209800   0.000000
C   0.698500  -1.209800   0.000000
H   2.481000   0.000000   0.000000
H   1.240500   2.148600   0.000000
H  -1.240500   2.148600   0.000000
H  -2.481000   0.000000   0.000000
H  -1.240500  -2.148600   0.000000
H   1.240500  -2.148600   0.000000"""


def test_arh_converger_converges_dft_without_fallback(tmp_path):
    log = _run(tmp_path, "benzene-arh", (DFT_ROUTE, BENZENE))
    assert "Native ARH (augmented Roothaan-Hall trust region)" in log
    assert "trust region collapsed" not in log, log[-3000:]
    assert "escalating" not in log, log[-3000:]



# The AO (density, Fock) history and the model response contraction are
# AO-indexed, so they are checked once more with spherical d functions
# (cc-pVDZ: 24 spherical vs 25 Cartesian functions for both molecules).  Exact
# Hessian, hybrid ARH and ARH alone must reach the same energy; for UHF
# Hartree-Fock the stored pairs must again satisfy G_i = F' X_i.
SPHERICAL = {
    "uhf-hf": ("uhf/cc-pvdz charge=1 mult=2", H2O, ""),
    "roks-bhhlyp": ("roks/bhhlyp/cc-pvdz mult=3", CH2, GRID),
}
SPHERICAL_AO = re.compile(r"Spherical AO functions =\s*(\d+)\s+Cartesian-equivalent =\s*(\d+)")


@pytest.mark.parametrize("case", sorted(SPHERICAL))
def test_arh_spherical_basis_same_energy(tmp_path, case):
    route, geom, grid = SPHERICAL[case]

    def deck(scf):
        return (f"{route} guess(type=huckel) "
                f"scf(conv=1e-07,maxit=200,{scf}) {grid}").rstrip(), geom

    exact = _run(tmp_path, "exact", deck("converger_type=trah,trh_arh=0"))
    hybrid = _run(tmp_path, "hybrid", deck("converger_type=trah,trh_arh=20,verbose=3"))
    alone = _run(tmp_path, "alone", deck("converger_type=arh,escalation=diis"))

    for log in (exact, hybrid, alone):
        assert "AO angular type: spherical harmonics" in log, log[:3000]
        assert SPHERICAL_AO.findall(log)[-1] == ("24", "25"), log[:3000]
    assert "Hessian: augmented Roothaan-Hall model, memory =  20" in hybrid
    assert "Native ARH (augmented Roothaan-Hall trust region)" in alone

    e_exact = float(ENERGY.findall(exact)[-1])
    assert float(ENERGY.findall(hybrid)[-1]) == pytest.approx(e_exact, abs=1e-8)
    assert float(ENERGY.findall(alone)[-1]) == pytest.approx(e_exact, abs=1e-8)
    assert "escalating to" not in alone

    residuals = [float(r) for r in LINEARITY.findall(hybrid)]
    assert residuals, hybrid[-3000:]
    if not grid:
        assert min(residuals) < 1e-10, residuals
        assert max(residuals) < 1e-3, residuals
