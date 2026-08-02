"""``[casscf] converger = trah``: the shared trust-region core over a
matrix-free CASSCF Hessian-vector product.

Three contracts, in order of how much they are worth:

  (a) the Hessian-vector product IS the Hessian.  ``_trah_hess_vec`` never
      assembles anything; this compares it against the symmetrized
      finite-difference Hessian ``_fd_orbital_hessian`` builds column by column,
      and separately checks that the closed-form Baker-Campbell-Hausdorff
      correction it subtracts is exactly the antisymmetric part of ``dg/dkappa``
      (coefficient 1/4, fitted, not assumed);
  (b) the native driver and the NumPy reference reach the same energy;
  (c) ``trah`` reproduces the CASSCF energy of the existing convergers.

Only (a) is a unit test of new numerics.  The rest is the same
"every converger finds the same solution" contract
``tests/test_casscf_convergers.py`` already applies to ``ah``/``diis``/``auto``.
"""
import re

import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False)) and hasattr(oqp, "fci_ao_integrals")


_H2O = "\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692"
_LIH = "\nLi 0 0 0\nH 0 0 3.0"
_CAS44 = {"active_electrons": "4", "active_orbitals": "4", "frozen_core": "3",
          "max_det": "5000"}
_CAS22 = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1",
          "max_det": "5000"}
_TOL = 1.0e-8


def _run(workdir, project, *, system, basis, cas, converger=None, hessian=None,
         native=True, maxmacro="60"):
    from oqp.library import casscf as cs
    from oqp.pyoqp import Runner

    casscf = {
        "max_macro_iterations": maxmacro, "optimizer": "newton", "root": "0",
        "gradient_norm_tol": "1.0e-7", "energy_decrease_tol": "1.0e-10",
        "step_norm_tol": "1.0e-8", "max_rotation_norm": "2.0e-1",
        "level_shift": "1.0e-3", "canonicalize": "true",
    }
    if converger:
        casscf["converger"] = converger
    if hessian:
        casscf["hessian"] = hessian
    input_dict = {
        "input": {"system": system, "charge": "0", "basis": basis,
                  "method": "casscf", "runtype": "energy"},
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": "1", "maxit": "100",
                "forced_attempt": "3", "save_molden": "False"},
        "cas": dict(cas),
        "ci": {"nroot": "1", "solver": "dense", "eig_tol": "1.0e-10",
               "integral_backend": "native", "target_spin": "any",
               "root_tracking": "energy"},
        "casscf": casscf,
        "tests": {"exception": "True"},
    }
    saved = cs._casscf_energy_backend
    if not native:
        cs._casscf_energy_backend = lambda: None
    try:
        runner = Runner(project=project, input_file=None,
                        log=str(workdir / f"{project}.log"), silent=1,
                        input_dict=input_dict)
        runner.run()
        return float(runner.mol.energies[0])
    finally:
        cs._casscf_energy_backend = saved


# ------------------------------------------------------- (a) the product IS H
def _hessian_probe(workdir, system, basis, cas, checks):
    """Run one CASSCF through the NumPy trah converger, capturing the Hessian
    comparison at the first macroiteration."""
    from oqp.library import casscf as cs
    from oqp.library import casscf_convergers as cv

    real = cv._trah_optimize

    def probe(C, raw_evaluate, evaluate, pairs, nbf, options, params, trace):
        if not checks:
            _obj, grad, _e, _c = evaluate(C)
            gmat = cv._trah_gradient_matrix(raw_evaluate.state["F"])
            npar = len(pairs)
            hh = 1.0e-4
            J = np.zeros((npar, npar))
            for k in range(npar):
                e = np.zeros(npar)
                e[k] = hh
                _, gp, _, _ = evaluate(cs._orbital_rotate(C, e, pairs, nbf))
                _, gm, _, _ = evaluate(cs._orbital_rotate(C, -e, pairs, nbf))
                J[:, k] = (gp - gm) / (2.0 * hh)
            H = 0.5 * (J + J.T)
            rng = np.random.RandomState(7)
            for _ in range(3):
                v = rng.standard_normal(npar)
                v /= np.linalg.norm(v)
                hv = cv._trah_hess_vec(C, evaluate, pairs, nbf, v, gmat)
                ref = H @ v
                K = cs._kappa_matrix(v, pairs, nbf)
                R = K @ gmat - gmat @ K
                corr = np.array([R[q, p] - R[p, q] for (p, q) in pairs])
                checks.append({
                    "rel": float(np.linalg.norm(hv - ref) / np.linalg.norm(ref)),
                    "uncorrected": float(np.linalg.norm(J @ v - ref)
                                         / np.linalg.norm(ref)),
                    # the coefficient of the BCH correction, fitted rather than
                    # assumed: (J - H).v projected onto the correction vector
                    "coeff": float(((J @ v - ref) @ corr) / (corr @ corr)),
                    "gnorm": float(np.linalg.norm(grad)),
                })
        return real(C, raw_evaluate, evaluate, pairs, nbf, options, params, trace)

    cv._trah_optimize = probe
    try:
        _run(workdir, "probe", system=system, basis=basis, cas=cas,
             converger="trah", native=False)
    finally:
        cv._trah_optimize = real


def test_hess_vec_reproduces_the_assembled_hessian(tmp_path):
    """H.v from two CI solves == the symmetrized 2*n_par-solve FD Hessian . v.

    The point of the whole converger: the same curvature information without
    ever building the matrix.  n_par = 54 here, so the reference costs 108
    gradient evaluations and each product costs 2."""
    if not _backend_available():
        pytest.skip("native OQP backend not built")
    checks = []
    _hessian_probe(tmp_path, _H2O, "6-31g", _CAS44, checks)
    assert checks, "the trah converger never ran"
    for c in checks:
        assert c["rel"] < 1.0e-6, c


def test_bch_correction_coefficient_is_one_quarter(tmp_path):
    """``dg/dkappa`` is NOT the Hessian, and the difference is not noise.

    ``g`` is the gradient in the displaced frame, so ``J = dg/dkappa`` differs
    from the second derivative by the BCH term of
    ``exp(K(kappa)) exp(K(delta))`` -- exactly antisymmetric, and of order |g|.
    That term is what ``_fd_orbital_hessian`` removes by symmetrizing; the
    matrix-free product removes it in closed form instead.  This asserts both
    halves: that the uncorrected product is measurably wrong, and that the
    correction's coefficient is 1/4 (fitted from the data, so a sign or factor
    slip cannot pass)."""
    if not _backend_available():
        pytest.skip("native OQP backend not built")
    checks = []
    _hessian_probe(tmp_path, _H2O, "6-31g", _CAS44, checks)
    assert checks
    for c in checks:
        # the correction is real: dropping it costs ~4 orders of magnitude
        assert c["uncorrected"] > 100.0 * c["rel"], c
        assert c["coeff"] == pytest.approx(0.25, abs=1.0e-5), c


def test_hess_vec_correction_vanishes_without_virtual_orbitals(tmp_path):
    """H2O/STO-3G CAS(4,4) has no virtual orbitals, so the BCH commutator has
    no non-redundant component with a non-zero gradient and the uncorrected
    product is already the Hessian.  Recorded because it is why the correction
    is invisible on the smallest CASSCF example in the tree."""
    if not _backend_available():
        pytest.skip("native OQP backend not built")
    checks = []
    _hessian_probe(tmp_path, _H2O, "sto-3g", _CAS44, checks)
    assert checks
    for c in checks:
        assert c["uncorrected"] < 1.0e-6, c


# ------------------------------------------------- (b) native vs NumPy pin
@pytest.mark.parametrize("system,basis,cas", [
    (_H2O, "6-31g", _CAS44),
    (_LIH, "6-31g", _CAS22),
])
def test_native_trah_matches_the_numpy_reference(tmp_path, system, basis, cas):
    if not _backend_available():
        pytest.skip("native OQP backend not built")
    native = _run(tmp_path, "trah_native", system=system, basis=basis, cas=cas,
                  converger="trah", native=True)
    numpy = _run(tmp_path, "trah_numpy", system=system, basis=basis, cas=cas,
                 converger="trah", native=False)
    # Validate on ENERGY, never on iteration count: a trust-region trajectory
    # is not reproducible run to run at this tolerance.
    assert native == pytest.approx(numpy, abs=_TOL)


# ------------------------------------------- (c) same solution as the others
@pytest.mark.parametrize("hessian", [None, "analytic"])
def test_trah_reproduces_the_default_converger(tmp_path, hessian):
    """Matrix-free (``hessian=fd``, the default) and assembled
    (``hessian=analytic``) both land on the two-phase converger's solution."""
    if not _backend_available():
        pytest.skip("native OQP backend not built")
    ref = _run(tmp_path, "twophase", system=_H2O, basis="6-31g", cas=_CAS44)
    trah = _run(tmp_path, f"trah_{hessian or 'fd'}", system=_H2O, basis="6-31g",
                cas=_CAS44, converger="trah", hessian=hessian)
    assert trah == pytest.approx(ref, abs=_TOL)


def test_trah_costs_far_fewer_ci_solves_than_ah(tmp_path):
    """The whole point, as a number: ``ah`` assembles the FD Hessian
    (2*n_par CI solves per macroiteration), ``trah`` never assembles anything."""
    if not _backend_available():
        pytest.skip("native OQP backend not built")
    _run(tmp_path, "cost_ah", system=_H2O, basis="6-31g", cas=_CAS44,
         converger="ah")
    _run(tmp_path, "cost_trah", system=_H2O, basis="6-31g", cas=_CAS44,
         converger="trah")

    def evals(name):
        text = (tmp_path / f"{name}.log").read_text()
        m = re.search(r"converger CI evaluations:\s*(\d+)", text)
        assert m, text[-2000:]
        return int(m.group(1))

    n_ah, n_trah = evals("cost_ah"), evals("cost_trah")
    assert n_trah * 3 < n_ah, (n_trah, n_ah)


def test_trah_is_its_own_converger_not_an_ah_alias():
    """``trah`` used to resolve to the ``ah`` code.  It now has its own."""
    from oqp.library.casscf import _native_converger_codes
    from oqp.library.casscf_convergers import _CONV_ALIASES

    assert _native_converger_codes({"converger": "trah"}) == (4, 0, True)
    assert _native_converger_codes({"converger": "ah"}) == (1, 0, True)
    assert _CONV_ALIASES["trah"] == "trah"
    assert _CONV_ALIASES["augmented-hessian"] == "ah"
