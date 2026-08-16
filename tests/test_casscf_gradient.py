"""Analytic state-specific CASSCF nuclear gradient (``source/modules/casscf_gradient.F90``).

The gradient is fully variational -- no orbital-response or CI-response term --
so the pins here are (a) that it IS the derivative of the CASSCF energy, and
(b) that the conditions making it valid are enforced rather than assumed.

Finite-difference convention
----------------------------
The reference gradient below is a FIVE-POINT, O(h^4) stencil, not the usual
two-point central difference.  That is not gold-plating: along a bond-stretching
coordinate the third energy derivative is large enough that the O(h^2)
truncation error of a two-point difference reaches ~1e-5 at h = 2e-3 Bohr, which
is the same size as a real defect and therefore useless as a test.  With the
O(h^4) formula the residual is ~1e-8, which is a statement about the gradient.

Each finite-difference test also checks that the energies it differentiated are
SMOOTH (the residual of a cubic fit through the five sampled points).  CASSCF
has multiple solutions, and a displaced point can converge to a different
branch or -- with an orbital-energy-ordered active space -- to a different
active space entirely.  When that happens the finite difference is not a
derivative of anything and the comparison has to be rejected rather than
reported as a gradient failure.  The geometries used here were chosen because
they are smooth; the smoothness assertion is what keeps a future change to the
CASSCF optimizer from turning this into a confusing false alarm.
"""
import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401
    except Exception:
        return False
    return (bool(getattr(oqp, "BACKEND_AVAILABLE", False))
            and hasattr(oqp, "casscf_gradient"))


# Linear H4/STO-3G CAS(2,2) with the symmetry deliberately broken.  At the
# EQUALLY SPACED geometry the optimizer converges to a different (higher)
# solution than it does at displaced geometries, so the energy is not smooth
# there and finite differences are meaningless; a slightly distorted chain has
# a single smooth branch.
_H4 = ("\nH 0.0  0.041 0.0"
       "\nH 0.017 0.0   0.795"
       "\nH -0.023 0.062 1.523"
       "\nH 0.009 -0.031 2.310")

_H2O = "\nO 0 0 0.1173\nH 0 0.7572 -0.4692\nH 0 -0.7572 -0.4692"


def _config(system, cas, ci, casscf, runtype="grad"):
    return {
        "input": {
            "system": system,
            "charge": "0",
            "basis": "sto-3g",
            "method": "casscf",
            "runtype": runtype,
        },
        "guess": {"type": "hcore"},
        "scf": {
            "type": "rhf",
            "multiplicity": "1",
            "conv": "1.0e-10",
            "maxit": "80",
            "forced_attempt": "3",
            "save_molden": "False",
        },
        "properties": {"grad": "0"},
        "cas": cas,
        "ci": ci,
        "casscf": casscf,
        "tests": {"exception": "True"},
    }


def _runner(tmp_path, project, config):
    from oqp.pyoqp import Runner
    return Runner(
        project=project,
        input_file=None,
        log=str(tmp_path / f"{project}.log"),
        input_dict=config,
        silent=1,
        usempi=False,
    )


_CAS22 = {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1",
          "max_det": "5000", "orbital_source": "rhf", "sort_orbitals": "energy"}
_CAS44 = {"active_electrons": "4", "active_orbitals": "4", "frozen_core": "3",
          "max_det": "5000", "orbital_source": "rhf", "sort_orbitals": "energy"}


def _ci(nroot):
    return {"nroot": str(nroot), "solver": "dense", "eig_tol": "1.0e-10",
            "integral_backend": "native", "target_spin": "any"}


def _casscf_opts(root, macro=100):
    return {"max_macro_iterations": str(macro), "optimizer": "newton",
            "root": str(root), "gradient_norm_tol": "1.0e-7",
            "energy_decrease_tol": "1.0e-10", "step_norm_tol": "1.0e-8",
            "max_rotation_norm": "2.0e-1", "level_shift": "1.0e-3",
            "canonicalize": "true"}


def _analytic_gradient(tmp_path, project, config):
    """Run the gradient and return ``(gradient (natom, 3), energy)``."""
    from oqp.library.single_point import SinglePoint, Gradient
    runner = _runner(tmp_path, project, config)
    mol = runner.mol
    sp = SinglePoint(mol)
    sp.energy()
    root = int(config["casscf"]["root"])
    energy = float(np.asarray(
        mol.data["OQP::CASSCF_ENERGIES"], dtype=float).reshape(-1)[root])
    grad = np.asarray(Gradient(mol).gradient(), dtype=float).reshape(-1, 3)
    return grad, energy, mol, sp


_RESTART_TAGS = ('OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::DM_A', 'OQP::DM_B',
                 'OQP::E_MO_A', 'OQP::E_MO_B')


def _fd_along(mol, sp, root, k, x0, e0, step):
    """O(h^4) derivative along coordinate ``k`` and the energy smoothness there."""
    snap = {}
    for tag in _RESTART_TAGS:
        try:
            snap[tag] = np.array(mol.data[tag], dtype=float, copy=True)
        except Exception:
            pass

    energies = {}
    for m in (-2, -1, 1, 2):
        x = x0.copy()
        x[k] += m * step
        for tag, val in snap.items():
            try:
                mol.data[tag][...] = val
            except Exception:
                pass
        mol.update_system(x)
        sp.energy()
        energies[m] = float(np.asarray(
            mol.data["OQP::CASSCF_ENERGIES"], dtype=float).reshape(-1)[root])

    for tag, val in snap.items():
        try:
            mol.data[tag][...] = val
        except Exception:
            pass
    mol.update_system(x0)

    deriv = (energies[-2] - 8 * energies[-1]
             + 8 * energies[1] - energies[2]) / (12 * step)
    xs = np.array([-2, -1, 0, 1, 2], dtype=float) * step
    ys = np.array([energies[-2], energies[-1], e0, energies[1], energies[2]],
                  dtype=float)
    smoothness = float(np.max(np.abs(ys - np.polyval(np.polyfit(xs, ys, 3), xs))))
    return deriv, smoothness


def test_casscf_gradient_matches_finite_difference_h4(tmp_path):
    """The analytic gradient is the derivative of the CASSCF energy."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    config = _config(_H4, _CAS22, _ci(1), _casscf_opts(0))
    grad, e0, mol, sp = _analytic_gradient(tmp_path, "casgrad_h4", config)
    x0 = np.asarray(mol.get_system(), dtype=float).reshape(-1).copy()

    # Differentiate along the two coordinates carrying the largest force, one
    # along the chain (large third derivative) and one perpendicular to it.
    flat = grad.reshape(-1)
    for k in (int(np.argmax(np.abs(flat))), 0):
        deriv, smoothness = _fd_along(mol, sp, 0, k, x0, e0, 1.0e-3)
        assert smoothness < 1.0e-8, (
            f"coordinate {k}: the displaced CASSCF energies are not smooth "
            f"(cubic-fit residual {smoothness:.2e}), so the finite difference "
            "is not a derivative -- the displaced points converged to a "
            "different CASSCF solution or active space")
        assert flat[k] == pytest.approx(deriv, abs=5.0e-7)


def test_casscf_gradient_is_translationally_invariant(tmp_path):
    """Sum of the forces vanishes: a rigid shift cannot change the energy.

    This is the sharpest cheap check on the Pulay and two-particle terms --
    a missing or mis-signed derivative-integral contribution breaks it, while
    the Hellmann-Feynman part alone would not.
    """
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    config = _config(_H2O, _CAS44, _ci(1), _casscf_opts(0, macro=50))
    grad, _e0, _mol, _sp = _analytic_gradient(tmp_path, "casgrad_trans", config)
    assert np.max(np.abs(grad.sum(axis=0))) < 1.0e-9


def test_casscf_gradient_follows_the_selected_root(tmp_path):
    """``[casscf] root`` selects which root is differentiated.

    An excited root is still a stationary point of the CASSCF energy -- the CI
    vector is an eigenvector whether or not it is the lowest one -- so the same
    state-specific expression applies with no Z-vector solve.  The pin is that
    the two roots give genuinely different energies AND different gradients.
    """
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    g0, e0, _m0, _s0 = _analytic_gradient(
        tmp_path, "casgrad_root0",
        _config(_H4, _CAS22, _ci(2), _casscf_opts(0)))
    g1, e1, mol1, sp1 = _analytic_gradient(
        tmp_path, "casgrad_root1",
        _config(_H4, _CAS22, _ci(2), _casscf_opts(1)))

    assert e1 > e0 + 0.05, "root 1 must be a genuinely higher state"
    assert np.max(np.abs(g1 - g0)) > 1.0e-3, (
        "the root-1 gradient is indistinguishable from the root-0 one; "
        "root selection is not reaching the gradient")

    # ...and the root-1 gradient is the derivative of the root-1 energy.
    x0 = np.asarray(mol1.get_system(), dtype=float).reshape(-1).copy()
    k = int(np.argmax(np.abs(g1.reshape(-1))))
    deriv, smoothness = _fd_along(mol1, sp1, 1, k, x0, e1, 1.0e-3)
    assert smoothness < 1.0e-8
    assert g1.reshape(-1)[k] == pytest.approx(deriv, abs=5.0e-7)


def test_casscf_gradient_supports_a_full_active_space(tmp_path):
    """Zero non-redundant rotations is stationary, not a gradient refusal."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    h2 = "\nH 0 0 0\nH 0 0 0.74"
    cas = {"active_electrons": "2", "active_orbitals": "2",
           "frozen_core": "0", "max_det": "5000",
           "orbital_source": "rhf", "sort_orbitals": "energy"}
    config = _config(h2, cas, _ci(1), _casscf_opts(0))
    grad, e0, mol, sp = _analytic_gradient(
        tmp_path, "casgrad_full_active", config)
    x0 = np.asarray(mol.get_system(), dtype=float).reshape(-1).copy()
    k = int(np.argmax(np.abs(grad.reshape(-1))))
    deriv, smoothness = _fd_along(mol, sp, 0, k, x0, e0, 1.0e-3)
    assert smoothness < 1.0e-8
    assert grad.reshape(-1)[k] == pytest.approx(deriv, abs=5.0e-7)


def test_casscf_gradient_accepts_excited_root_davidson_subspace(tmp_path):
    """Gradient eligibility does not inherit native-energy fallback rules."""
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    ci = _ci(2)
    ci.update({"solver": "davidson", "davidson_subspace": "4",
               "davidson_maxiter": "100"})
    config = _config(_H4, _CAS22, ci, _casscf_opts(1))
    grad, _energy, _mol, _sp = _analytic_gradient(
        tmp_path, "casgrad_root1_davidson", config)
    assert grad.shape == (4, 3)
    assert np.all(np.isfinite(grad))


def _preflight(method, runtype, root=0, state_average=None):
    """Typed config for the input checker (it reads native types, not strings)."""
    from oqp.utils import input_checker
    config = {
        "input": {"method": method, "runtype": runtype, "basis": "sto-3g",
                  "system": "\nH 0 0 0\nH 0 0 0.74\nH 0 0 1.48\nH 0 0 2.22"},
        "scf": {"type": "rhf", "multiplicity": 1},
        "tdhf": {"type": "rpa", "nstate": 1, "multiplicity": 1},
        "properties": {"grad": [0]},
        "cas": {"active_electrons": 2, "active_orbitals": 2, "frozen_core": 1,
                "max_det": 5000},
        "ci": {"nroot": 2, "solver": "dense"},
        "casscf": {"root": root},
        "optimize": {"istate": 0},
    }
    if state_average is not None:
        config["state_average"] = state_average
    return input_checker.check_input_values(config, raise_error=False, emit=False)


def test_state_specific_casscf_gradient_runtypes_pass_preflight():
    """runtype=grad and the gradient-driven optimizers are accepted for casscf."""
    for runtype in ("energy", "grad", "optimize", "ts", "mep", "irc"):
        report = _preflight("casscf", runtype)
        assert report.ok, f"runtype={runtype} rejected: {report.to_text()}"


def test_state_averaged_casscf_gradient_uses_the_numerical_path():
    """SA-CASSCF passes preflight for the numerical derivative on main."""
    sa = {"enabled": True, "nstate": 2}
    report = _preflight("sa-casscf", "grad", state_average=sa)
    assert report.ok, report.to_text()

    # The spelling method=casscf plus an enabled state average must take the
    # same numerical path, never the state-specific analytic expression.
    assert _preflight("casscf", "grad", state_average=sa).ok
    assert _preflight("sa-casscf", "energy", state_average=sa).ok


def test_state_averaged_casscf_is_refused_by_the_analytic_library(tmp_path):
    """The analytic entry point refuses an SA wavefunction when called directly.

    The public ``Gradient`` dispatcher sends this calculation to the numerical
    SA-CASSCF derivative; the analytic kernel still fails closed for a caller
    that bypasses that dispatch.
    """
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    from oqp.library.casscf_gradient import casscf_analytic_gradient
    from oqp.library.single_point import SinglePoint

    # runtype=energy so preflight accepts it; the gradient is then requested
    # programmatically, the way a driver would.
    config = _config(_H4, _CAS22, _ci(2), _casscf_opts(0), runtype="energy")
    config["input"]["method"] = "sa-casscf"
    config["state_average"] = {"enabled": "true", "nstate": "2"}
    runner = _runner(tmp_path, "casgrad_sa", config)
    SinglePoint(runner.mol).energy()

    with pytest.raises(ValueError, match="state-specific"):
        casscf_analytic_gradient(runner.mol)


def test_state_specific_gradient_selectors_are_public_slot_zero():
    """The physical [casscf] root never becomes a gradient-array index."""
    direct = _preflight("casscf", "grad", root=1)
    assert direct.ok, direct.to_text()

    # Rebuild the small configuration explicitly so only the selector changes.
    from oqp.utils import input_checker
    config = {
        "input": {"method": "casscf", "runtype": "grad", "basis": "sto-3g",
                  "system": _H4},
        "scf": {"type": "rhf", "multiplicity": 1},
        "tdhf": {"type": "rpa", "nstate": 1, "multiplicity": 1},
        "properties": {"grad": [1]},
        "optimize": {"istate": 0},
        "cas": {"active_electrons": 2, "active_orbitals": 2,
                "frozen_core": 1, "max_det": 5000},
        "ci": {"nroot": 2, "solver": "dense"},
        "casscf": {"root": 1},
    }
    bad_direct = input_checker.check_input_values(
        config, raise_error=False, emit=False)
    assert not bad_direct.ok
    assert "public gradient slot 0" in bad_direct.to_text()

    config["input"]["runtype"] = "optimize"
    config["properties"]["grad"] = [0]
    config["optimize"]["istate"] = 1
    bad_opt = input_checker.check_input_values(
        config, raise_error=False, emit=False)
    assert not bad_opt.ok
    assert "public gradient slot 0" in bad_opt.to_text()


def test_excited_root_optimizer_uses_matching_energy_and_gradient_slot(tmp_path):
    """The selected physical root occupies public energy and gradient slot 0."""
    if not _backend_available():
        pytest.skip("OpenQP backend not built; build liboqp to run this test")

    config = _config(_H4, _CAS22, _ci(2), _casscf_opts(1))
    runner = _runner(tmp_path, "casgrad_root1_opt", config)
    # Runner fills the optimizer defaults.  Select the derivative run type
    # after parsing so this focused one-step test does not need to duplicate
    # the complete [optimize] input section.
    runner.mol.config["input"]["runtype"] = "optimize"
    runner.mol.config["optimize"]["istate"] = 0

    from oqp.library.libscipy import StateSpecificOpt
    optimizer = StateSpecificOpt(runner.mol)
    energy, gradient = optimizer.one_step(
        np.asarray(runner.mol.get_system(), dtype=float).reshape(-1).copy())

    physical = np.asarray(
        runner.mol.data["OQP::CASSCF_ENERGIES"], dtype=float).reshape(-1)
    assert len(physical) >= 2
    assert len(runner.mol.energies) == 1
    assert energy == pytest.approx(physical[1], abs=1.0e-10)
    assert runner.mol.energies[0] == pytest.approx(physical[1], abs=1.0e-10)
    assert np.asarray(gradient).reshape(-1) == pytest.approx(
        np.asarray(runner.mol.grads[0]).reshape(-1), abs=1.0e-12)


def test_casscf_gradient_rejects_generalized_fock_asymmetry():
    """CI non-stationarity is rejected even when |g_orb| is zero."""
    from oqp.library.casscf_gradient import (
        _FOCK_ASYMMETRY_LIMIT, _enforce_stationarity,
    )

    _enforce_stationarity(0.0, 0.5 * _FOCK_ASYMMETRY_LIMIT, 1.0e-7)
    for bad in (2.0 * _FOCK_ASYMMETRY_LIMIT, np.inf, np.nan):
        with pytest.raises(ValueError, match="generalized Fock asymmetry"):
            _enforce_stationarity(0.0, bad, 1.0e-7)


def test_casscf_gradient_reports_stationarity_diagnostics(tmp_path):
    """The driver reports |g_orb| and the generalized-Fock asymmetry.

    The state-specific expression has no orbital-response term, so it is the
    energy derivative only at a stationary point.  Both diagnostics must come
    back small for a converged run -- they are what lets the Python side refuse
    a gradient taken somewhere the formula does not apply.
    """
    if not _backend_available():
        pytest.skip("native OQP backend not built; build liboqp to run this test")

    import re
    config = _config(_H2O, _CAS44, _ci(1), _casscf_opts(0, macro=50))
    log = str(tmp_path / "casgrad_diag.log")
    from oqp.pyoqp import Runner
    from oqp.library.single_point import SinglePoint, Gradient
    runner = Runner(project="casgrad_diag", input_file=None, log=log,
                    input_dict=config, silent=1, usempi=False)
    SinglePoint(runner.mol).energy()
    Gradient(runner.mol).gradient()

    text = open(log, encoding="utf-8", errors="replace").read()
    gnorm = re.search(r"Orbital-rotation gradient norm\s+(\S+)", text)
    fasym = re.search(r"Generalized Fock asymmetry\s+(\S+)", text)
    assert gnorm and fasym, "the gradient driver did not report its diagnostics"
    assert float(gnorm.group(1).replace("D", "E")) < 1.0e-6
    assert float(fasym.group(1).replace("D", "E")) < 1.0e-6
