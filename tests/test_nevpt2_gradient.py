"""Tests for the analytic strongly contracted NEVPT2 (SC-NEVPT2) nuclear gradient.

The gradient is analytic; finite differences appear here only as the reference
it is measured against, and only in the end-to-end tests.  The unit tests below
check the derivative algebra directly, so a regression is localized rather than
showing up as one wrong number at the end.

Layout
------
* derivative algebra -- the exact adjoints of the SC-NEVPT2 energy functional;
* structure -- the two invariance facts the Lagrangian is built on, including
  the negative control that the semicanonical response is not optional;
* end-to-end -- the nuclear gradient against multi-step 5-point central
  differences, translational/rotational behaviour, and route agreement;
* scope -- every configuration the derivative is NOT the derivative of;
* native contraction -- the Fortran derivative-integral entry point pinned
  against the independently validated CASSCF gradient.
"""
import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False)) and hasattr(
        oqp, "fci_ao_integrals")


def _native_gradient_available() -> bool:
    if not _backend_available():
        return False
    try:
        from oqp.library.fci import _lib_backend
    except Exception:
        return False
    backend = _lib_backend()
    return backend is not None and hasattr(backend[0], "nevpt2_gradient")


needs_backend = pytest.mark.skipif(
    not _backend_available(),
    reason="native OQP backend not built; build liboqp to run SC-NEVPT2 tests")
needs_native_gradient = pytest.mark.skipif(
    not _native_gradient_available(),
    reason="liboqp has no nevpt2_gradient entry point; rebuild liboqp")


_H4 = """
   H   0.000000000   0.000000000   0.000000000
   H   0.000000000   0.000000000   0.735000000
   H   0.000000000   0.000000000   1.560000000
   H   0.000000000   0.000000000   2.310000000"""

_H2O = """
   O   0.000000000   0.000000000   0.000000000
   H   0.000000000  -0.757000000   0.587000000
   H   0.000000000   0.757000000   0.587000000"""

#: Tight CASSCF convergence.  The Lagrangian assumes a stationary reference and
#: its error is FIRST order in the residual orbital-rotation gradient, so the
#: agreement these tests assert is only meaningful with |g_orb| driven down.
_TIGHT_CASSCF = {"max_macro_iterations": "300",
                 "gradient_norm_tol": "1.0e-10",
                 "energy_decrease_tol": "1.0e-13"}

_SC_NEVPT2 = {"reference": "casscf", "variant": "caspt2", "h0": "dyall",
              "contraction": "strong", "frozen": "0",
              "ipea_shift": "0.0", "level_shift": "0.0",
              "imaginary_shift": "0.0", "semi_canonical": "true"}


def _runner(tmp_path, project, *, system, basis, cas, pt2, runtype="energy",
            casscf=None, coords=None, extra=None):
    from oqp.pyoqp import Runner
    config = {
        "input": {"system": system, "charge": "0", "basis": basis,
                  "method": "caspt2", "runtype": runtype},
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": "1", "maxit": "100",
                "conv": "1.0e-11", "forced_attempt": "3",
                "save_molden": "False"},
        "properties": {"scf_prop": "", "grad": "0"},
        "casscf": dict(casscf or _TIGHT_CASSCF),
        "cas": cas,
        "ci": {"nroot": "1", "solver": "dense", "eig_tol": "1.0e-12",
               "integral_backend": "native", "integral_cutoff": "5.0e-12",
               "target_spin": "any"},
        "pt2": pt2,
        "tests": {"exception": "True"},
    }
    for section, block in (extra or {}).items():
        config.setdefault(section, {}).update(block)
    runner = Runner(project=project, input_file=None,
                    log=str(tmp_path / f"{project}.log"),
                    input_dict=config, silent=1, usempi=False)
    if coords is not None:
        runner.mol.update_system(np.asarray(coords, dtype=float))
    return runner


_H4_CAS22 = {"active_electrons": "2", "active_orbitals": "2",
             "frozen_core": "1", "max_det": "20000", "max_memory": "2048",
             "orbital_source": "rhf", "sort_orbitals": "energy"}
_H4_CAS44 = {"active_electrons": "4", "active_orbitals": "4",
             "frozen_core": "0", "max_det": "20000", "max_memory": "4096",
             "orbital_source": "rhf", "sort_orbitals": "energy"}


# --------------------------------------------------------------------------
# 1. derivative algebra: the exact adjoints of the SC-NEVPT2 energy functional
# --------------------------------------------------------------------------
def _random_case(seed=4242, nbf=8, ncore=2, nact=4):
    from oqp.library.fci import _determinants
    rng = np.random.default_rng(seed)
    h = rng.normal(size=(nbf, nbf))
    h = 0.5 * (h + h.T) - 3.0 * np.eye(nbf)
    g = rng.normal(size=(nbf,) * 4) * 0.1
    g = g + g.transpose(1, 0, 2, 3)
    g = g + g.transpose(0, 1, 3, 2)
    g = g + g.transpose(2, 3, 0, 1)
    eps = np.sort(rng.normal(size=nbf)) * 2.0
    nelec = (2, 2)
    dets = _determinants(nact, nelec)
    ci = rng.normal(size=len(dets))
    ci /= np.linalg.norm(ci)
    return h, g, eps, ncore, nact, nelec, ci, rng


def _five_point(f, t):
    return (f(-2 * t) - 8 * f(-t) + 8 * f(t) - f(2 * t)) / (12 * t)


@needs_backend
def test_ein_bar_is_the_exact_einsum_adjoint():
    """The generic adjoint helper on the two structurally awkward cases.

    A repeated index inside one operand is a diagonal READ, whose adjoint is a
    diagonal WRITE; an index private to the differentiated operand was summed
    away, so its adjoint broadcasts.  Both are checked against the definition
    ``dE/dA_ij = lim (E(A + t e_ij) - E(A)) / t`` evaluated exactly, which for a
    multilinear expression is just the expression with that element set to one.
    """
    from oqp.library.nevpt2_adjoint import _ein_bar

    rng = np.random.default_rng(7)
    n = 4
    # repeated index: out_mn = sum_j A[m,j,j,n]
    A = rng.normal(size=(n, n, n, n))
    ob = rng.normal(size=(n, n))
    bar = _ein_bar("mjjn->mn", ob, (A,), 0, A.shape)
    brute = np.zeros_like(A)
    for m in range(n):
        for j in range(n):
            for k in range(n):
                for q in range(n):
                    if j == k:
                        brute[m, j, k, q] = ob[m, q]
    assert np.allclose(bar, brute, atol=1e-14)

    # private summed index: out_i = sum_jk A[i,j] B[j,k]
    A2 = rng.normal(size=(n, n))
    B2 = rng.normal(size=(n, n))
    ob2 = rng.normal(size=n)
    bar_a = _ein_bar("ij,jk->i", ob2, (A2, B2), 0, A2.shape)
    brute_a = np.einsum("i,jk->ij", ob2, B2)
    assert np.allclose(bar_a, brute_a, atol=1e-13)


@needs_backend
@pytest.mark.parametrize("target", ["h1e", "eri", "eps"])
def test_energy_adjoints_match_directional_derivatives(target):
    """``dE2/dh``, ``dE2/dg`` and ``dE2/deps`` against 5-point directional FD."""
    from oqp.library.nevpt2_sc import sc_nevpt2_energy
    from oqp.library.nevpt2_adjoint import sc_nevpt2_energy_adjoints

    h, g, eps, ncore, nact, nelec, ci, rng = _random_case()
    e2, _comp, hbar, gbar, epsbar, _dmbars = sc_nevpt2_energy_adjoints(
        h, g, eps, ncore, nact, nelec, ci)
    assert e2 == pytest.approx(
        sc_nevpt2_energy(h, g, eps, ncore, nact, nelec, ci,
                         max_memory=8192)[0], abs=0.0)

    def energy(hh, gg, ee):
        return sc_nevpt2_energy(hh, gg, ee, ncore, nact, nelec, ci,
                                max_memory=8192)[0]

    if target == "h1e":
        d = rng.normal(size=h.shape)
        d = 0.5 * (d + d.T)
        analytic = float(np.sum(hbar * d))
        numeric = _five_point(lambda t: energy(h + t * d, g, eps), 1.0e-4)
    elif target == "eri":
        d = rng.normal(size=g.shape) * 0.05
        d = d + d.transpose(1, 0, 2, 3)
        d = d + d.transpose(0, 1, 3, 2)
        d = d + d.transpose(2, 3, 0, 1)
        analytic = float(np.sum(gbar * d))
        numeric = _five_point(lambda t: energy(h, g + t * d, eps), 1.0e-4)
    else:
        d = rng.normal(size=eps.shape)
        analytic = float(np.sum(epsbar * d))
        numeric = _five_point(lambda t: energy(h, g, eps + t * d), 1.0e-4)

    assert analytic == pytest.approx(numeric, rel=1.0e-7, abs=1.0e-7)


@needs_backend
def test_ci_derivative_matches_finite_differences():
    """``dE2/dc`` from the density-matrix adjoints, against a direct FD in c.

    The SC-NEVPT2 energy depends on the CI vector only through the active
    1-4 RDMs, so this closes the chain from the RDM adjoints to the CI
    parameters that the response equations are solved in.
    """
    from oqp.library.nevpt2_sc import sc_nevpt2_energy
    from oqp.library.nevpt2_adjoint import sc_nevpt2_energy_adjoints
    from oqp.library.nevpt2_gradient import _ci_bar_from_dm_bars
    from oqp.library.casscf_hessian import _excitation_matrices

    h, g, eps, ncore, nact, nelec, ci, rng = _random_case()
    _e2, _comp, _hb, _gb, _eb, dmbars = sc_nevpt2_energy_adjoints(
        h, g, eps, ncore, nact, nelec, ci)
    _dets, stack = _excitation_matrices(nact, nelec[0], nelec[1])
    ci_bar = _ci_bar_from_dm_bars(dmbars, stack, ci)

    d = rng.normal(size=ci.shape)
    d -= float(ci @ d) * ci                     # stay on the normalized sphere
    analytic = float(np.sum(ci_bar * d))
    numeric = _five_point(
        lambda t: sc_nevpt2_energy(h, g, eps, ncore, nact, nelec,
                                   ci + t * d, max_memory=8192)[0], 1.0e-4)
    assert analytic == pytest.approx(numeric, rel=1.0e-7, abs=1.0e-8)


# --------------------------------------------------------------------------
# 2. structure: the invariance facts the Lagrangian is built on
# --------------------------------------------------------------------------
def _rotate_integrals(h, g, U):
    return (U.T @ h @ U,
            np.einsum("pqrs,pa,qb,rc,sd->abcd", g, U, U, U, U, optimize=True))


def _effective_fock(h, g, D):
    return (h + np.einsum("rs,pqrs->pq", D, g, optimize=True)
            - 0.5 * np.einsum("rs,prsq->pq", D, g, optimize=True))


def _full_density(ci, dets, ncore, nact, nbf):
    from oqp.library.rdm import make_rdm1_spatial
    D = np.zeros((nbf, nbf))
    for i in range(ncore):
        D[i, i] = 2.0
    D[np.ix_(range(ncore, ncore + nact), range(ncore, ncore + nact))] = \
        make_rdm1_spatial(ci, dets, nact)
    return D


@needs_backend
def test_active_rotation_leaves_sc_nevpt2_invariant():
    """The active semicanonical gauge is FREE, so it needs no multiplier.

    Rotating the active orbitals and the CI vector together is a change of
    parameterization, not of the wavefunction.  This is what makes the
    intra-active orbital stationarity of the Lagrangian follow from CI
    stationarity instead of needing its own Lagrange multiplier.

    The orbital energies must be the ones the driver actually uses -- the
    generalized-Fock diagonal -- on BOTH sides, or the two evaluations are not
    the same functional and the comparison says nothing.
    """
    from oqp.library.nevpt2_sc import sc_nevpt2_energy
    from oqp.library.caspt2_dyall import _rotate_det_vector
    from oqp.library.fci import _determinants

    h, g, _eps_unused, ncore, nact, nelec, ci, _rng = _random_case()
    nbf = h.shape[0]
    dets = _determinants(nact, nelec)
    D = _full_density(ci, dets, ncore, nact, nbf)
    eps = np.diag(_effective_fock(h, g, D)).copy()
    e0 = sc_nevpt2_energy(h, g, eps, ncore, nact, nelec, ci, max_memory=8192)[0]

    for angle in (1.0e-3, 0.137):
        R = np.eye(nact)
        c, s = np.cos(angle), np.sin(angle)
        R[0, 0] = R[1, 1] = c
        R[0, 1], R[1, 0] = s, -s
        U = np.eye(nbf)
        U[np.ix_(range(ncore, ncore + nact), range(ncore, ncore + nact))] = R
        h_r, g_r = _rotate_integrals(h, g, U)
        ci_r = _rotate_det_vector(ci, R, nact, nelec)
        D_r = _full_density(ci_r, dets, ncore, nact, nbf)
        eps_r = np.diag(_effective_fock(h_r, g_r, D_r)).copy()
        e_r = sc_nevpt2_energy(h_r, g_r, eps_r, ncore, nact, nelec, ci_r,
                               max_memory=8192)[0]
        assert e_r == pytest.approx(e0, abs=1.0e-10)


@needs_backend
def test_semicanonical_response_is_not_optional():
    """Negative control: the energy DOES change under a virtual-virtual rotation.

    The strongly contracted denominators are diagonal, which is correct only in
    the semicanonical basis, so SC-NEVPT2 is not invariant here even with the
    orbital energies recomputed.  If this test ever passed as an invariance, the
    ``Lambda`` term in the gradient would be unnecessary -- and it is not: on
    this probe it moves the derivative by several percent.
    """
    from oqp.library.nevpt2_sc import sc_nevpt2_energy
    from oqp.library.fci import _determinants

    h, g, _eps_unused, ncore, nact, nelec, ci, _rng = _random_case()
    nbf = h.shape[0]
    dets = _determinants(nact, nelec)
    D = _full_density(ci, dets, ncore, nact, nbf)
    eps = np.diag(_effective_fock(h, g, D)).copy()
    e0 = sc_nevpt2_energy(h, g, eps, ncore, nact, nelec, ci, max_memory=8192)[0]

    lo = ncore + nact
    U = np.eye(nbf)
    angle = 1.0e-3
    c, s = np.cos(angle), np.sin(angle)
    U[lo, lo] = U[lo + 1, lo + 1] = c
    U[lo, lo + 1], U[lo + 1, lo] = s, -s
    h_r, g_r = _rotate_integrals(h, g, U)
    eps_r = np.diag(_effective_fock(h_r, g_r, D)).copy()
    e_r = sc_nevpt2_energy(h_r, g_r, eps_r, ncore, nact, nelec, ci,
                           max_memory=8192)[0]
    assert abs(e_r - e0) > 1.0e-8


# --------------------------------------------------------------------------
# 3. end-to-end: the nuclear gradient
# --------------------------------------------------------------------------
def _analytic_gradient(tmp_path, project, *, system, basis, cas, coords=None,
                       pt2=None):
    from oqp.library.single_point import SinglePoint
    from oqp.library.nevpt2_gradient import sc_nevpt2_analytic_gradient

    runner = _runner(tmp_path, project, system=system, basis=basis, cas=cas,
                     pt2=dict(pt2 or _SC_NEVPT2), coords=coords)
    SinglePoint(runner.mol).energy()
    grad = np.asarray(sc_nevpt2_analytic_gradient(runner.mol), dtype=float)
    return runner, grad.reshape(-1)


def _energy_at(tmp_path, project, *, system, basis, cas, coords, pt2=None):
    from oqp.library.single_point import SinglePoint
    runner = _runner(tmp_path, project, system=system, basis=basis, cas=cas,
                     pt2=dict(pt2 or _SC_NEVPT2), coords=coords)
    return float(SinglePoint(runner.mol).energy()[0])


@needs_native_gradient
def test_analytic_gradient_matches_five_point_finite_differences(tmp_path):
    """H4/STO-3G CAS(2,2): analytic derivative vs 5-point central differences.

    Only the four z components are differentiated: linear H4 lies on the z axis,
    so the eight transverse components are zero by symmetry and are checked
    separately (and much more cheaply) below.
    """
    runner, ana = _analytic_gradient(tmp_path, "scnevpt2_fd", system=_H4,
                                     basis="sto-3g", cas=_H4_CAS22)
    x0 = np.asarray(runner.mol.get_system(), dtype=float).reshape(-1).copy()

    for step in (2.0e-3, 1.0e-3):
        for i in (2, 5, 8, 11):
            vals = {}
            for k in (-2, -1, 1, 2):
                x = x0.copy()
                x[i] += k * step
                vals[k] = _energy_at(tmp_path, "scnevpt2_fd_disp", system=_H4,
                                     basis="sto-3g", cas=_H4_CAS22, coords=x)
            num = (vals[-2] - 8 * vals[-1] + 8 * vals[1] - vals[2]) / (12 * step)
            assert ana[i] == pytest.approx(num, abs=5.0e-8), (
                f"coordinate {i} at step {step}")


@needs_native_gradient
def test_gradient_is_translationally_invariant(tmp_path):
    """Sum over atoms vanishes: no net force on a free molecule."""
    _runner_obj, ana = _analytic_gradient(tmp_path, "scnevpt2_trans",
                                          system=_H2O, basis="sto-3g",
                                          cas={"active_electrons": "4",
                                               "active_orbitals": "4",
                                               "frozen_core": "3",
                                               "max_det": "20000",
                                               "max_memory": "4096",
                                               "orbital_source": "rhf",
                                               "sort_orbitals": "energy"})
    total = ana.reshape(-1, 3).sum(axis=0)
    assert np.allclose(total, 0.0, atol=1.0e-9)


@needs_native_gradient
def test_gradient_rotates_with_the_molecule(tmp_path):
    """A rigidly rotated geometry returns the rigidly rotated gradient.

    This is the rotational counterpart of the translation test and, unlike a
    torque check on a linear molecule, it is non-trivial for every geometry.
    """
    cas = {"active_electrons": "4", "active_orbitals": "4",
           "frozen_core": "3", "max_det": "20000", "max_memory": "4096",
           "orbital_source": "rhf", "sort_orbitals": "energy"}
    runner, ana = _analytic_gradient(tmp_path, "scnevpt2_rot0", system=_H2O,
                                     basis="sto-3g", cas=cas)
    x0 = np.asarray(runner.mol.get_system(), dtype=float).reshape(-1, 3)

    theta = 0.37
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    _runner2, ana_rot = _analytic_gradient(
        tmp_path, "scnevpt2_rot1", system=_H2O, basis="sto-3g", cas=cas,
        coords=(x0 @ R.T).ravel())
    assert np.allclose(ana_rot.reshape(-1, 3), ana.reshape(-1, 3) @ R.T,
                       atol=1.0e-8)


@needs_native_gradient
def test_analytic_and_numerical_routes_agree_through_the_public_dispatch(tmp_path):
    """``[pt2] gradient=analytic`` and ``=numerical`` give the same gradient.

    This exercises the production ``runtype=grad`` dispatch rather than the
    library entry point, so it also pins that the analytic route is what the
    optimizer, MEP and IRC drivers actually receive.
    """
    from oqp.library.single_point import SinglePoint, Gradient

    results = {}
    for route in ("analytic", "numerical"):
        pt2 = dict(_SC_NEVPT2, gradient=route, grad_step="1.0e-3")
        runner = _runner(tmp_path, f"scnevpt2_route_{route}", system=_H4,
                         basis="sto-3g", cas=_H4_CAS22, pt2=pt2,
                         runtype="grad")
        SinglePoint(runner.mol).energy()
        results[route] = np.asarray(
            Gradient(runner.mol).gradient(), dtype=float).reshape(-1)
    assert np.allclose(results["analytic"], results["numerical"], atol=2.0e-6)


@needs_native_gradient
def test_run_reports_its_own_stationarity_and_symmetry_diagnostics(tmp_path):
    """The self-checks the derivation depends on are reported, and pass.

    ``max|X - X^T|`` is the strongest single check available: a stationary
    Lagrangian has a symmetric generalized Fock, so a nonzero value would mean
    a response multiplier is wrong rather than merely imprecise.
    """
    runner, _ana = _analytic_gradient(tmp_path, "scnevpt2_diag", system=_H4,
                                      basis="6-31g", cas=_H4_CAS44)
    log = open(runner.mol.log).read()
    assert "analytic SC-NEVPT2 nuclear gradient" in log

    def _value(label):
        for line in log.splitlines():
            if label in line:
                return float(line.split()[-1])
        raise AssertionError(f"diagnostic {label!r} not reported")

    assert _value("CASSCF orbital gradient") < 1.0e-8
    assert _value("Lagrangian asymmetry") < 1.0e-8
    assert _value("response solve relative residual") < 1.0e-10


@needs_native_gradient
def test_gradient_is_continuous_along_a_bond_scan(tmp_path):
    """Root continuity: the analytic gradient tracks the energy along a scan.

    Each geometry is an independent calculation -- its own SCF, CASSCF,
    semicanonicalization and response solve -- so a root flip, a reordered
    active space or a discontinuity in the semicanonical basis would show up
    here as a gradient that no longer differentiates the sampled curve.

    The scan deliberately runs in the +z direction.  Compressing this bond
    instead crosses a CASSCF solution change near -0.028 Bohr (measured: the
    reference energy's slope jumps from -0.0135 to -0.031 Eh/Bohr across it),
    where the surface being differentiated is genuinely discontinuous and no
    finite-difference comparison is meaningful.  That is a property of the
    CAS(2,2) reference for stretched H4, not of the derivative: the analytic
    gradient reproduces a local five-point difference to ~1e-10 on each smooth
    branch either side of the crossing.
    """
    from oqp.library.single_point import SinglePoint

    base = _runner(tmp_path, "scnevpt2_scan0", system=_H4, basis="sto-3g",
                   cas=_H4_CAS22, pt2=dict(_SC_NEVPT2))
    x0 = np.asarray(base.mol.get_system(), dtype=float).reshape(-1).copy()

    energies, grads, offsets = [], [], (0.0, 0.03, 0.06, 0.09, 0.12)
    for k, d in enumerate(offsets):
        x = x0.copy()
        x[11] += d                                   # stretch the terminal H-H
        runner, ana = _analytic_gradient(tmp_path, f"scnevpt2_scan{k}",
                                         system=_H4, basis="sto-3g",
                                         cas=_H4_CAS22, coords=x)
        energies.append(float(runner.mol.energies[0]))
        grads.append(ana[11])

    # The five samples are equally spaced, so they carry their own O(dx^4)
    # five-point derivative at the centre -- for free, and far tighter than a
    # secant.  The analytic gradient there must reproduce it.
    dx = offsets[1] - offsets[0]
    five_point = (energies[0] - 8 * energies[1]
                  + 8 * energies[3] - energies[4]) / (12 * dx)
    # The 3.4e-7 residual is the formula's own dx^4 truncation at dx = 0.03,
    # not a disagreement; the same comparison at dx = 1e-3 closes to ~1e-10.
    assert grads[2] == pytest.approx(five_point, abs=2.0e-6)
    # ... and the gradient must vary smoothly across the scan: a root flip or a
    # discontinuity in the semicanonical basis would break the second
    # difference, which is O(dx^2) for a smooth curve.
    second = np.diff(grads, 2)
    assert np.all(np.abs(second) < 5.0e-3)
    assert float(np.max(np.abs(np.diff(second)))) < 1.0e-3
    assert all(np.isfinite(grads))


# --------------------------------------------------------------------------
# 4. scope: configurations this derivative is NOT the derivative of
# --------------------------------------------------------------------------
@needs_backend
@pytest.mark.parametrize("pt2_override,needle", [
    ({"contraction": "none"}, "STRONGLY CONTRACTED"),
    ({"reference": "casci"}, "CASSCF reference"),
    ({"level_shift": "0.1"}, "level_shift"),
    ({"h0": "fock", "contraction": "none"}, "STRONGLY CONTRACTED"),
])
def test_out_of_scope_calculations_are_refused_not_approximated(
        tmp_path, pt2_override, needle):
    """``gradient=analytic`` reports why rather than returning a wrong number."""
    from oqp.library.nevpt2_gradient import sc_nevpt2_gradient_route

    pt2 = dict(_SC_NEVPT2, gradient="analytic", **pt2_override)
    runner = _runner(tmp_path, "scnevpt2_scope", system=_H4, basis="sto-3g",
                     cas=_H4_CAS22, pt2=pt2)
    with pytest.raises(ValueError, match=needle):
        sc_nevpt2_gradient_route(runner.mol)


@needs_backend
@pytest.mark.parametrize("pt2_override", [
    {"contraction": "none"},
    {"reference": "casci"},
])
def test_auto_falls_back_to_central_differences_out_of_scope(
        tmp_path, pt2_override):
    """``gradient=auto`` degrades to the numerical route with a stated reason."""
    from oqp.library.nevpt2_gradient import sc_nevpt2_gradient_route

    pt2 = dict(_SC_NEVPT2, gradient="auto", **pt2_override)
    runner = _runner(tmp_path, "scnevpt2_auto", system=_H4, basis="sto-3g",
                     cas=_H4_CAS22, pt2=pt2)
    route, reason = sc_nevpt2_gradient_route(runner.mol)
    assert route == "numerical"
    assert reason


@needs_backend
def test_numerical_is_honoured_even_when_the_analytic_route_applies(tmp_path):
    from oqp.library.nevpt2_gradient import sc_nevpt2_gradient_route

    pt2 = dict(_SC_NEVPT2, gradient="numerical")
    runner = _runner(tmp_path, "scnevpt2_forced_num", system=_H4,
                     basis="sto-3g", cas=_H4_CAS22, pt2=pt2)
    route, reason = sc_nevpt2_gradient_route(runner.mol)
    assert route == "numerical"
    assert "gradient=numerical" in reason


@needs_backend
def test_state_averaged_reference_is_refused(tmp_path):
    from oqp.library.nevpt2_gradient import sc_nevpt2_gradient_route

    pt2 = dict(_SC_NEVPT2, gradient="analytic")
    runner = _runner(tmp_path, "scnevpt2_sa", system=_H4, basis="sto-3g",
                     cas=_H4_CAS22, pt2=pt2)
    # [state_average] on a caspt2 route is rejected by the input checker at
    # construction, so the state-averaged objective is staged on the built
    # config -- which is also how the PT2 driver itself sets it up for a
    # multi-root reference.
    runner.mol.config["state_average"] = {"enabled": "true"}
    with pytest.raises(ValueError, match="STATE-SPECIFIC"):
        sc_nevpt2_gradient_route(runner.mol)


@needs_backend
def test_degenerate_semicanonical_orbitals_are_refused():
    """The semicanonical multiplier divides by an intra-block orbital-energy gap.

    Exactly degenerate inactive or virtual orbitals leave the semicanonical
    basis -- and therefore the response of the rotation that produces it --
    undetermined.  Refusing names the offending pair instead of returning a
    large, plausible-looking number.
    """
    from oqp.library.nevpt2_gradient import _semicanonical_multipliers

    nbf, ncore, nact = 8, 2, 4
    eps = np.array([-1.0, -1.0, 0.1, 0.2, 0.3, 0.4, 1.0, 2.0])   # core pair tied
    F0 = np.eye(nbf) * 0.5
    with pytest.raises(ValueError, match="degenerate"):
        _semicanonical_multipliers(F0, eps, ncore, nact, nbf)


# --------------------------------------------------------------------------
# 5. the native derivative-integral contraction
# --------------------------------------------------------------------------
@needs_native_gradient
def test_native_contraction_reproduces_the_casscf_gradient(tmp_path):
    """Pin ``nevpt2_gradient.F90`` against the independently validated CASSCF one.

    The Fortran entry point is a pure quadrature of supplied AO densities, so
    feeding it the plain CASSCF relaxed densities -- no second-order or response
    contribution at all -- must reproduce what ``casscf_gradient`` computes from
    the same wavefunction through its own, differently factorized route.  That
    isolates every AO-side convention (packing, the energy-weighted density's
    sign, the four-fold grd2 density factor, the spherical-to-Cartesian
    expansion) from the perturbation theory above it.
    """
    import oqp
    from oqp.library.single_point import SinglePoint
    from oqp.library.casscf import (
        _full_rdms, _generalized_fock, _casscf_options)
    from oqp.library.casscf_gradient import casscf_analytic_gradient
    from oqp.library.fci import (
        _unpack_lower_triangle, contiguous_active_space,
        settings_from_casci_config)
    from oqp.library.nevpt2_gradient import _gradient_backend, _sym8
    from oqp.library.rdm import make_rdm1_spatial, make_rdm2_spatial
    from oqp.pyoqp import Runner

    config = {
        "input": {"system": _H4, "charge": "0", "basis": "sto-3g",
                  "method": "casscf", "runtype": "grad"},
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": "1", "maxit": "100",
                "conv": "1.0e-11", "forced_attempt": "3",
                "save_molden": "False"},
        "properties": {"scf_prop": "", "grad": "0"},
        "casscf": dict(_TIGHT_CASSCF),
        "cas": _H4_CAS22,
        "ci": {"nroot": "1", "solver": "dense", "eig_tol": "1.0e-12",
               "integral_backend": "native", "target_spin": "any"},
        "tests": {"exception": "True"},
    }
    runner = Runner(project="casscf_pin", input_file=None,
                    log=str(tmp_path / "casscf_pin.log"),
                    input_dict=config, silent=1, usempi=False)
    mol = runner.mol
    SinglePoint(mol).energy()
    reference = np.asarray(casscf_analytic_gradient(mol), dtype=float).reshape(-1)

    # Rebuild the same densities and push them through the NEVPT2 entry point.
    settings = settings_from_casci_config(mol.config)
    nbf = int(mol.data.get_basis()["nbf"])
    nelec = (int(mol.data["nelec_A"]), int(mol.data["nelec_B"]))
    ncore, nact, active_nelec, _plan = contiguous_active_space(
        nbf, nelec, settings, "CASSCF")
    oqp.fci_ao_integrals(mol)
    hcore_ao = _unpack_lower_triangle(
        np.asarray(mol.data["OQP::Hcore"], dtype=float), nbf)
    eri_ao = np.asarray(mol.data["OQP::AO_ERI"], dtype=float).reshape(
        (nbf,) * 4, order="F")
    C = np.asarray(mol.data["OQP::VEC_MO_A"], dtype=float).reshape((nbf, nbf)).T
    h1e = C.T @ hcore_ao @ C
    eri = np.einsum("pqrs,pa,qb,rc,sd->abcd", eri_ao, C, C, C, C, optimize=True)

    from oqp.library.casscf import _solve_active
    root = int(_casscf_options(mol.config).root)
    _e, coeffs, dets, _D, _G = _solve_active(
        h1e, eri, ncore, nact, active_nelec, float(mol.mol_energy.nenergy),
        settings, np.array([1.0]), [root])
    ci = np.asarray(coeffs[:, root], dtype=float)
    ci /= np.linalg.norm(ci)
    D_cas, G_cas = _full_rdms(make_rdm1_spatial(ci, dets, nact),
                              make_rdm2_spatial(ci, dets, nact),
                              ncore, nact, nbf)
    G_cas = _sym8(G_cas)
    X = _generalized_fock(D_cas, G_cas, h1e, eri)
    W = 0.5 * (X + X.T)

    D_ao = C @ D_cas @ C.T
    W_ao = C @ W @ C.T
    G_ao = _sym8(np.einsum("pqrs,ap,bq,cr,ds->abcd", G_cas, C, C, C, C,
                           optimize=True))
    lib, ffi = _gradient_backend()
    d_c = np.ascontiguousarray(D_ao)
    w_c = np.ascontiguousarray(W_ao)
    g_c = np.ascontiguousarray(G_ao)
    status = int(lib.nevpt2_gradient(
        mol.data._data, int(nbf),
        ffi.cast("double *", d_c.ctypes.data),
        ffi.cast("double *", w_c.ctypes.data),
        ffi.cast("double *", g_c.ctypes.data)))
    assert status == 0
    natom = int(mol.data["natom"])
    produced = np.asarray(mol.get_grad(), dtype=float).reshape(natom, 3).ravel()
    assert np.allclose(produced, reference, atol=1.0e-9)


@needs_native_gradient
def test_native_contraction_rejects_an_unsymmetric_two_particle_density(tmp_path):
    """The eight-fold symmetry grd2 assumes is verified, not trusted.

    A density lacking it would be contracted against the wrong integrals, so the
    failure has to be loud rather than a quietly wrong gradient.
    """
    import oqp
    from oqp.library.single_point import SinglePoint
    from oqp.library.nevpt2_gradient import _gradient_backend
    from oqp.pyoqp import Runner

    config = {
        "input": {"system": _H4, "charge": "0", "basis": "sto-3g",
                  "method": "casscf", "runtype": "energy"},
        "guess": {"type": "hcore"},
        "scf": {"type": "rhf", "multiplicity": "1", "maxit": "100",
                "conv": "1.0e-10", "forced_attempt": "3",
                "save_molden": "False"},
        "properties": {"scf_prop": ""},
        "casscf": dict(_TIGHT_CASSCF),
        "cas": _H4_CAS22,
        "ci": {"nroot": "1", "solver": "dense", "eig_tol": "1.0e-11",
               "integral_backend": "native", "target_spin": "any"},
        "tests": {"exception": "True"},
    }
    runner = Runner(project="symm_guard", input_file=None,
                    log=str(tmp_path / "symm_guard.log"),
                    input_dict=config, silent=1, usempi=False)
    mol = runner.mol
    SinglePoint(mol).energy()
    oqp.fci_ao_integrals(mol)
    nbf = int(mol.data.get_basis()["nbf"])

    rng = np.random.default_rng(11)
    d = np.ascontiguousarray(np.eye(nbf))
    w = np.ascontiguousarray(np.eye(nbf))
    g = np.ascontiguousarray(rng.normal(size=(nbf,) * 4))   # no symmetry at all
    lib, ffi = _gradient_backend()
    status = int(lib.nevpt2_gradient(
        mol.data._data, int(nbf),
        ffi.cast("double *", d.ctypes.data),
        ffi.cast("double *", w.ctypes.data),
        ffi.cast("double *", g.ctypes.data)))
    assert status == -3


# --------------------------------------------------------------------------
# 6. input plumbing
# --------------------------------------------------------------------------
def test_pt2_gradient_option_parses_its_documented_spellings():
    from oqp.library.caspt2_dyall import _caspt2_options

    assert _caspt2_options({"pt2": {}}).gradient == "auto"
    assert _caspt2_options({"pt2": {"gradient": "analytic"}}).gradient == "analytic"
    assert _caspt2_options({"pt2": {"gradient": "ANALYTICAL"}}).gradient == "analytic"
    assert _caspt2_options({"pt2": {"gradient": "numerical"}}).gradient == "numerical"
    assert _caspt2_options({"pt2": {"gradient": "fd"}}).gradient == "numerical"
    with pytest.raises(ValueError, match="pt2.gradient"):
        _caspt2_options({"pt2": {"gradient": "sometimes"}})


def test_pt2_gradient_is_a_recognized_input_keyword():
    """A [pt2] key the parser rejects cannot be set from an input file."""
    from oqp.utils import oqp_input

    spec = oqp_input.parse_canonical_oqp(
        'caspt2/sto-3g geom="h4.xyz" grad pt2(gradient=analytic)')
    legacy = oqp_input.lower_to_legacy(spec, source_dir=".")
    assert legacy["pt2"]["gradient"] == "analytic"
