"""End-to-end tests for the analytic CASPT2 / XMS-CASPT2 nuclear gradient.

The gradient itself contains no finite difference: every term is an explicit
derivative of the functional ``caspt2_dyall`` evaluates, with a multiplier
equation for each non-variational parameter (amplitudes, reference CI vector,
XMS state rotation, reference orbitals).  Finite differences appear here only as
the reference the analytic implementation is measured against, and in two
independent forms -- the production central-difference driver
(``wf_numgrad.wavefunction_numerical_gradient``) and a five-point stencil built
from fresh single-point runs, which shares no code with it.

The scope tests are as important as the numerical ones.  A PT2 option
combination the derivative does not cover must say so, not answer with a
plausible number: ``[pt2] gradient=analytic`` refuses, and the default
``gradient=auto`` falls back to central differences and records that it did.
"""
import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
    except Exception:
        return False
    if not bool(getattr(oqp, "BACKEND_AVAILABLE", False)):
        return False
    return hasattr(getattr(oqp, "lib", None), "caspt2_gradient")


needs_backend = pytest.mark.skipif(
    not _backend_available(),
    reason="native OQP backend with the caspt2_gradient entry point not built")

BOHR = 0.52917721092

#: Deliberately off-symmetry linear-ish H4, in Bohr: with a symmetric geometry
#: several gradient components vanish for reasons that have nothing to do with
#: the derivative being right.
H4_BOHR = np.array([0.0, 0.0, 0.0,
                    0.0, 0.0, 1.35,
                    0.18, 0.0, 2.80,
                    0.0, 0.12, 4.20])


def _make_runner(tmp_path, project, coords_bohr, atoms=("H", "H", "H", "H"), *,
                 method="caspt2", reference="casci", nroot=1,
                 target_roots=None, pt2_extra=None, runtype="energy",
                 grad=("0",), basis="sto-3g", active_electrons="2",
                 active_orbitals="2", frozen_core="1"):
    from oqp.pyoqp import Runner
    lines = []
    for sym, xyz in zip(atoms, np.asarray(coords_bohr, float).reshape(-1, 3)):
        a = xyz * BOHR
        lines.append(f"{sym} {a[0]:.12f} {a[1]:.12f} {a[2]:.12f}")
    pt2 = {"reference": reference, "frozen": "auto"}
    if target_roots is not None:
        pt2["target_roots"] = ",".join(str(r) for r in target_roots)
        pt2["nroot"] = str(len(target_roots))
    pt2.update(pt2_extra or {})
    return Runner(
        project=project, input_file=None, log=str(tmp_path / f"{project}.log"),
        input_dict={
            "input": {"system": "\n" + "\n".join(lines), "charge": "0",
                      "basis": basis, "method": method, "runtype": runtype},
            "guess": {"type": "hcore"},
            "scf": {"type": "rhf", "multiplicity": "1", "maxit": "120",
                    "conv": "1.0e-11", "save_molden": "False"},
            "cas": {"active_electrons": active_electrons,
                    "active_orbitals": active_orbitals,
                    "frozen_core": frozen_core, "orbital_source": "rhf",
                    "sort_orbitals": "energy", "max_det": "5000",
                    "max_memory": "512"},
            "ci": {"nroot": str(max(nroot, 1)), "solver": "dense",
                   "eig_tol": "1.0e-12", "integral_backend": "native",
                   "target_spin": "any"},
            # A loosely converged CASSCF reference puts its own convergence into
            # every displaced point, and the finite-difference side -- not the
            # analytic gradient -- then sets the floor: at the default
            # gradient_norm_tol the residual is 5e-7, at 1e-9 it is 2e-10.
            "casscf": {"gradient_norm_tol": "1.0e-9",
                       "energy_decrease_tol": "1.0e-12",
                       "max_macro_iterations": "200"},
            "pt2": pt2,
            "properties": {"scf_prop": "", "grad": ",".join(grad)},
            "symmetry": {"enabled": "False"},
            "tests": {"exception": "True"},
        },
        silent=1, usempi=False)


def _analytic(tmp_path, tag, coords, states=(0,), **kw):
    from oqp.library.caspt2_gradient import caspt2_analytic_gradient
    runner = _make_runner(tmp_path, tag, coords, **kw)
    runner.run(test_mod=True)
    return caspt2_analytic_gradient(runner.mol, list(states)), runner.mol


def _energies(tmp_path, tag, coords, **kw):
    runner = _make_runner(tmp_path, tag, coords, **kw)
    runner.run(test_mod=True)
    return np.asarray([float(e) for e in runner.mol.energies], dtype=float)


def _five_point(tmp_path, coords, state, h=2.0e-3, **kw):
    """Five-point stencil from fresh single-point runs (error O(h^4))."""
    coords = np.asarray(coords, float).reshape(-1)
    out = np.zeros(coords.size)
    for i in range(coords.size):
        for k, mult in ((2, -1.0), (1, 8.0), (-1, -8.0), (-2, 1.0)):
            disp = coords.copy()
            disp[i] += k * h
            out[i] += mult * _energies(tmp_path, f"fd{i}_{k}", disp, **kw)[state]
    return (out / (12.0 * h)).reshape(-1, 3)


@needs_backend
def test_single_state_caspt2_matches_five_point_finite_differences(tmp_path):
    """CASPT2/STO-3G on an off-symmetry H4: the analytic gradient reproduces a
    five-point difference of independently computed total energies."""
    grads, _mol = _analytic(tmp_path, "an", H4_BOHR)
    fd = _five_point(tmp_path, H4_BOHR, 0)
    assert np.max(np.abs(grads[0] - fd)) < 1.0e-7


@needs_backend
@pytest.mark.parametrize("shift", [
    {"level_shift": "0.15"},
    {"imaginary_shift": "0.20"},
    {"edshft": "0.05"},
])
def test_denominator_shifts_are_differentiated_exactly(tmp_path, shift):
    """Every regularization the option validator admits is carried exactly.

    A shifted amplitude is not stationary in the Hylleraas functional, so the
    shift enters the derivative through the divided-difference derivative of the
    denominator FUNCTION.  Dropping it (the usual approximation) would show up
    here as a first-order error, not a rounding one.
    """
    tag = "_".join(shift)
    grads, _mol = _analytic(tmp_path, f"an_{tag}", H4_BOHR, pt2_extra=shift)
    fd = _five_point(tmp_path, H4_BOHR, 0, pt2_extra=shift)
    assert np.max(np.abs(grads[0] - fd)) < 1.0e-7


@needs_backend
def test_xms_caspt2_state_gradients_match_finite_differences(tmp_path):
    """Both XMS-CASPT2 roots, including the state-rotation response."""
    grads, _mol = _analytic(tmp_path, "xms", H4_BOHR, states=(0, 1),
                            method="xms-caspt2", nroot=2, target_roots=(0, 1),
                            grad=("0", "1"))
    for state in (0, 1):
        fd = _five_point(tmp_path, H4_BOHR, state, method="xms-caspt2",
                         nroot=2, target_roots=(0, 1))
        assert np.max(np.abs(grads[state] - fd)) < 1.0e-7


@needs_backend
@pytest.mark.parametrize("method,states", [
    ("mrmp2", (0,)),
    ("mcqdpt2", (0, 1)),
    ("xmcqdpt2", (0, 1)),
])
def test_qdpt_family_on_its_default_engine(tmp_path, method, states):
    """The GAMESS-convention QDPT labels, on the DIRECT engine they default to.

    The gradient is always reconstructed on the dense path, so this also
    exercises the cross-check that refuses when the reconstruction and the
    reported energy disagree -- the direct and dense engines are different
    reductions of the same quantity.
    """
    kw = {}
    if len(states) > 1:
        kw = dict(nroot=len(states), target_roots=states,
                  grad=tuple(str(s) for s in states))
    grads, _mol = _analytic(tmp_path, f"q_{method}", H4_BOHR, states=states,
                            method=method, **kw)
    for state in states:
        fd = _five_point(tmp_path, H4_BOHR, state, method=method,
                         **{k: v for k, v in kw.items() if k != "grad"})
        assert np.max(np.abs(grads[state] - fd)) < 1.0e-7


@needs_backend
def test_casscf_reference_gradient_matches_finite_differences(tmp_path):
    """A state-specific CASSCF reference brings in the CASSCF orbital Hessian
    and the orbital-CI coupling, which an RHF-orbital CASCI reference does not."""
    grads, _mol = _analytic(tmp_path, "cas_an", H4_BOHR, reference="casscf")
    fd = _five_point(tmp_path, H4_BOHR, 0, h=1.0e-3, reference="casscf")
    assert np.max(np.abs(grads[0] - fd)) < 1.0e-7


@needs_backend
def test_state_averaged_casscf_reference_gradient(tmp_path):
    """A STATE-AVERAGED CASSCF reference: the constraint is the SA orbital
    gradient and the orbital-CI coupling runs over every averaged root.

    This is the combination with the most response machinery live at once --
    SA-CASSCF orbital relaxation, per-root CI relaxation, the XMS state rotation
    and the effective-Hamiltonian mixing vector.
    """
    grads, _mol = _analytic(tmp_path, "sa_an", H4_BOHR, states=(0, 1),
                            method="xms-caspt2", reference="casscf", nroot=2,
                            target_roots=(0, 1), grad=("0", "1"))
    for state in (0, 1):
        fd = _five_point(tmp_path, H4_BOHR, state, h=1.0e-3,
                         method="xms-caspt2", reference="casscf", nroot=2,
                         target_roots=(0, 1))
        assert np.max(np.abs(grads[state] - fd)) < 1.0e-7


@needs_backend
def test_gradient_is_translationally_and_rotationally_invariant(tmp_path):
    """Sum rules the derivative-integral contraction must satisfy exactly."""
    grads, _mol = _analytic(tmp_path, "inv", H4_BOHR)
    xyz = H4_BOHR.reshape(-1, 3)
    g = grads[0]
    assert np.max(np.abs(g.sum(axis=0))) < 1.0e-10
    com = xyz.mean(axis=0)
    assert np.max(np.abs(np.cross(xyz - com, g).sum(axis=0))) < 1.0e-9


@needs_backend
def test_frozen_core_is_unfolded_exactly(tmp_path):
    """LiH with the default PT2 frozen core: the fold is an exact algebraic
    transformation, so undoing it must leave the derivative unchanged."""
    coords = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 3.0])
    kw = dict(atoms=("Li", "H"))
    grads, mol = _analytic(tmp_path, "lih", coords, **kw)
    from oqp.library import caspt2_gradient as cg
    state = cg._build_state(mol)
    assert state.nfrozen == 1
    fd = _five_point(tmp_path, coords, 0, **kw)
    assert np.max(np.abs(grads[0] - fd)) < 1.0e-7
    assert np.max(np.abs(grads[0].sum(axis=0))) < 1.0e-10


@needs_backend
def test_d_shell_basis_exercises_the_cartesian_expansion(tmp_path):
    """LiH/cc-pVDZ: the first case with `d` functions.

    Every other case here is `s`/`p` only, where the spherical and Cartesian
    shells have the same size and the expansion in the derivative-ERI compute
    type is an identity. A `d` shell is 5 spherical against 6 Cartesian, so this
    is what actually tests `build_cart_density` on the factorized two-particle
    density.
    """
    coords = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 3.0])
    kw = dict(atoms=("Li", "H"), basis="cc-pvdz")
    grads, _mol = _analytic(tmp_path, "dsh", coords, **kw)
    fd = _five_point(tmp_path, coords, 0, **kw)
    assert np.max(np.abs(grads[0] - fd)) < 1.0e-7
    assert np.max(np.abs(grads[0].sum(axis=0))) < 1.0e-10


@needs_backend
@pytest.mark.parametrize("reference", ["casci", "casscf"])
def test_split_frozen_core_carries_its_own_response(tmp_path, reference):
    """BeH2: two inactive orbitals, of which the PT2 frozen core takes one.

    That split is not free.  The frozen orbitals are the LOWEST eigenvectors of
    the closed+active Fock restricted to the inactive span; the span is fixed by
    the inter-block constraints, but the Fock moves with the geometry, so the
    split moves with it and the energy is not invariant under the resulting
    within-inactive rotation.  Leaving that response out costs 8.5e-5 in the
    gradient and breaks rotational invariance at 2e-8 -- both invisible whenever
    nfrozen is 0 or equals ncore, which is every other case in this file.
    """
    beh2 = np.array([0.0, 0.0, 0.0,
                     0.0, 0.10, 2.55,
                     0.08, 0.0, -2.48])
    kw = dict(atoms=("Be", "H", "H"), basis="6-31g", reference=reference)
    grads, mol = _analytic(tmp_path, f"beh2_{reference}", beh2,
                           active_electrons="2", frozen_core="2", **kw)
    from oqp.library import caspt2_gradient as cg
    state = cg._build_state(mol)
    assert 0 < state.nfrozen < state.setup.ncore, "case does not split the core"

    fd = _five_point(tmp_path, beh2, 0, h=1.0e-3,
                     active_electrons="2", frozen_core="2", **kw)
    assert np.max(np.abs(grads[0] - fd)) < 1.0e-7
    xyz = beh2.reshape(-1, 3)
    com = xyz.mean(axis=0)
    assert np.max(np.abs(grads[0].sum(axis=0))) < 1.0e-10
    assert np.max(np.abs(np.cross(xyz - com, grads[0]).sum(axis=0))) < 1.0e-9


@needs_backend
def test_lagrangian_is_stationary_in_every_constrained_parameter(tmp_path):
    """Internal consistency: after the multipliers are solved, the Lagrangian
    has zero gradient along every constrained orbital rotation.  A missing term
    in the CI, XMS or orbital response leaves a residual here."""
    from oqp.library import caspt2_gradient as cg
    runner = _make_runner(tmp_path, "stat", H4_BOHR, method="xms-caspt2",
                          nroot=2, target_roots=(0, 1))
    runner.run(test_mod=True)
    state = cg._build_state(runner.mol)
    for target in range(state.nstate):
        relaxed = cg.relaxed_densities(state, target)
        assert relaxed.stationarity < 1.0e-9
        assert relaxed.orbital_residual < 1.0e-8


# --------------------------------------------------------------------- scope
@needs_backend
@pytest.mark.parametrize("pt2_extra,method,needle", [
    ({"h0": "dyall"}, "caspt2", "zeroth-order"),
    ({"h0": "dyall", "contraction": "strong"}, "caspt2", "contraction"),
    ({"ipea_shift": "0.30"}, "caspt2", "ipea_shift"),
    ({}, "ms-caspt2", "MULTI-SET"),
])
def test_out_of_scope_variants_fail_closed(tmp_path, pt2_extra, method, needle):
    """Every unsupported combination refuses with its own reason."""
    from oqp.library.caspt2_gradient import (
        CASPT2GradientNotImplemented, caspt2_analytic_gradient,
    )
    kw = {}
    if method == "ms-caspt2":
        kw = dict(nroot=2, target_roots=(0, 1))
    runner = _make_runner(tmp_path, "scope", H4_BOHR, method=method,
                          pt2_extra=pt2_extra, **kw)
    runner.run(test_mod=True)
    with pytest.raises(CASPT2GradientNotImplemented) as excinfo:
        caspt2_analytic_gradient(runner.mol, [0])
    assert needle in str(excinfo.value)


@needs_backend
def test_auto_falls_back_and_analytic_refuses(tmp_path, monkeypatch):
    """[pt2] gradient=auto routes an out-of-scope variant to central
    differences; gradient=analytic raises instead of quietly answering a
    different question.

    The numerical driver is replaced by a sentinel: what is under test is the
    ROUTING decision, and running 24 displaced PT2 energies to observe it would
    only make the test slow and couple it to an unrelated driver.
    """
    from oqp.library import wf_numgrad
    from oqp.library.caspt2_gradient import CASPT2GradientNotImplemented
    from oqp.library.single_point import Gradient

    called = {}

    def _sentinel(mol, grad_list, sp=None):
        called["yes"] = True
        return np.zeros((1, 4, 3))

    monkeypatch.setattr(wf_numgrad, "wavefunction_numerical_gradient", _sentinel)

    runner = _make_runner(tmp_path, "auto", H4_BOHR,
                          pt2_extra={"h0": "dyall", "gradient": "auto"})
    runner.run(test_mod=True)
    grads = np.asarray(Gradient(runner.mol).gradient(), dtype=float)
    assert called.get("yes"), "gradient=auto did not fall back to the numerical driver"
    assert grads.shape[-2:] == (4, 3)

    runner2 = _make_runner(tmp_path, "strict", H4_BOHR,
                           pt2_extra={"h0": "dyall", "gradient": "analytic"})
    runner2.run(test_mod=True)
    with pytest.raises(CASPT2GradientNotImplemented):
        Gradient(runner2.mol).gradient()


@needs_backend
def test_in_scope_variant_uses_the_analytic_route(tmp_path, monkeypatch):
    """The converse: a supported variant must NOT reach the numerical driver."""
    from oqp.library import wf_numgrad
    from oqp.library.single_point import Gradient

    def _sentinel(mol, grad_list, sp=None):          # pragma: no cover
        raise AssertionError("a supported variant fell back to central differences")

    monkeypatch.setattr(wf_numgrad, "wavefunction_numerical_gradient", _sentinel)
    runner = _make_runner(tmp_path, "inscope", H4_BOHR,
                          pt2_extra={"gradient": "auto"})
    runner.run(test_mod=True)
    grads = np.asarray(Gradient(runner.mol).gradient(), dtype=float)
    assert grads.shape[-2:] == (4, 3)
    assert np.max(np.abs(grads)) > 1.0e-3


@needs_backend
def test_imported_orbitals_are_refused(tmp_path):
    """CASCI orbitals read from a file are not a differentiable function of the
    geometry, so there is no orbital response to solve for.

    The gate is checked directly: reaching it through a full run would first
    hit the energy path's own "orbital_file is required" validation, which is a
    different (and correct) refusal about a different thing.
    """
    from oqp.library import caspt2_gradient as cg

    runner = _make_runner(tmp_path, "json_orb", H4_BOHR)
    runner.run(test_mod=True)
    setup = cg._caspt2_setup(runner.mol, run_reference=False)
    setup.orbital_source = "json:/tmp/orbitals.json"
    with pytest.raises(cg.CASPT2GradientNotImplemented) as excinfo:
        cg._gate_variant(setup.options, setup)
    assert "orbital_source" in str(excinfo.value)


def test_gradient_route_option_is_validated():
    """A misspelled route is an error, not a silent default."""
    from oqp.library.single_point import PT2_GRAD_METHODS
    from oqp.library.wf_numgrad import PT2_NUMGRAD_METHODS
    # The analytic dispatch set and the numerical one describe the same family;
    # a method reaching only one of them would silently lose its gradient.
    assert PT2_GRAD_METHODS == PT2_NUMGRAD_METHODS
