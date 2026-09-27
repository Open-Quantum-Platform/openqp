"""Saddle convergence and pair-potential checks without a native backend."""
import importlib.util
import sys
import types
from pathlib import Path
import numpy as np
import pytest

LIB = Path(__file__).resolve().parents[1] / 'pyoqp/oqp/library'
def _modules():
    names = ['oqp', 'oqp.library'] + ['oqp.library.' + n for n in
             ('oqp_coords', 'oqp_transit', 'oqp_model_hessian', 'oqp_gpr', 'oqp_engine', 'oqp_neb', 'oqp_idpp')]
    saved = {n: sys.modules.get(n) for n in names}
    modules = []
    try:
        for n in names[:2]:
            package = types.ModuleType(n); package.__path__ = []
            sys.modules[n] = package
        for n in names[2:]:
            spec = importlib.util.spec_from_file_location(n, LIB / (n.split('.')[-1]+'.py'))
            module = importlib.util.module_from_spec(spec)
            sys.modules[n] = module; spec.loader.exec_module(module)
            modules.append(module)
    finally:
        for n, old in saved.items():
            if old is None: sys.modules.pop(n, None)
            else: sys.modules[n] = old
    return modules
COORDS, TRANSIT, MODEL, GPR, ENGINE, NEB, IDPP = _modules()


def curved_double_well(x):
    a, b, z = x
    residual = b - 0.4*(1-a*a)
    return (a*a-1)**2 + 2*residual**2 + 3*z*z, np.array([
        4*a*(a*a-1) + 3.2*a*residual, 4*residual, 6*z])


@pytest.mark.parametrize('guess', [[0., 0., 0.], [0.2, 0.25, 0.05]])
def test_transit_converges_to_curved_saddle_using_only_energy_gradient(guess):
    ends = (np.array([-1., 0., 0.]), np.array([1., 0., 0.]))
    engine = ENGINE.OQPEngine([1], guess, mode='ts', coordsys='cartesian',
                             maxiter=100, transit_endpoints=ends)
    evaluations = []
    def eg(x):
        e, g = curved_double_well(x); evaluations.append((x.copy(), g.copy()))
        return e, g
    def stop():
        if np.max(np.abs(evaluations[-1][1])) < 1e-7: raise StopIteration
    result = engine.run(eg, stop)
    assert np.allclose(result, [0, .4, 0], atol=1e-6)
    assert np.max(np.abs(evaluations[-1][1])) < 1e-7
    assert len(evaluations) < 60
    # Own endpoint copies; repeated calls and caller mutation cannot alter them.
    ends[0][:] = 999
    assert engine.transit_endpoints[0][0] == -1


def test_circle_tangent_is_perpendicular_to_radius():
    cart = COORDS.CartesianCoordinates(1)
    point = np.array([0., 1., 0.])
    tangent = TRANSIT.transit_tangent(cart, point, [-1, 0, 0], [1, 0, 0])
    assert np.allclose(tangent, [1, 0, 0])
    assert tangent @ point == pytest.approx(0)


@pytest.mark.parametrize('endpoints', [([0,0,0],[0,0,0]), ([0,0,0],[np.nan,0,0]), ([0,0],[1,0])])
def test_invalid_transit_endpoints_rejected(endpoints):
    with pytest.raises(ValueError):
        ENGINE.OQPEngine([1], [0,0,0], mode='ts', transit_endpoints=endpoints)


def test_idpp_gradient_matches_finite_difference_and_zero_total_force():
    x = np.array([[0.,0.,0.], [1.2,.1,0.], [.4,1.3,.2]])
    target = np.array([1., 1.5, 1.7])
    e, g = IDPP.idpp_energy_gradient(x, target)
    numerical = np.empty(x.size)
    for k in range(x.size):
        plus, minus = x.reshape(-1).copy(), x.reshape(-1).copy()
        plus[k] += 1e-6; minus[k] -= 1e-6
        numerical[k] = (IDPP.idpp_energy_gradient(plus,target)[0] -
                        IDPP.idpp_energy_gradient(minus,target)[0])/2e-6
    assert np.allclose(g, numerical, atol=1e-9)
    assert np.allclose(g.reshape(-1,3).sum(axis=0), 0, atol=1e-14)
    assert IDPP.idpp_energy_gradient(x+12,target)[0] == pytest.approx(e)


@pytest.mark.parametrize("nimage", [4, 6, 7, 8])
def test_idpp_resolves_crossing_without_changing_endpoints_or_input(nimage):
    r = np.array([-1.,0,0, 1.,0,0])
    p = -r
    images = [(1-f)*r + f*p for f in np.linspace(0,1,nimage)]
    original = [x.copy() for x in images]
    result = IDPP.idpp_interpolate(images, maxiter=1500)
    assert result['converged']
    assert np.array_equal(result['images'][0],r)
    assert np.array_equal(result['images'][-1],p)
    assert all(np.array_equal(x,y) for x,y in zip(images,original))
    assert min(np.linalg.norm(x.reshape(-1,3)[0]-x.reshape(-1,3)[1])
               for x in result['images']) > 1.9
    # Verify the continuous straight segments, not just sampled distances.
    for left, right in zip(result['images'][:-1], result['images'][1:]):
        a = left[:3] - left[3:]
        d = right[:3] - right[3:] - a
        f = np.clip(-a @ d / (d @ d), 0., 1.)
        assert np.linalg.norm(a + f*d) > .2


def test_idpp_rejects_coincident_endpoints():
    with pytest.raises(ValueError, match='endpoints'):
        IDPP.idpp_interpolate([np.zeros(6)]*3)


def test_midpoint_preserves_symmetric_ammonia_inversion():
    r = np.array([[0,0,.7],[1.8,0,0],[-.9,1.56,0],[-.9,-1.56,0]])
    p = r.copy(); p[0,2] *= -1
    midpoint, kind = TRANSIT.internal_midpoint([7,1,1,1],r,p)
    assert np.isfinite(midpoint).all()
    assert kind in {'redundant_internal', 'cartesian_fallback'}
    assert abs(midpoint.reshape(-1,3)[0,2]) < 1e-8


def test_transit_rejects_calculated_hessian():
    with pytest.raises(ValueError, match='model Hessian'):
        ENGINE.OQPEngine([1], [0,0,0], mode='ts', initial_hessian=np.eye(3),
                         transit_endpoints=([-1,0,0],[1,0,0]))


def ammonia_example_geometries():
    folder = LIB.parents[2] / 'examples/OPT'
    # XYZ examples use Angstrom; the engine takes Bohr.
    return [np.loadtxt(folder / ('NH3_QST_' + name + '.xyz'),
                       skiprows=2, usecols=(1, 2, 3)).reshape(-1) / 0.529177210903
            for name in ('reactant', 'product', 'guess')]


@pytest.mark.parametrize('coordsys', ['auto', 'dlc', 'ric'])
def test_ammonia_qst3_indistinguishable_internal_endpoints(coordsys):
    r, p, guess = ammonia_example_geometries()
    engine = ENGINE.OQPEngine(
        [7, 1, 1, 1], guess, mode='ts', coordsys=coordsys,
        project_global_rigid_modes=True, transit_endpoints=(r, p))
    # Exercise the search before checking its coordinates: the old engine
    # raises ValueError on this valid inversion path at the first step.
    midpoint = (r + p) / 2
    direction = (p - r) / np.linalg.norm(p - r)
    half_length = np.linalg.norm(p - r) / 2
    gradients = []

    def double_well(x):
        displacement = x - midpoint
        height = displacement @ direction
        transverse = displacement - height * direction
        energy = (height**2 - half_length**2)**2 + transverse @ transverse
        gradient = (4 * height * (height**2 - half_length**2) * direction
                    + 2 * transverse)
        gradients.append(gradient)
        return energy, gradient

    def stop():
        if np.max(np.abs(gradients[-1])) < 1e-8:
            raise StopIteration

    result = engine.run(double_well, stop)
    assert np.allclose(result, midpoint, atol=1e-7)
    assert len(gradients) < 30
    assert isinstance(engine.coords, COORDS.CartesianCoordinates)
    assert engine.coordsys.endswith('->CART(fallback)')
    assert engine.H.shape == (guess.size, guess.size)


@pytest.mark.parametrize('transit', [False, True])
def test_ammonia_distinct_internal_endpoints_keep_dlc(transit):
    r, p, guess = ammonia_example_geometries()
    p[2] *= 1.2
    engine = ENGINE.OQPEngine(
        [7, 1, 1, 1], guess, mode='ts', project_global_rigid_modes=True,
        transit_endpoints=(r, p) if transit else None)
    assert isinstance(engine.coords, COORDS.DelocalizedInternalCoordinates)
    assert engine.coordsys == 'DLC'
    if transit:
        tangent = TRANSIT.transit_tangent(engine.coords, guess, r, p)
        assert np.linalg.norm(tangent) == pytest.approx(1)


@pytest.mark.parametrize('guess,accepted', [([0.2, 0.25, 0.05], True),
                                           ([-0.3, 0.1, -0.1], False)])
def test_gpr_qst_uses_actual_gradients_and_converges_to_same_saddle(guess, accepted):
    messages = []
    engine = ENGINE.OQPEngine(
        [1], guess, mode='ts', coordsys='cartesian', maxiter=100,
        transit_endpoints=([-1., 0., 0.], [1., 0., 0.]),
        hessian_update='gpr', logger=messages.append)
    evaluations = []

    def eg(x):
        e, g = curved_double_well(x)
        evaluations.append((x.copy(), e, g.copy()))
        return e, g

    def stop():
        if np.max(np.abs(evaluations[-1][2])) < 1e-7:
            raise StopIteration

    result = engine.run(eg, stop)
    assert np.allclose(result, [0, .4, 0], atol=1e-6)
    assert (engine.gpr_accepted_steps > 0) == accepted
    assert engine.gpr_accepted_steps + engine.gpr_rejected_steps == len(evaluations)-1
    assert all(any(np.array_equal(x, sample[0]) and e == sample[1]
                   and np.array_equal(g, sample[2]) for sample in evaluations)
               for x, e, g in engine.gpr.samples)
    assert any('accepted' in line for line in messages) == accepted
    if not accepted:
        assert any('uncertain_curvature' in line for line in messages)
    assert len(evaluations) < 60


def test_gpr_failure_retains_finite_quasi_newton_step(monkeypatch):
    engine = ENGINE.OQPEngine([1], [.2, .1, .1], mode='ts', coordsys='cart',
                             hessian_update='gpr')
    def reject(prior):
        engine.gpr.status = 'ill_conditioned_fit'
        return None
    monkeypatch.setattr(engine.gpr, 'correction', reject)
    e, g = curved_double_well(engine.x)
    step = engine._take_step(e, g)
    assert np.isfinite(step).all() and np.isfinite(engine.H).all()
    assert engine.gpr_rejected_steps == 1
    engine.rebase_previous_objective(e, g)
    assert not engine.gpr.samples


@pytest.mark.parametrize('coordsys', ['cart', 'dlc', 'ric', 'tric'])
def test_lindh_and_gpr_curvature_transform_together(coordsys):
    x = np.array([0., 0., 0., 1.8, 0., 0., -.45, 1.7, .1])
    engine = ENGINE.OQPEngine([8, 1, 1], x, coordsys=coordsys,
                             model_hessian='lindh', hessian_update='gpr')
    initial = engine.H.copy()
    b = engine.coords.b_matrix(x)
    prior = b.T @ initial @ b
    direction = np.arange(1., 10.)
    direction /= np.linalg.norm(direction)
    exact = prior + .3 * np.outer(direction, direction)
    # Analytic quadratic observations, in a single sampled displacement
    # direction. No electronic-structure Hessian or callback is involved.
    for offset in [-.15, -.1, -.05]:
        dx = offset * direction
        engine.gpr.add(x + dx, .5 * dx @ exact @ dx, exact @ dx)
    delta = .3 * np.outer(direction, direction)
    expected = engine._transform_initial_hessian(prior + delta, engine._guess_hessian(),
                                                 np.zeros(x.size))
    engine._apply_gpr_curvature(0., np.zeros(x.size), b)
    assert engine.gpr.status == 'accepted'
    assert np.allclose(engine.H, expected, atol=2e-3)
    assert not np.allclose(engine.H, initial, atol=1e-6)
    assert np.isfinite(engine._take_step(0., np.zeros(x.size))).all()


def test_lindh_gpr_coordinate_recovery_discards_history_and_uses_selected_model():
    x = np.array([0., 0., 0., 1.8, 0., 0., -.45, 1.7, .1])
    engine = ENGINE.OQPEngine([8, 1, 1], x, coordsys='dlc', model_hessian='lindh',
                             hessian_update='gpr')
    engine.gpr.add(x, 0., np.ones(x.size))
    step = engine._bounded_cartesian_recovery(np.ones(x.size))
    assert np.isfinite(step).all()
    assert isinstance(engine.coords, COORDS.CartesianCoordinates)
    assert not engine.gpr.samples
    assert np.allclose(engine.H, MODEL.lindh_cartesian_hessian([8,1,1], x))


@pytest.mark.parametrize('coordsys', ['cart', 'dlc', 'ric', 'tric'])
def test_nonstationary_gpr_transformation_keeps_coordinate_curvature_consistent(coordsys, monkeypatch):
    x = np.array([0., 0., 0., 1.8, 0., 0., -.45, 1.7, .1])
    engine = ENGINE.OQPEngine([8, 1, 1], x, coordsys=coordsys,
                             model_hessian='lindh', hessian_update='gpr')
    original = engine.H.copy()
    b = engine.coords.b_matrix(x)
    pinv = np.linalg.pinv(b, rcond=1e-6)
    gradient = np.arange(9) / 100.
    direction = np.arange(1., 10.) / 10.
    delta = .1 * np.outer(direction, direction)
    def correct(prior):
        assert np.allclose(prior, b.T @ original @ b + engine._coordinate_curvature(gradient))
        engine.gpr.status = 'accepted'
        return delta
    monkeypatch.setattr(engine.gpr, 'correction', correct)
    engine._apply_gpr_curvature(0., gradient, b)
    assert engine.gpr_accepted_steps == 1
    assert np.allclose(engine.H, original + pinv.T @ delta @ pinv, atol=1e-10)


@pytest.mark.parametrize('mode,transit', [('min', False), ('ts', False), ('ts', True)])
def test_default_model_uses_lindh_with_or_without_qst(mode, transit):
    r, p, guess = ammonia_example_geometries()
    messages = []
    engine = ENGINE.OQPEngine([7, 1, 1, 1], guess, mode=mode,
        project_global_rigid_modes=True, logger=messages.append,
        transit_endpoints=(r, p) if transit else None)
    assert engine.requested_model_hessian == 'auto'
    assert engine.model_hessian == 'lindh'
    assert engine.hessian_update == 'auto' and engine.gpr is None
    assert np.isfinite(engine.H).all()
    assert any('auto -> lindh' in line for line in messages)
    explicit = ENGINE.OQPEngine([7, 1, 1, 1], guess, mode=mode,
        project_global_rigid_modes=True, model_hessian='lindh',
        transit_endpoints=(r, p) if transit else None)
    assert np.allclose(engine.H, explicit.H)


@pytest.mark.parametrize('options', [
    dict(atoms=[35, 1]), dict(mode='meci'), dict(project_global_rigid_modes=False),
    dict(frozen_distances=[(1, 2)]), dict(initial_hessian=np.eye(6)),
])
def test_automatic_model_preserves_unsupported_existing_paths(options):
    kwargs = dict(atoms=[1, 1], x0=[0, 0, 0, 1.4, 0, 0], coordsys='cart',
                  project_global_rigid_modes=True)
    kwargs.update(options)
    engine = ENGINE.OQPEngine(**kwargs)
    assert engine.model_hessian == 'constant'
    assert engine.gpr is None
    assert np.isfinite(engine.H).all()


def test_explicit_constant_overrides_automatic_lindh():
    r, _, _ = ammonia_example_geometries()
    engine = ENGINE.OQPEngine([7, 1, 1, 1], r, model_hessian='constant',
                             project_global_rigid_modes=True)
    assert engine.model_hessian == 'constant'
    assert np.allclose(engine.H, engine.coords.guess_hessian(r))


@pytest.mark.parametrize('early_gradient', ['finite', 'zero', 'nonfinite'])
def test_qst_cartesian_recovery_invalidates_internal_tangent(monkeypatch, early_gradient):
    r, p, guess = ammonia_example_geometries()
    p[2] *= 1.2  # Distinct internal endpoints retain a six-dimensional DLC.
    engine = ENGINE.OQPEngine([7,1,1,1], guess, mode='ts', coordsys='dlc',
        project_global_rigid_modes=True, transit_endpoints=(r,p))
    assert isinstance(engine.coords, COORDS.DelocalizedInternalCoordinates)
    monkeypatch.setattr(engine.coords, 'back_transform',
                        lambda x, dq: (np.full_like(x, np.nan), True))
    g = np.linspace(-.1, .1, len(guess))
    if early_gradient == 'finite':
        new = engine._take_step(0., g)  # Populates tangent, then forces recovery.
    else:
        engine._transit_direction = TRANSIT.transit_tangent(engine.coords, guess, r, p)
        engine._followed = engine._transit_direction.copy()
        engine._transit_steps = 1
        new = engine._bounded_cartesian_recovery(
            np.zeros_like(g) if early_gradient == 'zero' else np.full_like(g, np.nan))
    assert isinstance(engine.coords, COORDS.CartesianCoordinates)
    assert engine._transit_direction is None
    assert engine._followed is None or engine._followed.shape == guess.shape
    assert engine._prev is None and np.isfinite(new).all()
    assert np.sqrt(np.mean((new-guess)**2)) <= engine.trust * (1+1e-12)
    # Repeated recovery and the next regular QST step must use Cartesian modes.
    assert np.isfinite(engine._bounded_cartesian_recovery(g)).all()
    assert np.isfinite(engine._take_step(0., g)).all()
    assert engine._transit_direction.shape == guess.shape


def test_idpp_rejects_false_force_convergence_with_crossing_between_images(monkeypatch):
    r = np.array([-1., 0, 0, 1., 0, 0])
    images = [(1-f)*r - f*r for f in np.linspace(0, 1, 6)]
    # Simulate projected-force convergence that still leaves a crossing band.
    monkeypatch.setattr(NEB.NEB, 'run', lambda *args, **kwargs:
                        {'converged': True, 'images': [x.copy() for x in images]})
    result = IDPP.idpp_interpolate(images)
    assert not result['converged']


def test_idpp_does_not_perturb_a_non_crossing_zero_objective_band():
    r = np.array([-1., 0, 0, 1., 0, 0])
    translation = np.array([0., 1., .3, 0., 1., .3])
    images = [r + f*translation for f in np.linspace(0, 1, 6)]
    result = IDPP.idpp_interpolate(images)
    assert result['converged']
    assert all(np.array_equal(a,b) for a,b in zip(result['images'], images))
