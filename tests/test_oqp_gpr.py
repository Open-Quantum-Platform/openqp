"""Derivative-GP model curvature checked without an electronic backend."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import cho_factor, cho_solve

PATH = Path(__file__).resolve().parents[1] / 'pyoqp/oqp/library/oqp_gpr.py'
SPEC = importlib.util.spec_from_file_location('_test_oqp_gpr', PATH)
GP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GP)


def kernel(x, y):
    return np.exp(-0.5 * np.sum((x - y)**2))


def test_derivative_covariance_against_scalar_kernel_differences():
    points = np.array([[-0.3, 0.1], [0.2, 0.4], [0.5, -0.2]])
    matrix = GP.derivative_covariance(points)
    step = 1e-4
    eye = np.eye(2) * step
    for i, x in enumerate(points):
        for j, y in enumerate(points):
            assert matrix[3*i, 3*j] == pytest.approx(kernel(x, y))
            for a in range(2):
                derivative = (kernel(x, y+eye[a]) - kernel(x, y-eye[a])) / (2*step)
                assert matrix[3*i, 3*j+1+a] == pytest.approx(derivative, abs=1e-8)
                for b in range(2):
                    mixed = (kernel(x+eye[a], y+eye[b]) - kernel(x+eye[a], y-eye[b])
                             - kernel(x-eye[a], y+eye[b]) + kernel(x-eye[a], y-eye[b])) / (4*step**2)
                    assert matrix[3*i+1+a, 3*j+1+b] == pytest.approx(mixed, abs=2e-8)
    assert np.allclose(matrix, matrix.T)
    assert np.linalg.eigvalsh(matrix)[0] > 0


def test_posterior_hessian_matches_finite_difference_of_posterior_energy():
    points = np.array([[-0.3, 0.1], [0.2, 0.4], [0.5, -0.2]])
    targets = np.arange(9) / 10
    matrix = GP.derivative_covariance(points) + np.eye(9) * 1e-8
    weights = cho_solve(cho_factor(matrix), targets)
    query = np.array([0.05, 0.03])
    analytic = GP.hessian_covariance(query, points) @ weights
    def energy(q):
        values = []
        for p in points:
            k = kernel(q, p)
            values.extend([k, *((q-p)*k)])
        return np.array(values) @ weights
    step = 1e-4
    eye = np.eye(2) * step
    for a in range(2):
        for b in range(2):
            finite = (energy(query+eye[a]+eye[b]) - energy(query+eye[a]-eye[b])
                      - energy(query-eye[a]+eye[b]) + energy(query-eye[a]-eye[b])) / (4*step**2)
            assert analytic[a, b] == pytest.approx(finite, rel=2e-5, abs=2e-5)


def quadratic_samples(model, hessian, points, origin=None):
    origin = np.zeros(len(points[0])) if origin is None else origin
    for p in points:
        x = p - origin
        model.add(p, 0.5*x @ hessian @ x, hessian @ x)


def test_recovers_sampled_negative_curvature_and_preserves_unobserved_space():
    model = GP.GPRCurvature(length_scale=0.5)
    exact = np.diag([-2., 5., 3.])
    prior = np.diag([0.5, 0.7, 0.9])
    points = [np.array([x, 0., 0.]) for x in (-0.15, -0.1, -0.05, 0.)]
    quadratic_samples(model, exact, points)
    delta = model.correction(prior)
    assert model.status == 'accepted'
    assert model.rank == 1
    assert abs((prior+delta)[0, 0] - exact[0, 0]) < 0.02
    assert np.allclose(delta[1:], 0, atol=1e-13)
    assert np.allclose(delta[:, 1:], 0, atol=1e-13)


def test_exact_quadratic_prior_is_unchanged():
    model = GP.GPRCurvature()
    exact = np.array([[-1., 0.3], [0.3, 2.]])
    quadratic_samples(model, exact, np.array([[-.2, 0], [0, .2], [.1, .1], [0, 0]]))
    assert np.max(np.abs(model.correction(exact))) < 1e-10


def test_covariant_under_cartesian_rotation_and_energy_origin_shift():
    h = np.diag([-2., 3., 1.])
    prior = np.eye(3) * 0.5
    points = np.array([[-.1, 0, 0], [0, -.1, 0], [.1, .05, 0], [0, 0, 0]])
    rot, _ = np.linalg.qr(np.array([[1., 2., 3.], [2., 5., 1.], [3., 1., 4.]]))
    a, b = GP.GPRCurvature(), GP.GPRCurvature()
    shift = np.array([5., -3., 7.])
    for x in points:
        e, g = 0.5*x @ h @ x, h @ x
        a.add(x, e, g)
        b.add(rot @ x + shift, e-123., rot @ g)
    assert np.allclose(b.correction(rot @ prior @ rot.T),
                       rot @ a.correction(prior) @ rot.T, atol=1e-6)


def test_history_ownership_bounds_duplicates_and_changed_objective():
    model = GP.GPRCurvature(history=3)
    for i in range(6):
        x = np.array([i/100.])
        g = 2*x
        model.add(x, float(x @ x), g)
        x[:] = 123.; g[:] = 999.
    assert len(model.samples) == 3
    assert model.samples[-1][0][0] == 0.05
    model.add([.05], .0025, [.1])
    assert len(model.samples) == 3
    model.add([.05], 50., [.1])
    assert len(model.samples) == 1
    assert model.correction(np.eye(1)) is None
    model.add([.04], np.nan, [.08])
    assert model.samples == []
    assert model.status == 'nonfinite_observation'


def test_distant_history_is_rejected_and_dimension_change_resets():
    model = GP.GPRCurvature(length_scale=.01)
    quadratic_samples(model, np.eye(1), np.array([[0.], [1.], [2.]]))
    assert model.correction(np.eye(1)) is None
    assert model.status == 'insufficient_local_history'
    model.add([0., 0.], 0., [0., 0.])
    assert len(model.samples) == 1


@pytest.mark.parametrize('kwargs', [dict(history=2), dict(history=21),
    dict(history=3.5), dict(history=True), dict(history=np.inf), dict(history=np.nan),
    dict(length_scale=True), dict(length_scale=0),
    dict(length_scale=np.inf), dict(length_scale=np.nan)])
def test_invalid_controls(kwargs):
    with pytest.raises(ValueError):
        GP.GPRCurvature(**kwargs)


def test_sparse_samples_reject_uncertain_curvature_despite_exact_observations():
    model = GP.GPRCurvature(length_scale=.5)
    quadratic_samples(model, np.array([[-2.]]), np.array([[-1.], [1.], [0.]]))
    assert model.correction(np.array([[.5]])) is None
    assert model.status == 'uncertain_curvature'
    assert model.variance_fraction > .05


def test_well_sampled_transverse_directions_cannot_hide_uncertain_reaction_curvature():
    model = GP.GPRCurvature(length_scale=.5)
    exact = np.diag([1., 1., -1.])
    points = np.array([[-.1, 0, 0], [.1, 0, 0], [0, -.1, 0],
                       [0, .1, 0], [0, 0, .7], [0, 0, 0]])
    quadratic_samples(model, exact, points)
    assert model.correction(np.eye(3)*.5) is None
    assert model.status == 'uncertain_curvature'
    assert model.variance_fraction > .05
