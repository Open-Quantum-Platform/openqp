"""Geometry-only checks for the modified Lindh initial model."""
import importlib.util
import sys
import types
from pathlib import Path

import numpy as np
import pytest


def _load_model():
    lib = Path(__file__).resolve().parents[1] / 'pyoqp' / 'oqp' / 'library'
    names = ('oqp', 'oqp.library', 'oqp.library.oqp_coords', 'oqp.library.oqp_model_hessian')
    saved = {name: sys.modules.get(name) for name in names}
    try:
        for name in names[:2]:
            package = types.ModuleType(name)
            package.__path__ = []
            sys.modules[name] = package
        for name in names[2:]:
            spec = importlib.util.spec_from_file_location(name, lib / (name.rsplit('.', 1)[1] + '.py'))
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
        return module
    finally:
        for name, previous in saved.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous


MODEL = _load_model()
WATER = ([8, 1, 1], np.array([[0., 0., 0.], [1.8, 0., 0.], [-.45, 1.74, 0.]]))
AMMONIA = ([7, 1, 1, 1], np.array([[0., 0., .5], [1.8, 0., 0.], [-.9, 1.56, 0.], [-.9, -1.56, 0.]]))
CHAIN = ([6, 6, 6, 6], np.array([[0., 0., 0.], [2.8, 0., 0.], [3.8, 2.6, 0.], [5.6, 2.8, 2.1]]))


@pytest.mark.parametrize('atoms, xyz', [WATER, AMMONIA, CHAIN,
    (AMMONIA[0], AMMONIA[1] * [1, 1, 0])])
def test_raw_matches_finite_difference_of_frozen_internal_model(atoms, xyz):
    """At its origin, 1/2 sum k (q-q0)^2 has exactly B.T K B curvature."""
    z, x, distances = MODEL._validate_geometry(atoms, xyz)
    terms = list(MODEL._lindh_terms(z, x, distances))
    q0 = [primitive.value(x) for primitive, _ in terms]

    def energy(flat):
        e = 0.
        for (primitive, k), initial in zip(terms, q0):
            dq = primitive.value(flat.reshape(-1, 3)) - initial
            if primitive.kind == 'dihedral':
                dq = (dq + np.pi) % (2 * np.pi) - np.pi
            e += .5 * k * dq**2
        return e

    flat = x.ravel()
    step = np.eye(flat.size) * 2.e-4
    numeric = np.empty((flat.size, flat.size))
    for i in range(flat.size):
        for j in range(flat.size):
            numeric[i, j] = (energy(flat + step[i] + step[j]) - energy(flat + step[i] - step[j])
                             - energy(flat - step[i] + step[j]) + energy(flat - step[i] - step[j])) / (4 * (2.e-4)**2)
    np.testing.assert_allclose(MODEL.lindh_cartesian_hessian(atoms, x, eigenvalue_floor=0), numeric, atol=2.e-7, rtol=2.e-6)


@pytest.mark.parametrize('floor', [0., .05])
def test_rotation_translation_and_permutation_covariance(floor):
    atoms, x = AMMONIA
    rotation, _ = np.linalg.qr(np.array([[.2, .7, .3], [-.5, .4, .8], [.8, .1, -.6]]))
    h = MODEL.lindh_cartesian_hessian(atoms, x, eigenvalue_floor=floor)
    transform = np.kron(np.eye(len(atoms)), rotation.T)
    rotated = MODEL.lindh_cartesian_hessian(atoms, x @ rotation + np.array([2., -3., .4]), eigenvalue_floor=floor)
    np.testing.assert_allclose(rotated, transform @ h @ transform.T, atol=2.e-8)
    order = np.array([2, 0, 3, 1])
    indices = (3 * order[:, None] + np.arange(3)).ravel()
    permuted = MODEL.lindh_cartesian_hessian(np.array(atoms)[order], x[order], eigenvalue_floor=floor)
    np.testing.assert_allclose(permuted, h[np.ix_(indices, indices)], atol=2.e-8)


def test_raw_rigid_modes_and_positive_floor():
    atoms, x = WATER
    raw = MODEL.lindh_cartesian_hessian(atoms, x, eigenvalue_floor=0)
    eig = np.linalg.eigvalsh(raw)
    assert np.count_nonzero(eig > 1.e-8) == 3
    regularized = MODEL.lindh_cartesian_hessian(atoms, x)
    np.testing.assert_allclose(np.linalg.eigvalsh(regularized), np.maximum(eig, .05), atol=1.e-12)
    assert np.linalg.norm(raw @ np.tile([.3, -.2, .4], 3)) < 1.e-10


def test_planar_inversion_retains_positive_curvature():
    atoms, x = AMMONIA
    planar = x.copy()
    planar[:, 2] = 0.
    direction = np.zeros_like(planar)
    direction[:, 2] = [-3, 1, 1, 1]
    direction = direction.ravel() / np.linalg.norm(direction)
    h = MODEL.lindh_cartesian_hessian(atoms, planar)
    assert direction @ h @ direction >= .05 - 1.e-12


def test_linear_and_disconnected_geometries_are_bounded():
    for atoms, x in [([8, 6, 8], [[-2.2, 0, 0], [0, 0, 0], [2.2, 0, 0]]),
                     ([1, 1], [[0, 0, 0], [100, 0, 0]]), ([18], [[0, 0, 0]])]:
        h = MODEL.lindh_cartesian_hessian(atoms, x)
        assert np.all(np.isfinite(h))
        assert np.linalg.eigvalsh(h)[0] >= .05 - 1.e-12


def test_stretched_bond_has_smaller_model_force_constant():
    near = MODEL.lindh_cartesian_hessian([1, 1], [[0, 0, 0], [1.4, 0, 0]], eigenvalue_floor=0)
    far = MODEL.lindh_cartesian_hessian([1, 1], [[0, 0, 0], [2.4, 0, 0]], eigenvalue_floor=0)
    assert 0 < np.trace(far) < np.trace(near)


@pytest.mark.parametrize('atoms, x', [([], []), ([19], [0, 0, 0]), ([1.5], [0, 0, 0]),
    ([0], [0, 0, 0]), ([True], [0, 0, 0]), ([1], [0, 0]), ([1], [np.nan, 0, 0]),
    ([1, 1], [0, 0, 0, 0, 0, 0]), ([1, 1], [-1.e308, 0, 0, 1.e308, 0, 0]), ([1], [1j, 0, 0])])
def test_invalid_geometry_is_rejected(atoms, x):
    with pytest.raises(ValueError):
        MODEL.lindh_cartesian_hessian(atoms, x)


@pytest.mark.parametrize('floor', [-1., np.nan, np.inf, [0.1], 'bad', True])
def test_invalid_regularization_is_rejected(floor):
    with pytest.raises(ValueError):
        MODEL.lindh_cartesian_hessian(*WATER, eigenvalue_floor=floor)


@pytest.mark.parametrize('case, positive_eigenvalues', [
    (WATER, [0.2663695187801505, 1.193155456619044, 1.3772612264442028]),
    (AMMONIA, [0.19517377542873132, 0.30955532410619313, 0.30990257106713187,
               0.8238539663760722, 1.7458610335268605, 1.7462106002831685]),
    (CHAIN, [0.006220364520052757, 0.09205331265771964, 0.10604617538587659,
             0.8702794355486857, 1.052819483380577, 1.2764928653780643]),
])
def test_raw_spectrum_against_executable_reference(case, positive_eigenvalues):
    # pysisyphus a4ce10dd6d7fdcb3d813f1c730eb365d29041999, lindh_guess,
    # with exactly the same screened primitive set. Its analytic B matrix is
    # independent of OpenQP's derivative implementation. Hartree/bohr**2.
    h = MODEL.lindh_cartesian_hessian(*case, eigenvalue_floor=0)
    np.testing.assert_allclose(np.linalg.eigvalsh(h)[6:], positive_eigenvalues,
                               rtol=1.e-8, atol=1.e-10)
