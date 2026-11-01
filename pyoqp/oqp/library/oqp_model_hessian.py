"""Geometry-dependent model curvature; no electronic Hessian is evaluated.

The modified Lindh model uses exp[alpha * (r_cov**2 - r**2)] screening
and stretch/bend/torsion factors 0.45/0.15/0.005 (atomic units), following
Lindh et al., Chem. Phys. Lett. 241, 423 (1995), DOI
10.1016/0009-2614(95)00646-L. Covalent-radius sums replace the original
period-pair reference distances, as in pysisyphus. This is explicitly a
modified Lindh model, not a reproduction of the original parameterization.

Only the first three periods (H--Ar) are supported. The returned Cartesian
matrix is B.T K B with an optional eigenvalue floor. The floor is optimizer
regularization, not part of the Lindh model and not molecular curvature.
"""
from itertools import combinations

import numpy as np

from .oqp_coords import Angle, Bond, Dihedral

# Cordero covalent radii in bohr rounded to match the tested modified-Lindh
# reference convention; H uses 0.4 Angstrom, rather than Cordero's 0.31.
# pysisyphus a4ce10dd6d7fdcb3d813f1c730eb365d29041999, elem_data.py.
_RADII = np.array([
    0.0, 0.7561, 0.5291, 2.4188, 1.8141, 1.5874, 1.4362, 1.3417,
    1.2472, 1.0771, 1.0960, 3.1369, 2.6645, 2.2866, 2.0976, 2.0220,
    1.9842, 1.9275, 2.0031,
])
# Discard negligible screened pair interactions. Linear angular coordinates
# are omitted because their scalar derivatives are singular; the explicitly
# separate spectral floor supplies bounded curvature in omitted directions.
_RHO_CUTOFF = 1.0e-4
_MIN_SINE = 0.1


def _validate_geometry(atoms, x):
    z = np.asarray(atoms)
    if (z.ndim != 1 or z.size == 0 or z.dtype.kind not in 'iuf'
            or not np.all(np.isfinite(z)) or np.any(z != np.floor(z))
            or np.any(z < 1) or np.any(z > 18)):
        raise ValueError('Modified Lindh model requires integer atomic numbers H--Ar')
    xyz = np.asarray(x)
    if xyz.dtype.kind not in 'iuf':
        raise ValueError('Lindh geometry must contain real numeric coordinates')
    xyz = xyz.astype(float)
    if xyz.shape not in ((3 * z.size,), (z.size, 3)) or not np.all(np.isfinite(xyz)):
        raise ValueError('Lindh geometry must contain finite 3N Cartesian coordinates in bohr')
    xyz = xyz.reshape(-1, 3)
    # Overflow in a coordinate difference or norm is an invalid geometry,
    # not permission to quietly replace its model by the regularization.
    with np.errstate(over='ignore', invalid='ignore'):
        delta = xyz[:, None, :] - xyz[None, :, :]
        distances = np.linalg.norm(delta, axis=-1)
    pairs = distances[np.triu_indices(z.size, 1)]
    if np.any(~np.isfinite(pairs)) or np.any(pairs < 1.0e-6):
        raise ValueError('Lindh geometry contains coincident atoms or unrepresentable distances')
    return z.astype(int), xyz, distances


def _lindh_terms(z, xyz, distances):
    """Yield screened (primitive, force constant) pairs in atomic units."""
    first = z <= 2
    alpha = np.where(first[:, None] & first[None, :], 1.0,
                     np.where(first[:, None] | first[None, :], 0.3949, 0.28))
    radii_sum = _RADII[z, None] + _RADII[z][None, :]
    with np.errstate(over='ignore', under='ignore'):
        rho = np.exp(alpha * (radii_sum**2 - distances**2))
    np.fill_diagonal(rho, 0.0)
    neighbors = [np.flatnonzero(row >= _RHO_CUTOFF).tolist() for row in rho]
    bends = set()
    for j, adjacent in enumerate(neighbors):
        for i in adjacent:
            if i < j:
                yield Bond(i, j), 0.45 * rho[i, j]
        for i, k in combinations(adjacent, 2):
            u = (xyz[i] - xyz[j]) / distances[i, j]
            v = (xyz[k] - xyz[j]) / distances[k, j]
            if np.linalg.norm(np.cross(u, v)) < _MIN_SINE:
                continue
            bends.add((min(i, k), j, max(i, k)))
            yield Angle(i, j, k), 0.15 * rho[i, j] * rho[j, k]
    for j, adjacent in enumerate(neighbors):
        for k in adjacent:
            if j >= k:
                continue
            for i in adjacent:
                if i == k or (min(i, k), j, max(i, k)) not in bends:
                    continue
                for l in neighbors[k]:
                    if l in (i, j) or (min(j, l), k, max(j, l)) not in bends:
                        continue
                    yield Dihedral(i, j, k, l), 0.005 * rho[i, j] * rho[j, k] * rho[k, l]


def lindh_cartesian_hessian(atoms, x, *, eigenvalue_floor=0.05):
    """Return a regularized modified-Lindh model in Hartree/bohr**2.

    ``atoms`` are integer atomic numbers; ``x`` is flat or (N, 3), in bohr.
    Set ``eigenvalue_floor=0`` to obtain the raw positive-semidefinite B.T K B
    model. The default floor (0.05 Hartree/bohr**2) also regularizes translations,
    rotations, and rank deficiencies such as planar inversion or linear bends.
    It leaves all model eigenvalues above the floor unchanged. The optimizer
    must still remove rigid displacements and choose its uphill TS direction.
    Neither matrix is an electronic Hessian or suitable for frequencies.

    Pair screening and exclusion of nearly linear angles make this an initial
    curvature estimate, not a globally smooth potential-energy Hessian.
    """
    if (not isinstance(eigenvalue_floor, (int, float, np.integer, np.floating))
            or isinstance(eigenvalue_floor, (bool, np.bool_))
            or not np.isfinite(eigenvalue_floor) or eigenvalue_floor < 0):
        raise ValueError('Lindh eigenvalue floor must be finite and nonnegative')
    z, xyz, distances = _validate_geometry(atoms, x)
    hessian = np.zeros((xyz.size, xyz.size))
    for primitive, force_constant in _lindh_terms(z, xyz, distances):
        derivatives = primitive.derivatives(xyz)
        indices = np.array([3 * atom + axis for atom, _ in derivatives for axis in range(3)])
        row = np.concatenate([derivative for _, derivative in derivatives])
        hessian[np.ix_(indices, indices)] += force_constant * np.outer(row, row)
    if not np.all(np.isfinite(hessian)):
        raise ValueError('Modified Lindh model produced nonfinite curvature')
    if eigenvalue_floor:
        values, vectors = np.linalg.eigh(hessian)
        # dot avoids spurious Accelerate matmul floating-point warnings.
        hessian = np.dot(vectors * np.maximum(values, eigenvalue_floor), vectors.T)
    hessian = 0.5 * hessian + 0.5 * hessian.T
    if not np.all(np.isfinite(hessian)):
        raise ValueError('Lindh regularization produced nonfinite curvature')
    return hessian
