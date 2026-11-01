"""Image-dependent pair-potential initialization of molecular NEB paths.

Smidstrup et al., J. Chem. Phys. 140, 214106 (2014), Eqs. 2--4.
Coordinates are Bohr. This surrogate objective uses no electronic calls.
"""
from __future__ import annotations
import numpy as np
from oqp.library.oqp_neb import NEB


def idpp_energy_gradient(coordinates, target):
    """All-pair objective with distance-dependent d**-4 weights and gradient."""
    xyz = np.asarray(coordinates, dtype=float).reshape(-1, 3)
    i, j = np.triu_indices(len(xyz), 1)
    target = np.asarray(target, dtype=float)
    if (target.shape != i.shape or not np.isfinite(xyz).all()
            or not np.isfinite(target).all() or np.any(target <= 0)):
        raise ValueError("IDPP requires finite coordinates and positive pair targets")
    delta = xyz[i] - xyz[j]
    distance = np.linalg.norm(delta, axis=1)
    if np.any(distance < 1e-10):
        raise ValueError("IDPP cannot evaluate coincident atoms")
    error = distance - target
    energy = np.sum(error**2 / distance**4)
    radial = 2 * error / distance**4 - 4 * error**2 / distance**5
    pair_gradient = (radial / distance)[:, None] * delta
    gradient = np.zeros_like(xyz)
    np.add.at(gradient, i, pair_gradient)
    np.add.at(gradient, j, -pair_gradient)
    return float(energy), gradient.reshape(-1)


def _close_segments(images, targets, i, j):
    """Yield near-collisions along straight segments, including their interiors."""
    for k, (left, right) in enumerate(zip(images[:-1], images[1:])):
        left, right = left.reshape(-1, 3), right.reshape(-1, 3)
        relative = left[i] - left[j]
        motion = right[i] - right[j] - relative
        norm2 = np.einsum('ij,ij->i', motion, motion)
        fraction = np.clip(np.divide(-np.einsum('ij,ij->i', relative, motion),
                                    norm2, out=np.zeros_like(norm2), where=norm2 > 1e-24), 0, 1)
        nearest = np.linalg.norm(relative + fraction[:, None] * motion, axis=1)
        target = (1-fraction) * targets[k] + fraction * targets[k+1]
        for pair in np.flatnonzero(nearest < 0.1 * target):
            yield pair, motion[pair]


def _transverse_direction(motion):
    axis = np.eye(3)[int(np.argmin(np.abs(motion)))]
    normal = np.cross(motion, axis)
    if np.linalg.norm(normal) < 1e-12:
        normal = axis
    return normal / np.linalg.norm(normal)


def idpp_interpolate(images, maxiter=1000, fmax_tol=1e-3):
    """Optimize an owned band with fixed endpoints and no electronic Hessian.

    Pairs approaching coincidence at or between images receive a small
    deterministic transverse displacement before optimizing the singular
    pair potential. For a collinear exchange the Cartesian axis selects one
    of the equivalent transverse directions. Endpoints are never displaced.
    """
    band = NEB([np.asarray(x, dtype=float).reshape(-1).copy() for x in images], k_spring=0.1, climbing=False)
    natom = band.images[0].size // 3
    if natom < 2 or band.images[0].size != 3 * natom:
        raise ValueError("IDPP requires at least two atoms in Cartesian coordinates")
    i, j = np.triu_indices(natom, 1)
    def distances(x):
        xyz = x.reshape(-1, 3)
        return np.linalg.norm(xyz[i] - xyz[j], axis=1)
    r, p = band.images[0], band.images[-1]
    dr, dp = distances(r), distances(p)
    if np.any(dr < 1e-10) or np.any(dp < 1e-10):
        raise ValueError("IDPP endpoints contain coincident atoms")
    targets = [(1-f)*dr + f*dp for f in np.linspace(0, 1, band.nimage)]
    # A collinear exchange can cross between sampled images and have zero
    # perpendicular NEB force. Seed a smooth, consistently directed transverse
    # displacement for every affected pair, even when no image is close.
    crossing_pairs = dict(_close_segments(band.images, targets, i, j))
    for pair, motion in crossing_pairs.items():
        normal = _transverse_direction(motion)
        for k in range(1, band.nimage-1):
            xyz = band.images[k].reshape(-1, 3)
            shift = 0.1 * targets[k][pair] * np.sin(np.pi*k/(band.nimage-1)) * normal
            xyz[i[pair]] += shift
            xyz[j[pair]] -= shift
    # Resolve remaining image-local collisions, including interacting pairs.
    relative_motion = (p-r).reshape(-1, 3)
    for k in range(1, band.nimage-1):
        xyz = band.images[k].reshape(-1, 3)
        for _ in range(10):
            close = distances(band.images[k]) < 0.1 * targets[k]
            if not np.any(close):
                break
            for pair in np.flatnonzero(close):
                a, b = i[pair], j[pair]
                motion = relative_motion[a] - relative_motion[b]
                normal = _transverse_direction(motion)
                shift = 0.1 * targets[k][pair] * normal
                xyz[a] += shift
                xyz[b] -= shift
        if np.any(distances(band.images[k]) < 1e-10):
            raise ValueError("IDPP could not separate coincident intermediate atoms")

    class PairBand(NEB):
        def _evaluate(self, energy_gradient, indices):
            for k in indices:
                self.energies[k], self.gradients[k] = idpp_energy_gradient(self.images[k], targets[k])

    surrogate = PairBand(band.images, k_spring=0.1, climbing=False)
    # PairBand evaluates its own image-indexed objective, never this callback.
    def no_electronic_calls(_):
        raise AssertionError("IDPP must not call an electronic backend")
    result = surrogate.run(no_electronic_calls, maxiter=maxiter,
                         fmax_tol=fmax_tol, frms_tol=fmax_tol,
                         dt=0.05, dt_max=0.5, maxmove=0.1)
    # A zero projected band force alone cannot certify removal of an atom
    # exchange. Never pass a still-crossing band to the electronic NEB driver.
    if any(True for _ in _close_segments(result["images"], targets, i, j)):
        result["converged"] = False
    return result
