"""Local energy/gradient Gaussian-process corrections to model curvature.

The squared-exponential derivative covariance follows Rasmussen and Williams,
Gaussian Processes for Machine Learning (2006), section 9.4. This is a bounded,
displacement-subspace Hessian correction, not the GPR-dimer algorithm and not a
calculated molecular Hessian. No electronic-structure callback is used here.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular


def derivative_covariance(points):
    """Unit-length, unit-variance RBF covariance for interleaved E, gradient."""
    points = np.asarray(points, dtype=float)
    n, rank = points.shape
    width = rank + 1
    covariance = np.empty((n * width, n * width))
    identity = np.eye(rank)
    for i in range(n):
        for j in range(i + 1):
            delta = points[i] - points[j]
            kernel = np.exp(-0.5 * (delta @ delta))
            block = np.empty((width, width))
            block[0, 0] = 1.0
            block[0, 1:] = delta
            block[1:, 0] = -delta
            block[1:, 1:] = identity - np.outer(delta, delta)
            block *= kernel
            a, b = slice(i * width, (i + 1) * width), slice(j * width, (j + 1) * width)
            covariance[a, b] = block
            covariance[b, a] = block.T
    return covariance


def hessian_covariance(point, points):
    """Covariance of the query Hessian with interleaved E, gradient samples."""
    rank = len(point)
    identity = np.eye(rank)
    result = np.empty((rank, rank, len(points) * (rank + 1)))
    for i, sample in enumerate(points):
        delta = np.asarray(point) - sample
        kernel = np.exp(-0.5 * (delta @ delta))
        second = np.outer(delta, delta) - identity
        offset = i * (rank + 1)
        result[:, :, offset] = second * kernel
        for j in range(rank):
            third = (second * delta[j] - np.outer(identity[:, j], delta)
                     - np.outer(delta, identity[:, j]))
            result[:, :, offset + j + 1] = third * kernel
    return result


class GPRCurvature:
    """Keep a small owned history of Cartesian E/G samples (Hartree, Bohr).

    The current quadratic model is the GP mean. Only its residual curvature
    within the span of recent displacements is changed. Residual energy and
    gradient are reanchored to zero at the current sample; the electronic
    energy and gradient are never replaced by GP predictions. Convergence and
    the observed energy change in the trust ratio use electronic evaluations.
    """

    def __init__(self, history=8, length_scale=0.5):
        try:
            valid_history = (not isinstance(history, (bool, np.bool_))
                             and int(history) == history and 3 <= history <= 20)
        except (TypeError, ValueError, OverflowError):
            valid_history = False
        if not valid_history:
            raise ValueError("gpr_history must be an integer from 3 to 20")
        try:
            valid_scale = (not isinstance(length_scale, (bool, np.bool_))
                           and np.isfinite(length_scale) and 1e-3 <= length_scale <= 10)
        except (TypeError, ValueError):
            valid_scale = False
        if not valid_scale:
            raise ValueError("gpr_length_scale must be finite and between 0.001 and 10 Bohr")
        self.history = int(history)
        self.length_scale = float(length_scale)
        self.samples = []
        self.status = "insufficient_history"
        self.rank = 0
        self.variance_fraction = None

    def clear(self, reason="reset"):
        self.samples.clear()
        self.status = reason
        self.rank = 0
        self.variance_fraction = None

    def add(self, x, energy, gradient):
        x = np.asarray(x, dtype=float).reshape(-1)
        gradient = np.asarray(gradient, dtype=float).reshape(-1)
        if (x.shape != gradient.shape or not x.size or not np.isfinite(energy)
                or not np.isfinite(x).all() or not np.isfinite(gradient).all()):
            self.clear("nonfinite_observation")
            return
        if self.samples and self.samples[-1][0].shape != x.shape:
            self.clear("dimension_changed")
        # A repeated point must not give a nearly singular duplicate row. A
        # changed objective at the same point invalidates the old observations.
        for old_x, old_e, old_g in self.samples:
            if np.linalg.norm(x - old_x) < 1e-8 * self.length_scale:
                if abs(energy - old_e) > 1e-8 or np.linalg.norm(gradient - old_g) > 1e-6:
                    self.clear("objective_changed")
                else:
                    self.samples = [s for s in self.samples if s[0] is not old_x]
                break
        self.samples.append((x.copy(), float(energy), gradient.copy()))
        self.samples = self.samples[-self.history:]

    def correction(self, prior):
        """Return a Cartesian model-curvature correction, or None on rejection."""
        self.rank = 0
        self.variance_fraction = None
        self.status = "insufficient_history"
        if len(self.samples) < 3:
            return None
        x, energy, gradient = self.samples[-1]
        prior = np.asarray(prior, dtype=float)
        if prior.shape != (x.size, x.size) or not np.isfinite(prior).all():
            self.status = "invalid_prior"
            return None
        prior = 0.5 * (prior + prior.T)
        ell = self.length_scale
        nearby = [s for s in self.samples if np.linalg.norm(s[0] - x) <= 3 * ell]
        if len(nearby) < 3:
            self.status = "insufficient_local_history"
            return None
        displacement = np.array([s[0] - x for s in nearby])
        try:
            _, singular, vectors = np.linalg.svd(displacement, full_matrices=False)
            keep = singular > max(1e-6 * ell, 1e-5 * singular[0])
            basis = vectors[keep].T
            self.rank = basis.shape[1]
            if not self.rank:
                self.status = "dependent_displacements"
                return None
            points = displacement @ basis / ell
            e_residual = np.array([s[1] for s in nearby]) - energy
            e_residual -= displacement @ gradient
            e_residual -= 0.5 * np.einsum('ij,ij->i', displacement @ prior, displacement)
            g_residual = np.array([s[2] for s in nearby]) - gradient
            g_residual -= displacement @ prior
            observations = np.column_stack((e_residual, ell * g_residual @ basis)).reshape(-1)
            amplitude = max(float(np.max(np.abs(observations))), 1e-8)
            covariance = derivative_covariance(points)
            covariance.flat[::len(covariance) + 1] += 1e-10
            factor = cho_factor(covariance, lower=True, check_finite=True)
            weights = cho_solve(factor, observations / amplitude, check_finite=True)
            cross = hessian_covariance(np.zeros(self.rank), points)
            rows = cross.reshape(self.rank**2, -1)
            triangular = solve_triangular(factor[0], rows.T, lower=True)
            # Use the worst posterior/prior variance over symmetric Hessian
            # components. An average hides a poorly sampled reaction direction
            # among well sampled transverse directions. Orthonormal svec has
            # sqrt(2) off-diagonal weights and prior covariance 2I + u u.T,
            # where u is the vectorized identity; whitening is analytic.
            indices = np.triu_indices(self.rank)
            diagonal = (indices[0] == indices[1]).astype(float)
            symmetric_rows = indices[0] * self.rank + indices[1]
            weights_svec = np.where(diagonal, 1., np.sqrt(2.))
            observed = triangular[:, symmetric_rows].T * weights_svec[:, None]
            whitened = observed / np.sqrt(2.)
            whitened += ((1 / np.sqrt(self.rank + 2.) - 1 / np.sqrt(2.))
                         / self.rank * np.outer(diagonal, diagonal @ observed))
            posterior = np.eye(len(diagonal)) - whitened @ whitened.T
            self.variance_fraction = float(max(np.linalg.eigvalsh(posterior)[-1], 0.))
            if self.variance_fraction > 0.05:
                self.status = "uncertain_curvature"
                return None
            reduced = (rows @ weights).reshape(self.rank, self.rank) * amplitude / ell**2
            delta = basis @ (0.5 * (reduced + reduced.T)) @ basis.T
            if (not np.isfinite(delta).all()
                    or np.linalg.norm(delta, ord=2) > 10 * max(np.linalg.norm(prior, ord=2), 0.5)):
                self.status = "excessive_curvature"
                return None
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            self.status = "ill_conditioned_fit"
            return None
        self.status = "accepted"
        return delta
