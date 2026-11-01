"""Synchronous-transit geometry and tangent construction, in atomic units.

Peng and Schlegel, Isr. J. Chem. 33, 449 (1993), equations 1 and 8.
The nuclear Hessian is never evaluated here.
"""
from __future__ import annotations
import numpy as np
from oqp.library.oqp_coords import RedundantInternalCoordinates


def transit_tangent(coords, current, reactant, product):
    """Circle tangent through R, X, P; LST limit for collinear points."""
    q = coords.q(current)
    r = coords.q_displacement(coords.q(reactant), q)
    p = coords.q_displacement(coords.q(product), q)
    rr, pp = float(r @ r), float(p @ p)
    chord = coords.q_displacement(coords.q(product), coords.q(reactant))
    if rr > 1e-16 and pp > 1e-16:
        tangent = p / pp - r / rr
    else:
        tangent = chord.copy()
    norm = np.linalg.norm(tangent)
    if not np.isfinite(norm) or norm < 1e-12:
        tangent = chord.copy()
        norm = np.linalg.norm(tangent)
    if not np.isfinite(norm) or norm < 1e-12:
        raise ValueError("Synchronous transit requires distinct endpoint geometries")
    tangent /= norm
    if tangent @ chord < 0:
        tangent = -tangent
    return tangent


def internal_midpoint(atoms, reactant, product):
    """Interpolate redundant internals; return explicit Cartesian fallback flag."""
    r = np.asarray(reactant, dtype=float).reshape(-1)
    p = np.asarray(product, dtype=float).reshape(-1)
    if r.shape != p.shape or not np.isfinite([r, p]).all():
        raise ValueError("Endpoint coordinates must be finite and have equal shape")
    midpoint = (r + p) / 2
    if np.linalg.norm(r-p) < 1e-10:
        raise ValueError("Synchronous transit requires distinct endpoint geometries")
    ric = RedundantInternalCoordinates.from_geometry(atoms, midpoint)
    if ric is not None:
        target = ric.q(r) + 0.5 * ric.q_displacement(ric.q(p), ric.q(r))
        displacement = ric.q_displacement(target, ric.q(midpoint))
        result, ok = ric.back_transform(midpoint, displacement)
        if ok and np.isfinite(result).all():
            return result, "redundant_internal"
    return midpoint, "cartesian_fallback"
