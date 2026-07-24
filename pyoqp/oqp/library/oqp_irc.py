"""OQP intrinsic reaction coordinate (IRC) for OpenQP.

Backend-agnostic (NumPy only) Gonzalez--Schlegel mass-weighted steepest-descent
reaction path (Gonzalez & Schlegel, J. Phys. Chem. 94, 5523 (1990)).

The IRC is the steepest-descent path in mass-weighted Cartesian coordinates
q = M^{1/2} x.  Each step of mass-weighted arc length ``s`` is taken with the
Gonzalez--Schlegel constrained scheme:

  1. pivot point  q* = q_k - (s/2) g_hat   (g_hat = mass-weighted gradient dir);
  2. find q_{k+1} on the hypersphere |q - q*| = s/2 that minimizes E, i.e. the
     point whose gradient is radial (no tangential component).

The first point is displaced from the transition state along the
imaginary-frequency mode (mass-weighted), in the requested forward/backward
direction.  The caller supplies ``energy_gradient(x_flat) -> (E, g_flat)`` in
atomic units, so this module is testable on analytic surfaces without the OQP
core.
"""

from __future__ import annotations

import numpy as np


class IRC:
    def __init__(self, mass, x_ts, init_direction_mw, step=0.10, sign=1,
                 logger=None):
        # mass: per-atom masses (natom,) in amu; x_ts: flat Cartesian (Bohr).
        mass = np.asarray(mass, dtype=float).reshape(-1)
        if mass.size == 0 or not np.all(np.isfinite(mass)) or np.any(mass <= 0.0):
            raise ValueError("IRC masses must be finite and positive")
        self.sqm = np.sqrt(np.repeat(mass, 3))          # M^{1/2}, length 3N
        self.x_ts = np.asarray(x_ts, dtype=float).reshape(-1)
        if self.x_ts.shape != self.sqm.shape or not np.all(np.isfinite(self.x_ts)):
            raise ValueError("IRC coordinates must be finite and match 3*natom")
        d = np.asarray(init_direction_mw, dtype=float).reshape(-1)
        dnorm = np.linalg.norm(d)
        if d.shape != self.sqm.shape or not np.all(np.isfinite(d)) or dnorm <= 0.0:
            raise ValueError("IRC initial direction must be a finite nonzero 3*natom vector")
        self.dir = d / dnorm
        self.step = float(step)
        if not np.isfinite(self.step) or self.step <= 0.0:
            raise ValueError("IRC step must be finite and positive")
        self.sign = 1 if sign >= 0 else -1
        self.logger = logger

    # mass-weighting helpers
    def _to_mw(self, x):
        return x * self.sqm

    def _to_cart(self, q):
        return q / self.sqm

    def _mw_grad(self, gx):
        return gx / self.sqm

    def run(self, energy_gradient, max_points=50, gtol=1.0e-4,
            micro_max=6, micro_tol=1.0e-4, on_point=None):
        if int(max_points) < 1:
            raise ValueError("IRC max_points must be at least 1")
        for name, value in {"gtol": gtol, "micro_tol": micro_tol}.items():
            value = float(value)
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError("IRC %s must be finite and positive" % name)
        if int(micro_max) < 1:
            raise ValueError("IRC micro_max must be at least 1")
        s = self.step
        path = []
        prev_e = None

        # First point: displace from the TS along the imaginary mode.
        q = self._to_mw(self.x_ts) + self.sign * s * self.dir

        for k in range(max_points):
            x = self._to_cart(q)
            e, gx = energy_gradient(x)
            g = self._mw_grad(np.asarray(gx, dtype=float).reshape(-1))
            gnorm = float(np.linalg.norm(g))
            gmax_cart = float(np.max(np.abs(np.asarray(gx).reshape(-1))))
            if gnorm < gtol:
                path.append({"x": x.copy(), "energy": float(e), "gmax": gmax_cart})
                if on_point is not None:
                    on_point(len(path) - 1, float(e), gmax_cart)
                return {"converged": True, "points": path, "npoint": len(path),
                        "reason": "gradient"}
            # The path has entered a minimum basin once the energy stops
            # decreasing (a fixed-arc-length IRC otherwise oscillates about it).
            # Do not retain the uphill overshoot as a false "minimum". Restore
            # the lower, previously accepted point so the returned path and the
            # caller's molecule/wavefunction remain synchronized there.
            if prev_e is not None and e > prev_e + 1.0e-9 and k >= 2:
                previous = path[-1]
                restored_e, restored_gx = energy_gradient(previous["x"])
                previous["energy"] = float(restored_e)
                previous["gmax"] = float(np.max(np.abs(
                    np.asarray(restored_gx, dtype=float).reshape(-1)
                )))
                return {"converged": True, "points": path,
                        "npoint": len(path), "reason": "minimum_bracket"}
            path.append({"x": x.copy(), "energy": float(e), "gmax": gmax_cart})
            if on_point is not None:
                on_point(len(path) - 1, float(e), gmax_cart)
            prev_e = e

            # Do not launch unrecorded micro-iterations after the final allowed
            # macro point.  This keeps the callback/molecule geometry identical
            # to the last point returned to the caller.
            if k == max_points - 1:
                break

            ghat = g / gnorm
            qstar = q - 0.5 * s * ghat
            # initial guess: full step along the gradient.
            qnew = q - s * ghat
            for _ in range(micro_max):
                xnew = self._to_cart(qnew)
                _, gxn = energy_gradient(xnew)
                gn = self._mw_grad(np.asarray(gxn, dtype=float).reshape(-1))
                rvec = qnew - qstar
                rnorm = np.linalg.norm(rvec)
                rhat = rvec / rnorm if rnorm > 1.0e-12 else rvec
                g_tan = gn - np.dot(gn, rhat) * rhat
                tnorm = np.linalg.norm(g_tan)
                if tnorm < micro_tol:
                    break
                # move ~s/2 tangentially, then project back onto the sphere.
                qnew = qnew - (0.5 * s / tnorm) * g_tan
                rvec = qnew - qstar
                qnew = qstar + 0.5 * s * rvec / np.linalg.norm(rvec)
            q = qnew

        return {"converged": False, "points": path, "npoint": max_points}


def _vibrational_basis(mass, coordinates):
    """Return an orthonormal mass-weighted vibrational basis.

    Translation and rotation vectors are constructed at the supplied geometry
    and removed with an SVD.  Its numerical rank naturally gives three
    external modes for an atom, five for a linear molecule, and six otherwise.
    """
    mass = np.asarray(mass, dtype=float).reshape(-1)
    coordinates = np.asarray(coordinates, dtype=float).reshape(-1, 3)
    if coordinates.shape != (mass.size, 3):
        raise ValueError("IRC coordinates must have shape (natom, 3)")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("IRC coordinates must contain only finite values")

    centered = coordinates - np.average(coordinates, axis=0, weights=mass)
    sqm = np.sqrt(mass)
    external = np.zeros((3 * mass.size, 6), dtype=float)
    for atom, (weight, position) in enumerate(zip(sqm, centered)):
        external[3 * atom:3 * atom + 3, :3] = np.eye(3) * weight
        for axis in range(3):
            unit = np.zeros(3)
            unit[axis] = 1.0
            external[3 * atom:3 * atom + 3, 3 + axis] = (
                weight * np.cross(unit, position)
            )

    left, singular, _ = np.linalg.svd(external, full_matrices=True)
    scale = float(singular[0]) if singular.size else 1.0
    rank = int(np.count_nonzero(singular > 1.0e-10 * max(scale, 1.0)))
    return left[:, rank:]


def imaginary_mode(mass, hessian, coordinates=None, negative_tol=1.0e-8):
    """Return the unique imaginary mode of a first-order saddle Hessian.

    ``mass`` in amu (natom,), ``hessian`` in Hartree/Bohr**2 (3N,3N), and
    optional ``coordinates`` in Bohr.  With coordinates, translations and
    rotations are projected before diagonalization so their near-zero
    numerical curvature cannot be mistaken for imaginary vibrations. Returns
    the normalized mass-weighted eigenvector and curvature. Numerical negative
    noise above ``-negative_tol`` is ignored. A valid IRC start must contain
    exactly one significant vibrational negative mode; minima and higher-order
    saddles are rejected.
    """
    mass = np.asarray(mass, dtype=float).reshape(-1)
    if mass.size == 0 or not np.all(np.isfinite(mass)) or np.any(mass <= 0.0):
        raise ValueError("IRC masses must be finite and positive")
    mr = np.repeat(mass, 3)
    hessian = np.asarray(hessian, dtype=float)
    expected = (mr.size, mr.size)
    if hessian.shape != expected or not np.all(np.isfinite(hessian)):
        raise ValueError("IRC Hessian must be a finite matrix with shape %s" % (expected,))
    mw_hess = hessian / np.sqrt(np.outer(mr, mr))
    mw_hess = 0.5 * (mw_hess + mw_hess.T)
    if coordinates is None:
        w, v = np.linalg.eigh(mw_hess)
    else:
        basis = _vibrational_basis(mass, coordinates)
        if basis.shape[1] == 0:
            raise ValueError("IRC requires at least one vibrational degree of freedom")
        w, internal_modes = np.linalg.eigh(basis.T @ mw_hess @ basis)
        v = basis @ internal_modes
    negative = np.flatnonzero(w < -abs(float(negative_tol)))
    if len(negative) != 1:
        raise ValueError(
            "IRC requires exactly one significant negative Hessian mode; found %d"
            % len(negative)
        )
    mode_index = int(negative[0])
    return v[:, mode_index], float(w[mode_index])
