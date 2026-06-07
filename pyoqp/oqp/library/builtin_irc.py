"""Builtin intrinsic reaction coordinate (IRC) for OpenQP.

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
        self.sqm = np.sqrt(np.repeat(mass, 3))          # M^{1/2}, length 3N
        self.x_ts = np.asarray(x_ts, dtype=float).reshape(-1)
        d = np.asarray(init_direction_mw, dtype=float).reshape(-1)
        self.dir = d / np.linalg.norm(d)
        self.step = float(step)
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
            path.append({"x": x.copy(), "energy": float(e), "gmax": gmax_cart})
            if on_point is not None:
                on_point(k, float(e), gmax_cart)
            if gnorm < gtol:
                return {"converged": True, "points": path, "npoint": k + 1,
                        "reason": "gradient"}
            # The path has entered a minimum basin once the energy stops
            # decreasing (a fixed-arc-length IRC otherwise oscillates about it).
            if prev_e is not None and e > prev_e + 1.0e-9 and k >= 2:
                return {"converged": True, "points": path, "npoint": k + 1,
                        "reason": "minimum"}
            prev_e = e

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


def imaginary_mode(mass, hessian):
    """Return the mass-weighted lowest-eigenvalue mode of a Cartesian Hessian.

    ``mass`` in amu (natom,), ``hessian`` in Hartree/Bohr**2 (3N,3N).  Returns
    the (normalized) eigenvector of the mass-weighted Hessian with the most
    negative eigenvalue and that eigenvalue (negative at a true TS).
    """
    mass = np.asarray(mass, dtype=float).reshape(-1)
    mr = np.repeat(mass, 3)
    mw_hess = np.asarray(hessian, dtype=float) / np.sqrt(np.outer(mr, mr))
    w, v = np.linalg.eigh(0.5 * (mw_hess + mw_hess.T))
    return v[:, 0], float(w[0])
