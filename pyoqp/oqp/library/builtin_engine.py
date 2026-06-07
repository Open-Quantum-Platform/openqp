"""Step-determination engine for the builtin OpenQP optimizer.

Backend-agnostic (NumPy only).  The :class:`BuiltinEngine` owns the optimization
loop and drives an energy/gradient callback supplied by the caller, so it can be
exercised on analytic test potentials without the compiled OQP core.

Algorithm (v1):

* working coordinates from :mod:`oqp.library.builtin_coords` (redundant internals
  with Cartesian fall-back);
* restricted-step Rational Function Optimization (RFO) for minima and
  partitioned-RFO (P-RFO, eigenvector following) for transition states;
* Schlegel-type diagonal model Hessian updated with BFGS (minima) or the
  Bofill SR1/PSB mixture (TS);
* a predictive trust region in Cartesian RMSD, grown/shrunk on the ratio of
  actual to predicted energy change.

References: Banerjee, Adams, Simons, Shepard, J. Phys. Chem. 89, 52 (1985);
Bofill, J. Comput. Chem. 15, 1 (1994); Peng, Ayala, Schlegel, Frisch,
J. Comput. Chem. 17, 49 (1996).
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla

from oqp.library.builtin_coords import build_coordinates, CartesianCoordinates

# All dense linear algebra in this module is routed through LAPACK
# (scipy.linalg / numpy.linalg) and BLAS-3 GEMM (the ``@`` operators on the
# B-matrix and Hessian), so it inherits the threaded, vendor-tuned backends
# OpenQP already links against.  The matrices here are O(3N), so the wall-clock
# is dominated by the electronic-structure energy/gradient, not these solves.


def _sym_eigh(a):
    """Symmetric eigensolve via LAPACK (syevd driver)."""
    return sla.eigh(0.5 * (a + a.T))


class ConvergenceSignal(Exception):
    """Raised by the convergence callback to stop the loop cleanly."""


class BuiltinEngine:
    def __init__(self, atoms, x0, mode="min", trust=0.2, trust_min=5.0e-3,
                 trust_max=0.5, maxiter=100, follow_mode=0, coordsys="auto",
                 logger=None):
        self.atoms = np.asarray(atoms, dtype=int).reshape(-1)
        self.x = np.asarray(x0, dtype=float).reshape(-1)
        self.mode = mode
        self.trust = float(trust)
        self.trust_min = float(trust_min)
        self.trust_max = float(trust_max)
        self.maxiter = int(maxiter)
        self.follow_mode = int(follow_mode)
        self.logger = logger

        self.coords = build_coordinates(self.atoms, self.x, coordsys=coordsys)
        self.coordsys = type(self.coords).__name__
        self.H = self.coords.guess_hessian(self.x)
        self._prev = None

    # -- public driver ------------------------------------------------------ #
    def run(self, energy_gradient, on_converged=None):
        """Optimize until ``on_converged`` raises or ``maxiter`` is reached.

        ``energy_gradient(x_flat) -> (E, g_flat)`` returns energy (Hartree) and
        Cartesian gradient (Hartree/Bohr).  ``on_converged()`` is called after
        every evaluation and may raise to stop (mirrors OQP's ``StopIteration``
        convergence protocol).
        """
        try:
            for _ in range(self.maxiter):
                e, g_cart = energy_gradient(self.x)
                if on_converged is not None:
                    on_converged()
                self.x = self._take_step(e, np.asarray(g_cart, dtype=float).reshape(-1))
        except (StopIteration, ConvergenceSignal):
            pass
        return self.x

    # -- one macro-iteration ------------------------------------------------ #
    def _take_step(self, e, g_cart):
        b = self.coords.b_matrix(self.x)
        g_q = self.coords.grad_to_q(self.x, g_cart)

        if self._prev is not None:
            s = self._prev["dq"]
            y = g_q - self._prev["g_q"]
            self.H = self._update_hessian(self.H, s, y)
            actual = e - self._prev["e"]
            pred = self._prev["pred"]
            self._update_trust(actual, pred, self._prev["cart_step"])

        dq = self._rfo_step(self.H, g_q)
        dq = self._restrict_to_trust(b, dq)

        x_new, ok = self.coords.back_transform(self.x, dq)
        if not ok:
            # Shrink and retry once; otherwise fall back to a plain scaled step.
            self.trust = max(self.trust * 0.5, self.trust_min)
            dq = self._restrict_to_trust(b, dq)
            x_new, ok = self.coords.back_transform(self.x, dq)
            if not ok:
                x_new = self.x + (b.T @ dq) * 0.5

        dq_taken = self.coords.q_displacement(self.coords.q(x_new),
                                              self.coords.q(self.x))
        cart_step = self.coords.cart_rmsd(self.x, x_new)
        pred = float(g_q @ dq_taken + 0.5 * dq_taken @ self.H @ dq_taken)
        self._prev = {"e": e, "g_q": g_q, "dq": dq_taken,
                      "pred": pred, "cart_step": cart_step}
        return x_new

    # -- RFO / P-RFO -------------------------------------------------------- #
    def _rfo_step(self, h, g):
        if self.mode == "ts":
            return self._prfo_step(h, g)
        return self._simple_rfo_step(h, g)

    @staticmethod
    def _simple_rfo_step(h, g):
        """Minimization RFO: lowest-eigenvalue augmented-Hessian step."""
        n = len(g)
        aug = np.zeros((n + 1, n + 1))
        aug[:n, :n] = h
        aug[:n, n] = g
        aug[n, :n] = g
        w, v = _sym_eigh(aug)
        vec = v[:, 0]
        if abs(vec[-1]) < 1.0e-10:
            # Degenerate scaling: fall back to a damped steepest-descent step.
            return -g / (np.linalg.norm(g) + 1.0)
        return vec[:n] / vec[n]

    def _prfo_step(self, h, g):
        """Partitioned RFO for TS: maximize along one mode, minimize the rest."""
        w, v = _sym_eigh(h)
        f = v.T @ g  # gradient in the Hessian eigenbasis
        t = self._select_follow_mode(w, v)

        step = np.zeros_like(g)
        # Maximized mode: take the upper root of the 2x2 augmented problem.
        lam = w[t]
        nu_max = 0.5 * (lam + np.sqrt(lam * lam + 4.0 * f[t] * f[t]))
        denom = lam - nu_max
        if abs(denom) > 1.0e-10:
            step += (-f[t] / denom) * v[:, t]

        # Minimized subspace: RFO over the remaining modes.
        idx = [i for i in range(len(w)) if i != t]
        if idx:
            lam_min = w[idx]
            f_min = f[idx]
            m = len(idx)
            aug = np.zeros((m + 1, m + 1))
            aug[:m, :m] = np.diag(lam_min)
            aug[:m, m] = f_min
            aug[m, :m] = f_min
            wa, va = _sym_eigh(aug)
            nu_min = wa[0]
            for jj, i in enumerate(idx):
                d = lam_min[jj] - nu_min
                if abs(d) > 1.0e-10:
                    step += (-f_min[jj] / d) * v[:, i]
        self._followed = v[:, t].copy()
        return step

    def _select_follow_mode(self, w, v):
        """Pick the mode to maximize: lowest curvature, with mode following."""
        prev = getattr(self, "_followed", None)
        if prev is not None:
            overlaps = np.abs(v.T @ prev)
            return int(np.argmax(overlaps))
        order = np.argsort(w)
        return int(order[self.follow_mode])

    # -- trust region & step restriction ------------------------------------ #
    def _restrict_to_trust(self, b, dq):
        # Estimate the Cartesian motion of this internal step and scale so its
        # RMSD does not exceed the trust radius.
        ginv, _ = self.coords._g_inverse(b) if hasattr(self.coords, "_g_inverse") \
            else (np.eye(b.shape[0]), 0)
        dx = b.T @ (ginv @ dq) if hasattr(self.coords, "_g_inverse") else dq
        rmsd = float(np.sqrt(np.mean(dx ** 2)))
        if rmsd > self.trust and rmsd > 1.0e-14:
            dq = dq * (self.trust / rmsd)
        return dq

    def _update_trust(self, actual, pred, cart_step):
        if abs(pred) < 1.0e-14:
            return
        ratio = actual / pred
        if ratio < 0.25:
            self.trust = max(self.trust * 0.25, self.trust_min)
        elif ratio > 0.75 and cart_step >= 0.8 * self.trust:
            self.trust = min(self.trust * 2.0, self.trust_max)

    # -- Hessian updates ---------------------------------------------------- #
    def _update_hessian(self, h, s, y):
        if self.mode == "ts":
            return self._bofill_update(h, s, y)
        return self._bfgs_update(h, s, y)

    @staticmethod
    def _bfgs_update(h, s, y):
        sy = float(s @ y)
        hs = h @ s
        shs = float(s @ hs)
        if sy > 1.0e-12 and shs > 1.0e-12:
            h = h + np.outer(y, y) / sy - np.outer(hs, hs) / shs
        return h

    @staticmethod
    def _bofill_update(h, s, y):
        ss = float(s @ s)
        if ss < 1.0e-14:
            return h
        dg = y - h @ s            # gradient-difference residual
        dgs = float(dg @ s)
        dgdg = float(dg @ dg)
        if dgdg < 1.0e-14:
            return h
        phi = (dgs * dgs) / (ss * dgdg)
        # SR1 (Murtagh-Sargent) part.
        if abs(dgs) > 1.0e-12:
            h_sr1 = h + np.outer(dg, dg) / dgs
        else:
            h_sr1 = h
        # PSB part.
        h_psb = (h + (np.outer(dg, s) + np.outer(s, dg)) / ss
                 - dgs * np.outer(s, s) / (ss * ss))
        return phi * h_sr1 + (1.0 - phi) * h_psb
