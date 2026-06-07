"""OQP nudged elastic band (NEB) for OpenQP.

Backend-agnostic (NumPy only) chain-of-states reaction-path optimizer.  The band
of images is optimized in Cartesian coordinates with:

* the improved tangent estimate (Henkelman & Jonsson, J. Chem. Phys. 113, 9978
  (2000));
* a parallel spring force on the segment-length difference and the perpendicular
  component of the true force;
* an optional climbing image (Henkelman, Uberuaga, Jonsson, J. Chem. Phys. 113,
  9901 (2000)) that drives the highest-energy image to the saddle; and
* the FIRE optimizer (Bitzek et al., Phys. Rev. Lett. 97, 170201 (2006)), which
  is robust for the non-conservative NEB force where a quasi-Newton RFO step is
  not appropriate.

The endpoints (reactant/product) are held fixed.  The caller supplies an
``energy_gradient(x_flat) -> (E, g_flat)`` callback (atomic units), so this
module is testable on analytic potentials without the compiled OQP core.
"""

from __future__ import annotations

import numpy as np


class NEB:
    def __init__(self, images, k_spring=0.05, climbing=False,
                 climb_fmax=0.05, logger=None):
        # images: list of flat Cartesian coordinate arrays (Bohr); endpoints
        # at index 0 and -1 are held fixed.
        self.images = [np.asarray(x, dtype=float).reshape(-1) for x in images]
        self.nimage = len(self.images)
        if self.nimage < 3:
            raise ValueError("NEB needs at least 3 images (2 endpoints + 1)")
        self.k = float(k_spring)
        self.climbing = bool(climbing)
        # Relax-then-climb: the climbing image is only switched on once the
        # regular-NEB perpendicular force first drops below this threshold, so
        # it does not destabilize a still-tangled band.
        self.climb_fmax = float(climb_fmax)
        self._climb_active = False
        self.logger = logger
        self.energies = [None] * self.nimage
        self.gradients = [None] * self.nimage

    # -- forces ------------------------------------------------------------- #
    @staticmethod
    def _tangent(xm, xi, xp, em, ei, ep):
        """Improved tangent (Henkelman-Jonsson 2000)."""
        tp = xp - xi
        tm = xi - xm
        if ep > ei > em:
            tau = tp
        elif ep < ei < em:
            tau = tm
        else:
            d1 = abs(ep - ei)
            d2 = abs(em - ei)
            dmax, dmin = max(d1, d2), min(d1, d2)
            if ep > em:
                tau = tp * dmax + tm * dmin
            else:
                tau = tp * dmin + tm * dmax
        n = np.linalg.norm(tau)
        return tau / n if n > 1.0e-12 else tau

    def _band_forces(self, climb_on):
        """Return the NEB force on each interior image and the max |F_perp|."""
        forces = [np.zeros_like(x) for x in self.images]
        fmax = 0.0
        # index of the climbing image: highest energy interior image.
        climb_idx = -1
        if climb_on:
            interior_e = [self.energies[i] for i in range(1, self.nimage - 1)]
            climb_idx = 1 + int(np.argmax(interior_e))

        for i in range(1, self.nimage - 1):
            xm, xi, xp = self.images[i - 1], self.images[i], self.images[i + 1]
            em, ei, ep = self.energies[i - 1], self.energies[i], self.energies[i + 1]
            tau = self._tangent(xm, xi, xp, em, ei, ep)
            g = self.gradients[i]
            g_par = np.dot(g, tau) * tau
            g_perp = g - g_par

            if climb_on and i == climb_idx:
                # Climbing image: invert the parallel force, no spring.  Its
                # convergence must be judged by the full climbing force, which
                # vanishes only at the saddle -- using g_perp alone would report
                # convergence while a large parallel gradient still pushes the
                # image off the saddle (e.g. a symmetric path with g_perp == 0).
                forces[i] = -g + 2.0 * g_par
                fmax = max(fmax, float(np.max(np.abs(forces[i]))))
            else:
                seg_p = np.linalg.norm(xp - xi)
                seg_m = np.linalg.norm(xi - xm)
                f_spring = self.k * (seg_p - seg_m) * tau
                forces[i] = -g_perp + f_spring
                fmax = max(fmax, float(np.max(np.abs(g_perp))))
        return forces, fmax

    # -- driver ------------------------------------------------------------- #
    def run(self, energy_gradient, fmax_tol=2.0e-3, maxiter=200,
            dt=0.5, dt_max=1.0, maxmove=0.2, on_iteration=None):
        # Evaluate the fixed endpoints once.
        for i in (0, self.nimage - 1):
            self.energies[i], self.gradients[i] = energy_gradient(self.images[i])

        # FIRE state over the concatenated interior images.
        interior = list(range(1, self.nimage - 1))
        v = {i: np.zeros_like(self.images[i]) for i in interior}
        alpha0, f_alpha, f_inc, f_dec, n_min = 0.1, 0.99, 1.1, 0.5, 5
        dt0 = dt
        alpha = alpha0
        n_pos = 0

        for it in range(maxiter):
            for i in interior:
                self.energies[i], self.gradients[i] = energy_gradient(self.images[i])

            forces, fmax = self._band_forces(self._climb_active)
            # Latch the climbing image on once the band is sufficiently relaxed,
            # and reset the FIRE state so relaxation-phase momentum does not
            # overshoot the saddle when the parallel force flips sign.
            if self.climbing and not self._climb_active and fmax < self.climb_fmax:
                self._climb_active = True
                for i in interior:
                    v[i] = np.zeros_like(v[i])
                dt = dt0
                alpha = alpha0
                n_pos = 0
                forces, fmax = self._band_forces(self._climb_active)

            if on_iteration is not None:
                on_iteration(it, fmax, list(self.energies))
            if fmax < fmax_tol:
                return {"converged": True, "iters": it + 1, "fmax": fmax,
                        "energies": list(self.energies)}

            # FIRE update (global power over the whole band).
            power = sum(float(np.dot(forces[i], v[i])) for i in interior)
            if power > 0.0:
                n_pos += 1
                vnorm = np.sqrt(sum(float(v[i] @ v[i]) for i in interior))
                fnorm = np.sqrt(sum(float(forces[i] @ forces[i]) for i in interior))
                scale = (alpha * vnorm / fnorm) if fnorm > 1.0e-12 else 0.0
                for i in interior:
                    v[i] = (1.0 - alpha) * v[i] + scale * forces[i]
                if n_pos > n_min:
                    dt = min(dt * f_inc, dt_max)
                    alpha *= f_alpha
            else:
                n_pos = 0
                dt *= f_dec
                alpha = alpha0
                for i in interior:
                    v[i] = np.zeros_like(v[i])

            for i in interior:
                v[i] = v[i] + dt * forces[i]
                step = dt * v[i]
                # Cap the per-image displacement so a large climbing-image force
                # cannot fling an image into a geometry the SCF backend chokes on.
                snorm = np.linalg.norm(step)
                if snorm > maxmove:
                    step = step * (maxmove / snorm)
                self.images[i] = self.images[i] + step

        return {"converged": False, "iters": maxiter, "fmax": fmax,
                "energies": list(self.energies)}
