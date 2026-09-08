"""Ewald electrostatics for periodic full-ESPF QM/MM (orthorhombic cells).

Implements the lattice-summed electrostatics of Bonfrate, Ferre and
Huix-Rotllant, JCTC 2024, 20, 4338 (eqs 7-8 of the paper; standard Ewald
splitting instead of PME) in atomic units (bohr, e, Hartree), with tin-foil
boundary conditions and a uniform neutralising background for non-neutral
charge sets (a position-independent constant, so no force):

* ``mm_potential(r_qm)``      Phi^MM at the QM(+link) centres from all MM charges
                              and their images, plus its gradient with respect
                              to the QM centre (for the QM coupling force);
* ``mm_forces(r_qm, q_qm)``   force on every MM atom from the QM charges q_qm
                              interacting with all MM images;
* ``qm_image_matrix(r_qm)``   psi_img(A,B): Ewald pair potential between QM
                              centre A and the *images* of QM centre B (the
                              in-cell 1/r is removed; the diagonal is the
                              self-image term), and its gradient.

The pair potential is  psi(r) = sum_n erfc(beta|r+nL|)/|r+nL|
                              + (4 pi/V) sum_{k/=0} exp(-k^2/4beta^2)/k^2 cos(k.r)
                              - pi/(beta^2 V),
with beta chosen so the real-space sum is converged at half the shortest box
edge (minimum image only) and the reciprocal sum truncated at a relative
accuracy ``tol``.
"""

import numpy as np
from scipy.special import erfc


class EwaldQMMM:
    def __init__(self, box_bohr, tol=1e-9):
        box = np.asarray(box_bohr, dtype=float).reshape(3)
        if np.any(box <= 0.0):
            raise ValueError("EwaldQMMM needs an orthorhombic box with positive edges")
        self.box = box
        self.volume = float(np.prod(box))
        self.tol = float(tol)
        self.rc = 0.5 * float(box.min())
        # erfc(beta*rc) ~ tol  ->  beta = x/rc with erfc(x)=tol
        x = 4.5 if tol <= 1e-9 else 4.0
        self.beta = x / self.rc
        # exp(-k^2/(4 beta^2)) < tol  ->  k_max = 2 beta sqrt(-ln tol)
        kmax = 2.0 * self.beta * np.sqrt(-np.log(tol))
        mmax = np.ceil(kmax * box / (2.0 * np.pi)).astype(int)
        grids = [np.arange(-m, m + 1) for m in mmax]
        m = np.array(np.meshgrid(*grids, indexing="ij")).reshape(3, -1).T
        m = m[np.any(m != 0, axis=1)]
        k = 2.0 * np.pi * m / box                       # (nk, 3)
        k2 = np.einsum("ij,ij->i", k, k)
        keep = k2 <= kmax * kmax
        self.k = k[keep]
        k2 = k2[keep]
        self.w = (4.0 * np.pi / self.volume) * np.exp(-k2 / (4.0 * self.beta ** 2)) / k2
        self.self_term = float(self.w.sum()) - 2.0 * self.beta / np.sqrt(np.pi)
        self.background = -np.pi / (self.beta ** 2 * self.volume)

    # ------------------------------------------------------------------ #
    def _min_image(self, d):
        return d - self.box * np.round(d / self.box)

    def max_damping_width(self):
        """Largest Gaussian MM-charge width (bohr) for which the erf-damping
        correction -erfc(mu r)/r, evaluated for the minimum image inside the
        real-space cutoff only, is converged to the lattice-sum tolerance:
        erfc(mu rc) <= tol  with  mu = 1/(sqrt(2) w)."""
        x = 4.5 if self.tol <= 1e-9 else 4.0        # erfc(x) ~ tol, as for beta
        return self.rc / (np.sqrt(2.0) * x)

    def check_damping(self, mu):
        """Raise unless the damping parameter ``mu`` (1/bohr) is short-ranged
        on the real-space cutoff scale (see max_damping_width)."""
        if mu is None:
            return
        mu = float(mu)
        if not np.isfinite(mu) or mu <= 0.0:
            raise ValueError(f"damping parameter mu must be finite and positive; got {mu!r}")
        w = 1.0 / (np.sqrt(2.0) * mu)
        w_max = self.max_damping_width()
        if w > w_max:
            raise ValueError(
                f"[qmmm] mm_charge_width = {w / 1.8897259886:.3f} A is too large "
                f"for this box: the erf-damping correction is summed for the "
                f"minimum image inside the real-space cutoff rc = "
                f"{self.rc / 1.8897259886:.2f} A only, which needs a width <= "
                f"{w_max / 1.8897259886:.3f} A (erfc(mu rc) <= {self.tol:.0e}).")

    def _real_pair(self, d, mu=None):
        """erfc(beta r)/r and its radial derivative for displacement(s) d.
        With ``mu`` the MM charge is a Gaussian of width 1/(sqrt(2) mu) instead
        of a point: the pair potential 1/r becomes erf(mu r)/r, i.e. the
        short-range term -erfc(mu r)/r is added to the Ewald real-space part."""
        r = np.linalg.norm(d, axis=-1)
        r = np.where(r < 1e-12, 1e-12, r)
        e = erfc(self.beta * r)
        g = np.exp(-(self.beta * r) ** 2)
        psi = e / r
        dpsi = -(e / r ** 2 + 2.0 * self.beta / np.sqrt(np.pi) * g / r)   # d psi / d r
        if mu is not None:
            em = erfc(mu * r)
            gm = np.exp(-(mu * r) ** 2)
            psi = psi - em / r
            dpsi = dpsi + (em / r ** 2 + 2.0 * mu / np.sqrt(np.pi) * gm / r)
        return r, psi, dpsi

    #: atoms per block in the reciprocal-space sums; bounds the (n, nk) work
    #: arrays to a few hundred MB for solvated proteins (a 16 A water box fits
    #: in one block, so small systems are summed exactly as before)
    CHUNK = 2048

    def _structure_factor(self, pos, q):
        s = np.zeros(len(self.k), dtype=complex)
        for i0 in range(0, len(pos), self.CHUNK):
            phase = pos[i0:i0 + self.CHUNK] @ self.k.T                 # (n, nk)
            s += (q[i0:i0 + self.CHUNK, None] * np.exp(1j * phase)).sum(axis=0)
        return s                                                       # (nk,)

    # ------------------------------------------------------------------ #
    def mm_potential(self, r_qm, r_mm, q_mm, mu=None):
        """Phi^MM at each QM centre (Hartree/e) and d Phi^MM / d r_qm (Hartree/e/bohr)."""
        r_qm = np.atleast_2d(np.asarray(r_qm, dtype=float))
        r_mm = np.atleast_2d(np.asarray(r_mm, dtype=float))
        q_mm = np.asarray(q_mm, dtype=float)
        phi = np.zeros(len(r_qm))
        grad = np.zeros((len(r_qm), 3))
        # real space (minimum image, converged within rc by construction of beta)
        for a in range(len(r_qm)):
            d = self._min_image(r_qm[a] - r_mm)
            r, psi, dpsi = self._real_pair(d, mu)
            mask = r < self.rc
            phi[a] = np.sum(q_mm[mask] * psi[mask])
            grad[a] = np.sum((q_mm[mask] * dpsi[mask] / r[mask])[:, None] * d[mask], axis=0)
        # reciprocal space
        s_mm = self._structure_factor(r_mm, q_mm)                 # (nk,)
        phase = r_qm @ self.k.T                                    # (na, nk)
        z = s_mm[None, :] * np.exp(-1j * phase)                    # S e^{-ik.r_A}
        phi += (self.w[None, :] * z.real).sum(axis=1)
        grad += (self.w[None, :] * z.imag) @ self.k
        # neutralising background (constant)
        phi += self.background * q_mm.sum()
        return phi, grad

    def mm_forces(self, r_qm, q_qm, r_mm, q_mm, mu=None):
        """Force (Hartree/bohr) on every MM atom from the interaction
        E = sum_A q_A Phi^MM_A with the QM charges q_qm."""
        r_qm = np.atleast_2d(np.asarray(r_qm, dtype=float))
        r_mm = np.atleast_2d(np.asarray(r_mm, dtype=float))
        q_qm = np.asarray(q_qm, dtype=float)
        q_mm = np.asarray(q_mm, dtype=float)
        f = np.zeros((len(r_mm), 3))
        for a in range(len(r_qm)):
            d = self._min_image(r_qm[a] - r_mm)                    # r_A - r_i
            r, psi, dpsi = self._real_pair(d, mu)
            mask = r < self.rc
            # dE/dr_i = q_A q_i dpsi * (r_i - r_A)/r = -q_A q_i dpsi d/r ; F = -dE/dr_i
            f[mask] += (q_qm[a] * q_mm[mask] * dpsi[mask] / r[mask])[:, None] * d[mask]
        s_qm = self._structure_factor(r_qm, q_qm)
        for i0 in range(0, len(r_mm), self.CHUNK):
            sl = slice(i0, i0 + self.CHUNK)
            phase = r_mm[sl] @ self.k.T
            z = s_qm[None, :] * np.exp(-1j * phase)                # S_QM e^{-ik.r_i}
            # dE/dr_i = q_i sum_k w_k k Im[S_QM e^{-ik r_i}]  -> F = -that
            f[sl] -= q_mm[sl, None] * ((self.w[None, :] * z.imag) @ self.k)
        return f

    def qm_image_matrix(self, r_qm):
        """psi_img (na, na) and its gradient d psi_img(A,B)/d r_A (na, na, 3)."""
        r_qm = np.atleast_2d(np.asarray(r_qm, dtype=float))
        na = len(r_qm)
        psi = np.zeros((na, na))
        dpsi = np.zeros((na, na, 3))
        for a in range(na):
            for b in range(na):
                if a == b:
                    psi[a, a] = self.self_term + self.background
                    continue
                d = self._min_image(r_qm[a] - r_qm[b])
                d0 = r_qm[a] - r_qm[b]                              # in-cell direct pair
                r, p, dp = self._real_pair(d)
                r0 = np.linalg.norm(d0)
                kr = self.k @ d
                p_rec = float(np.sum(self.w * np.cos(kr)))
                g_rec = -(self.w * np.sin(kr)) @ self.k
                psi[a, b] = p + p_rec - 1.0 / r0 + self.background
                dpsi[a, b] = dp * d / r + g_rec + d0 / r0 ** 3
        return psi, dpsi

    def qm_image_energy_force(self, r_qm, q_qm):
        """E_img = 1/2 sum_AB q_A psi_img(A,B) q_B and the force on each QM centre
        at fixed charges."""
        q = np.asarray(q_qm, dtype=float)
        psi, dpsi = self.qm_image_matrix(r_qm)
        e = 0.5 * float(q @ psi @ q)
        f = -np.einsum("a,b,abc->ac", q, q, dpsi)
        return e, f, psi
