"""Reference/regression test for the nuclear-repulsion Hessian contract.

The native ``grd1::hess_nn`` kernel accumulates the second derivatives of the
nuclear-repulsion energy into the OpenQP ``(3N, 3N)`` Cartesian Hessian
(atom-major order ``x1, y1, z1, x2, ...``). This test pins the closed-form
expression the Fortran routine implements by checking it against a central
finite difference of the nuclear-repulsion gradient. It is dependency-light
(NumPy only) so it runs in CI without the compiled native library.

The Fortran link test ``tests/fortran/test_hess_nn.F90`` performs the same
check against the production routines.
"""

import unittest

import numpy as np


def nuclear_repulsion_energy(xyz, zeff):
    """E_nn = sum_{k>l} Zk*Zl / r_kl, with xyz shape (natom, 3)."""
    e = 0.0
    nat = len(zeff)
    for k in range(nat):
        for l in range(k):
            r = np.linalg.norm(xyz[k] - xyz[l])
            e += zeff[k] * zeff[l] / r
    return e


def nuclear_repulsion_gradient(xyz, zeff):
    """Analytic gradient, shape (natom, 3); mirror of grd1::grad_nn."""
    nat = len(zeff)
    grad = np.zeros((nat, 3))
    for k in range(nat):
        for l in range(nat):
            if k == l:
                continue
            p = xyz[k] - xyz[l]
            r3 = np.linalg.norm(p) ** 3
            grad[k] += -zeff[k] * zeff[l] * p / r3
    return grad


def nuclear_repulsion_hessian(xyz, zeff):
    """Analytic (3N, 3N) Hessian; mirror of grd1::hess_nn.

    Pair (k,l) 3x3 block:  Zk*Zl * (3 p_a p_b / r^5 - delta_ab / r^3).
    """
    nat = len(zeff)
    hess = np.zeros((3 * nat, 3 * nat))
    eye = np.eye(3)
    for k in range(1, nat):
        for l in range(k):
            p = xyz[k] - xyz[l]
            r = np.linalg.norm(p)
            zz = zeff[k] * zeff[l]
            blk = zz * (3.0 * np.outer(p, p) / r ** 5 - eye / r ** 3)
            ks, le = slice(3 * k, 3 * k + 3), slice(3 * l, 3 * l + 3)
            hess[ks, ks] += blk
            hess[le, le] += blk
            hess[ks, le] -= blk
            hess[le, ks] -= blk
    return hess


class NuclearRepulsionHessianTest(unittest.TestCase):
    def setUp(self):
        # arbitrary water-like geometry in bohr
        self.xyz = np.array(
            [
                [0.0, 0.0, 0.20],
                [0.0, 1.45, -0.95],
                [0.0, -1.45, -0.95],
            ]
        )
        self.zeff = np.array([8.0, 1.0, 1.0])

    def test_gradient_matches_finite_difference(self):
        analytic = nuclear_repulsion_gradient(self.xyz, self.zeff).reshape(-1)
        fd = self._fd_gradient()
        self.assertLess(np.max(np.abs(analytic - fd)), 1e-7)

    def test_hessian_matches_finite_difference(self):
        analytic = nuclear_repulsion_hessian(self.xyz, self.zeff)
        fd = self._fd_hessian()
        self.assertLess(np.max(np.abs(analytic - fd)), 1e-6)

    def test_hessian_is_symmetric(self):
        h = nuclear_repulsion_hessian(self.xyz, self.zeff)
        self.assertLess(np.max(np.abs(h - h.T)), 1e-12)

    def test_translational_invariance(self):
        # Sum of force-constant blocks over one atom index must vanish
        # (rigid translation is a zero-frequency mode of E_nn).
        h = nuclear_repulsion_hessian(self.xyz, self.zeff)
        nat = len(self.zeff)
        acc = np.zeros((3, 3))
        for k in range(nat):
            for l in range(nat):
                acc += h[3 * k:3 * k + 3, 3 * l:3 * l + 3]
        self.assertLess(np.max(np.abs(acc)), 1e-10)

    def _fd_gradient(self, h=1e-6):
        nat = len(self.zeff)
        g = np.zeros(3 * nat)
        x = self.xyz.copy()
        for i in range(3 * nat):
            k, a = divmod(i, 3)
            x[k, a] += h
            ep = nuclear_repulsion_energy(x, self.zeff)
            x[k, a] -= 2 * h
            em = nuclear_repulsion_energy(x, self.zeff)
            x[k, a] += h
            g[i] = (ep - em) / (2 * h)
        return g

    def _fd_hessian(self, h=1e-5):
        nat = len(self.zeff)
        n = 3 * nat
        hess = np.zeros((n, n))
        x = self.xyz.copy()
        for i in range(n):
            k, a = divmod(i, 3)
            x[k, a] += h
            gp = nuclear_repulsion_gradient(x, self.zeff).reshape(-1)
            x[k, a] -= 2 * h
            gm = nuclear_repulsion_gradient(x, self.zeff).reshape(-1)
            x[k, a] += h
            hess[:, i] = (gp - gm) / (2 * h)
        return hess


if __name__ == "__main__":
    unittest.main()
