"""Unit tests for the Ewald electrostatics of periodic full-ESPF QM/MM
(oqp.library.qmmm_ewald).  Pure numpy/scipy; no compiled backend needed.

* the Ewald potential of a neutral MM charge set at a QM point equals the
  explicit lattice sum over many images;
* every analytic derivative (potential gradient at the QM centre, force on the
  MM atoms, force from the QM-image energy) matches a central finite
  difference of the corresponding energy.
"""
import unittest

import numpy as np

try:
    from oqp.library.qmmm_ewald import EwaldQMMM
    _HAVE = True
except Exception:  # pragma: no cover
    _HAVE = False


@unittest.skipUnless(_HAVE, "oqp.library.qmmm_ewald not importable")
class TestEwaldQMMM(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(7)
        self.box = np.array([14.0, 16.0, 15.0])
        self.ew = EwaldQMMM(self.box, tol=1e-10)
        # neutral MM set with zero dipole about the box centre (charge q at
        # c+d and at c-d): a cubic block of images then converges to the
        # tin-foil Ewald value without a shape-dependent surface term.
        d = rng.uniform(-0.3, 0.3, (3, 3)) * self.box
        q = rng.normal(size=3); q -= q.mean()
        c = 0.5 * self.box
        self.r_mm = np.vstack([c + d, c - d])
        self.q_mm = np.concatenate([q, q])
        # 3 QM centres
        self.r_qm = np.array([[7.0, 8.0, 7.5], [8.4, 8.2, 7.1], [6.9, 9.3, 8.0]])
        self.q_qm = np.array([0.3, -0.5, 0.2])

    def test_damping_width_bound(self):
        ang = 1.8897259886
        ew = EwaldQMMM(np.array([16.0, 16.0, 16.0]) * ang)     # the example box
        w_max = ew.max_damping_width()
        # a 0.7 A Gaussian (the Tinker default) is fine in a 16 A box (bound
        # 1.26 A); the bound scales with the real-space cutoff (half the
        # shortest edge), so the 14-bohr test box of setUp rejects 0.7 A
        self.assertGreater(w_max, 0.7 * ang)
        self.assertLess(self.ew.max_damping_width(), 0.7 * ang)
        self.assertAlmostEqual(w_max, ew.rc / (np.sqrt(2.0) * 4.5), places=12)
        ew.check_damping(None)
        ew.check_damping(1.0 / (np.sqrt(2.0) * 0.7 * ang))
        with self.assertRaisesRegex(ValueError, 'mm_charge_width'):
            ew.check_damping(1.0 / (np.sqrt(2.0) * 1.05 * w_max))

    def test_chunked_reciprocal_sums_are_exact(self):
        """The block-wise reciprocal sums (CHUNK atoms at a time) reproduce the
        single-block result to round-off for a set larger than one block."""
        rng = np.random.default_rng(3)
        n = 700
        r_mm = rng.uniform(0, 1, (n, 3)) * self.box
        q_mm = rng.normal(size=n); q_mm -= q_mm.mean()
        ew = EwaldQMMM(self.box, tol=1e-8)
        ref_phi, ref_grad = ew.mm_potential(self.r_qm, r_mm, q_mm)
        ref_f = ew.mm_forces(self.r_qm, self.q_qm, r_mm, q_mm)
        ew.CHUNK = 128
        phi, grad = ew.mm_potential(self.r_qm, r_mm, q_mm)
        f = ew.mm_forces(self.r_qm, self.q_qm, r_mm, q_mm)
        np.testing.assert_allclose(phi, ref_phi, rtol=0, atol=1e-12)
        np.testing.assert_allclose(grad, ref_grad, rtol=0, atol=1e-12)
        np.testing.assert_allclose(f, ref_f, rtol=0, atol=1e-12)

    def test_potential_matches_explicit_lattice_sum(self):
        # The Ewald (tin-foil, zero-average) potential and a direct sum over a
        # block of images differ by a constant fixed by the second moment of the
        # charge distribution; potential *differences* between points agree.
        phi, _ = self.ew.mm_potential(self.r_qm, self.r_mm, self.q_mm)
        n = 12
        shifts = np.array(np.meshgrid(*[np.arange(-n, n + 1)] * 3, indexing="ij")).reshape(3, -1).T * self.box
        direct = []
        for a in range(len(self.r_qm)):
            d = self.r_qm[a] - self.r_mm[:, None, :] - shifts[None, :, :]
            r = np.linalg.norm(d, axis=-1)
            direct.append(float(np.sum(self.q_mm[:, None] / r)))
        direct = np.array(direct)
        for a in range(1, len(self.r_qm)):
            self.assertAlmostEqual(phi[a] - phi[0], direct[a] - direct[0], places=5)

    def test_qm_gradient_of_potential(self):
        h = 1e-4
        _, grad = self.ew.mm_potential(self.r_qm, self.r_mm, self.q_mm)
        for a in range(len(self.r_qm)):
            for c in range(3):
                rp = self.r_qm.copy(); rp[a, c] += h
                rm = self.r_qm.copy(); rm[a, c] -= h
                fd = (self.ew.mm_potential(rp, self.r_mm, self.q_mm)[0][a]
                      - self.ew.mm_potential(rm, self.r_mm, self.q_mm)[0][a]) / (2 * h)
                self.assertAlmostEqual(grad[a, c], fd, places=7)

    def test_mm_forces(self):
        h = 1e-4
        f = self.ew.mm_forces(self.r_qm, self.q_qm, self.r_mm, self.q_mm)

        def energy(r_mm):
            phi, _ = self.ew.mm_potential(self.r_qm, r_mm, self.q_mm)
            return float(self.q_qm @ phi)

        for i in range(len(self.r_mm)):
            for c in range(3):
                rp = self.r_mm.copy(); rp[i, c] += h
                rm = self.r_mm.copy(); rm[i, c] -= h
                fd = -(energy(rp) - energy(rm)) / (2 * h)
                self.assertAlmostEqual(f[i, c], fd, places=7)

    def test_damped_mm_potential_and_forces(self):
        """Gaussian-smeared MM charges (mu): gradient and MM forces stay consistent."""
        h = 1e-4; mu = 0.5
        phi, grad = self.ew.mm_potential(self.r_qm, self.r_mm, self.q_mm, mu=mu)
        phi0, _ = self.ew.mm_potential(self.r_qm, self.r_mm, self.q_mm)
        self.assertFalse(np.allclose(phi, phi0))                # damping changes the potential
        for a in range(len(self.r_qm)):
            for c in range(3):
                rp = self.r_qm.copy(); rp[a, c] += h
                rm = self.r_qm.copy(); rm[a, c] -= h
                fd = (self.ew.mm_potential(rp, self.r_mm, self.q_mm, mu=mu)[0][a]
                      - self.ew.mm_potential(rm, self.r_mm, self.q_mm, mu=mu)[0][a]) / (2 * h)
                self.assertAlmostEqual(grad[a, c], fd, places=7)
        f = self.ew.mm_forces(self.r_qm, self.q_qm, self.r_mm, self.q_mm, mu=mu)

        def energy(r_mm):
            return float(self.q_qm @ self.ew.mm_potential(self.r_qm, r_mm, self.q_mm, mu=mu)[0])

        for i in range(len(self.r_mm)):
            for c in range(3):
                rp = self.r_mm.copy(); rp[i, c] += h
                rm = self.r_mm.copy(); rm[i, c] -= h
                self.assertAlmostEqual(f[i, c], -(energy(rp) - energy(rm)) / (2 * h), places=7)

    def test_qm_image_force(self):
        h = 1e-4
        e0, f, psi = self.ew.qm_image_energy_force(self.r_qm, self.q_qm)
        self.assertTrue(np.allclose(psi, psi.T))
        for a in range(len(self.r_qm)):
            for c in range(3):
                rp = self.r_qm.copy(); rp[a, c] += h
                rm = self.r_qm.copy(); rm[a, c] -= h
                ep = self.ew.qm_image_energy_force(rp, self.q_qm)[0]
                em = self.ew.qm_image_energy_force(rm, self.q_qm)[0]
                self.assertAlmostEqual(f[a, c], -(ep - em) / (2 * h), places=7)


if __name__ == "__main__":
    unittest.main()
