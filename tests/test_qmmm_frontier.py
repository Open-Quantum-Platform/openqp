"""Unit tests for QM/MM frontier (M1) charge redistribution.

These exercise the pure-Python RCD/RC/Z1 redistribution used by the QM/MM
driver to treat the electrostatics of a covalent QM/MM boundary: deleting the
MM host (M1) charge and redistributing it to virtual point charges at the M1-M2
bond midpoints.  They load the module by path so they run without the compiled
OpenQP backend or OpenMM.

The key property tested is gradient consistency: because every virtual charge
sits at a bond midpoint (a linear function of two real atom positions), the
Coulomb force on it must redistribute onto its hosts by the chain rule, so the
analytic coupling force matches a finite-difference of the coupling energy.
"""

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _coulomb_energy(qA, rA, qs, rs):
    """E = sum_A sum_s qA * qs / |rA - rs| (arbitrary consistent units)."""
    E = 0.0
    for a in range(len(qA)):
        d = rA[a] - rs
        r = np.linalg.norm(d, axis=1)
        E += qA[a] * np.sum(qs / r)
    return E


def _analytic_mm_forces(qA, rA, qs, rs, scatter, n_mm, row_of):
    """Force on each real MM atom, f = -dE/dR, via site forces + scatter."""
    f_site = np.zeros_like(rs)
    for a in range(len(qA)):
        d = rA[a] - rs
        r = np.linalg.norm(d, axis=1)
        coeff = qA[a] * qs / r ** 3
        f_site -= coeff[:, None] * d
    f_mm = np.zeros((n_mm, 3))
    for s, contribs in enumerate(scatter):
        for atom, w in contribs:
            f_mm[row_of[atom]] += w * f_site[s]
    return f_mm


class TestFrontierCharge(unittest.TestCase):
    def setUp(self):
        self.c = load_module(
            "qmmm_connectivity_frontier_under_test",
            "pyoqp/oqp/library/qmmm_connectivity.py",
        )

    @staticmethod
    def _charge_of(charges):
        return lambda a: charges[a]

    # -- charge conservation -------------------------------------------

    def test_rcd_conserves_total_charge(self):
        charges = {10: 0.5, 11: -0.3, 12: -0.2, 13: 0.11}   # M1=10, M2={11,12}
        deleted, dq, virt = self.c.redistribute_frontier_charges(
            [(10, [11, 12])], self._charge_of(charges), "rcd")
        orig = sum(charges.values())
        new = (sum(charges[a] for a in charges if a not in deleted)
               + sum(dq.values()) + sum(v.charge for v in virt))
        self.assertAlmostEqual(new, orig, places=12)

    def test_rc_conserves_total_charge(self):
        charges = {10: 0.5, 11: -0.3, 12: -0.2}
        deleted, dq, virt = self.c.redistribute_frontier_charges(
            [(10, [11, 12])], self._charge_of(charges), "rc")
        self.assertEqual(dq, {})   # RC does not touch M2 charges
        new = (sum(charges[a] for a in charges if a not in deleted)
               + sum(v.charge for v in virt))
        self.assertAlmostEqual(new, sum(charges.values()), places=12)

    def test_z1_deletes_only(self):
        charges = {10: 0.5, 11: -0.3, 12: -0.2}
        deleted, dq, virt = self.c.redistribute_frontier_charges(
            [(10, [11, 12])], self._charge_of(charges), "z1")
        self.assertEqual(deleted, {10})
        self.assertEqual(dq, {})
        self.assertEqual(virt, [])

    def test_none_is_noop(self):
        deleted, dq, virt = self.c.redistribute_frontier_charges(
            [(10, [11, 12])], self._charge_of({10: 1.0}), "none")
        self.assertEqual((deleted, dq, virt), (set(), {}, []))

    def test_isolated_host_falls_back_to_deletion(self):
        # M1 with no MM neighbours: cannot redistribute -> delete only.
        deleted, dq, virt = self.c.redistribute_frontier_charges(
            [(10, [])], self._charge_of({10: 0.4}), "rcd")
        self.assertEqual((deleted, dq, virt), ({10}, {}, []))

    def test_unknown_scheme_raises(self):
        with self.assertRaises(ValueError):
            self.c.redistribute_frontier_charges(
                [(10, [11])], self._charge_of({10: 1.0, 11: 0.0}), "bogus")

    # -- dipole conservation -------------------------------------------

    def _dipole_about_m1(self, pos, m1, dq, virt):
        mu = np.zeros(3)
        for v in virt:
            r = v.weights[0] * pos[v.hosts[0]] + v.weights[1] * pos[v.hosts[1]]
            mu += v.charge * (r - pos[m1])
        for a, dqa in dq.items():
            mu += dqa * (pos[a] - pos[m1])
        return mu

    def test_rcd_conserves_dipole_about_m1_but_rc_does_not(self):
        # M1 deleted (dipole of a point charge about itself is zero), so the
        # redistributed set must also have zero dipole about M1 for RCD.
        pos = {10: np.array([0.0, 0.0, 0.0]),
               11: np.array([1.5, 0.0, 0.0]),
               12: np.array([-0.4, 1.3, 0.2])}
        charges = {10: 0.7, 11: -0.2, 12: -0.2}
        _, dq, virt = self.c.redistribute_frontier_charges(
            [(10, [11, 12])], self._charge_of(charges), "rcd")
        mu_rcd = self._dipole_about_m1(pos, 10, dq, virt)
        self.assertTrue(np.allclose(mu_rcd, 0.0, atol=1e-12), mu_rcd)

        _, dq_rc, virt_rc = self.c.redistribute_frontier_charges(
            [(10, [11, 12])], self._charge_of(charges), "rc")
        mu_rc = self._dipole_about_m1(pos, 10, dq_rc, virt_rc)
        self.assertFalse(np.allclose(mu_rc, 0.0, atol=1e-9), mu_rc)

    def test_virtual_charges_sit_at_bond_midpoints(self):
        _, _, virt = self.c.redistribute_frontier_charges(
            [(10, [11, 12])], self._charge_of({10: 0.6, 11: 0.0, 12: 0.0}), "rcd")
        for v in virt:
            self.assertEqual(v.weights, (0.5, 0.5))
            self.assertEqual(v.hosts[0], 10)

    # -- site assembly --------------------------------------------------

    def test_assemble_identity_without_frontier(self):
        mm_idx = [3, 4, 5]
        q = [0.1, -0.2, 0.3]
        xyz = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
        qs, rs, scatter = self.c.assemble_embedding_sites(
            mm_idx, q, xyz, set(), {}, [])
        self.assertTrue(np.allclose(qs, q))
        self.assertTrue(np.allclose(rs, xyz))
        self.assertEqual([s for s in scatter], [((3, 1.0),), ((4, 1.0),), ((5, 1.0),)])

    def test_assemble_drops_deleted_and_adds_virtuals(self):
        mm_idx = [10, 11, 12]
        q = [0.6, -0.2, -0.2]
        xyz = np.array([[0, 0, 0], [1.5, 0, 0], [0, 1.5, 0]], float)
        deleted, dq, virt = self.c.redistribute_frontier_charges(
            [(10, [11, 12])], self._charge_of({10: 0.6, 11: -0.2, 12: -0.2}), "rcd")
        qs, rs, scatter = self.c.assemble_embedding_sites(
            mm_idx, q, xyz, deleted, dq, virt)
        # M1 (index 10) has no direct site; 2 real (M2) + 2 virtual midpoints.
        self.assertEqual(len(qs), 4)
        # a midpoint site is built from M1 and one M2 with 0.5/0.5 weights
        mids = [s for s in scatter if len(s) == 2]
        self.assertEqual(len(mids), 2)
        for s in mids:
            self.assertEqual(s[0][1], 0.5)
            self.assertEqual(s[1][1], 0.5)

    # -- assembled-charge conservation (through assemble, not redistribute) --

    def _assembled_total(self, mm_idx, mm_q, mm_xyz, frontier, scheme):
        charge_of = lambda a: mm_q[mm_idx.index(a)]
        deleted, dq, virt = self.c.redistribute_frontier_charges(
            frontier, charge_of, scheme)
        qs, _, _ = self.c.assemble_embedding_sites(
            mm_idx, mm_q, mm_xyz, deleted, dq, virt)
        return float(np.sum(qs))

    def test_assembled_charge_conserved_isolated_host(self):
        mm_idx = [10, 11, 12, 13]
        mm_q = [0.5, -0.3, -0.2, 0.11]
        xyz = np.array([[0, 0, 0], [1.5, 0, 0], [0, 1.5, 0], [3, 0, 0]], float)
        for scheme in ("rc", "rcd"):
            self.assertAlmostEqual(
                self._assembled_total(mm_idx, mm_q, xyz, [(10, [11, 12])], scheme),
                sum(mm_q), places=12, msg=scheme)

    def test_assembled_charge_conserved_adjacent_hosts(self):
        # Two frontier hosts bonded to each other (adjacent cuts): each is in the
        # other's M2 list, so an RCD delta_q lands on a deleted atom. The residual
        # -charge site must keep total charge conserved.
        mm_idx = [10, 11, 12, 13]
        mm_q = [0.5, 0.3, -0.2, -0.1]
        xyz = np.array([[0, 0, 0], [1.5, 0, 0], [-1.4, 0.3, 0], [2.9, 0.2, 0]], float)
        frontier = [(10, [11, 12]), (11, [10, 13])]
        for scheme in ("rc", "rcd"):
            self.assertAlmostEqual(
                self._assembled_total(mm_idx, mm_q, xyz, frontier, scheme),
                sum(mm_q), places=12, msg=scheme)

    # -- gradient consistency (the decisive test) -----------------------

    def _scenario(self):
        # 2 QM centres (fixed ESP charges/positions); MM = M1(0)+M2a(1)+M2b(2)
        # + 2 bulk (3,4). M1 deleted & redistributed by RCD.
        qA = np.array([0.3, -0.25])
        rA = np.array([[0.2, 0.1, -0.3], [-0.7, 0.4, 0.5]])
        mm_idx = [0, 1, 2, 3, 4]
        mm_q = [0.55, -0.15, -0.18, 0.2, -0.31]
        mm_xyz = np.array([
            [1.4, 0.0, 0.0],    # M1
            [2.6, 0.3, -0.2],   # M2a
            [1.1, 1.5, 0.4],    # M2b
            [3.5, -1.0, 1.2],   # bulk
            [-2.0, -1.5, 0.8],  # bulk
        ], float)
        return qA, rA, mm_idx, mm_q, mm_xyz

    def test_coupling_force_matches_finite_difference(self):
        qA, rA, mm_idx, mm_q, mm_xyz = self._scenario()
        row_of = {m: j for j, m in enumerate(mm_idx)}
        charge_of = self._charge_of({i: mm_q[j] for j, i in enumerate(mm_idx)})

        # isolated host (0 -> M2 {1,2}), and two adjacent hosts (0,1 bonded)
        for frontier in ([(0, [1, 2])], [(0, [1, 3]), (1, [0, 4])]):
          for scheme in ("none", "z1", "rc", "rcd"):
            deleted, dq, virt = self.c.redistribute_frontier_charges(
                frontier, charge_of, scheme)

            def energy(xyz):
                qs, rs, _ = self.c.assemble_embedding_sites(
                    mm_idx, mm_q, xyz, deleted, dq, virt)
                return _coulomb_energy(qA, rA, qs, rs)

            qs, rs, scatter = self.c.assemble_embedding_sites(
                mm_idx, mm_q, mm_xyz, deleted, dq, virt)
            f_mm = _analytic_mm_forces(qA, rA, qs, rs, scatter, len(mm_idx), row_of)

            h = 1e-6
            fd = np.zeros_like(mm_xyz)
            for j in range(len(mm_idx)):
                for c in range(3):
                    xp = mm_xyz.copy(); xp[j, c] += h
                    xm = mm_xyz.copy(); xm[j, c] -= h
                    fd[j, c] = -(energy(xp) - energy(xm)) / (2 * h)
            self.assertTrue(
                np.allclose(f_mm, fd, atol=1e-6),
                f"scheme={scheme}: analytic vs FD mismatch\n{f_mm - fd}")


if __name__ == "__main__":
    unittest.main()
