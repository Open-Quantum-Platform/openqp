"""Unit tests for QM/MM dangling-bond (link-atom) connectivity.

These exercise the pure-Python connectivity core used by the QM/MM driver to
cap covalent bonds cut by the QM/MM boundary.  They load the module by path so
they run without the compiled OpenQP backend or OpenMM.
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


class TestQMMMConnectivity(unittest.TestCase):
    def setUp(self):
        self.c = load_module(
            "qmmm_connectivity_under_test",
            "pyoqp/oqp/library/qmmm_connectivity.py",
        )

    # -- detection ------------------------------------------------------

    def test_single_severed_bond_creates_one_link_atom(self):
        # QM = {0(C),1(O),2(N),3(H),4(H)}; bond 0(C)-5(C_MM) is cut.
        z = {0: 6, 1: 8, 2: 7, 3: 1, 4: 1, 5: 6, 6: 1}
        bonds = [(0, 1), (0, 5), (1, 2), (2, 3), (2, 4), (5, 6)]
        links = self.c.detect_link_atoms(bonds, [0, 1, 2, 3, 4], lambda i: z[i])
        self.assertEqual(len(links), 1)
        la = links[0]
        self.assertEqual((la.qm_index, la.mm_index, la.host_row), (0, 5, 0))
        g_expected = (
            self.c.COVALENT_RADII[1] + self.c.COVALENT_RADII[6]
        ) / (self.c.COVALENT_RADII[6] + self.c.COVALENT_RADII[6])
        self.assertAlmostEqual(la.g, g_expected, places=12)

    def test_whole_molecule_qm_region_has_no_link_atoms(self):
        # Two intact waters; QM = first water only -> no bond crosses boundary.
        bonds = [(0, 1), (0, 2), (3, 4), (3, 5)]
        z = lambda i: 8 if i % 3 == 0 else 1
        self.assertEqual(self.c.detect_link_atoms(bonds, [0, 1, 2], z), [])

    def test_multiple_link_atoms_on_one_frontier_atom(self):
        # A QM carbon bonded to two MM carbons -> two link atoms.
        links = self.c.detect_link_atoms(
            [(0, 1), (0, 2)], [0], lambda i: 6
        )
        self.assertEqual(len(links), 2)
        self.assertEqual(
            sorted((l.qm_index, l.mm_index) for l in links), [(0, 1), (0, 2)]
        )

    def test_detection_is_order_independent_and_deterministic(self):
        z = {0: 6, 1: 6, 2: 6}
        a = self.c.detect_link_atoms([(0, 1), (0, 2)], [0], lambda i: z[i])
        b = self.c.detect_link_atoms([(2, 0), (1, 0)], [0], lambda i: z[i])
        self.assertEqual(
            [(l.qm_index, l.mm_index) for l in a],
            [(l.qm_index, l.mm_index) for l in b],
        )

    def test_unphysical_partition_raises(self):
        # Cutting an H-H bond gives g = 1 (link atom coincides with MM host).
        with self.assertRaises(ValueError):
            self.c.detect_link_atoms([(0, 1)], [0], lambda i: 1)

    # -- geometry -------------------------------------------------------

    def test_link_atom_lies_on_qm_mm_axis(self):
        qm = np.array([0.0, 0.0, 0.0])
        mm = np.array([1.5, 0.0, 0.0])
        g = 0.71
        pos = self.c.link_atom_position(qm, mm, g)
        self.assertTrue(np.allclose(pos, [1.5 * g, 0.0, 0.0]))

    # -- gradient projection -------------------------------------------

    def test_gradient_projection_conserves_total(self):
        grad = np.array([0.3, -1.2, 0.7])
        g = 0.65
        to_qm, to_mm = self.c.project_link_gradient(grad, g)
        # Translational invariance: the two shares reconstruct the full grad.
        self.assertTrue(np.allclose(to_qm + to_mm, grad))
        self.assertTrue(np.allclose(to_qm, (1.0 - g) * grad))
        self.assertTrue(np.allclose(to_mm, g * grad))

    def test_g_factor_is_between_zero_and_one_for_typical_bonds(self):
        for qm_z, mm_z in [(6, 6), (6, 7), (7, 6), (8, 6), (6, 16)]:
            g = self.c.link_g_factor(qm_z, mm_z)
            self.assertTrue(0.0 < g < 1.0, (qm_z, mm_z, g))


if __name__ == "__main__":
    unittest.main()
