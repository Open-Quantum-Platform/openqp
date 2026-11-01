"""Default frontier-charge treatment at a covalent QM/MM cut (issue #551).

The link hydrogen sits on the cut QM-M1 bond, 0.4-0.5 A from the MM host
M1.  With frontier_scheme=none the full M1 charge stays in the embedding as a
point charge inside the link-atom density, the ESPF energy is unbounded below,
and optimization and MD collapse the QM frontier atom onto M1 (rhodopsin
Lys296 CB|CA in the issue; the ACE CH3|C cut of the alanine dipeptide falls
from -41.5 to -249 Eh in ten RHF/6-31G steps).  These tests pin that the
default removes every embedding charge from the link-atom neighbourhood, and
that an explicit 'none' is reported.

Loaded by path, like test_qmmm_frontier.py, so they run without the compiled
OpenQP backend or OpenMM.
"""

import importlib.util
import re
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

_spec = importlib.util.spec_from_file_location(
    "qmmm_frontier_default_under_test",
    ROOT / "pyoqp/oqp/library/qmmm_connectivity.py",
)
c = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = c
_spec.loader.exec_module(c)

# ACE CH3|C cut of the alanine dipeptide (examples/QMMM/ala.pdb, Angstrom):
# QM frontier CH3 carbon, M1 = ACE C (AMBER +0.5972), M2 = ACE O and ALA N.
# The link H is ~0.40 A from M1 here.
QM_C = np.array([-0.028, 0.101, -0.091])
M1 = np.array([-0.023, 0.157, 1.419])
M2 = {2: np.array([1.014, 0.404, 2.021]), 6: np.array([-1.183, -0.077, 2.026])}
Q = {1: 0.5972, 2: -0.5679, 6: -0.4157}


def _sites(scheme):
    deleted, dq, virtuals = c.redistribute_frontier_charges(
        [(1, sorted(M2))], Q.__getitem__, scheme)
    idx = np.array([1, 2, 6])
    pos = np.array([M1, M2[2], M2[6]])
    charges, positions, _ = c.assemble_embedding_sites(
        idx, np.array([Q[i] for i in idx]), pos, deleted, dq, virtuals)
    return charges, positions


def _closest_charge_to_link(scheme):
    g = c.link_g_factor(6, 6)
    link = c.link_atom_position(QM_C, M1, g)
    charges, positions = _sites(scheme)
    d = np.linalg.norm(positions - link, axis=1)
    keep = np.abs(charges) > 1e-12
    return float(d[keep].min())


class TestFrontierDefault(unittest.TestCase):
    def test_default_redistributes_the_host_charge(self):
        self.assertEqual(c.DEFAULT_FRONTIER_SCHEME, "rcd")

    def test_input_schema_uses_the_same_default(self):
        src = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()
        m = re.search(r"'frontier_scheme':\s*\{[^}]*'default':\s*'([a-z0-9]+)'", src)
        self.assertIsNotNone(m)
        self.assertEqual(m.group(1), c.DEFAULT_FRONTIER_SCHEME)

    def test_none_leaves_the_host_charge_inside_the_link_atom(self):
        # The cause of #551: the M1 charge 0.4-0.5 A from the link hydrogen.
        self.assertLess(_closest_charge_to_link("none"), 0.5)

    def test_default_keeps_charges_away_from_the_link_atom(self):
        # The nearest remaining sites are the M1-M2 bond midpoints, which
        # are farther from the link atom than M1 is.
        d = _closest_charge_to_link(c.DEFAULT_FRONTIER_SCHEME)
        # Measured on this geometry: 0.40 A for 'none', 0.885 A for 'rcd'.
        self.assertGreater(d, 0.8)
        self.assertGreater(d, 2.0 * _closest_charge_to_link("none"))

    def test_default_conserves_the_embedding_charge(self):
        q_none, _ = _sites("none")
        q_def, _ = _sites(c.DEFAULT_FRONTIER_SCHEME)
        self.assertAlmostEqual(q_def.sum(), q_none.sum(), places=12)

    def test_explicit_none_at_a_covalent_cut_is_reported(self):
        msg = c.frontier_scheme_warning(1, "none")
        self.assertIsNotNone(msg)
        self.assertIn("#551", msg)
        self.assertIn(c.DEFAULT_FRONTIER_SCHEME, msg)
        self.assertIsNotNone(c.frontier_scheme_warning(2, " None "))

    def test_no_report_without_cut_or_with_redistribution(self):
        self.assertIsNone(c.frontier_scheme_warning(0, "none"))
        for scheme in ("rcd", "rc", "z1"):
            self.assertIsNone(c.frontier_scheme_warning(1, scheme))


if __name__ == "__main__":
    unittest.main()
