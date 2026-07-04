"""Integration test: QM/MM frontier (M1) charge redistribution on a real force
field and covalent boundary.

Builds the alanine dipeptide (ACE-ALA-NH2) with the AMBER-14 force field via
OpenMM, cuts the QM/MM boundary through the ALA C-CA backbone bond (QM = the
C-terminal amide), and checks that:

  * the cut is detected as exactly one hydrogen link atom, and
  * RCD / RC frontier-charge redistribution conserves the total MM charge that
    the QM density is embedded in (Z1 deletes and so does not).

This exercises the same pure connectivity core the QM/MM driver uses, on real
AMBER partial charges rather than synthetic ones. It is skipped when OpenMM (or
the bundled amber14-all.xml) is unavailable, so it does not require the compiled
OpenQP backend.
"""

import importlib.util
import os
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u
    _HAVE_OPENMM = True
except Exception:  # pragma: no cover - optional dependency
    _HAVE_OPENMM = False


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@unittest.skipUnless(_HAVE_OPENMM, "OpenMM not installed")
class TestFrontierChargeOnAmberAlanine(unittest.TestCase):
    def setUp(self):
        self.c = load_module(
            "qmmm_connectivity_amber_under_test",
            "pyoqp/oqp/library/qmmm_connectivity.py",
        )
        pdb_path = ROOT / "examples" / "QMMM" / "ala.pdb"
        if not pdb_path.exists():
            self.skipTest("examples/QMMM/ala.pdb missing")
        self.pdb = app.PDBFile(str(pdb_path))
        try:
            ff = app.ForceField("amber14-all.xml")
            # Keep the System referenced: its Force objects become invalid once
            # the owning System is garbage-collected.
            self.system = ff.createSystem(self.pdb.topology)
        except Exception as err:  # pragma: no cover
            self.skipTest(f"amber14-all.xml unavailable: {err}")
        self.nb = next(f for f in self.system.getForces()
                       if isinstance(f, mm.NonbondedForce))
        self.charge = {
            i: self.nb.getParticleParameters(i)[0].value_in_unit(u.elementary_charge)
            for i in range(self.nb.getNumParticles())
        }
        self.z = {a.index: a.element.atomic_number for a in self.pdb.topology.atoms()}
        self.bonds = [(b[0].index, b[1].index) for b in self.pdb.topology.bonds()]
        self.name = {a.index: (a.residue.name, a.name)
                     for a in self.pdb.topology.atoms()}
        # QM = C-terminal amide (ALA C, ALA O + all NH2 atoms): cuts only ALA CA-C.
        self.qm = [i for i, n in self.name.items()
                   if n[0] == "NH2" or n in (("ALA", "C"), ("ALA", "O"))]

    def _frontier(self, links):
        qm_set = set(self.qm)
        nbrs = {}
        for i, j in self.bonds:
            nbrs.setdefault(i, set()).add(j)
            nbrs.setdefault(j, set()).add(i)
        hosts = {l.mm_index: sorted(x for x in nbrs[l.mm_index] if x not in qm_set)
                 for l in links}
        return [(m1, hosts[m1]) for m1 in sorted(hosts)]

    def test_single_link_atom_on_backbone_cut(self):
        links = self.c.detect_link_atoms(self.bonds, self.qm, lambda i: self.z[i])
        self.assertEqual(len(links), 1)
        la = links[0]
        self.assertEqual(self.name[la.qm_index], ("ALA", "C"))
        self.assertEqual(self.name[la.mm_index], ("ALA", "CA"))
        self.assertTrue(0.0 < la.g < 1.0)

    def test_rcd_and_rc_conserve_amber_charge(self):
        links = self.c.detect_link_atoms(self.bonds, self.qm, lambda i: self.z[i])
        frontier = self._frontier(links)
        mm_idx = [i for i in range(self.nb.getNumParticles()) if i not in set(self.qm)]
        mm_q = [self.charge[i] for i in mm_idx]
        mm_xyz = np.zeros((len(mm_idx), 3))
        raw = sum(mm_q)
        for scheme in ("rcd", "rc"):
            deleted, dq, virt = self.c.redistribute_frontier_charges(
                frontier, lambda a: self.charge[a], scheme)
            qs, _, _ = self.c.assemble_embedding_sites(
                mm_idx, mm_q, mm_xyz, deleted, dq, virt)
            self.assertAlmostEqual(float(qs.sum()), raw, places=9, msg=scheme)

    def test_z1_deletes_host_charge(self):
        links = self.c.detect_link_atoms(self.bonds, self.qm, lambda i: self.z[i])
        frontier = self._frontier(links)
        (m1, _), = frontier
        mm_idx = [i for i in range(self.nb.getNumParticles()) if i not in set(self.qm)]
        mm_q = [self.charge[i] for i in mm_idx]
        deleted, dq, virt = self.c.redistribute_frontier_charges(
            frontier, lambda a: self.charge[a], "z1")
        qs, _, _ = self.c.assemble_embedding_sites(
            mm_idx, mm_q, np.zeros((len(mm_idx), 3)), deleted, dq, virt)
        self.assertAlmostEqual(float(qs.sum()), sum(mm_q) - self.charge[m1], places=9)


if __name__ == "__main__":
    unittest.main()
