"""Regression test: the QM/MM driver normalizes qm_atoms into topology order.

The QM geometry handed to the engine is built in topology (ascending atom.index)
order, so the returned QM gradient / coupling force / ESP charges are in that
order. The force-scatter loops and the link-atom host_row projection index those
arrays by position in self.qm_atoms, so qm_atoms must be sorted ascending or the
QM forces land on the wrong atoms. The driver sorts qm_atoms at construction;
this test pins that invariant.

Requires OpenMM and the compiled OpenQP backend (the driver imports it at module
load), so it is skipped in environments without them.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

try:
    import openmm.app as app
    from oqp.library.qmmm_driver import OpenQpQMMM
    _HAVE = True
except Exception:  # pragma: no cover - optional deps / uncompiled backend
    _HAVE = False


@unittest.skipUnless(_HAVE, "OpenMM or compiled OpenQP backend unavailable")
class TestQMAtomsOrder(unittest.TestCase):
    def test_qm_atoms_sorted_into_topology_order(self):
        pdb_path = ROOT / "examples" / "QMMM" / "ala.pdb"
        if not pdb_path.exists():
            self.skipTest("examples/QMMM/ala.pdb missing")
        pdb = app.PDBFile(str(pdb_path))
        try:
            ff = app.ForceField("amber14-all.xml")
        except Exception as err:  # pragma: no cover
            self.skipTest(f"amber14-all.xml unavailable: {err}")
        # deliberately shuffled selection (amide atoms, out of order)
        shuffled = [16, 8, 18, 9, 17]
        drv = OpenQpQMMM(
            pdb.positions, pdb.topology, ff, shuffled, oqp_cfg={},
            Cutoff=app.NoCutoff, Embedding="electrostatic",
        )
        self.assertEqual(list(drv.qm_atoms), sorted(shuffled))

    def test_full_espf_mm_system_has_no_qm_mm_exception_charges(self):
        """Full ESPF carries the whole QM-MM electrostatics; the MM system must
        not keep the scaled 1-4 Coulomb pairs across the covalent boundary
        (OpenMM stores exception charge products independently of the
        particle charges that forces_mm zeroes).  The split scheme keeps them."""
        import openmm as mm
        import openmm.unit as unit
        pdb_path = ROOT / "examples" / "QMMM" / "ala.pdb"
        if not pdb_path.exists():
            self.skipTest("examples/QMMM/ala.pdb missing")
        pdb = app.PDBFile(str(pdb_path))
        try:
            ff = app.ForceField("amber14-all.xml")
        except Exception as err:  # pragma: no cover
            self.skipTest(f"amber14-all.xml unavailable: {err}")
        qm = [8, 9, 16, 17, 18]
        seen = {}
        for embedding in ("electrostatic", "split"):
            drv = OpenQpQMMM(pdb.positions, pdb.topology, ff, qm, oqp_cfg={},
                             Cutoff=app.NoCutoff, Embedding=embedding)
            nb = next(f for f in drv.mm_systems["sys0"].getForces()
                      if isinstance(f, mm.NonbondedForce))
            qm_set = set(qm)
            cross = []
            for i in range(nb.getNumExceptions()):
                p1, p2, cp, _, _ = nb.getExceptionParameters(i)
                if (p1 in qm_set) != (p2 in qm_set):
                    cross.append(cp.value_in_unit(unit.elementary_charge ** 2))
            seen[embedding] = sum(1 for c in cross if c != 0.0)
        self.assertEqual(seen["electrostatic"], 0)
        self.assertEqual(seen["split"], 13)        # the alanine partition has 13 QM-MM 1-4 pairs

    def test_compute_force_restores_topology_order(self):
        """compute_force() receives a caller-supplied qm_atoms (QMMM_MD passes
        its own config-order list) and must re-sort it before scattering forces.
        Ported from PR #270 (mock-based, no OpenMM system needed)."""
        from unittest import mock
        drv = object.__new__(OpenQpQMMM)
        drv.Embedding = "electrostatic"
        drv.espf_full = True
        drv.Cutoff = app.NoCutoff                      # non-periodic: no Ewald branch
        drv.electrostatic_potential = mock.Mock(return_value=(None, None))
        drv.forces_qm_openqp = mock.Mock(return_value=(0.0, None, object()))
        drv._assemble_force_espf = mock.Mock(return_value=(0.0, None))

        shuffled = [16, 8, 18, 9, 17]
        drv.compute_force(None, None, None, shuffled)

        self.assertEqual(list(drv.qm_atoms), sorted(shuffled))
        drv._assemble_force_espf.assert_called_once()


if __name__ == "__main__":
    unittest.main()
