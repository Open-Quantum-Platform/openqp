"""The ORCA-style active/frozen atom selection shared by the QM/MM drivers.

These exercise ``oqp.library.qmmm_active`` directly, so they need OpenMM (for
the topology and the PDB reader) but not the compiled OpenQP backend, and they
therefore run in CI's Python leg.  The driver-side wiring is covered by
tests/test_qmmm_optimization.py.
"""
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

try:
    from openmm import Vec3, app, unit

    from oqp.library.qmmm_active import (
        active_from_pdb_file, freeze_constrained_partners, parse_atom_selection,
        resolve_active_set, selection_requested)
    _HAVE = True
except Exception:  # pragma: no cover - OpenMM missing
    _HAVE = False


def _topology():
    """Two QM carbons, two MM atoms of the same residue, and two waters at 3 and
    8 A -- the same layout the optimiser tests use."""
    top = app.Topology()
    chain = top.addChain()
    pos = []
    res = top.addResidue("QMA", chain)
    for name, p in (("C1", (0, 0, 0)), ("C2", (1.2, 0, 0)),
                    ("CA", (-1.5, 0, 0)), ("CB", (-4.5, 0, 0))):
        top.addAtom(name, app.element.carbon, res)
        pos.append(p)
    for x in (3.0, 8.0):
        water = top.addResidue("HOH", chain)
        for name, element, dy in (("O", app.element.oxygen, 0.0),
                                  ("H1", app.element.hydrogen, 0.6),
                                  ("H2", app.element.hydrogen, -0.6)):
            top.addAtom(name, element, water)
            pos.append((x, dy, 0.0))
    return top, np.array(pos, dtype=float)


@unittest.skipUnless(_HAVE, "OpenMM unavailable")
class TestSelectionSyntax(unittest.TestCase):
    def test_orca_and_openqp_spellings_agree(self):
        """ORCA writes {0:5 16}, we wrote 0-5,16; a deck may use either."""
        atoms = list(_topology()[0].atoms())
        expected = {0, 1, 2, 3, 4, 8}
        for spec in ("{0:3 4 8}", "0-3,4,8", "0:3 4 8", "0-3 4 8", "0,1,2,3,4,8",
                     "{0:3};4;8"):
            self.assertEqual(parse_atom_selection(spec, atoms, "active_atoms"),
                             expected, spec)

    def test_name_groups_match_every_matching_atom(self):
        atoms = list(_topology()[0].atoms())
        self.assertEqual(parse_atom_selection("name:O", atoms, "active_atoms"), {4, 7})
        self.assertEqual(parse_atom_selection("name:H1,H2", atoms, "active_atoms"),
                         {5, 6, 8, 9})
        # a name group and indices may be mixed in one selection
        self.assertEqual(parse_atom_selection("name:O 0-1", atoms, "active_atoms"),
                         {0, 1, 4, 7})

    def test_blank_selection_is_not_a_selection(self):
        atoms = list(_topology()[0].atoms())
        self.assertEqual(parse_atom_selection("", atoms, "active_atoms"), set())
        for cfg in ({}, {"active_atoms": "", "frozen_atoms": "", "active_radius": 0.0,
                         "active_from_pdb": False},
                    {"active_radius": "0.0"}):
            self.assertFalse(selection_requested(cfg), cfg)
        for cfg in ({"active_atoms": "0-3"}, {"frozen_atoms": "7"},
                    {"active_radius": 3.0}, {"active_from_pdb": True}):
            self.assertTrue(selection_requested(cfg), cfg)


@unittest.skipUnless(_HAVE, "OpenMM unavailable")
class TestResolveActiveSet(unittest.TestCase):
    def test_dynamics_propagates_everything_unless_asked(self):
        """default_all=True is what the MD and NAMD drivers pass: a deck that
        says nothing keeps moving every atom, as it always did."""
        top, pos = _topology()
        active, frozen = resolve_active_set({}, top, pos, [0, 1], default_all=True)
        np.testing.assert_array_equal(active, np.arange(10))
        self.assertEqual(frozen, set())

    def test_optimisation_moves_the_qm_region_unless_asked(self):
        top, pos = _topology()
        active, frozen = resolve_active_set({}, top, pos, [0, 1], default_all=False)
        np.testing.assert_array_equal(active, [0, 1])
        self.assertEqual(frozen, set())

    def test_frozen_alone_holds_only_those_atoms(self):
        """The headline use: hold a backbone and let everything else move."""
        top, pos = _topology()
        active, frozen = resolve_active_set(
            {"frozen_atoms": "7-9"}, top, pos, [0, 1], default_all=True)
        np.testing.assert_array_equal(active, [0, 1, 2, 3, 4, 5, 6])
        self.assertEqual(frozen, {7, 8, 9})

    def test_radius_takes_whole_mm_residues(self):
        top, pos = _topology()
        active, _ = resolve_active_set(
            {"active_radius": 2.0}, top, pos, [0, 1], default_all=False)
        np.testing.assert_array_equal(active, [0, 1, 2, 4, 5, 6])

    def test_the_qm_region_may_not_be_frozen(self):
        top, pos = _topology()
        with self.assertRaises(ValueError) as err:
            resolve_active_set({"frozen_atoms": "1"}, top, pos, [0, 1], default_all=True)
        self.assertIn("names QM atoms", str(err.exception))

    def test_a_negative_radius_is_rejected(self):
        top, pos = _topology()
        with self.assertRaises(ValueError) as err:
            resolve_active_set({"active_radius": -1.0}, top, pos, [0, 1], default_all=True)
        self.assertIn("finite distance >= 0 angstrom", str(err.exception))


@unittest.skipUnless(_HAVE, "OpenMM unavailable")
class TestActiveFromPdb(unittest.TestCase):
    """ORCA's Use_Active_InfoFromPDB: B-factor 1 marks an active atom.

    OpenMM's PDBFile drops that column, which is why qmmm_active reads the file
    again through pdbstructure -- if that ever regressed, every atom would come
    back with a B-factor of 0 and the selection would be empty.
    """

    PDB = ROOT / "examples" / "QMMM" / "formaldehyde_water_active.pdb"

    def test_the_example_pdb_marks_the_qm_region_and_the_nearest_water(self):
        atoms = list(app.PDBFile(str(self.PDB)).topology.atoms())
        self.assertEqual(len(atoms), 19)
        self.assertEqual(active_from_pdb_file(str(self.PDB), atoms), set(range(7)))

    def test_it_agrees_with_writing_the_same_selection_out(self):
        """The deck's comment offers active_atoms=0-6 as the explicit spelling
        of what the B-factor column says; they must select the same atoms."""
        pdb = app.PDBFile(str(self.PDB))
        top = pdb.topology
        pos = np.array(pdb.positions.value_in_unit(unit.angstrom))
        from_pdb, _ = resolve_active_set(
            {"active_from_pdb": True}, top, pos, [0, 1, 2, 3],
            default_all=True, pdb_path=str(self.PDB))
        explicit, _ = resolve_active_set(
            {"active_atoms": "0-6"}, top, pos, [0, 1, 2, 3], default_all=True)
        orca_spelling, _ = resolve_active_set(
            {"active_atoms": "{0:6}"}, top, pos, [0, 1, 2, 3], default_all=True)
        np.testing.assert_array_equal(from_pdb, np.arange(7))
        np.testing.assert_array_equal(explicit, np.arange(7))
        np.testing.assert_array_equal(orca_spelling, np.arange(7))

    def test_an_unmarked_pdb_is_an_error_not_an_empty_selection(self):
        plain = ROOT / "examples" / "QMMM" / "formaldehyde_water.pdb"
        atoms = list(app.PDBFile(str(plain)).topology.atoms())
        with self.assertRaises(ValueError) as err:
            active_from_pdb_file(str(plain), atoms)
        self.assertIn("is 0 everywhere", str(err.exception))


@unittest.skipUnless(_HAVE, "OpenMM unavailable")
class TestFrozenAtomsAndConstraints(unittest.TestCase):
    def test_a_constraint_never_ties_a_moving_atom_to_a_fixed_one(self):
        pairs = [(4, 5), (4, 6), (5, 6)]          # O-H, O-H, H-H of one water
        active, frozen = freeze_constrained_partners(pairs, [0, 1, 5, 6], {4})
        self.assertEqual((active, frozen), ({0, 1}, {4, 5, 6}))
        active, frozen = freeze_constrained_partners(pairs, [0, 1, 4], set())
        self.assertEqual((active, frozen), ({0, 1, 4, 5, 6}, set()))
        # nothing straddles the boundary: the sets come back untouched
        active, frozen = freeze_constrained_partners(pairs, [0, 1, 4, 5, 6], set())
        self.assertEqual((active, frozen), ({0, 1, 4, 5, 6}, set()))

    def test_the_move_mask_holds_frozen_rows_through_velocity_verlet(self):
        """The NAMD drivers hold an atom by multiplying its update by 0.  A held
        row must not move and must not acquire velocity, whatever force acts on
        it -- a frozen atom still feels the field, it just does not respond."""
        mask = np.array([[1.0], [1.0], [0.0], [0.0]])
        r = np.arange(12, dtype=float).reshape(4, 3)
        v = np.ones((4, 3))
        accel = np.full((4, 3), 0.7)            # a real force on every atom
        dt = 0.5
        r_new = r + (v * dt + 0.5 * accel * dt ** 2) * mask
        v_new = v + 0.5 * (accel + accel) * dt * mask
        np.testing.assert_array_equal(r_new[2:], r[2:])
        np.testing.assert_array_equal(v_new[2:], v[2:])
        self.assertTrue(np.all(r_new[:2] != r[:2]))
        self.assertTrue(np.all(v_new[:2] != v[:2]))


if __name__ == "__main__":
    unittest.main()
