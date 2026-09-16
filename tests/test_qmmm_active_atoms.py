"""The ORCA-style active/frozen atom selection shared by the QM/MM drivers.

These exercise ``oqp.library.qmmm_active`` directly, so they need OpenMM (for
the topology and the PDB reader) but not the compiled OpenQP backend, and they
therefore run in CI's Python leg.  The driver-side wiring is covered by
tests/test_qmmm_optimization.py.
"""
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PDB_ACTIVE = ROOT / "examples" / "QMMM" / "formaldehyde_water_active.pdb"

try:
    from openmm import Vec3, app, unit

    from oqp.library.qmmm_active import (
        active_from_pdb_file, freeze_constrained_partners, held_atoms,
        parse_atom_selection, resolve_active_set, selection_requested)
    _HAVE = True
except Exception:  # pragma: no cover - OpenMM missing
    _HAVE = False

try:  # the driver classes additionally need the compiled OpenQP backend
    from oqp.library.namd import NAMD_QMMM
    from oqp.library.qmmm_md import QMMM_MD
    _HAVE_DRIVERS = True
except Exception:  # pragma: no cover - backend missing
    _HAVE_DRIVERS = False


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


@unittest.skipUnless(_HAVE, "OpenMM unavailable")
class TestHeldIsTheComplementOfActive(unittest.TestCase):
    """What a driver must hold is everything outside the active set.

    ``frozen_atoms`` is only what the deck named explicitly, so a driver that
    holds *that* set propagates every atom an ``active_atoms`` /
    ``active_radius`` / ``active_from_pdb`` selection failed to select -- i.e.
    the selection silently does nothing.
    """

    def test_an_active_only_selection_still_holds_everything_else(self):
        top, pos = _topology()
        active, frozen = resolve_active_set(
            {"active_atoms": "0-3"}, top, pos, [0, 1], default_all=True)
        self.assertEqual(frozen, set())            # nothing was frozen BY NAME
        self.assertEqual(held_atoms(top, active), {4, 5, 6, 7, 8, 9})

    def test_a_radius_selection_holds_the_residues_it_did_not_reach(self):
        top, pos = _topology()
        active, frozen = resolve_active_set(
            {"active_radius": 2.0}, top, pos, [0, 1], default_all=True)
        self.assertEqual(frozen, set())
        self.assertEqual(held_atoms(top, active), {3, 7, 8, 9})

    def test_frozen_atoms_alone_holds_exactly_those(self):
        top, pos = _topology()
        active, _ = resolve_active_set(
            {"frozen_atoms": "7-9"}, top, pos, [0, 1], default_all=True)
        self.assertEqual(held_atoms(top, active), {7, 8, 9})

    def test_nothing_is_held_when_every_atom_is_active(self):
        top, pos = _topology()
        active, _ = resolve_active_set({}, top, pos, [0, 1], default_all=True)
        self.assertEqual(held_atoms(top, active), set())


@unittest.skipUnless(_HAVE, "OpenMM unavailable")
class TestAtomZeroIsARealSelection(unittest.TestCase):
    """Index 0 must never be mistaken for "nothing selected" -- `str(0 or "")`
    is `""`, which would silently drop a valid request."""

    def test_selection_requested_sees_index_zero_as_text_or_number(self):
        for cfg in ({"frozen_atoms": "0"}, {"frozen_atoms": 0},
                    {"active_atoms": "0"}, {"active_atoms": 0}):
            self.assertTrue(selection_requested(cfg), cfg)

    def test_parsing_accepts_index_zero_in_either_type(self):
        atoms = list(_topology()[0].atoms())
        for spec in ("0", 0):
            self.assertEqual(parse_atom_selection(spec, atoms, "frozen_atoms"), {0})

    def test_a_zero_radius_is_still_no_shell(self):
        for cfg in ({"active_radius": 0.0}, {"active_radius": "0.0"},
                    {"active_radius": 0}, {}):
            self.assertFalse(selection_requested(cfg), cfg)
        self.assertTrue(selection_requested({"active_radius": 2.5}))


class TestPositionUpdateStaysBitIdentical(unittest.TestCase):
    """With nothing held, the masked velocity-Verlet update must reproduce the
    old expression bit-for-bit.  Multiplying by an exact 1.0 is exact, but the
    terms have to keep their original association to stay so."""

    def test_masking_each_term_is_bit_identical_but_regrouping_is_not(self):
        r = np.array([[10.0, 10.0, 10.0]])
        v = np.full((1, 3), 1.0e-4)
        a = np.full((1, 3), 1.0e-5)
        dt = 20.0
        mask = np.ones((1, 1))
        old = r + v * dt + 0.5 * a * dt ** 2
        masked = r + v * dt * mask + 0.5 * a * dt ** 2 * mask
        regrouped = r + (v * dt + 0.5 * a * dt ** 2) * mask
        self.assertTrue(np.array_equal(old, masked))
        # control: the grouping this replaced really did move the number, so the
        # assertion above is not vacuous
        self.assertFalse(np.array_equal(old, regrouped))

    def test_a_held_row_does_not_move_and_gains_no_velocity(self):
        mask = np.array([[1.0], [0.0]])
        r = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 5.0]])
        v = np.array([[0.2, 0.0, 0.0], [0.0, 0.0, 0.0]])
        a = np.full((2, 3), 0.7)                   # a real force on both atoms
        dt = 0.5
        r_new = r + v * dt * mask + 0.5 * a * dt ** 2 * mask
        v_new = v + 0.5 * (a + a) * dt * mask
        np.testing.assert_array_equal(r_new[1], r[1])
        np.testing.assert_array_equal(v_new[1], 0.0)
        self.assertTrue(np.any(r_new[0] != r[0]))


@unittest.skipUnless(_HAVE and _HAVE_DRIVERS, "OpenMM or OpenQP backend unavailable")
class TestNamdHeldAtomBookkeeping(unittest.TestCase):
    def _bare(self, mask):
        o = object.__new__(NAMD_QMMM)
        o._move_mask = np.asarray(mask, dtype=float).reshape(-1, 1)
        o._all_atoms_move = bool(o._move_mask.all())
        o.natom_all = len(o._move_mask)
        o.m_all = np.array([1.0, 2.0, 3.0, 4.0])[:o.natom_all]
        return o

    def test_held_rows_are_zeroed_and_the_moving_com_is_put_at_rest(self):
        o = self._bare([1, 1, 0, 0])
        o.v_all = np.array([[1.0, 0.0, 0.0], [-0.4, 0.2, 0.0],
                            [5.0, 5.0, 5.0], [7.0, 0.0, 0.0]])
        o._zero_held_velocities()
        np.testing.assert_array_equal(o.v_all[2:], 0.0)
        moving_mass = o.m_all * o._move_mask[:, 0]
        momentum = (moving_mass[:, None] * o.v_all).sum(axis=0)
        np.testing.assert_allclose(momentum, 0.0, atol=1e-14)

    def test_it_leaves_an_all_moving_system_exactly_alone(self):
        o = self._bare([1, 1, 1, 1])
        v = np.array([[1.0, 0.0, 0.0], [-0.4, 0.2, 0.0],
                      [5.0, 5.0, 5.0], [7.0, 0.0, 0.0]])
        o.v_all = v.copy()
        o._zero_held_velocities()
        np.testing.assert_array_equal(o.v_all, v)

    def test_a_held_atom_does_not_shrink_the_adaptive_timestep(self):
        """dt_adaptive sizes the step from the largest predicted displacement.
        A held atom has no displacement, so the force on it must not slow the
        trajectory down."""
        o = self._bare([1, 0])
        o.dt_adaptive = True
        o.dt_max, o.dt_min, o.dx_max = 20.0, 1.0, 0.01
        v = np.zeros((2, 3))
        accel = np.zeros((2, 3))
        accel[1, 0] = 10.0                      # a large force on the HELD atom
        # control: unmasked, that force really does drag the timestep down, so
        # the assertion below is not vacuous
        self.assertLess(o._adaptive_dt(v, accel), o.dt_max)
        self.assertEqual(
            o._adaptive_dt(v * o._move_mask, accel * o._move_mask), o.dt_max)

    def test_a_restart_identity_treats_atom_zero_as_a_selection(self):
        base = NAMD_QMMM._qmmm_identity_config({"pdb_file": "x.pdb"})
        zero = NAMD_QMMM._qmmm_identity_config(
            {"pdb_file": "x.pdb", "frozen_atoms": "0"})
        self.assertEqual(zero.get("frozen_atoms"), "0")
        self.assertNotEqual(base, zero)

    def test_unset_selection_keys_leave_the_restart_identity_untouched(self):
        """Checkpoints written before these keys existed must keep validating."""
        base = NAMD_QMMM._qmmm_identity_config({"pdb_file": "x.pdb"})
        defaults = NAMD_QMMM._qmmm_identity_config({
            "pdb_file": "x.pdb", "active_atoms": "", "frozen_atoms": "",
            "active_radius": 0.0, "active_from_pdb": False})
        self.assertEqual(base, defaults)


@unittest.skipUnless(_HAVE and _HAVE_DRIVERS, "OpenMM or OpenQP backend unavailable")
class TestOpenMMMdHoldsTheRightAtoms(unittest.TestCase):
    """runtype=md holds an atom by giving it zero mass, so it must be handed the
    complement of the active set rather than the explicitly frozen set."""

    def _bare(self, cfg):
        o = object.__new__(QMMM_MD)
        o.pdb = app.PDBFile(str(PDB_ACTIVE))
        o._pdb_path = str(PDB_ACTIVE)
        o.qm_atoms = np.array([0, 1, 2, 3], dtype=int)
        o._selection_cfg = cfg
        o.rigidwater = False
        o.oqp_driver = SimpleNamespace(_box_lengths_bohr=lambda: None)
        return o

    def test_an_active_only_selection_holds_the_unselected_waters(self):
        self.assertEqual(self._bare({"active_atoms": "0-6"})._resolve_frozen_atoms(),
                         set(range(7, 19)))

    def test_active_from_pdb_holds_the_atoms_the_b_factors_did_not_mark(self):
        self.assertEqual(self._bare({"active_from_pdb": True})._resolve_frozen_atoms(),
                         set(range(7, 19)))

    def test_frozen_atoms_alone_holds_exactly_those(self):
        self.assertEqual(self._bare({"frozen_atoms": "7-9"})._resolve_frozen_atoms(),
                         {7, 8, 9})

    def test_no_selection_holds_nothing(self):
        self.assertEqual(self._bare({})._resolve_frozen_atoms(), set())


if __name__ == "__main__":
    unittest.main()
