"""Gates around staging the per-pair irrep table for the Davidson guess.

Three defects, all reachable only once ``[symmetry] enabled`` started defaulting
to true, and all of which failed SILENTLY -- the guess fell back to the
historical seeds and nothing in the output said so. That silence is why they are
pinned here: the failure mode of this mechanism is "the answer looks fine and is
a different state".

1. ``label_mo`` is a DISPLAY flag. It must not decide whether the irreps are
   computed, because the guess coverage needs them.
2. Quintet MRSF runs a transposed response space (beta occupied to alpha
   virtual). Staging the singlet/triplet space produces a table of the wrong
   length, which every Fortran consumer then discards at its size guard.
3. The state-label formatter must read the validated schema key
   ``[tdhf] multiplicity``, not the test-only ``mult`` alias.

These tests drive the real methods against injected data rather than grepping
the source. An earlier version of this file asserted on source text -- that a
call site still reads ``required=True``, that a branch still mentions
``_td_multiplicity() == 5`` -- which would have passed on a dead branch, a
wrong slice, or a comment. Those assertions are gone.
"""

import sys
import types
import unittest
from pathlib import Path

import numpy as np

from test_symmetry_metadata import load_molecule_module

ROOT = Path(__file__).resolve().parents[1]


def _ensure_oqp_library():
    """Make ``oqp.library.symmetry`` importable under the stub environment.

    ``load_molecule_module`` installs a stub ``oqp`` package whose ``__path__``
    is ``pyoqp/`` rather than ``pyoqp/oqp/``, so ``oqp.library`` does not
    resolve and the methods under test bail into their own ``except`` with
    "No module named 'oqp.library'" -- which would look exactly like the
    defect. Point a real path at it; ``symmetry.py`` needs only numpy.
    """
    load_molecule_module()
    if 'oqp.library' not in sys.modules:
        pkg = types.ModuleType('oqp.library')
        pkg.__path__ = [str(ROOT / 'pyoqp' / 'oqp' / 'library')]
        sys.modules['oqp.library'] = pkg


# C2v, in the order stage_response_symmetry indexes it (1-based into this dict).
C2V = {
    'a1': [1, 1, 1, 1],
    'a2': [1, 1, -1, -1],
    'b1': [1, -1, 1, -1],
    'b2': [1, -1, -1, 1],
}


class FakeData(dict):
    """Stand-in for Molecule.data: a tag store that records what was staged."""


def _molecule(config, metadata=None, data=None):
    _ensure_oqp_library()
    molecule_module = load_molecule_module()
    mol = molecule_module.Molecule.__new__(molecule_module.Molecule)
    mol.config = config
    mol.symmetry_metadata = metadata if metadata is not None else {}
    mol.data = data if data is not None else FakeData()
    return mol


def _staging_molecule(multiplicity, na=5, nb=3, nbf=8):
    """A molecule ready for stage_response_symmetry, with MO labels supplied.

    na/nb/nbf are chosen so the two candidate response spaces have DIFFERENT
    sizes -- na*(nbf-nb) = 25 and nb*(nbf-na) = 9 -- so the assertion can tell
    which one was staged rather than merely that something was staged.
    """
    labels = ['a1', 'b2', 'a1', 'b1', 'a1', 'b2', 'a1', 'b1'][:nbf]
    meta = {
        'status': 'enabled',
        'label_mo': True,
        'detection': {'character_table': C2V, 'operations': []},
        'mo_labels': {
            'status': 'ok',
            'alpha': {'labels': labels},
            'beta': {'labels': labels},
        },
        'use_response_symmetry': False,
    }
    data = FakeData({'nelec_A': np.array([na]), 'nelec_B': np.array([nb])})
    cfg = {'tdhf': {'type': 'mrsf', 'multiplicity': multiplicity}}
    return _molecule(cfg, meta, data), na, nb, nbf


class QuintetResponseSpaceTests(unittest.TestCase):
    """Quintet MRSF stages beta-occupied to alpha-virtual."""

    def test_quintet_stages_the_transposed_space(self):
        mol, na, nb, nbf = _staging_molecule(5)
        self.assertTrue(mol.stage_response_symmetry())
        staged = np.asarray(mol.data['OQP::sym_pair_irrep']).size
        self.assertEqual(staged, nb * (nbf - na),
                         'quintet must stage beta-occ x alpha-vir')
        self.assertNotEqual(staged, na * (nbf - nb))

    def test_triplet_keeps_the_original_space(self):
        mol, na, nb, nbf = _staging_molecule(3)
        self.assertTrue(mol.stage_response_symmetry())
        staged = np.asarray(mol.data['OQP::sym_pair_irrep']).size
        self.assertEqual(staged, na * (nbf - nb))

    def test_pair_labels_are_the_products_of_the_right_orbital_sets(self):
        """Not just the length -- the contents must come from beta occ x alpha
        vir, which a transposed-but-same-size bug would get wrong."""
        from oqp.library.symmetry import product_irrep
        mol, na, nb, nbf = _staging_molecule(5)
        mol.stage_response_symmetry()
        staged = np.asarray(mol.data['OQP::sym_pair_irrep']).ravel()
        labels = mol.symmetry_metadata['mo_labels']['alpha']['labels']
        irreps = list(C2V.keys())
        expected = [irreps.index(product_irrep([occ, vir], C2V)) + 1
                    for vir in labels[na:] for occ in labels[:nb]]
        self.assertEqual(staged.tolist(), expected)

    def test_projection_flag_and_status_follow_use_response_symmetry(self):
        mol, *_ = _staging_molecule(3)
        mol.stage_response_symmetry()
        self.assertEqual(int(np.asarray(
            mol.data['OQP::sym_response_project_enable']).ravel()[0]), 0)
        self.assertEqual(
            mol.symmetry_metadata['response_symmetry']['status'],
            'pair_table_staged')

        mol.symmetry_metadata['use_response_symmetry'] = True
        mol.stage_response_symmetry()
        self.assertEqual(int(np.asarray(
            mol.data['OQP::sym_response_project_enable']).ravel()[0]), 1)
        self.assertEqual(
            mol.symmetry_metadata['response_symmetry']['status'], 'active')

    def test_a_bail_clears_both_the_table_and_its_metadata(self):
        """A stale table is worse than none: it is well formed, just wrong."""
        mol, *_ = _staging_molecule(3)
        mol.stage_response_symmetry()
        self.assertIn('OQP::sym_pair_irrep', mol.data)

        del mol.symmetry_metadata['detection']          # next step: no detection
        self.assertFalse(mol.stage_response_symmetry())
        self.assertNotIn('OQP::sym_pair_irrep', mol.data)
        self.assertIsNone(mol.symmetry_metadata.get('response_symmetry'))


class LabelMoIsADisplayFlagTests(unittest.TestCase):
    """label_mo=false must not disable the irrep computation."""

    def _mol(self):
        meta = {
            'status': 'enabled',
            'label_mo': False,
            'detection': {'character_table': C2V, 'operations': []},
        }
        mol = _molecule({'symmetry': {'enabled': 'true', 'label_mo': False}}, meta)
        # Stop the method just past the label_mo gate, so reaching this stub is
        # the observable proof that the gate was bypassed.
        mol._symmetry_labeling_inputs = lambda: (None, None, 0, 'skipped_probe')
        return mol

    def test_display_off_suppresses_an_ordinary_call(self):
        mol = self._mol()
        self.assertIsNone(mol.label_molecular_orbitals())
        self.assertNotIn('mo_labels', mol.symmetry_metadata)

    def test_required_bypasses_the_display_flag(self):
        mol = self._mol()
        self.assertIsNone(mol.label_molecular_orbitals(required=True))
        self.assertEqual(mol.symmetry_metadata['mo_labels']['status'],
                         'skipped_probe')

    def test_required_is_not_the_default(self):
        import inspect
        molecule_module = load_molecule_module()
        sig = inspect.signature(molecule_module.Molecule.label_molecular_orbitals)
        self.assertIs(sig.parameters['required'].default, False)

    def test_staging_gets_labels_even_with_the_display_off(self):
        """The end-to-end consequence: the table is still staged."""
        mol, na, nb, nbf = _staging_molecule(3)
        mol.symmetry_metadata['label_mo'] = False
        self.assertTrue(mol.stage_response_symmetry())
        self.assertEqual(np.asarray(mol.data['OQP::sym_pair_irrep']).size,
                         na * (nbf - nb))


class FullGroupTagsAreClearedTests(unittest.TestCase):
    """A declined or failed full tier must not leave its blocks behind.

    Their mere presence makes the Fortran side take the blocked path, which
    then returns on a size mismatch -- skipping the abelian symmetrization
    entirely, and silently.
    """

    def test_clearing_removes_both_the_tag_and_the_metadata(self):
        meta = {'reduction_maps_full': {'n_operations': 24}}
        data = FakeData({'OQP::sym_op_blocks': np.zeros(4)})
        mol = _molecule({}, meta, data)
        mol._clear_full_group_tags(meta)
        self.assertNotIn('OQP::sym_op_blocks', mol.data)
        self.assertNotIn('reduction_maps_full', meta)

    def test_clearing_is_safe_when_nothing_was_staged(self):
        mol = _molecule({}, {}, FakeData())
        mol._clear_full_group_tags({})  # must not raise


class StaleIntegralStagingTests(unittest.TestCase):
    """A staging attempt that BAILS must not leave the previous one active.

    Every early return in the staging body -- no detection, frame not
    reoriented, no basis, no shell purity, AO-count mismatch, failed overlap
    invariance -- used to return before any invalidation, so the prior step's
    petite maps stayed in the tag store with ``sym_petite_enable`` still 1.
    Fortran would then apply them to the NEW geometry: wrong integrals, not a
    missed optimisation.
    """

    ALL_TAGS = (
        'OQP::sym_shell_map', 'OQP::sym_ao_target', 'OQP::sym_ao_sign',
        'OQP::sym_atom_weight', 'OQP::sym_op_blocks', 'OQP::sym_petite_enable',
    )

    def _staged(self):
        """A molecule carrying a complete, previously-staged petite setup."""
        data = FakeData({t: np.ones(2) for t in self.ALL_TAGS})
        meta = {
            'use_integral_symmetry': True,
            'reduction_maps': {'n_operations': 8},
            'reduction_maps_full': {'n_operations': 24},
            # No 'detection' -> the FIRST gate in the staging body returns.
        }
        return _molecule({'symmetry': {'use_integral_symmetry': 'true'}}, meta, data)

    def test_a_bail_on_missing_detection_clears_everything(self):
        mol = self._staged()
        self.assertFalse(mol.stage_integral_symmetry_maps())
        for tag in self.ALL_TAGS:
            self.assertNotIn(tag, mol.data, f'{tag} survived a failed staging')
        self.assertNotIn('reduction_maps', mol.symmetry_metadata)
        self.assertNotIn('reduction_maps_full', mol.symmetry_metadata)

    def test_turning_integral_symmetry_off_also_clears(self):
        """The flag gate returns even earlier than the staging body."""
        mol = self._staged()
        mol.symmetry_metadata['use_integral_symmetry'] = False
        self.assertFalse(mol.stage_integral_symmetry_maps())
        for tag in self.ALL_TAGS:
            self.assertNotIn(tag, mol.data, f'{tag} survived with the flag off')

    def test_petite_enable_is_deleted_not_left_asserting_active(self):
        """Absence and 0 are equivalent to the consumers; a surviving 1 is not."""
        mol = self._staged()
        mol.stage_integral_symmetry_maps()
        self.assertNotIn('OQP::sym_petite_enable', mol.data)

    def test_active_status_is_dropped_so_metadata_matches_the_tags(self):
        """`status == 'active'` drives _petite_is_staged() and
        symmetrize_gradient(); leaving it projects a gradient with the previous
        geometry's operations after the reduction was turned off."""
        mol = self._staged()
        mol.symmetry_metadata['integral_symmetry'] = {'status': 'active',
                                                      'n_operations': 8}
        mol.symmetry_metadata['use_integral_symmetry'] = False
        mol.stage_integral_symmetry_maps()
        self.assertNotEqual(
            mol.symmetry_metadata.get('integral_symmetry', {}).get('status'),
            'active')

    def test_reoriented_status_survives_invalidation(self):
        """The staging body REQUIRES it; clearing it would disable the
        reduction outright, so only 'active' may be dropped."""
        mol = self._staged()
        mol.symmetry_metadata['integral_symmetry'] = {'status': 'reoriented'}
        mol._clear_integral_symmetry_tags(mol.symmetry_metadata)
        self.assertEqual(
            mol.symmetry_metadata['integral_symmetry']['status'], 'reoriented')


class StaleMoLabelsTests(unittest.TestCase):
    """MO labels describe one geometry and must not be reused across steps.

    Reachable whenever the ordinary post-SCF relabelling does not run -- most
    obviously with label_mo=false, where label_molecular_orbitals() returns
    early and step N's labels persist into step N+1. The staged pair table is
    then well formed and describes the wrong orbitals.
    """

    def _mol(self, labels_key, current_geom):
        meta = {
            'status': 'enabled',
            'label_mo': False,
            'detection': {'character_table': C2V, 'operations': []},
            'mo_labels': {'status': 'ok', 'geometry_key': labels_key,
                          'alpha': {'labels': ['a1', 'b2']}},
        }
        mol = _molecule({'tdhf': {'type': 'mrsf'}}, meta)
        mol.get_system = lambda: current_geom
        return mol

    def test_labels_from_the_same_geometry_are_reused(self):
        mol = self._mol(None, [0.0, 0.0, 0.0])
        key = mol._symmetry_geometry_key()
        mol.symmetry_metadata['mo_labels']['geometry_key'] = key
        self.assertIsNotNone(mol._usable_mo_labels(mol.symmetry_metadata))

    def test_labels_from_a_different_geometry_are_rejected(self):
        mol = self._mol('a-key-from-some-other-structure', [0.0, 0.0, 0.0])
        self.assertIsNone(mol._usable_mo_labels(mol.symmetry_metadata))

    def test_a_moved_atom_changes_the_key(self):
        a = self._mol(None, [0.0, 0.0, 0.0])._symmetry_geometry_key()
        b = self._mol(None, [0.0, 0.0, 0.001])._symmetry_geometry_key()
        self.assertIsNotNone(a)
        self.assertNotEqual(a, b)

    def test_unknown_geometry_keeps_the_cached_labels(self):
        """None means 'cannot tell'. Recomputing unconditionally would pay the
        dense O(|G| n_AO^3) relabelling on every call."""
        mol = self._mol('whatever', [0.0, 0.0, 0.0])
        mol.get_system = lambda: (_ for _ in ()).throw(RuntimeError('no geometry'))
        self.assertIsNone(mol._symmetry_geometry_key())
        self.assertIsNotNone(mol._usable_mo_labels(mol.symmetry_metadata))


class TdMultiplicityKeyTests(unittest.TestCase):
    """The schema key is 'multiplicity'; 'mult' is only an alias."""

    def test_schema_key_is_preferred(self):
        self.assertEqual(_molecule({'tdhf': {'multiplicity': 3}})._td_multiplicity(), 3)

    def test_alias_still_works_for_hand_written_configs(self):
        self.assertEqual(_molecule({'tdhf': {'mult': 3}})._td_multiplicity(), 3)

    def test_schema_key_wins_over_the_alias(self):
        mol = _molecule({'tdhf': {'multiplicity': 5, 'mult': 1}})
        self.assertEqual(mol._td_multiplicity(), 5)

    def test_string_values_from_the_input_parser_are_accepted(self):
        self.assertEqual(_molecule({'tdhf': {'multiplicity': '5'}})._td_multiplicity(), 5)

    def test_absent_or_unparsable_is_none_not_an_exception(self):
        # This feeds a label formatter that must never abort a converged run.
        self.assertIsNone(_molecule({})._td_multiplicity())
        self.assertIsNone(_molecule({'tdhf': {}})._td_multiplicity())
        self.assertIsNone(_molecule({'tdhf': {'multiplicity': ''}})._td_multiplicity())
        self.assertIsNone(_molecule({'tdhf': {'multiplicity': 'x'}})._td_multiplicity())

    def test_the_schema_default_is_reachable_through_the_helper(self):
        """Guards the actual defect: reading only 'mult' made every validated
        config fall through, so `terms` was never built."""
        mol = _molecule({'tdhf': {'type': 'mrsf', 'multiplicity': 3}})
        self.assertEqual(f"{mol._td_multiplicity()}A1", '3A1')


class StaleModeLabelsTests(unittest.TestCase):
    """Mode labels must never outlive the modes they describe."""

    def test_labels_are_dropped_when_labelling_is_skipped(self):
        meta = {
            'status': 'enabled',
            'label_modes': False,
            'mode_labels': {'status': 'ok', 'labels': ['a1', 'a1', 'b2']},
        }
        mol = _molecule({}, meta)
        self.assertIsNone(mol.label_normal_modes())
        self.assertNotIn('mode_labels', meta)

    def test_labels_are_dropped_when_detection_is_unavailable(self):
        meta = {
            'status': 'enabled',
            'label_modes': True,
            'mode_labels': {'status': 'ok', 'labels': ['a1', 'a1', 'b2']},
        }
        mol = _molecule({}, meta)
        self.assertIsNone(mol.label_normal_modes())
        self.assertNotIn('mode_labels', meta)


if __name__ == '__main__':
    unittest.main()
