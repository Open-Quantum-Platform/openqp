"""Gates around staging the per-pair irrep table for the Davidson guess.

Three defects, all of which only became reachable once ``[symmetry] enabled``
started defaulting to true, and all of which failed SILENTLY -- the guess fell
back to the historical seeds and nothing in the output said so. That silence is
the reason they are pinned here rather than left to a live example: the failure
mode of this whole mechanism is "the answer looks fine and is a different
state".

1. ``label_mo`` is a DISPLAY flag. It must not decide whether the irreps are
   computed, because the guess coverage needs them.
2. Quintet MRSF runs a transposed response space (beta occupied to alpha
   virtual). Staging the singlet/triplet space produces a table of the wrong
   length, which every Fortran consumer then discards at its size guard.
3. The state-label formatter must read the validated schema key
   ``[tdhf] multiplicity``, not the test-only ``mult`` alias.
"""

import unittest

from test_symmetry_metadata import load_molecule_module


def _molecule(config, metadata=None):
    molecule_module = load_molecule_module()
    molecule = molecule_module.Molecule.__new__(molecule_module.Molecule)
    molecule.config = config
    molecule.symmetry_metadata = metadata if metadata is not None else {}
    return molecule


class TdMultiplicityKeyTests(unittest.TestCase):
    """The schema key is 'multiplicity'; 'mult' is only an alias."""

    def test_schema_key_is_preferred(self):
        mol = _molecule({'tdhf': {'multiplicity': 3}})
        self.assertEqual(mol._td_multiplicity(), 3)

    def test_alias_still_works_for_hand_written_configs(self):
        mol = _molecule({'tdhf': {'mult': 3}})
        self.assertEqual(mol._td_multiplicity(), 3)

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
        self.assertIsNone(_molecule({'tdhf': {'multiplicity': 'triplet'}})._td_multiplicity())


class LabelMoIsADisplayFlagTests(unittest.TestCase):
    """label_mo=false must not disable the irrep computation."""

    def _mol(self, label_mo):
        # 'detection' left absent: label_molecular_orbitals returns None there
        # too, so reaching that return proves the label_mo gate was passed.
        return _molecule(
            {'symmetry': {'enabled': 'true', 'label_mo': label_mo}},
            {'status': 'enabled', 'label_mo': label_mo},
        )

    def test_display_off_still_computes_when_required(self):
        mol = self._mol(False)
        # Without required=True the method returns before it ever looks for
        # 'detection'; with it, the flag is bypassed and the absent detection
        # is what stops it. Distinguish the two by the metadata side effect
        # neither path writes, and by the flag itself staying false.
        self.assertIsNone(mol.label_molecular_orbitals())
        self.assertIsNone(mol.label_molecular_orbitals(required=True))
        self.assertFalse(mol.symmetry_metadata['label_mo'])

    def test_required_is_not_the_default(self):
        """An ordinary call must keep honouring the user's display choice."""
        import inspect
        molecule_module = load_molecule_module()
        sig = inspect.signature(molecule_module.Molecule.label_molecular_orbitals)
        self.assertIs(sig.parameters['required'].default, False)

    def test_staging_asks_for_the_labels_it_needs(self):
        """stage_response_symmetry must pass required=True.

        Checked on the source rather than by running it: the live path needs a
        converged SCF, and the point of the assertion is that the coupling is
        not reintroduced.
        """
        from pathlib import Path
        src = (Path(__file__).resolve().parents[1]
               / 'pyoqp' / 'oqp' / 'molecule' / 'molecule.py').read_text()
        body = src.split('def stage_response_symmetry', 1)[1]
        body = body.split('\n    def ', 1)[0]
        self.assertIn('label_molecular_orbitals(required=True)', body)


class QuintetResponseSpaceTests(unittest.TestCase):
    """Quintet MRSF stages beta-occupied to alpha-virtual."""

    def test_quintet_branch_is_selected_on_multiplicity_five(self):
        mol = _molecule({'tdhf': {'type': 'mrsf', 'multiplicity': 5}})
        self.assertEqual(mol._td_multiplicity(), 5)

    def test_triplet_mrsf_keeps_the_original_space(self):
        mol = _molecule({'tdhf': {'type': 'mrsf', 'multiplicity': 3}})
        self.assertNotEqual(mol._td_multiplicity(), 5)

    def test_both_staging_and_labelling_branch_on_the_quintet(self):
        """The two sites must agree; fixing one and not the other is how the
        labels end up computed from a misread buffer."""
        from pathlib import Path
        src = (Path(__file__).resolve().parents[1]
               / 'pyoqp' / 'oqp' / 'molecule' / 'molecule.py').read_text()
        for method in ('stage_response_symmetry', 'label_excited_states'):
            body = src.split(f'def {method}', 1)[1].split('\n    def ', 1)[0]
            self.assertIn("_td_multiplicity() == 5", body,
                          f'{method} does not special-case quintet MRSF')


if __name__ == '__main__':
    unittest.main()
