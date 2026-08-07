"""Child jobs must not inherit symmetry state from a parent or a guess file.

Finite-difference drivers (numerical Hessian, IR/Raman properties, NACME) run
each displaced geometry as its own job and assemble the results in the parent's
frame.  If a child reorients into the standard frame of its own geometry, the
gradient it returns is expressed in a different frame from the one the parent
files it into -- and a displaced geometry generally has lower symmetry, hence a
different standard frame, than the reference.

Water/6-31G HF numerical frequencies measured with the reduction enabled on the
parent, in cm^-1:

    before   -2675.60    697.64   2728.85     (spurious imaginary mode)
    after     1828.86   3906.50   4001.54     (= the symmetry-off result)

Two independent leaks produced that, so both are pinned here.
"""

import unittest

import numpy as np

from oqp.molecule.molecule import Molecule


class ChildConfigDropsIntegralSymmetry(unittest.TestCase):
    """Leak 1: children deep-copy the parent config, including [symmetry]."""

    def test_helper_forces_the_flag_off(self):
        from oqp.library.single_point import _no_integral_symmetry_in_child

        config = {'input': {'runtype': 'grad'},
                  'symmetry': {'enabled': 'true',
                               'use_integral_symmetry': 'true'}}
        _no_integral_symmetry_in_child(config)
        self.assertEqual(config['symmetry']['use_integral_symmetry'], 'false')
        # Untouched: labelling is metadata-only and costs nothing.
        self.assertEqual(config['symmetry']['enabled'], 'true')

    def test_helper_creates_the_section_when_absent(self):
        from oqp.library.single_point import _no_integral_symmetry_in_child

        config = {'input': {'runtype': 'grad'}}
        _no_integral_symmetry_in_child(config)
        self.assertEqual(config['symmetry']['use_integral_symmetry'], 'false')

    def test_every_finite_difference_child_is_covered(self):
        import inspect

        from oqp.library import single_point

        source = inspect.getsource(single_point)
        # One call per child-config site: vibrational properties, the numerical
        # Hessian gradient child, and the NACME child.
        self.assertGreaterEqual(
            source.count('_no_integral_symmetry_in_child(config)'), 3)


class GuessFileDoesNotOverrideSymmetryConfig(unittest.TestCase):
    """Leak 2: put_data replaced this job's metadata with the producer's.

    This is the one that actually bit: the child's own input file said
    ``use_integral_symmetry=false`` and it reoriented anyway, because the
    parent's guess JSON carried ``True`` and won.
    """

    def _molecule(self, use_integral_symmetry):
        mol = Molecule.__new__(Molecule)
        mol.symmetry_metadata = {
            'status': 'enabled',
            'enabled': True,
            'use_integral_symmetry': use_integral_symmetry,
            'use_response_symmetry': False,
            'tolerance': 1.0e-5,
            'point_group': 'c1',
        }
        mol._system = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.4])
        mol.get_system = lambda: mol._system
        return mol

    def test_config_switches_come_from_this_job(self):
        mol = self._molecule(use_integral_symmetry=False)
        merged = mol._merge_restored_symmetry_metadata(
            {'use_integral_symmetry': True, 'use_response_symmetry': True,
             'enabled': True, 'status': 'enabled'},
            stored_coord=mol._system)
        self.assertFalse(merged['use_integral_symmetry'])
        self.assertFalse(merged['use_response_symmetry'])

    def test_geometry_derived_entries_are_dropped_at_a_new_geometry(self):
        mol = self._molecule(use_integral_symmetry=False)
        restored = {
            'point_group': 'c2v',
            'detection': {'point_group': 'c2v', 'operations': [1, 2, 3, 4]},
            'integral_symmetry': {'status': 'reoriented'},
            'mo_labels': {'alpha': ['a1']},
        }
        displaced = mol._system + np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.005])
        merged = mol._merge_restored_symmetry_metadata(
            restored, stored_coord=displaced)
        for key in ('detection', 'integral_symmetry', 'mo_labels'):
            self.assertNotIn(key, merged,
                             f'{key} describes the producer geometry')

    def test_geometry_derived_entries_survive_at_the_same_geometry(self):
        mol = self._molecule(use_integral_symmetry=False)
        restored = {'detection': {'point_group': 'c2v'},
                    'point_group': 'c2v'}
        merged = mol._merge_restored_symmetry_metadata(
            restored, stored_coord=mol._system)
        self.assertEqual(merged['detection']['point_group'], 'c2v')

    def test_missing_coordinates_are_treated_as_a_different_geometry(self):
        mol = self._molecule(use_integral_symmetry=False)
        merged = mol._merge_restored_symmetry_metadata(
            {'detection': {'point_group': 'c2v'}}, stored_coord=None)
        self.assertNotIn('detection', merged)


if __name__ == '__main__':
    unittest.main()
