"""Reorientation must not survive a path that does not end with the reduction.

`reorient_for_integral_symmetry` moves the molecule into the symmetry standard
frame before the basis exists. Staging the petite maps can then still decline,
for reasons only visible once the basis is built. Those runs were left in the
rotated frame with nothing to show for it.

Two things made that worse than a frame convention:

  - The standard frame is centred on the charge-weighted centroid, so every
    molecule is also TRANSLATED -- including C1, where the rotation is exactly
    the identity and there is no symmetry to exploit at all.
  - `input_coords` was captured in the reorientation loop and never used, so
    both failure exits returned False with the geometry already moved.
"""

import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


class FailureExitsRestoreTheGeometry(unittest.TestCase):
    """The captured coordinates must actually be used on both exits."""

    def _reorient_source(self):
        source = (ROOT / 'pyoqp/oqp/molecule/molecule.py').read_text()
        start = source.index('def reorient_for_integral_symmetry')
        end = source.index('def stage_integral_symmetry_maps')
        return source[start:end]

    def test_the_non_converged_exit_restores(self):
        body = self._reorient_source()
        bail = body.index("'skipped_orientation_not_converged'")
        window = body[:bail]
        self.assertIn('self.update_system(input_coords.ravel())', window,
                      'the non-converged exit must put the molecule back')

    def test_the_exception_handler_restores(self):
        body = self._reorient_source()
        handler = body.rindex('except Exception as exc:')
        self.assertIn('self.update_system(input_coords.ravel())',
                      body[handler:],
                      'the error exit must put the molecule back')

    def test_the_coordinates_survive_a_staging_bail(self):
        """Every bail in staging replaces meta['integral_symmetry'] wholesale.

        Keeping the coordinates inside that block silently loses them -- which
        is exactly what happened first time: the restore never fired for any
        spherical basis, because skipped_basis_mismatch had already wiped it.
        """
        body = self._reorient_source()
        self.assertIn("meta['_reorient_input_coords']", body)
        block = body[body.index("'status': 'reoriented'"):]
        self.assertNotIn("'input_coords'", block[:block.index('return True')],
                         'must live outside the block staging overwrites')

    def test_the_caller_pops_the_key_so_it_never_reaches_save_data(self):
        source = (ROOT / 'pyoqp/oqp/library/single_point.py').read_text()
        self.assertIn("meta.pop('_reorient_input_coords', None)", source)

    def test_the_caller_refreshes_the_guess_after_restoring(self):
        """The 1e integrals were built in the rotated frame."""
        source = (ROOT / 'pyoqp/oqp/library/single_point.py').read_text()
        restore = source.index("meta.pop('_reorient_input_coords', None)")
        window = source[restore:restore + 900]
        self.assertIn('self.mol.update_system(', window)
        self.assertIn('self._prep_guess()', window)
        self.assertLess(window.index('self.mol.update_system('),
                        window.index('self._prep_guess()'))


class PropDoesNotGetReoriented(unittest.TestCase):
    """`prop` consumes an externally supplied geometry in the caller's frame.

    load_previous_data writes OQP::xyz_old from it, and get_basis_overlap
    overlaps that against the current -- by then rotated and translated --
    structure. The cross-geometry MO overlap, the phase/reorder alignment built
    on it, and NACME are then all wrong. That is PR #319's finite-difference
    failure mode one level up.
    """

    def test_the_allow_list_is_energy_and_grad_only(self):
        source = (ROOT / 'pyoqp/oqp/molecule/molecule.py').read_text()
        self.assertIn("if runtype not in ('energy', 'grad'):", source)

    def test_prop_is_not_admitted(self):
        source = (ROOT / 'pyoqp/oqp/molecule/molecule.py').read_text()
        allow = source.index('if runtype not in (')
        line = source[allow:source.index('\n', allow)]
        self.assertNotIn("'prop'", line)
        self.assertNotIn("'properties'", line)


class ReorientationIsATranslationToo(unittest.TestCase):
    """Not a rotation of symmetric molecules only -- everything moves."""

    def test_a_c1_geometry_is_still_translated(self):
        import importlib.util
        import sys

        spec = importlib.util.spec_from_file_location(
            'symmetry_detect_frame_under_test',
            ROOT / 'pyoqp/oqp/library/symmetry_detect.py')
        detect = importlib.util.module_from_spec(spec)
        sys.modules['symmetry_detect_frame_under_test'] = detect
        spec.loader.exec_module(detect)

        # A genuine C1 case, deliberately off-origin. Five distinct elements,
        # non-planar -- note three atoms are always coplanar, so anything
        # smaller is at least Cs and would not test what this claims to.
        geometry = np.array([[1.0, 2.0, 3.0],
                             [1.6, 2.7, 3.4],
                             [2.1, 1.3, 2.6],
                             [0.4, 1.4, 4.1],
                             [0.6, 3.0, 2.2]])
        result = detect.detect_point_group([6, 1, 9, 17, 35], geometry,
                                           tolerance=1e-5)

        self.assertEqual(result['point_group'], 'c1')
        rotation = np.asarray(result['orientation'], dtype=float)
        origin = np.asarray(result['origin'], dtype=float)
        self.assertLess(float(np.max(np.abs(rotation - np.eye(3)))), 1e-9,
                        'C1 should not be rotated')
        # ~2.8 bohr for this geometry.
        self.assertGreater(float(np.max(np.abs(origin))), 1.0,
                           'but it IS translated -- which is why the restore '
                           'matters even with no symmetry to exploit')


if __name__ == '__main__':
    unittest.main()
