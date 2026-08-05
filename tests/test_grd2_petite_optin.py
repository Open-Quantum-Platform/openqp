"""grd2 must require a per-caller opt-in for the symmetry petite reduction.

int2 has always required one (`enable_petite`, used by exactly one caller --
the SCF Fock build). grd2 loaded the shell map unconditionally, so every
grd2_driver caller inherited the reduction from the global tag alone.

That is only valid when the density being contracted is totally symmetric.
It is not for the CPHF probe densities routed through fock_deriv_contract --
a single occupied-virtual outer product -- nor for hf_hessian's resp_grad,
which additionally runs at a displaced geometry while the staged shell map
describes the reference frame.

The defect is latent today only because runtype 'hess' fails the allow-list in
Molecule.reorient_for_integral_symmetry. Measured with that allow-list
temporarily widened to admit 'hess', water/6-31G HF analytical Hessian:

    use_integral_symmetry=false   1829.11   3906.40   4001.41 cm^-1
    use_integral_symmetry=true  -44999.00 -35788.33 -32424.62   (before)
    use_integral_symmetry=true    1829.11   3906.40   4001.41   (after)

So this must land before anything widens that allow-list.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# The opt-in set is exactly the set whose skeleton output is projected
# afterwards by Molecule.symmetrize_gradient.
OPTED_IN = (
    'source/modules/hf_gradient.F90',
    'source/modules/tdhf_gradient.F90',
    'source/modules/tdhf_sf_gradient.F90',
    'source/modules/tdhf_mrsf_gradient.F90',
)

# Callers that contract a density which is NOT totally symmetric.
NOT_OPTED_IN = (
    'source/modules/fock_deriv.F90',
    'source/modules/hf_hessian.F90',
)


def call_arguments(path, callee):
    """Text of every `call <callee>(...)` in a file, continuations joined."""
    source = ROOT.joinpath(path).read_text()
    joined = re.sub(r'&\s*\n\s*', ' ', source)
    return re.findall(r'call\s+%s\s*\((?:[^()]|\([^()]*\))*\)' % callee,
                      joined, flags=re.IGNORECASE)


class Grd2RequiresAnOptIn(unittest.TestCase):
    def test_the_shell_map_load_is_gated(self):
        source = ROOT.joinpath('source/integrals/grd2.F90').read_text()
        self.assertRegex(
            source,
            r'if \(present\(petite\)\) then\s*\n\s*'
            r'if \(petite\) call load_petite_shell_map',
            'load_petite_shell_map must be reached only through the opt-in')

    def test_locals_are_initialised_on_the_off_path(self):
        """load_petite_shell_map used to be what nulled these.

        sym_nops is read at the quartet loop and sym_map dereferenced inside
        petite_quartet_weight, so gating the call without initialising them
        would read uninitialised memory -- an intermittent wrong gradient.
        """
        source = ROOT.joinpath('source/integrals/grd2.F90').read_text()
        gate = source.index('if (present(petite)) then')
        preamble = source[:gate]
        self.assertIn('sym_map => null()', preamble)
        self.assertIn('sym_nops = 0', preamble)

    def test_both_cam_passes_forward_the_opt_in(self):
        """grd2_driver splits into two grd2_driver_gen calls for CAM/LC.

        Dropping the argument on one pass would give a range-separated
        excited-state gradient the reduction on only half of itself.
        """
        calls = call_arguments('source/integrals/grd2.F90', 'grd2_driver_gen')
        self.assertEqual(len(calls), 3, 'expected the CAM two-pass split plus '
                                        'the plain path')
        for call in calls:
            self.assertIn('petite', call, call)


class OnlySymmetricDensitiesOptIn(unittest.TestCase):
    def test_the_four_gradient_sites_opt_in(self):
        for path in OPTED_IN:
            with self.subTest(path=path):
                calls = call_arguments(path, 'grd2_driver')
                self.assertEqual(len(calls), 1, path)
                self.assertRegex(calls[0], r'petite\s*=\s*\.true\.', path)

    def test_cphf_and_displaced_geometry_callers_do_not(self):
        for path in NOT_OPTED_IN:
            with self.subTest(path=path):
                for call in call_arguments(path, 'grd2_driver'):
                    self.assertNotIn('petite', call,
                                     f'{path}: this density is not totally '
                                     f'symmetric, it must not opt in')

    def test_no_caller_outside_these_files_opts_in(self):
        """A new caller must think about its density, not inherit a default."""
        opted = []
        for path in ROOT.joinpath('source').rglob('*.F90'):
            rel = str(path.relative_to(ROOT))
            if rel == 'source/integrals/grd2.F90':
                continue
            for call in call_arguments(rel, 'grd2_driver'):
                if 'petite' in call:
                    opted.append(rel)
        self.assertCountEqual(opted, list(OPTED_IN))


class PetiteStateIsResetPerReference(unittest.TestCase):
    def test_reference_starts_from_the_reduction_off(self):
        """The enable flag lives in the tag store and survives across jobs."""
        source = ROOT.joinpath('pyoqp/oqp/library/single_point.py').read_text()
        reset = source.index('self._set_petite_enabled(False)\n        if symmetry_on:')
        reorient = source.index('self.mol.reorient_for_integral_symmetry()')
        self.assertLess(reset, reorient,
                        'the reset must precede reorientation, so a stale '
                        'enable cannot survive into a job whose maps are '
                        'never staged')


if __name__ == '__main__':
    unittest.main()
