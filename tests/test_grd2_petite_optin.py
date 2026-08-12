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
    """A staged reduction must not survive into the next job on the same
    Molecule: the enable flag and the maps live in the tag store, which
    outlives one `reference()` call.

    This asserted the literal text ``self._set_petite_enabled(False)`` followed
    by ``if symmetry_on:``, and its position relative to the reorientation
    call. That went red the moment the invalidation was strengthened -- main
    replaced the flag reset with ``clear_integral_symmetry_state()``, which
    DELETES ``OQP::sym_petite_enable`` along with the maps instead of zeroing
    it, so absence and 0 cannot disagree. The stronger implementation failed a
    test written for the weaker one, which is a statement about source text
    rather than about the invariant.

    It now drives the real `reference()` with the collaborators stubbed and
    records the call order, so it tracks the invariant it is named for and is
    indifferent to which mechanism enforces it.
    """

    def test_density_guard_persists_the_fallback_for_xc_and_gradients(self):
        """A local JK fallback must also disable later reduced consumers."""
        source = ROOT.joinpath('source/scf_addons.F90').read_text()
        fock_start = source.index('subroutine fock_jk(')
        routec = source.index('routec_try_fock_jk(', fock_start)
        guard_start = source.index('if (present(petite)) then', fock_start)
        guard = source[guard_start:routec]
        self.assertLess(guard_start, routec,
                        'density guard must run before Route-C can return')
        self.assertIn(
            'call tagarray_get_data(infos%dat, OQP_sym_petite,', guard)
        self.assertIn('global_petite(1) = 0_8', guard)

    def test_dense_input_frame_operators_keep_conservative_screening(self):
        source = ROOT.joinpath('source/integrals/int2.F90').read_text()
        self.assertIn('logical :: sym_dense = .false.', source)
        self.assertIn(
            'if (status == TA_OK) this%sym_dense = size(blocks) > 0',
            source,
        )
        self.assertIn(
            'if (this%petite .and. this%sym_dense) '
            'test = test*real(this%sym_nops, dp)',
            source,
        )
        self.assertIn('eri_data%weighted_cutoff = this%sym_dense', source)

    class _Stop(Exception):
        pass

    def test_the_stale_reduction_is_dropped_before_anything_can_consume_it(self):
        import types

        from oqp.library import single_point as sp

        calls = []
        obj = sp.SinglePoint.__new__(sp.SinglePoint)
        obj.method = 'hf'
        obj.init_scf = 'no'
        obj._init_convergence = lambda: calls.append('init_convergence')
        obj._prep_guess = lambda: calls.append('prep_guess')
        obj.swapmo = lambda: calls.append('swapmo')
        obj._set_petite_enabled = lambda flag: calls.append(f'petite_enabled({flag})')

        def run_scf():
            calls.append('run_scf')
            raise self._Stop()

        obj._run_scf = run_scf

        mol = types.SimpleNamespace()
        # A molecule carrying a reduction staged by a PREVIOUS job.
        mol.symmetry_metadata = {'status': 'enabled',
                                 'use_integral_symmetry': True}
        mol.has_staged_integral_symmetry = lambda: True
        mol.clear_integral_symmetry_state = lambda: calls.append('clear_staged')
        mol.reorient_for_integral_symmetry = lambda: calls.append('reorient')
        mol.stage_integral_symmetry_maps = lambda: calls.append('stage')
        mol.update_system = lambda x: calls.append('restore_geometry')
        mol._detect_symmetry_metadata = lambda: calls.append('redetect')
        obj.mol = mol

        original = sp.dump_log
        sp.dump_log = lambda *a, **k: None
        try:
            obj.reference()
        except self._Stop:
            pass
        finally:
            sp.dump_log = original

        # Whichever mechanism is used, the previous job's staged state has to
        # be dropped, and dropped before reorientation -- otherwise a job whose
        # maps are never staged runs with a stale enable against a shell map
        # built for another geometry.
        dropped = [c for c in calls
                   if c == 'clear_staged' or c == 'petite_enabled(False)']
        self.assertTrue(dropped,
                        f'nothing invalidated the staged reduction: {calls}')
        self.assertLess(
            calls.index(dropped[0]), calls.index('reorient'),
            f'invalidation must precede reorientation, got {calls}')
        self.assertLess(
            calls.index(dropped[0]), calls.index('prep_guess'),
            f'invalidation must precede the guess/basis stage, got {calls}')


if __name__ == '__main__':
    unittest.main()
