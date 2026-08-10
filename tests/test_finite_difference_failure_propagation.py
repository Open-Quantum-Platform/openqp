"""A detected finite-difference failure must reach every rank, and a scratch
gradient must only be reused by the run that produced it.

These pin the three defects the bot review found in the first round of fixes on
this branch. Each one is written so that reverting the fix makes it fail: the
assertions are on observed behaviour, never on the presence of a line of source
text. (The earlier files on this branch assert on source text, and a mutation
that inverted the reorientation-restore guard passed all of them.)
"""

import os
import types
import unittest

import numpy as np

from oqp.library.single_point import _displacement_run_signature


class ScratchGradientsBelongToOneRun(unittest.TestCase):
    """`[hess] restart=true` must not reuse another run's gradients.

    The tag used to be the displaced coordinate index and sign only, so the
    same deck rerun at a PERTURBED geometry produced bit-identical scratch
    names. `grad_wrapper` matches by name and never compares coordinates, so
    every old gradient came back as 'loaded'. Water/6-31G HF, populating the
    scratch at one geometry and rerunning at another, in cm^-1:

        reused scratch   1798.99   4047.08   4141.15
        actual answer    1384.73   4886.70   5092.94

    -- exit status 0, no warning.
    """

    GEOM = [0.0, 0.0, -0.0776, -1.0076, 1.0076, -1.1612, 1.0076, -1.0076, -1.1612]
    ATOMS = [8.0, 1.0, 1.0]

    @staticmethod
    def MODEL(state=0, **model):
        """The shape Molecule._hessian_request_signature returns."""
        cfg = {'input': {'basis': '6-31g', 'method': 'hf'}, 'scf': {'type': 'rhf'}}
        cfg['input'].update(model)
        return {'version': 2, 'state': int(state), 'model_config': cfg}

    def signature(self, **over):
        kw = dict(origin_coord=self.GEOM, atoms=self.ATOMS, dx=0.01,
                  model_signature=self.MODEL())
        kw.update(over)
        return _displacement_run_signature(**kw)

    def test_the_same_run_keeps_the_same_signature(self):
        """Restart has to still work -- an unchanged rerun must match."""
        self.assertEqual(self.signature(), self.signature())

    def test_a_perturbed_geometry_changes_the_signature(self):
        moved = list(self.GEOM)
        moved[2] += 0.05  # the displacement that produced the wrong Hessian
        self.assertNotEqual(self.signature(), self.signature(origin_coord=moved))

    def test_a_geometry_change_far_below_the_step_size_still_changes_it(self):
        """A perturbation smaller than dx is exactly the dangerous case: the
        orbit representatives are unchanged, so every tag used to collide."""
        moved = list(self.GEOM)
        moved[2] += 1.0e-9
        self.assertNotEqual(self.signature(), self.signature(origin_coord=moved))

    def test_a_different_step_size_changes_the_signature(self):
        self.assertNotEqual(self.signature(), self.signature(dx=0.02))

    def test_a_different_target_state_changes_the_signature(self):
        self.assertNotEqual(self.signature(),
                            self.signature(model_signature=self.MODEL(state=1)))

    def test_a_different_molecule_changes_the_signature(self):
        self.assertNotEqual(self.signature(),
                            self.signature(atoms=[8.0, 1.0, 9.0]))

    def test_a_different_hamiltonian_changes_the_signature(self):
        """Geometry is not the only axis. Change only the functional, keep the
        directory and geometry, rerun with restart=true: without the model in
        the signature every tag still matched and the driver assembled a
        Hessian from the old Hamiltonian's gradients."""
        self.assertNotEqual(
            self.signature(),
            self.signature(model_signature=self.MODEL(functional='bhhlyp')))

    def test_a_different_basis_changes_the_signature(self):
        self.assertNotEqual(
            self.signature(),
            self.signature(model_signature=self.MODEL(basis='cc-pvdz')))


class ReorientationRollbackLeavesADetectionThatDescribesUs(unittest.TestCase):
    """After rolling back the standard-frame move, the stored detection must
    describe the geometry the molecule is actually in.

    `reorient_for_integral_symmetry` calls `attach_detection_metadata` once per
    attempt, so `meta['detection']` tracks the last ROTATED frame. On the
    non-converged exit the coordinates go back but the detection did not, and
    the caller's restore block cannot repair it -- that block keys off
    `_reorient_input_coords`, which only exists on the success path. The
    MO/state/mode labellers then run with operations for a frame the molecule
    has left.
    """

    # Water, C2v, rotated 45 deg about z so the standard frame is genuinely
    # different from the input frame and the rollback has something to undo.
    COORD = np.array([0.0, 0.0, -0.0776,
                      -1.0076, 1.0076, -1.1612,
                      1.0076, -1.0076, -1.1612])
    ATOMS = [8.0, 1.0, 1.0]

    def _molecule(self):
        from oqp.molecule.molecule import Molecule
        mol = Molecule.__new__(Molecule)
        mol._system = self.COORD.copy()
        mol.config = {'input': {'runtype': 'grad'}}
        mol.get_system = lambda: mol._system
        mol.get_atoms = lambda: np.array(self.ATOMS)
        mol.update_system = lambda x: setattr(
            mol, '_system', np.asarray(x, dtype=float).copy())
        mol.symmetry_metadata = {
            'status': 'enabled', 'enabled': True, 'use_integral_symmetry': True,
            'tolerance': 1e-5, 'strict': False,
            'requested_point_group': 'auto', 'requested_subgroup': 'auto',
        }
        return mol

    def test_a_non_converged_reorientation_redetects(self):
        import oqp.library.symmetry_detect as sd

        real = sd.attach_detection_metadata

        def never_converges(meta, atoms, coords):
            """Real detection, but never reports the frame as already standard,
            so all four attempts are used up and the rollback fires."""
            real(meta, atoms, coords)
            rot = np.asarray(meta['detection']['orientation'], dtype=float)
            if np.max(np.abs(rot - np.eye(3))) < 1e-12:
                meta['detection']['orientation'] = [[0.0, 1.0, 0.0],
                                                    [-1.0, 0.0, 0.0],
                                                    [0.0, 0.0, 1.0]]
            return meta

        mol = self._molecule()
        real(mol.symmetry_metadata, np.array(self.ATOMS), self.COORD.reshape(-1, 3))

        sd.attach_detection_metadata = never_converges
        try:
            staged = mol.reorient_for_integral_symmetry()
        finally:
            sd.attach_detection_metadata = real

        self.assertFalse(staged)
        self.assertEqual(
            mol.symmetry_metadata['integral_symmetry']['status'],
            'skipped_orientation_not_converged')

        left = np.asarray(mol.get_system(), dtype=float)
        np.testing.assert_array_equal(
            left, self.COORD, 'the rollback must restore the geometry exactly')

        # The real assertion: detect fresh on the geometry we were left in and
        # require the stored detection to agree with it.
        fresh = {}
        real(fresh, np.array(self.ATOMS), left.reshape(-1, 3))
        np.testing.assert_allclose(
            np.asarray(mol.symmetry_metadata['detection']['orientation'], float),
            np.asarray(fresh['detection']['orientation'], float),
            atol=1e-12,
            err_msg='detection still describes the frame that was rolled back')


class GradWrapperReportsAFailedChild(unittest.TestCase):
    """`grad_wrapper` must report 'failed' whenever the child did not produce
    a fresh gradient for THIS run -- not merely when no file exists.

    Failure is signalled by the read raising. A gradient left by an earlier run
    therefore made a failed child look successful: with a previous run's
    scratch present, a numerical Hessian in which every one of 18 displacements
    failed still exited 0 and printed a full set of frequencies assembled from
    the old files.

    These drive the real `grad_wrapper`, with the child execution stubbed out so
    the test costs nothing. Reverting either the pre-launch removal or the
    widened exception handler turns them red.
    """

    def _key_dict(self, tmp, restart=False):
        return {
            'idx': 0, 'tag': 'sig0c0p',
            'atoms': np.array([8.0, 1.0, 1.0]),
            'coord': np.zeros(9),
            'dir_hess': tmp,
            'project_name': 'probe',
            'config': {'input': {}, 'guess': {}, 'properties': {},
                       'hess': {'state': 0, 'temperature': [298.15]},
                       'tests': {}},
            'guess_file': os.path.join(tmp, 'probe.json'),
            'state': 0, 'restart': restart,
        }

    def _run(self, tmp, restart=False):
        """Call grad_wrapper with the child neutered, as a failed child would
        leave things."""
        from oqp.library import single_point

        os.environ.setdefault('OMP_NUM_THREADS', '1')
        original = single_point._run_oqp_external
        single_point._run_oqp_external = lambda *a, **k: None
        try:
            return single_point.grad_wrapper(self._key_dict(tmp, restart))
        finally:
            single_point._run_oqp_external = original

    def test_a_child_that_produced_nothing_is_failed(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            _, _, status, _ = self._run(tmp)
            self.assertEqual(status, 'failed')

    def test_a_gradient_from_an_earlier_run_does_not_mask_the_failure(self):
        """The reproduced defect: the old file is readable, so without the
        pre-launch removal the worker reports 'computed' and the driver
        assembles a Hessian from stale gradients."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            stale = os.path.join(tmp, 'probe.sig0c0p.grad_0')
            np.savetxt(stale, np.arange(9, dtype=float))
            self.assertTrue(os.path.exists(stale))

            _, grad, status, _ = self._run(tmp)

            self.assertEqual(status, 'failed',
                             'a leftover gradient must not read as success')
            np.testing.assert_array_equal(
                grad, np.zeros(9),
                'the stale values must not reach the Hessian')

    def test_a_truncated_gradient_is_failed_not_reshaped(self):
        """A child killed mid-write leaves a short file. np.loadtxt parses it
        happily, so only a size check catches it -- and the handler has to
        cover ValueError, which FileNotFoundError alone did not."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            key = self._key_dict(tmp, restart=True)
            dat = os.path.join(tmp, 'probe.sig0c0p.grad_0')
            np.savetxt(dat, np.arange(4, dtype=float))  # 4 values, 9 expected

            from oqp.library import single_point
            os.environ.setdefault('OMP_NUM_THREADS', '1')
            _, grad, status, _ = single_point.grad_wrapper(key)

            self.assertEqual(status, 'failed')
            np.testing.assert_array_equal(grad, np.zeros(9))

    def test_a_complete_gradient_still_loads_on_restart(self):
        """The guard must not break the feature it protects."""
        import tempfile

        from oqp.library import single_point

        with tempfile.TemporaryDirectory() as tmp:
            dat = os.path.join(tmp, 'probe.sig0c0p.grad_0')
            expected = np.arange(9, dtype=float)
            np.savetxt(dat, expected)

            os.environ.setdefault('OMP_NUM_THREADS', '1')
            _, grad, status, _ = single_point.grad_wrapper(
                self._key_dict(tmp, restart=True))

            self.assertEqual(status, 'loaded')
            np.testing.assert_allclose(grad, expected)


class TheScratchTagCarriesTheSignature(unittest.TestCase):
    """The signature has to reach the filenames, not just exist as a helper.

    Asserting only on `_displacement_run_signature` would pass with the tag
    integration reverted -- the failure mode this branch's earlier tests had.
    """

    def test_the_tag_changes_when_the_signature_changes(self):
        from oqp.library.single_point import _displacement_run_signature

        geom = [0.0, 0.0, -0.0776, -1.0076, 1.0076, -1.1612,
                1.0076, -1.0076, -1.1612]
        atoms = [8.0, 1.0, 1.0]
        model = {'version': 2, 'state': 0,
                 'model_config': {'input': {'basis': '6-31g', 'method': 'hf'}}}
        other = {'version': 2, 'state': 0,
                 'model_config': {'input': {'basis': '6-31g',
                                            'functional': 'bhhlyp'}}}

        a = _displacement_run_signature(geom, atoms, 0.01, model)
        b = _displacement_run_signature(geom, atoms, 0.01, other)
        self.assertNotEqual(
            a, b, 'changing only the Hamiltonian must change the signature')

        # And the tag really is built from it, exactly as numerical_hess does.
        self.assertNotEqual(f'{a}c0p', f'{b}c0p')
        self.assertTrue(f'{a}c0p'.startswith(a))


if __name__ == '__main__':
    unittest.main()


class GuessMetadataIsOnlyImportedWhenItIsComparable(unittest.TestCase):
    """Matching coordinates do not by themselves make a producer's symmetry
    results usable.

    Two further conditions matter. Detection resolves a point group and an
    operation list from the geometry *under the settings in force*, so a
    producer running a looser `tolerance` describes the same coordinates with
    operations the reader would have rejected -- and with
    `use_response_symmetry=true` the reader then blocks the response with them.
    And results derived from the producer's WAVEFUNCTION (orbital, state and
    mode labels) are not geometry-derived at all: `stage_response_symmetry`
    consumes `mo_labels` verbatim whenever it is present with status 'ok', so an
    imported set sizes the blocking from the producer's orbital count.
    """

    COORD = [0.0, 0.0, 0.0, 0.0, 0.0, 1.4]

    def _molecule(self, **over):
        from oqp.molecule.molecule import Molecule
        mol = Molecule.__new__(Molecule)
        mol.tag = []
        mol.config = {}
        mol.config_tag = {}
        mol._state_tracking_fresh = True
        mol._system = np.array(self.COORD)
        mol.get_system = lambda: mol._system
        meta = {
            'status': 'enabled', 'enabled': True,
            'use_integral_symmetry': False, 'use_response_symmetry': True,
            'tolerance': 1.0e-5,
            'requested_point_group': 'auto', 'requested_subgroup': 'auto',
            'point_group': 'c1', 'subgroup': 'c1',
            'requested_matches_detected': True,
        }
        meta.update(over)
        mol.symmetry_metadata = meta
        return mol

    def _producer(self, **over):
        block = {
            'tolerance': 1.0e-5,
            'requested_point_group': 'auto', 'requested_subgroup': 'auto',
            'point_group': 'c2v', 'subgroup': 'c2v',
            'detection': {'point_group': 'c2v', 'operations': [1, 2, 3, 4]},
            'requested_matches_detected': False,
            'mo_labels': {'status': 'ok',
                          'alpha': {'labels': ['a1', 'b2', 'a1']}},
            'state_labels': {'status': 'ok'},
        }
        block.update(over)
        return block

    def test_matching_settings_still_import_geometry_results(self):
        """The feature must keep working when the two jobs really do agree."""
        mol = self._molecule()
        mol.put_data({'coord': self.COORD,
                      'symmetry_metadata': self._producer()})
        self.assertEqual(
            mol.symmetry_metadata['detection']['point_group'], 'c2v')

    def test_a_looser_producer_tolerance_blocks_the_import(self):
        mol = self._molecule()
        mol.put_data({'coord': self.COORD,
                      'symmetry_metadata': self._producer(tolerance=1.0e-3)})
        self.assertNotIn(
            'detection', mol.symmetry_metadata,
            'operations found under a tolerance this job rejects must not be '
            'imported')
        self.assertEqual(mol.symmetry_metadata['point_group'], 'c1')

    def test_a_different_requested_subgroup_blocks_the_import(self):
        mol = self._molecule()
        mol.put_data({'coord': self.COORD,
                      'symmetry_metadata': self._producer(
                          requested_subgroup='c2')})
        self.assertNotIn('detection', mol.symmetry_metadata)

    def test_an_older_guess_without_the_settings_blocks_the_import(self):
        """'unknown' has to read as 'not mine' -- files written before the
        settings were recorded must not be trusted."""
        producer = self._producer()
        del producer['tolerance']
        mol = self._molecule()
        mol.put_data({'coord': self.COORD, 'symmetry_metadata': producer})
        self.assertNotIn('detection', mol.symmetry_metadata)

    def test_wavefunction_labels_are_never_imported(self):
        """Even with the geometry AND every detection setting matching."""
        mol = self._molecule()
        mol.put_data({'coord': self.COORD,
                      'symmetry_metadata': self._producer()})
        for key in ('mo_labels', 'state_labels'):
            self.assertNotIn(
                key, mol.symmetry_metadata,
                f'{key} comes from the producer wavefunction, not its geometry')

    def test_the_readers_own_consistency_flag_survives(self):
        mol = self._molecule()
        mol.put_data({'coord': self.COORD,
                      'symmetry_metadata': self._producer()})
        self.assertTrue(
            mol.symmetry_metadata['requested_matches_detected'],
            'this flag describes the reader, not the producer')


class SymmetryUniqueDeclinesOnABrokenElectronicSolution(unittest.TestCase):
    """The orbit reduction is a statement about the nuclei; the identity it
    relies on also needs the electronic solution to respect those operations.

    A symmetry-broken unrestricted branch sits at a symmetric geometry with a
    density that does not transform, so an operation maps the branch being
    differentiated onto a different one.
    """

    def _hessian(self, labels):
        from oqp.library.single_point import Hessian

        hess = Hessian.__new__(Hessian)
        mol = types.SimpleNamespace()
        mol.config = {'hess': {'symmetry_unique': True}}
        mol.symmetry_metadata = {'status': 'enabled', 'tolerance': 1.0e-5}
        mol.get_atoms = lambda: np.array([8.0, 1.0, 1.0])
        mol.data = {'nelec_A': np.array([5]), 'nelec_B': np.array([5])}
        mol.label_molecular_orbitals = lambda: labels
        mol._parse_bool_like = staticmethod(
            lambda v: str(v).strip().lower() in ('true', 't', '1', 'yes'))
        mol._parse_bool_like = lambda v: str(v).strip().lower() in (
            'true', 't', '1', 'yes')
        hess.mol = mol
        hess.state = 0
        return hess

    GEOM = np.array([0.0, 0.0, -0.0776,
                     -1.0076, 1.0076, -1.1612,
                     1.0076, -1.0076, -1.1612])

    def test_a_mixed_occupied_orbital_declines_the_reduction(self):
        symmetric = {'status': 'ok',
                     'alpha': {'labels': ['a1'] * 5 + ['b2'] * 8}}
        broken = {'status': 'ok',
                  'alpha': {'labels': ['a1', 'a1', 'mixed', 'b1', 'a1']
                            + ['b2'] * 8}}

        good = self._hessian(symmetric)._symmetry_unique_displacements(
            self.GEOM, 0.01)
        self.assertIsNotNone(
            good[0], 'a symmetric solution must still get the reduction')

        bad = self._hessian(broken)._symmetry_unique_displacements(
            self.GEOM, 0.01)
        self.assertEqual(
            bad, (None, None),
            'a symmetry-broken reference must fall back to the full 6N set')

    def test_the_decline_reason_is_recorded(self):
        broken = {'status': 'ok',
                  'alpha': {'labels': ['a1', 'mixed', 'a1', 'b1', 'a1']}}
        hess = self._hessian(broken)
        hess._symmetry_unique_displacements(self.GEOM, 0.01)
        self.assertEqual(
            hess.mol.symmetry_metadata['hess_symmetry_unique']['status'],
            'symmetry_broken_electronic_solution')

    def test_unavailable_labels_decline_rather_than_assume(self):
        hess = self._hessian({'status': 'skipped_no_basis'})
        self.assertEqual(
            hess._symmetry_unique_displacements(self.GEOM, 0.01), (None, None))
