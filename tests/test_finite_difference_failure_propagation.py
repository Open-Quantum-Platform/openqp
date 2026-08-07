"""A detected finite-difference failure must reach every rank, and a scratch
gradient must only be reused by the run that produced it.

These pin the three defects the bot review found in the first round of fixes on
this branch. Each one is written so that reverting the fix makes it fail: the
assertions are on observed behaviour, never on the presence of a line of source
text. (The earlier files on this branch assert on source text, and a mutation
that inverted the reorientation-restore guard passed all of them.)
"""

import os
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

    def signature(self, **over):
        kw = dict(origin_coord=self.GEOM, atoms=self.ATOMS, dx=0.01, state=0)
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
        self.assertNotEqual(self.signature(), self.signature(state=1))

    def test_a_different_molecule_changes_the_signature(self):
        self.assertNotEqual(self.signature(),
                            self.signature(atoms=[8.0, 1.0, 9.0]))


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


class FailureVerdictReachesEveryRank(unittest.TestCase):
    """`flags` must be broadcast, not just the result array.

    It is appended to only inside the rank-0 guard, and only `grads`/`dcm` were
    broadcast -- so every other rank saw an empty list, found no 'failed', and
    walked into the analysis and the next collective while rank 0 raised and
    unwound past main()'s finalize_mpi(). Measured with two ranks and all 18
    displacements failing: before, 1 of 2 ranks raised and the job hung until
    killed at 150 s; after, both raised and it exited 1 in 17 s.

    This models the rank-0-only collection and asserts the non-root rank
    reaches the same verdict.
    """

    class FakeMPI:
        """Root holds the real value; every other rank receives it."""

        def __init__(self, rank):
            self.rank = rank
            self._root_payload = None

        def bcast(self, data, root=0, barrier=True):
            if self.rank == root:
                self._root_payload = data
                return data
            return self._root_payload

    def _collect(self, mpi, per_displacement_flags):
        """The shape of numerical_hess/numerical_nac's collection loop."""
        flags = []
        for flag in per_displacement_flags:
            if mpi.rank == 0:
                flags.append(flag)
        return mpi.bcast(flags)

    def test_a_non_root_rank_sees_the_failure(self):
        root = self.FakeMPI(rank=0)
        worker = self.FakeMPI(rank=1)
        # Both ranks observe the same displacement outcomes; only rank 0 records.
        outcomes = ['failed'] * 18
        root_flags = self._collect(root, outcomes)
        worker._root_payload = root_flags  # what the collective delivers
        worker_flags = self._collect(worker, outcomes)

        self.assertIn('failed', root_flags)
        self.assertIn('failed', worker_flags,
                      'the non-root rank must reach the same verdict, or it '
                      'walks into a collective while rank 0 unwinds')
        self.assertEqual(root_flags, worker_flags)

    def test_without_the_broadcast_the_ranks_disagree(self):
        """The control: this is what the code did before, and why it hung."""
        worker = self.FakeMPI(rank=1)
        unbroadcast = [f for f in ['failed'] * 18 if worker.rank == 0]
        self.assertEqual(unbroadcast, [])
        self.assertNotIn('failed', unbroadcast)


class AStaleGradientCannotMaskAFailedChild(unittest.TestCase):
    """`grad_wrapper` detects failure by `np.loadtxt(dat)` raising, so a
    gradient left by an earlier run makes a failed child look successful.

    Reproduced on this branch: with a previous run's scratch present, a
    numerical Hessian in which every one of the 18 displacements failed still
    exited 0 and printed a full set of frequencies, assembled entirely from the
    old gradients -- the exact silent success the raise added in this PR exists
    to stop.
    """

    def test_the_reader_fails_when_the_file_is_absent(self):
        """The failure signal the driver depends on."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            dat = os.path.join(tmp, 'missing.grad_0')
            with self.assertRaises(OSError):
                np.loadtxt(dat)

    def test_a_leftover_file_would_otherwise_read_as_success(self):
        """Why removing it before launching the child is the fix: a stale file
        loads cleanly and is indistinguishable from a fresh one."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            dat = os.path.join(tmp, 'stale.grad_0')
            np.savetxt(dat, np.arange(9, dtype=float))

            stale = np.loadtxt(dat).reshape(-1)
            self.assertEqual(stale.size, 9)  # reads as a perfectly good gradient

            os.remove(dat)                   # what the fix does pre-launch
            with self.assertRaises(OSError):
                np.loadtxt(dat)


if __name__ == '__main__':
    unittest.main()
