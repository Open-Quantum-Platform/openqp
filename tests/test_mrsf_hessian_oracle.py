import unittest

import numpy as np

from pyoqp.oqp.library.mrsf_hessian import (
    MRSFRootTrackingError,
    assemble_root_tracked_central_fd_hessian,
    compare_fd_step_hessians,
    physical_root_label,
    assess_root_tracking,
)


class MRSFHessianRootTrackingTests(unittest.TestCase):
    def test_physical_root_labels_follow_openqp_mrsf_convention(self):
        self.assertEqual(physical_root_label(0), "high-spin reference")
        self.assertEqual(physical_root_label(1), "physical S0")
        self.assertEqual(physical_root_label(2), "physical S1")

    def test_assess_root_tracking_accepts_same_character_target_root(self):
        reference = np.eye(3)
        displaced = np.eye(3)

        diagnostic = assess_root_tracking(
            requested_state=2,
            reference_vectors=reference,
            displaced_vectors=displaced,
            reference_energies=[-75.0, -74.9, -74.8],
            displaced_energies=[-75.0, -74.899, -74.801],
        )

        self.assertTrue(diagnostic.ok)
        self.assertEqual(diagnostic.requested_state, 2)
        self.assertEqual(diagnostic.matched_state, 2)
        self.assertGreaterEqual(diagnostic.overlap, 0.99)
        self.assertFalse(diagnostic.root_flip)

    def test_assess_root_tracking_fails_when_vector_overlap_selects_different_root(self):
        reference = np.eye(3)
        displaced = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
        ])

        with self.assertRaisesRegex(MRSFRootTrackingError, "root flip.*requested root 2.*matched root 1"):
            assess_root_tracking(
                requested_state=2,
                reference_vectors=reference,
                displaced_vectors=displaced,
                reference_energies=[-75.0, -74.9, -74.8],
                displaced_energies=[-75.0, -74.801, -74.899],
            )

    def test_assess_root_tracking_fails_when_overlap_is_ambiguous(self):
        reference = np.eye(3)
        displaced = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 0.71, 0.70],
            [0.0, 0.70, 0.71],
        ])

        with self.assertRaisesRegex(MRSFRootTrackingError, "ambiguous.*overlap"):
            assess_root_tracking(
                requested_state=1,
                reference_vectors=reference,
                displaced_vectors=displaced,
                reference_energies=[-75.0, -74.9000, -74.8999],
                displaced_energies=[-75.0, -74.89995, -74.89994],
                min_overlap=0.90,
                min_overlap_gap=0.10,
                min_energy_gap=1.0e-3,
            )

    def test_assemble_root_tracked_central_fd_hessian_requires_clean_plus_and_minus_roots(self):
        reference_vectors = np.eye(3)
        gradients_plus = [np.array([1.0, 2.0]), np.array([3.0, 5.0])]
        gradients_minus = [np.array([0.0, 1.0]), np.array([2.0, 3.0])]
        displaced_vectors = [np.eye(3), np.eye(3), np.eye(3), np.eye(3)]
        displaced_energies = [
            [-75.0, -74.9, -74.8],
            [-75.0, -74.9, -74.8],
            [-75.0, -74.9, -74.8],
            [-75.0, -74.9, -74.8],
        ]

        hessian, diagnostics = assemble_root_tracked_central_fd_hessian(
            requested_state=2,
            reference_vectors=reference_vectors,
            reference_energies=[-75.0, -74.9, -74.8],
            gradients_plus=gradients_plus,
            gradients_minus=gradients_minus,
            displaced_vectors=displaced_vectors,
            displaced_energies=displaced_energies,
            dx=0.5,
        )

        np.testing.assert_allclose(hessian, [[1.0, 1.0], [1.0, 2.0]])
        self.assertEqual(len(diagnostics), 4)
        self.assertTrue(all(item.ok for item in diagnostics))

    def test_assemble_root_tracked_central_fd_hessian_rejects_root_flip_before_assembly(self):
        reference_vectors = np.eye(3)
        flipped = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
        ])

        with self.assertRaisesRegex(MRSFRootTrackingError, "root flip"):
            assemble_root_tracked_central_fd_hessian(
                requested_state=2,
                reference_vectors=reference_vectors,
                reference_energies=[-75.0, -74.9, -74.8],
                gradients_plus=[np.array([1.0, 2.0])],
                gradients_minus=[np.array([0.0, 1.0])],
                displaced_vectors=[flipped, np.eye(3)],
                displaced_energies=[[-75.0, -74.8, -74.9], [-75.0, -74.9, -74.8]],
                dx=0.5,
            )

    def test_compare_fd_step_hessians_reports_stable_symmetric_step(self):
        fine = np.array([[1.0, 0.1], [0.1, 2.0]])
        coarser = np.array([[1.00001, 0.10001], [0.10001, 2.00001]])

        diagnostic = compare_fd_step_hessians(
            dx=3.0e-4,
            hessian=fine,
            reference_dx=1.0e-3,
            reference_hessian=coarser,
            max_step_delta=1.0e-3,
            max_asymmetry=1.0e-10,
        )

        self.assertTrue(diagnostic.ok)
        self.assertLess(diagnostic.max_abs_delta, 1.0e-3)
        self.assertEqual(diagnostic.shape, [2, 2])

    def test_compare_fd_step_hessians_rejects_asymmetric_or_step_unstable_matrix(self):
        bad = np.array([[1.0, 0.5], [0.1, 2.0]])
        ref = np.array([[1.0, 0.1], [0.1, 2.0]])

        diagnostic = compare_fd_step_hessians(
            dx=1.0e-4,
            hessian=bad,
            reference_dx=3.0e-4,
            reference_hessian=ref,
            max_step_delta=1.0e-3,
            max_asymmetry=1.0e-3,
        )

        self.assertFalse(diagnostic.ok)
        self.assertGreater(diagnostic.max_asymmetry, 1.0e-3)
        self.assertGreater(diagnostic.max_abs_delta, 1.0e-3)


if __name__ == "__main__":
    unittest.main()
