import unittest

import numpy as np

from pyoqp.oqp.library.mrsf_hessian import (
    MRSFRootTrackingError,
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


if __name__ == "__main__":
    unittest.main()
