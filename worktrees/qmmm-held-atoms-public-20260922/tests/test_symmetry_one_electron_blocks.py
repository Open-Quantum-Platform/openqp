import importlib.util
import sys
from pathlib import Path

import numpy as np
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load_symmetry_module():
    spec = importlib.util.spec_from_file_location(
        "openqp_symmetry_under_test",
        str(ROOT / "pyoqp/oqp/library/symmetry.py"),
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class TestOneElectronBlockDiagnostics(unittest.TestCase):
    def setUp(self):
        self.symmetry = load_symmetry_module()

    def test_zero_off_block_leakage_for_block_diagonal_matrix(self):
        labels = ["A1", "A1", "B1", "B1"]
        transform = np.eye(4)
        matrix = np.array(
            [
                [1.0, 0.2, 0.0, 0.0],
                [0.2, 2.0, 0.0, 0.0],
                [0.0, 0.0, 3.0, 0.4],
                [0.0, 0.0, 0.4, 4.0],
            ]
        )

        payload = self.symmetry.one_electron_block_diagnostics(
            matrix,
            transform,
            labels,
            tolerance=1.0e-8,
        )

        self.assertEqual(payload["max_off_block_abs"], 0.0)
        self.assertEqual(payload["off_block_element_count"], 8)
        self.assertTrue(payload["within_tolerance"])
        self.assertTrue(payload["transform_orthogonality"]["ok"])
        self.assertFalse(payload["use_integral_symmetry"])
        self.assertFalse(payload["use_response_symmetry"])

    def test_nonzero_off_block_leakage_is_detected(self):
        labels = ["A1", "A1", "B1", "B1"]
        transform = np.eye(4)
        matrix = np.array(
            [
                [1.0, 0.2, 0.001, 0.0],
                [0.2, 2.0, 0.0, 0.0],
                [0.001, 0.0, 3.0, 0.4],
                [0.0, 0.0, 0.4, 4.0],
            ]
        )

        payload = self.symmetry.one_electron_block_diagnostics(
            matrix,
            transform,
            labels,
            tolerance=1.0e-4,
        )

        self.assertAlmostEqual(payload["max_off_block_abs"], 0.001)
        self.assertFalse(payload["within_tolerance"])
        self.assertEqual(len(payload["max_off_block_indices"]), 1)

    def test_update_one_electron_block_diagnostics_updates_metadata(self):
        labels = ["A1", "A1", "B1", "B1"]
        transform = np.eye(4)
        metadata = {
            "status": "auto",
            "use_integral_symmetry": True,
            "use_response_symmetry": True,
        }
        matrices = {
            "overlap": np.diag([1.0, 1.0, 1.0, 1.0]),
            "hcore": np.diag([0.0, 0.0, 0.0, 0.0]),
        }

        updated = self.symmetry.update_one_electron_block_diagnostics(
            metadata,
            matrices,
            transform,
            labels,
            tolerance=1.0e-8,
        )

        self.assertIn("one_electron_block_diagnostics", updated)
        diagnostics = updated["one_electron_block_diagnostics"]
        self.assertIn("overlap", diagnostics)
        self.assertIn("hcore", diagnostics)
        self.assertEqual(updated["use_integral_symmetry"], False)
        self.assertEqual(updated["use_response_symmetry"], False)

        for name, entry in diagnostics.items():
            self.assertFalse(entry["use_integral_symmetry"], msg=f"{name} should keep reductions disabled")
            self.assertFalse(entry["use_response_symmetry"], msg=f"{name} should keep reductions disabled")

    def test_label_count_mismatch_rejected(self):
        labels = ["A1", "A1", "B1"]
        transform = np.eye(4)
        matrix = np.eye(4)

        with self.assertRaises(ValueError):
            self.symmetry.one_electron_block_diagnostics(
                matrix,
                transform,
                labels,
                tolerance=1.0e-8,
            )

    def test_reject_non_square_matrix(self):
        labels = ["A1", "A1", "B1", "B1"]
        transform = np.eye(4)
        matrix = np.array([[1.0, 0.0], [0.0, 1.0]])

        with self.assertRaises(ValueError):
            self.symmetry.one_electron_block_diagnostics(matrix, transform, labels, tolerance=1e-8)
