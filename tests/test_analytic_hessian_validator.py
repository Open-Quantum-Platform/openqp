import importlib.util
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_validator():
    spec = importlib.util.spec_from_file_location(
        "validate_analytic_hessian_under_test",
        ROOT / "tools" / "validate_analytic_hessian.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AnalyticHessianValidatorTests(unittest.TestCase):
    def test_compare_hessians_reports_max_rms_symmetry_and_largest_components(self):
        validator = load_validator()
        analytic = np.array(
            [
                [1.0, 2.0, 3.0],
                [2.2, 5.0, 6.0],
                [3.0, 6.0, 9.0],
            ]
        )
        reference = np.array(
            [
                [1.0, 2.0, 2.5],
                [2.0, 4.0, 6.0],
                [2.5, 6.0, 8.0],
            ]
        )

        summary = validator.compare_hessians(analytic, reference, top_n=2)

        self.assertEqual(summary["shape"], [3, 3])
        self.assertAlmostEqual(summary["max_abs_diff"], 1.0)
        self.assertAlmostEqual(summary["rms_diff"], np.sqrt(2.54 / 9.0))
        self.assertAlmostEqual(summary["max_asymmetry"], 0.2)
        self.assertEqual(len(summary["largest_components"]), 2)
        self.assertEqual(summary["largest_components"][0]["row"], 1)
        self.assertEqual(summary["largest_components"][0]["col"], 1)
        self.assertAlmostEqual(summary["largest_components"][0]["diff"], 1.0)

    def test_compare_hessians_rejects_shape_mismatch_before_reporting(self):
        validator = load_validator()

        with self.assertRaisesRegex(ValueError, "same shape"):
            validator.compare_hessians(np.zeros((2, 2)), np.zeros((3, 3)))

    def test_summary_json_is_stable_and_does_not_include_full_matrices(self):
        validator = load_validator()
        summary = validator.compare_hessians(np.eye(2), np.eye(2), top_n=1)

        payload = validator.summary_to_json(summary)

        self.assertIn('"max_abs_diff": 0.0', payload)
        self.assertNotIn('"analytic"', payload)
        self.assertNotIn('"reference"', payload)


if __name__ == "__main__":
    unittest.main()
