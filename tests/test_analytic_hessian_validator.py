import importlib.util
import json
import tempfile
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

    def test_build_validation_summary_includes_context_tolerances_and_pass_flag(self):
        validator = load_validator()
        analytic = np.array([[1.0, 1.002], [1.001, 2.0]])
        reference = np.array([[1.0, 1.0], [1.0, 2.0]])

        summary = validator.build_validation_summary(
            analytic,
            reference,
            method="hf",
            td_type="none",
            state=0,
            molecule="h2",
            basis="sto-3g",
            displacement=0.005,
            max_tolerance=0.003,
            rms_tolerance=0.002,
            top_n=1,
        )

        self.assertEqual(summary["method"], "hf")
        self.assertEqual(summary["td_type"], "none")
        self.assertEqual(summary["state"], 0)
        self.assertEqual(summary["molecule"], "h2")
        self.assertEqual(summary["basis"], "sto-3g")
        self.assertAlmostEqual(summary["displacement"], 0.005)
        self.assertAlmostEqual(summary["tolerances"]["max_abs_diff"], 0.003)
        self.assertAlmostEqual(summary["tolerances"]["rms_diff"], 0.002)
        self.assertTrue(summary["passed"])

    def test_build_validation_summary_fails_when_either_tolerance_is_exceeded(self):
        validator = load_validator()
        analytic = np.array([[1.0, 1.01], [1.0, 2.0]])
        reference = np.array([[1.0, 1.0], [1.0, 2.0]])

        summary = validator.build_validation_summary(
            analytic,
            reference,
            method="tdhf",
            td_type="rpa",
            state=1,
            molecule="h2o",
            basis="sto-3g",
            displacement=0.005,
            max_tolerance=0.003,
            rms_tolerance=0.02,
        )

        self.assertFalse(summary["passed"])
        self.assertIn("max_abs_diff", summary["failed_metrics"])

    def test_cli_writes_contextual_validation_summary_when_metadata_is_supplied(self):
        validator = load_validator()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            analytic = tmp / "analytic.txt"
            reference = tmp / "reference.txt"
            output = tmp / "summary.json"
            analytic.write_text("1.0 1.002\n1.001 2.0\n")
            reference.write_text("1.0 1.0\n1.0 2.0\n")

            status = validator.main(
                [
                    str(analytic),
                    str(reference),
                    "--method",
                    "hf",
                    "--td-type",
                    "none",
                    "--state",
                    "0",
                    "--molecule",
                    "h2",
                    "--basis",
                    "sto-3g",
                    "--displacement",
                    "0.005",
                    "--max-tolerance",
                    "0.003",
                    "--rms-tolerance",
                    "0.002",
                    "--output",
                    str(output),
                ]
            )

            payload = json.loads(output.read_text())
        self.assertEqual(status, 0)
        self.assertEqual(payload["method"], "hf")
        self.assertEqual(payload["basis"], "sto-3g")
        self.assertTrue(payload["passed"])


if __name__ == "__main__":
    unittest.main()
