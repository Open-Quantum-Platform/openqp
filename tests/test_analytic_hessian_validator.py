import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

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
        analytic = [
            [1.0, 2.0, 3.0],
            [2.2, 5.0, 6.0],
            [3.0, 6.0, 9.0],
        ]
        reference = [
            [1.0, 2.0, 2.5],
            [2.0, 4.0, 6.0],
            [2.5, 6.0, 8.0],
        ]

        summary = validator.compare_hessians(analytic, reference, top_n=2)

        self.assertEqual(summary["shape"], [3, 3])
        self.assertAlmostEqual(summary["max_abs_diff"], 1.0)
        self.assertAlmostEqual(summary["rms_diff"], (2.54 / 9.0) ** 0.5)
        self.assertAlmostEqual(summary["max_asymmetry"], 0.2)
        self.assertEqual(len(summary["largest_components"]), 2)
        self.assertEqual(summary["largest_components"][0]["row"], 1)
        self.assertEqual(summary["largest_components"][0]["col"], 1)
        self.assertAlmostEqual(summary["largest_components"][0]["diff"], 1.0)

    def test_compare_hessians_rejects_shape_mismatch_before_reporting(self):
        validator = load_validator()

        with self.assertRaisesRegex(ValueError, "same shape"):
            validator.compare_hessians([[0.0, 0.0], [0.0, 0.0]], [[0.0, 0.0, 0.0]])

    def test_compare_hessians_rejects_nan_or_infinite_components(self):
        validator = load_validator()

        with self.assertRaisesRegex(ValueError, r"analytic\[0\]\[1\] must be finite"):
            validator.compare_hessians([[0.0, float("nan")], [0.0, 0.0]], [[0.0, 0.0], [0.0, 0.0]])

        with self.assertRaisesRegex(ValueError, r"reference\[1\]\[0\] must be finite"):
            validator.compare_hessians([[0.0, 0.0], [0.0, 0.0]], [[0.0, 0.0], [float("inf"), 0.0]])

    def test_summary_json_is_stable_and_does_not_include_full_matrices(self):
        validator = load_validator()
        summary = validator.compare_hessians([[1.0, 0.0], [0.0, 1.0]], [[1.0, 0.0], [0.0, 1.0]], top_n=1)

        payload = validator.summary_to_json(summary)

        self.assertIn('"max_abs_diff": 0.0', payload)
        self.assertNotIn('"analytic"', payload)
        self.assertNotIn('"reference"', payload)

    def test_summary_metadata_marks_compact_schema_without_full_matrices(self):
        validator = load_validator()
        summary = validator.compare_hessians([[1.0, 0.0], [0.0, 1.0]], [[1.0, 0.0], [0.0, 1.0]])

        self.assertEqual(summary["schema_version"], "analytic_hessian_validation.v1")
        self.assertEqual(summary["report_type"], "metric_comparison")
        self.assertEqual(summary["matrix_payload"], "omitted")
        self.assertNotIn("analytic_matrix", summary)
        self.assertNotIn("reference_matrix", summary)

    def test_load_matrix_accepts_openqp_hess_json_files(self):
        validator = load_validator()
        with tempfile.TemporaryDirectory() as tmpdir:
            hess_json = Path(tmpdir) / "h2.hess.json"
            hess_json.write_text(json.dumps({"hessian": [[1.0, 0.0], [0.0, 2.0]], "freqs": []}))

            matrix = validator._load_matrix(hess_json)

        self.assertEqual(matrix, [[1.0, 0.0], [0.0, 2.0]])

    def test_load_matrix_rejects_json_without_hessian_field(self):
        validator = load_validator()
        with tempfile.TemporaryDirectory() as tmpdir:
            bad_json = Path(tmpdir) / "missing.hess.json"
            bad_json.write_text(json.dumps({"freqs": []}))

            with self.assertRaisesRegex(ValueError, "must contain a hessian field"):
                validator._load_matrix(bad_json)

    def test_build_validation_summary_includes_context_tolerances_and_pass_flag(self):
        validator = load_validator()
        analytic = [[1.0, 1.002], [1.001, 2.0]]
        reference = [[1.0, 1.0], [1.0, 2.0]]

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
        self.assertEqual(summary["schema_version"], "analytic_hessian_validation.v1")
        self.assertEqual(summary["report_type"], "contextual_validation")
        self.assertEqual(summary["matrix_payload"], "omitted")
        self.assertTrue(summary["passed"])

    def test_build_validation_summary_fails_when_either_tolerance_is_exceeded(self):
        validator = load_validator()
        analytic = [[1.0, 1.01], [1.0, 2.0]]
        reference = [[1.0, 1.0], [1.0, 2.0]]

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
        self.assertEqual(summary["failed_metric_details"][0]["metric"], "max_abs_diff")
        self.assertAlmostEqual(summary["failed_metric_details"][0]["observed"], 0.01)
        self.assertAlmostEqual(summary["failed_metric_details"][0]["tolerance"], 0.003)
        self.assertAlmostEqual(summary["failed_metric_details"][0]["excess"], 0.007)
        self.assertEqual(summary["failed_metric_details"][0]["worst_component"]["row"], 0)
        self.assertEqual(summary["failed_metric_details"][0]["worst_component"]["col"], 1)
        self.assertAlmostEqual(summary["failed_metric_details"][0]["worst_component"]["abs_diff"], 0.01)

    def test_build_validation_summary_can_fail_on_asymmetry_tolerance(self):
        validator = load_validator()
        analytic = [[1.0, 0.20], [0.0, 2.0]]
        reference = [[1.0, 0.10], [0.10, 2.0]]

        summary = validator.build_validation_summary(
            analytic,
            reference,
            method="hf",
            td_type="none",
            state=0,
            molecule="h2",
            basis="sto-3g",
            displacement=0.005,
            max_tolerance=0.2,
            rms_tolerance=0.2,
            asymmetry_tolerance=0.05,
        )

        self.assertFalse(summary["passed"])
        self.assertEqual(summary["tolerances"]["max_asymmetry"], 0.05)
        self.assertIn("max_asymmetry", summary["failed_metrics"])
        asymmetry_detail = [
            detail for detail in summary["failed_metric_details"] if detail["metric"] == "max_asymmetry"
        ][0]
        self.assertAlmostEqual(asymmetry_detail["observed"], 0.2)
        self.assertAlmostEqual(asymmetry_detail["tolerance"], 0.05)
        self.assertNotIn("worst_component", asymmetry_detail)

    def test_build_validation_summary_rejects_invalid_context_scalars(self):
        validator = load_validator()
        matrix = [[1.0, 0.0], [0.0, 1.0]]

        with self.assertRaisesRegex(ValueError, "displacement must be finite and positive"):
            validator.build_validation_summary(
                matrix,
                matrix,
                method="hf",
                td_type="none",
                state=0,
                molecule="h2",
                basis="sto-3g",
                displacement=0.0,
                max_tolerance=0.003,
                rms_tolerance=0.002,
            )

        with self.assertRaisesRegex(ValueError, "max_tolerance must be finite and positive"):
            validator.build_validation_summary(
                matrix,
                matrix,
                method="hf",
                td_type="none",
                state=0,
                molecule="h2",
                basis="sto-3g",
                displacement=0.005,
                max_tolerance=float("nan"),
                rms_tolerance=0.002,
            )

        with self.assertRaisesRegex(ValueError, "rms_tolerance must be finite and positive"):
            validator.build_validation_summary(
                matrix,
                matrix,
                method="hf",
                td_type="none",
                state=0,
                molecule="h2",
                basis="sto-3g",
                displacement=0.005,
                max_tolerance=0.003,
                rms_tolerance=-0.002,
            )

    def test_cli_requires_td_type_when_any_validation_metadata_is_supplied(self):
        validator = load_validator()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            analytic = tmp / "analytic.txt"
            reference = tmp / "reference.txt"
            analytic.write_text("1.0 0.0\n0.0 1.0\n")
            reference.write_text("1.0 0.0\n0.0 1.0\n")

            with self.assertRaises(SystemExit) as err:
                validator.main(
                    [
                        str(analytic),
                        str(reference),
                        "--method",
                        "hf",
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
                    ]
                )

        self.assertNotEqual(err.exception.code, 0)

    def test_cli_noncontext_summary_records_matrix_sources_without_full_matrices(self):
        validator = load_validator()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            analytic = tmp / "analytic.txt"
            reference = tmp / "reference.txt"
            output = tmp / "summary.json"
            analytic.write_text("1.0 0.0\n0.0 1.0\n")
            reference.write_text("1.0 0.0\n0.0 1.0\n")
            expected_analytic_sha = hashlib.sha256(analytic.read_bytes()).hexdigest()
            expected_reference_sha = hashlib.sha256(reference.read_bytes()).hexdigest()

            status = validator.main([str(analytic), str(reference), "--output", str(output)])

            payload = json.loads(output.read_text())
        self.assertEqual(status, 0)
        self.assertEqual(payload["matrix_sources"]["analytic"], str(analytic))
        self.assertEqual(payload["matrix_sources"]["reference"], str(reference))
        self.assertEqual(payload["matrix_source_sha256"]["analytic"], expected_analytic_sha)
        self.assertEqual(payload["matrix_source_sha256"]["reference"], expected_reference_sha)
        self.assertNotIn("analytic_matrix", payload)
        self.assertNotIn("reference_matrix", payload)

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

    def test_cli_returns_nonzero_when_contextual_validation_fails_tolerances(self):
        validator = load_validator()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            analytic = tmp / "analytic.txt"
            reference = tmp / "reference.txt"
            output = tmp / "summary.json"
            analytic.write_text("1.0 1.010\n1.000 2.0\n")
            reference.write_text("1.0 1.000\n1.000 2.0\n")

            status = validator.main(
                [
                    str(analytic),
                    str(reference),
                    "--method",
                    "tdhf",
                    "--td-type",
                    "rpa",
                    "--state",
                    "1",
                    "--molecule",
                    "h2o",
                    "--basis",
                    "sto-3g",
                    "--displacement",
                    "0.005",
                    "--max-tolerance",
                    "0.003",
                    "--rms-tolerance",
                    "0.020",
                    "--output",
                    str(output),
                ]
            )

            payload = json.loads(output.read_text())
        self.assertEqual(status, 1)
        self.assertFalse(payload["passed"])
        self.assertEqual(payload["failed_metrics"], ["max_abs_diff"])

    def test_cli_reports_invalid_context_scalar_as_usage_error(self):
        validator = load_validator()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            analytic = tmp / "analytic.txt"
            reference = tmp / "reference.txt"
            analytic.write_text("1.0 0.0\n0.0 1.0\n")
            reference.write_text("1.0 0.0\n0.0 1.0\n")

            with self.assertRaises(SystemExit) as err:
                validator.main(
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
                        "-0.005",
                        "--max-tolerance",
                        "0.003",
                        "--rms-tolerance",
                        "0.002",
                    ]
                )

        self.assertNotEqual(err.exception.code, 0)


if __name__ == "__main__":
    unittest.main()
