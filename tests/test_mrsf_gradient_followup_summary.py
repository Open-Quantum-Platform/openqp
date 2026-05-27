import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "tools" / "mrsf_gradient_followup_summary.py"


def load_module():
    spec = importlib.util.spec_from_file_location("mrsf_gradient_followup_summary", MODULE_PATH)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class MrsfGradientFollowupSummaryTests(unittest.TestCase):
    def write_compact_csv(self, rows):
        tmp = tempfile.NamedTemporaryFile("w", newline="", suffix=".csv", delete=False)
        fieldnames = [
            "method",
            "molecule",
            "root",
            "physical_state",
            "max_abs_diff_ha_per_bohr",
            "rms_diff_ha_per_bohr",
            "worst_component",
            "worst_analytic",
            "worst_fd",
            "trah_total",
            "failed_count",
            "elapsed_s",
        ]
        with tmp:
            writer = csv.DictWriter(tmp, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        return Path(tmp.name)

    def test_summary_ranks_target_post_ovov_failures_without_relabeling_reference(self):
        module = load_module()
        csv_path = self.write_compact_csv(
            [
                {
                    "method": "mrsf",
                    "molecule": "ch2o",
                    "root": "4",
                    "physical_state": "",
                    "max_abs_diff_ha_per_bohr": "0.02121502",
                    "rms_diff_ha_per_bohr": "0.02121502",
                    "worst_component": "a1_z",
                    "worst_analytic": "-0.25282498",
                    "worst_fd": "-0.27404",
                    "trah_total": "0",
                    "failed_count": "0",
                    "elapsed_s": "2.18",
                },
                {
                    "method": "mrsf",
                    "molecule": "h2o",
                    "root": "2",
                    "physical_state": "",
                    "max_abs_diff_ha_per_bohr": "0.0000545",
                    "rms_diff_ha_per_bohr": "0.0000545",
                    "worst_component": "a0_z",
                    "worst_analytic": "-0.1452205",
                    "worst_fd": "-0.145275",
                    "trah_total": "0",
                    "failed_count": "0",
                    "elapsed_s": "1.82",
                },
                {
                    "method": "mrsf",
                    "molecule": "nh3",
                    "root": "4",
                    "physical_state": "",
                    "max_abs_diff_ha_per_bohr": "0.0922304",
                    "rms_diff_ha_per_bohr": "0.0922304",
                    "worst_component": "a0_y",
                    "worst_analytic": "-0.0922604",
                    "worst_fd": "-0.00003",
                    "trah_total": "1",
                    "failed_count": "0",
                    "elapsed_s": "1.77",
                },
            ]
        )

        summary = module.summarize_compact_csv(csv_path, threshold=1.0e-3)

        self.assertEqual(3, summary["total_cases"])
        self.assertEqual(2, summary["failure_count"])
        self.assertEqual("nh3", summary["failures"][0]["molecule"])
        self.assertEqual("S3", summary["failures"][0]["physical_state"])
        self.assertEqual("trah_or_failed", summary["failures"][0]["classification"])
        self.assertEqual("ch2o", summary["failures"][1]["molecule"])
        self.assertEqual("S3", summary["failures"][1]["physical_state"])
        self.assertEqual("fd_mismatch_no_trah", summary["failures"][1]["classification"])
        self.assertEqual(["h2o root 2"], summary["clean_controls"])

    def test_component_summary_flags_state_character_changes_before_algebra_hints(self):
        module = load_module()
        tmp = tempfile.NamedTemporaryFile("w", newline="", suffix=".csv", delete=False)
        fieldnames = [
            "molecule",
            "method",
            "root",
            "physical_state",
            "component",
            "analytic_ha_per_bohr",
            "fd_ha_per_bohr",
            "diff_ha_per_bohr",
            "abs_diff_ha_per_bohr",
            "trah_count",
            "failed_any",
            "s2_grad",
            "s2_plus",
            "s2_minus",
        ]
        with tmp:
            writer = csv.DictWriter(tmp, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(
                {
                    "molecule": "nh3",
                    "method": "mrsf",
                    "root": "4",
                    "physical_state": "S3",
                    "component": "a0_y",
                    "analytic_ha_per_bohr": "-0.0922604",
                    "fd_ha_per_bohr": "-0.00003",
                    "diff_ha_per_bohr": "-0.0922304",
                    "abs_diff_ha_per_bohr": "0.0922304",
                    "trah_count": "0",
                    "failed_any": "False",
                    "s2_grad": "{3: 0.01, 4: 0.02}",
                    "s2_plus": "{3: 0.02, 4: 1.85}",
                    "s2_minus": "{3: 0.01, 4: 0.03}",
                }
            )
            writer.writerow(
                {
                    "molecule": "h2s",
                    "method": "mrsf",
                    "root": "5",
                    "physical_state": "S4",
                    "component": "a0_z",
                    "analytic_ha_per_bohr": "-0.14875174",
                    "fd_ha_per_bohr": "-0.22803",
                    "diff_ha_per_bohr": "0.07927826",
                    "abs_diff_ha_per_bohr": "0.07927826",
                    "trah_count": "0",
                    "failed_any": "False",
                    "s2_grad": "{5: 0.01}",
                    "s2_plus": "{5: 0.01}",
                    "s2_minus": "{5: 0.01}",
                }
            )

        summary = module.summarize_components_csv(Path(tmp.name), threshold=1.0e-3)

        self.assertEqual(2, summary["component_group_count"])
        self.assertEqual("nh3", summary["groups"][0]["molecule"])
        self.assertTrue(summary["groups"][0]["possible_state_character_change"])
        self.assertEqual("root_tracking_or_state_character_change", summary["groups"][0]["mechanism_hint"])
        self.assertEqual("h2s", summary["groups"][1]["molecule"])
        self.assertFalse(summary["groups"][1]["possible_state_character_change"])
        self.assertEqual("localized_z_component_z_vector_or_operator_mapping", summary["groups"][1]["mechanism_hint"])

    def test_component_summary_preserves_bad_component_details_for_next_diagnostic(self):
        module = load_module()
        rows = [
            {
                "method": "mrsf",
                "molecule": "hcn",
                "root": 6,
                "physical_state": "S5",
                "component": "a0_z",
                "axis": "z",
                "analytic_ha_per_bohr": 0.100,
                "fd_ha_per_bohr": 0.020,
                "diff_ha_per_bohr": 0.080,
                "abs_diff_ha_per_bohr": 0.080,
                "trah_count": 0,
                "failed_any": False,
                "s2_max_delta": 0.0,
                "s2_evidence": "present",
                "bad_component": True,
                "target_case": True,
            },
            {
                "method": "mrsf",
                "molecule": "hcn",
                "root": 6,
                "physical_state": "S5",
                "component": "a1_x",
                "axis": "x",
                "analytic_ha_per_bohr": 0.010,
                "fd_ha_per_bohr": 0.0104,
                "diff_ha_per_bohr": -0.0004,
                "abs_diff_ha_per_bohr": 0.0004,
                "trah_count": 0,
                "failed_any": False,
                "s2_max_delta": 0.0,
                "s2_evidence": "present",
                "bad_component": False,
                "target_case": True,
            },
        ]

        summary = module.summarize_component_rows(rows, threshold=1.0e-3)

        self.assertEqual(
            [
                {
                    "component": "a0_z",
                    "axis": "z",
                    "abs_diff_ha_per_bohr": 0.080,
                    "analytic_ha_per_bohr": 0.100,
                    "fd_ha_per_bohr": 0.020,
                    "s2_evidence": "present",
                }
            ],
            summary["groups"][0]["bad_components"],
        )
        self.assertEqual(summary["groups"][0]["bad_components"], summary["target_bad_groups"][0]["bad_components"])

    def test_component_summary_marks_empty_or_all_zero_s2_maps_as_unknown_evidence(self):
        module = load_module()
        rows = [
            {
                "method": "mrsf",
                "molecule": "ch2o",
                "root": 4,
                "physical_state": "S3",
                "component": "a1_z",
                "axis": "z",
                "analytic_ha_per_bohr": -0.25282498,
                "fd_ha_per_bohr": -0.27404,
                "diff_ha_per_bohr": 0.02121502,
                "abs_diff_ha_per_bohr": 0.02121502,
                "trah_count": 0,
                "failed_any": False,
                "s2_max_delta": 0.0,
                "s2_evidence": "unknown",
                "bad_component": True,
            }
        ]

        summary = module.summarize_component_rows(rows, threshold=1.0e-3)

        self.assertEqual("unknown", summary["groups"][0]["s2_evidence"])
        self.assertFalse(summary["groups"][0]["possible_state_character_change"])
        self.assertIn("state-character evidence missing", summary["groups"][0]["mechanism_hint"])

    def test_component_summary_counts_target_remaining_highroot_groups(self):
        module = load_module()
        rows = [
            {
                "method": "mrsf",
                "molecule": "ch2o",
                "root": 4,
                "physical_state": "S3",
                "component": "a1_z",
                "axis": "z",
                "analytic_ha_per_bohr": -0.25282498,
                "fd_ha_per_bohr": -0.27404,
                "diff_ha_per_bohr": 0.02121502,
                "abs_diff_ha_per_bohr": 0.02121502,
                "trah_count": 0,
                "failed_any": False,
                "s2_max_delta": 0.0,
                "s2_evidence": "unknown",
                "bad_component": True,
            },
            {
                "method": "mrsf",
                "molecule": "h2o",
                "root": 2,
                "physical_state": "S1",
                "component": "a0_z",
                "axis": "z",
                "analytic_ha_per_bohr": -0.1452205,
                "fd_ha_per_bohr": -0.145275,
                "diff_ha_per_bohr": 0.0000545,
                "abs_diff_ha_per_bohr": 0.0000545,
                "trah_count": 0,
                "failed_any": False,
                "s2_max_delta": 0.0,
                "s2_evidence": "present",
                "bad_component": False,
            },
        ]

        summary = module.summarize_component_rows(rows, threshold=1.0e-3)

        self.assertEqual(1, summary["target_group_count"])
        self.assertEqual(1, summary["target_bad_group_count"])
        self.assertEqual(1, len(summary["target_bad_groups"]))
        self.assertEqual("ch2o", summary["target_bad_groups"][0]["molecule"])
        self.assertEqual(4, summary["target_bad_groups"][0]["root"])
        self.assertNotIn("h2o", {group["molecule"] for group in summary["target_bad_groups"]})
        self.assertTrue(summary["groups"][0]["target_case"])
        self.assertFalse(summary["groups"][1]["target_case"])

    def test_cli_components_mode_writes_component_group_summary(self):
        module = load_module()
        tmp = tempfile.NamedTemporaryFile("w", newline="", suffix=".csv", delete=False)
        output = tempfile.NamedTemporaryFile("r", suffix=".json", delete=False)
        output.close()
        fieldnames = [
            "molecule",
            "method",
            "root",
            "physical_state",
            "component",
            "analytic_ha_per_bohr",
            "fd_ha_per_bohr",
            "diff_ha_per_bohr",
            "abs_diff_ha_per_bohr",
            "trah_count",
            "failed_any",
            "s2_grad",
            "s2_plus",
            "s2_minus",
        ]
        with tmp:
            writer = csv.DictWriter(tmp, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(
                {
                    "molecule": "hcn",
                    "method": "mrsf",
                    "root": "6",
                    "physical_state": "",
                    "component": "a0_z",
                    "analytic_ha_per_bohr": "0.100",
                    "fd_ha_per_bohr": "0.020",
                    "diff_ha_per_bohr": "0.080",
                    "abs_diff_ha_per_bohr": "0.080",
                    "trah_count": "0",
                    "failed_any": "False",
                    "s2_grad": "{6: 0.02}",
                    "s2_plus": "{6: 0.02}",
                    "s2_minus": "{6: 0.02}",
                }
            )

        status = module.main(["--components", tmp.name, "--output", output.name])

        self.assertEqual(0, status)
        written = output.name and Path(output.name).read_text()
        self.assertIn('"component_group_count": 1', written)
        self.assertIn('"physical_state": "S5"', written)
        self.assertIn("localized_z_component", written)


if __name__ == "__main__":
    unittest.main()
