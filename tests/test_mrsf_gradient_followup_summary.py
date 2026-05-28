import csv
import importlib.util
import json
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

    def test_multi_component_summary_preserves_source_csv_for_artifact_triage(self):
        module = load_module()
        first = [
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
                "target_case": True,
            }
        ]
        second = [
            {
                "method": "mrsf",
                "molecule": "h2s",
                "root": 5,
                "physical_state": "S4",
                "component": "a0_z",
                "axis": "z",
                "analytic_ha_per_bohr": -0.14875174,
                "fd_ha_per_bohr": -0.22803,
                "diff_ha_per_bohr": 0.07927826,
                "abs_diff_ha_per_bohr": 0.07927826,
                "trah_count": 0,
                "failed_any": False,
                "s2_max_delta": 0.0,
                "s2_evidence": "present",
                "bad_component": True,
                "target_case": True,
            }
        ]

        summary = module.summarize_component_datasets(
            [("first/components.csv", first), ("second/components.csv", second)],
            threshold=1.0e-3,
        )

        self.assertEqual(2, summary["dataset_count"])
        self.assertEqual(2, summary["target_bad_group_count"])
        self.assertEqual(
            {"first/components.csv", "second/components.csv"},
            {group["source_csv"] for group in summary["target_bad_groups"]},
        )
        self.assertEqual(
            {"first/components.csv", "second/components.csv"},
            {group["source_csv"] for group in summary["groups"]},
        )

    def write_components_csv(self, rows):
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
            for row in rows:
                writer.writerow(row)
        return Path(tmp.name)

    def test_cli_components_mode_writes_component_group_summary(self):
        module = load_module()
        csv_path = self.write_components_csv(
            [
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
            ]
        )
        output = tempfile.NamedTemporaryFile("r", suffix=".json", delete=False)
        output.close()

        status = module.main(["--components", str(csv_path), "--output", output.name])

        self.assertEqual(0, status)
        written = output.name and Path(output.name).read_text()
        self.assertIn('"component_group_count": 1', written)
        self.assertIn('"physical_state": "S5"', written)
        self.assertIn("localized_z_component", written)

    def test_parse_mrsf_log_extracts_state_character_table(self):
        module = load_module()
        log_text = """
  Spin-adapted spin-flip excitations
 State #   1  Energy =   -7.315353 eV
               <S^2> =    0.0000
 State #   4  Energy =    2.734056 eV
               <S^2> =    0.0000
     Summary table
 State      Energy       Excitation   Excitation(eV)  <S^2>         Transition dipole moment, a.u.        Oscillator
   1    -56.5040345278    -7.315353     0.000000      0.000     0.0000     0.0000     0.0000     0.0000      0.0000
   0    -56.2352002518     0.000000     7.315353        (ROHF/UHF Reference state)
   4    -56.1347255488     2.734056    10.049409      0.000    -0.5851     0.0000     0.0000     0.5851      0.0843
"""

        states = module.parse_mrsf_log_state_table(log_text)

        self.assertEqual("S3", states[4]["physical_state"])
        self.assertAlmostEqual(2.734056, states[4]["raw_mrsf_root_value_ev"])
        self.assertAlmostEqual(10.049409, states[4]["physical_excitation_energy_ev"])
        self.assertAlmostEqual(0.0843, states[4]["oscillator_strength"])
        self.assertEqual("ROHF/UHF Reference state", states[0]["state_type"])

    def test_root_continuity_summary_detects_near_degenerate_target_root(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            root_dir = Path(tmpdir)
            for rel, root4_energy in {
                "grad/grad.log": 2.734056,
                "e_a0_y_plus/e_a0_y_plus.log": 2.734051,
                "e_a0_y_minus/e_a0_y_minus.log": 2.734052,
            }.items():
                path = root_dir / rel
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(
                    f"""
  Spin-adapted spin-flip excitations
 State #   1  Energy =   -7.315353 eV
               <S^2> =    0.0000
 State #   3  Energy =    2.734050 eV
               <S^2> =    0.0000
 State #   4  Energy =    {root4_energy:.6f} eV
               <S^2> =    0.0000
     Summary table
 State      Energy       Excitation   Excitation(eV)  <S^2>         Transition dipole moment, a.u.        Oscillator
   1    -56.5040345278    -7.315353     0.000000      0.000     0.0000     0.0000     0.0000     0.0000      0.0000
   0    -56.2352002518     0.000000     7.315353        (ROHF/UHF Reference state)
   3    -56.1347257513     2.734050    10.049404      0.000    -0.0000    -0.5850    -0.0000     0.5850      0.0843
   4    -56.1347255488     {root4_energy:.6f}    10.049409      0.000    -0.5851     0.0000     0.0000     0.5851      0.0843
"""
                )

            summary = module.summarize_root_continuity_dir(root_dir, root=4)

        self.assertEqual(4, summary["root"])
        self.assertEqual("S3", summary["physical_state"])
        self.assertEqual(3, summary["log_count"])
        self.assertEqual("present", summary["s2_evidence"])
        self.assertLess(summary["min_neighbor_gap_ev"], 1.0e-4)
        self.assertTrue(summary["near_degenerate_target"])
        self.assertIn("root-continuity", summary["evidence_hint"])

    def test_cli_root_continuity_mode_writes_log_evidence_summary(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            root_dir = Path(tmpdir)
            log_path = root_dir / "grad" / "grad.log"
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(
                """
  Spin-adapted spin-flip excitations
 State #   1  Energy =   -7.315353 eV
               <S^2> =    0.0000
 State #   4  Energy =    2.734056 eV
               <S^2> =    0.0000
     Summary table
 State      Energy       Excitation   Excitation(eV)  <S^2>         Transition dipole moment, a.u.        Oscillator
   1    -56.5040345278    -7.315353     0.000000      0.000     0.0000     0.0000     0.0000     0.0000      0.0000
   0    -56.2352002518     0.000000     7.315353        (ROHF/UHF Reference state)
   4    -56.1347255488     2.734056    10.049409      0.000    -0.5851     0.0000     0.0000     0.5851      0.0843
"""
            )
            output = root_dir / "root_continuity.json"

            status = module.main(["--root-continuity", "--root", "4", str(root_dir), "--output", str(output)])
            written = output.read_text()

        self.assertEqual(0, status)
        self.assertIn('"physical_state": "S3"', written)
        self.assertIn('"s2_evidence": "present"', written)

    def test_cli_components_mode_accepts_multiple_csvs_for_artifact_ranking(self):
        module = load_module()
        first = self.write_components_csv(
            [
                {
                    "molecule": "ch2o",
                    "method": "mrsf",
                    "root": "4",
                    "physical_state": "",
                    "component": "a1_z",
                    "analytic_ha_per_bohr": "-0.25282498",
                    "fd_ha_per_bohr": "-0.27404",
                    "diff_ha_per_bohr": "0.02121502",
                    "abs_diff_ha_per_bohr": "0.02121502",
                    "trah_count": "0",
                    "failed_any": "False",
                    "s2_grad": "{}",
                    "s2_plus": "{}",
                    "s2_minus": "{}",
                }
            ]
        )
        second = self.write_components_csv(
            [
                {
                    "molecule": "h2s",
                    "method": "mrsf",
                    "root": "5",
                    "physical_state": "",
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
            ]
        )
        output = tempfile.NamedTemporaryFile("r", suffix=".json", delete=False)
        output.close()

        status = module.main(["--components", str(first), str(second), "--output", output.name])

        self.assertEqual(0, status)
        written = Path(output.name).read_text()
        self.assertIn('"dataset_count": 2', written)
        self.assertIn('"target_bad_group_count": 2', written)
        self.assertIn(str(first), written)
        self.assertIn(str(second), written)

    def test_root_continuity_targets_resolve_existing_dirs_and_record_missing_cases(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir) / "mrsf_s6_ovov_fix"
            source_csv = base / "components.csv"
            source_csv.parent.mkdir(parents=True, exist_ok=True)
            source_csv.write_text("placeholder\n")
            root_dir = base / "nh3" / "mrsf" / "root_4"
            log_path = root_dir / "grad" / "grad.log"
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(
                """
  Spin-adapted spin-flip excitations
 State #   1  Energy =   -7.315353 eV
               <S^2> =    0.0000
 State #   4  Energy =    2.734056 eV
               <S^2> =    0.0000
     Summary table
 State      Energy       Excitation   Excitation(eV)  <S^2>         Transition dipole moment, a.u.        Oscillator
   1    -56.5040345278    -7.315353     0.000000      0.000     0.0000     0.0000     0.0000     0.0000      0.0000
   0    -56.2352002518     0.000000     7.315353        (ROHF/UHF Reference state)
   4    -56.1347255488     2.734056    10.049409      0.000    -0.5851     0.0000     0.0000     0.5851      0.0843
"""
            )
            component_summary = Path(tmpdir) / "component_summary.json"
            component_summary.write_text(
                json.dumps(
                    {
                        "target_bad_groups": [
                            {"molecule": "nh3", "method": "mrsf", "root": 4, "source_csv": str(source_csv)},
                            {"molecule": "h2s", "method": "mrsf", "root": 5, "source_csv": str(source_csv)},
                        ]
                    }
                )
            )

            summary = module.summarize_root_continuity_targets(component_summary)

        self.assertEqual(2, summary["case_count"])
        self.assertEqual(1, summary["parsed_case_count"])
        self.assertEqual(1, summary["missing_case_count"])
        self.assertEqual("nh3", summary["cases"][0]["molecule"])
        self.assertEqual("present", summary["cases"][0]["s2_evidence"])
        self.assertEqual("h2s", summary["missing_cases"][0]["molecule"])
        self.assertIn("h2s/mrsf/root_5", summary["missing_cases"][0]["expected_root_dir"])

    def test_source_diagnostic_targets_prioritize_stable_no_trah_cases(self):
        module = load_module()
        component_summary = {
            "target_bad_groups": [
                {
                    "molecule": "nh3",
                    "method": "mrsf",
                    "root": 4,
                    "physical_state": "S3",
                    "source_csv": "artifact/components.csv",
                    "max_abs_diff_ha_per_bohr": 0.092,
                    "worst_component": "a0_y",
                    "bad_component_count": 1,
                    "bad_axes": ["y"],
                    "bad_components": [{"component": "a0_y", "axis": "y", "abs_diff_ha_per_bohr": 0.092}],
                    "mechanism_hint": "root_tracking_or_state_character_change",
                },
                {
                    "molecule": "h2s",
                    "method": "mrsf",
                    "root": 5,
                    "physical_state": "S4",
                    "source_csv": "artifact/components.csv",
                    "max_abs_diff_ha_per_bohr": 0.079,
                    "worst_component": "a0_z",
                    "bad_component_count": 1,
                    "bad_axes": ["z"],
                    "bad_components": [{"component": "a0_z", "axis": "z", "abs_diff_ha_per_bohr": 0.079}],
                    "mechanism_hint": "localized_z_component_z_vector_or_operator_mapping",
                },
                {
                    "molecule": "ch2o",
                    "method": "mrsf",
                    "root": 4,
                    "physical_state": "S3",
                    "source_csv": "artifact/components.csv",
                    "max_abs_diff_ha_per_bohr": 0.021,
                    "worst_component": "a1_z",
                    "bad_component_count": 1,
                    "bad_axes": ["z"],
                    "bad_components": [{"component": "a1_z", "axis": "z", "abs_diff_ha_per_bohr": 0.021}],
                    "mechanism_hint": "localized_z_component_z_vector_or_operator_mapping",
                },
            ]
        }
        continuity_summary = {
            "cases": [
                {"molecule": "nh3", "method": "mrsf", "root": 4, "near_degenerate_target": True, "trah_log_count": 0, "s2_evidence": "present"},
                {"molecule": "h2s", "method": "mrsf", "root": 5, "near_degenerate_target": False, "trah_log_count": 0, "s2_evidence": "present"},
                {"molecule": "ch2o", "method": "mrsf", "root": 4, "near_degenerate_target": False, "trah_log_count": 0, "s2_evidence": "present"},
            ],
            "missing_cases": [],
        }

        summary = module.summarize_source_diagnostic_targets(component_summary, continuity_summary)

        self.assertEqual(2, summary["stable_source_candidate_count"])
        self.assertEqual(1, summary["deferred_root_continuity_risk_count"])
        self.assertEqual("h2s", summary["stable_source_candidates"][0]["molecule"])
        self.assertEqual(5, summary["stable_source_candidates"][0]["root"])
        self.assertEqual("source_diagnostic_candidate", summary["stable_source_candidates"][0]["diagnostic_status"])
        self.assertIn("Z-vector/operator", summary["stable_source_candidates"][0]["recommended_next_check"])
        self.assertEqual("nh3", summary["deferred_root_continuity_risks"][0]["molecule"])

    def test_cli_source_diagnostic_targets_writes_ranked_candidates(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            component_summary = Path(tmpdir) / "component_summary.json"
            continuity_summary = Path(tmpdir) / "continuity_summary.json"
            output = Path(tmpdir) / "source_targets.json"
            component_summary.write_text(
                json.dumps(
                    {
                        "target_bad_groups": [
                            {
                                "molecule": "h2s",
                                "method": "mrsf",
                                "root": 5,
                                "physical_state": "S4",
                                "source_csv": "artifact/components.csv",
                                "max_abs_diff_ha_per_bohr": 0.079,
                                "bad_axes": ["z"],
                                "bad_components": [{"component": "a0_z", "axis": "z", "abs_diff_ha_per_bohr": 0.079}],
                            }
                        ]
                    }
                )
            )
            continuity_summary.write_text(
                json.dumps(
                    {
                        "cases": [
                            {
                                "molecule": "h2s",
                                "method": "mrsf",
                                "root": 5,
                                "source_csv": "artifact/components.csv",
                                "near_degenerate_target": False,
                                "trah_log_count": 0,
                                "s2_evidence": "present",
                            }
                        ],
                        "missing_cases": [],
                    }
                )
            )

            status = module.main([
                "--source-diagnostic-targets",
                str(component_summary),
                str(continuity_summary),
                "--output",
                str(output),
            ])

            written = output.read_text()
        self.assertEqual(0, status)
        self.assertIn('"stable_source_candidate_count": 1', written)
        self.assertIn('"diagnostic_status": "source_diagnostic_candidate"', written)

    def test_source_diagnostic_plan_selects_top_candidate_and_source_guards(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            source_root = Path(tmpdir)
            module_dir = source_root / "source" / "modules"
            module_dir.mkdir(parents=True)
            (module_dir / "tdhf_mrsf_gradient.F90").write_text("call utddft_xc_gradient\n")
            (module_dir / "tdhf_mrsf_z_vector.F90").write_text("call mrsf_z_vector\n")
            source_targets = {
                "stable_source_candidates": [
                    {
                        "molecule": "h2s",
                        "method": "mrsf",
                        "root": 5,
                        "physical_state": "S4",
                        "max_abs_diff_ha_per_bohr": 0.07927826,
                        "bad_axes": ["z"],
                        "bad_components": [
                            {"component": "a0_z", "axis": "z", "abs_diff_ha_per_bohr": 0.07927826}
                        ],
                        "root_dir": "artifact/h2s/mrsf/root_5",
                    }
                ]
            }

            plan = module.summarize_source_diagnostic_plan(source_targets, source_root=source_root)

        self.assertEqual("h2s", plan["selected_candidate"]["molecule"])
        self.assertEqual(5, plan["selected_candidate"]["root"])
        self.assertEqual("S4", plan["selected_candidate"]["physical_state"])
        self.assertEqual("localized_z_component", plan["diagnostic_family"])
        self.assertTrue(all(item["exists"] for item in plan["source_files_to_inspect"]))
        self.assertIn("finite-difference", " ".join(plan["validation_required_before_fix_claim"]))
        self.assertIn("no production algebra edit", plan["scope_guard"])

    def test_source_diagnostic_evidence_records_static_source_signals_without_fix_claims(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            source_root = Path(tmpdir)
            module_dir = source_root / "source" / "modules"
            module_dir.mkdir(parents=True)
            (module_dir / "tdhf_mrsf_gradient.F90").write_text(
                """
call utddft_xc_gradient(basis=basis, &
     molGrid=molGrid, &
     dedft=infos%atoms%grad)
df1 = df1 + sgnk*qfspcp2*db2
spcscale = [infos%tddft%spc_coco, &
            infos%tddft%spc_ovov, &
            infos%tddft%spc_coov]
"""
            )
            (module_dir / "tdhf_mrsf_z_vector.F90").write_text(
                """
call mrsfcbc(infos, mo_a, mo_a, wrk1, fmrst1(1,:,:,:))
fmrst1(1,7,:,:) = td_abxc
td_mrsf_den(1:7,:,:) = fmrst1(1,1:7,:,:)
"""
            )
            plan = {
                "selected_candidate": {
                    "molecule": "h2s",
                    "root": 5,
                    "physical_state": "S4",
                    "bad_components": [{"component": "a0_z", "axis": "z", "abs_diff_ha_per_bohr": 0.079}],
                },
                "diagnostic_family": "localized_z_component",
            }

            evidence = module.summarize_source_diagnostic_evidence(plan, source_root=source_root)

        self.assertEqual("h2s", evidence["selected_candidate"]["molecule"])
        self.assertEqual("static_source_diagnostic_only", evidence["evidence_scope"])
        signals = evidence["source_signals"]
        self.assertTrue(signals["ovov_gradient_sign_uses_post_pr153_plus"])
        self.assertTrue(signals["z_vector_channel7_overwrites_mrsfcbc_with_td_abxc"])
        self.assertTrue(signals["z_vector_mrsfcbc_uses_rohf_same_mo"])
        self.assertTrue(signals["gradient_spcscale_order_present"])
        self.assertTrue(signals["z_vector_td_mrsf_den_consumes_seven_channels"])
        self.assertFalse(signals["gradient_xc_call_has_explicit_xa_xb_handoff"])
        locations = evidence["source_signal_locations"]
        self.assertEqual([5], locations["ovov_gradient_sign_post_pr153_plus_lines"])
        self.assertEqual([2], locations["gradient_xc_call_lines"])
        self.assertEqual([2], locations["z_vector_mrsfcbc_rohf_same_mo_lines"])
        self.assertEqual([3], locations["z_vector_channel7_td_abxc_overwrite_lines"])
        self.assertEqual([6], locations["gradient_spcscale_order_lines"])
        self.assertEqual([4], locations["z_vector_td_mrsf_den_handoff_lines"])
        snippets = evidence["source_signal_snippets"]
        self.assertEqual(
            [{"line": 5, "text": "df1 = df1 + sgnk*qfspcp2*db2"}],
            snippets["ovov_gradient_sign_post_pr153_plus"],
        )
        self.assertEqual(
            [{"line": 3, "text": "fmrst1(1,7,:,:) = td_abxc"}],
            snippets["z_vector_channel7_td_abxc_overwrite"],
        )
        self.assertEqual(
            [{"line": 6, "text": "spcscale = [infos%tddft%spc_coco, &"}],
            snippets["gradient_spcscale_order"],
        )
        self.assertEqual(
            [{"line": 4, "text": "td_mrsf_den(1:7,:,:) = fmrst1(1,1:7,:,:)"}],
            snippets["z_vector_td_mrsf_den_handoff"],
        )
        next_validation = evidence["next_validation_plan"]
        self.assertEqual("h2s", next_validation["molecule"])
        self.assertEqual(5, next_validation["root"])
        self.assertEqual("S4", next_validation["physical_state"])
        self.assertEqual("localized_z_component", next_validation["diagnostic_family"])
        self.assertEqual(["a0_z"], next_validation["components_to_validate"])
        self.assertTrue(next_validation["requires_no_fix_control"])
        self.assertTrue(next_validation["requires_finite_difference_rerun"])
        self.assertIn("no production algebra edit", evidence["scope_guard"])
        self.assertIn("finite-difference", " ".join(evidence["validation_required_before_fix_claim"]))

    def test_cli_source_diagnostic_evidence_writes_static_source_signals(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source_root = tmp / "repo"
            module_dir = source_root / "source" / "modules"
            module_dir.mkdir(parents=True)
            (module_dir / "tdhf_mrsf_gradient.F90").write_text("df1 = df1 + sgnk*qfspcp2*db2\n")
            (module_dir / "tdhf_mrsf_z_vector.F90").write_text("fmrst1(1,7,:,:) = td_abxc\n")
            plan_path = tmp / "plan.json"
            output = tmp / "evidence.json"
            plan_path.write_text(json.dumps({"selected_candidate": {"molecule": "h2s", "root": 5}}))

            status = module.main([
                "--source-diagnostic-evidence",
                str(plan_path),
                "--source-root",
                str(source_root),
                "--output",
                str(output),
            ])
            written = output.read_text()

        self.assertEqual(0, status)
        self.assertIn('"evidence_scope": "static_source_diagnostic_only"', written)
        self.assertIn('"z_vector_channel7_overwrites_mrsfcbc_with_td_abxc": true', written)

    def test_source_validation_manifest_records_existing_evidence_and_missing_controls(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            root_dir = Path(tmpdir) / "h2s" / "mrsf" / "root_5"
            (root_dir / "grad").mkdir(parents=True)
            (root_dir / "e_a0_z_plus").mkdir(parents=True)
            (root_dir / "e_a0_z_minus").mkdir(parents=True)
            (root_dir / "grad" / "grad.log").write_text("no TRAH here\n")
            (root_dir / "e_a0_z_plus" / "plus.log").write_text("no TRAH here\n")
            control_dir = root_dir / "validation_controls"
            control_dir.mkdir()
            (control_dir / "fd_rerun_a0_z_components.csv").write_text("component,abs_diff_ha_per_bohr\n")
            evidence = {
                "next_validation_plan": {
                    "molecule": "h2s",
                    "method": "mrsf",
                    "root": 5,
                    "physical_state": "S4",
                    "diagnostic_family": "localized_z_component",
                    "root_dir": str(root_dir),
                    "components_to_validate": ["a0_z"],
                    "bad_components_to_validate": [
                        {
                            "component": "a0_z",
                            "axis": "z",
                            "abs_diff_ha_per_bohr": 0.07927826,
                            "analytic_ha_per_bohr": -0.14875174,
                            "fd_ha_per_bohr": -0.22803,
                        }
                    ],
                }
            }

            manifest = module.summarize_source_validation_manifest(evidence)

        self.assertEqual("h2s", manifest["selected_case"]["molecule"])
        self.assertEqual(["a0_z"], manifest["components_to_validate"])
        self.assertEqual(
            [
                {
                    "component": "a0_z",
                    "axis": "z",
                    "abs_diff_ha_per_bohr": 0.07927826,
                    "analytic_ha_per_bohr": -0.14875174,
                    "fd_ha_per_bohr": -0.22803,
                }
            ],
            manifest["bad_components_to_validate"],
        )
        self.assertTrue(manifest["existing_evidence"]["root_dir_exists"])
        self.assertEqual(2, manifest["existing_evidence"]["log_count"])
        self.assertFalse(manifest["existing_evidence"]["trah_detected"])
        self.assertIn("finite_difference_rerun", manifest["required_controls"][0]["control"])
        self.assertEqual("partial", manifest["required_controls"][0]["status"])
        self.assertEqual(
            [
                {
                    "component": "a0_z",
                    "fd_component_csv": str(root_dir / "validation_controls" / "fd_rerun_a0_z_components.csv"),
                    "fd_component_csv_exists": True,
                    "fd_summary_json": str(root_dir / "validation_controls" / "fd_rerun_a0_z_summary.json"),
                    "fd_summary_json_exists": False,
                    "no_fix_control_json": str(root_dir / "validation_controls" / "no_fix_a0_z_control.json"),
                    "no_fix_control_json_exists": False,
                }
            ],
            manifest["control_artifact_plan"],
        )
        readiness = manifest["validation_readiness"]
        self.assertFalse(readiness["ready_for_source_edit"])
        self.assertEqual("blocked_missing_or_partial_controls", readiness["status"])
        self.assertIn("finite_difference_rerun_for_selected_components", readiness["blocking_controls"])
        self.assertIn("no_fix_or_pre_change_control_same_case", readiness["blocking_controls"])
        self.assertIn("diagnostic manifest only", manifest["scope_guard"])

    def test_validation_control_inputs_package_existing_inputs_without_running_jobs(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            root_dir = Path(tmpdir) / "h2s" / "mrsf" / "root_5"
            (root_dir / "grad").mkdir(parents=True)
            (root_dir / "e_a0_z_plus").mkdir(parents=True)
            (root_dir / "e_a0_z_minus").mkdir(parents=True)
            (root_dir / "grad" / "grad.inp").write_text("gradient input\n")
            (root_dir / "e_a0_z_plus" / "e_a0_z_plus.inp").write_text("plus input\n")
            (root_dir / "e_a0_z_minus" / "e_a0_z_minus.inp").write_text("minus input\n")
            manifest = {
                "selected_case": {
                    "molecule": "h2s",
                    "method": "mrsf",
                    "root": 5,
                    "physical_state": "S4",
                    "root_dir": str(root_dir),
                },
                "bad_components_to_validate": [
                    {
                        "component": "a0_z",
                        "axis": "z",
                        "abs_diff_ha_per_bohr": 0.07927826,
                        "analytic_ha_per_bohr": -0.14875174,
                        "fd_ha_per_bohr": -0.22803,
                    }
                ],
                "control_artifact_plan": [
                    {
                        "component": "a0_z",
                        "fd_component_csv": str(root_dir / "validation_controls" / "fd_rerun_a0_z_components.csv"),
                        "fd_summary_json": str(root_dir / "validation_controls" / "fd_rerun_a0_z_summary.json"),
                        "no_fix_control_json": str(root_dir / "validation_controls" / "no_fix_a0_z_control.json"),
                    }
                ],
            }

            controls = module.summarize_validation_control_inputs(manifest)

        self.assertEqual("validation_control_inputs_only", controls["control_scope"])
        self.assertFalse(controls["jobs_launched"])
        self.assertEqual("h2s", controls["selected_case"]["molecule"])
        self.assertEqual(1, controls["component_count"])
        component = controls["components"][0]
        self.assertEqual("a0_z", component["component"])
        self.assertEqual("z", component["axis"])
        self.assertEqual(0.07927826, component["abs_diff_ha_per_bohr"])
        self.assertTrue(component["existing_input_files"]["gradient_input_exists"])
        self.assertTrue(component["existing_input_files"]["plus_input_exists"])
        self.assertTrue(component["existing_input_files"]["minus_input_exists"])
        self.assertEqual(
            str(root_dir / "validation_controls" / "fd_rerun_a0_z_components.csv"),
            component["planned_outputs"]["fd_component_csv"],
        )
        self.assertEqual("ready_to_generate_control_scripts", controls["next_action"])
        self.assertIn("no OpenQP jobs launched", controls["scope_guard"])

    def test_cli_validation_control_inputs_writes_no_run_package(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            root_dir = tmp / "h2s" / "mrsf" / "root_5"
            (root_dir / "grad").mkdir(parents=True)
            (root_dir / "e_a0_z_plus").mkdir(parents=True)
            (root_dir / "e_a0_z_minus").mkdir(parents=True)
            (root_dir / "grad" / "grad.inp").write_text("gradient input\n")
            (root_dir / "e_a0_z_plus" / "e_a0_z_plus.inp").write_text("plus input\n")
            (root_dir / "e_a0_z_minus" / "e_a0_z_minus.inp").write_text("minus input\n")
            manifest_path = tmp / "manifest.json"
            output = tmp / "controls.json"
            manifest_path.write_text(json.dumps({
                "selected_case": {"molecule": "h2s", "method": "mrsf", "root": 5, "root_dir": str(root_dir)},
                "bad_components_to_validate": [{"component": "a0_z", "axis": "z"}],
                "control_artifact_plan": [{"component": "a0_z", "fd_component_csv": str(root_dir / "validation_controls" / "fd_rerun_a0_z_components.csv")}],
            }))

            status = module.main(["--validation-control-inputs", str(manifest_path), "--output", str(output)])
            written = output.read_text()

        self.assertEqual(0, status)
        self.assertIn('"control_scope": "validation_control_inputs_only"', written)
        self.assertIn('"jobs_launched": false', written)

    def test_validation_control_scripts_plan_commands_without_launching_jobs(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source_root = tmp / "repo"
            module_dir = source_root / "source" / "modules"
            module_dir.mkdir(parents=True)
            (module_dir / "tdhf_mrsf_gradient.F90").write_text("df1 = df1 + sgnk*qfspcp2*db2\n")
            (module_dir / "tdhf_mrsf_z_vector.F90").write_text("fmrst1(1,7,:,:) = td_abxc\n")
            root_dir = tmp / "h2s" / "mrsf" / "root_5"
            grad_dir = root_dir / "grad"
            plus_dir = root_dir / "e_a0_z_plus"
            minus_dir = root_dir / "e_a0_z_minus"
            grad_dir.mkdir(parents=True)
            plus_dir.mkdir(parents=True)
            minus_dir.mkdir(parents=True)
            controls = {
                "selected_case": {"molecule": "h2s", "method": "mrsf", "root": 5, "physical_state": "S4"},
                "components": [
                    {
                        "component": "a0_z",
                        "axis": "z",
                        "existing_input_files": {
                            "gradient_input": str(grad_dir / "grad.inp"),
                            "gradient_input_exists": True,
                            "plus_input": str(plus_dir / "e_a0_z_plus.inp"),
                            "plus_input_exists": True,
                            "minus_input": str(minus_dir / "e_a0_z_minus.inp"),
                            "minus_input_exists": True,
                        },
                        "planned_outputs": {
                            "fd_component_csv": str(root_dir / "validation_controls" / "fd_rerun_a0_z_components.csv"),
                            "fd_summary_json": str(root_dir / "validation_controls" / "fd_rerun_a0_z_summary.json"),
                            "no_fix_control_json": str(root_dir / "validation_controls" / "no_fix_a0_z_control.json"),
                        },
                    }
                ],
            }

            scripts = module.summarize_validation_control_scripts(controls, source_root=source_root)

        self.assertEqual("validation_control_scripts_plan_only", scripts["control_scope"])
        self.assertFalse(scripts["jobs_launched"])
        self.assertFalse(scripts["scripts_written"])
        self.assertEqual("manual_review_before_launch", scripts["next_action"])
        self.assertFalse(scripts["launch_allowed"])
        self.assertEqual(
            [
                "manual_review_before_launch",
                "finite_difference_jobs_not_started_by_this_planner",
                "no_fix_control_not_started_by_this_planner",
            ],
            scripts["launch_blockers"],
        )
        self.assertIn("confirm current branch/source hash", scripts["manual_review_checklist"])
        self.assertIn("confirm no-fix/pre-change control source ref", scripts["manual_review_checklist"])
        self.assertIn("confirm root-continuity/no-TRAH evidence", scripts["manual_review_checklist"])
        self.assertEqual(
            {
                "manual_review_required": True,
                "approved_to_launch": False,
                "reviewed_by": None,
                "review_note": "review-only plan; no validation-control jobs may be launched from this artifact",
            },
            scripts["manual_review_status"],
        )
        self.assertEqual(1, scripts["component_count"])
        component = scripts["components"][0]
        self.assertEqual("a0_z", component["component"])
        self.assertTrue(component["script_path"].endswith("validation_controls/run_a0_z_controls.sh"))
        self.assertEqual(3, component["command_count"])
        self.assertEqual(
            {
                "launch_allowed": False,
                "jobs_launched": False,
                "execution_status": "not_started_review_only",
            },
            component["execution_guard"],
        )
        self.assertIn("OMP_NUM_THREADS=4", component["commands"][0])
        self.assertIn("openqp --nompi grad.inp", component["commands"][0])
        self.assertIn("openqp --nompi e_a0_z_plus.inp", component["commands"][1])
        self.assertIn("openqp --nompi e_a0_z_minus.inp", component["commands"][2])
        self.assertEqual(
            [
                {
                    "role": "gradient",
                    "input": str(grad_dir / "grad.inp"),
                    "workdir": str(grad_dir),
                    "log": "grad_a0_z_control.log",
                    "thread_count": 4,
                    "command": component["commands"][0],
                },
                {
                    "role": "plus_displacement",
                    "input": str(plus_dir / "e_a0_z_plus.inp"),
                    "workdir": str(plus_dir),
                    "log": "e_a0_z_plus_control.log",
                    "thread_count": 4,
                    "command": component["commands"][1],
                },
                {
                    "role": "minus_displacement",
                    "input": str(minus_dir / "e_a0_z_minus.inp"),
                    "workdir": str(minus_dir),
                    "log": "e_a0_z_minus_control.log",
                    "thread_count": 4,
                    "command": component["commands"][2],
                },
            ],
            component["command_manifest"],
        )
        snapshot = scripts["source_snapshot"]
        self.assertEqual("source_file_hashes_for_manual_review", snapshot["snapshot_scope"])
        self.assertTrue(snapshot["source_files"][0]["exists"])
        self.assertEqual("source/modules/tdhf_mrsf_gradient.F90", snapshot["source_files"][0]["path"])
        self.assertEqual(64, len(snapshot["source_files"][0]["sha256"] or ""))
        self.assertTrue(snapshot["all_source_files_present"])
        self.assertIn("no shell scripts written", scripts["scope_guard"])

    def test_validation_control_results_summarizes_completed_artifacts_without_fix_claim(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root_dir = Path(tmp) / "h2s" / "mrsf" / "root_5"
            control_dir = root_dir / "validation_controls"
            control_dir.mkdir(parents=True)
            fd_csv = control_dir / "fd_rerun_a0_z_components.csv"
            fd_csv.write_text("component,axis,abs_diff_ha_per_bohr\na0_z,z,0.07927826\n")
            fd_summary = control_dir / "fd_rerun_a0_z_summary.json"
            fd_summary.write_text(
                json.dumps(
                    {
                        "component": "a0_z",
                        "jobs_launched": True,
                        "max_abs_diff_ha_per_bohr": 0.07927826,
                        "prior_abs_diff_ha_per_bohr": 0.07927825,
                        "reproduces_prior_residual": True,
                        "production_gradient_algebra_edited": False,
                        "root_continuity_no_trah_status": "present",
                        "trah_detected": False,
                        "next_action": "debug_source_level_z_vector_or_density_handoff",
                    }
                )
            )
            no_fix = control_dir / "no_fix_a0_z_control.json"
            no_fix.write_text(
                json.dumps(
                    {
                        "component": "a0_z",
                        "jobs_launched": True,
                        "production_gradient_algebra_edited": False,
                        "status": "current_baseline_control_recorded",
                    }
                )
            )
            manifest = {
                "selected_case": {
                    "molecule": "h2s",
                    "method": "mrsf",
                    "root": 5,
                    "physical_state": "S4",
                    "root_dir": str(root_dir),
                },
                "components_to_validate": ["a0_z"],
                "control_artifact_plan": [
                    {
                        "component": "a0_z",
                        "fd_component_csv": str(fd_csv),
                        "fd_summary_json": str(fd_summary),
                        "no_fix_control_json": str(no_fix),
                    }
                ],
                "validation_readiness": {
                    "ready_for_source_edit": True,
                    "status": "ready_for_source_edit",
                    "blocking_controls": [],
                },
            }

            results = module.summarize_validation_control_results(manifest)

        self.assertEqual("validation_control_results_summary", results["control_scope"])
        self.assertEqual("h2s root 5 / physical S4", results["selected"])
        self.assertEqual(1, results["component_count"])
        self.assertEqual(0, results["missing_artifact_count"])
        self.assertTrue(results["controls_complete"])
        self.assertFalse(results["ready_for_production_fix_claim"])
        self.assertEqual("source_diagnostic_hypothesis_required_before_source_edit", results["next_action"])
        component = results["components"][0]
        self.assertEqual("a0_z", component["component"])
        self.assertTrue(component["fd_rerun_present"])
        self.assertTrue(component["no_fix_control_present"])
        self.assertTrue(component["reproduces_prior_residual"])
        self.assertEqual(0.07927826, component["max_abs_diff_ha_per_bohr"])
        self.assertEqual("current_baseline_control_recorded", component["no_fix_control_status"])
        self.assertIn("not a production fix claim", results["scope_guard"])

    def test_source_hypothesis_diagnostic_records_ranked_hypotheses_without_source_edit(self):
        module = load_module()
        source_evidence = {
            "selected_candidate": {
                "molecule": "h2s",
                "method": "mrsf",
                "root": 5,
                "physical_state": "S4",
            },
            "diagnostic_family": "localized_z_component",
            "source_signals": {
                "z_vector_channel7_overwrites_mrsfcbc_with_td_abxc": True,
                "gradient_xc_call_has_explicit_xa_xb_handoff": False,
                "ovov_gradient_sign_uses_post_pr153_plus": True,
            },
            "source_signal_locations": {
                "z_vector_channel7_td_abxc_overwrite_lines": [842],
                "gradient_xc_call_lines": [156],
            },
            "source_signal_snippets": {
                "z_vector_channel7_td_abxc_overwrite": [{"line": 842, "text": "fmrst1(1,7,:,:) = td_abxc"}],
                "gradient_xc_call": [{"line": 156, "text": "call utddft_xc_gradient(...)"}],
            },
        }
        validation_results = {
            "selected": "h2s root 5 / physical S4",
            "ready_for_source_hypothesis_diagnostic": True,
            "controls_complete": True,
            "trah_detected": False,
            "ready_for_production_fix_claim": False,
            "components": [
                {
                    "component": "a0_z",
                    "max_abs_diff_ha_per_bohr": 0.07927826,
                    "reproduces_prior_residual": True,
                    "root_continuity_no_trah_status": "present",
                }
            ],
        }

        diagnostic = module.summarize_source_hypothesis_diagnostic(source_evidence, validation_results)

        self.assertEqual("source_hypothesis_diagnostic_only", diagnostic["diagnostic_scope"])
        self.assertEqual("h2s root 5 / physical S4", diagnostic["selected"])
        self.assertTrue(diagnostic["controls_complete"])
        self.assertFalse(diagnostic["production_gradient_algebra_edited"])
        self.assertFalse(diagnostic["ready_for_production_fix_claim"])
        self.assertEqual("a0_z", diagnostic["residual_components"][0]["component"])
        self.assertEqual("channel7_density_provenance", diagnostic["ranked_source_hypotheses"][0]["hypothesis_id"])
        self.assertIn("fmrst1(1,7", diagnostic["ranked_source_hypotheses"][0]["source_snippets"][0]["text"])
        self.assertEqual("source_hypothesis_validation_required", diagnostic["next_action"])
        self.assertIn("no production algebra edit", diagnostic["scope_guard"])

    def test_source_level_validation_compares_channel7_provenance_before_xc_handoff(self):
        module = load_module()
        source_hypothesis = {
            "selected": "h2s root 5 / physical S4",
            "controls_complete": True,
            "trah_detected": False,
            "ready_for_production_fix_claim": False,
            "residual_components": [
                {
                    "component": "a0_z",
                    "max_abs_diff_ha_per_bohr": 0.07927826,
                    "reproduces_prior_residual": True,
                    "root_continuity_no_trah_status": "present",
                }
            ],
            "ranked_source_hypotheses": [
                {
                    "hypothesis_id": "channel7_density_provenance",
                    "rank": 1,
                    "source_locations": [842],
                    "source_snippets": [{"line": 842, "text": "fmrst1(1,7,:,:) = td_abxc"}],
                },
                {
                    "hypothesis_id": "mrsf_xc_density_handoff",
                    "rank": 2,
                    "source_locations": [156],
                    "source_snippets": [{"line": 156, "text": "call utddft_xc_gradient(basis=basis, &"}],
                },
            ],
        }

        validation = module.summarize_source_level_validation(source_hypothesis)

        self.assertEqual("source_level_validation_only", validation["validation_scope"])
        self.assertEqual("h2s root 5 / physical S4", validation["selected"])
        self.assertEqual(["channel7_density_provenance", "mrsf_xc_density_handoff"], validation["compared_hypotheses"])
        self.assertEqual("channel7_density_provenance", validation["primary_next_source_test"])
        self.assertEqual("mrsf_xc_density_handoff", validation["secondary_next_source_test"])
        self.assertEqual("validated_for_source_trial", validation["source_level_validation_status"])
        self.assertFalse(validation["production_gradient_algebra_edited"])
        self.assertFalse(validation["ready_for_production_fix_claim"])
        self.assertIn("fmrst1(1,7", validation["validation_points"][0]["source_snippets"][0]["text"])
        self.assertIn("do not pass td_abxc directly", validation["validation_points"][1]["guardrail"])
        self.assertEqual("prepare_one_variable_channel7_source_trial_or_instrumentation", validation["next_action"])
        self.assertIn("diagnostic-only", validation["scope_guard"])

    def test_channel7_source_trial_plan_is_review_only_and_one_variable(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            source_root = Path(tmp)
            gradient = source_root / "source" / "modules" / "tdhf_mrsf_gradient.F90"
            z_vector = source_root / "source" / "modules" / "tdhf_mrsf_z_vector.F90"
            gradient.parent.mkdir(parents=True)
            gradient.write_text("call utddft_xc_gradient(basis=basis, &\n")
            z_vector.write_text("call mrsfcbc(..., mo_a, mo_a, ball)\nfmrst1(1,7,:,:) = td_abxc\n")
            source_level_validation = {
                "selected": "h2s root 5 / physical S4",
                "source_level_validation_status": "validated_for_source_trial",
                "primary_next_source_test": "channel7_density_provenance",
                "secondary_next_source_test": "mrsf_xc_density_handoff",
                "residual_components": [
                    {
                        "component": "a0_z",
                        "max_abs_diff_ha_per_bohr": 0.07927826,
                        "reproduces_prior_residual": True,
                        "root_continuity_no_trah_status": "present",
                    }
                ],
                "validation_points": [
                    {
                        "hypothesis_id": "channel7_density_provenance",
                        "source_locations": [842],
                        "source_snippets": [{"line": 842, "text": "fmrst1(1,7,:,:) = td_abxc"}],
                    },
                    {
                        "hypothesis_id": "mrsf_xc_density_handoff",
                        "source_locations": [156],
                        "source_snippets": [{"line": 156, "text": "call utddft_xc_gradient(basis=basis, &"}],
                    },
                ],
            }

            plan = module.summarize_channel7_source_trial_plan(source_level_validation, source_root=source_root)

        self.assertEqual("channel7_source_trial_plan_only", plan["trial_scope"])
        self.assertEqual("h2s root 5 / physical S4", plan["selected"])
        self.assertEqual("channel7_density_provenance", plan["one_variable_under_test"])
        self.assertTrue(plan["ready_for_manual_source_trial"])
        self.assertFalse(plan["production_gradient_algebra_edited"])
        self.assertFalse(plan["source_files_modified_by_planner"])
        self.assertFalse(plan["ready_for_production_fix_claim"])
        self.assertIn("mrsf_xc_density_handoff", plan["forbidden_bundled_changes"])
        self.assertIn("ovov_sign_baseline_control", plan["forbidden_bundled_changes"])
        self.assertIn("a0_z", plan["residual_components_to_recheck"])
        self.assertEqual("source_file_hashes_for_manual_review", plan["source_snapshot"]["snapshot_scope"])
        self.assertTrue(plan["source_snapshot"]["all_source_files_present"])
        self.assertEqual("manual_review_before_source_edit", plan["manual_review_status"]["next_gate"])
        self.assertEqual("run_one_variable_channel7_trial_then_repeat_fd_control", plan["next_action"])
        self.assertIn("no source files modified", plan["scope_guard"])

    def test_source_trial_outcome_defers_channel7_after_negative_result(self):
        module = load_module()
        source_level_validation = {
            "selected": "h2s root 5 / physical S4",
            "primary_next_source_test": "channel7_density_provenance",
            "secondary_next_source_test": "mrsf_xc_density_handoff",
            "validation_points": [
                {"hypothesis_id": "channel7_density_provenance", "source_locations": [842]},
                {"hypothesis_id": "mrsf_xc_density_handoff", "source_locations": [156]},
            ],
        }
        trial_results = {
            "trial_scope": "channel7_source_trial_results",
            "selected": "h2s root 5 / physical S4",
            "one_variable_under_test": "channel7_density_provenance",
            "component": "a0_z",
            "abs_diff_ha_per_bohr": 0.07927826468642385,
            "control_no_fix_abs_diff_ha_per_bohr": 0.07927826468642385,
            "delta_vs_no_fix_abs_diff_ha_per_bohr": 0.0,
            "moved_toward_fd_control": False,
            "residual_removed": False,
            "trah_detected": False,
            "production_gradient_algebra_edited": False,
        }

        outcome = module.summarize_source_trial_outcome(source_level_validation, trial_results)

        self.assertEqual("source_trial_outcome_diagnostic_only", outcome["outcome_scope"])
        self.assertEqual("h2s root 5 / physical S4", outcome["selected"])
        self.assertEqual("channel7_density_provenance", outcome["completed_source_test"])
        self.assertEqual("negative_no_change", outcome["completed_source_test_status"])
        self.assertEqual("mrsf_xc_density_handoff", outcome["next_source_test"])
        self.assertIn("channel7_density_provenance", outcome["deferred_hypotheses"])
        self.assertFalse(outcome["production_gradient_algebra_edited"])
        self.assertFalse(outcome["ready_for_production_fix_claim"])
        self.assertEqual("plan_mrsf_xc_density_handoff_diagnostic", outcome["next_action"])
        self.assertIn("diagnostic-only", outcome["scope_guard"])

    def test_xc_density_handoff_diagnostic_plans_safe_static_run_after_negative_channel7(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            source_root = Path(tmpdir)
            gradient = source_root / "source/modules/tdhf_mrsf_gradient.F90"
            z_vector = source_root / "source/modules/tdhf_mrsf_z_vector.F90"
            gradient.parent.mkdir(parents=True, exist_ok=True)
            z_vector.parent.mkdir(parents=True, exist_ok=True)
            gradient.write_text(
                "call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)\n"
                "call tagarray_get_data(infos%dat, OQP_td_p, td_p)\n"
                "call tagarray_get_data(infos%dat, OQP_td_mrsf_density, td_mrsf_density)\n"
                "call utddft_xc_gradient(basis=basis, da=d(:,:,1), db=d(:,:,2), pa=p(:,:,1:1), pb=p(:,:,2:2), infos=infos)\n"
            )
            z_vector.write_text("call infos%dat%reserve_data(OQP_td_mrsf_density, TA_TYPE_REAL64, nbf*nbf*7)\n")
            source_level_validation = {
                "selected": "h2s root 5 / physical S4",
                "secondary_next_source_test": "mrsf_xc_density_handoff",
                "validation_points": [
                    {
                        "hypothesis_id": "mrsf_xc_density_handoff",
                        "source_locations": [4],
                        "source_snippets": [{"line": 4, "text": "call utddft_xc_gradient(... pa=p, pb=p ...)"}],
                    }
                ],
                "residual_components": [{"component": "a0_z", "max_abs_diff_ha_per_bohr": 0.07927826468642385}],
            }
            source_trial_outcome = {
                "completed_source_test": "channel7_density_provenance",
                "completed_source_test_status": "negative_no_change",
                "next_source_test": "mrsf_xc_density_handoff",
                "deferred_hypotheses": ["channel7_density_provenance"],
                "trah_detected": False,
            }

            diagnostic = module.summarize_mrsf_xc_density_handoff_diagnostic(
                source_level_validation,
                source_trial_outcome,
                source_root=source_root,
            )

        self.assertEqual("mrsf_xc_density_handoff_diagnostic_only", diagnostic["diagnostic_scope"])
        self.assertEqual("h2s root 5 / physical S4", diagnostic["selected"])
        self.assertEqual("mrsf_xc_density_handoff", diagnostic["one_variable_under_test"])
        self.assertEqual("static_source_trace_completed", diagnostic["diagnostic_run_status"])
        self.assertTrue(diagnostic["channel7_density_provenance_deferred"])
        self.assertFalse(diagnostic["quantum_jobs_launched"])
        self.assertFalse(diagnostic["source_files_modified_by_diagnostic"])
        self.assertFalse(diagnostic["production_gradient_algebra_edited"])
        self.assertFalse(diagnostic["ready_for_production_fix_claim"])
        self.assertIn("direct_td_abxc_as_xa_xb", diagnostic["forbidden_source_trials"])
        self.assertEqual("instrument_mrsf_spin_resolved_xc_density_candidate", diagnostic["next_action"])
        self.assertEqual(["a0_z"], diagnostic["residual_components_to_recheck"])
        self.assertIn("OQP_td_mrsf_density", diagnostic["source_observations"]["mrsf_density_tags_seen"])
        self.assertTrue(diagnostic["source_snapshot"]["all_source_files_present"])

    def test_xc_density_instrumentation_plan_is_review_only_and_requires_real_candidate(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            source_root = Path(tmpdir)
            gradient = source_root / "source/modules/tdhf_mrsf_gradient.F90"
            z_vector = source_root / "source/modules/tdhf_mrsf_z_vector.F90"
            gradient.parent.mkdir(parents=True, exist_ok=True)
            z_vector.parent.mkdir(parents=True, exist_ok=True)
            gradient.write_text(
                "call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)\n"
                "call tagarray_get_data(infos%dat, OQP_td_p, td_p)\n"
                "call tagarray_get_data(infos%dat, OQP_td_mrsf_density, td_mrsf_density)\n"
                "call utddft_xc_gradient(basis=basis, da=d(:,:,1), db=d(:,:,2), pa=p(:,:,1:1), pb=p(:,:,2:2), infos=infos)\n"
            )
            z_vector.write_text("call infos%dat%reserve_data(OQP_td_mrsf_density, TA_TYPE_REAL64, nbf*nbf*7)\n")
            diagnostic = {
                "diagnostic_scope": "mrsf_xc_density_handoff_diagnostic_only",
                "selected": "h2s root 5 / physical S4",
                "diagnostic_run_status": "static_source_trace_completed",
                "one_variable_under_test": "mrsf_xc_density_handoff",
                "residual_components_to_recheck": ["a0_z"],
                "source_observations": {
                    "gradient_xc_call_has_explicit_xa_xb_handoff": False,
                    "gradient_xc_call_lines": [4],
                    "mrsf_density_tags_seen": ["OQP_td_abxc", "OQP_td_p", "OQP_td_mrsf_density"],
                },
                "source_trial_outcome_status": "negative_no_change",
                "channel7_density_provenance_deferred": True,
            }

            plan = module.summarize_mrsf_xc_density_instrumentation_plan(diagnostic, source_root=source_root)

        self.assertEqual("mrsf_xc_density_instrumentation_plan_only", plan["instrumentation_scope"])
        self.assertEqual("h2s root 5 / physical S4", plan["selected"])
        self.assertEqual("mrsf_xc_density_handoff", plan["one_variable_under_test"])
        self.assertTrue(plan["ready_for_manual_instrumentation_review"])
        self.assertFalse(plan["source_files_modified_by_planner"])
        self.assertFalse(plan["jobs_launched"])
        self.assertFalse(plan["production_gradient_algebra_edited"])
        self.assertFalse(plan["ready_for_production_fix_claim"])
        self.assertEqual(["a0_z"], plan["residual_components_to_recheck"])
        self.assertIn("manual_review_before_source_instrumentation", plan["manual_review_status"]["blockers"])
        self.assertIn("direct_td_abxc_as_xa_xb", plan["forbidden_bundled_changes"])
        self.assertEqual("review_mrsf_spin_resolved_xc_density_candidate_before_source_trial", plan["next_action"])
        self.assertEqual("missing_real_spin_resolved_mrsf_xc_density_candidate", plan["candidate_density_requirements"]["current_blocker"])
        self.assertTrue(plan["source_snapshot"]["all_source_files_present"])


if __name__ == "__main__":
    unittest.main()
