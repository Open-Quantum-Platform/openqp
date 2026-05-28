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


if __name__ == "__main__":
    unittest.main()
