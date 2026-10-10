import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# A real MRSF-TDDFT run log (examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_ENERGY.inp)
# captured from a built OpenQP, including the OQP_TIMER lines emitted by the
# Fortran log_oqp_timer subroutine. This pins the parser to the *actual* emitter
# output, which a Python-formatter round-trip cannot do on its own.
MRSF_LOG_FIXTURE = ROOT / "tests" / "fixtures" / "mrsf_tddft_h2o_energy.log"


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class TDHFBenchmarkTimerTests(unittest.TestCase):
    def setUp(self):
        self.timers = load_module("tdhf_benchmark_timers_under_test", "python/openqp_gpu/tdhf_benchmark_timers.py")

    def test_davidson_timer_manifest_has_stable_gpu_benchmark_labels(self):
        manifest = self.timers.davidson_timer_manifest()
        labels = [timer.label for timer in manifest]

        self.assertEqual(labels[0], "tdhf.response.total")
        self.assertIn("tdhf.davidson.total", labels)
        self.assertIn("tdhf.davidson.sigma_build", labels)
        self.assertIn("tdhf.davidson.metc_contract", labels)
        self.assertIn("tdhf.davidson.eri_buffer", labels)
        self.assertTrue(all(timer.private_gpu_group for timer in manifest))

    def test_format_timer_line_is_machine_parseable_and_unit_explicit(self):
        line = self.timers.format_timer_line(
            "tdhf.davidson.sigma_build",
            elapsed_seconds=1.25,
            metadata={"branch": "perf/tdhf-davidson-timers", "nstate": 4},
        )

        self.assertEqual(
            line,
            "OQP_TIMER label=tdhf.davidson.sigma_build seconds=1.250000 branch=perf/tdhf-davidson-timers nstate=4",
        )

    def test_parse_timer_line_round_trips_formatted_timer(self):
        line = self.timers.format_timer_line(
            "tdhf.davidson.metc_contract",
            elapsed_seconds=0.0032,
            metadata={"kernel": "cpu_baseline"},
        )
        parsed = self.timers.parse_timer_line(line)

        self.assertEqual(parsed["label"], "tdhf.davidson.metc_contract")
        self.assertAlmostEqual(parsed["seconds"], 0.0032)
        self.assertEqual(parsed["kernel"], "cpu_baseline")

    def test_parse_timer_lines_extracts_only_timer_records_from_log_text(self):
        log_text = "\n".join(
            [
                "normal OpenQP output",
                self.timers.format_timer_line("tdhf.davidson.total", 2.5, {"iter": 3}),
                "unrelated warning line",
                self.timers.format_timer_line("tdhf.response.total", 3.0),
            ]
        )

        records = self.timers.parse_timer_lines(log_text)

        self.assertEqual([record["label"] for record in records], ["tdhf.davidson.total", "tdhf.response.total"])
        self.assertEqual(records[0]["iter"], "3")
        self.assertAlmostEqual(records[1]["seconds"], 3.0)

    def test_summarize_timer_records_groups_counts_and_total_seconds_by_label(self):
        records = [
            self.timers.parse_timer_line(self.timers.format_timer_line("tdhf.davidson.total", 2.5, {"iter": 1})),
            self.timers.parse_timer_line(self.timers.format_timer_line("tdhf.davidson.total", 3.5, {"iter": 2})),
            self.timers.parse_timer_line(self.timers.format_timer_line("tdhf.response.total", 9.0)),
        ]

        summary = self.timers.summarize_timer_records(records)

        self.assertEqual(summary["tdhf.davidson.total"], {"count": 2, "seconds_total": 6.0, "seconds_mean": 3.0})
        self.assertEqual(summary["tdhf.response.total"], {"count": 1, "seconds_total": 9.0, "seconds_mean": 9.0})

    def test_format_timer_summary_csv_emits_stable_rows_for_overleaf_data_snapshots(self):
        summary = {
            "tdhf.davidson.total": {"count": 2, "seconds_total": 6.0, "seconds_mean": 3.0},
            "tdhf.response.total": {"count": 1, "seconds_total": 9.0, "seconds_mean": 9.0},
        }

        csv_text = self.timers.format_timer_summary_csv(summary)

        self.assertEqual(
            csv_text,
            "label,count,seconds_total,seconds_mean\n"
            "tdhf.response.total,1,9.000000,9.000000\n"
            "tdhf.davidson.total,2,6.000000,3.000000\n",
        )

    def test_parses_oqp_timer_lines_from_real_mrsf_run_log(self):
        # Regression: the Fortran emitter must produce lines the parser accepts.
        # A fixed-width seconds descriptor would pad with blanks
        # ("seconds=    0.001884") and split into two tokens, so this asserts on
        # a genuine run log rather than a Python-formatted line.
        log_text = MRSF_LOG_FIXTURE.read_text()

        records = self.timers.parse_timer_lines(log_text)
        labels = [record["label"] for record in records]

        # H2O MRSF-TDDFT energy converges in 3 Davidson iterations: one
        # sigma_build + metc_contract per iteration, plus one davidson.total
        # and one response.total at the end.
        self.assertEqual(labels.count("tdhf.davidson.sigma_build"), 3)
        self.assertEqual(labels.count("tdhf.davidson.metc_contract"), 3)
        self.assertEqual(labels.count("tdhf.davidson.total"), 1)
        self.assertEqual(labels.count("tdhf.response.total"), 1)

        # Every seconds value parsed cleanly to a positive float.
        self.assertTrue(all(isinstance(r["seconds"], float) and r["seconds"] >= 0.0 for r in records))

        summary = self.timers.summarize_timer_records(records)
        self.assertEqual(summary["tdhf.davidson.sigma_build"]["count"], 3)
        # response.total brackets the whole solve, so it is the largest bucket.
        self.assertGreaterEqual(
            summary["tdhf.response.total"]["seconds_total"],
            summary["tdhf.davidson.total"]["seconds_total"],
        )

    def test_parses_timer_lines_with_slurm_label_prefix(self):
        # Launchers (e.g. `srun --label`) may prepend a rank tag to each stdout
        # line. The parser must still recover the records from a tee'd log.
        log_text = MRSF_LOG_FIXTURE.read_text()
        prefixed = "\n".join(
            f"0: {line}" if "OQP_TIMER" in line.split() else line
            for line in log_text.splitlines()
        )

        plain_records = self.timers.parse_timer_lines(log_text)
        prefixed_records = self.timers.parse_timer_lines(prefixed)

        self.assertEqual(len(prefixed_records), len(plain_records))
        self.assertEqual(
            [r["label"] for r in prefixed_records],
            [r["label"] for r in plain_records],
        )
        # A prose line that merely mentions OQP_TIMER as a substring must not parse.
        self.assertEqual(self.timers.parse_timer_lines("see OQP_TIMER_NOTES below"), [])


if __name__ == "__main__":
    unittest.main()
