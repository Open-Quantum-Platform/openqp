import importlib.util
import os
import sys
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]


def load_oqp_tester_with_stubs():
    oqp_stub = types.ModuleType("oqp")
    pyoqp_stub = types.ModuleType("oqp.pyoqp")

    class RunnerShouldNotRun:
        def __init__(self, *args, **kwargs):
            raise AssertionError("Runner should not be constructed for skipped optional-backend tests")

    setattr(pyoqp_stub, "Runner", RunnerShouldNotRun)
    sys.modules["oqp"] = oqp_stub
    sys.modules["oqp.pyoqp"] = pyoqp_stub
    spec = importlib.util.spec_from_file_location(
        "oqp_tester_optional_backend_under_test",
        ROOT / "pyoqp/oqp/utils/oqp_tester.py",
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class MPIStub:
    use_mpi = False
    rank = 0
    size = 1


class OQPTesterOptionalBackendTests(unittest.TestCase):
    def test_dftbplus_examples_are_skipped_when_executable_is_missing(self):
        tester_mod = load_oqp_tester_with_stubs()
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            inp = tmp_path / "missing_dftb.inp"
            inp.write_text(
                """
[input]
method=dftb
runtype=energy
system=
 H 0 0 0
 H 0 0 0.74

[dftb]
executable=/definitely/missing/dftb+
sk_path=/tmp/sk
""".strip()
            )
            tester = tester_mod.OQPTester(
                base_test_dir=str(tmp_path),
                output_dir=str(tmp_path / "out"),
                total_cpus=1,
                omp_threads=1,
                mpi_manager=MPIStub(),
            )
            result = tester.run_single_test(str(inp))

        self.assertEqual(result["status"], "SKIPPED")
        self.assertIn("DFTB+ executable not found", result["message"])
        self.assertIn("install DFTB+", result["message"])

    def test_report_counts_skipped_tests_without_failing_suite(self):
        tester_mod = load_oqp_tester_with_stubs()
        with TemporaryDirectory() as tmp:
            tester = tester_mod.OQPTester(
                base_test_dir=tmp,
                output_dir=str(Path(tmp) / "out"),
                total_cpus=1,
                omp_threads=1,
                mpi_manager=MPIStub(),
            )
            tester.start_time = 1.0
            tester.end_time = 2.0
            tester.results = [{
                "project": "missing_dftb",
                "input_file": str(Path(tmp) / "missing_dftb.inp"),
                "log_file": str(Path(tmp) / "missing_dftb.log"),
                "status": "SKIPPED",
                "message": "DFTB+ executable not found; install DFTB+ for full test coverage.",
                "execution_time": 0.0,
            }]
            Path(tester.output_dir).mkdir(parents=True, exist_ok=True)
            report = tester.generate_report()

        self.assertIn("Skipped: 1", report)
        self.assertEqual(tester.status, 0)
    def test_invalid_test_path_generates_empty_report_instead_of_type_error(self):
        tester_mod = load_oqp_tester_with_stubs()
        with TemporaryDirectory() as tmp:
            tester = tester_mod.OQPTester(
                base_test_dir=tmp,
                output_dir=str(Path(tmp) / "out"),
                total_cpus=1,
                omp_threads=1,
                mpi_manager=MPIStub(),
            )
            report = tester.run("missing-test-name")

        self.assertIn("Total tests: 0", report)
        self.assertEqual(tester.status, 0)

    def test_harness_exception_still_generates_report(self):
        tester_mod = load_oqp_tester_with_stubs()
        with TemporaryDirectory() as tmp:
            tester = tester_mod.OQPTester(
                base_test_dir=tmp,
                output_dir=str(Path(tmp) / "out"),
                total_cpus=1,
                omp_threads=1,
                mpi_manager=MPIStub(),
            )

            def fail_run_tests(_test_path):
                raise RuntimeError("worker pool failed")

            tester.run_tests = fail_run_tests
            report = tester.run("all")
            report_path = Path(tester.output_dir) / "test_report.txt"
            detailed_report = report_path.read_text()

        self.assertIn("Errors: 1", report)
        self.assertIn("PyOQP test harness error: RuntimeError - worker pool failed", detailed_report)
        self.assertEqual(tester.status, 1)

    def test_run_tests_all_defaults_to_one_worker_for_stability(self):
        tester_mod = load_oqp_tester_with_stubs()
        with TemporaryDirectory() as tmp:
            tester = tester_mod.OQPTester(
                base_test_dir=tmp,
                output_dir=str(Path(tmp) / "out"),
                total_cpus=8,
                omp_threads=2,
                mpi_manager=MPIStub(),
            )
            self.assertEqual(tester.max_workers, 4)
            tester.run_tests("all")

        self.assertEqual(tester.max_workers, 1)

    def test_run_tests_all_respects_worker_override(self):
        tester_mod = load_oqp_tester_with_stubs()
        with TemporaryDirectory() as tmp:
            tester = tester_mod.OQPTester(
                base_test_dir=tmp,
                output_dir=str(Path(tmp) / "out"),
                total_cpus=8,
                omp_threads=2,
                mpi_manager=MPIStub(),
            )
            old_value = os.environ.get("OQP_TEST_MAX_WORKERS")
            os.environ["OQP_TEST_MAX_WORKERS"] = "3"
            try:
                tester.run_tests("all")
            finally:
                if old_value is None:
                    os.environ.pop("OQP_TEST_MAX_WORKERS", None)
                else:
                    os.environ["OQP_TEST_MAX_WORKERS"] = old_value

        self.assertEqual(tester.max_workers, 3)


if __name__ == "__main__":
    unittest.main()
