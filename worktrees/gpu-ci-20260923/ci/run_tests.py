"""Run explicit CI suites; missing dependencies/devices must never pass as skips."""
import argparse
import os
from pathlib import Path
import sys
import unittest

import xmlrunner

ROOT = Path(__file__).resolve().parents[1]
PYTHON_TESTS = (
    "test_gpu_workspace.py",
    "test_gpu_workspace_legacy.py",
    "test_gpu_metc_persistent_buffers.py",
    "test_gpu_package_import.py",
    "test_gpu_metc_contraction_equivalence.py",
    "test_gpu_benchmark_timers.py",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", choices=("host", "python", "cuda"))
    args = parser.parse_args()
    os.chdir(ROOT)
    patterns = PYTHON_TESTS
    if args.suite == "host":
        patterns += ("test_gpu_workspace_runtime.py",)
    elif args.suite == "cuda":
        os.environ["OPENQP_GPU_METC_REQUIRE"] = "1"
        patterns = ("test_gpu_metc_regression.py",)
    suite = unittest.TestSuite()
    for pattern in patterns:
        tests = unittest.defaultTestLoader.discover(str(ROOT / "tests"), pattern=pattern)
        if not tests.countTestCases():
            raise RuntimeError(f"No tests collected from {pattern}")
        suite.addTests(tests)
    result = xmlrunner.XMLTestRunner(output="reports", verbosity=2).run(suite)
    if result.skipped:
        print(f"CI forbids skipped tests: {result.skipped}", file=sys.stderr)
    return 0 if result.wasSuccessful() and result.testsRun and not result.skipped else 1


if __name__ == "__main__":
    sys.exit(main())
