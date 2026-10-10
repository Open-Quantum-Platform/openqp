"""Validation of schema-backed explicit DFT-D4 damping parameters."""

import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_input_checker():
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        use_mpi = False
        size = 1

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    name = "input_checker_d4_under_test"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/utils/input_checker.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class D4InputCheckerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.checker = load_input_checker()

    def check(self, section):
        report = self.checker.CheckReport()
        self.checker._check_d4({"d4": section}, report)
        return report

    def test_functional_defaults_require_no_explicit_values(self):
        self.assertTrue(self.check({}).ok)

    def test_complete_finite_parameter_set_is_accepted(self):
        report = self.check({
            "s6": "1.0", "s8": "0.95948085", "s9": "1.0",
            "a1": "0.38574991", "a2": "4.80688534", "alp": "16.0",
        })
        self.assertTrue(report.ok, report.to_text())

    def test_partial_or_nonfinite_parameter_sets_are_rejected(self):
        partial = self.check({"s6": "1.0"})
        self.assertFalse(partial.ok)
        self.assertIn("all six", partial.to_text())

        nonfinite = self.check({
            "s6": "1.0", "s8": "nan", "s9": "1.0",
            "a1": "0.38574991", "a2": "4.80688534", "alp": "16.0",
        })
        self.assertFalse(nonfinite.ok)
        self.assertIn("d4.s8", nonfinite.to_text())


if __name__ == "__main__":
    unittest.main()
