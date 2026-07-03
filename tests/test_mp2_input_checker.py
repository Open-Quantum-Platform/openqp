import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_input_checker(name):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/utils/input_checker.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_minimal_oqp_stubs():
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        use_mpi = False
        size = 1

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils


class TestMP2InputChecker(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()

    def test_mp2_rejects_dft_functional(self):
        input_checker = load_input_checker("input_checker_mp2_functional_under_test")
        config = {
            "input": {"method": "mp2", "functional": "pbe"},
            "mp2": {"variant": "mp2"},
        }

        report = input_checker.CheckReport()
        input_checker._check_mp2(config, report)

        self.assertFalse(report.ok)
        self.assertIn("input.functional", report.to_text())

    def test_mp2_rejects_derivative_runtype(self):
        input_checker = load_input_checker("input_checker_mp2_runtype_under_test")
        config = {
            "input": {"method": "mp2", "runtype": "optimize"},
        }

        report = input_checker.CheckReport()
        input_checker._check_runtype(config, report)

        self.assertFalse(report.ok)
        self.assertIn("energy-only", report.to_text())

    def test_mp2_accepts_named_and_custom_variants(self):
        input_checker = load_input_checker("input_checker_mp2_variants_under_test")
        for variant in ("mp2", "scs-mp2", "sos-mp2", "os-mp2", "ss-mp2", "scs-mi-mp2"):
            with self.subTest(variant=variant):
                report = input_checker.CheckReport()
                input_checker._check_mp2(
                    {"input": {"method": "mp2"}, "mp2": {"variant": variant}},
                    report,
                )
                self.assertTrue(report.ok, report.to_text())

        report = input_checker.CheckReport()
        input_checker._check_mp2(
            {
                "input": {"method": "mp2"},
                "mp2": {
                    "variant": "custom",
                    "same_spin_scale": 0.5,
                    "opposite_spin_scale": 1.1,
                },
            },
            report,
        )
        self.assertTrue(report.ok, report.to_text())


if __name__ == "__main__":
    unittest.main()
