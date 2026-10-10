import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_minimal_oqp_stubs():
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        pass

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils


class TestAdvancedGuessConfig(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()

    def test_input_checker_accepts_native_guess_modes(self):
        input_checker = load_module(
            "input_checker_native_guess_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )

        for guess_type in ("hcore", "huckel", "modhuckel", "minao", "sap"):
            with self.subTest(guess_type=guess_type):
                config = {"guess": {"type": guess_type}}
                report = input_checker.CheckReport()
                input_checker._check_guess(config, report)
                self.assertTrue(report.ok, report.to_text())

    def test_input_checker_rejects_removed_sad_guess(self):
        input_checker = load_module(
            "input_checker_sad_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {"guess": {"type": "sad"}}

        report = input_checker.CheckReport()
        input_checker._check_guess(config, report)

        self.assertFalse(report.ok)
        self.assertIn("Unknown guess type", report.to_text())


if __name__ == "__main__":
    unittest.main()
