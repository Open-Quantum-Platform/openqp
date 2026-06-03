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
        use_mpi = False
        size = 1

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils

    bse = types.ModuleType("basis_set_exchange")

    def get_basis(name, elements=None):
        max_l = {"cc-pV5Z": 4, "cc-pVTZ": 3}[name]
        return {
            "elements": {
                "8": {"electron_shells": [{"angular_momentum": [max_l]}]},
                "1": {"electron_shells": [{"angular_momentum": [0]}]},
            }
        }

    bse.get_basis = get_basis
    sys.modules["basis_set_exchange"] = bse


class TestAnalyticalHessianBasisLGate(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()
        self.input_checker = load_module(
            "input_checker_hessian_l_gate_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )

    def test_analytical_hessian_rejects_g_or_higher_basis_functions(self):
        config = {
            "input": {
                "method": "hf",
                "runtype": "hess",
                "basis": "cc-pV5Z",
                "system": "\nO 0.0 0.0 0.0\nH 0.0 0.0 0.9",
            },
            "scf": {"type": "rhf"},
            "hess": {"type": "analytical", "state": 0},
        }
        report = self.input_checker.CheckReport()

        self.input_checker._check_hess(config, report)

        text = report.to_text()
        self.assertFalse(report.ok, text)
        self.assertIn("max L=4", text)
        self.assertIn("max L <= 3", text)

    def test_analytical_hessian_accepts_f_or_lower_basis_functions(self):
        config = {
            "input": {
                "method": "hf",
                "runtype": "hess",
                "basis": "cc-pVTZ",
                "system": "\nO 0.0 0.0 0.0\nH 0.0 0.0 0.9",
            },
            "scf": {"type": "rhf"},
            "hess": {"type": "analytical", "state": 0},
        }
        report = self.input_checker.CheckReport()

        self.input_checker._check_hess(config, report)

        self.assertTrue(report.ok, report.to_text())


if __name__ == "__main__":
    unittest.main()
