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


def tdhf_config(td_type, td_mult, runtype="energy"):
    return {
        "input": {"method": "tdhf", "runtype": runtype, "functional": "bhhlyp"},
        "scf": {"type": "rhf", "multiplicity": 1},
        "tdhf": {"type": td_type, "multiplicity": td_mult, "nstate": 3},
    }


class TestConventionalTDHFMultiplicity(unittest.TestCase):
    """The closed-shell RPA/TDA response has no triplet path.

    tdhf_energy always includes the Coulomb term and the alpha-plus-beta XC
    kernel, so [tdhf] multiplicity=3 used to return the singlet roots under a
    triplet label.  The input checker must reject it instead.
    """

    def setUp(self):
        install_minimal_oqp_stubs()
        self.input_checker = load_input_checker("input_checker_tdhf_mult_under_test")

    def test_rpa_and_tda_reject_triplet_multiplicity(self):
        for td_type in ("rpa", "tda"):
            for runtype in ("energy", "grad"):
                with self.subTest(td_type=td_type, runtype=runtype):
                    report = self.input_checker.CheckReport()
                    self.input_checker._check_tdhf(
                        tdhf_config(td_type, 3, runtype), report
                    )
                    self.assertFalse(report.ok)
                    self.assertIn("tdhf.multiplicity", report.to_text())
                    self.assertIn("singlet excited states only", report.to_text())

    def test_rpa_and_tda_accept_singlet_multiplicity(self):
        for td_type in ("rpa", "tda"):
            with self.subTest(td_type=td_type):
                report = self.input_checker.CheckReport()
                self.input_checker._check_tdhf(tdhf_config(td_type, 1), report)
                self.assertTrue(report.ok, report.to_text())

    def test_mrsf_triplet_is_still_accepted(self):
        config = tdhf_config("mrsf", 3)
        config["scf"] = {"type": "rohf", "multiplicity": 3}
        report = self.input_checker.CheckReport()
        self.input_checker._check_tdhf(config, report)
        self.assertTrue(report.ok, report.to_text())


if __name__ == "__main__":
    unittest.main()
