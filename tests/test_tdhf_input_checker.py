import importlib.util
import sys
import types
import unittest
from pathlib import Path

from oqp_module_isolation import isolate_modules, stub_module


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
    stub_module("oqp")
    stub_module("oqp.utils")
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
    """The closed-shell RPA/TDA response has singlet and triplet branches
    (no Coulomb term, alpha-minus-beta kernel) in tdhf_energy, tdhf_z_vector
    and tdhf_gradient: triplet energies, gradients and optimizations are
    accepted on RHF; Hessians/NACs, other multiplicities and open-shell
    references are rejected.
    """

    def setUp(self):
        self.addCleanup(isolate_modules())
        install_minimal_oqp_stubs()
        self.input_checker = load_input_checker("input_checker_tdhf_mult_under_test")

    def test_rpa_and_tda_accept_triplet_energy_on_rhf(self):
        for td_type in ("rpa", "tda"):
            with self.subTest(td_type=td_type):
                report = self.input_checker.CheckReport()
                self.input_checker._check_tdhf(tdhf_config(td_type, 3, "energy"), report)
                self.assertTrue(report.ok, report.to_text())

    def test_rpa_and_tda_accept_triplet_gradient_and_optimization(self):
        for td_type in ("rpa", "tda"):
            for runtype in ("grad", "optimize"):
                with self.subTest(td_type=td_type, runtype=runtype):
                    report = self.input_checker.CheckReport()
                    self.input_checker._check_tdhf(tdhf_config(td_type, 3, runtype), report)
                    self.assertTrue(report.ok, report.to_text())

    def test_triplet_gradient_rejects_meta_gga(self):
        for runtype in ("grad", "optimize"):
            with self.subTest(runtype=runtype):
                config = tdhf_config("rpa", 3, runtype)
                config["input"]["functional"] = "m06-2x"
                report = self.input_checker.CheckReport()
                self.input_checker._check_tdhf(config, report)
                self.assertFalse(report.ok)
                self.assertIn("meta-GGA", report.to_text())
        config = tdhf_config("rpa", 3, "energy")
        config["input"]["functional"] = "m06-2x"
        report = self.input_checker.CheckReport()
        self.input_checker._check_tdhf(config, report)
        self.assertTrue(report.ok, report.to_text())

    def test_rpa_and_tda_reject_triplet_hessian_and_nac(self):
        for td_type in ("rpa", "tda"):
            for runtype in ("hess", "nac"):
                with self.subTest(td_type=td_type, runtype=runtype):
                    report = self.input_checker.CheckReport()
                    self.input_checker._check_tdhf(tdhf_config(td_type, 3, runtype), report)
                    self.assertFalse(report.ok)
                    self.assertIn("Hessians, NACs, crossings", report.to_text())

    def test_rpa_and_tda_reject_triplet_on_open_shell_reference(self):
        config = tdhf_config("rpa", 3, "energy")
        config["scf"] = {"type": "uhf", "multiplicity": 1}
        report = self.input_checker.CheckReport()
        self.input_checker._check_tdhf(config, report)
        self.assertFalse(report.ok)
        self.assertIn("open-shell", report.to_text())

    def test_rpa_and_tda_reject_other_multiplicities(self):
        for td_type in ("rpa", "tda"):
            with self.subTest(td_type=td_type):
                report = self.input_checker.CheckReport()
                self.input_checker._check_tdhf(tdhf_config(td_type, 5, "energy"), report)
                self.assertFalse(report.ok)
                self.assertIn("tdhf.multiplicity", report.to_text())

    def test_rpa_and_tda_accept_singlet_multiplicity(self):
        for td_type in ("rpa", "tda"):
            with self.subTest(td_type=td_type):
                report = self.input_checker.CheckReport()
                self.input_checker._check_tdhf(tdhf_config(td_type, 1), report)
                self.assertTrue(report.ok, report.to_text())

    def test_rpa_and_tda_mecp_reject_a_triplet_surface(self):
        # MECPOpt switches tdhf.multiplicity to imult/jmult, so a singlet
        # tdhf.multiplicity alone does not keep the triplet surface away.
        for td_type in ("rpa", "tda"):
            for imult, jmult in ((1, 3), (3, 1), (3, 3)):
                with self.subTest(td_type=td_type, imult=imult, jmult=jmult):
                    config = tdhf_config(td_type, 1, "mecp")
                    config["optimize"] = {"imult": imult, "jmult": jmult}
                    report = self.input_checker.CheckReport()
                    self.input_checker._check_tdhf(config, report)
                    self.assertFalse(report.ok)
                    self.assertIn("MECP cannot evaluate", report.to_text())

    def test_rpa_mecp_malformed_multiplicity_does_not_raise(self):
        config = tdhf_config("rpa", 1, "mecp")
        config["optimize"] = {"imult": "abc", "jmult": 1}
        report = self.input_checker.CheckReport()
        self.input_checker._check_tdhf(config, report)
        self.assertNotIn("MECP cannot evaluate", report.to_text())

    def test_mrsf_mecp_is_still_accepted(self):
        config = tdhf_config("mrsf", 1, "mecp")
        config["scf"] = {"type": "rohf", "multiplicity": 3}
        config["optimize"] = {"imult": 1, "jmult": 3}
        report = self.input_checker.CheckReport()
        self.input_checker._check_tdhf(config, report)
        self.assertTrue(report.ok, report.to_text())

    def test_mrsf_triplet_is_still_accepted(self):
        config = tdhf_config("mrsf", 3)
        config["scf"] = {"type": "rohf", "multiplicity": 3}
        report = self.input_checker.CheckReport()
        self.input_checker._check_tdhf(config, report)
        self.assertTrue(report.ok, report.to_text())


if __name__ == "__main__":
    unittest.main()
