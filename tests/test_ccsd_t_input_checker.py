"""Input-validation coverage for the closed-shell CCSD(T) method.

The Fortran kernels are validated numerically elsewhere; these tests pin the
guard rails that stop an unsupported job reaching them -- a DFT reference, an
open-shell reference, or a derivative runtype the module does not implement.
"""
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


class TestCCSDTInputChecker(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()
        self.checker = load_input_checker("input_checker_cc_under_test")

    def _report(self):
        return self.checker.CheckReport()

    def test_cc_methods_are_recognised(self):
        self.assertIn("ccsd", self.checker.METHODS)
        self.assertIn("ccsd(t)", self.checker.METHODS)

    def test_cc_rejects_dft_functional(self):
        config = {"input": {"method": "ccsd(t)", "functional": "pbe"},
                  "scf": {"type": "rhf"}}
        report = self._report()
        self.checker._check_cc(config, report)

        self.assertFalse(report.ok)
        self.assertIn("input.functional", report.to_text())

    def test_cc_accepts_uhf_with_a_note(self):
        """UHF goes through the spin-orbital solver: allowed, but costly enough
        that the user is told so rather than left to wonder."""
        config = {"input": {"method": "ccsd", "functional": ""},
                  "scf": {"type": "uhf"}}
        report = self._report()
        self.checker._check_cc(config, report)

        self.assertTrue(report.ok)
        self.assertIn("spin-orbital", report.to_text())

    def test_cc_accepts_rohf(self):
        """ROHF is supported once f_ov is carried through the equations and the
        correlation energy; it goes through the same spin-orbital solver."""
        config = {"input": {"method": "ccsd(t)", "functional": ""},
                  "scf": {"type": "rohf"}}
        report = self._report()
        self.checker._check_cc(config, report)

        self.assertTrue(report.ok)
        self.assertIn("spin-orbital", report.to_text())

    def test_cc_rejects_unknown_reference(self):
        config = {"input": {"method": "ccsd", "functional": ""},
                  "scf": {"type": "gvb"}}
        report = self._report()
        self.checker._check_cc(config, report)

        self.assertFalse(report.ok)
        self.assertIn("scf.type", report.to_text())

    def test_cc_rejects_negative_frozen_core(self):
        config = {"input": {"method": "ccsd(t)", "functional": ""},
                  "scf": {"type": "rhf"}, "cc": {"nfzc": -1}}
        report = self._report()
        self.checker._check_cc(config, report)

        self.assertFalse(report.ok)
        self.assertIn("cc.nfzc", report.to_text())

    def test_cc_rejects_non_positive_maxit(self):
        config = {"input": {"method": "ccsd(t)", "functional": ""},
                  "scf": {"type": "rhf"}, "cc": {"maxit": 0}}
        report = self._report()
        self.checker._check_cc(config, report)

        self.assertFalse(report.ok)
        self.assertIn("cc.maxit", report.to_text())

    def test_cc_accepts_closed_shell_hf_reference(self):
        config = {"input": {"method": "ccsd(t)", "functional": ""},
                  "scf": {"type": "rhf"}, "cc": {"nfzc": 1, "maxit": 50}}
        report = self._report()
        self.checker._check_cc(config, report)

        self.assertTrue(report.ok, msg=report.to_text())

    def test_cc_rejects_derivative_runtype(self):
        config = {"input": {"method": "ccsd(t)", "runtype": "grad"}}
        report = self._report()
        self.checker._check_runtype(config, report)

        self.assertFalse(report.ok)
        self.assertIn("input.runtype", report.to_text())

    def test_cc_allows_energy_runtype(self):
        config = {"input": {"method": "ccsd(t)", "runtype": "energy"}}
        report = self._report()
        self.checker._check_runtype(config, report)

        self.assertTrue(report.ok, msg=report.to_text())


    def test_cc_rejects_a_non_positive_convergence_threshold(self):
        """conv=0 can never be met: both solvers require the amplitude RMS and
        the energy change to fall strictly below it, so the job would burn
        every iteration and then abort."""
        for bad in ("0", "0.0", "-1e-7", "nonsense"):
            config = {"input": {"method": "ccsd(t)", "functional": ""},
                      "scf": {"type": "rhf"}, "cc": {"conv": bad}}
            report = self._report()
            self.checker._check_cc(config, report)

            self.assertFalse(report.ok, bad)
            self.assertIn("cc.conv", report.to_text(), bad)

    def test_cc_accepts_a_positive_convergence_threshold(self):
        config = {"input": {"method": "ccsd(t)", "functional": ""},
                  "scf": {"type": "rhf"}, "cc": {"conv": "1e-8"}}
        report = self._report()
        self.checker._check_cc(config, report)

        self.assertTrue(report.ok)


if __name__ == "__main__":
    unittest.main()
