"""The input checker refuses TRAH for MOM/rstctmo calculations.

TRAH minimises over orbital rotations without re-diagonalising, so neither
occupation rule can act; from a non-aufbau stationary point it descends to a
lower state and reports convergence.  scf_driver aborts on the combination,
and the checker reports it before any work starts.
"""

import importlib.util
import sys
import types
import unittest
from pathlib import Path

from oqp_module_isolation import isolate_modules, stub_module


ROOT = Path(__file__).resolve().parents[1]
CHECKER_MODULE = "input_checker_mom_trah_under_test"


def load_input_checker():
    stub_module("oqp")
    stub_module("oqp.utils")
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        use_mpi = False
        size = 1

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    name = CHECKER_MODULE
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/utils/input_checker.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class MomTrahInputCheckerTests(unittest.TestCase):
    def setUp(self):
        # The stubs live only for this test; the conftest guard fails any
        # test that leaves them in sys.modules for the next one.
        self.addCleanup(isolate_modules(CHECKER_MODULE))
        self.checker = load_input_checker()

    def check(self, **scf):
        report = self.checker.CheckReport()
        config = {"input": {"method": "hf"},
                  "scf": {"type": "uhf", "multiplicity": 1, **scf}}
        self.checker._check_scf(config, report)
        return report

    def keys(self, report):
        return {issue.path for issue in report.errors}

    def test_trah_with_occupation_rule_is_an_error(self):
        for rule in ("mom", "rstctmo"):
            for value in (True, "true"):
                report = self.check(**{rule: value, "converger_type": "trah"})
                self.assertIn("scf.converger_type", self.keys(report), (rule, value))
                self.assertIn("TRAH cannot keep", report.to_text())
                report = self.check(**{rule: value, "init_converger": "trah"})
                self.assertIn("scf.init_converger", self.keys(report), (rule, value))

    def test_occupation_rule_with_diis_or_soscf_is_accepted(self):
        for rule in ("mom", "rstctmo"):
            for conv in ("diis", "soscf"):
                report = self.check(**{rule: True, "converger_type": conv})
                self.assertNotIn("scf.converger_type", self.keys(report), (rule, conv))

    def test_trah_without_occupation_rule_is_accepted(self):
        report = self.check(converger_type="trah", mom=False, rstctmo=False)
        self.assertNotIn("scf.converger_type", self.keys(report))


if __name__ == "__main__":
    unittest.main()
