"""Input validation for the deterministic E-DIIS/A-DIIS solver limit."""

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
    name = "input_checker_simplex_limit_under_test"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp/oqp/utils/input_checker.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class SCFSimplexLimitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.checker = load_input_checker()

    def check(self, diis_type, maxdiis):
        report = self.checker.CheckReport()
        self.checker._check_scf(
            {"scf": {"type": "rhf", "multiplicity": 1,
                     "diis_type": diis_type, "maxdiis": maxdiis}},
            report,
        )
        return report

    def test_ediis_family_rejects_histories_above_exact_limit(self):
        for diis_type in ("ediis", "adiis", "vdiis"):
            with self.subTest(diis_type=diis_type):
                report = self.check(diis_type, 14)
                self.assertFalse(report.ok)
                self.assertIn("scf.maxdiis", report.to_text())
                self.assertIn("at most 13", report.to_text())

    def test_limit_and_large_cdiis_history_remain_valid(self):
        for diis_type, maxdiis in (("ediis", 13), ("adiis", 13),
                                   ("vdiis", 13), ("cdiis", 20)):
            with self.subTest(diis_type=diis_type, maxdiis=maxdiis):
                self.assertTrue(self.check(diis_type, maxdiis).ok)


if __name__ == "__main__":
    unittest.main()
