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

    def check(self, diis_type, maxdiis, vshift=0.0, converger="diis",
              init_scf="no", init_converger="diis", escalation="",
              alternative_scf="trah"):
        report = self.checker.CheckReport()
        self.checker._check_scf(
            {"scf": {"type": "rhf", "multiplicity": 1,
                     "diis_type": diis_type, "maxdiis": maxdiis,
                     "vshift": vshift, "converger_type": converger,
                     "init_scf": init_scf,
                     "init_converger": init_converger,
                     "escalation": escalation,
                     "alternative_scf": alternative_scf}},
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

    def test_level_shifted_cdiis_uses_the_simplex_limit(self):
        report = self.check("cdiis", 14, vshift=0.1)
        self.assertFalse(report.ok)
        self.assertIn("level-shifted DIIS", report.to_text())
        self.assertTrue(self.check("cdiis", 13, vshift=0.1).ok)

    def test_level_shift_overrides_every_diis_subtype(self):
        report = self.check("none", 14, vshift=0.1)
        self.assertFalse(report.ok)
        self.assertIn("scf.maxdiis", report.to_text())

    def test_ml_converger_is_limited_because_it_can_select_ediis(self):
        report = self.check("cdiis", 14, converger="ml")
        self.assertFalse(report.ok)
        self.assertIn("ML-selected DIIS", report.to_text())
        self.assertTrue(self.check("cdiis", 13, converger="ml").ok)

    def test_dormant_diis_settings_do_not_limit_soscf_or_trah(self):
        for converger in ("soscf", "trah"):
            with self.subTest(converger=converger):
                self.assertTrue(
                    self.check("ediis", 14, converger=converger).ok
                )

    def test_active_initial_diis_stage_is_limited(self):
        report = self.check(
            "adiis", 14, converger="soscf", init_scf="rhf",
            init_converger="diis",
        )
        self.assertFalse(report.ok)
        self.assertTrue(
            self.check(
                "adiis", 14, converger="soscf", init_scf="rhf",
                init_converger="soscf",
            ).ok
        )

    def test_explicit_diis_escalation_is_limited(self):
        report = self.check(
            "ediis", 14, converger="soscf", escalation="diis,trah"
        )
        self.assertFalse(report.ok)
        self.assertIn("scf.maxdiis", report.to_text())

    def test_diis_alternative_is_limited_only_without_override(self):
        report = self.check(
            "adiis", 14, converger="soscf", alternative_scf="diis"
        )
        self.assertFalse(report.ok)
        self.assertTrue(
            self.check(
                "adiis", 14, converger="soscf", alternative_scf="diis",
                escalation="trah",
            ).ok
        )


if __name__ == "__main__":
    unittest.main()
