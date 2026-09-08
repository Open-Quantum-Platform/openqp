"""The periodic/embedding [qmmm] controls belong to the OpenQpQMMM driver
(runtype=md / namd); the legacy single-point QM/MM path must reject them
instead of silently running a NoCutoff point-charge job."""
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_checker():
    """Load the checker module the way tests/test_d4_input_checker.py does:
    a stub oqp.utils.mpi_utils so the module imports without the runtime."""
    import sys, types
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    if "oqp.utils.mpi_utils" not in sys.modules:
        mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

        class MPIManager:
            use_mpi = False
            size = 1
        mpi_utils.MPIManager = MPIManager
        sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    name = "input_checker_qmmm_sp_under_test"
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


class TestSinglePointQMMMOptions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chk = _load_checker()

    def _errors(self, qmmm, runtype="energy", md=None):
        cfg = {"input": {"runtype": runtype, "qmmm_flag": True, "method": "hf", "basis": "6-31g",
                         "system": "ala.pdb 9 10 17 18 19", "charge": 0},
               "scf": {"type": "rhf", "multiplicity": 1},
               "qmmm": dict(qmmm)}
        if md is not None:
            cfg["md"] = dict(md)
        report = self.chk.CheckReport()
        self.chk._check_qmmm_driver_options(cfg, report)
        return [d for d in report.diagnostics
                if d.severity == "ERROR" and d.path in ("qmmm.cutoff", "qmmm.embedding")]

    def test_single_point_rejects_driver_only_controls(self):
        for qmmm in ({"cutoff": "PME"}, {"mm_charge_width": "0.7"}, {"ewald_tol": "1e-6"},
                     {"lj_switch": True}, {"h_lj": "true"}):
            self.assertTrue(self._errors(qmmm), qmmm)

    def test_single_point_accepts_plain_nocutoff(self):
        self.assertFalse(self._errors({"cutoff": "NoCutoff", "mm_charge_width": "0", "lj_switch": False}))

    def test_soc_namd_rejects_periodic_cutoff(self):
        for cutoff in ("PME", "Ewald", "CutoffPeriodic"):
            errs = self._errors({"cutoff": cutoff}, runtype="namd", md={"soc": True})
            self.assertEqual(len(errs), 1, cutoff)
            self.assertIn("SOC-NAMD", errs[0].message)
        # same-spin FSSH keeps the periodic box; SOC-NAMD keeps a cluster
        self.assertFalse(self._errors({"cutoff": "PME"}, runtype="namd", md={"soc": False}))
        self.assertFalse(self._errors({"cutoff": "NoCutoff"}, runtype="namd", md={"soc": True}))

    def test_periodic_dynamics_rejects_split_embedding(self):
        for rt in ("md", "namd"):
            for cutoff in ("PME", "Ewald", "CutoffPeriodic"):
                errs = self._errors({"cutoff": cutoff, "embedding": "split"}, runtype=rt)
                self.assertEqual([e.path for e in errs], ["qmmm.embedding"], (rt, cutoff))
            self.assertFalse(self._errors({"cutoff": "NoCutoff", "embedding": "split"}, runtype=rt))
            self.assertFalse(self._errors({"cutoff": "PME", "embedding": "electrostatic"}, runtype=rt))
            self.assertFalse(self._errors({"cutoff": "PME", "embedding": "mechanical"}, runtype=rt))

    def test_periodic_tdhf_dynamics_warns_about_zvector_convergence(self):
        def warns(cfg_extra, runtype="namd"):
            cfg = {"input": {"runtype": runtype, "qmmm_flag": True, "method": "tdhf", "basis": "6-31g",
                             "system": "ala.pdb 9 10 17 18 19", "charge": 0},
                   "scf": {"type": "rohf", "multiplicity": 3}, "qmmm": {"cutoff": "PME"}}
            cfg.update(cfg_extra)
            report = self.chk.CheckReport()
            self.chk._check_qmmm_driver_options(cfg, report)
            return [d for d in report.diagnostics if d.path == "tdhf.zvconv"]
        self.assertEqual([d.severity for d in warns({})], ["WARNING"])              # default 1e-6
        self.assertEqual([d.severity for d in warns({"tdhf": {"zvconv": "1e-6"}})], ["WARNING"])
        self.assertFalse(warns({"tdhf": {"zvconv": "1e-8"}}))
        self.assertFalse(warns({"qmmm": {"cutoff": "NoCutoff"}}))
        self.assertFalse(warns({"input": {"runtype": "namd", "qmmm_flag": True, "method": "hf"}}))

    def test_md_and_namd_keep_the_controls(self):
        for rt in ("md", "namd"):
            self.assertFalse(self._errors({"cutoff": "PME", "mm_charge_width": "0.7"}, runtype=rt))


if __name__ == "__main__":
    unittest.main()
