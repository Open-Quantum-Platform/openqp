import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_OPT = ROOT / "examples" / "OPT"


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_minimal_oqp_stubs():
    oqp = sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    oqp_utils = sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        pass

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    return oqp, oqp_utils


def install_runfunc_stubs():
    install_minimal_oqp_stubs()
    oqp_library = sys.modules.setdefault("oqp.library", types.ModuleType("oqp.library"))
    oqp_library.__path__ = []

    single_point = types.ModuleType("oqp.library.single_point")
    for class_name in (
        "SinglePoint", "Gradient", "Hessian", "LastStep",
        "BasisOverlap", "NACME", "NAC",
    ):
        setattr(single_point, class_name, type(class_name, (), {}))
    sys.modules["oqp.library.single_point"] = single_point

    libscipy = types.ModuleType("oqp.library.libscipy")
    libscipy.StateSpecificOpt = type("StateSpecificOpt", (), {})
    libscipy.MECIOpt = type("MECIOpt", (), {})
    libscipy.MECPOpt = type("MECPOpt", (), {})
    libscipy.MEP = type("MEP", (), {})
    libscipy.QMMMOpt = type("QMMMOpt", (), {})
    sys.modules["oqp.library.libscipy"] = libscipy

    libgeometric = types.ModuleType("oqp.library.libgeometric")

    class GeometricOpt:
        def __init__(self, mol):
            self.mol = mol

    class GeometricMECIOpt:
        def __init__(self, mol):
            self.mol = mol

    class GeometricMECPOpt:
        def __init__(self, mol):
            self.mol = mol

    class GeometricTSOpt:
        def __init__(self, mol):
            self.mol = mol

    class GeometricIRCOpt:
        def __init__(self, mol):
            self.mol = mol

    class GeometricNEBOpt:
        def __init__(self, mol):
            self.mol = mol

    libgeometric.GeometricOpt = GeometricOpt
    libgeometric.GeometricMECIOpt = GeometricMECIOpt
    libgeometric.GeometricMECPOpt = GeometricMECPOpt
    libgeometric.GeometricTSOpt = GeometricTSOpt
    libgeometric.GeometricIRCOpt = GeometricIRCOpt
    setattr(libgeometric, "GeometricNEBOpt", GeometricNEBOpt)
    sys.modules["oqp.library.libgeometric"] = libgeometric

    # runfunc also imports the oqp optimizer backend; stub it too.
    liboqp = types.ModuleType("oqp.library.liboqp")
    for _cls in ("OQPOpt", "OQPTSOpt", "OQPMECIOpt", "OQPMECPOpt",
                 "OQPTCIOpt", "OQPNEBOpt", "OQPIRCOpt", "OQPMEPOpt"):
        setattr(liboqp, _cls, type(_cls, (), {}))
    sys.modules["oqp.library.liboqp"] = liboqp

    return GeometricOpt, GeometricMECIOpt, GeometricMECPOpt, GeometricTSOpt, GeometricIRCOpt


def input_section(text, section_name):
    marker = f"[{section_name}]"
    lower = text.lower()
    if marker not in lower:
        return ""
    return lower.split(marker, 1)[1].split("[", 1)[0]


def section_value(section_text, key):
    for line in section_text.splitlines():
        stripped = line.strip()
        if stripped.startswith(f"{key}="):
            return stripped.split("=", 1)[1].strip()
    return None


class TestOptimizationExampleCaps(unittest.TestCase):
    def test_geometry_path_examples_have_bounded_iterations(self):
        geometry_runtypes = {"optimize", "meci", "mecp", "ts", "irc", "neb", "mep"}
        failures = []
        checked = 0
        for path in sorted(EXAMPLES_OPT.glob("*.inp")):
            text = path.read_text()
            input_sec = input_section(text, "input")
            runtype = section_value(input_sec, "runtype")
            if runtype not in geometry_runtypes:
                continue
            checked += 1
            optimize = input_section(text, "optimize")
            scf = input_section(text, "scf")
            tdhf = input_section(text, "tdhf")
            opt_maxit = section_value(optimize, "maxit")
            if opt_maxit is None:
                failures.append(f"{path.name}: missing [optimize] maxit")
            elif int(float(opt_maxit)) > 10:
                failures.append(f"{path.name}: [optimize] maxit={opt_maxit} > 10")
            scf_maxit = section_value(scf, "maxit")
            if scf_maxit is not None and int(float(scf_maxit)) > 30:
                failures.append(f"{path.name}: [scf] maxit={scf_maxit} > 30")
            tdhf_maxit = section_value(tdhf, "maxit")
            if "type=mrsf" in tdhf and tdhf_maxit is None:
                failures.append(f"{path.name}: MRSF [tdhf] missing maxit")
            elif tdhf_maxit is not None and int(float(tdhf_maxit)) > 30:
                failures.append(f"{path.name}: [tdhf] maxit={tdhf_maxit} > 30")
        self.assertGreater(checked, 0)
        self.assertEqual(failures, [])


class TestGeometricOptimizerConfig(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()

    def test_geometric_optimizer_examples_are_available(self):
        expected = {
            "H2O_RHF-DFT_OPTIMIZE_GEOMETRIC.inp": "runtype=optimize",
            "H2O_RHF-DFT_OPTIMIZE_GEOMETRIC.json": None,
            "C2H4_BHHLYP-MRSFTDDFT_MECI_GEOMETRIC.inp": "runtype=meci",
            "C2H4_BHHLYP-MRSFTDDFT_MECI_GEOMETRIC.json": None,
            "C2H4_BHHLYP-MRSFTDDFT_MECP_GEOMETRIC.inp": "runtype=mecp",
            "C2H4_BHHLYP-MRSFTDDFT_MECP_GEOMETRIC.json": None,
            "HCN_RHF-DFT_TS_GEOMETRIC.inp": "runtype=ts",
            "HCN_RHF-DFT_TS_GEOMETRIC.json": None,
            "HCN_BHHLYP-MRSFTDDFT_TS_GEOMETRIC.inp": "runtype=ts",
            "HCN_BHHLYP-MRSFTDDFT_TS_GEOMETRIC.json": None,
            "HCN_RHF-DFT_IRC_GEOMETRIC.inp": "runtype=irc",
            "HCN_RHF-DFT_IRC_GEOMETRIC.json": None,
            "HCN_RHF-DFT_CONSTRAINED_GEOMETRIC.inp": "runtype=optimize",
            "HCN_RHF-DFT_CONSTRAINED_GEOMETRIC.constraints": "$freeze",
            "HCN_RHF-DFT_CONSTRAINED_GEOMETRIC.json": None,
        }

        missing = sorted(name for name in expected if not (EXAMPLES_OPT / name).is_file())
        self.assertEqual(missing, [])
        for name, runtype in expected.items():
            if not name.endswith(".inp"):
                continue
            text = (EXAMPLES_OPT / name).read_text()
            self.assertIn(runtype, text)
            self.assertIn("lib=geometric", text)
            self.assertIn("[geometric]", text)

    def test_geometric_config_supports_constraints_file_options(self):
        text = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()

        self.assertIn("'constraints_file': {'type': str, 'default': ''}", text)
        self.assertIn("'enforce': {'type': float, 'default': '0.0'}", text)
        self.assertIn("'conmethod': {'type': int, 'default': '0'}", text)

    def test_geometric_runner_passes_constraints_to_optimizer(self):
        install_runfunc_stubs()
        file_utils = types.ModuleType("oqp.utils.file_utils")
        file_utils.dump_log = lambda *args, **kwargs: None
        sys.modules["oqp.utils.file_utils"] = file_utils
        libgeometric = load_module(
            "libgeometric_constraints_under_test",
            "pyoqp/oqp/library/libgeometric.py",
        )

        runner = libgeometric._GeometricRunner.__new__(libgeometric._GeometricRunner)
        runner.geometric_config = {
            "constraints_file": "constraints.txt",
            "enforce": 1e-4,
            "conmethod": 1,
        }

        keywords = runner._optimizer_keywords()

        self.assertEqual(keywords["constraints"], "constraints.txt")
        self.assertEqual(keywords["enforce"], 1e-4)
        self.assertEqual(keywords["conmethod"], 1)

    def test_input_checker_accepts_geometric_for_state_specific_optimize(self):
        input_checker = load_module(
            "input_checker_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "optimize", "method": "hf"},
            "optimize": {"lib": "geometric", "istate": 0},
        }

        report = input_checker.CheckReport()
        input_checker._check_optimize(config, report)

        self.assertTrue(report.ok, report.to_text())

    def test_input_checker_accepts_geometric_for_meci(self):
        input_checker = load_module(
            "input_checker_geometric_meci_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "meci", "method": "tdhf"},
            "optimize": {"lib": "geometric", "istate": 1, "jstate": 2},
        }

        report = input_checker.CheckReport()
        input_checker._check_optimize(config, report)

        self.assertTrue(report.ok, report.to_text())

    def test_get_optimizer_dispatches_geometric_optimize(self):
        GeometricOpt, _, _, _, _ = install_runfunc_stubs()
        runfunc = load_module(
            "runfunc_under_test",
            "pyoqp/oqp/library/runfunc.py",
        )

        mol = types.SimpleNamespace(
            config={
                "input": {"runtype": "optimize"},
                "optimize": {"lib": "geometric"},
            }
        )

        optimizer = runfunc.get_optimizer(mol)

        self.assertIsInstance(optimizer, GeometricOpt)
        self.assertIs(optimizer.mol, mol)

    def test_get_optimizer_dispatches_geometric_meci(self):
        _, GeometricMECIOpt, _, _, _ = install_runfunc_stubs()
        runfunc = load_module(
            "runfunc_geometric_meci_under_test",
            "pyoqp/oqp/library/runfunc.py",
        )

        mol = types.SimpleNamespace(
            config={
                "input": {"runtype": "meci"},
                "optimize": {"lib": "geometric"},
            }
        )

        optimizer = runfunc.get_optimizer(mol)

        self.assertIsInstance(optimizer, GeometricMECIOpt)
        self.assertIs(optimizer.mol, mol)

    def test_input_checker_accepts_geometric_for_mecp(self):
        input_checker = load_module(
            "input_checker_geometric_mecp_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "mecp", "method": "tdhf"},
            "optimize": {
                "lib": "geometric",
                "istate": 1,
                "jstate": 1,
                "imult": 1,
                "jmult": 3,
            },
        }

        report = input_checker.CheckReport()
        input_checker._check_optimize(config, report)

        self.assertTrue(report.ok, report.to_text())

    def test_input_checker_accepts_geometric_for_ts(self):
        input_checker = load_module(
            "input_checker_geometric_ts_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "ts", "method": "hf"},
            "optimize": {"lib": "geometric", "istate": 0},
        }

        report = input_checker.CheckReport()
        input_checker._check_optimize(config, report)

        self.assertTrue(report.ok, report.to_text())

    def test_get_optimizer_dispatches_geometric_mecp(self):
        _, _, GeometricMECPOpt, _, _ = install_runfunc_stubs()
        runfunc = load_module(
            "runfunc_geometric_mecp_under_test",
            "pyoqp/oqp/library/runfunc.py",
        )

        mol = types.SimpleNamespace(
            config={
                "input": {"runtype": "mecp"},
                "optimize": {"lib": "geometric"},
            }
        )

        optimizer = runfunc.get_optimizer(mol)

        self.assertIsInstance(optimizer, GeometricMECPOpt)
        self.assertIs(optimizer.mol, mol)

    def test_get_optimizer_dispatches_geometric_ts(self):
        _, _, _, GeometricTSOpt, _ = install_runfunc_stubs()
        runfunc = load_module(
            "runfunc_geometric_ts_under_test",
            "pyoqp/oqp/library/runfunc.py",
        )

        mol = types.SimpleNamespace(
            config={
                "input": {"runtype": "ts"},
                "optimize": {"lib": "geometric"},
            }
        )

        optimizer = runfunc.get_optimizer(mol)

        self.assertIsInstance(optimizer, GeometricTSOpt)
        self.assertIs(optimizer.mol, mol)

    def test_full_input_checker_accepts_default_lib_for_irc(self):
        input_checker = load_module(
            "input_checker_full_irc_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {
                "runtype": "irc",
                "method": "hf",
                "basis": "3-21g",
                "system": "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.8",
            },
            "guess": {},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {},
            "properties": {},
            "optimize": {"istate": 0},
        }

        report = input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())

    def test_optimize_lib_default_is_oqp(self):
        text = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()

        self.assertIn("'lib': {'type': str, 'default': 'oqp'}", text)

    def test_geometric_ts_uses_initial_hessian_by_default(self):
        text = (ROOT / "pyoqp/oqp/library/libgeometric.py").read_text()
        ts_block = text.split("class GeometricTSOpt", 1)[1].split("class GeometricIRCOpt", 1)[0]

        self.assertIn('if self.hessian == "never":', ts_block)
        self.assertIn('self.hessian = "first"', ts_block)

    def test_input_checker_accepts_geometric_for_irc(self):
        input_checker = load_module(
            "input_checker_geometric_irc_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "irc", "method": "hf"},
            "optimize": {"lib": "geometric", "istate": 0},
        }

        report = input_checker.CheckReport()
        input_checker._check_optimize(config, report)

        self.assertTrue(report.ok, report.to_text())

    def test_get_optimizer_dispatches_geometric_irc(self):
        _, _, _, _, GeometricIRCOpt = install_runfunc_stubs()
        runfunc = load_module(
            "runfunc_geometric_irc_under_test",
            "pyoqp/oqp/library/runfunc.py",
        )

        mol = types.SimpleNamespace(
            config={
                "input": {"runtype": "irc"},
                "optimize": {"lib": "geometric"},
            }
        )

        optimizer = runfunc.get_optimizer(mol)

        self.assertIsInstance(optimizer, GeometricIRCOpt)
        self.assertIs(optimizer.mol, mol)


if __name__ == "__main__":
    unittest.main()
