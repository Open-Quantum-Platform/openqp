import importlib.util
import sys
import tempfile
import types
import unittest
from contextlib import contextmanager
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


@contextmanager
def load_oqp_tester(runner_class, module_name):
    """Load OQPTester with a controlled Runner and no compiled extension."""
    module_names = (
        "oqp",
        "oqp.utils",
        "oqp.pyoqp",
        "oqp.runtime",
        module_name,
    )
    missing = object()
    saved_modules = {
        name: sys.modules.get(name, missing)
        for name in module_names
    }

    try:
        oqp = types.ModuleType("oqp")
        oqp.__path__ = []
        sys.modules["oqp"] = oqp

        oqp_utils = types.ModuleType("oqp.utils")
        oqp_utils.__path__ = []
        sys.modules["oqp.utils"] = oqp_utils

        pyoqp = types.ModuleType("oqp.pyoqp")
        pyoqp.Runner = runner_class
        sys.modules["oqp.pyoqp"] = pyoqp

        runtime = types.ModuleType("oqp.runtime")
        runtime.resolve_oqp_root = lambda: (str(ROOT), "test")
        sys.modules["oqp.runtime"] = runtime

        yield load_module(module_name, "pyoqp/oqp/utils/oqp_tester.py")
    finally:
        for name, module in saved_modules.items():
            if module is missing:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


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

    def test_constrained_optimization_is_now_a_native_example(self):
        path = EXAMPLES_OPT / "HCN_RHF-DFT_CONSTRAINED_OQP.inp"
        self.assertTrue(path.is_file())
        self.assertTrue(path.with_suffix(".json").is_file())
        text = path.read_text()
        self.assertIn("runtype=optimize", text)
        self.assertIn("lib=oqp", text)
        self.assertIn("freeze=distance(1,2)", text)
        self.assertNotIn("[geometric]", text)
        self.assertEqual(list(EXAMPLES_OPT.glob("*GEOMETRIC.inp")), [])

    def test_primary_geometry_examples_use_the_native_backend(self):
        expected = {
            "H2O_RHF-DFT_OPTIMIZE.inp": "runtype=optimize",
            "H2O_RHF-DFT_OPTIMIZE_OQP.inp": "runtype=optimize",
            "HCN_RHF-DFT_CONSTRAINED_OQP.inp": "runtype=optimize",
            "C2H4_BHHLYP-MRSFTDDFT_MECI.inp": "runtype=meci",
            "C2H4_BHHLYP-MRSFTDDFT_MECP_OQP.inp": "runtype=mecp",
            "C2H4_BHHLYP-MRSFTDDFT_TCI_OQP.inp": "runtype=tci",
            "C2H4_BHHLYP-MRSFTDDFT_MEP.inp": "runtype=mep",
            "HCN_RHF-DFT_TS_OQP.inp": "runtype=ts",
            "HCN_BHHLYP-MRSFTDDFT_TS_OQP.inp": "runtype=ts",
            "HCN_RHF-DFT_IRC_OQP.inp": "runtype=irc",
            "HCN_RHF-DFT_NEB_OQP.inp": "runtype=neb",
        }

        for name, runtype in expected.items():
            path = EXAMPLES_OPT / name
            self.assertTrue(path.is_file(), name)
            self.assertTrue(path.with_suffix(".json").is_file(), path.with_suffix(".json").name)
            text = path.read_text()
            self.assertIn(runtype, text)
            self.assertNotIn("lib=geometric", text)
            self.assertNotIn("[geometric]", text)

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

    def test_baeka_states_are_authoritative_over_legacy_pair(self):
        input_checker = load_module(
            "input_checker_baeka_states_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {"runtype": "meci", "method": "tdhf"},
            "tdhf": {"nstate": 4},
            "optimize": {
                "lib": "oqp", "meci_search": "baeka",
                "states": [2, 3, 4],
                # These legacy pair values are intentionally contradictory;
                # the multistate list is the BaekA source of truth.
                "istate": 5, "jstate": 2,
            },
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

    def test_tdhf_irc_state_is_positive_and_within_nstate(self):
        input_checker = load_module(
            "input_checker_full_irc_state_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        base = {
            "input": {"runtype": "irc", "method": "tdhf", "basis": "3-21g",
                      "system": "\nH 0 0 0\nH 0 0 0.7"},
            "guess": {},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "multiplicity": 3, "nstate": 2},
            "properties": {},
            "optimize": {"lib": "oqp", "istate": 0},
            "hess": {"type": "numerical", "state": 1, "nproc": 1},
        }

        report = input_checker.check_input_values(base, raise_error=False, emit=False)
        self.assertIn("TDHF optimization cannot target state 0", report.to_text())

        base["optimize"]["istate"] = 3
        report = input_checker.check_input_values(base, raise_error=False, emit=False)
        self.assertIn("Requested state index exceeds", report.to_text())

        base["optimize"]["istate"] = 1
        base["hess"]["restart"] = True
        report = input_checker.check_input_values(base, raise_error=False, emit=False)
        self.assertIn("restart artifacts are not yet signed", report.to_text())

        base["hess"]["restart"] = False
        base["input"]["qmmm_flag"] = True
        report = input_checker.check_input_values(base, raise_error=False, emit=False)
        self.assertIn("not connected to the active QM/MM force backend", report.to_text())

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


class TestOQPTesterCollection(unittest.TestCase):
    REPRESENTATIVE_LEGACY_DECKS = (
        "HF/H2O_RHF-HF_ENERGY.inp",
        "MP2/h2o_ump2_6-31g.inp",
        "MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.inp",
        "NMR/H2O_RHF-NMR.inp",
        "OPT/H2O_RHF-DFT_OPTIMIZE_OQP.inp",
        "OPT/C2H4_BHHLYP-MRSFTDDFT_TCI_OQP.inp",
        "OPT/C2H4_BHHLYP-MRSFTDDFT_BAEKA_OQP.inp",
    )

    @staticmethod
    def _write(path, text=""):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def test_collection_keeps_only_representative_paired_legacy_decks(self):
        class NoopRunner:
            pass

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            for relative in self.REPRESENTATIVE_LEGACY_DECKS:
                legacy = root / relative
                self._write(legacy)
                self._write(legacy.with_suffix(".oqp"))

            ordinary_legacy = root / "HF/H2O_RHF-HF_GRADIENT.inp"
            self._write(ordinary_legacy)
            self._write(ordinary_legacy.with_suffix(".oqp"))
            self._write(root / "HF/legacy_only.inp")
            self._write(root / "HF/ignored.resolved.oqp")

            with load_oqp_tester(
                NoopRunner, "oqp_tester_collection_under_test"
            ) as tester_module:
                tester = tester_module.OQPTester.__new__(tester_module.OQPTester)
                tester.base_test_dir = str(root)
                selected = {
                    Path(path).relative_to(root).as_posix()
                    for path in tester._get_input_files(str(root))
                }

        expected = set(self.REPRESENTATIVE_LEGACY_DECKS)
        expected.update(
            str(Path(relative).with_suffix(".oqp"))
            for relative in self.REPRESENTATIVE_LEGACY_DECKS
        )
        expected.update({
            "HF/H2O_RHF-HF_GRADIENT.oqp",
            "HF/legacy_only.inp",
        })
        self.assertEqual(selected, expected)

    def test_paired_legacy_run_has_distinct_project_and_log(self):
        class RecordingRunner:
            calls = []

            def __init__(self, **kwargs):
                self.calls.append(kwargs)

            def run(self, test_mod=False):
                pass

            def test(self):
                return "matched", 0

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            legacy = root / "HF/H2O_RHF-HF_ENERGY.inp"
            semantic = legacy.with_suffix(".oqp")
            self._write(legacy)
            self._write(semantic)
            output_dir = root / "output"
            output_dir.mkdir()

            with load_oqp_tester(
                RecordingRunner, "oqp_tester_project_names_under_test"
            ) as tester_module:
                tester = tester_module.OQPTester.__new__(tester_module.OQPTester)
                tester.output_dir = str(output_dir)
                tester.mpi_manager = types.SimpleNamespace(use_mpi=0, rank=0)
                legacy_result = tester.run_single_test(str(legacy))
                semantic_result = tester.run_single_test(str(semantic))

        self.assertEqual(legacy_result["project"], "H2O_RHF-HF_ENERGY__legacy")
        self.assertEqual(
            Path(legacy_result["log_file"]).name,
            "H2O_RHF-HF_ENERGY__legacy.log",
        )
        self.assertEqual(semantic_result["project"], "H2O_RHF-HF_ENERGY")
        self.assertEqual(
            Path(semantic_result["log_file"]).name,
            "H2O_RHF-HF_ENERGY.log",
        )
        self.assertEqual(
            [call["project"] for call in RecordingRunner.calls],
            ["H2O_RHF-HF_ENERGY__legacy", "H2O_RHF-HF_ENERGY"],
        )

    def test_explicit_source_example_path_keeps_legacy_matrix(self):
        class NoopRunner:
            pass

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            installed_examples = root / "installed/share/examples"
            source_examples = root / "checkout/examples"
            legacy = source_examples / "HF/H2O_RHF-HF_ENERGY.inp"
            self._write(legacy)
            self._write(legacy.with_suffix(".oqp"))

            with load_oqp_tester(
                NoopRunner, "oqp_tester_explicit_source_path_under_test"
            ) as tester_module:
                tester = tester_module.OQPTester.__new__(tester_module.OQPTester)
                tester.base_test_dir = str(installed_examples)
                selected = {
                    Path(path).relative_to(source_examples).as_posix()
                    for path in tester._get_input_files(str(source_examples))
                }

        self.assertEqual(
            selected,
            {
                "HF/H2O_RHF-HF_ENERGY.inp",
                "HF/H2O_RHF-HF_ENERGY.oqp",
            },
        )

    def test_controlled_runner_modules_are_restored(self):
        class NoopRunner:
            pass

        module_name = "oqp_tester_module_restoration_under_test"
        tracked = ("oqp", "oqp.utils", "oqp.pyoqp", "oqp.runtime", module_name)
        missing = object()
        before = {name: sys.modules.get(name, missing) for name in tracked}

        with load_oqp_tester(NoopRunner, module_name):
            self.assertIs(sys.modules["oqp.pyoqp"].Runner, NoopRunner)
            self.assertIn(module_name, sys.modules)

        for name, module in before.items():
            if module is missing:
                self.assertNotIn(name, sys.modules)
            else:
                self.assertIs(sys.modules.get(name), module)


class TestOptionalGeometricExample(unittest.TestCase):
    def _run_with_error(self, error):
        class FailingRunner:
            def __init__(self, **kwargs):
                pass

            def run(self, test_mod=False):
                raise error

        with load_oqp_tester(
            FailingRunner,
            f"oqp_tester_geometric_{type(error).__name__}_{id(error)}",
        ) as tester_module:
            with tempfile.TemporaryDirectory() as output_dir:
                tester = tester_module.OQPTester.__new__(tester_module.OQPTester)
                tester.output_dir = output_dir
                tester.mpi_manager = types.SimpleNamespace(use_mpi=0, rank=0)
                return tester.run_single_test(
                    str(EXAMPLES_OPT / "HCN_RHF-DFT_CONSTRAINED_OQP.inp")
                )

    def test_missing_optional_geometric_backend_is_still_a_recognized_adapter_error(self):
        result = self._run_with_error(ModuleNotFoundError(
            "geomeTRIC is required for [optimize] lib=geometric. "
            "Install it with `pip install geometric`."
        ))

        self.assertEqual(result["status"], "SKIPPED")
        self.assertIn("optional geomeTRIC optimizer", result["message"])

    def test_unrelated_import_error_remains_visible(self):
        result = self._run_with_error(
            ModuleNotFoundError("No module named 'unrelated_backend'")
        )

        self.assertEqual(result["status"], "ERROR")
        self.assertIn("unrelated_backend", result["message"])


if __name__ == "__main__":
    unittest.main()
