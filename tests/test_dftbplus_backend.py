import importlib.util
import sys
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "dftbplus"


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_minimal_oqp_stubs():
    oqp_stub = sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    oqp_stub.ffi = types.SimpleNamespace()
    oqp_stub.lib = types.SimpleNamespace()
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    if "numpy" not in sys.modules:
        sys.modules["numpy"] = types.ModuleType("numpy")
    sys.modules.setdefault("oqp.periodic_table", types.ModuleType("oqp.periodic_table"))
    sys.modules["oqp.periodic_table"].MASSES = {1: 1.00782503223, 8: 15.99491461957}
    sys.modules["oqp.periodic_table"].SYMBOL_MAP = {"H": 1, "O": 8}
    constants = types.ModuleType("oqp.utils.constants")
    constants.ANGSTROM_TO_BOHR = 1.8897259886
    sys.modules["oqp.utils.constants"] = constants
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        use_mpi = False
        size = 1

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils


class DFTBPlusParserTests(unittest.TestCase):
    def setUp(self):
        self.dftbplus = load_module("dftbplus_under_test", "pyoqp/oqp/library/dftbplus.py")

    def test_parse_detailed_out_energy(self):
        result = self.dftbplus.parse_detailed_out(FIXTURES / "detailed.out")
        self.assertAlmostEqual(result.energy, -4.0611433278)

    def test_parse_results_tag_energy_and_gradient(self):
        result = self.dftbplus.parse_results_tag(FIXTURES / "results.tag")
        self.assertAlmostEqual(result.energy, -4.0611433278)
        self.assertEqual(len(result.gradient), 3)
        self.assertEqual(result.gradient[0], [-0.001, -0.0, 0.001])
        self.assertEqual(result.gradient[2], [-0.0, -0.003, -0.0])

    def test_parse_results_tag_with_padded_dftbplus_headers(self):
        text = """mermin_energy       :real:0:
 -0.670550865464077E+000
forces              :real:2:3,2
 -0.0 -0.0 -0.208705952983491E-002
 -0.0 -0.0  0.208705952983491E-002
"""
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "results.tag"
            path.write_text(text)
            result = self.dftbplus.parse_results_tag(path)
        self.assertAlmostEqual(result.energy, -0.670550865464077)
        self.assertEqual(len(result.gradient), 2)
        self.assertAlmostEqual(result.gradient[0][2], 0.00208705952983491)
        self.assertAlmostEqual(result.gradient[1][2], -0.00208705952983491)


class DFTBPlusSchemaTests(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()

    def test_method_dftb_and_dftb_section_are_validated(self):
        schema_text = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()
        input_checker = load_module("input_checker_dftb_under_test", "pyoqp/oqp/utils/input_checker.py")
        self.assertIn("'dftb': {", schema_text)
        self.assertIn("'sk_path': {'type': str", schema_text)

        config = {
            "input": {"method": "dftb", "runtype": "grad", "system": "\nH 0 0 0\nH 0 0 0.74"},
            "dftb": {"sk_path": "/tmp/3ob-3-1", "scc": True},
        }
        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

    def test_method_dftb_rejects_unsupported_runtype(self):
        input_checker = load_module("input_checker_dftb_runtype_under_test", "pyoqp/oqp/utils/input_checker.py")
        config = {"input": {"method": "dftb", "runtype": "hess", "system": "H 0 0 0"}}
        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("DFTB+ backend currently supports runtype=energy or grad", report.to_text())

    def test_method_dftb_allows_ground_state_gradient_index(self):
        input_checker = load_module("input_checker_dftb_grad_under_test", "pyoqp/oqp/utils/input_checker.py")
        config = {
            "input": {"method": "dftb", "runtype": "grad", "system": "\nH 0 0 0\nH 0 0 0.74"},
            "properties": {"grad": [0]},
            "dftb": {"sk_path": "/tmp/3ob-3-1"},
        }
        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())


class DFTBPlusWriterRunnerTests(unittest.TestCase):
    def setUp(self):
        self.dftbplus = load_module("dftbplus_writer_under_test", "pyoqp/oqp/library/dftbplus.py")

    def test_write_dftb_input_contains_geometry_driver_and_slakos(self):
        atoms = [8, 1, 1]
        coords_bohr = [0.0, 0.0, 0.0, 0.0, 1.43, 1.1, 0.0, -1.43, 1.1]
        config = {"input": {"charge": 0}, "dftb": {"sk_path": "/opt/dftb/3ob-3-1", "scc": True, "max_scc_iterations": 77}}
        with TemporaryDirectory() as tmp:
            self.dftbplus.write_dftbplus_input(tmp, atoms, coords_bohr, config, gradient=True)
            text = Path(tmp, "dftb_in.hsd").read_text()
        self.assertIn("Driver = ConjugateGradient { MaxSteps = 0 }", text)
        self.assertIn('Prefix = "/opt/dftb/3ob-3-1/"', text)
        self.assertIn("MaxSCCIterations = 77", text)
        self.assertIn('O = "p"', text)
        self.assertIn('H = "s"', text)
        self.assertIn("O H", text)

    def test_runner_missing_executable_has_clear_error(self):
        runner = self.dftbplus.DFTBPlusRunner({"dftb": {"executable": "/definitely/missing/dftb+", "sk_path": "/tmp/sk"}})
        with self.assertRaisesRegex(self.dftbplus.DFTBPlusError, "DFTB\\+ executable not found"):
            runner.run([1, 1], [0.0, 0.0, 0.0, 0.0, 0.0, 1.4], gradient=False)

    def test_runner_missing_parameters_has_clear_error(self):
        runner = self.dftbplus.DFTBPlusRunner({"dftb": {"executable": sys.executable, "sk_path": "/definitely/missing/sk"}})
        with self.assertRaisesRegex(self.dftbplus.DFTBPlusError, "DFTB\\+ parameter directory not found"):
            runner.run([1, 1], [0.0, 0.0, 0.0, 0.0, 0.0, 1.4], gradient=False)

    def test_runner_keep_workdir_preserves_generated_inputs_and_parses_outputs(self):
        fake_code = """#!/usr/bin/env python3
from pathlib import Path
Path('results.tag').write_text('total_energy:real:0:\\n-1.25\\nforces:real:2:2,3\\n0.1 0.0 -0.1\\n0.0 0.2 0.0\\n')
"""
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fake_dftb = tmp_path / "fake-dftbplus.py"
            fake_dftb.write_text(fake_code)
            fake_dftb.chmod(0o755)
            sk_path = tmp_path / "sk"
            sk_path.mkdir()
            runner = self.dftbplus.DFTBPlusRunner({"dftb": {"executable": str(fake_dftb), "sk_path": str(sk_path), "keep_workdir": True}})
            result = runner.run([1, 1], [0.0, 0.0, 0.0, 0.0, 0.0, 1.4], gradient=True)
            self.assertAlmostEqual(result.energy, -1.25)
            self.assertEqual(result.gradient, [[-0.1, -0.0, 0.1], [-0.0, -0.2, -0.0]])
            self.assertIsNotNone(result.workdir)
            kept_dir = Path(result.workdir)
            self.assertTrue(kept_dir.is_dir())
            self.assertTrue((kept_dir / "dftb_in.hsd").exists())
            self.assertTrue((kept_dir / "results.tag").exists())


if __name__ == "__main__":
    unittest.main()
