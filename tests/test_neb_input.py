import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_minimal_oqp_stubs():
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils


class TestNEBInputSchema(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()

    def test_neb_config_schema_defines_product_and_image_count(self):
        text = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()

        self.assertIn("'neb': {", text)
        self.assertIn("'product': {'type': str, 'default': ''}", text)
        self.assertIn("'nimage': {'type': int, 'default': '5'}", text)

    def test_full_input_checker_accepts_hf_s0_geometric_neb(self):
        input_checker = load_module(
            "input_checker_neb_valid_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            product = Path(tmpdir) / "product.xyz"
            product.write_text("2\nproduct\nH 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
            config = {
                "input": {
                    "runtype": "neb",
                    "method": "hf",
                    "basis": "3-21g",
                    "system": "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7",
                },
                "guess": {},
                "scf": {"type": "rhf", "multiplicity": 1},
                "tdhf": {},
                "properties": {},
                "optimize": {"lib": "geometric", "istate": 0},
                "neb": {"product": str(product), "nimage": 3},
            }

            report = input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())

    def test_neb_requires_geometric_optimizer_product_and_three_images(self):
        input_checker = load_module(
            "input_checker_neb_invalid_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {
                "runtype": "neb",
                "method": "hf",
                "basis": "3-21g",
                "system": "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7",
            },
            "guess": {},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {},
            "properties": {},
            "optimize": {"lib": "scipy", "istate": 0},
            "neb": {"product": "", "nimage": 2},
        }

        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        text = report.to_text()

        self.assertFalse(report.ok)
        self.assertIn("optimize.lib", text)
        self.assertIn("neb.product", text)
        self.assertIn("neb.nimage", text)


if __name__ == "__main__":
    unittest.main()
