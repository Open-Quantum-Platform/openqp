import importlib.util
import os
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


def install_geometric_stubs():
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.library", types.ModuleType("oqp.library"))
    libscipy = types.ModuleType("oqp.library.libscipy")
    libscipy.MECIOpt = object
    libscipy.MECPOpt = object
    libscipy.StateSpecificOpt = object
    sys.modules["oqp.library.libscipy"] = libscipy
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    file_utils = types.ModuleType("oqp.utils.file_utils")
    file_utils.dump_log = lambda *args, **kwargs: None
    sys.modules["oqp.utils.file_utils"] = file_utils


class FakeMol:
    def __init__(self, log_path):
        self.log_path = str(log_path)
        self.config = {
            "neb": {"product": "product.xyz", "nimage": 4},
            "geometric": {"prefix": "neb-test"},
        }


class TestNEBDispatchScaffold(unittest.TestCase):
    def setUp(self):
        install_geometric_stubs()

    def test_geometric_neb_scaffold_plans_isolated_image_directories(self):
        libgeometric = load_module(
            "libgeometric_neb_scaffold_under_test",
            "pyoqp/oqp/library/libgeometric.py",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            optimizer = libgeometric.GeometricNEBOpt(FakeMol(tmpdir))
            image_dirs = optimizer.plan_image_directories()

        self.assertEqual(
            [Path(path).name for path in image_dirs],
            ["image_000", "image_001", "image_002", "image_003"],
        )
        self.assertTrue(all(str(Path(tmpdir) / "neb") in path for path in image_dirs))

    def test_geometric_neb_resolves_repo_relative_product_paths_without_doubling_input_dir(self):
        libgeometric = load_module(
            "libgeometric_neb_product_path_under_test",
            "pyoqp/oqp/library/libgeometric.py",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            product = root / "examples" / "OPT" / "product.xyz"
            product.parent.mkdir(parents=True)
            product.write_text("2\nproduct\nH 0 0 0\nH 0 0 1\n")
            fake = FakeMol(root / "examples" / "OPT")
            fake.input_file = str(root / "examples" / "OPT" / "neb.inp")
            fake.config["neb"]["product"] = "examples/OPT/product.xyz"
            optimizer = libgeometric.GeometricNEBOpt(fake)

            cwd = os.getcwd()
            try:
                os.chdir(root)
                resolved = optimizer._resolve_product_xyz_path()
            finally:
                os.chdir(cwd)

        self.assertEqual(Path(resolved).resolve(), product.resolve())

    def test_dispatch_maps_geometric_neb_to_scaffold(self):
        runfunc_text = (ROOT / "pyoqp/oqp/library/runfunc.py").read_text()
        pyoqp_text = (ROOT / "pyoqp/oqp/pyoqp.py").read_text()

        self.assertIn("GeometricNEBOpt", runfunc_text)
        self.assertIn("'neb': GeometricNEBOpt", runfunc_text)
        self.assertIn("'neb': compute_geom", pyoqp_text)


if __name__ == "__main__":
    unittest.main()
