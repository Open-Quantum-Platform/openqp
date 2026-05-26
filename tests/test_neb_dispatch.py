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


def install_neb_runtime_stubs(calls):
    geometric = types.ModuleType("geometric")
    sys.modules["geometric"] = geometric

    molecule_module = types.ModuleType("geometric.molecule")
    setattr(molecule_module, "Elements", {1: "H", 6: "C", 7: "N", 8: "O"})

    class Molecule:
        def __init__(self):
            self.elem = []
            self.xyzs = []
            self.comms = []

        def build_topology(self):
            calls.append(("build_topology", len(self.xyzs)))

    setattr(molecule_module, "Molecule", Molecule)
    sys.modules["geometric.molecule"] = molecule_module

    params_module = types.ModuleType("geometric.params")

    class NEBParams:
        plain = 0

        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)
            calls.append(("params", kwargs["images"]))

    setattr(params_module, "NEBParams", NEBParams)
    sys.modules["geometric.params"] = params_module

    neb_module = types.ModuleType("geometric.neb")

    class ElasticBand:
        def __init__(self, molecule, engine, tmpdir, params, plain):
            self.molecule = molecule
            self.engine = engine
            self.tmpdir = tmpdir
            self.params = params
            calls.append(("elastic_band", len(molecule.xyzs), tmpdir))

        def SaveToDisk(self, fout):
            calls.append(("save", fout))

    def OptimizeChain(chain, engine, params):
        calls.append(("optimize_chain", params.images))
        return chain, 3

    setattr(neb_module, "ElasticBand", ElasticBand)
    setattr(neb_module, "OptimizeChain", OptimizeChain)
    sys.modules["geometric.neb"] = neb_module


class FakeMol:
    def __init__(self, log_path):
        self.log_path = str(log_path)
        self.input_file = str(Path(log_path) / "input.inp")
        self.config = {
            "neb": {"product": "product.xyz", "nimage": 4},
            "geometric": {"prefix": "neb-test"},
        }

    def get_atoms(self):
        return [1, 1]

    def get_system(self):
        return [[0.0, 0.0, 0.0], [0.0, 0.0, 1.88972612546]]


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

    def test_geometric_neb_optimize_calls_geometric_neb_runner(self):
        neb_utils = load_module("oqp.library.neb_utils", "pyoqp/oqp/library/neb_utils.py")
        sys.modules["oqp.library.neb_utils"] = neb_utils
        calls = []
        install_neb_runtime_stubs(calls)
        libgeometric = load_module(
            "libgeometric_neb_runner_under_test",
            "pyoqp/oqp/library/libgeometric.py",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "product.xyz").write_text(
                "2\nproduct\nH 0.0 0.0 0.2\nH 0.0 0.0 1.4\n"
            )
            optimizer = libgeometric.GeometricNEBOpt(FakeMol(tmpdir))
            optimizer.natom = 2
            optimizer.maxit = 9
            optimizer.max_grad = 0.001
            optimizer.rmsd_grad = 0.0005

            final_chain = optimizer.optimize()

        self.assertIsNotNone(final_chain)
        self.assertIn(("params", 4), calls)
        self.assertIn(("optimize_chain", 4), calls)
        self.assertIn(("save", "chain_final.xyz"), calls)


if __name__ == "__main__":
    unittest.main()
