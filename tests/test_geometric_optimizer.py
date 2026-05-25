import importlib.util
import sys
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
    sys.modules["oqp.library.libscipy"] = libscipy

    libdlfind = types.ModuleType("oqp.library.libdlfind")
    libdlfind.DLFindMin = type("DLFindMin", (), {})
    libdlfind.DLFindTS = type("DLFindTS", (), {})
    libdlfind.DLFindMECI = type("DLFindMECI", (), {})
    sys.modules["oqp.library.libdlfind"] = libdlfind

    libgeometric = types.ModuleType("oqp.library.libgeometric")

    class GeometricOpt:
        def __init__(self, mol):
            self.mol = mol

    class GeometricMECIOpt:
        def __init__(self, mol):
            self.mol = mol

    libgeometric.GeometricOpt = GeometricOpt
    libgeometric.GeometricMECIOpt = GeometricMECIOpt
    sys.modules["oqp.library.libgeometric"] = libgeometric
    return GeometricOpt, GeometricMECIOpt


class TestGeometricOptimizerConfig(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()

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
        GeometricOpt, _ = install_runfunc_stubs()
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
        _, GeometricMECIOpt = install_runfunc_stubs()
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


if __name__ == "__main__":
    unittest.main()
