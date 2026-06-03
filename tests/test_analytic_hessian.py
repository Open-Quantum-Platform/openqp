import importlib.util
import sys
import types
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

STUB_MODULES = [
    "oqp",
    "oqp.molecule",
    "oqp.utils",
    "oqp.utils.constants",
    "oqp.utils.mpi_utils",
    "oqp.utils.matrix",
    "oqp.utils.file_utils",
    "oqp.library",
    "oqp.library.frequency",
    "oqp.periodic_table",
]


def snapshot_modules(names=STUB_MODULES):
    return {name: sys.modules.get(name) for name in names}


def restore_modules(snapshot):
    for name, module in snapshot.items():
        if module is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = module


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class AnalyticHessianNativeDispatchTests(unittest.TestCase):
    def setUp(self):
        self._module_snapshot = snapshot_modules()
        oqp_stub = types.ModuleType("oqp")
        self.native_calls = []

        def hf_hessian(mol):
            self.native_calls.append(mol)

        setattr(oqp_stub, "hf_hessian", hf_hessian)
        sys.modules["oqp"] = oqp_stub
        molecule_mod = types.ModuleType("oqp.molecule")
        setattr(molecule_mod, "Molecule", object)
        sys.modules["oqp.molecule"] = molecule_mod
        sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
        mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

        class MPIManager:
            use_mpi = False
            rank = 0
            size = 1

            def barrier(self):
                return None

        setattr(mpi_utils, "MPIManager", MPIManager)
        setattr(mpi_utils, "MPIPool", object)
        sys.modules["oqp.utils.mpi_utils"] = mpi_utils
        matrix = types.ModuleType("oqp.utils.matrix")
        setattr(matrix, "DampingParam", object)
        setattr(matrix, "DispersionModel", object)
        sys.modules["oqp.utils.matrix"] = matrix
        library = types.ModuleType("oqp.library")
        sys.modules["oqp.library"] = library
        frequency = types.ModuleType("oqp.library.frequency")
        setattr(frequency, "normal_mode", lambda *args, **kwargs: (np.array([]), np.array([]), np.array([])))
        setattr(frequency, "thermal_analysis", lambda *args, **kwargs: {})
        sys.modules["oqp.library.frequency"] = frequency
        file_utils = types.ModuleType("oqp.utils.file_utils")
        setattr(file_utils, "dump_log", lambda *args, **kwargs: None)
        setattr(file_utils, "dump_data", lambda *args, **kwargs: None)
        setattr(file_utils, "write_config", lambda *args, **kwargs: None)
        setattr(file_utils, "write_xyz", lambda *args, **kwargs: None)
        sys.modules["oqp.utils.file_utils"] = file_utils
        self.single_point = load_module("single_point_analytic_hess_dispatch", "pyoqp/oqp/library/single_point.py")

    def tearDown(self):
        sys.modules.pop("single_point_analytic_hess_dispatch", None)
        restore_modules(self._module_snapshot)

    def test_hf_analytical_hessian_reads_native_fortran_hessian_without_external_backend(self):
        class Mol:
            config = {
                "guess": {"save_mol": False},
                "properties": {"export": False, "title": ""},
                "tests": {"exception": True},
                "hess": {"type": "analytical", "state": 0, "read": False, "restart": False, "temperature": [298.15], "clean": True},
                "input": {"method": "hf"},
                "scf": {"multiplicity": 1},
                "tdhf": {"type": "rpa", "multiplicity": 1},
            }
            data = {"natom": 2, "OQP::hf_hessian": np.eye(6)}

            def set_hessian_result(self, raw_hessian):
                self.hessian = np.asarray(raw_hessian, dtype=float)
                self.hessian_metadata = {"max_asymmetry": 0.0, "symmetrized": False}
                return self.hessian

        mol = Mol()
        hessian = self.single_point.Hessian(mol)

        result, flags = hessian.analytical_ground_state_hess()

        self.assertEqual(self.native_calls, [mol])
        self.assertEqual(flags, ["computed", "native_openqp"])
        self.assertEqual(result.shape, (6, 6))
        self.assertEqual(mol.hessian_metadata["backend"], "native_openqp")
        self.assertTrue(mol.hessian_metadata["native_openqp_kernel"])
        self.assertTrue(mol.hessian_metadata["native_openqp_cphf_solver_exercised"])
        self.assertTrue(mol.hessian_metadata["native_openqp_final_assembly"])
        self.assertTrue(mol.hessian_metadata["no_external_hessian_backend"])
        self.assertNotIn("reference_backend", mol.hessian_metadata)

    def test_sf_analytical_hessian_routes_separately_from_mrsf_private_path(self):
        class Mol:
            config = {
                "guess": {"save_mol": False},
                "properties": {"export": False, "title": ""},
                "tests": {"exception": True},
                "hess": {"type": "analytical", "state": 1, "read": False, "restart": False, "temperature": [298.15], "clean": True},
                "input": {"method": "tdhf"},
                "scf": {"multiplicity": 3},
                "tdhf": {"type": "sf", "multiplicity": 1},
            }

        hessian = self.single_point.Hessian(Mol())
        hessian.analytical_sf_hess = lambda: ("sf-route", ["stubbed"])

        self.assertEqual(hessian.analytical_hess(), ("sf-route", ["stubbed"]))


class AnalyticHessianInputValidationTests(unittest.TestCase):
    def setUp(self):
        self._module_snapshot = snapshot_modules()
        sys.modules.setdefault("oqp", types.ModuleType("oqp"))
        sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
        mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

        class MPIManager:
            use_mpi = False
            size = 1

        setattr(mpi_utils, "MPIManager", MPIManager)
        sys.modules["oqp.utils.mpi_utils"] = mpi_utils
        self.input_checker = load_module("input_checker_analytic_hess", "pyoqp/oqp/utils/input_checker.py")

    def tearDown(self):
        sys.modules.pop("input_checker_analytic_hess", None)
        restore_modules(self._module_snapshot)

    def test_hf_analytical_hessian_is_allowed_by_capability_matrix(self):
        config = {
            "input": {"method": "hf", "runtype": "hess", "system": "\nH 0 0 0\nH 0 0 0.74", "basis": "sto-3g"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "multiplicity": 1, "nstate": 1},
            "hess": {"type": "analytical", "state": 0, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())

    def test_rohf_hf_analytical_hessian_is_supported(self):
        # ROHF (and UHF) HF/DFT analytic Hessians are implemented and
        # finite-difference validated, so the capability matrix accepts them.
        config = {
            "input": {"method": "hf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "rpa", "nstate": 1, "multiplicity": 1},
            "hess": {"type": "analytical", "state": 0, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())

    def test_hf_analytical_hessian_allows_openqp_library_basis_mapping(self):
        config = {
            "input": {
                "method": "hf",
                "runtype": "hess",
                "system": "\nO 0 0 0 o1\nH 0 0 0.9 h1\nH 0 0.7 -0.3 h2",
                "basis": "library",
                "library": "\no1 sto-3g\nh1 sto-3g\nh2 sto-3g",
            },
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 1, "multiplicity": 1},
            "hess": {"type": "analytical", "state": 0, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())

    def test_tddft_analytical_hessian_is_rejected_until_scaffold_exists(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 3, "multiplicity": 1},
            "hess": {"type": "analytical", "state": 1, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertFalse(report.ok)
        self.assertIn("TDDFT analytic Hessian is not implemented", report.to_text())

    def test_mrsf_analytical_hessian_is_rejected_explicitly_not_silently_numerical(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3, "multiplicity": 3},
            "hess": {"type": "analytical", "state": 1, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertFalse(report.ok)
        self.assertIn("MRSF-TDDFT analytic Hessian is not implemented", report.to_text())

    def test_sf_analytical_hessian_has_sf_specific_rejection_message(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "sf", "nstate": 3, "multiplicity": 3},
            "hess": {"type": "analytical", "state": 1, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)
        text = report.to_text()

        self.assertFalse(report.ok)
        self.assertIn("SF-TDDFT analytic Hessian is not implemented", text)
        self.assertNotIn("MRSF gradient/Z-vector", text)

    def test_hf_analytic_hessian_example_documents_keyword(self):
        example = ROOT / "examples/HESS/H2O_RHF-DFT_ANA_HESS.inp"
        text = example.read_text()

        self.assertIn("runtype=hess", text)
        self.assertIn("type=analytical", text)
        self.assertIn("state=0", text)


if __name__ == "__main__":
    unittest.main()
