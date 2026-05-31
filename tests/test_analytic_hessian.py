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


class FakeHessian:
    def kernel(self):
        h = np.zeros((2, 2, 3, 3))
        h[0, 0, 0, 0] = 1.0
        h[0, 1, 0, 2] = 2.0
        h[1, 0, 2, 0] = 4.0
        h[1, 1, 1, 1] = 3.0
        return h


class FakeMF:
    def kernel(self):
        return -1.0

    def Hessian(self):
        return FakeHessian()


class AnalyticHessianExternalRuntimeTests(unittest.TestCase):
    def setUp(self):
        self._module_snapshot = snapshot_modules()
        oqp_stub = sys.modules.setdefault("oqp", types.ModuleType("oqp"))
        setattr(oqp_stub, "ffi", types.SimpleNamespace())
        setattr(oqp_stub, "lib", types.SimpleNamespace())
        sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
        constants = types.ModuleType("oqp.utils.constants")
        setattr(constants, "ANGSTROM_TO_BOHR", 1.8897259886)
        sys.modules["oqp.utils.constants"] = constants
        periodic_table = types.ModuleType("oqp.periodic_table")
        setattr(periodic_table, "SYMBOL_MAP", {1: "H", 8: "O", "H": 1, "O": 8})
        setattr(periodic_table, "ELEMENTS_NAME", {"H": "H", "O": "O"})
        sys.modules["oqp.periodic_table"] = periodic_table
        self.external = load_module("external_analytic_hess_under_test", "pyoqp/oqp/library/external.py")

    def tearDown(self):
        sys.modules.pop("external_analytic_hess_under_test", None)
        restore_modules(self._module_snapshot)

    def test_pyscf_hessian_is_flattened_symmetrized_and_labeled_for_openqp(self):
        class Mol:
            usempi = False
            project_name = "h2"
            config = {"input": {"method": "hf"}, "scf": {"type": "rhf"}}

        mol = Mol()
        hessian, flags = self.external.analytic_hessian_from_pyscf(mol, mf_factory=lambda _mol: FakeMF())

        self.assertEqual(flags, ["computed", "external_pyscf"])
        self.assertEqual(hessian.shape, (6, 6))
        self.assertAlmostEqual(hessian[0, 0], 1.0)
        self.assertAlmostEqual(hessian[0, 5], 3.0)
        self.assertAlmostEqual(hessian[5, 0], 3.0)
        self.assertAlmostEqual(hessian[4, 4], 3.0)
        self.assertEqual(mol.hessian_metadata["backend"], "external_pyscf")
        self.assertFalse(mol.hessian_metadata["native_openqp_kernel"])
        self.assertTrue(mol.hessian_metadata["no_numerical_fallback"])
        self.assertEqual(mol.hessian_metadata["shape"], [6, 6])
        self.assertAlmostEqual(mol.hessian_metadata["max_asymmetry_before_symmetrization"], 2.0)
        self.assertEqual(mol.hessian_metadata["compact_validation_summary"]["matrix_payload"], "omitted")
        self.assertNotIn("raw_hessian", mol.hessian_metadata["compact_validation_summary"])

    def test_mrsf_pyscf_hessian_dispatch_fails_explicitly(self):
        class Mol:
            usempi = False
            project_name = "h2o_mrsf"
            config = {"input": {"method": "tdhf"}, "tdhf": {"type": "mrsf"}}

        with self.assertRaisesRegex(NotImplementedError, "MRSF-TDDFT analytic Hessian"):
            self.external.analytic_hessian_from_pyscf(Mol(), mf_factory=lambda _mol: FakeMF())


    def test_design_note_documents_all_planned_hessian_tiers(self):
        design_note = ROOT / "docs/analytic_hessian_design.md"

        text = design_note.read_text()

        self.assertIn("## Units and data contract", text)
        self.assertIn("## HF analytic Hessian", text)
        self.assertIn("## DFT analytic Hessian", text)
        self.assertIn("## TDDFT analytic Hessian", text)
        self.assertIn("## MRSF-TDDFT analytic Hessian", text)
        self.assertIn("no silent numerical fallback", text.lower())


class AnalyticHessianNativeDispatchTests(unittest.TestCase):
    def setUp(self):
        self._module_snapshot = snapshot_modules()
        oqp_stub = types.ModuleType("oqp")
        self.native_calls = []

        def hf_hessian(mol):
            self.native_calls.append(mol)
            raw_hess = np.arange(36, dtype=float).reshape(6, 6)
            mol.set_hessian_result(raw_hess)

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

    def test_hf_analytical_hessian_uses_guarded_external_pyscf_bridge_not_native_scaffold(self):
        external_calls = []
        external_mod = types.ModuleType("oqp.library.external")

        def analytic_hessian_from_pyscf(mol):
            external_calls.append(mol)
            mol.hessian_metadata = {
                "backend": "external_pyscf",
                "native_openqp_kernel": False,
                "no_numerical_fallback": True,
                "shape": [6, 6],
                "max_asymmetry_before_symmetrization": 0.0,
                "compact_validation_summary": {"matrix_payload": "omitted"},
            }
            return np.eye(6), ["computed", "external_pyscf"]

        setattr(external_mod, "analytic_hessian_from_pyscf", analytic_hessian_from_pyscf)
        sys.modules["oqp.library.external"] = external_mod

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
            data = {"natom": 2}

        mol = Mol()
        hessian = self.single_point.Hessian(mol)

        result, flags = hessian.analytical_ground_state_hess()

        self.assertEqual(external_calls, [mol])
        self.assertEqual(self.native_calls, [])
        self.assertEqual(flags, ["computed", "external_pyscf"])
        self.assertEqual(result.shape, (6, 6))
        self.assertEqual(mol.hessian_metadata["backend"], "external_pyscf")
        self.assertFalse(mol.hessian_metadata["native_openqp_kernel"])

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

    def test_mrsf_analytical_hessian_guard_reports_root_mapping_and_spin(self):
        class Mol:
            config = {
                "guess": {"save_mol": False},
                "properties": {"export": False, "title": ""},
                "tests": {"exception": True},
                "hess": {"type": "analytical", "state": 2, "read": False, "restart": False, "temperature": [298.15], "clean": True},
                "input": {"method": "tdhf"},
                "scf": {"multiplicity": 3},
                "tdhf": {"type": "mrsf", "multiplicity": 1},
            }
            data = {"OQP::td_s2": np.array([2.0, 0.02, 1.04]), "natom": 3}
            hessian_metadata = {}

            @staticmethod
            def get_system():
                return np.array([
                    [0.0000000, 0.0000000, -0.0410615540],
                    [-0.5331943, 0.5331943, -0.6144692230],
                    [0.5331943, -0.5331943, -0.6144692230],
                ])

            @staticmethod
            def get_atoms():
                return np.array([8.0, 1.0, 1.0])

        mol = Mol()
        hessian = self.single_point.Hessian(mol)

        with self.assertRaisesRegex(
            NotImplementedError,
            r"MRSF-TDDFT analytic Hessian.*root 2 maps to physical S1.*<S\^2>=1\.04",
        ):
            hessian.analytical_hess()

        self.assertEqual(mol.hessian_metadata["backend"], "native_mrsf_analytical_partial")
        self.assertIn("nuclear_repulsion", mol.hessian_metadata["completed_terms"])
        self.assertIn("one_electron_second_derivative_contraction", mol.hessian_metadata["completed_terms"])
        self.assertIn("one_electron_integral_derivative_backend", mol.hessian_metadata["missing_terms"])
        self.assertIn("electronic_response", mol.hessian_metadata["missing_terms"])
        self.assertEqual(mol.hessian_metadata["completed_term_shapes"]["nuclear_repulsion"], [9, 9])
        self.assertEqual(mol.hessian_metadata["completed_term_status"]["one_electron_second_derivative_contraction"], "assembly_kernel_validated_requires_integral_tensor")

    def test_mrsf_numerical_hessian_runtime_guard_requires_root_tracking_oracle(self):
        class Mol:
            project_name = "h2o_mrsf"
            log_path = "/tmp"
            config = {
                "guess": {"save_mol": False},
                "properties": {"export": False, "title": ""},
                "tests": {"exception": True},
                "hess": {"type": "numerical", "state": 2, "read": False, "restart": False, "temperature": [298.15], "clean": True, "nproc": 1, "dx": 1.0e-3},
                "input": {"method": "tdhf"},
                "scf": {"multiplicity": 3},
                "tdhf": {"type": "mrsf", "multiplicity": 1, "nstate": 3},
            }
            data = {"natom": 3}

        hessian = self.single_point.Hessian(Mol())

        with self.assertRaisesRegex(
            NotImplementedError,
            r"MRSF-TDDFT numerical Hessian requires the Gate 3B root-tracked finite-difference oracle.*state-index-only",
        ):
            hessian.numerical_hess()

    def test_mrsf_root_tracked_oracle_mode_dispatches_to_explicit_path(self):
        class Mol:
            project_name = "h2o_mrsf"
            log_path = "/tmp"
            config = {
                "guess": {"save_mol": False},
                "properties": {"export": False, "title": ""},
                "tests": {"exception": True},
                "hess": {"type": "mrsf_numerical_oracle", "state": 2, "read": False, "restart": False, "temperature": [298.15], "clean": True, "nproc": 1, "dx": 1.0e-3},
                "input": {"method": "tdhf"},
                "scf": {"multiplicity": 3},
                "tdhf": {"type": "mrsf", "multiplicity": 1, "nstate": 3},
            }
            data = {"natom": 3}

        hessian = self.single_point.Hessian(Mol())

        self.assertEqual(hessian.hess_func.__name__, "mrsf_numerical_oracle_hess")


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

    def test_mrsf_analytical_hessian_partial_scaffold_is_allowed_for_term_validation(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3, "multiplicity": 3},
            "hess": {"type": "analytical", "state": 1, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())

    def test_mrsf_analytical_hessian_rejects_high_spin_reference_root(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3, "multiplicity": 3},
            "hess": {"type": "analytical", "state": 0, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertFalse(report.ok)
        self.assertIn("requires a positive physical root state", report.to_text())

    def test_mrsf_numerical_hessian_is_rejected_until_root_tracking_oracle_exists(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3, "multiplicity": 1},
            "hess": {"type": "numerical", "state": 2, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertFalse(report.ok)
        text = report.to_text()
        self.assertIn("MRSF-TDDFT numerical Hessian requires the Gate 3B root-tracked finite-difference oracle", text)
        self.assertIn("state-index-only finite differences", text)

    def test_mrsf_root_tracked_oracle_hessian_mode_is_allowed_by_input_checker(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3, "multiplicity": 1},
            "hess": {"type": "mrsf_numerical_oracle", "state": 2, "nproc": 1, "temperature": [298.15], "dx": 3.0e-4},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())

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
