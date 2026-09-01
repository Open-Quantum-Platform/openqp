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
    "oqp.library.openqp_dftb",
    "oqp.utils.state_labels",
    "oqp.utils.qmmm",
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
        openqp_dftb = types.ModuleType("oqp.library.openqp_dftb")
        setattr(openqp_dftb, "OpenQPDFTBAdapter", object)
        sys.modules["oqp.library.openqp_dftb"] = openqp_dftb
        state_tracking = types.ModuleType("oqp.library.state_tracking")
        setattr(state_tracking, "diagonal_phase_tracking", lambda *args, **kwargs: None)
        setattr(state_tracking, "maximum_overlap_assignment", lambda *args, **kwargs: None)
        sys.modules["oqp.library.state_tracking"] = state_tracking
        nac_utils = types.ModuleType("oqp.library.nac_utils")
        for name in (
            "canonical_state_overlap",
            "hst_derivative_coupling",
            "interstate_coupling",
            "load_numerical_nac_cache",
            "write_numerical_nac_cache_marker",
        ):
            setattr(nac_utils, name, lambda *args, **kwargs: None)
        sys.modules["oqp.library.nac_utils"] = nac_utils
        tb_backends = types.ModuleType("oqp.utils.tb_backends")
        setattr(tb_backends, "is_tb_method", lambda *_args, **_kwargs: False)
        setattr(tb_backends, "make_tb_adapter", lambda *_args, **_kwargs: None)
        setattr(tb_backends, "tb_config", lambda *_args, **_kwargs: {})
        sys.modules["oqp.utils.tb_backends"] = tb_backends
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
        state_labels = types.ModuleType("oqp.utils.state_labels")
        setattr(state_labels, "is_mrsf", lambda *_args, **_kwargs: False)
        setattr(state_labels, "public_state_label", lambda state, *_args, **_kwargs: f"S{state}")
        sys.modules["oqp.utils.state_labels"] = state_labels
        sys.modules["oqp.utils.qmmm"] = types.ModuleType("oqp.utils.qmmm")
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

    def test_matrix_only_hessian_skips_vibrational_analysis_and_cache(self):
        class Mol:
            config = {
                "guess": {"save_mol": False},
                "properties": {"export": False, "title": ""},
                "tests": {"exception": True},
                "hess": {"type": "numerical", "state": 0, "read": False,
                         "restart": False, "temperature": [298.15], "clean": True,
                         "dx": 0.01, "nproc": 1},
                "input": {"method": "hf"},
                "scf": {"multiplicity": 1},
                "tdhf": {"type": "rpa", "multiplicity": 1},
            }
            energies = np.array([-1.0])

            def save_freqs(self, _state):
                raise AssertionError("matrix-only Hessian must not write a frequency cache")

        expected = np.diag(np.arange(1.0, 7.0))
        hessian = self.single_point.Hessian(Mol())
        hessian.hess_func = lambda: (expected.copy(), ["computed"])
        self.single_point.normal_mode = lambda *_args, **_kwargs: (
            (_ for _ in ()).throw(AssertionError("normal-mode analysis must be skipped"))
        )

        result = hessian.hessian(analysis=False)

        self.assertTrue(np.array_equal(result, expected))
        self.assertTrue(np.array_equal(hessian.mol.hessian, expected))


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

    def test_pure_tdhf_analytical_hessian_is_supported(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 3, "multiplicity": 1},
            "hess": {"type": "analytical", "state": 1, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())

    def test_excited_state_analytic_hessian_functional_gate(self):
        base = {
            "input": {"method": "tdhf"},
            "scf": {"type": "rhf"},
            "tdhf": {"type": "rpa", "nstate": 2},
            "hess": {"state": 1},
        }

        for functional in (
            "", "SVWN", "svwn5", "LDA", "BLYP", "PBE", "PBEPBE",
            "b3lyp5", "B3LYPV5",
        ):
            config = {section: values.copy() for section, values in base.items()}
            config["input"]["functional"] = functional
            status, reason = self.input_checker.analytic_hessian_capability(config)
            with self.subTest(functional=functional):
                self.assertEqual(status, "supported", reason)

        for functional in ("B3LYP", "M06-L", "CAM-B3LYP", "TETER"):
            config = {section: values.copy() for section, values in base.items()}
            config["input"]["functional"] = functional
            status, reason = self.input_checker.analytic_hessian_capability(config)
            with self.subTest(functional=functional):
                self.assertEqual(status, "unsupported_feature")
                self.assertIn("LDA/GGA and global-hybrid paths", reason)

    def test_excited_state_analytic_hessian_rejects_triplet_rpa_during_input_check(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess",
                      "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3",
                      "basis": "sto-3g"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 3, "multiplicity": 3},
            "hess": {"type": "analytical", "state": 1, "nproc": 1,
                     "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(
            config, raise_error=False, emit=False,
        )

        self.assertFalse(report.ok)
        self.assertIn("singlet targets only", report.to_text())

    def test_excited_state_analytic_hessian_rejects_higher_roots_until_indefinite_solver(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess",
                      "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3",
                      "basis": "sto-3g"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 3, "multiplicity": 1},
            "hess": {"type": "analytical", "state": 2, "nproc": 1,
                     "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(
            config, raise_error=False, emit=False,
        )

        self.assertFalse(report.ok)
        self.assertIn("only the lowest excited root", report.to_text())

    def test_excited_state_analytic_hessian_requires_two_computed_roots(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess",
                      "system": "\nH 0 0 -0.37\nH 0 0 0.37",
                      "basis": "sto-3g"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 1, "multiplicity": 1},
            "hess": {"type": "analytical", "state": 1, "nproc": 1,
                     "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(
            config, raise_error=False, emit=False,
        )

        self.assertFalse(report.ok)
        self.assertIn("tdhf.nstate>=2", report.to_text())

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

    def test_mrsf_native_ts_analytical_initial_hessian_is_rejected_in_preflight(self):
        config = {
            "input": {"method": "tdhf", "runtype": "ts",
                      "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3",
                      "basis": "sto-3g"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3, "multiplicity": 3},
            "optimize": {"lib": "oqp", "istate": 1},
            "oqp": {"init_hessian": "analytical"},
            "hess": {"type": "analytical", "state": 1, "nproc": 1,
                     "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(
            config, raise_error=False, emit=False,
        )

        self.assertFalse(report.ok)
        self.assertIn("Analytical native TS initialization is unavailable", report.to_text())
        self.assertIn("MRSF-TDDFT analytic Hessian is not implemented", report.to_text())

    def test_native_ts_analytical_initial_hessian_applies_basis_l_gate(self):
        config = {
            "input": {"method": "hf", "runtype": "ts",
                      "system": "\nH 0 0 0\nH 0 0 0.74", "basis": "mock-g"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 1, "multiplicity": 1},
            "optimize": {"lib": "oqp", "istate": 0},
            "oqp": {"init_hessian": "analytical"},
            "hess": {"state": 0},
        }
        self.input_checker._basis_max_angular_momentum = lambda _config: 4

        report = self.input_checker.check_input_values(
            config, raise_error=False, emit=False,
        )

        self.assertFalse(report.ok)
        self.assertIn("basis angular momentum only up to L=3", report.to_text())

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
