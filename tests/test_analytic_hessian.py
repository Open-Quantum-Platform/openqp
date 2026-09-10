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

    def test_vibrational_intensity_opt_out_keeps_harmonic_analysis_separate(self):
        class Mol:
            config = {
                "guess": {"save_mol": False},
                "properties": {"export": False, "title": ""},
                "tests": {"exception": True},
                "hess": {"type": "analytical", "state": 1, "read": False,
                         "restart": False, "temperature": [298.15], "clean": True,
                         "vibrational_intensities": False},
                "input": {"method": "tdhf"},
                "scf": {"multiplicity": 3},
                "tdhf": {"type": "mrsf", "multiplicity": 1},
            }

        hessian = self.single_point.Hessian(Mol())
        hessian._compute_vibrational_intensities(np.ones((1, 3)))

        self.assertEqual(
            hessian.mol.vibrational_intensity_metadata["status"],
            "not_computed",
        )
        self.assertIn(
            "vibrational_intensities=False",
            hessian.mol.vibrational_intensity_metadata["reason"],
        )
        self.assertEqual(hessian.mol.infrared_intensities.size, 0)
        self.assertEqual(hessian.mol.raman_activities.size, 0)


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

    def _hess_grid_report(self, functional, pruned, rad_npts, ang_npts,
                          hess_type="analytical", method="tdhf", state=1):
        config = {
            "input": {"method": method, "runtype": "hess",
                      "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3",
                      "basis": "sto-3g", "functional": functional},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 2, "multiplicity": 1},
            "hess": {"type": hess_type, "state": state, "nproc": 1,
                     "temperature": [298.15]},
            "dftgrid": {"pruned": pruned, "rad_npts": rad_npts,
                        "ang_npts": ang_npts},
        }
        return self.input_checker.check_input_values(
            config, raise_error=False, emit=False
        )

    @staticmethod
    def _grid_warnings(report):
        return [d for d in report.diagnostics
                if d.severity == "WARNING" and d.path == "dftgrid"]

    def test_analytic_tddft_hessian_warns_on_a_grid_too_coarse_for_it(self):
        # The analytic TDDFT Hessian needs a finer, unpruned grid than the rest
        # of the derivative stack. Measured on H2O/STO-3G SVWN S1: the analytic
        # frequencies move 6.03 cm-1 between the default pruned SG2 96x302 grid
        # and an unpruned 128x590 grid, while the finite-difference Hessian is
        # identical to 0.01 cm-1 on both. The warning must fire wherever that
        # error is still of that size, and must stay silent once the grid is
        # good enough -- otherwise it is either useless or noise.
        for pruned, rad, ang, why in (
            ("SG2", 96, 302, "the shipped default grid: 6.0 cm-1 error"),
            ("", 96, 302, "unpruned but still coarse: 1.5 cm-1 error"),
            ("SG1", 200, 974, "pruned, however fine the nominal counts"),
        ):
            report = self._hess_grid_report("svwn", pruned, rad, ang)
            with self.subTest(pruned=pruned, rad=rad, ang=ang):
                self.assertTrue(self._grid_warnings(report), why)
                # A warning, never an error: the number is usable and converges.
                self.assertTrue(report.ok, report.to_text())

    def test_analytic_tddft_hessian_is_quiet_on_an_adequate_grid(self):
        for pruned, rad, ang in (("", 128, 590), ("", 155, 974)):
            report = self._hess_grid_report("svwn", pruned, rad, ang)
            with self.subTest(rad=rad, ang=ang):
                self.assertFalse(self._grid_warnings(report))

    def test_grid_warning_is_scoped_to_the_analytic_dft_hessian(self):
        # Pure TDHF has no quadrature at all, and the numerical Hessian is
        # already converged on the default grid, so neither may be warned about.
        self.assertFalse(self._grid_warnings(
            self._hess_grid_report("", "SG2", 96, 302)))
        self.assertFalse(self._grid_warnings(
            self._hess_grid_report("svwn", "SG2", 96, 302,
                                   hess_type="numerical")))
        # The ground-state HF/DFT analytic Hessian is a different, older kernel
        # whose grid behaviour was not measured here, so it must not pick up an
        # excited-state warning it knows nothing about.
        self.assertFalse(self._grid_warnings(
            self._hess_grid_report("pbe", "SG2", 96, 302,
                                   method="hf", state=0)))

    def test_only_real_pruning_schemes_count_as_pruned(self):
        # source/dftlib/dft.F90 pruning is `select case` over SG0/SG1/SG2/SG3
        # with no `case default`, so any other spelling leaves the grid
        # unpruned and must not be reported as pruned.
        for pruned in ("SG0", "SG1", "SG2", "SG3"):
            with self.subTest(pruned=pruned, expect="warn"):
                self.assertTrue(self._grid_warnings(
                    self._hess_grid_report("svwn", pruned, 128, 590)))
        for pruned in ("", "none", "off", "false", "no"):
            with self.subTest(pruned=pruned, expect="quiet"):
                self.assertFalse(self._grid_warnings(
                    self._hess_grid_report("svwn", pruned, 128, 590)))

    def test_cam_mode_is_rejected_for_the_excited_state_analytic_hessian(self):
        # [dftgrid] cam_flag switches on range separation independently of the
        # functional name, and the native gate aborts on it, so a name-only
        # check lets functional=pbe + cam_flag=true validate and then die in
        # Fortran after the SCF and response have run.
        base = {
            "input": {"method": "tdhf", "functional": "pbe"},
            "scf": {"type": "rhf"},
            "tdhf": {"type": "rpa", "nstate": 2},
            "hess": {"state": 1},
        }
        for cam, expected in ((False, "supported"),
                              (True, "unsupported_feature"),
                              ("true", "unsupported_feature")):
            config = {k: v.copy() for k, v in base.items()}
            config["dftgrid"] = {"cam_flag": cam}
            status, reason = self.input_checker.analytic_hessian_capability(config)
            with self.subTest(cam_flag=cam):
                self.assertEqual(status, expected, reason)

    def test_cam_rejection_does_not_touch_paths_that_support_it(self):
        # Pure TDHF has no XC at all, and the ground-state analytic Hessian
        # supports CAM (tests/test_cam_hessian.py), so neither may be rejected.
        for method, functional, state in (("tdhf", "", 1), ("hf", "pbe", 0)):
            config = {
                "input": {"method": method, "functional": functional},
                "scf": {"type": "rhf"},
                "tdhf": {"type": "rpa", "nstate": 2},
                "hess": {"state": state},
                "dftgrid": {"cam_flag": True},
            }
            status, reason = self.input_checker.analytic_hessian_capability(config)
            with self.subTest(method=method, functional=functional):
                self.assertEqual(status, "supported", reason)

    def _hess_errors_with_ranks(self, ranks, method="tdhf", hess_type="analytical"):
        real = self.input_checker.MPIManager
        self.input_checker.MPIManager = lambda: types.SimpleNamespace(
            size=ranks, use_mpi=int(ranks > 1), rank=0
        )
        try:
            config = {
                "input": {"method": method, "functional": "svwn"},
                "scf": {"type": "rhf", "multiplicity": 1},
                "tdhf": {"type": "rpa", "multiplicity": 1, "nstate": 2},
                "hess": {"type": hess_type,
                         "state": 1 if method == "tdhf" else 0, "nproc": 1},
                "dftgrid": {"pruned": "", "rad_npts": 128, "ang_npts": 590},
            }
            report = self.input_checker.CheckReport()
            self.input_checker._check_hess(config, report)
            return [d for d in report.diagnostics
                    if d.severity == "ERROR" and "one MPI rank" in d.message]
        finally:
            self.input_checker.MPIManager = real

    def test_multi_rank_excited_state_analytic_hessian_is_rejected(self):
        # tdhf_hessian_is_applicable requires mpi_size == 1 and aborts
        # otherwise, so validation must catch this rather than letting the run
        # die in Fortran after the SCF and response are already done.
        self.assertFalse(self._hess_errors_with_ranks(1))
        self.assertTrue(self._hess_errors_with_ranks(4))
        # The ground-state Hessian and the numerical path have no such limit.
        self.assertFalse(self._hess_errors_with_ranks(4, method="hf"))
        self.assertFalse(self._hess_errors_with_ranks(4, hess_type="numerical"))

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

    def test_mrsf_tdhf_analytical_hessian_is_supported_without_fallback(self):
        config = {
            "input": {"method": "tdhf", "runtype": "hess", "system": "\nO 0 0 0\nH 0 0 0.9\nH 0 0.7 -0.3", "basis": "sto-3g"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3, "multiplicity": 3},
            "hess": {"type": "analytical", "state": 1, "nproc": 1, "temperature": [298.15]},
        }

        report = self.input_checker.check_input_values(config, raise_error=False, emit=False)

        self.assertTrue(report.ok, report.to_text())

    def test_mrsf_tddft_semilocal_analytical_hessian_is_supported(self):
        config = {
            "input": {"method": "tdhf", "functional": "B3LYP"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "multiplicity": 1},
            "hess": {"state": 1},
        }
        status, reason = self.input_checker.analytic_hessian_capability(config)
        self.assertEqual(status, "supported", reason)
        self.assertIn("MRSF-TDDFT", reason)

        for functional in ("M06-L", "CAM-B3LYP"):
            config["input"]["functional"] = functional
            status, reason = self.input_checker.analytic_hessian_capability(config)
            self.assertEqual(status, "unsupported_feature")
            self.assertIn("remain fail-closed", reason)

    def test_mrsf_native_ts_analytical_initial_hessian_is_supported(self):
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

        self.assertTrue(report.ok, report.to_text())
        self.assertNotIn("not implemented", report.to_text())

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

    def test_analytic_rpa_hessian_examples_compute_an_isolation_root(self):
        examples = sorted((ROOT / "examples/HESS").glob("*_RPA_ANA_HESS.inp"))

        self.assertTrue(examples)
        for example in examples:
            nstate_lines = [
                line for line in example.read_text().splitlines()
                if line.strip().lower().startswith("nstate=")
            ]
            with self.subTest(example=example.name):
                self.assertEqual(len(nstate_lines), 1)
                self.assertGreaterEqual(int(nstate_lines[0].split("=", 1)[1]), 2)


if __name__ == "__main__":
    unittest.main()
