import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


class TestAnalyticHessianBindings(unittest.TestCase):
    def test_c_header_declares_planned_hessian_entry_points(self):
        header = read("include/oqp.h")

        for symbol in [
            "hf_hessian",
            "tdhf_hessian",
            "tdhf_sf_hessian",
        ]:
            self.assertIn(f"void {symbol}(struct oqp_handle_t *inf);", header)
        self.assertNotIn("tdhf_mrsf_hessian", header)

    def test_python_dispatch_mentions_native_hessian_entry_points(self):
        source = read("pyoqp/oqp/library/single_point.py")

        for symbol in [
            "oqp.hf_hessian",
            "oqp.tdhf_hessian",
            "oqp.tdhf_sf_hessian",
        ]:
            self.assertIn(symbol, source)

    def test_hf_dispatch_runs_native_cphf_prepass_before_pyscf_oracle_without_numerical_fallback(self):
        source = read("pyoqp/oqp/library/single_point.py")

        self.assertIn("native_hess_func", source)
        self.assertIn("no numerical fallback", source.lower())
        self.assertIn("native OpenQP CPHF prepass", source)
        self.assertIn("external PySCF final Hessian oracle", source)
        self.assertIn("analytic_hessian_from_pyscf", source)
        self.assertIn("native_hess_func = self.native_hess_func['hf']", source)
        self.assertIn("native_hess_func(self.mol)", source)
        self.assertIn("Native Fortran CPHF Hessian response log", source)

    def test_public_sf_dispatch_is_separate_from_private_mrsf_path(self):
        source = read("pyoqp/oqp/library/single_point.py")

        self.assertIn("def analytical_sf_hess", source)
        self.assertIn("td_type == 'sf'", source)
        self.assertIn("return self.analytical_sf_hess()", source)

    def test_hf_hessian_fortran_runs_native_cphf_solver_diagnostic_without_claiming_final_assembly(self):
        source = read("source/modules/hf_hessian.F90")

        self.assertIn("module hf_hessian_mod", source)
        self.assertIn('bind(C, name="hf_hessian")', source)
        self.assertIn("subroutine hf_hessian_C", source)
        self.assertIn("subroutine hf_hessian", source)
        self.assertIn("use cphf_mod, only: cphf_solve", source)
        self.assertIn("Native OpenQP HF/DFT Hessian CPHF response prepass", source)
        self.assertIn("call cphf_solve", source)
        self.assertIn("Final analytic Hessian assembly remains guarded", source)
        self.assertNotIn("implementation is not available yet", source)

    def test_tdhf_hessian_fortran_scaffolds_export_c_abi_without_claiming_support(self):
        expected = {
            "tdhf_hessian.F90": [
                "module tdhf_hessian_mod",
                'bind(C, name="tdhf_hessian")',
                "subroutine tdhf_hessian_C",
                "subroutine tdhf_hessian",
                "Analytic TDDFT Hessian kernel scaffold reached",
            ],
            "tdhf_sf_hessian.F90": [
                "module tdhf_sf_hessian_mod",
                'bind(C, name="tdhf_sf_hessian")',
                "subroutine tdhf_sf_hessian_C",
                "subroutine tdhf_sf_hessian",
                "Analytic SF-TDDFT Hessian kernel scaffold reached",
            ],
        }

        for filename, needles in expected.items():
            with self.subTest(filename=filename):
                source = read(f"source/modules/{filename}")
                for needle in needles:
                    self.assertIn(needle, source)
                self.assertIn("WITH_ABORT", source)
        self.assertFalse((ROOT / "source/modules/tdhf_mrsf_hessian.F90").exists())

    def test_molecule_has_single_hessian_storage_helper_with_asymmetry_metadata(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("def set_hessian_result", source)
        self.assertIn("hessian_metadata", source)
        self.assertIn("max_asymmetry", source)
        self.assertIn("0.5 * (hessian + hessian.T)", source)
        self.assertIn("def get_hess(self):", source)
        self.assertNotIn("def get_hess(self):\n        \"\"\"\n        Get hessian results\n        \"\"\"\n\n        return []", source)

    def test_saved_hessian_json_uses_inertia_not_modes_for_inertia_field(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("'inertia': self.inertia.tolist()", source)
        self.assertNotIn("'inertia': self.modes.tolist()", source)

    def test_saved_hessian_json_has_labeled_frequency_mode_vectors(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("'frequency_modes'", source)
        self.assertIn("'frequencies_cm-1'", source)
        self.assertIn("'normal_mode_eigenvectors'", source)
        self.assertIn("'normal_mode_eigenvectors_units'", source)

    def test_read_hessian_json_restores_hessian_metadata(self):
        source = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("self.hessian_metadata = data.get('hessian_metadata', {})", source)

    def test_native_cphf_exposes_reusable_static_polarizability_kernel(self):
        source = read("source/modules/cphf.F90")
        header = read("include/oqp.h")

        self.assertIn("public :: cphf_static_polarizability", source)
        self.assertIn("subroutine cphf_static_polarizability(infos, alpha)", source)
        self.assertIn("call cphf_static_polarizability(infos, alpha)", source)
        self.assertIn("bind(C, name=\"cphf_static_polarizability\")", source)
        self.assertIn("void cphf_static_polarizability(struct oqp_handle_t *inf, double *alpha);", header)
        self.assertIn("cphf_polarizability_selftest", source)

    def test_native_dipole_and_vibrational_intensity_entry_points_replace_pyscf_bridge(self):
        header = read("include/oqp.h")
        electric = read("source/modules/electric_moments.F90")
        vib = read("source/modules/vibrational_intensities.F90")
        single_point = read("pyoqp/oqp/library/single_point.py")
        external = read("pyoqp/oqp/library/external.py")

        self.assertIn("void electric_dipole_au(struct oqp_handle_t *inf, double *dipole);", header)
        self.assertIn("bind(C, name=\"electric_dipole_au\")", electric)
        self.assertIn("bind(C, name=\"vibrational_intensities_native\")", vib)
        self.assertIn("void vibrational_intensities_native(struct oqp_handle_t *inf, int64_t nmode, int64_t ncoord,", header)
        self.assertIn("oqp.vibrational_intensities_native", single_point)
        self.assertNotIn("vibrational_intensities_from_pyscf", single_point)
        self.assertNotIn("external_pyscf_finite_difference", single_point)
        self.assertNotIn("def vibrational_intensities_from_pyscf", external)

    def test_hessian_frequency_log_prints_normal_mode_eigenvectors(self):
        single_point = read("pyoqp/oqp/library/single_point.py")
        file_utils = read("pyoqp/oqp/utils/file_utils.py")

        self.assertIn("section='freq_modes'", single_point)
        self.assertIn("Normal mode eigenvectors", file_utils)
        self.assertIn("Frequencies --", file_utils)
        self.assertIn("Atom AN", file_utils)
        self.assertNotIn("for mode in block", file_utils)

    def test_cphf_solver_logs_iteration_residuals_and_timing(self):
        source = read("source/modules/cphf.F90")

        self.assertIn("use io_constants, only: iw", source)
        self.assertIn("CPHF/CPKS iterative solver", source)
        self.assertIn("INITIAL CPHF ERROR", source)
        self.assertIn("CPHF ITER", source)
        self.assertIn("CPHF wall time", source)
        self.assertIn("call cpu_time", source)

    def test_pyscf_bridge_enables_cphf_iteration_logging(self):
        source = read("pyoqp/oqp/library/external.py")

        self.assertIn("from pyscf.scf import cphf as pyscf_cphf", source)
        self.assertIn("_enable_pyscf_cphf_iteration_logging", source)
        self.assertIn("PyOQP: %s solver iterations begin", source)
        self.assertIn("kwargs.setdefault(\"verbose\", cphf_log)", source)
        self.assertIn("cphf_log.timer(f\"PyOQP: {label} solver total\"", source)
        self.assertIn("module.solve = original_solve", source)


if __name__ == "__main__":
    unittest.main()
