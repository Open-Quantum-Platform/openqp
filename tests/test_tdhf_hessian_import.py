"""Checks for the initial TDHF/TDDFT Hessian implementation."""

import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source" / "modules" / "tdhf_hessian_components.F90"
FIXED_DENSITY = ROOT / "source" / "modules" / "tdhf_hessian_fixed_density.F90"
RESPONSE = ROOT / "source" / "modules" / "tdhf_hessian_response.F90"
DESIGN = ROOT / "docs" / "TDDFT_HESSIAN_IMPORT.md"


class TdhfHessianImportTests(unittest.TestCase):
    def test_final_assembly_symmetrizes_response_rows_once(self):
        source = SOURCE.read_text()
        assignment = re.search(
            r"h_total\s*=\s*h_fixed\s*\+\s*h_xc"
            r"\s*&?\s*\+\s*0\.5_dp\s*\*\s*\(\s*response_rows\s*\+\s*"
            r"transpose\(response_rows\)\s*\)",
            source,
            re.IGNORECASE | re.DOTALL,
        )
        self.assertIsNotNone(
            assignment,
            "The excitation Hessian must add 1/2(G + G^T) exactly at final assembly.",
        )

    def test_final_assembly_does_not_add_ground_state_twice(self):
        source = SOURCE.read_text().lower()
        assignment = source.split("h_total =", 1)[1]
        self.assertNotIn("h_ground", assignment)
        self.assertIn("because h_fixed already includes", source)

    def test_directional_row_asymmetry_is_recorded_before_symmetrization(self):
        source = SOURCE.read_text()
        self.assertRegex(
            source,
            r"row_asymmetry\s*=\s*maxval\(abs\(response_rows\s*-\s*"
            r"transpose\(response_rows\)\)\)",
        )

    def test_shape_mismatch_returns_an_explicit_error_status(self):
        source = SOURCE.read_text().lower()
        self.assertIn("integer, intent(out) :: status", source)
        self.assertRegex(source, r"status\s*=\s*-1\s*\n\s*return")
        self.assertRegex(source, r"status\s*=\s*0")

    def test_initial_applicability_is_fail_closed(self):
        source = SOURCE.read_text().lower()
        required = (
            "scf_type == 1",
            "td_multiplicity == 1",
            ".not. tamm_dancoff",
            ".not. needs_tau",
            "mpi_size == 1",
        )
        for condition in required:
            self.assertIn(condition, source)

    def test_fixed_density_contraction_uses_the_gradient_densities(self):
        source = FIXED_DENSITY.read_text().lower()
        for tag in ("oqp_dm_a", "oqp_wao", "oqp_td_p", "oqp_td_xpy", "oqp_td_xmy"):
            self.assertIn(tag, source)
        self.assertIn("overlap_density = overlap_density - 2.0_dp*wao", source)
        self.assertIn("one_density = dmat_a + 2.0_dp*td_p(:,1)", source)
        self.assertIn("call grd2_hess_driver", source)

    def test_fixed_density_contraction_rejects_unported_xc_terms(self):
        source = FIXED_DENSITY.read_text().lower()
        self.assertRegex(source, r"infos%control%hamilton\s*>=\s*20")
        self.assertRegex(source, r"xc quadrature\s*'//\s*&\s*'second derivatives")

    def test_design_records_the_total_energy_decomposition(self):
        design = DESIGN.read_text()
        self.assertIn("E_I(\\mathbf R) = E_0(\\mathbf R) + \\omega_I(\\mathbf R)", design)
        self.assertIn("H^{1e}_{KL} + H^{2e}_{KL} + H^{xc}_{KL}", design)
        self.assertIn("directional-row asymmetry", design)

    def test_fortran_component_assembly(self):
        compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
        if compiler is None:
            self.skipTest("GNU Fortran compiler is not available")

        with tempfile.TemporaryDirectory(prefix="oqp-tdhf-hessian-") as tmp:
            exe = Path(tmp) / "test_tdhf_hessian_components"
            subprocess.run(
                [
                    compiler,
                    str(ROOT / "source" / "precision.F90"),
                    str(SOURCE),
                    str(ROOT / "tests" / "fortran" / "test_tdhf_hessian_components.F90"),
                    "-o",
                    str(exe),
                ],
                cwd=tmp,
                check=True,
            )
            subprocess.run([str(exe)], cwd=tmp, check=True)

    def test_fortran_projected_amplitude_response(self):
        compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
        if compiler is None:
            self.skipTest("GNU Fortran compiler is not available")

        with tempfile.TemporaryDirectory(prefix="oqp-tdhf-response-") as tmp:
            exe = Path(tmp) / "test_tdhf_hessian_response"
            subprocess.run(
                [
                    compiler,
                    str(ROOT / "source" / "precision.F90"),
                    str(RESPONSE),
                    str(ROOT / "tests" / "fortran" / "test_tdhf_hessian_response.F90"),
                    "-o",
                    str(exe),
                ],
                cwd=tmp,
                check=True,
            )
            subprocess.run([str(exe)], cwd=tmp, check=True)


if __name__ == "__main__":
    unittest.main()
