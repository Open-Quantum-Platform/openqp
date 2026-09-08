"""Checks for the initial TDHF/TDDFT Hessian implementation."""

import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source" / "modules" / "tdhf_hessian_components.F90"
FIXED_DENSITY = ROOT / "source" / "modules" / "tdhf_hessian_fixed_density.F90"
RESPONSE = ROOT / "source" / "modules" / "tdhf_hessian_response.F90"
Z_RHS = ROOT / "source" / "modules" / "tdhf_hessian_z_rhs.F90"
GRD2 = ROOT / "source" / "integrals" / "grd2.F90"
GRD2_RYS = ROOT / "source" / "integrals" / "grd2_rys.F90"
XC = ROOT / "source" / "modules" / "tdhf_hessian_xc.F90"
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

    def test_fixed_density_contraction_routes_xc_terms_separately(self):
        source = FIXED_DENSITY.read_text().lower()
        xc_source = XC.read_text().lower()
        self.assertRegex(source, r"infos%control%hamilton\s*>=\s*20")
        self.assertIn("xc quadrature part is assembled separately", source)
        self.assertIn("subroutine build_tdhf_xc_fixed_hessian", xc_source)
        self.assertIn("subroutine add_tdhf_xc_response_rows", xc_source)

    def test_explicit_channel_uses_one_blocked_derivative_eri_traversal(self):
        z_rhs = Z_RHS.read_text().lower()
        routine = z_rhs.split(
            "subroutine explicit_channel_derivative_matrix", 1
        )[1].split("end subroutine explicit_channel_derivative_matrix", 1)[0]
        self.assertEqual(routine.count("call grd2_operator_driver"), 1)
        self.assertNotIn("call grd2_driver", routine)
        self.assertIn("subroutine grd2_operator_driver", GRD2.read_text().lower())
        rys = GRD2_RYS.read_text().lower()
        self.assertIn("subroutine grd2_rys_compute_operator", rys)
        self.assertIn("call consumer%accumulate", rys)

    def test_blocked_channel_scatter_matches_central_difference_oracle(self):
        rng = np.random.default_rng(391)
        nbf = 6
        ids = (0, 2, 3, 5)
        i, j, k, l = ids
        coulscale = 0.83
        hfscale = 0.37
        value = -0.19

        def plus_quartet(matrix):
            return value * (
                32.0 * coulscale * matrix[i, j] * matrix[k, l]
                - 8.0
                * hfscale
                * (
                    matrix[k, i] * matrix[l, j]
                    + matrix[l, i] * matrix[k, j]
                )
            )

        def minus_quartet(matrix):
            return -2.0 * hfscale * value * (
                (matrix[i, l] - matrix[l, i])
                * (matrix[j, k] - matrix[k, j])
                + (matrix[j, l] - matrix[l, j])
                * (matrix[i, k] - matrix[k, i])
            )

        base_plus = rng.normal(size=(nbf, nbf))
        base_plus = 0.5 * (base_plus + base_plus.T)
        probe_plus = rng.normal(size=(nbf, nbf))
        probe_plus = 0.5 * (probe_plus + probe_plus.T)
        operator_plus = np.zeros((nbf, nbf))
        operator_plus[k, l] += 16.0 * coulscale * base_plus[i, j] * value
        operator_plus[i, j] += 16.0 * coulscale * base_plus[k, l] * value
        operator_plus[l, j] -= 4.0 * hfscale * base_plus[k, i] * value
        operator_plus[k, i] -= 4.0 * hfscale * base_plus[l, j] * value
        operator_plus[k, j] -= 4.0 * hfscale * base_plus[l, i] * value
        operator_plus[l, i] -= 4.0 * hfscale * base_plus[k, j] * value
        oracle_plus = 0.25 * (
            plus_quartet(base_plus + probe_plus)
            - plus_quartet(base_plus - probe_plus)
        )
        np.testing.assert_allclose(np.sum(operator_plus * probe_plus), oracle_plus)

        base_minus = rng.normal(size=(nbf, nbf))
        probe_minus = rng.normal(size=(nbf, nbf))
        operator_minus = np.zeros((nbf, nbf))

        def add_antisymmetric(row, col, coefficient):
            operator_minus[row, col] += coefficient
            operator_minus[col, row] -= coefficient

        add_antisymmetric(j, k, -hfscale * (base_minus[i, l] - base_minus[l, i]) * value)
        add_antisymmetric(i, l, -hfscale * (base_minus[j, k] - base_minus[k, j]) * value)
        add_antisymmetric(i, k, -hfscale * (base_minus[j, l] - base_minus[l, j]) * value)
        add_antisymmetric(j, l, -hfscale * (base_minus[i, k] - base_minus[k, i]) * value)
        oracle_minus = 0.25 * (
            minus_quartet(base_minus + probe_minus)
            - minus_quartet(base_minus - probe_minus)
        )
        np.testing.assert_allclose(np.sum(operator_minus * probe_minus), oracle_minus)

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
