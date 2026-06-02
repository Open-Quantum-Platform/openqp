"""Source-level stability guards for z-vector iterative solvers.

These tests intentionally avoid importing the compiled OpenQP runtime.  They pin
native Fortran safety invariants that prevent NaN/Inf propagation in the
performance-critical z-vector solver paths.
"""

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PCG_SRC = ROOT / "source" / "pcg.F90"
RHF_ZVEC_SRC = ROOT / "source" / "modules" / "tdhf_z_vector.F90"
SF_ZVEC_SRC = ROOT / "source" / "modules" / "tdhf_sf_z_vector.F90"
MRSF_ZVEC_SRC = ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90"
TDHF_LIB_SRC = ROOT / "source" / "tdhf_lib.F90"


class ZVectorSolverStabilityTests(unittest.TestCase):
    def test_pcg_init_builds_true_residual_for_nonzero_initial_guess(self):
        """PCG initialization must compute A*x0 before forming r=b-A*x0."""
        src = PCG_SRC.read_text()
        init = re.search(r"subroutine pcg_init\(this, b, update, precond, dat, x0, tol\).*?end subroutine", src, re.S | re.I)
        if init is None:
            self.fail("Could not locate pcg_init implementation")
        block = init.group(0)

        update_call = "call this%update(this%Ap, this%x, this%dat)"
        residual = "this%r(:) = this%b - this%Ap"
        self.assertIn(update_call, block)
        self.assertLess(block.index(update_call), block.index(residual))
        self.assertIn("if (any(.not. ieee_is_finite(this%Ap)))", block)

    def test_pcg_step_guards_breakdown_denominators_and_nonfinite_updates(self):
        """PCG must not divide by zero/tiny denominators or propagate NaN/Inf."""
        src = PCG_SRC.read_text()

        self.assertIn("use, intrinsic :: ieee_arithmetic", src)
        self.assertRegex(src, r"public\s+PCG_BREAKDOWN")
        self.assertRegex(src, r"PCG_BREAKDOWN\s*=\s*3")
        self.assertIn("pcg_safe_positive_denominator", src)

        step = re.search(r"subroutine pcg_step\(this\).*?end subroutine", src, re.S | re.I)
        if step is None:
            self.fail("Could not locate pcg_step implementation")
        block = step.group(0)

        self.assertNotRegex(
            block,
            r"alpha\s*=\s*rz\s*/\s*dot_product\(p,\s*Ap\)",
            "pcg_step must not divide directly by dot_product(p, Ap); guard it first.",
        )
        self.assertRegex(block, r"pap\s*=\s*dot_product\(p,\s*Ap\)")
        self.assertIn("pap = dot_product(p, Ap)", block)
        self.assertIn("if (.not. pcg_safe_positive_denominator(pap)", block)
        self.assertIn("ieee_is_finite(alpha)", block)
        self.assertIn("rz_new = dot_product(r, y)", block)
        self.assertIn("pcg_safe_positive_denominator(rz)", block)
        self.assertIn("if (.not. pcg_safe_positive_denominator(rz_new)", block)
        self.assertIn("if (.not. ieee_is_finite(beta)", block)
        self.assertIn("if (.not. ieee_is_finite(error)", block)
        self.assertIn("any(.not. ieee_is_finite(Ap))", block)
        self.assertIn("any(.not. ieee_is_finite(y)", block)
        self.assertIn("if (.not. ieee_is_finite(error)", block)
        self.assertIn("if (.not. all(ieee_is_finite(x)))", block)

    def test_pcg_optimize_only_writes_optional_outputs_when_present(self):
        """Reusable PCG driver must not assign absent optional err/cgiters outputs."""
        src = PCG_SRC.read_text()
        helper = re.search(r"subroutine pcg_optimize\(.*?end subroutine", src, re.S | re.I)
        if helper is None:
            self.fail("Could not locate pcg_optimize implementation")
        block = helper.group(0)

        self.assertIn("optional, intent(out) :: err", block)
        self.assertIn("optional, intent(out) :: cgiters", block)
        self.assertIn("if (present(err)) err = pcg%error", block)
        self.assertNotIn("\n      err = pcg%error", block)

    def test_rhf_zvector_preconditioner_clamps_near_zero_denominators(self):
        """RHF/RPA/TDA z-vector preconditioner must be finite and floor guarded."""
        src = RHF_ZVEC_SRC.read_text()

        self.assertIn("use, intrinsic :: ieee_arithmetic", src)
        self.assertIn("ZVEC_PRECOND_FLOOR", src)
        self.assertIn("sanitize_zvector_preconditioner", src)
        self.assertNotIn("xminv = 1.0d0/xm", src)
        self.assertRegex(src, r"call\s+sanitize_zvector_preconditioner\(xm,\s*xminv,\s*iw\)")

        helper = re.search(
            r"subroutine\s+sanitize_zvector_preconditioner\(xm,\s*xminv,\s*log_unit\).*?end subroutine",
            src,
            re.S | re.I,
        )
        if helper is None:
            self.fail("Missing z-vector preconditioner sanitizer helper")
        block = helper.group(0)
        self.assertRegex(block, r"abs\(denom\)\s*<\s*ZVEC_PRECOND_FLOOR")
        self.assertRegex(block, r"ieee_is_finite\(denom\)")
        self.assertRegex(block, r"1\.0_dp\s*/\s*denom")
        self.assertIn("regularized", block.lower())

    def test_sf_zvector_loop_guards_alpha_breakdown_and_nonfinite_updates(self):
        """SF z-vector loop must not divide by tiny/non-finite p^T A p or save bad vectors."""
        src = SF_ZVEC_SRC.read_text()

        self.assertIn("use, intrinsic :: ieee_arithmetic", src)
        self.assertIn("SF_ZVEC_DENOMINATOR_FLOOR", src)
        self.assertIn("sanitize_sf_zvector_preconditioner", src)
        self.assertRegex(src, r"call\s+sanitize_sf_zvector_preconditioner\(xm,\s*xminv,\s*iw\)")

        loop = re.search(r"do iter = 1, infos%control%maxit_zv.*?end do", src, re.S | re.I)
        if loop is None:
            self.fail("Could not locate SF z-vector PCG loop")
        block = loop.group(0)

        self.assertNotIn("alpha = 1.0_dp/dot_product(pk, lhs)", block)
        self.assertIn("pap = dot_product(pk, lhs)", block)
        self.assertIn("if (.not. ieee_is_finite(pap) .or. abs(pap) < SF_ZVEC_DENOMINATOR_FLOOR)", block)
        self.assertIn("alpha = 1.0_dp / pap", block)
        self.assertIn("if (.not. ieee_is_finite(alpha))", block)
        self.assertIn("if (any(.not. ieee_is_finite(lhs)) .or. any(.not. ieee_is_finite(pk)))", block)
        self.assertIn("if (any(.not. ieee_is_finite(xk)) .or. any(.not. ieee_is_finite(errv)))", block)
        self.assertIn("if (.not. ieee_is_finite(error))", block)
        self.assertIn("if (any(.not. ieee_is_finite(pk)))", block)
        self.assertIn("Z-Vector breakdown", src)

    def test_sf_zvector_breakdown_skips_density_and_w_updates_from_bad_solution(self):
        """SF z-vector must not build response densities/W with a non-finite breakdown solution."""
        src = SF_ZVEC_SRC.read_text()
        breakdown = src.index("if (zvector_breakdown) then")
        first_solution_use = min(
            src.index("call sfropcal"),
            src.index("call sfrowcal"),
        )
        self.assertLess(breakdown, first_solution_use)
        guard_block = src[breakdown:first_solution_use]
        self.assertIn("call int2_driver%clean()", guard_block)
        self.assertIn("if (dft) call dftclean(infos)", guard_block)
        self.assertIn("call measure_time", guard_block)
        self.assertRegex(guard_block, r"return\b")

    def test_rhf_zvector_pcg_breakdown_skips_density_and_w_updates_from_bad_solution(self):
        """RHF/TDA z-vector must not build P/W densities after PCG breakdown."""
        src = RHF_ZVEC_SRC.read_text()
        breakdown = src.index("case (PCG_BREAKDOWN)")
        first_solution_use = min(
            src.index("call iatogen"),
            src.index("call compute_w_mo"),
        )
        self.assertLess(breakdown, first_solution_use)
        guard_block = src[breakdown:first_solution_use]
        self.assertIn("call int2_data%clean()", guard_block)
        self.assertIn("call int2_driver%clean()", guard_block)
        self.assertIn("if (dft) call dftclean(infos)", guard_block)
        self.assertIn("call pcg%clean()", guard_block)
        self.assertIn("call measure_time", guard_block)
        self.assertRegex(guard_block, r"return\b")

    def test_mrsf_gmres_recomputes_true_residual_after_solution_update(self):
        """MRSF GMRES must not accept only the Givens residual estimate after restart updates."""
        src = MRSF_ZVEC_SRC.read_text()
        solve = re.search(r"subroutine gmres_solve\(.*?end subroutine gmres_solve", src, re.S | re.I)
        if solve is None:
            self.fail("Could not locate MRSF gmres_solve implementation")
        block = solve.group(0)

        self.assertIn("true_residual", block)
        self.assertIn("call recompute_gmres_true_residual", block)
        self.assertRegex(block, r"r\s*=\s*b\s*-\s*Ax")
        self.assertIn("true_residual = sqrt(dot_product(r, r))", block)
        self.assertIn("error_out = true_residual", block)
        self.assertIn("converged = true_residual < tol", block)

    def test_mrsf_gmres_restart_convergence_uses_true_residual_before_preconditioned_beta(self):
        """Restart-level convergence must use ||b-Ax||, not the preconditioned Arnoldi seed norm."""
        src = MRSF_ZVEC_SRC.read_text()
        solve = re.search(r"subroutine gmres_solve\(.*?end subroutine gmres_solve", src, re.S | re.I)
        if solve is None:
            self.fail("Could not locate MRSF gmres_solve implementation")
        block = solve.group(0)

        restart_setup = re.search(
            r"! Compute initial residual r = b - A\*x.*?! Apply preconditioner to residual",
            block,
            re.S | re.I,
        )
        if restart_setup is None:
            self.fail("Could not locate GMRES restart residual setup")
        setup = restart_setup.group(0)

        self.assertIn("true_residual = sqrt(dot_product(r, r))", setup)
        self.assertIn("if (.not. ieee_is_finite(true_residual))", setup)
        self.assertIn("if (true_residual < tol) then", setup)
        self.assertIn("error_out = true_residual", setup)

        pre_arnoldi = block[block.index("! Compute initial residual r = b - A*x"):block.index("! Arnoldi process")]
        self.assertLess(pre_arnoldi.index("true_residual = sqrt(dot_product(r, r))"), pre_arnoldi.index("call apply_precond"))
        self.assertNotIn("if (error < tol) then", pre_arnoldi)

    def test_mrsf_gmres_distinguishes_happy_breakdown_from_instability(self):
        """A tiny Arnoldi norm can be a happy breakdown; GMRES should update x and check true residual."""
        src = MRSF_ZVEC_SRC.read_text()
        solve = re.search(r"subroutine gmres_solve\(.*?end subroutine gmres_solve", src, re.S | re.I)
        if solve is None:
            self.fail("Could not locate MRSF gmres_solve implementation")
        block = solve.group(0)

        self.assertIn("happy_breakdown", block)
        self.assertIn("happy_breakdown = .false.", block)
        self.assertIn("GMRES: happy breakdown", block)
        self.assertIn("happy_breakdown = .true.", block)

        breakdown = block.index("happy_breakdown = .true.")
        solution_update = block.index("! Update solution: x = x + V*y")
        true_residual = block.index("call recompute_gmres_true_residual")
        self.assertLess(breakdown, solution_update)
        self.assertLess(solution_update, true_residual)

        breakdown_branch = block[block.rindex("else if", 0, breakdown):block.index("else\n          V(:,j+1)", breakdown)]
        self.assertNotIn("unstable = .true.", breakdown_branch)
        self.assertIn("H(j+1,j) = 0.0_dp", breakdown_branch)
        happy_block = block[breakdown:solution_update]
        self.assertRegex(happy_block, r"exit\b")

    def test_mrsf_gmres_reorthogonalizes_only_when_overlap_drift_exceeds_threshold(self):
        """GMRES should avoid unconditional second MGS passes but reorthogonalize drifted vectors."""
        src = MRSF_ZVEC_SRC.read_text()
        solve = re.search(r"subroutine gmres_solve\(.*?end subroutine gmres_solve", src, re.S | re.I)
        if solve is None:
            self.fail("Could not locate MRSF gmres_solve implementation")
        block = solve.group(0)

        self.assertIn("gmres_reorth_threshold", block)
        self.assertIn("max_abs_overlap", block)
        self.assertRegex(block, r"max_abs_overlap\s*=\s*max\(max_abs_overlap,\s*abs\(temp\)\)")
        self.assertIn("if (max_abs_overlap > gmres_reorth_threshold) then", block)
        reorth = block[block.index("if (max_abs_overlap > gmres_reorth_threshold) then"):]
        self.assertIn("! GMRES MGS reorthogonalization pass", reorth)
        self.assertIn("H(i,j) = H(i,j) + temp", reorth)
        self.assertIn("V(:,j+1) = V(:,j+1) - temp * V(:,i)", reorth)
        self.assertIn("if (.not. ieee_is_finite(temp) .or. .not. ieee_is_finite(H(i,j)))", reorth)

    def test_tdhf_davidson_basis_append_reorthogonalizes_when_overlap_drifts(self):
        """TDHF Davidson basis expansion should reorthogonalize vectors with residual overlap drift."""
        src = TDHF_LIB_SRC.read_text()
        helper = re.search(r"subroutine rpanewb\(.*?end subroutine rpanewb", src, re.S | re.I)
        if helper is None:
            self.fail("Could not locate TDHF Davidson rpanewb basis helper")
        block = helper.group(0)

        self.assertIn("reorth_threshold", block)
        self.assertIn("max_abs_overlap", block)
        self.assertRegex(block, r"max_abs_overlap\s*=\s*max\(max_abs_overlap,\s*abs\(bq\)\)")
        self.assertIn("if (max_abs_overlap > reorth_threshold) then", block)
        reorth = block[block.index("if (max_abs_overlap > reorth_threshold) then"):]
        self.assertIn("! Davidson MGS reorthogonalization pass", reorth)
        self.assertIn("q(:,k) = q(:,k) - bq*bvec(:,istat)", reorth)
        self.assertIn("if (.not. ieee_is_finite(bq))", reorth)
        self.assertIn("if (.not. ieee_is_finite(fnorm)) cycle", reorth)

    def test_mrsf_gmres_back_substitution_rejects_nonfinite_rhs_and_seed_solution(self):
        """GMRES triangular solve must not seed y(n) from non-finite RHS or produce NaN/Inf."""
        src = MRSF_ZVEC_SRC.read_text()
        helper = re.search(r"subroutine back_substitution\(.*?end subroutine back_substitution", src, re.S | re.I)
        if helper is None:
            self.fail("Could not locate MRSF GMRES back_substitution helper")
        block = helper.group(0)

        self.assertIn("if (.not. ieee_is_finite(b(n)))", block)
        self.assertIn("x(n) = b(n) / A(n,n)", block)
        seeded = block.index("x(n) = b(n) / A(n,n)")
        seed_check = block.index("if (.not. ieee_is_finite(x(n)))")
        self.assertLess(seeded, seed_check)
        self.assertIn("if (.not. ieee_is_finite(b(i)))", block)
        self.assertNotIn("rhs = b(i)\n        if (.not. ieee_is_finite(rhs))", block)

    def test_mrsf_operator_rejects_nonfinite_input_before_expensive_response_work(self):
        """MRSF GMRES operator applications must fail closed before integral/DFT work on NaN vectors."""
        src = MRSF_ZVEC_SRC.read_text()
        self.assertIn("ieee_value", src)
        self.assertIn("ieee_quiet_nan", src)
        operator = re.search(
            r"subroutine\s+apply_z_operator\(.*?end subroutine apply_z_operator",
            src,
            re.S | re.I,
        )
        if operator is None:
            self.fail("Could not locate MRSF z-vector operator helper")
        block = operator.group(0)
        self.assertIn("if (any(.not. ieee_is_finite(x_in))) then", block)
        self.assertIn("x_out = ieee_value(0.0_dp, ieee_quiet_nan)", block)
        self.assertIn("MRSF z-vector operator rejected non-finite input", block)
        self.assertLess(block.index("if (any(.not. ieee_is_finite(x_in)))"), block.index("allocate(int2_data)"))

    def test_mrsf_zvector_preconditioner_sanitizes_nonfinite_values_before_solver_choice(self):
        """MRSF z-vector CG/GMRES must not seed either solver with NaN/Inf preconditioners."""
        src = MRSF_ZVEC_SRC.read_text()
        self.assertIn("sanitize_mrsf_zvector_preconditioner", src)
        self.assertRegex(src, r"(?s)call\s+sfromcal\(xm,\s*xminv,.*?call\s+sanitize_mrsf_zvector_preconditioner\(xm,\s*xminv,\s*iw\)")

        helper = re.search(
            r"subroutine\s+sanitize_mrsf_zvector_preconditioner\(xm,\s*xminv,\s*log_unit\).*?end subroutine",
            src,
            re.S | re.I,
        )
        if helper is None:
            self.fail("Missing MRSF z-vector preconditioner sanitizer helper")
        block = helper.group(0)
        self.assertIn("if (.not. ieee_is_finite(denom) .or. abs(denom) < MRSF_ZVEC_DENOMINATOR_FLOOR)", block)
        self.assertIn("xminv(i) = 1.0_dp / denom", block)
        self.assertIn("if (regularized > 0) then", block)

    def test_mrsf_default_cg_guards_pap_denominator_and_breakdown_solution(self):
        """Default fast MRSF CG path must not divide by tiny p^T A p or use bad xk."""
        src = MRSF_ZVEC_SRC.read_text()
        self.assertIn("MRSF_ZVEC_DENOMINATOR_FLOOR", src)
        self.assertIn("mrsf_zvector_breakdown", src)

        loop = re.search(r"do iter = 1, infos%control%maxit_zv.*?end do", src, re.S | re.I)
        if loop is None:
            self.fail("Could not locate MRSF default CG z-vector loop")
        block = loop.group(0)

        self.assertNotIn("alpha = 1.0_dp/dot_product(pk, lhs)", block)
        self.assertIn("pap = dot_product(pk, lhs)", block)
        self.assertIn("if (.not. ieee_is_finite(pap) .or. abs(pap) < MRSF_ZVEC_DENOMINATOR_FLOOR)", block)
        self.assertIn("alpha = 1.0_dp / pap", block)
        self.assertIn("if (.not. ieee_is_finite(alpha))", block)
        self.assertIn("if (any(.not. ieee_is_finite(lhs)) .or. any(.not. ieee_is_finite(pk)))", block)
        self.assertIn("if (any(.not. ieee_is_finite(xk)) .or. any(.not. ieee_is_finite(errv)))", block)
        self.assertIn("if (.not. ieee_is_finite(error))", block)

        breakdown = src.index("if (mrsf_zvector_breakdown) then")
        first_solution_use = min(src.index("call sfropcal"), src.index("call mrsfqropcal"))
        self.assertLess(breakdown, first_solution_use)
        guard_block = src[breakdown:first_solution_use]
        self.assertIn("call int2_data%clean()", guard_block)
        self.assertIn("call int2_driver%clean()", guard_block)
        self.assertIn("if (dft) call dftclean(infos)", guard_block)
        self.assertRegex(guard_block, r"return\b")


if __name__ == "__main__":
    unittest.main()
