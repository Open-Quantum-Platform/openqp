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


class ZVectorSolverStabilityTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
