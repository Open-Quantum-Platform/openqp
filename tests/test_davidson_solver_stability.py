"""Source-level stability guards for TDHF/SF/MRSF Davidson expansions.

These tests avoid importing the compiled OpenQP runtime. They pin native Fortran
contracts that prevent Davidson residual/preconditioner code from silently
propagating NaN/Inf values or appending non-finite Krylov vectors.
"""

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TDHF_LIB_SRC = ROOT / "source" / "tdhf_lib.F90"
TDHF_SF_LIB_SRC = ROOT / "source" / "tdhf_sf_lib.F90"
TDHF_ENERGY_SRC = ROOT / "source" / "modules" / "tdhf_energy.F90"


class DavidsonSolverStabilityTests(unittest.TestCase):
    def test_rpa_residual_preconditioner_uses_finite_floor_guard(self):
        """RPA/TDA Davidson q vectors must not divide by tiny or non-finite gaps."""
        src = TDHF_LIB_SRC.read_text()

        self.assertIn("use, intrinsic :: ieee_arithmetic", src)
        self.assertIn("DAVIDSON_DENOMINATOR_FLOOR", src)
        self.assertIn("davidson_safe_denominator", src)

        solve = re.search(r"subroutine rparesvec\(.*?end subroutine rparesvec", src, re.S | re.I)
        if solve is None:
            self.fail("Could not locate rparesvec implementation")
        block = solve.group(0)

        self.assertNotIn("w_r(:,ivec)/(ee(ivec)-abd(:))", block)
        self.assertNotIn("1.0d0/(ee(ivec)-abd(:))", block)
        self.assertIn("call apply_davidson_preconditioner", block)
        self.assertIn("if (.not. ieee_is_finite(errors(ivec)))", block)
        self.assertIn("if (any(.not. ieee_is_finite(q(:,ivec,1))))", block)

        helper = re.search(
            r"subroutine apply_davidson_preconditioner\(.*?end subroutine apply_davidson_preconditioner",
            src,
            re.S | re.I,
        )
        if helper is None:
            self.fail("Missing Davidson denominator sanitizer helper")
        helper_block = helper.group(0)
        self.assertIn("davidson_safe_denominator(denom)", helper_block)
        self.assertIn("merge(DAVIDSON_DENOMINATOR_FLOOR, -DAVIDSON_DENOMINATOR_FLOOR", helper_block)
        self.assertIn("ieee_is_finite", helper_block)

    def test_rpanewb_rejects_nonfinite_candidate_vectors(self):
        """Davidson basis expansion must not normalize or append NaN/Inf vectors."""
        src = TDHF_LIB_SRC.read_text()
        solve = re.search(r"subroutine rpanewb\(.*?end subroutine rpanewb", src, re.S | re.I)
        if solve is None:
            self.fail("Could not locate rpanewb implementation")
        block = solve.group(0)

        self.assertIn("any(.not. ieee_is_finite(q(:,k)))", block)
        self.assertIn("if (.not. ieee_is_finite(fnorm))", block)
        self.assertIn("bvec(:,nvec) = q(:,k)/fnorm", block)

    def test_sf_davidson_preconditioner_guards_nonfinite_denominators(self):
        """SF/MRSF Davidson q vectors must use the same finite floor guard."""
        src = TDHF_SF_LIB_SRC.read_text()

        self.assertIn("use, intrinsic :: ieee_arithmetic", src)
        self.assertIn("SF_DAVIDSON_DENOMINATOR_FLOOR", src)
        self.assertIn("sf_davidson_safe_denominator", src)

        solve = re.search(r"subroutine sfqvec\(.*?end subroutine sfqvec", src, re.S | re.I)
        if solve is None:
            self.fail("Could not locate sfqvec implementation")
        block = solve.group(0)

        self.assertIn("sf_davidson_safe_denominator(val1)", block)
        self.assertIn("merge(SF_DAVIDSON_DENOMINATOR_FLOOR, -SF_DAVIDSON_DENOMINATOR_FLOOR", block)
        self.assertIn("ieee_is_finite(q(ii,ist))", block)
        self.assertNotIn("if( val2<1.0D-12 )then", block)

    def test_tdhf_davidson_driver_stops_before_outputs_on_nonfinite_residuals(self):
        """TDHF Davidson driver must not print/store transition data after NaN residuals."""
        src = TDHF_ENERGY_SRC.read_text()

        self.assertIn("use, intrinsic :: ieee_arithmetic", src)
        self.assertIn("TD-DFT Davidson breakdown: non-finite residual", src)

        residual_check = re.search(
            r"call rparesvec\(.*?call rpaprint",
            src,
            re.S | re.I,
        )
        if residual_check is None:
            self.fail("Could not locate TDHF Davidson residual check")
        loop_block = residual_check.group(0)
        self.assertIn("if (any(.not. ieee_is_finite(errors(imax+1:ndsr)))) then", loop_block)
        self.assertIn("ierr = 4", loop_block)
        self.assertLess(loop_block.index("if (any(.not. ieee_is_finite(errors(imax+1:ndsr)))) then"), loop_block.index("mxerr = maxval(errors(imax+1:ndsr))"))

        breakdown = src.index("case (4)")
        first_output = min(
            src.index("call get_td_transition_dipole"),
            src.index("call td_print_results"),
            src.index("call infos%dat%reserve_data(OQP_td_t"),
        )
        self.assertLess(breakdown, first_output)
        guard_block = src[breakdown:first_output]
        self.assertIn("infos%mol_energy%Davidson_converged=.false.", guard_block)
        self.assertIn("call int2_driver%clean()", guard_block)
        self.assertIn("if (dft) call dftclean(infos)", guard_block)
        self.assertIn("call measure_time", guard_block)
        self.assertRegex(guard_block, r"return\b")


if __name__ == "__main__":
    unittest.main()
