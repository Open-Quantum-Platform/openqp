"""Source-level stability guards for MRSF/SF GMRES and triangular back-solve."""

import re
import unittest
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "source" / "modules" / "tdhf_mrsf_z_vector.F90"


class TestMrsfGmresStability(unittest.TestCase):
    def setUp(self):
        self.src = SRC.read_text()

    def test_gmres_solver_checks_denominator_and_nonfinite_inputs(self):
        self.assertIn("use, intrinsic :: ieee_arithmetic", self.src)
        self.assertIn("GMRES_DENOMINATOR_FLOOR", self.src)

        gmres = re.search(r"subroutine gmres_solve\(.*?end subroutine gmres_solve", self.src, re.S | re.I)
        if gmres is None:
            self.fail("Could not locate gmres_solve")
        block = gmres.group(0)

        self.assertIn("if (unstable) then", block)
        self.assertIn("if (.not. ieee_is_finite(beta)", block)
        self.assertIn("if (.not. ieee_is_finite(H(j+1,j))", block)
        self.assertIn("call back_substitution", block)
        self.assertIn("error_out = huge(1.0_dp)", block)

        # New stability hardening checks
        self.assertIn("any(.not. ieee_is_finite(w))", block)
        self.assertIn("any(.not. ieee_is_finite(V(:,j)))", block)
        self.assertIn("any(.not. ieee_is_finite(V(:,j+1)))", block)
        self.assertIn("if (.not. ieee_is_finite(temp) .or. .not. ieee_is_finite(H(i,j)))", block)
        self.assertIn("if (.not. ieee_is_finite(H(i,j)))", block)
        self.assertIn("if (.not. ieee_is_finite(error_initial))", block)

    def test_gmres_restarts_honor_safe_dimension_bounds(self):
        gmres = re.search(r"subroutine gmres_solve\(.*?end subroutine gmres_solve", self.src, re.S | re.I)
        if gmres is None:
            self.fail("Could not locate gmres_solve")
        block = gmres.group(0)

        # Restart/basis-size guards prevent degenerate allocations and divide-by-zero in
        # residual-ratio reporting.
        self.assertIn("if (n <= 0 .or. restart <= 0) then", block)
        self.assertIn("if (error_initial == 0.0_dp)", block)
        self.assertIn("if (unstable .or. converged) then", block)
        self.assertIn("Relative reduction : not available", block)

    def test_givens_rotation_uses_overflow_safe_hypot_scaling(self):
        gmres = re.search(r"subroutine gmres_solve\(.*?end subroutine gmres_solve", self.src, re.S | re.I)
        if gmres is None:
            self.fail("Could not locate gmres_solve")
        block = gmres.group(0)

        # Direct sqrt(a*a + b*b) can overflow before normalization even when
        # the final cosine/sine are representable. GMRES must scale first.
        self.assertIn("scale = max(abs(a), abs(b))", block)
        self.assertIn("r = scale * sqrt((a/scale)**2 + (b/scale)**2)", block)
        self.assertNotIn("sqrt(a*a + b*b)", block)

    def test_back_substitution_avoids_zero_or_nonfinite_pivots(self):
        bs = re.search(r"subroutine back_substitution\(A, b, x, n, unstable\).*?end subroutine back_substitution", self.src, re.S | re.I)
        if bs is None:
            self.fail("Could not locate back_substitution")
        block = bs.group(0)

        self.assertIn("logical, intent(out) :: unstable", block)
        self.assertIn("if (n <= 0) then", block)
        self.assertIn("abs(A(n,n)) < GMRES_DENOMINATOR_FLOOR", block)
        self.assertIn("if (.not. ieee_is_finite(A(i,j))", block)
        self.assertIn("if (.not. ieee_is_finite(rhs)", block)
        self.assertIn("if (.not. ieee_is_finite(x(i)))", block)
        self.assertIn("unstable = .true.", block)


if __name__ == "__main__":
    unittest.main()
