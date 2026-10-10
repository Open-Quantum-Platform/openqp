import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class TestMP2GradientReviewGuards(unittest.TestCase):
    def test_mp2_opens_log_before_cphf_and_aborts_on_nonconvergence(self):
        source = (ROOT / "source/modules/mp2_gradient.F90").read_text()

        open_log = "open(unit=iw, file=infos%log_filename, position='append')"
        solve = "call cphf_solve(infos, 1, rhs, u"
        self.assertEqual(source.count(open_log), 1)
        self.assertLess(source.index(open_log), source.index(solve))
        self.assertIn("converged=cphf_converged", source)
        self.assertIn("if (.not. cphf_converged) then", source)
        self.assertIn("MP2 analytic gradient CPHF response did not converge", source)

    def test_closed_shell_cphf_exposes_all_rhs_convergence_status(self):
        source = (ROOT / "source/modules/cphf.F90").read_text()

        self.assertIn(
            "subroutine cphf_solve(infos, nrhs, bvec, uvec, tol, maxit, converged)",
            source,
        )
        self.assertIn("rhs_converged = pcg%errcode == PCG_CONVERGED", source)
        self.assertIn("all_converged = all_converged .and. rhs_converged", source)
        self.assertIn("if (present(converged)) converged = all_converged", source)


if __name__ == "__main__":
    unittest.main()
