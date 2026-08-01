"""Static safety contract for the native ROHF pair-adjoint MINRES path."""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def _subroutine(source, name):
    match = re.search(
        rf"subroutine\s+{name}\b.*?end\s+subroutine\s+{name}",
        source,
        re.IGNORECASE | re.DOTALL,
    )
    if match is None:
        raise AssertionError(f"missing Fortran subroutine {name}")
    return match.group(0).lower()


class CPHFROHFMinresContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = (ROOT / "source" / "modules" / "cphf.F90").read_text()
        cls.solve = _subroutine(cls.source, "cphf_solve_rohf")

    def test_pair_adjoint_minres_uses_an_spd_preconditioner(self):
        solve_loop = self.solve.split("do irhs = 1, nrhs", 1)[1]
        minres_branch = solve_loop.split("if (use_minres) then", 1)[1].split(
            "else", 1
        )[0]
        self.assertIn("precond=cphf_precond_rohf_minres", minres_branch)
        self.assertNotIn("precond=cphf_precond_rohf,", minres_branch)

        precond = _subroutine(self.source, "cphf_precond_rohf_minres")
        self.assertIn("y = abs(p%xminv)*x", precond)
        self.assertIn("symmetric positive definite", self.source.lower())

    def test_historical_pcg_path_retains_the_signed_preconditioner(self):
        solve_loop = self.solve.split("do irhs = 1, nrhs", 1)[1]
        pcg_branch = solve_loop.split("if (use_minres) then", 1)[1].split(
            "else", 1
        )[1]
        self.assertIn("precond=cphf_precond_rohf,", pcg_branch)
        signed = _subroutine(self.source, "cphf_precond_rohf")
        self.assertIn("y = p%xminv*x", signed)
        self.assertNotIn("abs(p%xminv)", signed)

    def test_convergence_is_certified_with_the_true_unpreconditioned_residual(self):
        self.assertRegex(
            self.solve,
            r"(?s)call\s+cphf_apbx_rohf\(ax,\s*uvec\(:,irhs\),\s*c_loc\(cgdata\)\)"
            r".*?residual_norm\s*=\s*norm2\(bvec\(:,irhs\)\s*-\s*ax\)",
        )
        solved = self.solve.split("solved =", 1)[1].split("write(iw", 1)[0]
        self.assertIn("residual_norm <= sqrt(abs(cnv))", solved)
        self.assertNotIn("minres%error", solved)

    def test_squared_residual_fails_closed_without_overflow(self):
        guard = self.solve.split(
            "residual_sq = huge(1.0_dp)", 1
        )[1].split("solved =", 1)[0]
        self.assertIn("ieee_is_finite(residual_norm)", guard)
        self.assertIn("huge(1.0_dp)/residual_norm", guard)
        self.assertIn("residual_sq = residual_norm*residual_norm", guard)
        self.assertNotIn("residual_norm**2", self.solve)

    def test_report_exposes_the_minres_exit_status(self):
        report = self.solve.split(
            'write(iw,\'(\" rohf z-vector minres rhs\"', 1
        )[1]
        report = report.split("call flush(iw)", 1)[0]
        self.assertIn('" stopped after"', report)
        self.assertIn("status=", report)
        self.assertIn("int(minres%errcode)", report)

    def test_callback_data_lifetime_covers_solver_and_true_residual(self):
        self.assertIn("type(cphf_cg_data_rohf), target :: cgdata", self.solve)
        self.assertIn("precond=cphf_precond_rohf_minres, dat=cgdata", self.solve)
        self.assertLess(
            self.solve.index("call minres%clean()"),
            self.solve.index("deallocate(famo, fbmo, xminv, fao, w2, w3, ax)"),
        )


if __name__ == "__main__":
    unittest.main()
