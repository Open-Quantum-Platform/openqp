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
        cls.apbx = _subroutine(cls.source, "cphf_apbx_rohf")
        cls.apbx_batch = _subroutine(cls.source, "cphf_apbx_rohf_batch")
        minres_source = (ROOT / "source" / "minres.F90").read_text()
        minres_init = re.search(
            r"subroutine\s+minres_init\b.*?end\s+subroutine",
            minres_source,
            re.IGNORECASE | re.DOTALL,
        )
        if minres_init is None:
            raise AssertionError("missing Fortran subroutine minres_init")
        cls.minres_init = minres_init.group(0).lower()

    def test_pair_adjoint_minres_uses_an_spd_preconditioner(self):
        minres_branch = self.solve.split(
            "! keep one scalar paige-saunders recurrence", 1
        )[1].split("\n    else", 1)[0]
        self.assertIn("precond=cphf_precond_rohf_minres", minres_branch)
        self.assertNotIn("precond=cphf_precond_rohf,", minres_branch)

        precond = _subroutine(self.source, "cphf_precond_rohf_minres")
        self.assertIn("y = abs(p%xminv)*x", precond)
        self.assertIn("symmetric positive definite", self.source.lower())

    def test_historical_pcg_path_retains_the_signed_preconditioner(self):
        pcg_branch = self.solve.split(
            "! keep one scalar paige-saunders recurrence", 1
        )[1].split("\n    else", 1)[1]
        self.assertIn("precond=cphf_precond_rohf,", pcg_branch)
        signed = _subroutine(self.source, "cphf_precond_rohf")
        self.assertIn("y = p%xminv*x", signed)
        self.assertNotIn("abs(p%xminv)", signed)

    def test_convergence_is_certified_with_the_true_unpreconditioned_residual(self):
        self.assertRegex(
            self.solve,
            r"(?s)do\s+irhs\s*=\s*1,\s*nrhs.*?"
            r"call\s+cphf_apbx_rohf\(ax,\s*uvec\(:,irhs\),\s*c_loc\(cgdata\)\)"
            r".*?residual_norm\s*=\s*norm2\(bvec\(:,irhs\)\s*-\s*ax\)",
        )
        certification = self.solve.split(
            "! krylov iterations deliberately use", 1
        )[1].split("else", 1)[0]
        self.assertNotIn("cphf_apbx_rohf_batch", certification)
        self.assertIn("legacy scalar operator", certification)
        solved = self.solve.split("solved =", 1)[1].split("write(iw", 1)[0]
        self.assertIn("residual_norm <= sqrt(abs(rhs_cnv(irhs)))", solved)
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
        self.assertIn("int(minres_batch(irhs)%errcode)", report)

    def test_callback_data_lifetime_covers_solver_and_true_residual(self):
        self.assertIn("type(cphf_cg_data_rohf), target :: cgdata", self.solve)
        self.assertIn("precond=cphf_precond_rohf_minres, dat=cgdata", self.solve)
        scalar_cert = self.solve.index(
            "call cphf_apbx_rohf(ax, uvec(:,irhs), c_loc(cgdata))"
        )
        self.assertLess(scalar_cert, self.solve.index("if (dft) call dftclean(infos)"))
        self.assertLess(
            scalar_cert,
            self.solve.index(
                "deallocate(famo, fbmo, xminv, fao, w2, w3, ax, rhs_cnv)"
            ),
        )
        self.assertLess(
            self.solve.index("call minres_batch(irhs)%clean()"),
            self.solve.index(
                "deallocate(famo, fbmo, xminv, fao, w2, w3, ax, rhs_cnv)"
            ),
        )

    def test_independent_minres_recursions_share_batched_physics_kernels(self):
        self.assertIn("call minres_batch(irhs)%prepare_step(ready)", self.solve)
        self.assertIn("call minres_batch(irhs)%finish_step()", self.solve)
        self.assertIn(
            "call cphf_apbx_rohf_batch(ax_batch(:,1:nactive)", self.solve
        )
        self.assertIn("call p%int2_driver%run(", self.apbx_batch)
        self.assertIn("call int2_driver_batch%set_screening()", self.solve)
        self.assertIn("nmtx=nvec", self.apbx_batch)
        self.assertIn("call utddft_fxc(", self.apbx_batch)
        self.assertIn("batching is conservative", self.apbx_batch)
        self.assertIn("active width may safely shrink from 3 to 2 to 1", self.apbx_batch)
        self.assertNotIn("allocate(", self.apbx_batch)

    def test_unrestricted_eri_consumer_accepts_adjacent_spin_batches(self):
        int2 = (ROOT / "source" / "integrals" / "int2.F90").read_text().lower()
        start = int2.split("subroutine int2_urohf_data_t_parallel_start", 1)[1]
        start = start.split("end subroutine", 1)[0]
        update = _subroutine(int2, "int2_urohf_data_t_update")
        self.assertIn("this%nfocks = size(this%d,2)", start)
        self.assertIn("mod(this%nfocks,2)", start)
        self.assertIn("show_message", start)
        self.assertIn("do ifock = 1, this%nfocks, 2", update)
        self.assertIn("ifock:ifock+1", update)

    def test_xc_consumer_reallocates_for_each_active_width(self):
        fxc = (
            ROOT / "source" / "dftlib" / "dft_gridint_fxc.F90"
        ).read_text().lower()
        parallel = fxc.split("subroutine parallel_start(self, xce, nthreads)", 1)[1]
        parallel = parallel.split("end subroutine", 1)[0]
        utd = fxc.split("subroutine utddft_fxc(", 1)[1]
        utd = utd.split("end subroutine", 1)[0]
        self.assertIn("call self%clean()", parallel)
        self.assertIn("self%nmtx", parallel)
        self.assertIn("type(xc_consumer_tde_t) :: dat", utd)
        self.assertIn("dat%nmtx = nmtx", utd)
        self.assertIn("call dat%clean()", utd)

    def test_rohf_operator_reuses_solver_owned_workspace(self):
        """Krylov Hessian actions must not allocate their eleven work arrays."""
        self.assertNotIn("allocate(", self.apbx)
        for name in (
            "xa_work",
            "xb_work",
            "x2a_work",
            "x2b_work",
            "dm_work",
            "v_work",
            "kmat_work",
            "dm_tri_work",
            "pfock_work",
        ):
            self.assertIn(f"cgdata%{name} =>", self.solve)
            self.assertIn(f"p%{name}", self.apbx)
        # Persistent packed Fock storage must retain the zero initialization
        # previously supplied by ALLOCATE(..., SOURCE=0).
        self.assertIn("pfock = 0.0_dp", self.apbx)

    def test_implicit_zero_guess_skips_the_initial_operator_action(self):
        residual_setup = self.minres_init.split("! r1 = b - a x", 1)[1].split(
            "! y = m^-1 r1", 1
        )[0]
        explicit_guess, zero_guess = residual_setup.rsplit("else", 1)
        self.assertIn("if (present(x0)) then", explicit_guess)
        self.assertIn("if (present(ax0)) then", explicit_guess)
        self.assertIn("this%av = ax0", explicit_guess)
        self.assertIn("call this%update(this%av, this%x, this%dat)", explicit_guess)
        self.assertIn("all(ieee_is_finite(this%av))", explicit_guess)
        self.assertIn("this%r1 = this%b - this%av", explicit_guess)
        self.assertIn("this%r1 = this%b", zero_guess)
        self.assertNotIn("this%update", zero_guess)

    def test_rohf_response_fock_transforms_only_the_required_vo_blocks(self):
        response = self.apbx.split("call get_response_packed", 1)[1].split(
            "call rohf_pack_trial", 1
        )[0]
        # V Co, followed by Cv^T (V Co), for alpha and beta respectively.
        for nocc, nvir, out in (
            ("nocca", "nvira", "x2a"),
            ("noccb", "nvirb", "x2b"),
        ):
            self.assertRegex(
                response,
                rf"(?s)call\s+dgemm\('n','n',\s*nbf,\s*{nocc},\s*nbf,"
                rf".*?p%mo,\s*nbf,\s*0\.0_dp,\s*work2,\s*nbf\)",
            )
            self.assertRegex(
                response,
                rf"(?s)call\s+dgemm\('t','n',\s*{nvir},\s*{nocc},\s*nbf,"
                rf".*?work2,\s*nbf,\s*1\.0_dp,\s*{out},\s*{nvir}\)",
            )
        self.assertNotRegex(
            response, r"call\s+dgemm\([^\n]*nbf,\s*nbf,\s*nbf"
        )
        self.assertNotIn("work3(nocca+1:nbf", response)
        self.assertNotIn("work3(noccb+1:nbf", response)


if __name__ == "__main__":
    unittest.main()
