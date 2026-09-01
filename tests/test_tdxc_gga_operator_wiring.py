"""Focused contracts for restricted-GGA response-operator derivatives."""

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GXC = (ROOT / "source" / "dftlib" / "dft_gridint_gxc.F90").read_text().lower()
RHS = (ROOT / "source" / "modules" / "tdhf_hessian_rhs.F90").read_text().lower()
TDHESS = (ROOT / "source" / "modules" / "tdhf_hessian.F90").read_text().lower()
Z_RHS = (ROOT / "source" / "modules" / "tdhf_hessian_z_rhs.F90").read_text().lower()
TDXC_GRAD = (ROOT / "source" / "dftlib" / "dft_gridint_tdxc_grad.F90").read_text().lower()
TDHESS_XC = (ROOT / "source" / "modules" / "tdhf_hessian_xc.F90").read_text().lower()
TDHF_GRAD = (ROOT / "source" / "modules" / "tdhf_gradient.F90").read_text().lower()


def test_restricted_gga_x3_uses_total_density_algebra():
    assert "use mod_dft_gridint_tdxc_deriv" in GXC
    assert "call restricted_gga_total_kernel(xce, i, kernel)" in GXC
    assert "call gga_x3_coefficients(kernel, grad_r_total" in GXC
    assert "grad_r_total = 2.0_fp*drho(1:3,i)" in GXC
    assert "grad_p_total = 2.0_fp*drrho(1:3,1,i,j)" in GXC
    assert "2.0_fp*rrho(1,i,j)" in GXC


def test_tdhxkg_probe_directions_have_closed_shell_factors():
    assert "cr(2) = [0.5_fp, 0.5_fp]" in GXC
    assert "cs(3) = [0.25_fp, 0.25_fp, 0.25_fp]" in GXC
    for field in ("fr", "fs", "frr", "frs", "fss", "frrr", "frrs", "frss", "fsss"):
        assert re.search(rf"kernel%{field}\s*=", GXC)


def test_lda_operation_order_is_preserved_in_its_own_branch():
    lda = GXC.split("if (xce%funtyp == oqp_funtyp_lda) then", 1)[1].split("else if", 1)[0]
    assert "tmp(:,i,1) = 0.5_fp*g_r(1)*aov(:,i)" in lda
    assert "gga_x3_coefficients" not in lda


def test_response_paths_reach_x2_skeleton_and_x3_response_terms():
    # KORD=2 skeleton derivative: polarized derivative of the existing analytic
    # TDDFT XC gradient.  KORD=3 response half: the GGA-wired Gxc consumer.
    assert "call tddft_xc_gradient" in Z_RHS
    assert "call tddft_gxc" in RHS
    assert "deri_full_p(:,:,kk)=deri_full_p(:,:,kk)+fx(:,:,3*kk-2)" in RHS


def test_kxc_ground_response_crosses_the_one_spin_grid_boundary():
    assert "pvs+0.5_dp*dground(:,:,kk)" in RHS
    assert "dx(:,:,3*kk)=0.5_dp*dground(:,:,kk)" in RHS


def test_zero_orbital_connection_shortcut_is_never_applied_to_dft():
    assert "zero_orbital_connection=infos%control%hamilton<20" in TDHESS
    assert "zero_orbital_connection = infos%control%hamilton < 20" in RHS


def test_gxc_keeps_cross_sigma_and_self_sigma_in_distinct_arrays():
    # Two unrestricted-spin passes use vector assignments; the restricted
    # pass uses scalar sums.  None may overwrite sigma with the self product.
    assert GXC.count("ssigma = [2*dsaa, 2*dsbb, (dsab+dsba)]") == 2
    assert "ssigma = 2*sum(drrho(1:3,1,i,j)*drrho(1:3,1,i,j))" in GXC
    assert sum(line.strip() == "sigma = [2*dsaa, 2*dsbb, (dsab+dsba)]"
               for line in GXC.splitlines()) == 2


def test_restricted_xc_gradient_accumulates_into_the_total_gradient():
    body = TDXC_GRAD.split("subroutine tddft_xc_gradient", 1)[1].split(
        "end subroutine", 1
    )[0]
    assert "intent(inout) :: dedft" in body
    assert "dedft = 0.0_fp" not in body
    # Hessian XC terms now use direct grid consumers, not a molecular-gradient
    # finite-difference helper that could accidentally reset the total.
    assert "subroutine xc_gradient" not in TDHESS_XC
    assert "tddft_lda_fixed_hessian" in TDHESS_XC
    assert "tddft_gga_fixed_hessian" in TDHESS_XC


def test_production_tddft_gradient_preserves_the_validated_fixed_grid_path():
    # The standard restricted TDDFT gradient already matches an independent
    # PySCF reference.  Moving-grid probes are used by the MRSF correction and
    # Hessian response machinery, but must not be added to this production call.
    body = TDHF_GRAD.split("subroutine tdhf_gradient(infos)", 1)[1].split(
        "end subroutine", 1
    )[0]
    assert "include_weight_derivative" not in body
    assert "include_ground_state" not in body
    assert "call derexc_blk" not in body
    assert "threshold=1.0d-14" in body
    assert "dedft = dedft + dat%nucgrad(:,:,1)" in TDXC_GRAD


def test_restricted_owner_motion_alone_gets_closed_shell_factor():
    assert "merge(1.0_fp,2.0_fp,xce%hasbeta)*sum(tmpgrad, dim=1)" in TDXC_GRAD
    assert "2.0_fp*dat%nucgrad(:,:,1)" not in TDXC_GRAD
