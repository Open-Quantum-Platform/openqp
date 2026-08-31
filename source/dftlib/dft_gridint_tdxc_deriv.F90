! Restricted LDA XC derivatives used by the analytic TDDFT Hessian.
!
! All densities in this module are total (alpha+beta) densities.  This is the
! convention used by the restricted OpenQP xc_der1/xc_der2_contr/
! xc_der3_contr path; the spin-slot factors in the GAMESS TDHXK3 routines must
! therefore not be applied a second time.
module mod_dft_gridint_tdxc_deriv

  use precision, only: fp

  implicit none
  private

  public :: lda_tdxc_fixed_value
  public :: lda_tdxc_fixed_der1
  public :: lda_tdxc_fixed_der2
  public :: lda_tdxc_response_probe
  public :: lda_tdxc_polarized_kxc
  public :: lda_fourth_from_third

contains

  ! v = d e_xc/d rho, f = d2 e_xc/d rho2.  The factor four is the
  ! closed-shell singlet transition-density factor of GAMESS EGRID(2).
  pure elemental function lda_tdxc_fixed_value(v, f, p, x) result(q)
    real(fp), intent(in) :: v, f, p, x
    real(fp) :: q
    q = v*p + 4.0_fp*f*x*x
  end function lda_tdxc_fixed_value

  ! First nuclear derivative of lda_tdxc_fixed_value.  r is the ground-state
  ! density, p the effective difference density, and x the transition density.
  pure elemental function lda_tdxc_fixed_der1(v, f, k, p, x, &
                                               r_a, p_a, x_a) result(q_a)
    real(fp), intent(in) :: v, f, k, p, x
    real(fp), intent(in) :: r_a, p_a, x_a
    real(fp) :: q_a
    q_a = f*r_a*p + v*p_a &
        + 4.0_fp*(k*r_a*x*x + 2.0_fp*f*x*x_a)
  end function lda_tdxc_fixed_der1

  ! Exact mixed second derivative.  l is d4 e_xc/d rho4.  The nuclear
  ! derivatives include AO translation, moving-grid weights outside this
  ! point function, and (where requested) density response.
  pure elemental function lda_tdxc_fixed_der2(v, f, k, l, p, x, &
      r_a, r_b, r_ab, p_a, p_b, p_ab, x_a, x_b, x_ab) result(q_ab)
    real(fp), intent(in) :: v, f, k, l, p, x
    real(fp), intent(in) :: r_a, r_b, r_ab
    real(fp), intent(in) :: p_a, p_b, p_ab
    real(fp), intent(in) :: x_a, x_b, x_ab
    real(fp) :: q_ab

    q_ab = (k*r_a*r_b + f*r_ab)*p &
         + f*(r_a*p_b + r_b*p_a) + v*p_ab
    q_ab = q_ab + 4.0_fp*( &
        (l*r_a*r_b + k*r_ab)*x*x &
      + 2.0_fp*k*(r_a*x*x_b + r_b*x*x_a) &
      + 2.0_fp*f*(x_a*x_b + x*x_ab))
  end function lda_tdxc_fixed_der2

  ! Density-response probe whose skeleton derivative is the XC response row:
  ! f rho(dD) rho(Peff) + 4 k rho(dD) rho(X)^2
  !                         + 8 f rho(X) rho(dX).
  pure elemental function lda_tdxc_response_probe(f, k, p, x, dp, dx) result(q)
    real(fp), intent(in) :: f, k, p, x, dp, dx
    real(fp) :: q
    q = f*dp*p + 4.0_fp*k*dp*x*x + 8.0_fp*f*x*dx
  end function lda_tdxc_response_probe

  ! Bilinear third-kernel action by polarization.  This avoids assumptions
  ! about packed spin components in a caller which already has quadratic Gxc.
  pure elemental function lda_tdxc_polarized_kxc(g_apb, g_a, g_b) result(g_ab)
    real(fp), intent(in) :: g_apb, g_a, g_b
    real(fp) :: g_ab
    g_ab = 0.5_fp*(g_apb - g_a - g_b)
  end function lda_tdxc_polarized_kxc

  ! Symmetric high-accuracy derivative of k(rho)=d3 e_xc/d rho3.  Supplying
  ! k at rho +/- h and rho +/- 2h keeps Libxc itself as the source of truth.
  pure elemental function lda_fourth_from_third(k_m2, k_m1, k_p1, k_p2, h) &
      result(l)
    real(fp), intent(in) :: k_m2, k_m1, k_p1, k_p2, h
    real(fp) :: l
    if (h <= 0.0_fp) then
      l = 0.0_fp
    else
      l = (k_m2 - 8.0_fp*k_m1 + 8.0_fp*k_p1 - k_p2)/(12.0_fp*h)
    end if
  end function lda_fourth_from_third

end module mod_dft_gridint_tdxc_deriv
