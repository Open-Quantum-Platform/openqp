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

  integer, parameter, public :: GGA_RHO_A = 1
  integer, parameter, public :: GGA_RHO_B = 2
  integer, parameter, public :: GGA_SIGMA_AA = 3
  integer, parameter, public :: GGA_SIGMA_AB = 4
  integer, parameter, public :: GGA_SIGMA_BB = 5

  type, public :: gga_tdxc_kernel_t
    real(fp) :: fr = 0.0_fp, fs = 0.0_fp
    real(fp) :: frr = 0.0_fp, frs = 0.0_fp, fss = 0.0_fp
    real(fp) :: frrr = 0.0_fp, frrs = 0.0_fp
    real(fp) :: frss = 0.0_fp, fsss = 0.0_fp
    real(fp) :: frrrr = 0.0_fp, frrrs = 0.0_fp
    real(fp) :: frrss = 0.0_fp, frsss = 0.0_fp, fssss = 0.0_fp
  end type gga_tdxc_kernel_t

  type, public :: gga_tdxc_variables_t
    real(fp) :: r = 0.0_fp, s = 0.0_fp
    real(fp) :: p = 0.0_fp, pg = 0.0_fp
    real(fp) :: v = 0.0_fp, gv = 0.0_fp, w = 0.0_fp
  end type gga_tdxc_variables_t

  public :: lda_tdxc_fixed_value
  public :: lda_tdxc_fixed_der1
  public :: lda_tdxc_fixed_der2
  public :: lda_tdxc_response_probe
  public :: lda_tdxc_polarized_kxc
  public :: lda_fourth_from_third
  public :: gga_total_kernel_from_spin
  public :: gga_fourth_from_third
  public :: gga_tdxc_fixed_value
  public :: gga_tdxc_fixed_derivatives
  public :: gga_x2_coefficients
  public :: gga_x2_coefficient_derivative
  public :: gga_x3_coefficients
  public :: gga_x3_coefficient_derivative

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

  ! Convert derivatives with respect to the Libxc spin variables
  ! (rho_a,rho_b,sigma_aa,sigma_ab,sigma_bb) to derivatives with respect to
  ! total closed-shell rho and sigma.  This is TDHXKG written as a tensor
  ! contraction: d(rho_a,rho_b)/d rho=1/2 and
  ! d(sigma_aa,sigma_ab,sigma_bb)/d sigma=1/4.
  pure subroutine gga_total_kernel_from_spin(d1, d2, d3, kernel)
    real(fp), intent(in) :: d1(5), d2(5,5), d3(5,5,5)
    type(gga_tdxc_kernel_t), intent(out) :: kernel
    real(fp), parameter :: cr(5) = [0.5_fp, 0.5_fp, 0.0_fp, 0.0_fp, 0.0_fp]
    real(fp), parameter :: cs(5) = [0.0_fp, 0.0_fp, 0.25_fp, 0.25_fp, 0.25_fp]
    integer :: i, j, l

    kernel = gga_tdxc_kernel_t()
    do i = 1, 5
      kernel%fr = kernel%fr + d1(i)*cr(i)
      kernel%fs = kernel%fs + d1(i)*cs(i)
      do j = 1, 5
        kernel%frr = kernel%frr + d2(i,j)*cr(i)*cr(j)
        kernel%frs = kernel%frs + d2(i,j)*cr(i)*cs(j)
        kernel%fss = kernel%fss + d2(i,j)*cs(i)*cs(j)
        do l = 1, 5
          kernel%frrr = kernel%frrr + d3(i,j,l)*cr(i)*cr(j)*cr(l)
          kernel%frrs = kernel%frrs + d3(i,j,l)*cr(i)*cr(j)*cs(l)
          kernel%frss = kernel%frss + d3(i,j,l)*cr(i)*cs(j)*cs(l)
          kernel%fsss = kernel%fsss + d3(i,j,l)*cs(i)*cs(j)*cs(l)
        end do
      end do
    end do
  end subroutine gga_total_kernel_from_spin

  ! GAMESS TDHXFD4: obtain the five fourth derivatives from central
  ! differences of total-density third derivatives.  The caller evaluates
  ! rho at rho*(1+-eps_r), and sigma by scaling grad(rho) with
  ! sqrt(1+-eps_s), so that sigma itself changes by (1+-eps_s).
  pure subroutine gga_fourth_from_third(rho, sigma, rho_plus, rho_minus, &
      sigma_plus, sigma_minus, eps_r, eps_s, rho_cut, kernel)
    real(fp), intent(in) :: rho, sigma, eps_r, eps_s, rho_cut
    type(gga_tdxc_kernel_t), intent(in) :: rho_plus, rho_minus
    type(gga_tdxc_kernel_t), intent(in) :: sigma_plus, sigma_minus
    type(gga_tdxc_kernel_t), intent(inout) :: kernel

    kernel%frrrr = 0.0_fp
    kernel%frrrs = 0.0_fp
    kernel%frrss = 0.0_fp
    kernel%frsss = 0.0_fp
    kernel%fssss = 0.0_fp
    if (rho < rho_cut) return
    if (eps_r > 0.0_fp .and. rho > 0.0_fp) &
      kernel%frrrr = (rho_plus%frrr-rho_minus%frrr)/(2.0_fp*eps_r*rho)
    if (sigma < rho_cut*rho_cut .or. eps_s <= 0.0_fp) return
    kernel%frrrs = (sigma_plus%frrr-sigma_minus%frrr)/(2.0_fp*eps_s*sigma)
    kernel%frrss = (sigma_plus%frrs-sigma_minus%frrs)/(2.0_fp*eps_s*sigma)
    kernel%frsss = (sigma_plus%frss-sigma_minus%frss)/(2.0_fp*eps_s*sigma)
    kernel%fssss = (sigma_plus%fsss-sigma_minus%fsss)/(2.0_fp*eps_s*sigma)
  end subroutine gga_fourth_from_third

  ! GAMESS TDHXD2G point function.  pg=grad(rho).grad(rho_P) and
  ! gv=grad(rho).grad(rho_v).
  pure elemental function gga_tdxc_fixed_value(kernel, u) result(q)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    type(gga_tdxc_variables_t), intent(in) :: u
    real(fp) :: q

    q = kernel%fr*u%p + 2.0_fp*kernel%fs*u%pg &
      + 4.0_fp*(kernel%frr*u%v*u%v &
      + 4.0_fp*kernel%frs*u%v*u%gv &
      + 4.0_fp*kernel%fss*u%gv*u%gv + 2.0_fp*kernel%fs*u%w)
  end function gga_tdxc_fixed_value

  ! GAMESS TDHXFAB: first and second derivatives of the seven-variable
  ! point function with u=(rho,sigma,rho_P,Gamma_P,rho_v,Gamma_v,|grad rho_v|2).
  pure subroutine gga_tdxc_fixed_derivatives(kernel, u, fa, fab)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    type(gga_tdxc_variables_t), intent(in) :: u
    real(fp), intent(out) :: fa(7), fab(7,7)

    fa = 0.0_fp
    fab = 0.0_fp
    fa(1) = kernel%frr*u%p + 2.0_fp*kernel%frs*u%pg &
      + 4.0_fp*kernel%frrr*u%v**2 + 16.0_fp*kernel%frrs*u%v*u%gv &
      + 16.0_fp*kernel%frss*u%gv**2 + 8.0_fp*kernel%frs*u%w
    fa(2) = kernel%frs*u%p + 2.0_fp*kernel%fss*u%pg &
      + 4.0_fp*kernel%frrs*u%v**2 + 16.0_fp*kernel%frss*u%v*u%gv &
      + 16.0_fp*kernel%fsss*u%gv**2 + 8.0_fp*kernel%fss*u%w
    fa(3) = kernel%fr
    fa(4) = 2.0_fp*kernel%fs
    fa(5) = 8.0_fp*kernel%frr*u%v + 16.0_fp*kernel%frs*u%gv
    fa(6) = 16.0_fp*kernel%frs*u%v + 32.0_fp*kernel%fss*u%gv
    fa(7) = 8.0_fp*kernel%fs

    fab(5,5) = 8.0_fp*kernel%frr
    fab(5,6) = 16.0_fp*kernel%frs
    fab(6,5) = fab(5,6)
    fab(6,6) = 32.0_fp*kernel%fss
    fab(3,1) = kernel%frr; fab(1,3) = fab(3,1)
    fab(3,2) = kernel%frs; fab(2,3) = fab(3,2)
    fab(4,1) = 2.0_fp*kernel%frs; fab(1,4) = fab(4,1)
    fab(4,2) = 2.0_fp*kernel%fss; fab(2,4) = fab(4,2)
    fab(7,1) = 8.0_fp*kernel%frs; fab(1,7) = fab(7,1)
    fab(7,2) = 8.0_fp*kernel%fss; fab(2,7) = fab(7,2)
    fab(5,1) = 8.0_fp*kernel%frrr*u%v + 16.0_fp*kernel%frrs*u%gv
    fab(1,5) = fab(5,1)
    fab(5,2) = 8.0_fp*kernel%frrs*u%v + 16.0_fp*kernel%frss*u%gv
    fab(2,5) = fab(5,2)
    fab(6,1) = 16.0_fp*kernel%frrs*u%v + 32.0_fp*kernel%frss*u%gv
    fab(1,6) = fab(6,1)
    fab(6,2) = 16.0_fp*kernel%frss*u%v + 32.0_fp*kernel%fsss*u%gv
    fab(2,6) = fab(6,2)
    fab(1,1) = kernel%frrr*u%p + 2.0_fp*kernel%frrs*u%pg &
      + 4.0_fp*kernel%frrrr*u%v**2 + 16.0_fp*kernel%frrrs*u%v*u%gv &
      + 16.0_fp*kernel%frrss*u%gv**2 + 8.0_fp*kernel%frrs*u%w
    fab(1,2) = kernel%frrs*u%p + 2.0_fp*kernel%frss*u%pg &
      + 4.0_fp*kernel%frrrs*u%v**2 + 16.0_fp*kernel%frrss*u%v*u%gv &
      + 16.0_fp*kernel%frsss*u%gv**2 + 8.0_fp*kernel%frss*u%w
    fab(2,1) = fab(1,2)
    fab(2,2) = kernel%frss*u%p + 2.0_fp*kernel%fsss*u%pg &
      + 4.0_fp*kernel%frrss*u%v**2 + 16.0_fp*kernel%frsss*u%v*u%gv &
      + 16.0_fp*kernel%fssss*u%gv**2 + 8.0_fp*kernel%fsss*u%w
  end subroutine gga_tdxc_fixed_derivatives

  pure subroutine gga_x2_coefficients(kernel, grad_r, p, grad_p, a, c, b)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    real(fp), intent(in) :: grad_r(3), p, grad_p(3)
    real(fp), intent(out) :: a, c, b(3)
    real(fp) :: gamma_p

    gamma_p = dot_product(grad_r, grad_p)
    a = kernel%frr*p + 2.0_fp*kernel%frs*gamma_p
    c = 2.0_fp*kernel%frs*p + 4.0_fp*kernel%fss*gamma_p
    b = c*grad_r + 2.0_fp*kernel%fs*grad_p
  end subroutine gga_x2_coefficients

  pure subroutine gga_x2_coefficient_derivative(kernel, grad_r, p, grad_p, &
      dr, dgrad_r, dp, dgrad_p, da, dc, db)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    real(fp), intent(in) :: grad_r(3), p, grad_p(3)
    real(fp), intent(in) :: dr, dgrad_r(3), dp, dgrad_p(3)
    real(fp), intent(out) :: da, dc, db(3)
    real(fp) :: gamma_p, dgamma_p, ds, dfrr, dfrs, dfss, dfs, c

    gamma_p = dot_product(grad_r, grad_p)
    dgamma_p = dot_product(dgrad_r, grad_p) + dot_product(grad_r, dgrad_p)
    ds = 2.0_fp*dot_product(grad_r, dgrad_r)
    dfrr = kernel%frrr*dr + kernel%frrs*ds
    dfrs = kernel%frrs*dr + kernel%frss*ds
    dfss = kernel%frss*dr + kernel%fsss*ds
    dfs = kernel%frs*dr + kernel%fss*ds
    c = 2.0_fp*kernel%frs*p + 4.0_fp*kernel%fss*gamma_p
    da = dfrr*p + kernel%frr*dp &
      + 2.0_fp*(dfrs*gamma_p + kernel%frs*dgamma_p)
    dc = 2.0_fp*(dfrs*p + kernel%frs*dp) &
      + 4.0_fp*(dfss*gamma_p + kernel%fss*dgamma_p)
    db = dc*grad_r + c*dgrad_r + 2.0_fp*dfs*grad_p &
      + 2.0_fp*kernel%fs*dgrad_p
  end subroutine gga_x2_coefficient_derivative

  pure subroutine gga_x3_coefficients(kernel, grad_r, p, grad_p, a3, c2, c3, b3)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    real(fp), intent(in) :: grad_r(3), p, grad_p(3)
    real(fp), intent(out) :: a3, c2, c3, b3(3)
    real(fp) :: gamma_p, grad_p_sq

    gamma_p = dot_product(grad_r, grad_p)
    grad_p_sq = dot_product(grad_p, grad_p)
    a3 = kernel%frrr*p*p + 4.0_fp*kernel%frrs*p*gamma_p &
      + 4.0_fp*kernel%frss*gamma_p**2 + 2.0_fp*kernel%frs*grad_p_sq
    c2 = 2.0_fp*kernel%frs*p + 4.0_fp*kernel%fss*gamma_p
    c3 = 2.0_fp*kernel%frrs*p*p + 8.0_fp*kernel%frss*p*gamma_p &
      + 8.0_fp*kernel%fsss*gamma_p**2 + 4.0_fp*kernel%fss*grad_p_sq
    b3 = c3*grad_r + 2.0_fp*c2*grad_p
  end subroutine gga_x3_coefficients

  pure subroutine gga_x3_coefficient_derivative(kernel, grad_r, p, grad_p, &
      dr, dgrad_r, dp, dgrad_p, da3, dc2, dc3, db3)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    real(fp), intent(in) :: grad_r(3), p, grad_p(3)
    real(fp), intent(in) :: dr, dgrad_r(3), dp, dgrad_p(3)
    real(fp), intent(out) :: da3, dc2, dc3, db3(3)
    real(fp) :: gamma_p, dgamma_p, grad_p_sq, dgrad_p_sq, ds
    real(fp) :: dfrrr, dfrrs, dfrss, dfsss, dfrs, dfss, c2, c3

    gamma_p = dot_product(grad_r, grad_p)
    dgamma_p = dot_product(dgrad_r, grad_p) + dot_product(grad_r, dgrad_p)
    grad_p_sq = dot_product(grad_p, grad_p)
    dgrad_p_sq = 2.0_fp*dot_product(grad_p, dgrad_p)
    ds = 2.0_fp*dot_product(grad_r, dgrad_r)
    dfrrr = kernel%frrrr*dr + kernel%frrrs*ds
    dfrrs = kernel%frrrs*dr + kernel%frrss*ds
    dfrss = kernel%frrss*dr + kernel%frsss*ds
    dfsss = kernel%frsss*dr + kernel%fssss*ds
    dfrs = kernel%frrs*dr + kernel%frss*ds
    dfss = kernel%frss*dr + kernel%fsss*ds
    c2 = 2.0_fp*kernel%frs*p + 4.0_fp*kernel%fss*gamma_p
    c3 = 2.0_fp*kernel%frrs*p*p + 8.0_fp*kernel%frss*p*gamma_p &
      + 8.0_fp*kernel%fsss*gamma_p**2 + 4.0_fp*kernel%fss*grad_p_sq
    da3 = dfrrr*p*p + 2.0_fp*kernel%frrr*p*dp &
      + 4.0_fp*(dfrrs*p*gamma_p + kernel%frrs*(dp*gamma_p+p*dgamma_p)) &
      + 4.0_fp*(dfrss*gamma_p**2 + 2.0_fp*kernel%frss*gamma_p*dgamma_p) &
      + 2.0_fp*(dfrs*grad_p_sq + kernel%frs*dgrad_p_sq)
    dc2 = 2.0_fp*(dfrs*p + kernel%frs*dp) &
      + 4.0_fp*(dfss*gamma_p + kernel%fss*dgamma_p)
    dc3 = 2.0_fp*dfrrs*p*p + 4.0_fp*kernel%frrs*p*dp &
      + 8.0_fp*(dfrss*p*gamma_p + kernel%frss*(dp*gamma_p+p*dgamma_p)) &
      + 8.0_fp*(dfsss*gamma_p**2 + 2.0_fp*kernel%fsss*gamma_p*dgamma_p) &
      + 4.0_fp*(dfss*grad_p_sq + kernel%fss*dgrad_p_sq)
    db3 = dc3*grad_r + c3*dgrad_r + 2.0_fp*dc2*grad_p + 2.0_fp*c2*dgrad_p
  end subroutine gga_x3_coefficient_derivative

end module mod_dft_gridint_tdxc_deriv
