! Self-test ABI for the TDHF/TDDFT Davidson residual step (rparesvec) and the
! convergence measure the driver takes from it (rpa_max_residual).
module rpa_residual_selftest_mod
  use tdhf_lib, only: rparesvec, rpa_max_residual
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t, c_bool
  implicit none

  public :: oqp_rpa_residual_selftest

contains

!> On entry w_l/w_r hold the (A+B)V and (A-B)V products of the Ritz vectors
!> v_l/v_r; on exit they hold the residuals. errors are the squared residual
!> norms, q(:,:,1:2) the left/right preconditioned corrections, mxerr the
!> driver's convergence measure and converged = (mxerr <= tol).
  subroutine oqp_rpa_residual_selftest(n, ndsr, imax, tamm_dancoff, tol, &
      w_l, w_r, v_l, v_r, ee, abd, errors, q, mxerr, converged) &
      bind(C, name='oqp_rpa_residual_selftest')
    integer(c_int64_t), value :: n, ndsr, imax
    logical(c_bool), value :: tamm_dancoff
    real(c_double), value :: tol
    real(c_double), intent(inout) :: w_l(n, ndsr), w_r(n, ndsr)
    real(c_double), intent(inout) :: v_l(n, ndsr), v_r(n, ndsr)
    real(c_double), intent(in) :: ee(ndsr), abd(n)
    real(c_double), intent(out) :: errors(ndsr), q(n, ndsr, 2), mxerr
    integer(c_int), intent(out) :: converged

    integer :: nd, im

    nd = int(ndsr)
    im = int(imax)
    q = 0.0_c_double
    call rparesvec(q, w_l, w_r, v_l, v_r, ee, abd, nd, errors, tol, im, &
                   logical(tamm_dancoff))
    mxerr = rpa_max_residual(errors, im, nd)
    converged = merge(1_c_int, 0_c_int, mxerr <= tol)
  end subroutine oqp_rpa_residual_selftest

end module rpa_residual_selftest_mod
