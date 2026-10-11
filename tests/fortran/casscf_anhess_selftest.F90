! Production context and Hessian assembly, with a test-only mode selector.
module casscf_anhess_selftest_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  use casscf_anhess_mod
  implicit none
  private
  public :: casscf_anhess_selftest
contains
  function casscf_anhess_selftest(n, nc, na, nalpha, nbeta, npar, nd, &
                                  nthreads, mode, dets, pairs, h1e, eri, ci, &
                                  hess, selected) result(status) &
      bind(C, name="casscf_anhess_selftest")
    integer(c_int32_t), value :: n, nc, na, nalpha, nbeta, npar, nthreads, mode
    integer(c_int64_t), value :: nd
    integer(c_int64_t), intent(in) :: dets(0:*)
    integer(c_int32_t), intent(in) :: pairs(0:*)
    real(c_double), intent(in) :: h1e(0:*), eri(0:*), ci(0:*)
    real(c_double), intent(inout) :: hess(0:*)
    integer(c_int32_t), intent(out) :: selected
    integer(c_int64_t) :: status, rc
    type(cas_anhess_ctx_t) :: ctx
    integer :: roots(0:0), pass
    integer(c_int64_t) :: tiny_dets(0:0)
    real(c_double) :: weights(0:0)
    real(c_double), allocatable :: first(:)

    selected = -1_c_int32_t
    roots = 0
    weights = 1.0_c_double
    tiny_dets = 0_c_int64_t
    allocate(first(0:int(npar,c_int64_t)**2-1))
    do pass = 1, 2
      ! Reinitialize one context from a smaller identity to the requested one.
      status = cas_anhess_init(ctx, 1, 0, 1, 0, 0, 1, 1, 1, &
                              1_c_int64_t, tiny_dets)
      if (status /= CAS_HESS_OK) exit
      select case (mode)
      case (-1)
        status = cas_anhess_init(ctx, int(n), int(nc), int(na), int(nalpha), &
                  int(nbeta), int(npar), 1, int(nthreads), nd, dets)
      case (0,1)
        status = cas_anhess_init(ctx, int(n), int(nc), int(na), int(nalpha), &
                  int(nbeta), int(npar), 1, int(nthreads), nd, dets, mode == 1)
      case default
        status = CAS_HESS_ERR_INPUT
      end select
      if (status /= CAS_HESS_OK) exit
      selected = merge(1_c_int32_t, 0_c_int32_t, ctx%use_mf)
      if (allocated(ctx%stack) .eqv. ctx%use_mf) then
        status = -100_c_int64_t
        exit
      end if
      status = cas_anhess_build(ctx, pairs, weights, roots, 1, dets, &
                                h1e, eri, ci, hess)
      if (status /= CAS_HESS_OK) exit
      if (pass == 1) then
        first = hess(0:int(npar,c_int64_t)**2-1)
      else
        if (any(first /= hess(0:int(npar,c_int64_t)**2-1))) then
          status = -101_c_int64_t
          exit
        end if
      end if
    end do
    ! Invalid reinitialization must dispose of the previous context too.
    rc = cas_anhess_init(ctx, 33, 0, 32, 0, 0, 1, 1, 1, &
                        1_c_int64_t, tiny_dets)
    if (rc /= CAS_HESS_ERR_INPUT .or. allocated(ctx%wmat)) status = -102_c_int64_t
    call cas_anhess_free(ctx)
    call cas_anhess_free(ctx)
  end function casscf_anhess_selftest
end module casscf_anhess_selftest_mod
