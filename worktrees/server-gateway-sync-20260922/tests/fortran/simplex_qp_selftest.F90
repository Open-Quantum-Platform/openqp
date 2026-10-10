module simplex_qp_selftest_mod
  use simplex_qp, only: solve_simplex_qp
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
  implicit none

contains

  subroutine oqp_simplex_qp_solve(n, h, g, x, value, status) &
      bind(C, name='oqp_simplex_qp_solve')
    integer(c_int64_t), value :: n
    real(c_double), intent(in) :: h(n, n), g(n)
    real(c_double), intent(out) :: x(n), value
    integer(c_int), intent(out) :: status
    integer :: native_status

    call solve_simplex_qp(h, g, x, value, native_status, &
                          preferred=int(n, kind(native_status)))
    status = int(native_status, c_int)
  end subroutine oqp_simplex_qp_solve

  subroutine oqp_simplex_qp_solve_avoid(n, h, g, nforbidden, forbidden, &
                                        forbid_vertices_before, x, value, status) &
      bind(C, name='oqp_simplex_qp_solve_avoid')
    integer(c_int64_t), value :: n, nforbidden, forbid_vertices_before
    real(c_double), intent(in) :: h(n, n), g(n), forbidden(n, nforbidden)
    real(c_double), intent(out) :: x(n), value
    integer(c_int), intent(out) :: status
    integer :: native_status

    call solve_simplex_qp( &
      h, g, x, value, native_status, forbidden=forbidden, &
      nforbidden=int(nforbidden, kind(native_status)), &
      forbid_vertices_before=int(forbid_vertices_before, kind(native_status)), &
      preferred=int(n, kind(native_status)))
    status = int(native_status, c_int)
  end subroutine oqp_simplex_qp_solve_avoid

end module simplex_qp_selftest_mod
