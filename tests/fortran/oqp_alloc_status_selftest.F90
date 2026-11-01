!> Self-test of the tagarray create status -> oqp_alloc return code mapping.
!>
!> A real tagarray allocation failure cannot be induced portably: GCC 15
!> libstdc++ on macOS traps (SIGILL) in a failing nothrow aligned operator
!> new[] instead of returning null.  This harness therefore feeds the named
!> status constants of the pinned tagarray_defines module through the mapping
!> that oqp_alloc itself uses, so a renamed or renumbered tagarray status, or
!> an allocation failure mapped to the invalid-input code, fails the test.
module oqp_alloc_status_selftest_m
  use iso_c_binding, only: c_int32_t, c_int64_t
  use tagarray_defines, only: TA_OK, TA_MEMORY_ALLOCATION_ERROR, &
      TA_DATA_INCORRECT_DIMENSIONS, TA_INCORRECT_TAG, &
      TA_CONTAINER_RECORD_EXISTS, TA_DATA_TOO_MANY_DIMENSIONS, &
      TA_UNKNOWN_ERROR, TA_NOT_IMPLEMENTED
  use c_interop, only: oqp_alloc_create_code
  implicit none
  private
  public :: oqp_alloc_status_selftest
contains
  !> Fill statuses(1:n) with tagarray create statuses in a fixed order and
  !> codes(1:n) with the oqp_alloc return code for each.  Order:
  !> OK, MEMORY_ALLOCATION_ERROR, DATA_INCORRECT_DIMENSIONS, INCORRECT_TAG,
  !> CONTAINER_RECORD_EXISTS, DATA_TOO_MANY_DIMENSIONS, UNKNOWN_ERROR,
  !> NOT_IMPLEMENTED.  n = -1 when nmax is too small.
  subroutine oqp_alloc_status_selftest(nmax, statuses, codes, n) &
      bind(C, name='oqp_alloc_status_selftest')
    integer(c_int64_t), value, intent(in) :: nmax
    integer(c_int32_t), intent(out) :: statuses(nmax)
    integer(c_int64_t), intent(out) :: codes(nmax)
    integer(c_int64_t), intent(out) :: n
    integer(c_int32_t), parameter :: cases(8) = [ &
        TA_OK, TA_MEMORY_ALLOCATION_ERROR, TA_DATA_INCORRECT_DIMENSIONS, &
        TA_INCORRECT_TAG, TA_CONTAINER_RECORD_EXISTS, &
        TA_DATA_TOO_MANY_DIMENSIONS, TA_UNKNOWN_ERROR, TA_NOT_IMPLEMENTED]
    integer(c_int64_t) :: i

    n = -1
    if (nmax < size(cases, kind=c_int64_t)) return
    n = size(cases, kind=c_int64_t)
    do i = 1, n
      statuses(i) = cases(i)
      codes(i) = oqp_alloc_create_code(cases(i))
    end do
  end subroutine oqp_alloc_status_selftest
end module oqp_alloc_status_selftest_m
