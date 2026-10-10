!> @brief Global maximum-overlap correspondence for orbital/state tracking.
!>
!> The assignment maximises sum_i |overlap(i, assignment(i))| subject to a
!> one-to-one row/column correspondence.  This avoids the order dependence of
!> row-wise greedy matching at orbital and state crossings.
module state_tracking_mod

  use iso_c_binding, only: c_double, c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision, only: dp

  implicit none
  private

  public :: maximum_overlap_assignment
  public :: oqp_maximum_overlap_assignment_C
  public :: oqp_diagonal_phase_tracking_C

contains

  !> Solve a square maximum-weight bipartite assignment.
  !>
  !> @param[in]  overlap     row/column overlap matrix
  !> @param[out] assignment  assigned column for every row (Fortran, 1-based)
  !> @param[out] signs       sign of the selected overlap; zero maps to +1
  !> @param[out] matched     absolute selected overlap for every row
  !> @param[out] margins     selected weight minus best unselected row weight
  !> @param[out] info        zero on success; negative on invalid input
  subroutine maximum_overlap_assignment(overlap, assignment, signs, matched, margins, info)
    real(kind=dp), intent(in) :: overlap(:,:)
    integer, intent(out) :: assignment(:)
    real(kind=dp), intent(out) :: signs(:), matched(:), margins(:)
    integer, intent(out) :: info

    real(kind=dp), allocatable :: cost(:,:), u(:), v(:), minv(:)
    integer, allocatable :: p(:), way(:)
    logical, allocatable :: used(:)
    real(kind=dp) :: alt, cur, delta, selected
    integer :: i, i0, j, j0, j1, n

    n = size(overlap, 1)
    info = 0
    if (n <= 0 .or. size(overlap, 2) /= n .or. &
        size(assignment) /= n .or. size(signs) /= n .or. &
        size(matched) /= n .or. size(margins) /= n) then
      info = -1
      return
    end if
    if (.not. all(ieee_is_finite(overlap))) then
      info = -2
      return
    end if

    allocate(cost(n,n), u(0:n), v(0:n), minv(n), p(0:n), way(0:n), used(0:n))
    cost = -abs(overlap)
    u = 0.0_dp
    v = 0.0_dp
    p = 0
    way = 0

    ! Hungarian shortest-augmenting-path algorithm.  Columns are visited in
    ! ascending order, which gives deterministic tie handling.
    do i = 1, n
      p(0) = i
      minv = huge(1.0_dp)
      used = .false.
      j0 = 0
      do
        used(j0) = .true.
        i0 = p(j0)
        delta = huge(1.0_dp)
        j1 = 0
        do j = 1, n
          if (used(j)) cycle
          cur = cost(i0,j) - u(i0) - v(j)
          if (cur < minv(j)) then
            minv(j) = cur
            way(j) = j0
          end if
          if (minv(j) < delta) then
            delta = minv(j)
            j1 = j
          end if
        end do
        do j = 0, n
          if (used(j)) then
            u(p(j)) = u(p(j)) + delta
            v(j) = v(j) - delta
          else if (j > 0) then
            minv(j) = minv(j) - delta
          end if
        end do
        j0 = j1
        if (p(j0) == 0) exit
      end do

      do
        j1 = way(j0)
        p(j0) = p(j1)
        j0 = j1
        if (j0 == 0) exit
      end do
    end do

    assignment = 0
    do j = 1, n
      assignment(p(j)) = j
    end do

    do i = 1, n
      selected = overlap(i, assignment(i))
      matched(i) = abs(selected)
      if (selected < 0.0_dp) then
        signs(i) = -1.0_dp
      else
        signs(i) = 1.0_dp
      end if

      if (n == 1) then
        margins(i) = matched(i)
      else
        alt = 0.0_dp
        do j = 1, n
          if (j /= assignment(i)) alt = max(alt, abs(overlap(i,j)))
        end do
        margins(i) = matched(i) - alt
      end if
    end do
  end subroutine maximum_overlap_assignment

  !> C entry point.  overlap_row_major is a C/Python row-major n-by-n matrix;
  !> assignment is returned zero-based for direct NumPy indexing.
  function oqp_maximum_overlap_assignment_C(n, overlap_row_major, assignment, &
      signs, matched, margins) result(info) &
      bind(C, name="oqp_maximum_overlap_assignment")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: overlap_row_major(*)
    integer(c_int), intent(out) :: assignment(*)
    real(c_double), intent(out) :: signs(*), matched(*), margins(*)
    integer(c_int) :: info

    real(kind=dp), allocatable :: overlap(:,:), signs_f(:), matched_f(:), margins_f(:)
    integer, allocatable :: assignment_f(:)
    integer :: i, info_f, j, nn

    info = -1_c_int
    if (n <= 0_c_int) return
    nn = int(n, kind=kind(1))
    allocate(overlap(nn,nn), assignment_f(nn), signs_f(nn), matched_f(nn), margins_f(nn))

    do i = 1, nn
      do j = 1, nn
        overlap(i,j) = overlap_row_major((i-1)*nn + j)
      end do
    end do
    call maximum_overlap_assignment(overlap, assignment_f, signs_f, matched_f, margins_f, info_f)
    info = int(info_f, kind=c_int)
    if (info_f /= 0) return

    do i = 1, nn
      assignment(i) = int(assignment_f(i) - 1, kind=c_int)
      signs(i) = real(signs_f(i), kind=c_double)
      matched(i) = real(matched_f(i), kind=c_double)
      margins(i) = real(margins_f(i), kind=c_double)
    end do
  end function oqp_maximum_overlap_assignment_C

  !> C entry point for adiabatic, phase-only transport.  Unlike the global
  !> assignment above, this deliberately preserves the energy-root index and
  !> takes the sign from overlap(i,i).  It is the only safe phase operation
  !> when the caller does not also transport the state permutation through
  !> amplitudes, energies, active surface, and forces (as in ordinary NAMD).
  function oqp_diagonal_phase_tracking_C(n, overlap_row_major, signs, matched, margins) &
      result(info) bind(C, name="oqp_diagonal_phase_tracking")
    integer(c_int), value, intent(in) :: n
    real(c_double), intent(in) :: overlap_row_major(*)
    real(c_double), intent(out) :: signs(*), matched(*), margins(*)
    integer(c_int) :: info

    integer :: i, j, nn
    real(kind=dp) :: diagonal, alternative

    info = -1_c_int
    if (n <= 0_c_int) return
    nn = int(n, kind=kind(1))
    do i = 1, nn
      diagonal = overlap_row_major((i-1)*nn + i)
      if (.not. ieee_is_finite(diagonal)) then
        info = -2_c_int
        return
      end if
      if (diagonal < 0.0_dp) then
        signs(i) = -1.0_dp
      else
        signs(i) = 1.0_dp
      end if
      matched(i) = abs(diagonal)
      alternative = 0.0_dp
      do j = 1, nn
        if (j == i) cycle
        if (.not. ieee_is_finite(overlap_row_major((i-1)*nn + j))) then
          info = -2_c_int
          return
        end if
        alternative = max(alternative, abs(overlap_row_major((i-1)*nn + j)))
      end do
      margins(i) = matched(i) - alternative
    end do
    info = 0_c_int
  end function oqp_diagonal_phase_tracking_C

end module state_tracking_mod
