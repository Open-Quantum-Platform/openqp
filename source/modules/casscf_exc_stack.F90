!> @brief Fortran engine for the spin-summed excitation-matrix stack of
!> pyoqp casscf_hessian.py.
!>
!> The analytic CASSCF Hessian evaluates its CI-relaxation half against
!>
!>     (E_tu)_{row,col} = <row| E_tu |col>,   E_tu = sum_sigma a+_{t sigma} a_{u sigma}
!>
!> over the determinant list `fci._determinants`.  The Python builder walks
!> every determinant and applies each annihilation/creation pair one element at
!> a time, which is the last pure-Python arithmetic left on this module's
!> execution path.  This replaces it (the Python builder stays as the numerical
!> pin).
!>
!> Each column is independent, so the whole stack is a single sweep over the
!> determinants with no accumulation between columns.  Lookups of the excited
!> determinant use the same sorted-key binary search idiom as
!> rdm_kernel.F90 -- and the same caution applies:
!>
!>   the determinant list arrives in **CI ordering, not key ordering**, so it
!>   must never be binary-searched directly.  A sorted copy carries a
!>   permutation back to the CI position, which is what indexes the output.
!>
!> Determinants are the fci.py integer keys (alpha bits 0..nact-1, beta bits
!> nact..2*nact-1; requires 2*nact <= 62).  The output is C-order
!> [nact,nact,ndet,ndet] with the last index fastest, matching the NumPy array
!> the caller allocates.  Note that array is nact^2 * ndet^2 doubles, which the
!> caller is expected to have size-guarded before calling.
module casscf_exc_stack_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: casscf_excitation_stack

contains

  !> Ascending sort of `keys` carrying `perm` alongside, so a search hit maps
  !> back to the original (CI-ordered) determinant position.
  subroutine sort_keys_perm(n, keys, perm)
    integer(i8), intent(in) :: n
    integer(i8), intent(inout) :: keys(0:), perm(0:)
    integer(i8) :: gap, i, j, key, pos
    gap = n / 2_i8
    do while (gap > 0_i8)
      do i = gap, n - 1_i8
        key = keys(i)
        pos = perm(i)
        j = i
        do while (j >= gap)
          if (keys(j - gap) <= key) exit
          keys(j) = keys(j - gap)
          perm(j) = perm(j - gap)
          j = j - gap
        end do
        keys(j) = key
        perm(j) = pos
      end do
      gap = gap / 2_i8
    end do
  end subroutine sort_keys_perm

  !> Position of `key` in the ascending array, or -1.
  pure function bsearch(n, keys, key) result(pos)
    integer(i8), intent(in) :: n, keys(0:), key
    integer(i8) :: pos, lo, hi, mid
    lo = 0_i8
    hi = n - 1_i8
    pos = -1_i8
    do while (lo <= hi)
      mid = (lo + hi) / 2_i8
      if (keys(mid) == key) then
        pos = mid
        return
      else if (keys(mid) < key) then
        lo = mid + 1_i8
      else
        hi = mid - 1_i8
      end if
    end do
  end function bsearch

  !> Spin-summed excitation matrices, C-order [nact,nact,ndet,ndet].
  !>
  !> @param[in]  nact   active orbitals (requires 2*nact <= 62)
  !> @param[in]  ndet   determinants
  !> @param[in]  dets   determinant keys in CI order
  !> @param[out] stack  (E_tu)_{row,col}
  !> Returns 0 on success, -1 if the bit encoding does not fit.
  function casscf_excitation_stack(nact, ndet, dets, stack) result(info) &
      bind(C, name="casscf_excitation_stack")
    integer(c_int32_t), value :: nact
    integer(i8), value :: ndet
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    integer(i8), intent(in) :: dets(0:*)
    real(dp), intent(inout) :: stack(0:*)
    integer(i8) :: info

    integer :: na, t, u, off, ioff, phase_u, phase_t
    integer :: offs(2)
    integer(i8) :: col, row, spos, det, det_u, det_tu, ubit, tbit, n2
    integer(i8), allocatable :: skeys(:), sperm(:)

    na = int(nact)
    info = 0_i8
    if (2 * na > 62) then
      info = -1_i8
      return
    end if
    if (ndet <= 0_i8 .or. na <= 0) return

    offs(1) = 0
    offs(2) = na
    n2 = ndet * ndet
    do col = 0_i8, int(na, i8)*int(na, i8)*n2 - 1_i8
      stack(col) = 0.0_dp
    end do

    ! The CI ordering is not key ordering, so sort a copy and keep the map back.
    allocate(skeys(0:ndet-1), sperm(0:ndet-1))
    do col = 0_i8, ndet - 1_i8
      skeys(col) = dets(col)
      sperm(col) = col
    end do
    call sort_keys_perm(ndet, skeys, sperm)

    do col = 0_i8, ndet - 1_i8
      det = dets(col)
      do ioff = 1, 2                ! alpha block, then beta block
        off = offs(ioff)
        do u = 0, na - 1
          ubit = ishft(1_i8, u + off)
          if (iand(det, ubit) == 0_i8) cycle
          phase_u = 1
          if (mod(popcnt(iand(det, ubit - 1_i8)), 2) /= 0) phase_u = -1
          det_u = ieor(det, ubit)
          do t = 0, na - 1
            tbit = ishft(1_i8, t + off)
            if (iand(det_u, tbit) /= 0_i8) cycle
            phase_t = 1
            if (mod(popcnt(iand(det_u, tbit - 1_i8)), 2) /= 0) phase_t = -1
            det_tu = ior(det_u, tbit)
            spos = bsearch(ndet, skeys, det_tu)
            if (spos < 0_i8) cycle
            row = sperm(spos)
            stack((int(t, i8)*int(na, i8) + int(u, i8))*n2 + row*ndet + col) = &
                stack((int(t, i8)*int(na, i8) + int(u, i8))*n2 + row*ndet + col) &
                + real(phase_u * phase_t, dp)
          end do
        end do
      end do
    end do

    deallocate(skeys, sperm)
  end function casscf_excitation_stack

end module casscf_exc_stack_mod
