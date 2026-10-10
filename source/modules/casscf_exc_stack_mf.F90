!> @brief Matrix-free (on-the-fly) evaluation of the spin-summed excitation
!> matrices needed for the analytic CASSCF orbital-rotation Hessian.
!>
!> The original `casscf_exc_stack.F90` builds the dense `[nact,nact,ndet,ndet]`
!> stack, O(nact² × ndet²) memory despite density ≈ 2/ndet (~0.04% at ndet=5000).
!> These routines evaluate the same products without ever storing the full
!> tensor, by walking the non-zero entries directly.
!>
!> Determinants are the fci.py integer keys (alpha bits 0..nact-1, beta bits
!> nact..2*nact-1; requires 2*nact <= 62).
module casscf_exc_stack_mf_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: mf_sort_keys_perm, mf_bsearch
  public :: casscf_exc_stack_apply_wmat

contains

  !> Shell sort for up to i8-sized arrays, carrying a permutation alongside.
  subroutine mf_sort_keys_perm(n, keys, perm)
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
  end subroutine mf_sort_keys_perm

  !> Binary search returning position or -1.
  pure subroutine mf_bsearch(n, keys, key, pos)
    integer(i8), intent(in) :: n, keys(0:), key
    integer(i8), intent(out) :: pos
    integer(i8) :: lo, hi, mid
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
  end subroutine mf_bsearch

  !> W_tua = (E_tu|c>)_a on the fly, replacing stack^T @ civec DGEMV.
  !>
  !> Enumerates non-zero entries of E_tu(a,b) and accumulates
  !> wmat[(t,u),a] += E_tu(a,b) * civec(b).
  !>
  !> @param[in]  nact   active orbitals
  !> @param[in]  ndet   determinants
  !> @param[in]  dets   determinant keys in CI order
  !> @param[in]  skeys  sorted determinant keys
  !> @param[in]  sperm  permutation back to CI ordering
  !> @param[in]  civec  CI vector
  !> @param[out] wmat   C-order [nact,nact,ndet]
  subroutine casscf_exc_stack_apply_wmat(nact, ndet, dets, skeys, sperm, &
                                         civec, wmat) &
      bind(C, name="casscf_exc_stack_apply_wmat")
    integer(c_int32_t), value :: nact
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(0:ndet-1), skeys(0:ndet-1), sperm(0:ndet-1)
    real(dp), intent(in) :: civec(0:ndet-1)
    real(dp), intent(inout) :: wmat(0:*)

    integer :: na
    integer(i8) :: n2, col, det, det_u, det_tu, ubit, tbit, row
    integer :: off, ioff, t, u, phase_u, phase_t, offs(2)
    real(dp) :: ci

    na = int(nact)
    n2 = int(na, i8) * int(na, i8)
    if (na <= 0 .or. ndet <= 0_i8) return
    if (2 * na > 62) return

    offs(1) = 0
    offs(2) = na

    ! Zero output (before any early return so wmat is always safe)
    wmat(0:n2*ndet - 1_i8) = 0.0_dp

    ! Bra-side loop: iterate over bra (row), enumerate de-excitations to find
    ! ket (col).  Each (tu,row) slot is written by exactly one thread (the owner
    ! of that row), so no atomic update is needed.
    !$omp parallel do default(shared) schedule(static) if(ndet >= 64_i8) &
    !$omp   private(row, det, ioff, off, t, tbit, phase_t, det_t, &
    !$omp           u, ubit, phase_u, det_col, col_sorted, col)
    do row = 0_i8, ndet - 1_i8
      det = dets(row)
      do ioff = 1, 2
        off = offs(ioff)
        ! Annihilate t from row: t is the creation index of E_tu
        do t = 0, na - 1
          tbit = ishft(1_i8, t + off)
          if (iand(det, tbit) == 0_i8) cycle
          phase_t = 1
          if (mod(popcnt(iand(det, tbit - 1_i8)), 2) /= 0) phase_t = -1
          det_t = ieor(det, tbit)
          ! Create u in det_t: u is the annihilation index of E_tu
          do u = 0, na - 1
            ubit = ishft(1_i8, u + off)
            if (iand(det_t, ubit) /= 0_i8) cycle
            phase_u = 1
            if (mod(popcnt(iand(det_t, ubit - 1_i8)), 2) /= 0) phase_u = -1
            ! |col> = E_ut|row> = +/- a^+_u a_t |row>
            det_col = ior(det_t, ubit)
            call mf_bsearch(ndet, skeys, det_col, col_sorted)
            if (col_sorted < 0_i8) cycle
            col = sperm(col_sorted)
            ! wmat(tu, bra=row) += E_tu(bra=row, ket=col) * civec(ket=col)
            ! Each thread owns its rows, so no atomic needed.
            wmat((int(t, i8)*int(na, i8) + int(u, i8))*ndet + row) = &
                wmat((int(t, i8)*int(na, i8) + int(u, i8))*ndet + row) &
                + real(phase_t * phase_u, dp) * civec(col)
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine casscf_exc_stack_apply_wmat

end module casscf_exc_stack_mf_mod
