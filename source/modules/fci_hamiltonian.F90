!> @brief Fortran engine for the determinant CI hot loops of pyoqp fci.py.
!>
!> Three bind(C) entry points replace the pure-Python inner loops of the
!> native CASCI/FCI solver (the Python fallbacks remain and pin the numerics):
!>
!>  * oqp_dsyevd            -- dense symmetric eigensolver through the SAME
!>    ILP64 LAPACK liboqp links (mathlib eigen::diag_symm_full, DSYEVD),
!>    immune to the LP64-numpy-vs-ILP64-core ABI hazard the Python layer
!>    guards against with its Jacobi fallback.  Input is a symmetric matrix
!>    (storage order immaterial); output eigenvectors are LAPACK column-major
!>    (the Python caller transposes its C-order view).
!>  * fci_dense_hamiltonian -- assemble the dense determinant Hamiltonian from
!>    the spin-orbital integrals with the exact enumeration semantics of
!>    fci.py::_iter_hamiltonian_elements (one-electron a+_p a_q; two-electron
!>    0.5 (pq|rs) a+_p a+_q a_s a_r; bit-count phases; cutoff screening),
!>    symmetrized 0.5(H+H^T), OpenMP over columns.
!>  * fci_hamiltonian_matvec / fci_hamiltonian_diag -- matrix-free application
!>    Y = H X (block of vectors, symmetrized as 0.5(H_raw+H_raw^T) exactly
!>    like the dense path) and the diagonal <D|H|D> for the Davidson
!>    preconditioner.
!>
!> Determinants are the fci.py integer keys (alpha bits 0..norb-1, beta bits
!> norb..2norb-1; requires nspin <= 63).  hspin is C-order [nspin,nspin],
!> gspin C-order [nspin,nspin,nspin,nspin]; both are active-space sized by the
!> caller, so the nspin^4 storage is small.  Vector blocks X/Y are C-order
!> [ndet, nvec] numpy arrays (vector k = X[:, k]); flat index (det, k) =
!> det*nvec + k on both sides.
module fci_hamiltonian_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  !> Sort-once reuse for an in-process driver (fci_driver.F90).  The bind(C)
  !> entry points below re-sort the determinant keys on every call, which is
  !> correct for a Python caller that hands over a fresh list each time but is
  !> pure waste inside a Davidson loop that applies H twenty times to the same
  !> determinant space.  These expose the same kernels with the sort hoisted
  !> out; the bind(C) wrappers are thin shells over them, so both callers run
  !> byte-for-byte the same enumeration.
  public :: fci_sort_dets, fci_dense_build, fci_diag_build, fci_matvec_apply
  public :: oqp_dsyevd_f

contains

  !> Dense symmetric eigensolve for an in-process Fortran caller (the bind(C)
  !> `oqp_dsyevd` is a shell over this).
  subroutine oqp_dsyevd_f(n, a, w, info)
    use eigen, only: diag_symm_full
    integer, intent(in) :: n
    real(dp), intent(inout) :: a(*)
    real(dp), intent(out) :: w(*)
    integer, intent(out) :: info
    info = 0
    call diag_symm_full(1, n, a, n, w, info)
  end subroutine oqp_dsyevd_f

  !> Sort the CI-ordered determinant list into ascending keys, carrying a
  !> permutation back to the CI position.  `skeys(k)` is the k-th smallest key
  !> and `sperm(k)` its CI row.
  subroutine fci_sort_dets(ndet, dets, skeys, sperm)
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(ndet)
    integer(i8), intent(out) :: skeys(ndet), sperm(ndet)
    integer(i8) :: col
    do col = 1, ndet
      sperm(col) = col
    end do
    call sort_keys(ndet, dets, sperm)
    do col = 1, ndet
      skeys(col) = dets(sperm(col))
    end do
  end subroutine fci_sort_dets

  function oqp_dsyevd(n, a, w) result(info) bind(C, name="oqp_dsyevd")
    integer(c_int32_t), value :: n
    real(dp), intent(inout) :: a(*)
    real(dp), intent(out) :: w(*)
    integer(i8) :: info
    integer :: ierr
    call oqp_dsyevd_f(int(n), a, w, ierr)
    info = int(ierr, i8)
  end function oqp_dsyevd

  !------------------------------------------------------------------ helpers
  pure function pop_below(det, pos) result(n)
    integer(i8), intent(in) :: det
    integer, intent(in) :: pos
    integer :: n
    n = popcnt(iand(det, shiftl(1_i8, pos) - 1_i8))
  end function pop_below

  subroutine sort_keys(n, keys, perm)
    !! iterative quicksort of perm so that keys(perm) ascends
    integer(i8), intent(in) :: n
    integer(i8), intent(in) :: keys(n)
    integer(i8), intent(inout) :: perm(n)
    integer(i8) :: stack_lo(64), stack_hi(64)
    integer :: sp
    integer(i8) :: lo, hi, i, j, tmp, pivot
    if (n < 2) return
    sp = 1; stack_lo(1) = 1; stack_hi(1) = n
    do while (sp > 0)
      lo = stack_lo(sp); hi = stack_hi(sp); sp = sp - 1
      do while (lo < hi)
        if (hi - lo < 24) then
          do i = lo + 1, hi
            tmp = perm(i)
            j = i - 1
            do while (j >= lo)
              if (keys(perm(j)) <= keys(tmp)) exit
              perm(j + 1) = perm(j)
              j = j - 1
            end do
            perm(j + 1) = tmp
          end do
          exit
        end if
        pivot = keys(perm(ishft(lo + hi, -1)))
        i = lo; j = hi
        do
          do while (keys(perm(i)) < pivot)
            i = i + 1
          end do
          do while (keys(perm(j)) > pivot)
            j = j - 1
          end do
          if (i >= j) exit
          tmp = perm(i); perm(i) = perm(j); perm(j) = tmp
          i = i + 1; j = j - 1
        end do
        ! recurse into the smaller side, loop on the larger
        if (j - lo < hi - j) then
          if (j + 1 < hi) then
            sp = sp + 1; stack_lo(sp) = j + 1; stack_hi(sp) = hi
          end if
          hi = j
        else
          if (lo < j) then
            sp = sp + 1; stack_lo(sp) = lo; stack_hi(sp) = j
          end if
          lo = j + 1
        end if
      end do
    end do
  end subroutine sort_keys

  function find_det(n, sorted_keys, perm, key) result(idx)
    !! 1-based ORIGINAL index of key, or 0 when absent
    integer(i8), intent(in) :: n
    integer(i8), intent(in) :: sorted_keys(n), perm(n)
    integer(i8), intent(in) :: key
    integer(i8) :: idx, lo, hi, mid
    lo = 1; hi = n
    do while (lo <= hi)
      mid = ishft(lo + hi, -1)
      if (sorted_keys(mid) < key) then
        lo = mid + 1
      else if (sorted_keys(mid) > key) then
        hi = mid - 1
      else
        idx = perm(mid)
        return
      end if
    end do
    idx = 0
  end function find_det

  !--------------------------------------------------- shared element walker
  subroutine walk_column(col, nspin, ndet, dets, skeys, sperm, hspin, gspin, &
                         cutoff, hcol, want_col, xrow, yacc, nvec, x, ldv)
    !! Enumerate <row|H|col> contributions for one column determinant.
    !! Dense mode  (want_col): hcol(row) += val.
    !! Matvec mode (.not. want_col): yacc(:,row) += val * xrow(:) where
    !! xrow = x(:,col) was loaded by the caller, AND the transpose term
    !! yacc(:,col) += val * x(:,row) is accumulated for exact 0.5(H+H^T)
    !! symmetrization by the caller (handled via the second accumulation).
    integer(i8), intent(in) :: col, ndet
    integer, intent(in) :: nspin, nvec
    integer(i8), intent(in) :: dets(ndet), skeys(ndet), sperm(ndet)
    real(dp), intent(in) :: hspin(*), gspin(*)
    real(dp), intent(in) :: cutoff
    real(dp), intent(inout) :: hcol(*)
    logical, intent(in) :: want_col
    real(dp), intent(in) :: xrow(*)
    real(dp), intent(inout) :: yacc(*)
    real(dp), intent(in) :: x(*)
    integer(i8), intent(in) :: ldv

    integer(i8) :: det, det_q, det_rs, det_qrs, det_pqrs, row
    integer :: occ(nspin), nocc
    integer :: p, q, r, s, io, jo
    real(dp) :: val, ph_q, ph_r, ph_s, ph_p, ph_cq
    integer :: k

    det = dets(col)
    nocc = 0
    do p = 0, nspin - 1
      if (btest(det, p)) then
        nocc = nocc + 1
        occ(nocc) = p
      end if
    end do

    ! one-electron a+_p a_q
    do io = 1, nocc
      q = occ(io)
      det_q = ibclr(det, q)
      ph_q = merge(1.0_dp, -1.0_dp, mod(pop_below(det, q), 2) == 0)
      do p = 0, nspin - 1
        val = hspin(int(p, i8) * nspin + q + 1)
        if (abs(val) <= cutoff) cycle
        if (btest(det_q, p)) cycle
        ph_p = merge(1.0_dp, -1.0_dp, mod(pop_below(det_q, p), 2) == 0)
        row = find_det(ndet, skeys, sperm, ibset(det_q, p))
        if (row == 0) cycle
        call put(row, val * ph_q * ph_p)
      end do
    end do

    ! two-electron 0.5 g(p,q,r,s) a+_p a+_q a_s a_r
    do io = 1, nocc
      r = occ(io)
      ph_r = merge(1.0_dp, -1.0_dp, mod(pop_below(det, r), 2) == 0)
      det_rs = ibclr(det, r)
      do jo = 1, nocc
        s = occ(jo)
        if (s == r) cycle
        ph_s = merge(1.0_dp, -1.0_dp, mod(pop_below(ibclr(det, r), s), 2) == 0)
        det_qrs = ibclr(ibclr(det, r), s)
        do q = 0, nspin - 1
          if (btest(det_qrs, q)) cycle
          ph_cq = merge(1.0_dp, -1.0_dp, mod(pop_below(det_qrs, q), 2) == 0)
          det_pqrs = ibset(det_qrs, q)
          do p = 0, nspin - 1
            val = gspin((((int(p, i8) * nspin + q) * nspin + r) * nspin + s) + 1)
            if (abs(val) <= cutoff) cycle
            if (btest(det_pqrs, p)) cycle
            ph_p = merge(1.0_dp, -1.0_dp, mod(pop_below(det_pqrs, p), 2) == 0)
            row = find_det(ndet, skeys, sperm, ibset(det_pqrs, p))
            if (row == 0) cycle
            call put(row, 0.5_dp * val * ph_r * ph_s * ph_cq * ph_p)
          end do
        end do
      end do
    end do

  contains

    subroutine put(row_, val_)
      integer(i8), intent(in) :: row_
      real(dp), intent(in) :: val_
      if (want_col) then
        hcol(row_) = hcol(row_) + val_
      else
        ! forward term: y(:,row) += val * x(:,col)   (uses preloaded xrow)
        ! transpose term: y(:,col) += val * x(:,row) (for exact symmetrization,
        ! both scaled by 0.5)
        do k = 1, nvec
          yacc((row_ - 1) * ldv + k) = yacc((row_ - 1) * ldv + k) &
                                     + 0.5_dp * val_ * xrow(k)
          yacc((col - 1) * ldv + k) = yacc((col - 1) * ldv + k) &
                                    + 0.5_dp * val_ * x((row_ - 1) * ldv + k)
        end do
      end if
    end subroutine put

  end subroutine walk_column

  !----------------------------------------------------------- dense build
  subroutine fci_dense_hamiltonian(nspin, ndet, dets, hspin, gspin, cutoff, &
                                   hmat, nthreads) bind(C, name="fci_dense_hamiltonian")
    integer(c_int32_t), value :: nspin, nthreads
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(*)
    real(dp), intent(in) :: hspin(*), gspin(*)
    real(dp), value :: cutoff
    real(dp), intent(inout) :: hmat(*)

    integer(i8), allocatable :: skeys(:), sperm(:)

    allocate(skeys(ndet), sperm(ndet))
    call fci_sort_dets(ndet, dets(1:ndet), skeys, sperm)
    call fci_dense_build(int(nspin), ndet, dets, skeys, sperm, hspin, gspin, &
                         cutoff, hmat, int(nthreads))
    deallocate(skeys, sperm)
  end subroutine fci_dense_hamiltonian

  !> Dense determinant Hamiltonian from a determinant list already sorted by
  !> `fci_sort_dets`.  `hmat` is [ndet,ndet]; the result is symmetric so its
  !> C-order and column-major views coincide.
  subroutine fci_dense_build(nspin, ndet, dets, skeys, sperm, hspin, gspin, &
                             cutoff, hmat, nthreads)
    integer, intent(in) :: nspin, nthreads
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(*), skeys(ndet), sperm(ndet)
    real(dp), intent(in) :: hspin(*), gspin(*)
    real(dp), intent(in) :: cutoff
    real(dp), intent(inout) :: hmat(*)

    integer(i8) :: col, i, j
    integer :: nthr
    real(dp) :: dummy(1)

    hmat(1:ndet * ndet) = 0.0_dp
    nthr = max(1, nthreads)        ! explicit request wins over the ambient
                                   ! (BLAS-clamped) OMP default

    !$omp parallel do num_threads(nthr) schedule(dynamic, 16) default(shared) private(col)
    do col = 1, ndet
      call walk_column(col, nspin, ndet, dets, skeys, sperm, hspin, gspin, &
                       cutoff, hmat((col - 1) * ndet + 1), .true., &
                       dummy, dummy, 0, dummy, 0_i8)
    end do
    !$omp end parallel do

    ! symmetrize 0.5(H + H^T), matching the Python builder exactly
    do j = 1, ndet
      do i = j + 1, ndet
        hmat((j - 1) * ndet + i) = 0.5_dp * (hmat((j - 1) * ndet + i) &
                                           + hmat((i - 1) * ndet + j))
        hmat((i - 1) * ndet + j) = hmat((j - 1) * ndet + i)
      end do
    end do
  end subroutine fci_dense_build

  !----------------------------------------------------------- diagonal
  subroutine fci_hamiltonian_diag(nspin, ndet, dets, hspin, gspin, diag) &
      bind(C, name="fci_hamiltonian_diag")
    integer(c_int32_t), value :: nspin
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(*)
    real(dp), intent(in) :: hspin(*), gspin(*)
    real(dp), intent(out) :: diag(*)
    call fci_diag_build(int(nspin), ndet, dets, hspin, gspin, diag)
  end subroutine fci_hamiltonian_diag

  !> <D|H|D> over the determinant list -- the Davidson preconditioner.
  subroutine fci_diag_build(nspin, ndet, dets, hspin, gspin, diag)
    integer, intent(in) :: nspin
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(*)
    real(dp), intent(in) :: hspin(*), gspin(*)
    real(dp), intent(out) :: diag(*)
    integer(i8) :: col
    integer :: occ(nspin), nocc, p, io, jo, ns
    integer(i8) :: det
    real(dp) :: e
    ns = nspin
    do col = 1, ndet
      det = dets(col)
      nocc = 0
      do p = 0, ns - 1
        if (btest(det, p)) then
          nocc = nocc + 1
          occ(nocc) = p
        end if
      end do
      e = 0.0_dp
      do io = 1, nocc
        e = e + hspin(int(occ(io), i8) * ns + occ(io) + 1)
      end do
      do io = 1, nocc
        do jo = 1, nocc
          if (io == jo) cycle
          ! 0.5 [ g(p,q,p,q) - g(p,q,q,p) ] over occupied pairs
          e = e + 0.5_dp * ( &
              gspin((((int(occ(io), i8) * ns + occ(jo)) * ns + occ(io)) * ns + occ(jo)) + 1) &
            - gspin((((int(occ(io), i8) * ns + occ(jo)) * ns + occ(jo)) * ns + occ(io)) + 1))
        end do
      end do
      diag(col) = e
    end do
  end subroutine fci_diag_build

  !----------------------------------------------------------- matvec block
  function fci_hamiltonian_matvec(nspin, ndet, dets, hspin, gspin, cutoff, &
                                  nvec, x, y, nthreads) &
      result(status) bind(C, name="fci_hamiltonian_matvec")
    integer(c_int32_t), value :: nspin, nvec, nthreads
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(*)
    real(dp), intent(in) :: hspin(*), gspin(*)
    real(dp), value :: cutoff
    real(dp), intent(in) :: x(*)
    real(dp), intent(inout) :: y(*)
    integer(i8) :: status

    integer(i8), allocatable :: skeys(:), sperm(:)

    allocate(skeys(ndet), sperm(ndet))
    call fci_sort_dets(ndet, dets(1:ndet), skeys, sperm)
    call fci_matvec_apply(int(nspin), ndet, dets, skeys, sperm, hspin, gspin, &
                          cutoff, int(nvec), x, y, int(nthreads))
    deallocate(skeys, sperm)
    status = 0_i8
  end function fci_hamiltonian_matvec

  !> Matrix-free block application Y = 0.5(H+H^T) X from a pre-sorted
  !> determinant list.  X and Y are C-order [ndet, nvec].
  !>
  !> The string-driven engine handles the canonical CAS product list -- which is
  !> every CASCI/CASSCF/FCI solve -- in `O(ndet * nact**2)` BLAS-3 work.  It
  !> declines any other list (a spin-filtered or restricted determinant set) and
  !> the determinant-pair walker below then runs exactly as before.  Setting
  !> `OQP_FCI_SIGMA=walk` forces the walker, which is what the two paths are
  !> compared with.
  subroutine fci_matvec_apply(nspin, ndet, dets, skeys, sperm, hspin, gspin, &
                              cutoff, nvec, x, y, nthreads)
    use fci_sigma_strings_mod, only: fci_sigma_strings
    integer, intent(in) :: nspin, nvec, nthreads
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(*), skeys(ndet), sperm(ndet)
    real(dp), intent(in) :: hspin(*), gspin(*)
    real(dp), intent(in) :: cutoff
    real(dp), intent(in) :: x(*)
    real(dp), intent(inout) :: y(*)

    real(dp), allocatable :: yloc(:)
    real(dp) :: xrow(nvec)
    integer(i8) :: col, i
    integer :: nthr, ns, k
    real(dp) :: dummy(1)
    character(len=8) :: sigma_mode
    integer :: env_len, env_stat

    sigma_mode = ' '
    call get_environment_variable("OQP_FCI_SIGMA", sigma_mode, env_len, env_stat)
    if (env_stat /= 0 .or. sigma_mode /= 'walk') then
      if (fci_sigma_strings(nspin, ndet, dets, hspin, gspin, nvec, x, y, &
                            nthreads) == 0) return
    end if

    ns = nspin
    y(1:ndet * nvec) = 0.0_dp
    nthr = max(1, nthreads)        ! explicit request wins over the ambient
                                   ! (BLAS-clamped) OMP default

    !$omp parallel num_threads(nthr) default(shared) private(yloc, xrow, col, i, k)
    allocate(yloc(ndet * nvec))
    yloc = 0.0_dp
    !$omp do schedule(dynamic, 16)
    do col = 1, ndet
      do k = 1, int(nvec)
        xrow(k) = x((col - 1) * nvec + k)
      end do
      call walk_column(col, ns, ndet, dets, skeys, sperm, hspin, gspin, &
                       cutoff, dummy, .false., xrow, yloc, int(nvec), x, &
                       int(nvec, i8))
    end do
    !$omp end do
    !$omp critical (fci_matvec_merge)
    do i = 1, ndet * nvec
      y(i) = y(i) + yloc(i)
    end do
    !$omp end critical (fci_matvec_merge)
    deallocate(yloc)
    !$omp end parallel
  end subroutine fci_matvec_apply

end module fci_hamiltonian_mod
