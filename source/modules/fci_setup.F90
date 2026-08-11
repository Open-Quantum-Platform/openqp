!> @brief Fortran engine for the active-space setup of pyoqp `fci.py`:
!> the spin-orbital expansion of the integrals and the frozen-core fold.
!>
!> Both sit between the AO->MO transformation (mo_transform.F90) and the
!> determinant Hamiltonian (fci_hamiltonian.F90), and both ran in Python on
!> every CI solve -- which for CASSCF means once per gradient evaluation, so
!> thousands of times per run.  The Python implementations stay in the tree as
!> the numerical pins and as the fallback when these symbols are unavailable.
!>
!> Two bind(C) entry points:
!>
!>  * fci_spin_orbital_integrals -- spatial (pq|rs) -> spin-orbital gspin
!>  * fci_fold_core              -- frozen-core energy and the folded active h
!>
!> Neither does any arithmetic that can lose precision.  The spin expansion is
!> a permuted copy: gspin[p,q,r,s] = (p r | q s) when spin(p)==spin(r) and
!> spin(q)==spin(s), zero otherwise, so exactly four of the sixteen spin blocks
!> are populated and each is the same spatial tensor with its second and third
!> indices swapped.  The fold is a sum of ERI elements over the core.  Both are
!> therefore expected to agree with the Python reference bit for bit, and the
!> tests assert exactly that rather than a tolerance.
!>
!> All arrays cross the boundary in the caller's C order (last index fastest).
!> The spin expansion is written directly in C-order index arithmetic (the
!> output is not symmetric under a full index reversal, so the Fortran
!> column-major view is NOT the same tensor and cannot be used here).
module fci_setup_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
!$ use omp_lib
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: fci_spin_orbital_integrals, fci_fold_core, fci_spin_square

contains

  !> Threads for one region, given the caller's request: the request when it
  !> is positive, otherwise whatever the environment already allows.
  !>
  !> This exists so the request can be expressed as a `num_threads` clause.
  !> The clause is scoped to its own region; omp_set_num_threads is not, and
  !> using it here leaked the FCI engine's thread cap into every subsequent
  !> parallel region in the process.
  function fci_region_threads(nthreads) result(n)
    ! c_int32_t, not the default integer: both callers take it straight off a
    ! bind(C) `value` dummy, and OpenQP's default integer is 8 bytes (ILP64).
    integer(c_int32_t), intent(in) :: nthreads
    integer :: n
    n = 1
!$  n = omp_get_max_threads()
    if (nthreads > 0) n = int(nthreads)
    n = max(1, n)
  end function fci_region_threads

  !--------------------------------------------------- determinant lookup
  !> Number of set bits below `pos`, i.e. the fermionic phase counter.
  pure function pop_below(det, pos) result(n)
    integer(i8), intent(in) :: det
    integer, intent(in) :: pos
    integer :: n
    n = popcnt(iand(det, shiftl(1_i8, pos) - 1_i8))
  end function pop_below

  !> Iterative quicksort of `perm` so that keys(perm) ascends.
  !>
  !> The determinant list arrives in CI ordering, NOT sorted by key, so it can
  !> never be binary-searched directly; the permutation carries each sorted
  !> entry back to its CI position.
  subroutine sort_keys(n, keys, perm)
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

  !> 1-based ORIGINAL (CI) index of `key`, or 0 when absent.
  pure function find_det(n, sorted_keys, perm, key) result(idx)
    integer(i8), intent(in) :: n
    integer(i8), intent(in) :: sorted_keys(n), perm(n)
    integer(i8), intent(in) :: key
    integer(i8) :: idx, lo, hi, mid
    lo = 1; hi = n
    idx = 0
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
  end function find_det

  !> Spin-orbital expansion of the spatial integrals.
  !>
  !> @param[in]  norb   spatial orbitals in the active space
  !> @param[in]  h1e    spatial one-electron integrals, C-order [norb,norb]
  !> @param[in]  eri    spatial (pq|rs), C-order [norb,norb,norb,norb]
  !> Both outputs are written in FULL -- every element, the structural zeros
  !> included -- so the caller passes uninitialized buffers (np.empty).
  !>
  !> @param[out] hspin  C-order [2*norb, 2*norb]
  !> @param[out] gspin  C-order [2*norb, 2*norb, 2*norb, 2*norb]
  !> @param[in]  nthreads  OpenMP threads for the block fill (0 = leave as is)
  subroutine fci_spin_orbital_integrals(norb, h1e, eri, hspin, gspin, nthreads) &
      bind(C, name="fci_spin_orbital_integrals")
    integer(c_int32_t), value :: norb, nthreads
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    real(dp), intent(in) :: h1e(0:*), eri(0:*)
    real(dp), intent(inout) :: hspin(0:*), gspin(0:*)

    integer :: n, ns, p, q, r, s, ps, qs
    integer(i8) :: ns2, ns3, ns4, base_p, base_pq, base_pqr
    integer(i8) :: n2, n3, esrc

    n = int(norb)
    if (n <= 0) return
    ns = 2 * n
    n2 = int(n, i8) * int(n, i8)
    n3 = n2 * int(n, i8)
    ns2 = int(ns, i8) * int(ns, i8)
    ns3 = ns2 * int(ns, i8)
    ns4 = ns3 * int(ns, i8)

    hspin(0:ns2 - 1_i8) = 0.0_dp
    do p = 0, n - 1
      do q = 0, n - 1
        hspin(int(p, i8) * int(ns, i8) + int(q, i8)) = h1e(int(p, i8)*int(n, i8) + int(q, i8))
        hspin(int(n + p, i8) * int(ns, i8) + int(n + q, i8)) = &
            h1e(int(p, i8)*int(n, i8) + int(q, i8))
      end do
    end do

    ! One bulk clear, then the four populated spin blocks scattered into it.
    !
    ! The clear looks like the wasteful half -- it writes 16*norb^4 doubles for
    ! 4*norb^4 doubles of payload -- but it is a single memset and runs at better
    ! than 100 GB/s, so at norb=8 it costs about 4 us of the kernel's 7.  The
    ! alternative, folding the zeros into the copy so the tensor is written
    ! exactly once, was MEASURED AND IS 2.5x SLOWER (16.7 us against 6.7 us at
    ! norb=8): the structural zeros come in runs of `norb`, so a single pass
    ! issues thousands of 64-byte memsets instead of one large one, and the
    ! per-call overhead swamps the bytes it saves.  Do not "optimise" this back
    ! into one pass without re-measuring.
    !
    ! What DID cost this kernel its deficit against numpy was doing the clear
    ! TWICE: the caller allocated with np.zeros, so the buffer was memset in C
    ! and then memset again here, while the numpy fallback it is racing writes
    ! the array only once.  The caller now passes np.empty and this is the only
    ! clear -- so every element of both outputs has to be written here.  Do not
    ! add an early return inside this nest without re-establishing that.
    !
    ! Layout is unchanged and no arithmetic happens in either version, so the
    ! result stays bit-identical to _spin_orbital_integrals_reference:
    ! gspin[p,q,r,s] = (p r | q s) when spin(p)==spin(r) and spin(q)==spin(s).
    gspin(0:ns4 - 1_i8) = 0.0_dp

    ! `nthreads` applies to THIS region only, through the num_threads clause.
    ! It used to be applied with omp_set_num_threads and never put back, and
    ! omp_set_num_threads is process-wide: the FCI caller asks for
    ! min(8, cpu_count) threads (fci._fci_lib_threads), so on any machine with
    ! more than eight cores the first CI solve silently pinned EVERY later
    ! OpenMP region in the process -- the SCF, the integral transforms, the
    ! analytic Hessian -- to eight threads.  It made a whole-job thread scan
    ! read as perfectly flat: H2O/cc-pVTZ CAS(6,6) took 12.0 s at
    ! OMP_NUM_THREADS=1 and 12.0 s at 128, with the Hessian kernel reporting
    ! nthr=8 in every one of them.
    !$omp parallel do collapse(2) default(shared) if(ns4 >= 1048576_i8) &
    !$omp   num_threads(fci_region_threads(nthreads)) &
    !$omp   private(p, q, r, s, ps, qs, base_p, base_pq, base_pqr, esrc)
    do p = 0, ns - 1
      do q = 0, ns - 1
        ps = p / n
        qs = q / n
        base_p = int(p, i8) * ns3
        base_pq = base_p + int(q, i8) * ns2
        ! spin(r) must equal spin(p), spin(s) must equal spin(q); both the
        ! destination run and the eri run are unit-stride in s
        do r = ps * n, ps * n + n - 1
          base_pqr = base_pq + int(r, i8) * int(ns, i8) + int(qs * n, i8)
          esrc = int(mod(p, n), i8) * n3 + int(mod(r, n), i8) * n2 &
                 + int(mod(q, n), i8) * int(n, i8)
          do s = 0, n - 1
            gspin(base_pqr + int(s, i8)) = eri(esrc + int(s, i8))
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine fci_spin_orbital_integrals

  !> Frozen-core fold: the inactive energy and the effective active h.
  !>
  !> The core orbitals are doubly occupied and frozen, so they contribute a
  !> constant
  !>
  !>     2 sum_i h[i,i] + sum_ij ( 2 (ii|jj) - (ij|ji) )
  !>
  !> and an effective one-electron operator on the active space
  !>
  !>     h_act[p,q] += sum_i ( 2 (pq|ii) - (pi|iq) ).
  !>
  !> Both index lists are explicit so the sequential and the explicitly
  !> selected CAS paths share one implementation.
  !>
  !> @param[in]  norb      full MO count
  !> @param[in]  nact      active orbitals
  !> @param[in]  ncore     frozen core orbitals
  !> @param[in]  active    active MO indices (0-based), [nact]
  !> @param[in]  core      core MO indices (0-based), [ncore]
  !> @param[in]  h1e       full MO one-electron integrals, C-order [norb,norb]
  !> @param[in]  eri       full MO (pq|rs), C-order [norb,norb,norb,norb]
  !> @param[out] h_act     folded active h, C-order [nact,nact]
  !> @param[inout] eref    on entry the caller's reference energy, on exit that
  !>                       energy plus the inactive contribution.  Accumulating
  !>                       into the caller's starting value rather than summing
  !>                       separately and adding keeps the rounding identical to
  !>                       the Python reference.
  subroutine fci_fold_core(norb, nact, ncore, active, core, h1e, eri, &
                           h_act, eref) bind(C, name="fci_fold_core")
    integer(c_int32_t), value :: norb, nact, ncore
    integer(c_int32_t), intent(in) :: active(0:*), core(0:*)
    real(dp), intent(in) :: h1e(0:*), eri(0:*)
    real(dp), intent(inout) :: h_act(0:*), eref(0:*)

    integer :: n, na, nc, p, q, i, j, pp, qq, ii, jj
    integer(i8) :: n2, n3, acc
    real(dp) :: ecore, term

    n = int(norb)
    na = int(nact)
    nc = int(ncore)
    n2 = int(n, i8) * int(n, i8)
    n3 = n2 * int(n, i8)

    do p = 0, na - 1
      pp = int(active(p))
      do q = 0, na - 1
        qq = int(active(q))
        h_act(int(p, i8)*int(na, i8) + int(q, i8)) = &
            h1e(int(pp, i8)*int(n, i8) + int(qq, i8))
      end do
    end do

    ecore = eref(0)
    if (nc <= 0) return

    do i = 0, nc - 1
      ii = int(core(i))
      ecore = ecore + 2.0_dp * h1e(int(ii, i8)*int(n, i8) + int(ii, i8))
    end do
    do i = 0, nc - 1
      ii = int(core(i))
      do j = 0, nc - 1
        jj = int(core(j))
        ! 2 (ii|jj) - (ij|ji).  The term is formed before it is accumulated so
        ! the rounding matches the Python reference's `ecore += 2*a - b`
        ! exactly; `ecore + 2*a - b` would associate the other way.
        term = 2.0_dp * eri(int(ii, i8)*n3 + int(ii, i8)*n2 &
                            + int(jj, i8)*int(n, i8) + int(jj, i8)) &
               - eri(int(ii, i8)*n3 + int(jj, i8)*n2 &
                     + int(jj, i8)*int(n, i8) + int(ii, i8))
        ecore = ecore + term
      end do
    end do

    do p = 0, na - 1
      pp = int(active(p))
      do q = 0, na - 1
        qq = int(active(q))
        acc = int(p, i8)*int(na, i8) + int(q, i8)
        do i = 0, nc - 1
          ii = int(core(i))
          ! 2 (pq|ii) - (pi|iq), formed before accumulation (see above)
          term = 2.0_dp * eri(int(pp, i8)*n3 + int(qq, i8)*n2 &
                              + int(ii, i8)*int(n, i8) + int(ii, i8)) &
                 - eri(int(pp, i8)*n3 + int(ii, i8)*n2 &
                       + int(ii, i8)*int(n, i8) + int(qq, i8))
          h_act(acc) = h_act(acc) + term
        end do
      end do
    end do

    eref(0) = ecore
  end subroutine fci_fold_core

  !> <S^2> for a block of CI vectors in the determinant basis.
  !>
  !>     <S^2> = ms(ms+1) + sum_{pq} <psi| a+_{p,beta} a_{p,alpha}
  !>                                       a+_{q,alpha} a_{q,beta} |psi>
  !>
  !> The Python reference walks that double sum per determinant per root with
  !> `bin(x).count("1")` phase counts and a dict lookup per term, which is
  !> O(nroot * ndet * norb^2) Python: 48 ms per root at CAS(8,8) and 843 ms per
  !> root at CAS(10,10).
  !>
  !> The determinant list is in CI order, so it is sorted once into a key array
  !> carrying a permutation back to CI position and binary-searched from there.
  !> The per-root loop shares that sort.
  !>
  !> @param[in]  norb     spatial orbitals (alpha bits 0..norb-1, beta above)
  !> @param[in]  ndet     determinant count
  !> @param[in]  nvec     number of CI vectors
  !> @param[in]  dets     determinant bit patterns in CI order, [ndet]
  !> @param[in]  civec    CI vectors, C-order [ndet, nvec] (column per root)
  !> @param[in]  ms2      2*ms = nalpha - nbeta
  !> @param[out] s2       <S^2> per root, [nvec]
  !> @return     0 on success, -1 when the scratch arrays cannot be allocated
  function fci_spin_square(norb, ndet, nvec, dets, civec, ms2, s2, nthreads) &
      result(status) bind(C, name="fci_spin_square")
    integer(c_int32_t), value :: norb, nvec, ms2, nthreads
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(0:*)
    real(dp), intent(in) :: civec(0:*)
    real(dp), intent(inout) :: s2(0:*)
    integer(i8) :: status

    integer :: n, nv, p, q, root, ierr, phase
    integer(i8) :: col, row, det, det_q, det_aq, det_p, mask
    integer(i8), allocatable :: keys(:), skeys(:), perm(:)
    real(dp), allocatable :: norms(:), cvec(:)
    real(dp) :: ms, base, flip, c_col, acc, term

    status = 0_i8
    n = int(norb)
    nv = int(nvec)
    if (n <= 0 .or. ndet <= 0 .or. nv <= 0) return

    allocate(keys(ndet), skeys(ndet), perm(ndet), norms(0:nv-1), &
             cvec(0:ndet*int(nv, i8) - 1_i8), stat=ierr)
    if (ierr /= 0) then
      status = -1_i8
      return
    end if

    ! The determinant list is in CI order, so it is sorted into a key array
    ! carrying a permutation back to the CI position; the sort is shared by
    ! every root.
    do col = 1_i8, ndet
      keys(col) = dets(col - 1_i8)
      perm(col) = col
    end do
    call sort_keys(ndet, keys, perm)
    do col = 1_i8, ndet
      skeys(col) = keys(perm(col))
    end do

    ms = 0.5_dp * real(ms2, dp)
    base = ms * (ms + 1.0_dp)

    ! The reference normalizes the CI vector before accumulating; doing the
    ! same here (rather than scaling the accumulated term afterwards) keeps the
    ! arithmetic in the same order, so the two agree to the last bit instead of
    ! differing by the sign of a 1e-17 spin-flip residual.
    do root = 0, nv - 1
      acc = 0.0_dp
      do col = 0_i8, ndet - 1_i8
        term = civec(col * int(nv, i8) + int(root, i8))
        acc = acc + term * term
      end do
      norms(root) = sqrt(acc)
      if (norms(root) > 0.0_dp) then
        do col = 0_i8, ndet - 1_i8
          cvec(col * int(nv, i8) + int(root, i8)) = &
              civec(col * int(nv, i8) + int(root, i8)) / norms(root)
        end do
      end if
    end do

    ! Region-local thread count -- see the note at the other call site.
    !$omp parallel do default(shared) if(ndet >= 512_i8) schedule(static) &
    !$omp   num_threads(fci_region_threads(nthreads)) &
    !$omp   private(root, col, det, det_q, det_aq, det_p, row, &
    !$omp           p, q, phase, mask, c_col, flip)
    do root = 0, nv - 1
      flip = 0.0_dp
      if (norms(root) > 0.0_dp) then
        do col = 0_i8, ndet - 1_i8
          c_col = cvec(col * int(nv, i8) + int(root, i8))
          if (c_col == 0.0_dp) cycle
          det = dets(col)
          do q = 0, n - 1
            ! a_{q,beta}
            mask = shiftl(1_i8, n + q)
            if (iand(det, mask) == 0_i8) cycle
            det_q = ieor(det, mask)
            ! a+_{q,alpha}
            mask = shiftl(1_i8, q)
            if (iand(det_q, mask) /= 0_i8) cycle
            det_aq = ior(det_q, mask)
            do p = 0, n - 1
              ! a_{p,alpha}
              mask = shiftl(1_i8, p)
              if (iand(det_aq, mask) == 0_i8) cycle
              det_p = ieor(det_aq, mask)
              ! a+_{p,beta}
              mask = shiftl(1_i8, n + p)
              if (iand(det_p, mask) /= 0_i8) cycle
              row = find_det(ndet, skeys, perm, ior(det_p, mask))
              if (row == 0_i8) cycle
              phase = pop_below(det, n + q) + pop_below(det_q, q) &
                      + pop_below(det_aq, p) + pop_below(det_p, n + p)
              if (mod(phase, 2) == 0) then
                flip = flip + cvec((row - 1_i8) * int(nv, i8) &
                                   + int(root, i8)) * c_col
              else
                flip = flip - cvec((row - 1_i8) * int(nv, i8) &
                                   + int(root, i8)) * c_col
              end if
            end do
          end do
        end do
      end if
      s2(root) = base + flip
    end do
    !$omp end parallel do

    deallocate(keys, skeys, perm, norms, cvec)
  end function fci_spin_square

end module fci_setup_mod
