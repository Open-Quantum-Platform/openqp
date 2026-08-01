!> @brief Fortran engine for the spin-orbital reduced density matrices of
!> pyoqp rdm.py.
!>
!> The determinant RDMs feed CASCI, CASSCF (generalized Fock / orbital
!> gradient), CASPT2 and NEVPT2, so they sit on the hot path of the whole
!> native wavefunction stack.  Two bind(C) entry points replace the Python
!> enumeration (the Python implementations remain as the numerical pin):
!>
!>  * rdm1_spinorb -- D1[p,q] = <a+_p a_q>
!>  * rdm2_spinorb -- D2[p,q,r,s] = <a+_p a+_q a_s a_r>
!>
!> D2 is built by factorising through the doubly-annihilated intermediate
!> space rather than walking creations back onto every determinant.  Writing
!>
!>     X[(a,b), t] = sum_d c_d <t| a_b a_a |d>
!>
!> the bra pair (p,q) and the ket pair (r,s) meet on that shared intermediate,
!> so the whole 2-RDM is one Gram matrix D2 = X X^T (real CI vectors), i.e. a
!> single DGEMM.  This mirrors the validated Python path exactly.
!>
!> The intermediates are collected by sorting the (nelec-2)-electron keys that
!> the double annihilation actually reaches and binary-searching them, the
!> same idiom fci_hamiltonian.F90 uses for its determinant lookups -- no hash
!> map, and the column order is deterministic (ascending key), so the DGEMM
!> reduction order is reproducible run to run.
!>
!> Determinants are the fci.py integer keys (alpha bits 0..norb-1, beta bits
!> norb..2norb-1; requires nspin <= 63).  Outputs are C-order: D1 is
!> [nspin,nspin] and D2 is [nspin,nspin,nspin,nspin] with the last index
!> fastest, matching the numpy arrays the caller allocates.
module rdm_kernel_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  public :: rdm1_spinorb, rdm2_spinorb, nevpt2_make_rdms

contains

  !> Population count of the low bits of a determinant key.
  pure function popcnt_below(det, bit) result(n)
    integer(i8), intent(in) :: det, bit
    integer :: n
    n = popcnt(iand(det, bit - 1_i8))
  end function popcnt_below

  !> Ascending in-place sort of keys with the matching permutation.
  subroutine sort_keys(n, keys)
    integer(i8), intent(in) :: n
    integer(i8), intent(inout) :: keys(0:)
    integer(i8) :: gap, i, j, key
    ! Shell sort: no scratch allocation, and n here is at most ndet*nelec^2.
    gap = n / 2_i8
    do while (gap > 0_i8)
      do i = gap, n - 1_i8
        key = keys(i)
        j = i
        do while (j >= gap)
          if (keys(j - gap) <= key) exit
          keys(j) = keys(j - gap)
          j = j - gap
        end do
        keys(j) = key
      end do
      gap = gap / 2_i8
    end do
  end subroutine sort_keys

  !> Ascending sort of `keys` carrying `perm` alongside, so the caller can map
  !> a search hit back to the original (unsorted) determinant position.  The
  !> determinant list handed in by fci.py is in CI ordering, not key order, so
  !> it must not be binary-searched directly.
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

  !> Index of `key` in the ascending array `keys(0:n-1)`, or -1.
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


  !> Spin-free single-excitation operators E_pq = sum_sigma a+_{p,sigma} a_{q,sigma}
  !> for one determinant basis, as (row, col, sign) triples per (p,q).
  !>
  !> E_pq never depends on the vector it is applied to, so the norb^2 operators
  !> are built once and reused; make_rdms applies them O(norb^6) times.
  subroutine build_epq(norb, ndet, dets, skeys, sperm, cnt, rows, cols, sgns)
    ! OpenQP is ILP64: the default integer is 8 bytes, so norb must be declared
    ! with the same 4-byte kind the bind(C) caller receives it in.
    integer(c_int32_t), intent(in) :: norb
    integer(i8), intent(in) :: ndet, dets(0:), skeys(0:), sperm(0:)
    integer(i8), intent(out) :: cnt(0:), rows(0:,0:), cols(0:,0:)
    real(dp), intent(out) :: sgns(0:,0:)

    integer(i8) :: k, det, det1, qbit, pbit, idx
    integer :: p, q, pq, shift
    real(dp) :: sgn

    cnt = 0_i8
    do p = 0, norb - 1
      do q = 0, norb - 1
        pq = p*norb + q
        do k = 0_i8, ndet - 1_i8
          det = dets(k)
          do shift = 0, norb, norb
            if (shift /= 0 .and. shift /= norb) cycle
            qbit = ishft(1_i8, q + shift)
            if (iand(det, qbit) == 0_i8) cycle
            sgn = 1.0_dp
            if (mod(popcnt_below(det, qbit), 2) /= 0) sgn = -1.0_dp
            det1 = ieor(det, qbit)
            pbit = ishft(1_i8, p + shift)
            if (iand(det1, pbit) /= 0_i8) cycle
            if (mod(popcnt_below(det1, pbit), 2) /= 0) sgn = -sgn
            idx = bsearch(ndet, skeys, ior(det1, pbit))
            if (idx < 0_i8) cycle
            rows(cnt(pq), pq) = sperm(idx)
            cols(cnt(pq), pq) = k
            sgns(cnt(pq), pq) = sgn
            cnt(pq) = cnt(pq) + 1_i8
          end do
        end do
      end do
    end do
  end subroutine build_epq

  !> out = E_pq . vec for one prebuilt operator.
  subroutine apply_epq(ndet, n, rows, cols, sgns, vec, out)
    integer(i8), intent(in) :: ndet, n, rows(0:), cols(0:)
    real(dp), intent(in) :: sgns(0:), vec(0:)
    real(dp), intent(out) :: out(0:)
    integer(i8) :: k
    out(0:ndet-1) = 0.0_dp
    do k = 0_i8, n - 1_i8
      out(rows(k)) = out(rows(k)) + sgns(k) * vec(cols(k))
    end do
  end subroutine apply_epq

  !> D1[p,q] = <a+_p a_q>, C-order [nspin,nspin].
  subroutine rdm1_spinorb(nspin, ndet, dets, civec, d1) &
      bind(C, name="rdm1_spinorb")
    integer(c_int32_t), value :: nspin
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(0:ndet-1)
    real(dp), intent(in) :: civec(0:ndet-1)
    real(dp), intent(inout) :: d1(0:nspin*nspin-1)

    integer(i8) :: col, det, det_q, det_pq, row, qbit, pbit
    integer(i8), allocatable :: skeys(:), sperm(:)
    integer :: p, q, phase
    real(dp) :: c

    allocate(skeys(0:ndet-1), sperm(0:ndet-1))
    do col = 0_i8, ndet - 1_i8
      skeys(col) = dets(col)
      sperm(col) = col
    end do
    call sort_keys_perm(ndet, skeys, sperm)

    d1 = 0.0_dp
    do col = 0_i8, ndet - 1_i8
      c = civec(col)
      if (c == 0.0_dp) cycle
      det = dets(col)
      do q = 0, nspin - 1
        qbit = ishft(1_i8, q)
        if (iand(det, qbit) == 0_i8) cycle
        phase = 1
        if (mod(popcnt_below(det, qbit), 2) /= 0) phase = -1
        det_q = ieor(det, qbit)
        do p = 0, nspin - 1
          pbit = ishft(1_i8, p)
          if (iand(det_q, pbit) /= 0_i8) cycle
          det_pq = ior(det_q, pbit)
          row = bsearch(ndet, skeys, det_pq)
          if (row < 0_i8) cycle
          row = sperm(row)
          if (mod(popcnt_below(det_q, pbit), 2) /= 0) then
            d1(p*nspin + q) = d1(p*nspin + q) - phase * civec(row) * c
          else
            d1(p*nspin + q) = d1(p*nspin + q) + phase * civec(row) * c
          end if
        end do
      end do
    end do
    deallocate(skeys, sperm)
  end subroutine rdm1_spinorb

  !> D2[p,q,r,s] = <a+_p a+_q a_s a_r>, C-order [nspin,nspin,nspin,nspin].
  !>
  !> Returns 0 on success, or -1 if the caller's `cap` on the number of
  !> reachable intermediates was too small (the Python fallback then runs).
  function rdm2_spinorb(nspin, ndet, dets, civec, cap, d2, nthreads) result(info) &
      bind(C, name="rdm2_spinorb")
    integer(c_int32_t), value :: nspin, nthreads
    integer(i8), value :: ndet, cap
    integer(i8), intent(in) :: dets(0:ndet-1)
    real(dp), intent(in) :: civec(0:ndet-1)
    real(dp), intent(inout) :: d2(0:nspin*nspin*nspin*nspin-1)
    integer(i8) :: info

    integer(i8), allocatable :: keys(:)
    real(dp), allocatable :: x(:,:), gram(:,:)
    integer(i8) :: col, det, det_a, det_ab, nkey, t, idx
    integer(i8) :: abit, bbit
    integer :: a, b, npair, pq, rs, p, q, r, s, phase
    real(dp) :: c

    info = 0_i8
    npair = nspin * nspin

    ! Pass 1: collect every reachable doubly-annihilated key.
    allocate(keys(0:cap-1))
    nkey = 0_i8
    do col = 0_i8, ndet - 1_i8
      if (civec(col) == 0.0_dp) cycle
      det = dets(col)
      do a = 0, nspin - 1
        abit = ishft(1_i8, a)
        if (iand(det, abit) == 0_i8) cycle
        det_a = ieor(det, abit)
        do b = 0, nspin - 1
          bbit = ishft(1_i8, b)
          if (iand(det_a, bbit) == 0_i8) cycle
          det_ab = ieor(det_a, bbit)
          if (nkey >= cap) then
            deallocate(keys)
            info = -1_i8
            return
          end if
          keys(nkey) = det_ab
          nkey = nkey + 1_i8
        end do
      end do
    end do

    if (nkey == 0_i8) then
      d2 = 0.0_dp
      deallocate(keys)
      return
    end if

    ! Sort and compact to unique ascending intermediates.
    call sort_keys(nkey, keys)
    idx = 0_i8
    do t = 1_i8, nkey - 1_i8
      if (keys(t) /= keys(idx)) then
        idx = idx + 1_i8
        keys(idx) = keys(t)
      end if
    end do
    nkey = idx + 1_i8

    ! Pass 2: accumulate X[(a,b), t].
    allocate(x(0:nkey-1, 0:npair-1))
    x = 0.0_dp
    do col = 0_i8, ndet - 1_i8
      c = civec(col)
      if (c == 0.0_dp) cycle
      det = dets(col)
      do a = 0, nspin - 1
        abit = ishft(1_i8, a)
        if (iand(det, abit) == 0_i8) cycle
        phase = 1
        if (mod(popcnt_below(det, abit), 2) /= 0) phase = -1
        det_a = ieor(det, abit)
        do b = 0, nspin - 1
          bbit = ishft(1_i8, b)
          if (iand(det_a, bbit) == 0_i8) cycle
          det_ab = ieor(det_a, bbit)
          t = bsearch(nkey, keys, det_ab)
          if (t < 0_i8) cycle
          if (mod(popcnt_below(det_a, bbit), 2) /= 0) then
            x(t, a*nspin + b) = x(t, a*nspin + b) - phase * c
          else
            x(t, a*nspin + b) = x(t, a*nspin + b) + phase * c
          end if
        end do
      end do
    end do
    deallocate(keys)

    ! D2 = X^T X over the intermediate index: one DGEMM.
    allocate(gram(0:npair-1, 0:npair-1))
    call dgemm('T', 'N', npair, npair, int(nkey), 1.0_dp, x, int(nkey), &
               x, int(nkey), 0.0_dp, gram, npair)
    deallocate(x)

    ! gram(pq, rs) with pq = p*nspin+q is exactly D2[p,q,r,s]; write it out in
    ! the caller's C order (last index fastest).
    do p = 0, nspin - 1
      do q = 0, nspin - 1
        pq = p*nspin + q
        do r = 0, nspin - 1
          do s = 0, nspin - 1
            rs = r*nspin + s
            d2(((p*nspin + q)*nspin + r)*nspin + s) = gram(pq, rs)
          end do
        end do
      end do
    end do
    deallocate(gram)
  end function rdm2_spinorb


  !> Spin-free dm1..dm4 in the PySCF make_dm1234 convention
  !> (dm_n[p,q,r,s,...] = <0| E_pq E_rs ... |0>), from the active-space
  !> determinant CI vector.  Mirrors nevpt2_sc.make_rdms exactly.
  !>
  !> `upto` selects how many are produced; the caller must have allocated every
  !> array it requests.  Returns 0 on success, -1 if allocation failed.
  function nevpt2_make_rdms(norb, ndet, dets, civec, upto, dm1, dm2, dm3, dm4) &
      result(info) bind(C, name="nevpt2_make_rdms")
    integer(c_int32_t), value :: norb, upto
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(0:ndet-1)
    real(dp), intent(in) :: civec(0:ndet-1)
    ! bind(C) receives bare C pointers, so these must be assumed-SIZE, not
    ! assumed-shape -- an assumed-shape dummy expects a descriptor and would
    ! read garbage bounds off the stack.
    real(dp), intent(inout) :: dm1(0:*), dm2(0:*), dm3(0:*), dm4(0:*)
    integer(i8) :: info

    integer(i8), allocatable :: skeys(:), sperm(:), cnt(:), rows(:,:), cols(:,:)
    real(dp), allocatable :: sgns(:,:), e1(:,:), bras(:,:), e2(:,:), kets(:,:), blk(:,:)
    integer(i8) :: k, maxnnz
    integer :: p, q, r, s, t, u, v, w, npair, pq, rs, tu, vw, tuvw, ierr

    info = 0_i8
    npair = norb * norb
    maxnnz = 2_i8 * ndet

    allocate(skeys(0:ndet-1), sperm(0:ndet-1), stat=ierr)
    if (ierr /= 0) then; info = -1_i8; return; end if
    do k = 0_i8, ndet - 1_i8
      skeys(k) = dets(k)
      sperm(k) = k
    end do
    call sort_keys_perm(ndet, skeys, sperm)

    allocate(cnt(0:npair-1), rows(0:maxnnz-1, 0:npair-1), &
             cols(0:maxnnz-1, 0:npair-1), sgns(0:maxnnz-1, 0:npair-1), stat=ierr)
    if (ierr /= 0) then; info = -1_i8; return; end if
    call build_epq(norb, ndet, dets, skeys, sperm, cnt, rows, cols, sgns)
    deallocate(skeys, sperm)

    ! E1(:, pq) = E_pq . ci
    allocate(e1(0:ndet-1, 0:npair-1), stat=ierr)
    if (ierr /= 0) then; info = -1_i8; return; end if
    do pq = 0, npair - 1
      call apply_epq(ndet, cnt(pq), rows(:,pq), cols(:,pq), sgns(:,pq), civec, e1(:,pq))
    end do

    do p = 0, norb - 1
      do q = 0, norb - 1
        dm1(p*norb + q) = dot_product(civec, e1(:, p*norb + q))
      end do
    end do
    if (upto == 1) return

    ! bras(:, p*norb+q) = E1[(q,p)], so every later GEMM is already in
    ! [p*norb+q, r*norb+s] order.
    allocate(bras(0:ndet-1, 0:npair-1), blk(0:npair-1, 0:npair-1), stat=ierr)
    if (ierr /= 0) then; info = -1_i8; return; end if
    do p = 0, norb - 1
      do q = 0, norb - 1
        bras(:, p*norb + q) = e1(:, q*norb + p)
      end do
    end do

    call dgemm('T', 'N', npair, npair, int(ndet), 1.0_dp, bras, int(ndet), &
               e1, int(ndet), 0.0_dp, blk, npair)
    do pq = 0, npair - 1
      do rs = 0, npair - 1
        dm2(int(pq, i8)*npair + rs) = blk(pq, rs)
      end do
    end do
    if (upto == 2) return

    ! E2(:, rs + npair*tu) = E_rs . E1[(t,u)]
    allocate(e2(0:ndet-1, 0:npair*npair-1), stat=ierr)
    if (ierr /= 0) then; info = -1_i8; return; end if
    do tu = 0, npair - 1
      do rs = 0, npair - 1
        call apply_epq(ndet, cnt(rs), rows(:,rs), cols(:,rs), sgns(:,rs), &
                       e1(:, tu), e2(:, rs + npair*tu))
      end do
    end do

    do tu = 0, npair - 1
      call dgemm('T', 'N', npair, npair, int(ndet), 1.0_dp, bras, int(ndet), &
                 e2(0, npair*tu), int(ndet), 0.0_dp, blk, npair)
      do pq = 0, npair - 1
        do rs = 0, npair - 1
          dm3((int(pq, i8)*npair + rs)*npair + tu) = blk(pq, rs)
        end do
      end do
    end do
    if (upto == 3) return

    ! dm4[p,q,r,s,t,u,v,w] = E1[(q,p)] . (E_rs E2[(t,u,v,w)]): one DGEMM per
    ! (t,u,v,w) block instead of norb^8 separate inner products.
    allocate(kets(0:ndet-1, 0:npair-1), stat=ierr)
    if (ierr /= 0) then; info = -1_i8; return; end if
    ! e2 stores E2[(r,s,t,u)] at column rs + npair*tu, so E2[(t,u,v,w)] lives
    ! at column tu + npair*vw -- while dm4's trailing [t,u,v,w] block wants
    ! tu*npair + vw.  The two orders are transposed; keep them separate.
    do tu = 0, npair - 1
      do vw = 0, npair - 1
        tuvw = tu + npair*vw
        do rs = 0, npair - 1
          call apply_epq(ndet, cnt(rs), rows(:,rs), cols(:,rs), sgns(:,rs), &
                         e2(:, tuvw), kets(:, rs))
        end do
        call dgemm('T', 'N', npair, npair, int(ndet), 1.0_dp, bras, int(ndet), &
                   kets, int(ndet), 0.0_dp, blk, npair)
        do pq = 0, npair - 1
          do rs = 0, npair - 1
            dm4((int(pq, i8)*npair + rs)*npair*npair + tu*npair + vw) = blk(pq, rs)
          end do
        end do
      end do
    end do
    deallocate(kets, e2, bras, blk, e1, cnt, rows, cols, sgns)
  end function nevpt2_make_rdms

end module rdm_kernel_mod
