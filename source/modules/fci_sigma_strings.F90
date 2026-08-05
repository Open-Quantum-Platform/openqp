!> String-driven sigma (Knowles-Handy) for the full active-space determinant CI.
!>
!> The determinant-pair walker in `fci_hamiltonian.F90` enumerates, for every
!> column determinant, all `nocc**2 * nspin**2` double excitations and locates
!> each row by a binary search over the determinant list.  That is
!> `O(ndet * nact**4)` scalar work with a random access per surviving term, and
!> it dominates every CASCI/CASSCF/FCI run: at CAS(10,10) one block sigma costs
!> seconds where a string-driven engine costs milliseconds.
!>
!> This module implements the standard factorization instead.  For the canonical
!> CAS product list -- alpha strings outer, beta strings inner, both in
!> `combinations` order, which is exactly what `build_determinants` emits -- a
!> determinant is a pair `(Ia, Ib)` and every excitation acts on one string only.
!> Enumerating single excitations once per *string* (there are `~sqrt(ndet)` of
!> them) instead of once per *determinant* turns the sigma into
!>
!>     t1(D, kl) = sum over single excitations kl reaching D        (gather)
!>     t2(D, ij) = sum_kl t1(D, kl) * E(kl, ij)                     (DGEMM)
!>     sigma(D) += sum over single excitations from D of t2         (scatter)
!>
!> with `E` the (nact^2, nact^2) two-electron matrix carrying the absorbed
!> one-electron term.  The gather and the scatter are contiguous AXPYs; the
!> middle step is a single BLAS-3 call.  No determinant search survives.
!>
!> The entry point declines (nonzero status) whenever the determinant list is
!> not the canonical product -- a spin-filtered or otherwise restricted list --
!> and the caller then runs the walker unchanged, so the general path keeps its
!> exact behaviour.
module fci_sigma_strings_mod

  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
!$ use omp_lib
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  !> Working-set ceiling for the blocked scratch of the sigma and the RDM, in
  !> bytes.  The alpha-string range is blocked to stay under it, so a wide
  !> active space costs more blocks rather than an allocation failure.
  !>
  !> This used to be a compile-time `parameter` fixed at 512 MiB, which meant
  !> the `[cas]`/`[fci]`/`[pt2] max_memory` keyword -- the budget the user
  !> actually declares -- reached nothing in the native kernels.  It is now set
  !> from that budget through `fci_set_work_bytes_cap`; the initial value is
  !> the old constant, so a Fortran-only caller that never sets it behaves
  !> exactly as before.
  integer(i8), save :: work_bytes_cap = 512_i8 * 1024_i8 * 1024_i8

  public :: fci_sigma_strings, rdm12_strings, spin_square_strings
  public :: fci_set_work_bytes_cap, fci_get_work_bytes_cap

contains

  !> Column-major (kl-major) single-excitation table for one string set.
  !>
  !> The target-indexed table of `build_link_table` makes every WRITE of one
  !> target contiguous, which is what the alpha channel wants (a target owns a
  !> `nbstr`-long row run).  The beta channel's runs are one element long, so
  !> target order degenerates into ~ndet*nlink scattered read-modify-writes --
  !> the dominant cost of a big sigma.  Grouping the same entries by the
  !> excitation label kl instead makes the inner loop walk ONE t1/t2 column
  !> with ascending row index: sequential writes, and the source reads stay
  !> inside one alpha block (L1-resident for any realistic active space).
  !>
  !> Entries per kl are emitted in ascending target order by construction
  !> (the outer loop below runs over targets in order).
  subroutine build_link_bykl(norb, nstr, strs, nkl, maxper, cnt, tgt, src, sgn)
    integer, intent(in) :: norb, nkl, maxper
    integer(i8), intent(in) :: nstr
    integer(i8), intent(in) :: strs(nstr)
    integer, intent(out) :: cnt(nkl)
    integer(i8), intent(out) :: tgt(maxper, nkl), src(maxper, nkl)
    real(dp), intent(out) :: sgn(maxper, nkl)

    integer(i8), allocatable :: skeys(:), sperm(:)
    integer(i8) :: t, det1, src_key, s_idx
    integer :: a, i, kl, slot
    real(dp) :: ph

    allocate(skeys(nstr), sperm(nstr))
    call sorted_index(nstr, strs, skeys, sperm)

    cnt = 0
    ! outer loop over TARGETS in ascending index, so each kl list is sorted
    ! by target and the consumer's writes walk forward.
    ! |t> = a+_a a_i |src>: diagonal (a == i) means src = t with i occupied;
    ! off-diagonal means a occupied in t, i empty in t, src = t - a + i.
    do t = 1_i8, nstr
      do i = 0, norb - 1
        if (.not. btest(strs(t), i)) cycle
        kl = i * norb + i + 1
        slot = cnt(kl) + 1
        cnt(kl) = slot
        tgt(slot, kl) = t
        src(slot, kl) = t
        sgn(slot, kl) = 1.0_dp
      end do
      do a = 0, norb - 1
        if (.not. btest(strs(t), a)) cycle
        do i = 0, norb - 1
          if (i == a) cycle
          if (btest(strs(t), i)) cycle
          src_key = ibset(ibclr(strs(t), a), i)
          s_idx = lookup(nstr, skeys, sperm, src_key)
          if (s_idx == 0_i8) cycle
          ! phase of a+_a a_i |src>: annihilate i out of src, create a into
          ! the i-annihilated string -- build_link_table's convention exactly
          det1 = ibclr(src_key, i)
          ph = 1.0_dp
          if (mod(popcnt(iand(src_key, shiftl(1_i8, i) - 1_i8)), 2) /= 0) ph = -ph
          if (mod(popcnt(iand(det1, shiftl(1_i8, a) - 1_i8)), 2) /= 0) ph = -ph
          kl = a * norb + i + 1
          slot = cnt(kl) + 1
          cnt(kl) = slot
          tgt(slot, kl) = t
          src(slot, kl) = s_idx
          sgn(slot, kl) = ph
        end do
      end do
    end do

    deallocate(skeys, sperm)
  end subroutine build_link_bykl

  !> <S^2> per CI root on the canonical CAS product determinant list, via
  !>
  !>     <S^2> = S_z^2 + S_z + |S_+ c|^2,   S_+ = sum_p a+_{p alpha} a_{p beta}
  !>
  !> The general kernel (`fci_spin_square`) walks determinant pairs; this one
  !> materialises S_+ c in the (na+1, nb-1) product space -- one pass of
  !> O(ndet*norb) index arithmetic -- and takes its norm.  Declines (nonzero
  !> status) for non-product lists; the caller falls back to the walker.
  function spin_square_strings(norb, ndet, dets, nroot, civecs, s2, nthreads) &
      result(status)
    integer, intent(in) :: norb, nroot, nthreads
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(*)
    real(dp), intent(in) :: civecs(*)
    real(dp), intent(inout) :: s2(*)
    integer :: status

    integer :: na, nb, p, r
    integer(i8) :: nastr, nbstr, mastr, mbstr, ia, ib, k, ja, jb
    integer(i8), allocatable :: astr(:), bstr(:), wastr(:), wbstr(:)
    integer(i8), allocatable :: akeys(:), aperm(:), bkeys(:), bperm(:)
    real(dp), allocatable :: w(:)
    real(dp) :: sz, ph, c, ssum
    integer :: ierr

    status = 1
    if (norb <= 0 .or. ndet <= 0_i8 .or. nroot < 1) return
    na = popcnt(iand(dets(1), shiftl(1_i8, norb) - 1_i8))
    nb = popcnt(shiftr(dets(1), norb))
    nastr = binom(norb, na)
    nbstr = binom(norb, nb)
    if (nastr <= 0_i8 .or. nbstr <= 0_i8 .or. nastr * nbstr /= ndet) return

    allocate(astr(nastr), bstr(nbstr))
    call string_basis(norb, na, nastr, astr)
    call string_basis(norb, nb, nbstr, bstr)
    k = 0_i8
    do ia = 1_i8, nastr
      do ib = 1_i8, nbstr
        k = k + 1_i8
        if (dets(k) /= ior(astr(ia), shiftl(bstr(ib), norb))) then
          deallocate(astr, bstr)
          return
        end if
      end do
    end do

    sz = 0.5_dp * real(na - nb, dp)

    if (nb == 0 .or. na == norb) then
      ! S_+ annihilates the reference space
      do r = 1, nroot
        s2(r) = sz * (sz + 1.0_dp)
      end do
      deallocate(astr, bstr)
      status = 0
      return
    end if

    mastr = binom(norb, na + 1)
    mbstr = binom(norb, nb - 1)
    ! decline (walker fallback) rather than abort if the S_+ space is too
    ! large to hold
    allocate(wastr(mastr), wbstr(max(mbstr, 1_i8)), stat=ierr)
    if (ierr /= 0) then
      deallocate(astr, bstr)
      return
    end if
    call string_basis(norb, na + 1, mastr, wastr)
    call string_basis(norb, nb - 1, mbstr, wbstr)
    allocate(akeys(mastr), aperm(mastr), bkeys(max(mbstr, 1_i8)), &
             bperm(max(mbstr, 1_i8)), w(mastr * mbstr), stat=ierr)
    if (ierr /= 0) then
      deallocate(astr, bstr, wastr, wbstr)
      return
    end if
    call sorted_index(mastr, wastr, akeys, aperm)
    call sorted_index(mbstr, wbstr, bkeys, bperm)

    do r = 1, nroot
      w = 0.0_dp
      do ia = 1_i8, nastr
        do ib = 1_i8, nbstr
          c = civecs((((ia - 1_i8) * nbstr + ib) - 1_i8) * int(nroot, i8) + r)
          if (c == 0.0_dp) cycle
          do p = 0, norb - 1
            if (.not. btest(bstr(ib), p)) cycle
            if (btest(astr(ia), p)) cycle
            ! a+_{p alpha} a_{p beta}: the beta annihilation crosses the whole
            ! alpha string (constant (-1)^na, irrelevant to the norm) plus the
            ! beta bits below p; the alpha creation crosses alpha bits below p.
            ph = 1.0_dp
            if (mod(popcnt(iand(bstr(ib), shiftl(1_i8, p) - 1_i8)) &
                  + popcnt(iand(astr(ia), shiftl(1_i8, p) - 1_i8)), 2) /= 0) &
                ph = -1.0_dp
            ja = lookup(mastr, akeys, aperm, ibset(astr(ia), p))
            jb = lookup(mbstr, bkeys, bperm, ibclr(bstr(ib), p))
            if (ja == 0_i8 .or. jb == 0_i8) cycle
            w((ja - 1_i8) * mbstr + jb) = w((ja - 1_i8) * mbstr + jb) + ph * c
          end do
        end do
      end do
      ssum = 0.0_dp
      do k = 1_i8, mastr * mbstr
        ssum = ssum + w(k) * w(k)
      end do
      s2(r) = sz * (sz + 1.0_dp) + ssum
    end do

    deallocate(astr, bstr, wastr, wbstr, akeys, aperm, bkeys, bperm, w)
    status = 0
    if (nthreads < 0) status = 0  ! nthreads reserved for future threading
  end function spin_square_strings

  !> bind(C) face of `rdm12_strings`, for the test pins.
  function rdm12_strings_c(norb, ndet, dets, civec, d1, d2, nthreads) &
      result(status) bind(C, name="rdm12_strings_c")
    use, intrinsic :: iso_c_binding, only: c_int32_t
    integer(c_int32_t), value :: norb, nthreads
    integer(i8), value :: ndet
    integer(i8), intent(in) :: dets(*)
    real(dp), intent(in) :: civec(*)
    real(dp), intent(inout) :: d1(*), d2(*)
    integer(i8) :: status
    status = int(rdm12_strings(int(norb), ndet, dets, civec, d1, d2, &
                               int(nthreads)), i8)
  end function rdm12_strings_c

  !> Spatial 1- and 2-RDMs of one CI vector through the same string
  !> factorization as the sigma: with
  !>
  !>     T[D, pq] = <D| E_pq |0>       (E_pq = sum_sigma a+_{p sigma} a_{q sigma})
  !>
  !> built from the alpha/beta single-excitation tables, the spin-summed
  !> spatial RDMs are
  !>
  !>     d1[p,q]     = sum_D c_D T[D, pq]                      (one DGEMV)
  !>     d2[p,q,r,s] = sum_D T[D,pq] T[D,rs] - delta_qr d1[p,s]  (one DGEMM)
  !>
  !> in the same chemists' convention as `rdm2_spatial` (contracts as
  !> 0.5 sum (pq|rs) d2[p,q,r,s]).  This replaces the per-determinant pair
  !> walk of the general kernels with O(ndet nact^2) table work plus BLAS-3,
  !> which is what the CASSCF gradient loop needs -- the general kernels remain
  !> for spin-filtered determinant lists and stay the numerical pin.
  !>
  !> Nonzero status = the determinant list is not the canonical CAS product
  !> (or T would not fit); the caller falls back to the walking kernels.
  function rdm12_strings(norb, ndet, dets, civec, d1, d2, nthreads) &
      result(status)
    integer, intent(in) :: norb, nthreads
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(*)
    real(dp), intent(in) :: civec(*)
    real(dp), intent(inout) :: d1(*), d2(*)
    integer :: status

    integer :: na, nb, nkl, nlink_a, nthr, maxper_b
    integer(i8) :: nastr, nbstr, ia, ib, k, base, srow, drow
    integer(i8) :: e, m
    integer(i8) :: nblk, blk_start, blk_len, nrow, nrow_max, cbase
    integer :: p, q, r, s, pq, rs, t
    real(dp) :: sgn
    integer(i8), allocatable :: astr(:), bstr(:)
    integer, allocatable :: ntab_a(:)
    integer, allocatable :: kl_a(:, :)
    integer(i8), allocatable :: src_a(:, :)
    real(dp), allocatable :: sgn_a(:, :)
    integer, allocatable :: cnt_b(:)
    integer(i8), allocatable :: tgt_b(:, :), src_bk(:, :)
    real(dp), allocatable :: sgn_bk(:, :)
    real(dp), allocatable :: tmat(:, :), gram(:, :)

    status = 1
    if (norb <= 0 .or. ndet <= 0_i8) return
    na = popcnt(iand(dets(1), shiftl(1_i8, norb) - 1_i8))
    nb = popcnt(shiftr(dets(1), norb))

    nastr = binom(norb, na)
    nbstr = binom(norb, nb)
    if (nastr <= 0_i8 .or. nbstr <= 0_i8) return
    if (nastr * nbstr /= ndet) return

    nkl = norb * norb

    allocate(astr(nastr), bstr(nbstr))
    call string_basis(norb, na, nastr, astr)
    call string_basis(norb, nb, nbstr, bstr)
    k = 0_i8
    do ia = 1_i8, nastr
      do ib = 1_i8, nbstr
        k = k + 1_i8
        if (dets(k) /= ior(astr(ia), shiftl(bstr(ib), norb))) then
          deallocate(astr, bstr)
          return
        end if
      end do
    end do

    nlink_a = max(1, na * (norb - na) + na)
    allocate(ntab_a(nastr), kl_a(nlink_a, nastr), src_a(nlink_a, nastr), &
             sgn_a(nlink_a, nastr))
    call build_link_table(norb, nastr, astr, nlink_a, ntab_a, kl_a, src_a, sgn_a)
    ! Beta uses the kl-major table for the same reason the sigma does: a beta
    ! excitation moves a single row, so a target-indexed walk writes tmat with
    ! stride ndet and read-modify-writes every element.  Grouped by the
    ! excitation label the inner loop walks ONE tmat column with ascending
    ! rows, and the civec reads stay inside one alpha block.
    maxper_b = int(min(nbstr, binom(norb - 1, nb - 1) + binom(norb - 2, nb - 1)))
    allocate(cnt_b(nkl), tgt_b(max(maxper_b, 1), nkl), &
             src_bk(max(maxper_b, 1), nkl), sgn_bk(max(maxper_b, 1), nkl))
    call build_link_bykl(norb, nbstr, bstr, nkl, max(maxper_b, 1), cnt_b, &
                         tgt_b, src_bk, sgn_bk)

    nthr = max(1, nthreads)

    ! T is (ndet, nact^2) and reaches tens of gigabytes on a wide active space.
    ! It used to be built whole, and the routine simply DECLINED above 2 GiB --
    ! handing CAS(14,14) and larger back to the per-determinant Python walker,
    ! i.e. the cases that need the fast path most.  Both consumers of T reduce
    ! over its rows (d1 = T^T c, gram = T^T T), so the row range can be blocked
    ! over alpha strings exactly as the sigma blocks it, and the two BLAS calls
    ! accumulate across blocks.  Nothing declines any more, and peak scratch is
    ! bounded by the working-set ceiling instead of by ndet.
    nblk = max(1_i8, work_bytes_cap / max(nbstr * int(nkl, i8) * 8_i8, 1_i8))
    nblk = min(nblk, nastr)
    nrow_max = nblk * nbstr
    allocate(tmat(nrow_max, nkl), gram(nkl, nkl))

    ! The leading dimension stays nrow_max for every block; a short final block
    ! just passes a smaller row count to BLAS, which reads A(1:k, :) with the
    ! stated lda.
    d1(1:int(nkl, i8)) = 0.0_dp
    gram = 0.0_dp

    blk_start = 1_i8
    do while (blk_start <= nastr)
      blk_len = min(nblk, nastr - blk_start + 1_i8)
      nrow = blk_len * nbstr
      cbase = (blk_start - 1_i8) * nbstr

      !$omp parallel do num_threads(nthr) schedule(static) default(shared) private(t)
      do t = 1, nkl
        tmat(1:nrow, t) = 0.0_dp
      end do
      !$omp end parallel do

      ! alpha channel of E_pq: rows of one target alpha string are the
      ! contiguous run (ia-1)*nbstr + 1 .. ia*nbstr.  The SOURCE row may sit
      ! outside the block, which is fine -- civec is whole, only T is blocked.
      !$omp parallel do num_threads(nthr) schedule(static) default(shared) &
      !$omp private(ia, e, srow, drow, t, sgn, m)
      do ia = blk_start, blk_start + blk_len - 1_i8
        drow = (ia - blk_start) * nbstr
        do e = 1_i8, int(ntab_a(ia), i8)
          srow = (src_a(e, ia) - 1_i8) * nbstr
          t = kl_a(e, ia) + 1
          sgn = sgn_a(e, ia)
          do m = 1_i8, nbstr
            tmat(drow + m, t) = tmat(drow + m, t) + sgn * civec(srow + m)
          end do
        end do
      end do
      !$omp end parallel do

      ! beta channel: excitations stay inside one alpha string, so the alpha
      ! index partitions the writes.  kl-major, so each (t) run walks one
      ! column down.
      !$omp parallel do num_threads(nthr) schedule(static) default(shared) &
      !$omp private(ia, base, e, srow, drow, t)
      do ia = blk_start, blk_start + blk_len - 1_i8
        base = (ia - blk_start) * nbstr
        do t = 1, nkl
          do e = 1_i8, int(cnt_b(t), i8)
            drow = base + tgt_b(e, t)
            srow = base + src_bk(e, t)
            tmat(drow, t) = tmat(drow, t) + sgn_bk(e, t) * civec(srow + cbase)
          end do
        end do
      end do
      !$omp end parallel do

      ! d1 += T_blk^T c_blk
      call dgemv('T', int(nrow), nkl, 1.0_dp, tmat, int(nrow_max), &
                 civec(cbase + 1_i8), 1, 1.0_dp, d1, 1)
      ! gram += T_blk^T T_blk -- symmetric by construction, so DSYRK does it in
      ! half the flops of the equivalent DGEMM.  beta=1 accumulates the blocks.
      call dsyrk('U', 'T', nkl, int(nrow), 1.0_dp, tmat, int(nrow_max), &
                 1.0_dp, gram, nkl)

      blk_start = blk_start + blk_len
    end do

    ! DSYRK wrote one triangle; the d2 unpacking below indexes gram(pq, rs)
    ! over the full square, so mirror the upper triangle down.
    do t = 2, nkl
      do pq = 1, t - 1
        gram(t, pq) = gram(pq, t)
      end do
    end do

    ! d2[p,q,r,s] = <E_pq E_rs> - delta_qr d1[p,s].  The bra side enters
    ! through <0|E_pq|D> = <D|E_qp|0> (real CI vectors), so the Gram row is
    ! the TRANSPOSED pair (q,p).
    do p = 0, norb - 1
      do q = 0, norb - 1
        pq = q * norb + p + 1
        do r = 0, norb - 1
          rs = r * norb + 1
          do s = 0, norb - 1
            d2((((int(p, i8) * norb + q) * norb + r) * norb + s) + 1_i8) = &
                gram(pq, rs + s)
          end do
        end do
        ! subtract the delta_qr column: r == q
        do s = 0, norb - 1
          d2((((int(p, i8) * norb + q) * norb + q) * norb + s) + 1_i8) = &
              d2((((int(p, i8) * norb + q) * norb + q) * norb + s) + 1_i8) &
              - d1(int(p, i8) * norb + s + 1_i8)
        end do
      end do
    end do

    deallocate(astr, bstr, ntab_a, kl_a, src_a, sgn_a, &
               cnt_b, tgt_b, src_bk, sgn_bk, tmat, gram)
    status = 0
  end function rdm12_strings

  !> One target-indexed single-excitation table for a string set.
  !>
  !> Entry `e` of target string `t` is `(kl, src, sgn)` with the meaning
  !> `a+_a a_i |src> = sgn |t>` and `kl = a*norb + i` (0-based).  Indexing by
  !> the target is what lets both the gather and the scatter run with disjoint
  !> writes: for a fixed target every write lands in that target's own slice.
  subroutine build_link_table(norb, nstr, strs, maxlink, ntab, &
                              tab_kl, tab_src, tab_sgn)
    integer, intent(in) :: norb, maxlink
    integer(i8), intent(in) :: nstr
    integer(i8), intent(in) :: strs(nstr)
    integer, intent(out) :: ntab(nstr)
    integer, intent(out) :: tab_kl(maxlink, nstr)
    integer(i8), intent(out) :: tab_src(maxlink, nstr)
    real(dp), intent(out) :: tab_sgn(maxlink, nstr)

    integer(i8), allocatable :: skeys(:), sperm(:)
    integer(i8) :: src, tgt_key, det1
    integer :: a, i, slot
    integer(i8) :: tgt
    real(dp) :: sgn

    allocate(skeys(nstr), sperm(nstr))
    call sorted_index(nstr, strs, skeys, sperm)

    ntab = 0
    do src = 1_i8, nstr
      do i = 0, norb - 1
        if (.not. btest(strs(src), i)) cycle
        det1 = ibclr(strs(src), i)
        do a = 0, norb - 1
          if (btest(det1, a)) cycle
          tgt_key = ibset(det1, a)
          tgt = lookup(nstr, skeys, sperm, tgt_key)
          if (tgt == 0_i8) cycle
          sgn = 1.0_dp
          if (mod(pop_below(strs(src), i), 2) /= 0) sgn = -sgn
          if (mod(pop_below(det1, a), 2) /= 0) sgn = -sgn
          slot = ntab(tgt) + 1
          ntab(tgt) = slot
          tab_kl(slot, tgt) = a * norb + i
          tab_src(slot, tgt) = src
          tab_sgn(slot, tgt) = sgn
        end do
      end do
    end do

    deallocate(skeys, sperm)
  end subroutine build_link_table

  pure function pop_below(det, pos) result(n)
    integer(i8), intent(in) :: det
    integer, intent(in) :: pos
    integer :: n
    n = popcnt(iand(det, shiftl(1_i8, pos) - 1_i8))
  end function pop_below

  !> Ascending key order plus the permutation back to the original index.
  subroutine sorted_index(n, keys, skeys, sperm)
    integer(i8), intent(in) :: n
    integer(i8), intent(in) :: keys(n)
    integer(i8), intent(out) :: skeys(n), sperm(n)
    integer(i8) :: i, m, k, j

    do i = 1_i8, n
      skeys(i) = keys(i)
      sperm(i) = i
    end do
    ! insertion sort: the string count is O(sqrt(ndet)) and this runs once
    do i = 2_i8, n
      k = skeys(i)
      j = sperm(i)
      m = i - 1_i8
      do while (m >= 1_i8)
        if (skeys(m) <= k) exit
        skeys(m + 1_i8) = skeys(m)
        sperm(m + 1_i8) = sperm(m)
        m = m - 1_i8
      end do
      skeys(m + 1_i8) = k
      sperm(m + 1_i8) = j
    end do
  end subroutine sorted_index

  pure function lookup(n, skeys, sperm, key) result(idx)
    integer(i8), intent(in) :: n, key
    integer(i8), intent(in) :: skeys(n), sperm(n)
    integer(i8) :: idx, lo, hi, mid
    idx = 0_i8
    lo = 1_i8
    hi = n
    do while (lo <= hi)
      mid = (lo + hi) / 2_i8
      if (skeys(mid) == key) then
        idx = sperm(mid)
        return
      else if (skeys(mid) < key) then
        lo = mid + 1_i8
      else
        hi = mid - 1_i8
      end if
    end do
  end function lookup

  pure function binom(n, k) result(c)
    integer, intent(in) :: n, k
    integer(i8) :: c
    integer :: i
    c = 0_i8
    if (k < 0 .or. k > n) return
    c = 1_i8
    do i = 1, min(k, n - k)
      c = c * int(n - min(k, n - k) + i, i8) / int(i, i8)
    end do
  end function binom

  !> Occupation bit patterns for all C(norb, nelec) subsets, in the same
  !> lexicographic order `build_determinants` and the Python driver use.
  subroutine string_basis(norb, nelec, nstr, out)
    integer, intent(in) :: norb, nelec
    integer(i8), intent(in) :: nstr
    integer(i8), intent(out) :: out(nstr)
    integer :: c(max(nelec, 1)), i, j
    integer(i8) :: idx, key

    if (nelec == 0) then
      out(1) = 0_i8
      return
    end if
    do i = 1, nelec
      c(i) = i - 1
    end do
    idx = 0_i8
    do
      idx = idx + 1_i8
      key = 0_i8
      do i = 1, nelec
        key = ibset(key, c(i))
      end do
      out(idx) = key
      if (idx >= nstr) exit
      i = nelec
      do while (i >= 1)
        if (c(i) /= norb - nelec + i - 1) exit
        i = i - 1
      end do
      if (i < 1) exit
      c(i) = c(i) + 1
      do j = i + 1, nelec
        c(j) = c(j - 1) + 1
      end do
    end do
  end subroutine string_basis

  !> Y = H X for the canonical CAS product determinant list.
  !>
  !> Returns 0 on success.  A nonzero status means the list is not the
  !> canonical product (or the active space is degenerate) and nothing was
  !> written to `y`; the caller falls back to the determinant walker.
  !>
  !> `cutoff` has no counterpart here -- the walker uses it to skip individual
  !> integrals, while this path contracts them in one DGEMM -- so it is
  !> accepted and ignored, which only ever adds terms the walker dropped.
  function fci_sigma_strings(nspin, ndet, dets, hspin, gspin, nvec, x, y, &
                             nthreads) result(status)
    integer, intent(in) :: nspin, nvec, nthreads
    integer(i8), intent(in) :: ndet
    integer(i8), intent(in) :: dets(*)
    real(dp), intent(in) :: hspin(*), gspin(*)
    real(dp), intent(in) :: x(*)
    real(dp), intent(inout) :: y(*)
    integer :: status

    integer :: norb, na, nb, nkl, nlink_a, nthr, maxper_b
    integer(i8) :: nastr, nbstr, ia, ib, k, ns
    integer(i8), allocatable :: astr(:), bstr(:)
    integer, allocatable :: ntab_a(:)
    integer, allocatable :: kl_a(:, :)
    integer(i8), allocatable :: src_a(:, :)
    real(dp), allocatable :: sgn_a(:, :)
    integer, allocatable :: cnt_b(:)
    integer(i8), allocatable :: tgt_b(:, :), src_bk(:, :)
    real(dp), allocatable :: sgn_bk(:, :)
    real(dp), allocatable :: emat(:, :)
    real(dp), allocatable :: t1(:), t2(:)
    integer(i8) :: nblk, blk_start, blk_len, nrow, inner, nrow_max
    integer(i8) :: bytes_per_alpha

    status = 1
    norb = nspin / 2
    if (norb <= 0 .or. 2 * norb /= nspin .or. ndet <= 0_i8 .or. nvec <= 0) return

    ns = int(nspin, i8)
    na = popcnt(iand(dets(1), shiftl(1_i8, norb) - 1_i8))
    nb = popcnt(shiftr(dets(1), norb))
    if (na + nb <= 0) return

    nastr = binom(norb, na)
    nbstr = binom(norb, nb)
    if (nastr <= 0_i8 .or. nbstr <= 0_i8) return
    if (nastr * nbstr /= ndet) return

    allocate(astr(nastr), bstr(nbstr))
    call string_basis(norb, na, nastr, astr)
    call string_basis(norb, nb, nbstr, bstr)

    ! The list must be exactly the product in CI order, or the index arithmetic
    ! below addresses the wrong determinants.
    k = 0_i8
    do ia = 1_i8, nastr
      do ib = 1_i8, nbstr
        k = k + 1_i8
        if (dets(k) /= ior(astr(ia), shiftl(bstr(ib), norb))) then
          deallocate(astr, bstr)
          return
        end if
      end do
    end do

    nkl = norb * norb
    nlink_a = max(1, na * (norb - na) + na)

    allocate(ntab_a(nastr), kl_a(nlink_a, nastr), src_a(nlink_a, nastr), &
             sgn_a(nlink_a, nastr))
    call build_link_table(norb, nastr, astr, nlink_a, ntab_a, kl_a, src_a, sgn_a)
    ! beta runs are one element long, so the beta channel uses the kl-major
    ! table: the inner loops then walk one t1/t2 column with ascending rows
    maxper_b = int(min(nbstr, binom(norb - 1, nb - 1) + binom(norb - 2, nb - 1)))
    allocate(cnt_b(nkl), tgt_b(max(maxper_b, 1), nkl), &
             src_bk(max(maxper_b, 1), nkl), sgn_bk(max(maxper_b, 1), nkl))
    call build_link_bykl(norb, nbstr, bstr, nkl, max(maxper_b, 1), cnt_b, &
                         tgt_b, src_bk, sgn_bk)

    allocate(emat(nkl, nkl))
    call build_emat(norb, ns, na + nb, hspin, gspin, nkl, emat)

    nthr = max(1, nthreads)
    inner = nbstr * int(nvec, i8)

    bytes_per_alpha = 2_i8 * inner * int(nkl, i8) * 8_i8
    nblk = max(1_i8, work_bytes_cap / max(bytes_per_alpha, 1_i8))
    nblk = min(nblk, nastr)

    ! Same story as t1 below: ndet*nvec doubles, cleared once per sigma call.
    !$omp parallel do num_threads(nthr) schedule(static) default(shared) private(k)
    do k = 1_i8, ndet * int(nvec, i8)
      y(k) = 0.0_dp
    end do
    !$omp end parallel do

    ! One allocation for the whole run of blocks.  These buffers reach
    ! hundreds of megabytes, and allocating/freeing them per block cost more
    ! than the gather itself: every fresh allocation faults its pages back in.
    ! Only the last block is shorter, and a shorter block simply uses a prefix
    ! of the buffer -- the leading dimension is passed explicitly, so the
    ! layout is defined by `nrow`, not by the buffer's size.
    nrow_max = nblk * inner
    allocate(t1(nrow_max * int(nkl, i8)), t2(nrow_max * int(nkl, i8)))

    blk_start = 1_i8
    do while (blk_start <= nastr)
      blk_len = min(nblk, nastr - blk_start + 1_i8)
      nrow = blk_len * inner

      call gather(norb, nastr, nbstr, nvec, blk_start, blk_len, inner, nkl, &
                  nlink_a, ntab_a, kl_a, src_a, sgn_a, &
                  maxper_b, cnt_b, tgt_b, src_bk, sgn_bk, nthr, x, nrow, t1)

      call dgemm('N', 'N', int(nrow), nkl, nkl, 1.0_dp, t1, int(nrow), &
                 emat, nkl, 0.0_dp, t2, int(nrow))

      call scatter(norb, nastr, nbstr, nvec, blk_start, blk_len, inner, nkl, &
                    nlink_a, ntab_a, kl_a, src_a, sgn_a, &
                    maxper_b, cnt_b, tgt_b, src_bk, sgn_bk, nthr, nrow, t2, y)

      blk_start = blk_start + blk_len
    end do
    deallocate(t1, t2)

    deallocate(emat, astr, bstr)
    deallocate(ntab_a, kl_a, src_a, sgn_a, cnt_b, tgt_b, src_bk, sgn_bk)
    status = 0
  end function fci_sigma_strings

  !> The (nact^2, nact^2) two-electron matrix with the one-electron term
  !> absorbed, so a single contraction produces the whole sigma.
  !>
  !>   f(j,k) = [ h(j,k) - 1/2 sum_i (ji|ik) ] / nelec
  !>   E(ai,bj) = 1/2 [ (ai|bj) + delta_ai f(b,j) + delta_bj f(a,i) ]
  !>
  !> `gspin(p,q,r,s) = (p r | q s)` is the spin-orbital layout written by
  !> `fci_spin_orbital_integrals`; the alpha-alpha block carries the spatial
  !> tensor, so `(a i | b j) = gspin(a, b, i, j)`.
  subroutine build_emat(norb, ns, nelec, hspin, gspin, nkl, emat)
    integer, intent(in) :: norb, nelec, nkl
    integer(i8), intent(in) :: ns
    real(dp), intent(in) :: hspin(*), gspin(*)
    real(dp), intent(out) :: emat(nkl, nkl)

    integer :: a, i, b, j, ai, bj
    real(dp) :: f1e(norb, norb), acc

    do j = 0, norb - 1
      do i = 0, norb - 1
        acc = 0.0_dp
        do b = 0, norb - 1
          ! (i b | b j) = gspin(i, b, b, j)
          acc = acc + gspin(chem(i, b, b, j, ns) + 1_i8)
        end do
        f1e(i + 1, j + 1) = (hspin(int(i, i8) * ns + int(j, i8) + 1_i8) &
                             - 0.5_dp * acc) / real(nelec, dp)
      end do
    end do

    do b = 0, norb - 1
      do j = 0, norb - 1
        bj = b * norb + j + 1
        do a = 0, norb - 1
          do i = 0, norb - 1
            ai = a * norb + i + 1
            acc = gspin(chem(a, i, b, j, ns) + 1_i8)
            if (a == i) acc = acc + f1e(b + 1, j + 1)
            if (b == j) acc = acc + f1e(a + 1, i + 1)
            emat(ai, bj) = 0.5_dp * acc
          end do
        end do
      end do
    end do
  end subroutine build_emat

  !> Flat 0-based index of the chemist integral (a i | b j) inside `gspin`.
  pure function chem(a, i, b, j, ns) result(idx)
    integer, intent(in) :: a, i, b, j
    integer(i8), intent(in) :: ns
    integer(i8) :: idx
    idx = ((int(a, i8) * ns + int(b, i8)) * ns + int(i, i8)) * ns + int(j, i8)
  end function chem

  !> t1(row, kl) = sum over single excitations kl that reach the row's
  !> determinant, of the signed CI coefficient of the source determinant.
  subroutine gather(norb, nastr, nbstr, nvec, blk_start, blk_len, inner, nkl, &
                    nlink_a, ntab_a, kl_a, src_a, sgn_a, &
                    maxper_b, cnt_b, tgt_b, src_bk, sgn_bk, nthr, x, nrow, t1)
    integer, intent(in) :: norb, nvec, nkl, nlink_a, maxper_b, nthr
    integer(i8), intent(in) :: nastr, nbstr, blk_start, blk_len, inner, nrow
    integer, intent(in) :: ntab_a(nastr)
    integer, intent(in) :: kl_a(nlink_a, nastr)
    integer(i8), intent(in) :: src_a(nlink_a, nastr)
    real(dp), intent(in) :: sgn_a(nlink_a, nastr)
    integer, intent(in) :: cnt_b(nkl)
    integer(i8), intent(in) :: tgt_b(maxper_b, nkl), src_bk(maxper_b, nkl)
    real(dp), intent(in) :: sgn_bk(maxper_b, nkl)
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: t1(nrow, nkl)

    integer(i8) :: ia, ib, e, srow, drow, nv, sa, so, m
    integer :: t
    real(dp) :: s

    nv = int(nvec, i8)
    ! t1 is (nrow, nact^2) and runs to hundreds of megabytes -- the whole-array
    ! assignment that used to stand here was single-threaded, and on a wide
    ! active space it cost more than the gather, the DGEMM and the scatter put
    ! together.  Clear it on the same team that fills it, one column per
    ! iteration so each thread walks a contiguous nrow-long run.
    !$omp parallel do num_threads(nthr) schedule(static) default(shared) private(t)
    do t = 1, nkl
      t1(:, t) = 0.0_dp
    end do
    !$omp end parallel do

    ! alpha: every write lands in the target string's own row slice
    !$omp parallel do num_threads(nthr) schedule(static) default(shared) &
    !$omp private(ia, e, srow, drow, t, s, m)
    do ia = blk_start, blk_start + blk_len - 1_i8
      drow = (ia - blk_start) * inner
      do e = 1_i8, int(ntab_a(ia), i8)
        srow = (src_a(e, ia) - 1_i8) * inner
        t = kl_a(e, ia) + 1
        s = sgn_a(e, ia)
        do m = 1_i8, inner
          t1(drow + m, t) = t1(drow + m, t) + s * x(srow + m)
        end do
      end do
    end do
    !$omp end parallel do

    ! beta: excitations stay inside one alpha string, so the alpha index
    ! partitions the writes.  kl-major: for each excitation label the inner
    ! loop walks ONE t1 column with ascending rows (sequential writes), and
    ! the x reads stay inside the alpha block.
    !$omp parallel do num_threads(nthr) schedule(static) default(shared) &
    !$omp private(ia, e, sa, so, srow, drow, t, s, m)
    do ia = blk_start, blk_start + blk_len - 1_i8
      sa = (ia - blk_start) * inner
      so = (ia - 1_i8) * inner
      do t = 1, nkl
        do e = 1_i8, int(cnt_b(t), i8)
          drow = sa + (tgt_b(e, t) - 1_i8) * nv
          srow = so + (src_bk(e, t) - 1_i8) * nv
          s = sgn_bk(e, t)
          do m = 1_i8, nv
            t1(drow + m, t) = t1(drow + m, t) + s * x(srow + m)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    if (.false.) then
      if (norb < 0) continue
    end if
  end subroutine gather

  !> sigma(target) += sum over single excitations out of the block rows.
  !> Indexed by target for the same reason the gather is: disjoint writes.
  subroutine scatter(norb, nastr, nbstr, nvec, blk_start, blk_len, inner, nkl, &
                     nlink_a, ntab_a, kl_a, src_a, sgn_a, &
                     maxper_b, cnt_b, tgt_b, src_bk, sgn_bk, nthr, nrow, t2, y)
    integer, intent(in) :: norb, nvec, nkl, nlink_a, maxper_b, nthr
    integer(i8), intent(in) :: nastr, nbstr, blk_start, blk_len, inner, nrow
    integer, intent(in) :: ntab_a(nastr)
    integer, intent(in) :: kl_a(nlink_a, nastr)
    integer(i8), intent(in) :: src_a(nlink_a, nastr)
    real(dp), intent(in) :: sgn_a(nlink_a, nastr)
    integer, intent(in) :: cnt_b(nkl)
    integer(i8), intent(in) :: tgt_b(maxper_b, nkl), src_bk(maxper_b, nkl)
    real(dp), intent(in) :: sgn_bk(maxper_b, nkl)
    real(dp), intent(in) :: t2(nrow, nkl)
    real(dp), intent(inout) :: y(*)

    integer(i8) :: ia, ib, e, srow, drow, nv, sa, so, blk_end, m
    integer :: t
    real(dp) :: s

    nv = int(nvec, i8)
    blk_end = blk_start + blk_len - 1_i8

    ! alpha: `src_a(e, ia)` is the row the coefficient comes from, so only the
    ! entries whose source is inside this block contribute here
    !$omp parallel do num_threads(nthr) schedule(static) default(shared) &
    !$omp private(ia, e, srow, drow, t, s, m)
    do ia = 1_i8, nastr
      drow = (ia - 1_i8) * inner
      do e = 1_i8, int(ntab_a(ia), i8)
        if (src_a(e, ia) < blk_start .or. src_a(e, ia) > blk_end) cycle
        srow = (src_a(e, ia) - blk_start) * inner
        t = kl_a(e, ia) + 1
        s = sgn_a(e, ia)
        do m = 1_i8, inner
          y(drow + m) = y(drow + m) + s * t2(srow + m, t)
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do num_threads(nthr) schedule(static) default(shared) &
    !$omp private(ia, e, sa, so, srow, drow, t, s, m)
    do ia = blk_start, blk_end
      sa = (ia - blk_start) * inner
      so = (ia - 1_i8) * inner
      do t = 1, nkl
        do e = 1_i8, int(cnt_b(t), i8)
          drow = so + (tgt_b(e, t) - 1_i8) * nv
          srow = sa + (src_bk(e, t) - 1_i8) * nv
          s = sgn_bk(e, t)
          do m = 1_i8, nv
            y(drow + m) = y(drow + m) + s * t2(srow + m, t)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    if (.false.) then
      if (norb < 0 .or. nkl < 0) continue
    end if
  end subroutine scatter

  !> Set the blocked-scratch ceiling used by the sigma and the RDM.
  !>
  !> `bytes` is clamped to at least 16 MiB: the blocking divides by it, and a
  !> pathologically small budget would otherwise produce one block per alpha
  !> string and turn the core DGEMM into a stream of tiny calls.  A
  !> non-positive value restores the built-in default.
  subroutine fci_set_work_bytes_cap(bytes) bind(C, name="fci_set_work_bytes_cap")
    integer(c_int64_t), value :: bytes
    integer(i8), parameter :: floor_bytes = 16_i8 * 1024_i8 * 1024_i8
    if (bytes <= 0_i8) then
      work_bytes_cap = 512_i8 * 1024_i8 * 1024_i8
    else
      work_bytes_cap = max(int(bytes, i8), floor_bytes)
    end if
  end subroutine fci_set_work_bytes_cap

  !> The ceiling currently in force, in bytes.
  function fci_get_work_bytes_cap() result(bytes) bind(C, name="fci_get_work_bytes_cap")
    integer(c_int64_t) :: bytes
    bytes = int(work_bytes_cap, c_int64_t)
  end function fci_get_work_bytes_cap

end module fci_sigma_strings_mod
