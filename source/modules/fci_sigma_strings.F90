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

  !> Working-set ceiling for the (t1, t2) pair, in bytes.  The alpha-string
  !> range is blocked to stay under it, so a wide active space costs more
  !> blocks rather than an allocation failure.
  integer(i8), parameter :: work_bytes_cap = 512_i8 * 1024_i8 * 1024_i8

  public :: fci_sigma_strings, rdm12_strings

contains

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

    integer :: na, nb, nkl, nlink_a, nlink_b, nthr
    integer(i8) :: nastr, nbstr, ia, ib, k, base, srow, drow
    integer(i8) :: e, m
    integer :: p, q, r, s, pq, rs, t
    real(dp) :: sgn
    integer(i8), allocatable :: astr(:), bstr(:)
    integer, allocatable :: ntab_a(:), ntab_b(:)
    integer, allocatable :: kl_a(:, :), kl_b(:, :)
    integer(i8), allocatable :: src_a(:, :), src_b(:, :)
    real(dp), allocatable :: sgn_a(:, :), sgn_b(:, :)
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
    ! T is (ndet, nact^2) dense; decline (Python-fallback semantics) rather
    ! than thrash when it would not fit comfortably.
    if (ndet * int(nkl, i8) * 8_i8 > work_bytes_cap * 4_i8) return

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
    nlink_b = max(1, nb * (norb - nb) + nb)
    allocate(ntab_a(nastr), kl_a(nlink_a, nastr), src_a(nlink_a, nastr), &
             sgn_a(nlink_a, nastr))
    allocate(ntab_b(nbstr), kl_b(nlink_b, nbstr), src_b(nlink_b, nbstr), &
             sgn_b(nlink_b, nbstr))
    call build_link_table(norb, nastr, astr, nlink_a, ntab_a, kl_a, src_a, sgn_a)
    call build_link_table(norb, nbstr, bstr, nlink_b, ntab_b, kl_b, src_b, sgn_b)

    allocate(tmat(ndet, nkl), gram(nkl, nkl))
    tmat = 0.0_dp
    nthr = max(1, nthreads)

    ! alpha channel of E_pq: rows of one target alpha string are the
    ! contiguous run (ia-1)*nbstr + 1 .. ia*nbstr
    !$omp parallel do num_threads(nthr) schedule(static) default(shared) &
    !$omp private(ia, e, srow, drow, t, sgn, m)
    do ia = 1_i8, nastr
      drow = (ia - 1_i8) * nbstr
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

    ! beta channel: excitations stay inside one alpha block
    !$omp parallel do num_threads(nthr) schedule(static) default(shared) &
    !$omp private(ia, ib, base, e, srow, drow, t, sgn)
    do ia = 1_i8, nastr
      base = (ia - 1_i8) * nbstr
      do ib = 1_i8, nbstr
        drow = base + ib
        do e = 1_i8, int(ntab_b(ib), i8)
          srow = base + src_b(e, ib)
          t = kl_b(e, ib) + 1
          tmat(drow, t) = tmat(drow, t) + sgn_b(e, ib) * civec(srow)
        end do
      end do
    end do
    !$omp end parallel do

    ! d1 = T^T c
    call dgemv('T', int(ndet), nkl, 1.0_dp, tmat, int(ndet), civec, 1, &
               0.0_dp, d1, 1)
    ! gram = T^T T
    call dgemm('T', 'N', nkl, nkl, int(ndet), 1.0_dp, tmat, int(ndet), &
               tmat, int(ndet), 0.0_dp, gram, nkl)

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
               ntab_b, kl_b, src_b, sgn_b, tmat, gram)
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

    integer :: norb, na, nb, nkl, nlink_a, nlink_b, nthr
    integer(i8) :: nastr, nbstr, ia, ib, k, ns
    integer(i8), allocatable :: astr(:), bstr(:)
    integer, allocatable :: ntab_a(:), ntab_b(:)
    integer, allocatable :: kl_a(:, :), kl_b(:, :)
    integer(i8), allocatable :: src_a(:, :), src_b(:, :)
    real(dp), allocatable :: sgn_a(:, :), sgn_b(:, :)
    real(dp), allocatable :: emat(:, :), t1(:, :), t2(:, :)
    integer(i8) :: nblk, blk_start, blk_len, nrow, inner
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
    nlink_b = max(1, nb * (norb - nb) + nb)

    allocate(ntab_a(nastr), kl_a(nlink_a, nastr), src_a(nlink_a, nastr), &
             sgn_a(nlink_a, nastr))
    allocate(ntab_b(nbstr), kl_b(nlink_b, nbstr), src_b(nlink_b, nbstr), &
             sgn_b(nlink_b, nbstr))
    call build_link_table(norb, nastr, astr, nlink_a, ntab_a, kl_a, src_a, sgn_a)
    call build_link_table(norb, nbstr, bstr, nlink_b, ntab_b, kl_b, src_b, sgn_b)

    allocate(emat(nkl, nkl))
    call build_emat(norb, ns, na + nb, hspin, gspin, nkl, emat)

    nthr = max(1, nthreads)
    inner = nbstr * int(nvec, i8)

    bytes_per_alpha = 2_i8 * inner * int(nkl, i8) * 8_i8
    nblk = max(1_i8, work_bytes_cap / max(bytes_per_alpha, 1_i8))
    nblk = min(nblk, nastr)

    y(1:ndet * int(nvec, i8)) = 0.0_dp

    blk_start = 1_i8
    do while (blk_start <= nastr)
      blk_len = min(nblk, nastr - blk_start + 1_i8)
      nrow = blk_len * inner
      allocate(t1(nrow, nkl), t2(nrow, nkl))

      call gather(norb, nastr, nbstr, nvec, blk_start, blk_len, inner, nkl, &
                  nlink_a, ntab_a, kl_a, src_a, sgn_a, &
                  nlink_b, ntab_b, kl_b, src_b, sgn_b, nthr, x, nrow, t1)

      call dgemm('N', 'N', int(nrow), nkl, nkl, 1.0_dp, t1, int(nrow), &
                 emat, nkl, 0.0_dp, t2, int(nrow))

      call scatter(norb, nastr, nbstr, nvec, blk_start, blk_len, inner, nkl, &
                   nlink_a, ntab_a, kl_a, src_a, sgn_a, &
                   nlink_b, ntab_b, kl_b, src_b, sgn_b, nthr, nrow, t2, y)

      deallocate(t1, t2)
      blk_start = blk_start + blk_len
    end do

    deallocate(emat, astr, bstr)
    deallocate(ntab_a, kl_a, src_a, sgn_a, ntab_b, kl_b, src_b, sgn_b)
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
                    nlink_b, ntab_b, kl_b, src_b, sgn_b, nthr, x, nrow, t1)
    integer, intent(in) :: norb, nvec, nkl, nlink_a, nlink_b, nthr
    integer(i8), intent(in) :: nastr, nbstr, blk_start, blk_len, inner, nrow
    integer, intent(in) :: ntab_a(nastr), ntab_b(nbstr)
    integer, intent(in) :: kl_a(nlink_a, nastr), kl_b(nlink_b, nbstr)
    integer(i8), intent(in) :: src_a(nlink_a, nastr), src_b(nlink_b, nbstr)
    real(dp), intent(in) :: sgn_a(nlink_a, nastr), sgn_b(nlink_b, nbstr)
    real(dp), intent(in) :: x(*)
    real(dp), intent(out) :: t1(nrow, nkl)

    integer(i8) :: ia, ib, e, srow, drow, nv, sa, so, m
    integer :: t
    real(dp) :: s

    nv = int(nvec, i8)
    t1 = 0.0_dp

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
    ! partitions the writes
    !$omp parallel do num_threads(nthr) schedule(static) default(shared) &
    !$omp private(ia, ib, e, sa, so, srow, drow, t, s, m)
    do ia = blk_start, blk_start + blk_len - 1_i8
      sa = (ia - blk_start) * inner
      so = (ia - 1_i8) * inner
      do ib = 1_i8, nbstr
        drow = sa + (ib - 1_i8) * nv
        do e = 1_i8, int(ntab_b(ib), i8)
          srow = so + (src_b(e, ib) - 1_i8) * nv
          t = kl_b(e, ib) + 1
          s = sgn_b(e, ib)
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
                     nlink_b, ntab_b, kl_b, src_b, sgn_b, nthr, nrow, t2, y)
    integer, intent(in) :: norb, nvec, nkl, nlink_a, nlink_b, nthr
    integer(i8), intent(in) :: nastr, nbstr, blk_start, blk_len, inner, nrow
    integer, intent(in) :: ntab_a(nastr), ntab_b(nbstr)
    integer, intent(in) :: kl_a(nlink_a, nastr), kl_b(nlink_b, nbstr)
    integer(i8), intent(in) :: src_a(nlink_a, nastr), src_b(nlink_b, nbstr)
    real(dp), intent(in) :: sgn_a(nlink_a, nastr), sgn_b(nlink_b, nbstr)
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
    !$omp private(ia, ib, e, sa, so, srow, drow, t, s, m)
    do ia = blk_start, blk_end
      sa = (ia - blk_start) * inner
      so = (ia - 1_i8) * inner
      do ib = 1_i8, nbstr
        drow = so + (ib - 1_i8) * nv
        do e = 1_i8, int(ntab_b(ib), i8)
          srow = sa + (src_b(e, ib) - 1_i8) * nv
          t = kl_b(e, ib) + 1
          s = sgn_b(e, ib)
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

end module fci_sigma_strings_mod
