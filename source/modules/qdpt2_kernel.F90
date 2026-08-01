!> @brief Matrix-free QDPT2 streaming kernel (diagonal-Fock zeroth order).
!>
!> Fortran/OpenMP port of the hot loop of pyoqp/oqp/library/qdpt2_direct.py:
!> for every reference determinant D (alpha/beta bit words over the frozen-
!> core-reduced orbital space) generate all EXTERNAL single and double
!> excitations, evaluate the Slater-Condon matrix element <Phi|H|D> with the
!> chemist-notation spatial ERI (pq|rs), and accumulate the per-state
!> first-order couplings V_I(Phi) += <Phi|H|D> * C(D, I) together with the
!> diagonal zeroth-order energy E0(Phi) = sum_occ eps.
!>
!> Accumulation strategy (v2): SORT/REDUCE, not hashing.  The first version
!> used per-thread open-addressing hash tables; measured at 114.5M streamed
!> terms the random-probe tables were DRAM-latency-bound (73 s serial,
!> negative thread scaling) while the NumPy lexsort merge ran in ~6 s.  This
!> version therefore mirrors the cache-friendly design in Fortran:
!>
!>   1. references are partitioned over threads; every reference produces an
!>      identical term count, so each thread writes into an EXACT slice of one
!>      shared term buffer (no locks, no over-allocation);
!>   2. each thread sorts its slice by the 128-bit (alpha,beta) key
!>      (index quicksort + gather) and reduces adjacent duplicates in place;
!>   3. a serial T-way merge of the sorted unique slices writes the final
!>      unique external space, summing V across threads on key ties.
!>
!> Conventions (must match qdpt2_direct.py exactly):
!>  * orbital bits 0..norb-1 per spin word; core = bits 0..ncore-1, then nact
!>    active bits, the rest virtual; "internal" targets (full core, empty
!>    virtual window in both spins) are excluded at generation;
!>  * single-excitation element = mean-field F_eff[a,i] from D's occupation
!>    (k=i self-term cancels exactly between direct and exchange);
!>  * same-spin double (i<j -> a<b): (ai|bj) - (aj|bi), sequential-operator
!>    bit-count phases; alpha-beta double: (ai|bj), within-word phases
!>    (cross-spin crossings cancel pairwise);
!>  * eri is the C-contiguous chemist numpy array eri[p,q,r,s] = (pq|rs);
!>    h1e C-contiguous [p,q]; cvec C-contiguous C[nsup, nstate];
!>    out_v C-contiguous [cap, nstate] (slot-major).
!>
!> Returns the number of unique external determinants compacted into the
!> leading output slots, or -1 when they exceed cap (caller sizes cap at the
!> exact streamed-term upper bound, so this cannot trigger in normal use).
!> Keys use bits 0..norb-1 only (norb <= 63).
module qdpt2_kernel_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

contains

  pure function popcnt_below(w, pos) result(n)
    integer(i8), intent(in) :: w
    integer, intent(in) :: pos
    integer :: n
    n = popcnt(iand(w, shiftl(1_i8, pos) - 1_i8))
  end function popcnt_below

  pure function between_parity(w, p1, p2) result(s)
    integer(i8), intent(in) :: w
    integer, intent(in) :: p1, p2
    real(dp) :: s
    integer :: lo, hi, n
    lo = min(p1, p2); hi = max(p1, p2)
    n = popcnt(iand(w, shiftl(1_i8, hi) - shiftl(1_i8, lo + 1)))
    s = merge(1.0_dp, -1.0_dp, mod(n, 2) == 0)
  end function between_parity

  pure function terms_per_reference(norb, na_occ, nb_occ) result(t)
    !! upper bound on generated terms per reference (internal targets are
    !! skipped at generation but counted here) -- must match the Python
    !! _terms_per_reference formula
    integer, intent(in) :: norb, na_occ, nb_occ
    integer(i8) :: t
    integer(i8) :: nva, nvb
    nva = norb - na_occ; nvb = norb - nb_occ
    t = int(na_occ, i8) * nva + int(nb_occ, i8) * nvb
    t = t + (int(na_occ, i8) * (na_occ - 1) / 2) * (nva * (nva - 1) / 2)
    t = t + (int(nb_occ, i8) * (nb_occ - 1) / 2) * (nvb * (nvb - 1) / 2)
    t = t + int(na_occ, i8) * nva * int(nb_occ, i8) * nvb
  end function terms_per_reference

  !------------------------------------------------------- two-key index sort
  subroutine sort_pairs(n, ka, kb, perm)
    !! iterative quicksort of perm so that (ka,kb)(perm) ascends lexicographically
    integer(i8), intent(in) :: n
    integer(i8), intent(in) :: ka(n), kb(n)
    integer(i8), intent(inout) :: perm(n)
    integer(i8) :: stack_lo(64), stack_hi(64)
    integer :: sp
    integer(i8) :: lo, hi, i, j, tmp, pa, pb

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
              if (.not. key_gt(perm(j), tmp)) exit
              perm(j + 1) = perm(j)
              j = j - 1
            end do
            perm(j + 1) = tmp
          end do
          exit
        end if
        tmp = perm(ishft(lo + hi, -1))
        pa = ka(tmp); pb = kb(tmp)
        i = lo; j = hi
        do
          do while (ka(perm(i)) < pa .or. (ka(perm(i)) == pa .and. kb(perm(i)) < pb))
            i = i + 1
          end do
          do while (ka(perm(j)) > pa .or. (ka(perm(j)) == pa .and. kb(perm(j)) > pb))
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

  contains

    pure function key_gt(ia, ib) result(g)
      integer(i8), intent(in) :: ia, ib
      logical :: g
      g = ka(ia) > ka(ib) .or. (ka(ia) == ka(ib) .and. kb(ia) > kb(ib))
    end function key_gt

  end subroutine sort_pairs

  !------------------------------------------------------------ main kernel
  function qdpt2_stream_kernel(norb, ncore, nact, nsup, nstate, nthreads, cap, &
                               sup_a, sup_b, cvec, h1e, eri, eps, &
                               out_ka, out_kb, out_e0, out_v) &
      result(n_uniq) bind(C, name="qdpt2_stream_kernel")
    integer(c_int32_t), value :: norb, ncore, nact, nstate, nthreads
    integer(i8), value :: nsup, cap
    integer(i8), intent(in) :: sup_a(*), sup_b(*)
    real(dp), intent(in) :: cvec(*), h1e(*), eri(*), eps(*)
    integer(i8), intent(inout) :: out_ka(*), out_kb(*)
    real(dp), intent(inout) :: out_e0(*), out_v(*)
    integer(i8) :: n_uniq

    integer(i8) :: cmask, vmask, tpr, total
    integer(i8), allocatable :: wka(:), wkb(:), perm(:)
    real(dp), allocatable :: we0(:), wv(:,:)
    integer(i8), allocatable :: t_off(:), t_cnt(:), t_uniq(:), head(:)
    integer :: nthr, t, tbest, norb_i, nstate_i, st
    integer(i8) :: d, lo_ref, hi_ref, pos, i, run_a, run_b

    norb_i = int(norb)
    nstate_i = int(nstate)
    cmask = shiftl(1_i8, int(ncore)) - 1_i8
    vmask = iand(not(shiftl(1_i8, int(ncore) + int(nact)) - 1_i8), &
                 shiftl(1_i8, norb_i) - 1_i8)

    tpr = terms_per_reference(norb_i, popcnt(sup_a(1)), popcnt(sup_b(1)))
    total = nsup * tpr

    ! honor the caller's explicit thread request: the ambient OMP default is
    ! clamped to 1 by the core's BLAS-thread control, and num_threads() below
    ! overrides it regardless
    nthr = max(1, int(nthreads))
    nthr = int(min(int(nthr, i8), max(1_i8, nsup / 4_i8)))

    allocate(wka(total), wkb(total), we0(total), wv(nstate_i, total))
    allocate(t_off(nthr + 1), t_cnt(nthr), t_uniq(nthr), head(nthr))
    do t = 1, nthr + 1
      t_off(t) = ((t - 1) * nsup / nthr) * tpr + 1
    end do

    !$omp parallel num_threads(nthr) default(shared) private(t, d, lo_ref, hi_ref, pos, perm, i, run_a, run_b)
    !$omp do schedule(static, 1)
    do t = 1, nthr
      lo_ref = (int(t, i8) - 1) * nsup / nthr + 1
      hi_ref = int(t, i8) * nsup / nthr
      pos = t_off(t) - 1
      do d = lo_ref, hi_ref
        call stream_one(int(d), norb_i, nstate_i, cmask, vmask, &
                        sup_a(d), sup_b(d), cvec, h1e, eri, eps, &
                        wka, wkb, we0, wv, pos)
      end do
      t_cnt(t) = pos - (t_off(t) - 1)

      ! sort this thread's slice by (ka, kb) and gather the payloads
      allocate(perm(t_cnt(t)))
      do i = 1, t_cnt(t)
        perm(i) = t_off(t) - 1 + i
      end do
      call sort_pairs(t_cnt(t), wka, wkb, perm)
      call gather_slice(t_off(t), t_cnt(t), perm, nstate_i, wka, wkb, we0, wv)
      deallocate(perm)

      ! reduce adjacent duplicates in place
      t_uniq(t) = 0
      i = t_off(t)
      do while (i < t_off(t) + t_cnt(t))
        run_a = wka(i); run_b = wkb(i)
        pos = t_off(t) + t_uniq(t)
        if (pos /= i) then
          wka(pos) = run_a; wkb(pos) = run_b; we0(pos) = we0(i)
          wv(:, pos) = wv(:, i)
        end if
        i = i + 1
        do while (i < t_off(t) + t_cnt(t))
          if (wka(i) /= run_a .or. wkb(i) /= run_b) exit
          wv(:, pos) = wv(:, pos) + wv(:, i)
          i = i + 1
        end do
        t_uniq(t) = t_uniq(t) + 1
      end do
    end do
    !$omp end do
    !$omp end parallel

    ! serial T-way merge of the sorted unique slices
    n_uniq = 0
    do t = 1, nthr
      head(t) = 0
    end do
    do
      tbest = 0
      do t = 1, nthr
        if (head(t) >= t_uniq(t)) cycle
        if (tbest == 0) then
          tbest = t
        else
          i = t_off(t) + head(t)
          d = t_off(tbest) + head(tbest)
          if (wka(i) < wka(d) .or. (wka(i) == wka(d) .and. wkb(i) < wkb(d))) tbest = t
        end if
      end do
      if (tbest == 0) exit
      if (n_uniq >= cap) then
        n_uniq = -1_i8
        deallocate(wka, wkb, we0, wv, t_off, t_cnt, t_uniq, head)
        return
      end if
      d = t_off(tbest) + head(tbest)
      n_uniq = n_uniq + 1
      out_ka(n_uniq) = wka(d); out_kb(n_uniq) = wkb(d)
      out_e0(n_uniq) = we0(d)
      do st = 1, nstate_i
        out_v((n_uniq - 1) * nstate_i + st) = wv(st, d)
      end do
      head(tbest) = head(tbest) + 1
      ! absorb equal keys from every other slice
      do t = 1, nthr
        do while (head(t) < t_uniq(t))
          i = t_off(t) + head(t)
          if (wka(i) /= out_ka(n_uniq) .or. wkb(i) /= out_kb(n_uniq)) exit
          do st = 1, nstate_i
            out_v((n_uniq - 1) * nstate_i + st) = &
                out_v((n_uniq - 1) * nstate_i + st) + wv(st, i)
          end do
          head(t) = head(t) + 1
        end do
      end do
    end do

    deallocate(wka, wkb, we0, wv, t_off, t_cnt, t_uniq, head)
  end function qdpt2_stream_kernel

  subroutine gather_slice(off, cnt, perm, nstate, wka, wkb, we0, wv)
    !! apply the sort permutation to the slice payloads (out-of-place gather)
    integer(i8), intent(in) :: off, cnt
    integer(i8), intent(in) :: perm(cnt)
    integer, intent(in) :: nstate
    integer(i8), intent(inout) :: wka(*), wkb(*)
    real(dp), intent(inout) :: we0(*)
    real(dp), intent(inout) :: wv(nstate, *)
    integer(i8), allocatable :: tk(:)
    real(dp), allocatable :: td(:), tv(:,:)
    integer(i8) :: i

    allocate(tk(cnt), td(cnt), tv(nstate, cnt))
    do i = 1, cnt
      tk(i) = wka(perm(i))
    end do
    wka(off:off + cnt - 1) = tk
    do i = 1, cnt
      tk(i) = wkb(perm(i))
    end do
    wkb(off:off + cnt - 1) = tk
    do i = 1, cnt
      td(i) = we0(perm(i))
      tv(:, i) = wv(:, perm(i))
    end do
    we0(off:off + cnt - 1) = td
    wv(:, off:off + cnt - 1) = tv
    deallocate(tk, td, tv)
  end subroutine gather_slice

  !----------------------------------------------------- per-reference stream
  subroutine stream_one(d, norb, nstate, cmask, vmask, wa, wb, cvec, h1e, eri, &
                        eps, wka, wkb, we0, wv, pos)
    integer, intent(in) :: d, norb, nstate
    integer(i8), intent(in) :: cmask, vmask, wa, wb
    real(dp), intent(in) :: cvec(*), h1e(*), eri(*), eps(*)
    integer(i8), intent(inout) :: wka(*), wkb(*)
    real(dp), intent(inout) :: we0(*)
    real(dp), intent(inout) :: wv(nstate, *)
    integer(i8), intent(inout) :: pos

    integer :: occ_a(norb), occ_b(norb), virt_a(norb), virt_b(norb)
    integer :: na_o, nb_o, na_v, nb_v
    real(dp) :: fa(norb, norb), fb(norb, norb), wgt(nstate)
    real(dp) :: e0_ref, elem, ph
    integer(i8) :: ka, kb, w1, w2, w3, w4
    integer :: io, jo, av, bv, i, j, a, b, st

    do st = 1, nstate
      wgt(st) = cvec((d - 1) * nstate + st)
    end do

    call decode(wa, norb, occ_a, na_o, virt_a, na_v)
    call decode(wb, norb, occ_b, nb_o, virt_b, nb_v)

    e0_ref = 0.0_dp
    do io = 1, na_o
      e0_ref = e0_ref + eps(occ_a(io) + 1)
    end do
    do io = 1, nb_o
      e0_ref = e0_ref + eps(occ_b(io) + 1)
    end do

    call mean_field(norb, occ_a, na_o, occ_b, nb_o, h1e, eri, fa)
    call mean_field(norb, occ_b, nb_o, occ_a, na_o, h1e, eri, fb)

    ! ---- alpha singles
    do io = 1, na_o
      i = occ_a(io)
      do av = 1, na_v
        a = virt_a(av)
        ka = ior(ieor(wa, shiftl(1_i8, i)), shiftl(1_i8, a))
        if (is_internal(ka, wb, cmask, vmask)) cycle
        elem = fa(a + 1, i + 1) * between_parity(wa, i, a)
        call push(ka, wb, e0_ref - eps(i + 1) + eps(a + 1), elem)
      end do
    end do

    ! ---- beta singles
    do io = 1, nb_o
      i = occ_b(io)
      do av = 1, nb_v
        a = virt_b(av)
        kb = ior(ieor(wb, shiftl(1_i8, i)), shiftl(1_i8, a))
        if (is_internal(wa, kb, cmask, vmask)) cycle
        elem = fb(a + 1, i + 1) * between_parity(wb, i, a)
        call push(wa, kb, e0_ref - eps(i + 1) + eps(a + 1), elem)
      end do
    end do

    ! ---- same-spin doubles (alpha then beta)
    call same_spin_doubles(.true.)
    call same_spin_doubles(.false.)

    ! ---- alpha-beta doubles
    do io = 1, na_o
      i = occ_a(io)
      do av = 1, na_v
        a = virt_a(av)
        ka = ior(ieor(wa, shiftl(1_i8, i)), shiftl(1_i8, a))
        do jo = 1, nb_o
          j = occ_b(jo)
          do bv = 1, nb_v
            b = virt_b(bv)
            kb = ior(ieor(wb, shiftl(1_i8, j)), shiftl(1_i8, b))
            if (is_internal(ka, kb, cmask, vmask)) cycle
            elem = eri(idx4(a, i, b, j, norb)) &
                 * between_parity(wa, i, a) * between_parity(wb, j, b)
            call push(ka, kb, &
                      e0_ref - eps(i + 1) + eps(a + 1) - eps(j + 1) + eps(b + 1), &
                      elem)
          end do
        end do
      end do
    end do

  contains

    subroutine push(pka, pkb, pe0, pelem)
      integer(i8), intent(in) :: pka, pkb
      real(dp), intent(in) :: pe0, pelem
      pos = pos + 1
      wka(pos) = pka; wkb(pos) = pkb; we0(pos) = pe0
      wv(:, pos) = pelem * wgt
    end subroutine push

    subroutine same_spin_doubles(alpha)
      logical, intent(in) :: alpha
      integer :: n_o, n_v
      integer :: occ(norb), virt(norb)
      integer(i8) :: w
      real(dp) :: s1, s2, s3, s4
      if (alpha) then
        n_o = na_o; n_v = na_v; occ = occ_a; virt = virt_a; w = wa
      else
        n_o = nb_o; n_v = nb_v; occ = occ_b; virt = virt_b; w = wb
      end if
      do io = 1, n_o - 1
        i = occ(io)
        do jo = io + 1, n_o
          j = occ(jo)
          w1 = ieor(w, shiftl(1_i8, i))
          s1 = merge(1.0_dp, -1.0_dp, mod(popcnt_below(w, i), 2) == 0)
          w2 = ieor(w1, shiftl(1_i8, j))
          s2 = merge(1.0_dp, -1.0_dp, mod(popcnt_below(w1, j), 2) == 0)
          do av = 1, n_v - 1
            a = virt(av)
            do bv = av + 1, n_v
              b = virt(bv)
              w3 = ior(w2, shiftl(1_i8, b))
              s3 = merge(1.0_dp, -1.0_dp, mod(popcnt_below(w2, b), 2) == 0)
              w4 = ior(w3, shiftl(1_i8, a))
              s4 = merge(1.0_dp, -1.0_dp, mod(popcnt_below(w3, a), 2) == 0)
              ph = s1 * s2 * s3 * s4
              elem = (eri(idx4(a, i, b, j, norb)) - eri(idx4(a, j, b, i, norb))) * ph
              if (alpha) then
                if (is_internal(w4, wb, cmask, vmask)) cycle
                call push(w4, wb, &
                          e0_ref - eps(i + 1) - eps(j + 1) + eps(a + 1) + eps(b + 1), elem)
              else
                if (is_internal(wa, w4, cmask, vmask)) cycle
                call push(wa, w4, &
                          e0_ref - eps(i + 1) - eps(j + 1) + eps(a + 1) + eps(b + 1), elem)
              end if
            end do
          end do
        end do
      end do
    end subroutine same_spin_doubles

  end subroutine stream_one

  pure function is_internal(ka, kb, cmask, vmask) result(res)
    integer(i8), intent(in) :: ka, kb, cmask, vmask
    logical :: res
    res = (iand(ka, cmask) == cmask) .and. (iand(kb, cmask) == cmask) &
        .and. (iand(ka, vmask) == 0_i8) .and. (iand(kb, vmask) == 0_i8)
  end function is_internal

  pure function idx4(p, q, r, s, norb) result(idx)
    integer, intent(in) :: p, q, r, s, norb
    integer(i8) :: idx
    idx = ((int(p, i8) * norb + q) * norb + r) * norb + s + 1_i8
  end function idx4

  subroutine decode(w, norb, occ, n_occ, virt, n_virt)
    integer(i8), intent(in) :: w
    integer, intent(in) :: norb
    integer, intent(out) :: occ(norb), virt(norb)
    integer, intent(out) :: n_occ, n_virt
    integer :: p
    n_occ = 0; n_virt = 0
    do p = 0, norb - 1
      if (btest(w, p)) then
        n_occ = n_occ + 1
        occ(n_occ) = p
      else
        n_virt = n_virt + 1
        virt(n_virt) = p
      end if
    end do
  end subroutine decode

  subroutine mean_field(norb, occ_same, n_same, occ_other, n_other, h1e, eri, f)
    integer, intent(in) :: norb, n_same, n_other
    integer, intent(in) :: occ_same(norb), occ_other(norb)
    real(dp), intent(in) :: h1e(*), eri(*)
    real(dp), intent(out) :: f(norb, norb)
    integer :: p, q, kk, k
    do q = 0, norb - 1
      do p = 0, norb - 1
        f(p + 1, q + 1) = h1e(int(p, i8) * norb + q + 1)
      end do
    end do
    do kk = 1, n_same
      k = occ_same(kk)
      do q = 0, norb - 1
        do p = 0, norb - 1
          f(p + 1, q + 1) = f(p + 1, q + 1) + eri(idx4(p, q, k, k, norb)) &
                                            - eri(idx4(p, k, k, q, norb))
        end do
      end do
    end do
    do kk = 1, n_other
      k = occ_other(kk)
      do q = 0, norb - 1
        do p = 0, norb - 1
          f(p + 1, q + 1) = f(p + 1, q + 1) + eri(idx4(p, q, k, k, norb))
        end do
      end do
    end do
  end subroutine mean_field

end module qdpt2_kernel_mod
