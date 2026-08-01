!> @brief Fortran engine for the fixed-CI half of the analytic CASSCF orbital
!> Hessian of pyoqp casscf_hessian.py.
!>
!> The Python assembly walks the non-redundant rotation pairs and, for each,
!> materialises the one-index derivative integrals as a full nbf^4 tensor and
!> contracts them against the 2-RDM:
!>
!>     T^(l) = one-index derivative of g along K^(l)          O(nbf^4) storage
!>     Z^(l) = (D t - t D)^T + 2 sum_bcd G[n,b,c,d] T^(l)[m,b,c,d]   O(nbf^5)
!>     B[:,l] = Z^(l)[p_k,q_k] - Z^(l)[q_k,p_k]
!>
!> so the loop costs npar * O(nbf^5).  With a large basis and a small active
!> space -- the regime the dense-spectrum analytic Hessian is actually good
!> for -- that dominates everything else in the module.
!>
!> Neither the nbf^5 nor the nbf^4 is necessary, because T^(l) is sparse: the
!> generator K^(l) = e_PQ - e_QP applied to each of the four slots leaves only
!> eight slabs of g,
!>
!>     T[a,b,c,d] = d(a,Q)g[P,b,c,d] - d(a,P)g[Q,b,c,d]
!>                + d(b,Q)g[a,P,c,d] - d(b,P)g[a,Q,c,d]
!>                + d(c,Q)g[a,b,P,d] - d(c,P)g[a,b,Q,d]
!>                + d(d,Q)g[a,b,c,P] - d(d,P)g[a,b,c,Q].
!>
!> Contracting slab by slab, the slot-1 term collapses onto a single
!> pair-independent intermediate and the other three become nbf^4 GEMMs:
!>
!>     Y[n,a]      = sum_bcd G[n,b,c,d] g[a,b,c,d]      ONCE, not per pair
!>     A[X,Y][m,n] = sum_cd  g[m,X,c,d] G[n,Y,c,d]
!>     B[X,Y][m,n] = sum_bd  g[m,b,X,d] G[n,b,Y,d]
!>     C[X,Y][m,n] = sum_bc  g[m,b,c,X] G[n,b,c,Y]
!>
!>     C1[m,n] = A[P,Q]-A[Q,P] + B[P,Q]-B[Q,P] + C[P,Q]-C[Q,P]
!>     C1[Q,n] += Y[n,P],   C1[P,n] -= Y[n,Q]
!>
!> so the per-pair cost drops from O(nbf^5) to O(nbf^4) and no nbf^4 temporary
!> is ever formed.  The A contraction needs no gather at all: g[:,X,:,:] is
!> already a valid column-major (nbf^2 x nbf) matrix with leading dimension
!> nbf^3.  B and C fix an interior index, so those two need their operands
!> gathered into contiguous buffers.
!>
!> The gathers, not the GEMMs, are what this loop costs.  Timing the two
!> separately at nbf=48 put 80% of the routine in its nbf^3 gathers and only
!> 20% in the six nbf^4 GEMMs -- the GEMMs already ran at ~110 GFlop/s,
!> while a gather that fixes an interior index of a 42 MB array reads one
!> nbf-long run per 18 kB stride and sustains ~5 GB/s.  Flops were never the
!> problem.
!>
!> So no gather is done for one index at a time.  Both indices are gathered in
!> RANGES, which is what the pair list already comes in:
!>
!>   * q.  _nonredundant_pairs emits p outermost and sweeps a whole space in
!>     q -- (active, inactive), (virtual, inactive), (virtual, active) -- so
!>     the pairs sharing a p are consecutive AND their q is a contiguous
!>     ascending range.  The six q-side operands are gathered for that whole
!>     range at once and keyed on it, so consecutive runs over the same space
!>     reuse them: two builds per call instead of four gathers per pair.
!>   * p.  p is likewise consecutive, so its four operands are gathered a
!>     chunk of p at a time.  Column kp of such a chunk is the (nbf^2 x nbf)
!>     matrix for p = pc0 + kp with leading dimension npc*nbf^2, so batching
!>     the gather costs the GEMM nothing but a leading dimension.
!>
!> Gathering a range is not merely amortised, it is cheaper per element: the
!> B operand becomes one nq*nbf-long run per stride instead of nq nbf-long
!> ones, and the C operand becomes nq-long runs instead of single elements.
!> Measured at nbf=48 the batched gathers move a double 22x faster than the
!> per-pair ones did.
!>
!> Batching q also fixes the GEMM shape.  Per pair the contraction is
!> (nbf x nbf^2).(nbf^2 x nbf), an inner-product shape whose k dominates both
!> free dimensions; over a q-range it becomes (nbf x nbf^2).(nbf^2 x nq*nbf)
!> for the same flops.  The A term, which needed no gather per pair, does need
!> one to be batched (its q operand is a strided set of columns of G), and that
!> gather is a plain nbf^2 block copy.
!>
!> The two families differ in which operand carries q, so they accumulate into
!> two buffers with the pair index on opposite sides:
!>
!>     resp[m, (q,n)] = A[P,Q]+B[P,Q]+C[P,Q]      q batched on the right
!>     resm[(q,m), n] = A[Q,P]+B[Q,P]+C[Q,P]      q batched on the left
!>     C1[m,n]        = resp[m,(q,n)] - resm[(q,m),n]
!>
!> Nothing here ASSUMES that ordering, it only exploits it: the runs are
!> scanned off the pair list, and a caller whose list is ordered differently
!> gets shorter runs, down to one pair per batch.
!>
!> The single contraction above (rather than the four the Python originally
!> evaluated) is justified separately: T inherits the full eight-fold symmetry
!> of the real MO ERIs, so all four coincide.  See casscf_hessian._z_matrix and
!> tests/test_casscf_zmatrix.py.
!>
!> The folded active derivative integrals f^(l), g^(l) that the CI-relaxation
!> half consumes are produced here too, evaluated directly from the eight slabs
!> rather than by slicing a materialised T.
!>
!> All arrays cross the boundary in the caller's C order (last index fastest).
module casscf_hess_bmat_mod
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_double
  implicit none
  private

  integer, parameter :: i8 = c_int64_t
  integer, parameter :: dp = c_double

  !> Budget for the six q-batched operand buffers together, in doubles
  !> (256 MiB).  They are O(nq * nbf^3) against inputs that are already
  !> 2 * nbf^4, so this only ever binds for a large inactive or active space.
  integer(i8), parameter :: q_budget = 33554432_i8
  !> Budget for the four p-batched operand buffers together, in doubles
  !> (128 MiB).  Measured: 128 MiB is the optimum at nbf=36 and nbf=48, and
  !> both 32 MiB and 512 MiB are slower -- too small stops amortising the
  !> strided reads, too large pushes the buffer out of cache for the GEMM
  !> that consumes it.
  integer(i8), parameter :: p_budget = 16777216_i8

  !> Threading threshold for a gather, in doubles moved.  The gathers are the
  !> only part of this routine BLAS is not already threading, and they are
  !> latency-bound, so spreading them wins ~1.2x at nbf=48.  Below this the
  !> parallel region costs more than the gather does: at nbf=24 the whole
  !> routine is 3.5 ms and threading every gather made it 0.85x.
  integer(i8), parameter :: omp_min_words = 262144_i8

  public :: casscf_hess_bmat

contains

  !> Fixed-CI Hessian columns B and the folded active derivative integrals.
  !>
  !> @param[in]  nbf    orbitals
  !> @param[in]  ncore  inactive orbitals
  !> @param[in]  nact   active orbitals
  !> @param[in]  npar   non-redundant rotation pairs
  !> @param[in]  pairs  pair list, C-order [npar,2] as (p,q)
  !> @param[in]  dmat   full spin-summed 1-RDM, C-order [nbf,nbf]
  !> @param[in]  rdm2   full spin-summed chemist 2-RDM, C-order [nbf,nbf,nbf,nbf]
  !> @param[in]  h1e    MO core Hamiltonian, C-order [nbf,nbf]
  !> @param[in]  eri    MO ERIs (chemist), C-order [nbf,nbf,nbf,nbf]
  !> @param[out] bmat   fixed-CI Hessian columns, C-order [npar,npar]
  !> @param[out] fder   folded 1e derivative integrals, C-order [npar,nact,nact]
  !> @param[out] gder   folded 2e derivative integrals, C-order
  !>                    [npar,nact,nact,nact,nact]
  subroutine casscf_hess_bmat(nbf, ncore, nact, npar, pairs, dmat, rdm2, &
                              h1e, eri, bmat, fder, gder) &
      bind(C, name="casscf_hess_bmat")
    integer(c_int32_t), value :: nbf, ncore, nact, npar
    ! bind(C) hands over bare pointers: array dummies must be assumed-SIZE.
    integer(c_int32_t), intent(in) :: pairs(0:*)
    real(dp), intent(in) :: dmat(0:*), rdm2(0:*), h1e(0:*), eri(0:*)
    real(dp), intent(inout) :: bmat(0:*), fder(0:*), gder(0:*)

    ! Default (ILP64) integer kind, matching the linked BLAS.
    integer :: n, nc, na, np, n2, l, p, q, i, j, m
    integer(i8) :: n3
    real(dp), allocatable :: ymat(:,:), cmat(:,:), tmat(:,:), w1(:,:), w2(:,:)
    real(dp), allocatable :: dm(:,:), zmat(:,:)
    ! p-batched operands: a chunk of consecutive p, so the strided reads are
    ! amortised the same way the q side amortises them.
    real(dp), allocatable :: gbp(:,:), rbp(:,:), gcp(:,:), rcp(:,:)
    integer :: pc0, npc, npmax, kp, ldp
    ! q-batched operands and the two accumulators they feed.
    real(dp), allocatable :: qar(:,:), qag(:,:), qbr(:,:), qbg(:,:)
    real(dp), allocatable :: qcr(:,:), qcg(:,:), resp(:,:), resm(:)
    integer :: q0, nrun, nqc, nqmax, maxrun, qb, qs, kk, ldm, last_qs, last_nq

    n = int(nbf)
    nc = int(ncore)
    na = int(nact)
    np = int(npar)
    n2 = n * n
    n3 = int(n, i8) * int(n, i8) * int(n, i8)
    if (n <= 0 .or. np <= 0) return

    ! Longest run of pairs that share a p and sweep a contiguous ascending q,
    ! which is what the q-batched buffers have to hold.  Scanned rather than
    ! assumed: a caller whose pair list is ordered differently simply gets
    ! shorter runs, down to the per-pair form at length one.
    maxrun = 1
    l = 0
    do while (l < np)
      p = int(pairs(2*l))
      q0 = int(pairs(2*l + 1))
      nrun = 1
      do while (l + nrun < np)
        if (int(pairs(2*(l + nrun))) /= p) exit
        if (int(pairs(2*(l + nrun) + 1)) /= q0 + nrun) exit
        nrun = nrun + 1
      end do
      maxrun = max(maxrun, nrun)
      l = l + nrun
    end do
    nqmax = int(max(1_i8, min(int(maxrun, i8), q_budget / (6_i8 * n3))))
    npmax = int(max(1_i8, min(int(n, i8), p_budget / (4_i8 * n3))))

    allocate(ymat(0:n-1, 0:n-1), cmat(0:n-1, 0:n-1), tmat(0:n-1, 0:n-1))
    allocate(w1(0:n-1, 0:n-1), w2(0:n-1, 0:n-1), zmat(0:n-1, 0:n-1))
    allocate(dm(0:n-1, 0:n-1))
    allocate(gbp(0:n2-1, 0:npmax*n-1), rbp(0:n2-1, 0:npmax*n-1))
    allocate(gcp(0:n2-1, 0:npmax*n-1), rcp(0:n2-1, 0:npmax*n-1))
    allocate(qar(0:n2-1, 0:nqmax*n-1), qag(0:n2-1, 0:nqmax*n-1))
    allocate(qbr(0:n2-1, 0:nqmax*n-1), qbg(0:n2-1, 0:nqmax*n-1))
    allocate(qcr(0:n2-1, 0:nqmax*n-1), qcg(0:n2-1, 0:nqmax*n-1))
    allocate(resp(0:n-1, 0:nqmax*n-1))
    allocate(resm(0:int(nqmax, i8)*int(n2, i8) - 1_i8))
    pc0 = -1
    npc = 0
    last_qs = -1
    last_nq = -1

    do i = 0, n - 1
      do j = 0, n - 1
        dm(i, j) = dmat(int(i, i8)*int(n, i8) + int(j, i8))
      end do
    end do

    ! ---- Y[n,a] = sum_bcd G[n,b,c,d] g[a,b,c,d], once for all pairs.
    ! Both buffers are already column-major ((b,c,d), leading index) matrices
    ! with leading dimension nbf^3, so this is one GEMM with no repacking.
    ! ymat(a,n) holds Y[n,a].
    call dgemm('T', 'N', n, n, int(n3), 1.0_dp, eri, int(n3), &
               rdm2, int(n3), 0.0_dp, ymat, n)

    l = 0
    do while (l < np)
      p = int(pairs(2*l))
      q0 = int(pairs(2*l + 1))
      nrun = 1
      do while (l + nrun < np)
        if (int(pairs(2*(l + nrun))) /= p) exit
        if (int(pairs(2*(l + nrun) + 1)) /= q0 + nrun) exit
        nrun = nrun + 1
      end do

      ! The four p-indexed operands are invariant across the run, and the same
      ! batched gathers the q side uses serve them: a chunk of consecutive p is
      ! gathered at once, then column kp of the chunk is the (nbf^2 x nbf)
      ! matrix for p = pc0 + kp with leading dimension npc*nbf^2.
      if (p < pc0 .or. p >= pc0 + npc) then
        pc0 = p
        npc = min(npmax, n - p)
        call gather_q_third(eri, pc0, npc, n, gbp)
        call gather_q_third(rdm2, pc0, npc, n, rbp)
        call gather_q_fourth(eri, pc0, npc, n, gcp)
        call gather_q_fourth(rdm2, pc0, npc, n, rcp)
      end if
      kp = p - pc0
      ldp = npc * n2

      qb = 0
      do while (qb < nrun)
        nqc = min(nqmax, nrun - qb)
        qs = q0 + qb
        ldm = nqc * n

        ! The six q-indexed operands, keyed on the range: consecutive runs
        ! sweep the same space, so this fires twice for a whole call.
        if (qs /= last_qs .or. nqc /= last_nq) then
          call gather_q_second(rdm2, qs, nqc, n, qar)
          call gather_q_second(eri, qs, nqc, n, qag)
          call gather_q_third(rdm2, qs, nqc, n, qbr)
          call gather_q_third(eri, qs, nqc, n, qbg)
          call gather_q_fourth(rdm2, qs, nqc, n, qcr)
          call gather_q_fourth(eri, qs, nqc, n, qcg)
          last_qs = qs
          last_nq = nqc
        end if

        ! resp(m, kk + nqc*n) = (A+B+C)[P,Q][m,n]: q rides on the right
        ! operand, so it lands in the GEMM's free column index.
        ! The A term needs no p gather: g[:,P,:,:] is already a column-major
        ! (nbf^2 x nbf) matrix, base P*nbf^2 and leading dimension nbf^3.
        call dgemm('T', 'N', n, ldm, n2, 1.0_dp, &
                   eri(int(p, i8)*int(n2, i8)), int(n3), qar, n2, 0.0_dp, resp, n)
        call dgemm('T', 'N', n, ldm, n2, 1.0_dp, &
                   gbp(0, kp), ldp, qbr, n2, 1.0_dp, resp, n)
        call dgemm('T', 'N', n, ldm, n2, 1.0_dp, &
                   gcp(0, kp), ldp, qcr, n2, 1.0_dp, resp, n)

        ! resm(kk + nqc*m, n) = (A+B+C)[Q,P][m,n]: here q is on the ERI side,
        ! which is the GEMM's left operand, so it lands in the row index and
        ! the leading dimension is nqc*nbf.
        call dgemm('T', 'N', ldm, n, n2, 1.0_dp, qag, n2, &
                   rdm2(int(p, i8)*int(n2, i8)), int(n3), 0.0_dp, resm, ldm)
        call dgemm('T', 'N', ldm, n, n2, 1.0_dp, qbg, n2, &
                   rbp(0, kp), ldp, 1.0_dp, resm, ldm)
        call dgemm('T', 'N', ldm, n, n2, 1.0_dp, qcg, n2, &
                   rcp(0, kp), ldp, 1.0_dp, resm, ldm)

        do kk = 0, nqc - 1
          do j = 0, n - 1
            do m = 0, n - 1
              cmat(m, j) = resp(m, kk + nqc*j) &
                         - resm(int(kk + nqc*m, i8) + int(ldm, i8)*int(j, i8))
            end do
          end do
          call finish_pair(l + qb + kk, p, qs + kk)
        end do

        qb = qb + nqc
      end do

      l = l + nrun
    end do

    deallocate(ymat, cmat, tmat, w1, w2, zmat, dm)
    deallocate(gbp, rbp, gcp, rcp)
    deallocate(qar, qag, qbr, qbg, qcr, qcg, resp, resm)

  contains

    !> Everything downstream of C1 for one pair: the slot-1 slabs, the one-index
    !> derivative of h, Z, the Hessian column and the folded active integrals.
    !> `cmat` holds C1 on entry, however it was formed.
    subroutine finish_pair(lp, pp, qq)
      integer, intent(in) :: lp, pp, qq
      integer :: ii, jj, kk2, mm, tt, uu, vv, ww
      real(dp) :: a

      ! ---- slot-1 slabs: the only pair dependence is which row they land on
      do jj = 0, n - 1
        cmat(qq, jj) = cmat(qq, jj) + ymat(pp, jj)
        cmat(pp, jj) = cmat(pp, jj) - ymat(qq, jj)
      end do

      ! ---- t = [h, K] for this pair (casscf_hessian._one_index_derivative_h)
      tmat = 0.0_dp
      do ii = 0, n - 1
        tmat(ii, qq) = tmat(ii, qq) + h1e(int(ii, i8)*int(n, i8) + int(pp, i8))
        tmat(ii, pp) = tmat(ii, pp) - h1e(int(ii, i8)*int(n, i8) + int(qq, i8))
      end do
      do jj = 0, n - 1
        tmat(pp, jj) = tmat(pp, jj) - h1e(int(qq, i8)*int(n, i8) + int(jj, i8))
        tmat(qq, jj) = tmat(qq, jj) + h1e(int(pp, i8)*int(n, i8) + int(jj, i8))
      end do

      ! ---- Z[m,n] = (D t - t D)[n,m] + 2 C1[m,n]
      call dgemm('N', 'N', n, n, n, 1.0_dp, dm, n, tmat, n, 0.0_dp, w1, n)
      call dgemm('N', 'N', n, n, n, 1.0_dp, tmat, n, dm, n, 0.0_dp, w2, n)
      do mm = 0, n - 1
        do jj = 0, n - 1
          zmat(mm, jj) = (w1(jj, mm) - w2(jj, mm)) + 2.0_dp * cmat(mm, jj)
        end do
      end do

      ! ---- B[k,l] = Z[p_k,q_k] - Z[q_k,p_k]
      do kk2 = 0, np - 1
        ii = int(pairs(2*kk2))
        jj = int(pairs(2*kk2 + 1))
        bmat(int(kk2, i8)*int(np, i8) + int(lp, i8)) = zmat(ii, jj) - zmat(jj, ii)
      end do

      ! ---- folded active derivative integrals, straight from the slabs
      do tt = 0, na - 1
        do uu = 0, na - 1
          do vv = 0, na - 1
            do ww = 0, na - 1
              gder((((int(lp, i8)*int(na, i8) + int(tt, i8))*int(na, i8) &
                      + int(uu, i8))*int(na, i8) + int(vv, i8))*int(na, i8) &
                      + int(ww, i8)) = &
                  tval(eri, n, pp, qq, nc + tt, nc + uu, nc + vv, nc + ww)
            end do
          end do
        end do
      end do
      do tt = 0, na - 1
        do uu = 0, na - 1
          a = tmat(nc + tt, nc + uu)
          do ii = 0, nc - 1
            a = a + 2.0_dp * tval(eri, n, pp, qq, nc + tt, nc + uu, ii, ii) &
                  - tval(eri, n, pp, qq, nc + tt, ii, ii, nc + uu)
          end do
          fder((int(lp, i8)*int(na, i8) + int(tt, i8))*int(na, i8) + int(uu, i8)) = a
        end do
      end do
    end subroutine finish_pair

    !> out((c,d), kk + nq*x) = src[x, qs+kk, c, d] for kk = 0..nq-1: the
    !> A-term operand over a contiguous q range.  Each (x, q) slab is nbf^2
    !> contiguous doubles, so this is a block copy.
    subroutine gather_q_second(src, qstart, nq, nn, out)
      real(dp), intent(in) :: src(0:*)
      integer, intent(in) :: qstart, nq, nn
      real(dp), intent(inout) :: out(0:,0:)
      integer :: xx, kq
      integer(i8) :: base, nsq
      nsq = int(nn, i8) * int(nn, i8)
      !$omp parallel do default(shared) private(kq, base) schedule(static) &
      !$omp   if (int(nq, i8)*int(nn, i8)**3 >= omp_min_words)
      do xx = 0, nn - 1
        do kq = 0, nq - 1
          base = int(xx, i8)*nsq*int(nn, i8) + int(qstart + kq, i8)*nsq
          out(:, kq + nq*xx) = src(base:base + nsq - 1_i8)
        end do
      end do
    end subroutine gather_q_second

    !> out((b,d), kk + nq*x) = src[x, b, qs+kk, d]: the B-term operand over a
    !> contiguous q range.  Batching q turns the per-pair nbf-long run into an
    !> nq*nbf-long one, which is where the gather cost goes.
    subroutine gather_q_third(src, qstart, nq, nn, out)
      real(dp), intent(in) :: src(0:*)
      integer, intent(in) :: qstart, nq, nn
      real(dp), intent(inout) :: out(0:,0:)
      integer :: xx, bb, kq, dd
      integer(i8) :: base
      !$omp parallel do default(shared) private(bb, kq, dd, base) schedule(static) &
      !$omp   if (int(nq, i8)*int(nn, i8)**3 >= omp_min_words)
      do xx = 0, nn - 1
        do bb = 0, nn - 1
          base = int(xx, i8)*int(nn, i8)**3 + int(bb, i8)*int(nn, i8)**2 &
               + int(qstart, i8)*int(nn, i8)
          do kq = 0, nq - 1
            do dd = 0, nn - 1
              out(bb*nn + dd, kq + nq*xx) = &
                  src(base + int(kq, i8)*int(nn, i8) + int(dd, i8))
            end do
          end do
        end do
      end do
    end subroutine gather_q_third

    !> out((b,c), kk + nq*x) = src[x, b, c, qs+kk]: the C-term operand over a
    !> contiguous q range.
    subroutine gather_q_fourth(src, qstart, nq, nn, out)
      real(dp), intent(in) :: src(0:*)
      integer, intent(in) :: qstart, nq, nn
      real(dp), intent(inout) :: out(0:,0:)
      integer :: xx, bb, kq, cc
      integer(i8) :: base
      !$omp parallel do default(shared) private(bb, kq, cc, base) schedule(static) &
      !$omp   if (int(nq, i8)*int(nn, i8)**3 >= omp_min_words)
      do xx = 0, nn - 1
        do bb = 0, nn - 1
          base = int(xx, i8)*int(nn, i8)**3 + int(bb, i8)*int(nn, i8)**2 &
               + int(qstart, i8)
          do cc = 0, nn - 1
            do kq = 0, nq - 1
              out(bb*nn + cc, kq + nq*xx) = &
                  src(base + int(cc, i8)*int(nn, i8) + int(kq, i8))
            end do
          end do
        end do
      end do
    end subroutine gather_q_fourth

    !> One element of the one-index derivative tensor T[a,b,c,d] for the pair
    !> (P,Q), evaluated from the eight slabs instead of a stored tensor.
    pure function tval(g, nn, pp, qq, a, b, c, d) result(val)
      real(dp), intent(in) :: g(0:*)
      integer, intent(in) :: nn, pp, qq, a, b, c, d
      real(dp) :: val
      val = 0.0_dp
      if (a == qq) val = val + gel(g, nn, pp, b, c, d)
      if (a == pp) val = val - gel(g, nn, qq, b, c, d)
      if (b == qq) val = val + gel(g, nn, a, pp, c, d)
      if (b == pp) val = val - gel(g, nn, a, qq, c, d)
      if (c == qq) val = val + gel(g, nn, a, b, pp, d)
      if (c == pp) val = val - gel(g, nn, a, b, qq, d)
      if (d == qq) val = val + gel(g, nn, a, b, c, pp)
      if (d == pp) val = val - gel(g, nn, a, b, c, qq)
    end function tval

    !> C-order element g[a,b,c,d] of an [n,n,n,n] buffer.
    pure function gel(g, nn, a, b, c, d) result(val)
      real(dp), intent(in) :: g(0:*)
      integer, intent(in) :: nn, a, b, c, d
      real(dp) :: val
      val = g(((int(a, i8)*int(nn, i8) + int(b, i8))*int(nn, i8) &
               + int(c, i8))*int(nn, i8) + int(d, i8))
    end function gel

  end subroutine casscf_hess_bmat

end module casscf_hess_bmat_mod
