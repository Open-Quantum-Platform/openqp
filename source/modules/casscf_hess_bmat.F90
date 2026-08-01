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
!> is ever formed.
!>
!> ---------------------------------------------------------------------------
!> WHAT THIS LOOP COSTS, measured rather than modelled (OQP_HESS_PROF=1 prints
!> the split; see `ptock` below).  Of the form this file replaced, on a
!> 152-core Xeon 8368, single thread:
!>
!>     nbf=58 ncore=2 nact=6 npar=412      nbf=76 ncore=2 nact=6 npar=556
!>       run_gemm    5057.5 ms  93.89 %      run_gemm   20098.1 ms  94.45 %
!>       p_gather     139.9 ms   2.60 %      p_gather     543.9 ms   2.56 %
!>       ymat_gemm    125.8 ms   2.33 %      ymat_gemm    475.0 ms   2.23 %
!>       q_gather      24.3 ms   0.45 %      q_gather      54.9 ms   0.26 %
!>       fp_z_gemm     30.0 ms   0.56 %      fp_z_gemm     87.3 ms   0.41 %
!>       everything else        0.17 %      everything else        0.09 %
!>
!> "Everything else" is the whole scalar tail: the tval/gel assembly of g_der
!> (O(nact^4) per pair), the zmat assembly (O(nbf^2)), the B extraction
!> (O(npar)) and the resp/resm difference.  At nbf=58 they are 0.03 %, 0.06 %,
!> 0.02 % and 0.03 % -- one part in nine hundred between them, and they get
!> *smaller* with nbf.  This kernel is GEMM, and only GEMM.  (An earlier
!> comment here claimed 80% of the routine was in its gathers.  That was true
!> of the per-pair gathers this file replaced; it has not been true since they
!> were batched, and it sent a later reader looking in the wrong place.)
!>
!> So the shape of those GEMMs is the whole performance story, and it is the
!> reason the routine used to be flat in the thread count.  Both indices are
!> batched, which is what the pair list already comes in:
!>
!>   * q.  _nonredundant_pairs emits p outermost and sweeps a whole space in
!>     q -- (active, inactive), (virtual, inactive), (virtual, active) -- so
!>     the pairs sharing a p are consecutive AND their q is a contiguous
!>     ascending range.
!>   * p.  Each of those three spaces repeats the SAME q-range for every p in
!>     a contiguous ascending range, so the pair list is three p-range x
!>     q-range rectangles -- and the first two merge, because (active,
!>     inactive) and (virtual, inactive) share both their q-range and a
!>     contiguous run of p.  A whole rectangle tile is done in one GEMM.
!>
!> Batching q alone fixes the k-dominance: per pair the contraction is
!> (nbf x nbf^2).(nbf^2 x nbf), an inner-product shape whose k dominates both
!> free dimensions.  Batching p as well is what makes it thread: the free
!> dimension M goes from nbf to npc*nbf, and a GEMM with M=58 cannot fill more
!> than a handful of cores no matter how many are present.  Measured on the
!> exact operands and leading dimensions at nbf=58 (32 threads, AVX-512
!> kernels): M=58 gives 226 GFlop/s and M=21*58 gives 778 GFlop/s for the same
!> flops per pair.
!>
!> Both operands of every GEMM are therefore gathered over a range, and the
!> gather is cheaper per element for it: the B operand becomes one nq*nbf-long
!> run per stride instead of nq nbf-long ones, and the C operand becomes
!> nq-long runs instead of single elements.  The A operand needs no strided
!> read at all, only a block copy.  Twelve operand buffers in total, six keyed
!> on the p range and six on the q range.  Both keys are remembered, so a range
!> that repeats across rectangles is gathered once: at the CAS(6,6)/cc-pVTZ
!> shape a whole call does two p-gather rounds and two q-gather rounds, not one
!> per pair.  That keying is why the p buffers cost six gathers rather than the
!> four the q-only form needed and still come out ahead.
!>
!> The two families differ in which operand carries q, so they accumulate into
!> two buffers with the tile indices on opposite sides:
!>
!>     resp[(kp,m), (kq,n)] = A[P,Q]+B[P,Q]+C[P,Q]
!>     resm[(kq,m), (kp,n)] = A[Q,P]+B[Q,P]+C[Q,P]
!>     C1[m,n]              = resp[(kp,m),(kq,n)] - resm[(kq,m),(kp,n)]
!>
!> Nothing here ASSUMES that ordering, it only exploits it: the rectangles are
!> scanned off the pair list and verified entry by entry, and a caller whose
!> list is ordered differently gets smaller tiles, down to one pair per tile.
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
  !> Swept on the same box at nbf=58: 16 MiB 267 ms, 64 MiB 165, 256 MiB 164,
  !> 1 GiB 166 -- flat once it holds the whole q-run, which 256 MiB does up to
  !> nact=12 at nbf=76.  Kept.
  integer(i8), parameter :: q_budget = 33554432_i8
  !> Budget for the six p-batched operand buffers together, in doubles
  !> (512 MiB).  This one sets the GEMM's free dimension M = npc*nbf, so it is
  !> a throughput knob and not merely a gather-amortisation knob, which is why
  !> it is the larger of the two.  Re-derived on a 152-core Xeon 8368 at 32
  !> threads, sweeping it against the whole routine:
  !>
  !>     nbf=58    16 MiB 834 ms | 64 MiB 233 | 128 MiB 203 | 256 MiB 187
  !>             | 512 MiB 167 | 1 GiB 167 | 2 GiB 169
  !>     nbf=76    64 MiB 1082 ms | 128 MiB 855 | 256 MiB 683 | 512 MiB 580
  !>             | 1 GiB 588 | 2 GiB 539
  !>
  !> 128 MiB -- the old value, calibrated on a laptop at nbf=36 and nbf=48 when
  !> the GEMM still took one p at a time -- costs 1.22x at nbf=58 and 1.47x at
  !> nbf=76.  512 MiB holds the widest rectangle whole at nbf=58 and is within
  !> 8% of unbounded at nbf=76, while staying a fixed cap rather than growing
  !> as 6*nbf^4 (six times the ERI array) the way "hold everything" would.
  integer(i8), parameter :: p_budget = 67108864_i8

  !> Threading threshold for a gather, in doubles moved.  The gathers are the
  !> only part of this routine BLAS is not already threading, and they are
  !> latency-bound, so spreading them wins ~7x at nbf=58 on 32 threads.  Below
  !> this the parallel region costs more than the gather does: at nbf=24 the
  !> whole routine is 3.5 ms and threading every gather made it 0.85x.
  integer(i8), parameter :: omp_min_words = 262144_i8

  public :: casscf_hess_bmat

  ! -------------------------------------------------------------------------
  ! Section profiler, off unless OQP_HESS_PROF=1 is set in the environment.
  ! This routine is 94% one section; that was not visible from the outside and
  ! cost a round of tuning aimed at the wrong three loops, so the split stays
  ! available.  When off, each probe is one predicted branch on a saved
  ! logical.  Not thread-safe, and not meant to be: it profiles a serial
  ! driver calling a threaded kernel.
  ! -------------------------------------------------------------------------
  integer, parameter :: nsec = 12
  character(len=11), parameter :: sec_name(nsec) = [ character(len=11) :: &
      'setup', 'ymat_gemm', 'p_gather', 'q_gather', 'run_gemm', 'cmat_diff', &
      'fp_slab1_t', 'fp_z_gemm', 'fp_z_asm', 'fp_bextract', 'fp_gder', 'fp_fder' ]
  real(dp), save :: prof_t(nsec) = 0.0_dp
  integer(i8), save :: prof_rate = 1_i8, prof_mark = 0_i8
  logical, save :: prof_on = .false., prof_ready = .false.

contains

  !> Charge the time since the previous probe to section `k`.
  subroutine ptock(k)
    integer, intent(in) :: k
    integer(i8) :: now
    if (.not. prof_on) return
    call system_clock(now)
    prof_t(k) = prof_t(k) + real(now - prof_mark, dp) / real(prof_rate, dp)
    prof_mark = now
  end subroutine ptock

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
    integer :: n, nc, na, np, n2, l, p, q0, i, j, m, ks
    integer(i8) :: n3
    real(dp), allocatable :: ymat(:,:), cmat(:,:), tmat(:,:), w1(:,:), w2(:,:)
    real(dp), allocatable :: dm(:,:), zmat(:,:)
    ! Operand buffers.  p?g / p?r are the ERI and 2-RDM gathered over a chunk
    ! of consecutive p, q?g / q?r the same over a range of q; the a/b/c suffix
    ! is which of the three interior slots the batched index occupies.
    real(dp), allocatable :: pag(:,:), pbg(:,:), pcg(:,:)
    real(dp), allocatable :: par(:,:), pbr(:,:), pcr(:,:)
    real(dp), allocatable :: qag(:,:), qbg(:,:), qcg(:,:)
    real(dp), allocatable :: qar(:,:), qbr(:,:), qcr(:,:)
    real(dp), allocatable :: resp(:), resm(:)
    integer :: nrun, nprow, maxrun, maxprow, nqmax, npmax
    integer :: ip0, iq0, npc, nqc, pc0, qs, kp, kq, ldp, ldq
    integer :: last_ps, last_np, last_qs, last_nq
    real(dp) :: tot

    if (.not. prof_ready) then
      block
        character(len=8) :: envv
        integer :: est
        call get_environment_variable("OQP_HESS_PROF", envv, status=est)
        prof_on = (est == 0 .and. envv(1:1) == '1')
      end block
      call system_clock(count_rate=prof_rate)
      prof_ready = .true.
    end if
    if (prof_on) then
      prof_t = 0.0_dp
      call system_clock(prof_mark)
    end if

    n = int(nbf)
    nc = int(ncore)
    na = int(nact)
    np = int(npar)
    n2 = n * n
    n3 = int(n, i8) * int(n, i8) * int(n, i8)
    if (n <= 0 .or. np <= 0) return

    ! Widest rectangle in the pair list, which is what the batched buffers have
    ! to hold.  Scanned rather than assumed: a caller whose pair list is ordered
    ! differently simply gets smaller tiles, down to one pair at (1,1).
    maxrun = 1
    maxprow = 1
    l = 0
    do while (l < np)
      call scan_block(l, p, q0, nrun, nprow)
      maxrun = max(maxrun, nrun)
      maxprow = max(maxprow, nprow)
      l = l + nprow * nrun
    end do
    nqmax = int(max(1_i8, min(int(maxrun, i8), q_budget / (6_i8 * n3))))
    npmax = int(max(1_i8, min(int(maxprow, i8), p_budget / (6_i8 * n3))))

    allocate(ymat(0:n-1, 0:n-1), cmat(0:n-1, 0:n-1), tmat(0:n-1, 0:n-1))
    allocate(w1(0:n-1, 0:n-1), w2(0:n-1, 0:n-1), zmat(0:n-1, 0:n-1))
    allocate(dm(0:n-1, 0:n-1))
    allocate(pag(0:n2-1, 0:npmax*n-1), pbg(0:n2-1, 0:npmax*n-1))
    allocate(pcg(0:n2-1, 0:npmax*n-1), par(0:n2-1, 0:npmax*n-1))
    allocate(pbr(0:n2-1, 0:npmax*n-1), pcr(0:n2-1, 0:npmax*n-1))
    allocate(qag(0:n2-1, 0:nqmax*n-1), qbg(0:n2-1, 0:nqmax*n-1))
    allocate(qcg(0:n2-1, 0:nqmax*n-1), qar(0:n2-1, 0:nqmax*n-1))
    allocate(qbr(0:n2-1, 0:nqmax*n-1), qcr(0:n2-1, 0:nqmax*n-1))
    allocate(resp(0:int(npmax, i8)*int(nqmax, i8)*int(n2, i8) - 1_i8))
    allocate(resm(0:int(npmax, i8)*int(nqmax, i8)*int(n2, i8) - 1_i8))
    last_ps = -1
    last_np = -1
    last_qs = -1
    last_nq = -1

    do i = 0, n - 1
      do j = 0, n - 1
        dm(i, j) = dmat(int(i, i8)*int(n, i8) + int(j, i8))
      end do
    end do
    call ptock(1)

    ! ---- Y[n,a] = sum_bcd G[n,b,c,d] g[a,b,c,d], once for all pairs.
    ! Both buffers are already column-major ((b,c,d), leading index) matrices
    ! with leading dimension nbf^3, so this is one GEMM with no repacking.
    ! ymat(a,n) holds Y[n,a].
    call dgemm('T', 'N', n, n, int(n3), 1.0_dp, eri, int(n3), &
               rdm2, int(n3), 0.0_dp, ymat, n)
    call ptock(2)

    l = 0
    do while (l < np)
      call scan_block(l, p, q0, nrun, nprow)

      ip0 = 0
      do while (ip0 < nprow)
        npc = min(npmax, nprow - ip0)
        pc0 = p + ip0
        ldp = npc * n

        ! The six p-indexed operands, keyed on the range so that the two
        ! rectangles sharing a p range gather it once.
        if (pc0 /= last_ps .or. npc /= last_np) then
          call gather_q_second(eri,  pc0, npc, n, pag)
          call gather_q_third (eri,  pc0, npc, n, pbg)
          call gather_q_fourth(eri,  pc0, npc, n, pcg)
          call gather_q_second(rdm2, pc0, npc, n, par)
          call gather_q_third (rdm2, pc0, npc, n, pbr)
          call gather_q_fourth(rdm2, pc0, npc, n, pcr)
          last_ps = pc0
          last_np = npc
          call ptock(3)
        end if

        iq0 = 0
        do while (iq0 < nrun)
          nqc = min(nqmax, nrun - iq0)
          qs = q0 + iq0
          ldq = nqc * n

          ! The six q-indexed operands, likewise keyed on the range.
          if (qs /= last_qs .or. nqc /= last_nq) then
            call gather_q_second(eri,  qs, nqc, n, qag)
            call gather_q_third (eri,  qs, nqc, n, qbg)
            call gather_q_fourth(eri,  qs, nqc, n, qcg)
            call gather_q_second(rdm2, qs, nqc, n, qar)
            call gather_q_third (rdm2, qs, nqc, n, qbr)
            call gather_q_fourth(rdm2, qs, nqc, n, qcr)
            last_qs = qs
            last_nq = nqc
            call ptock(4)
          end if

          ! resp((kp,m), (kq,n)) = (A+B+C)[P,Q][m,n] for every pair of the
          ! tile at once: p rides the left operand's free index and q the
          ! right operand's, so one GEMM covers npc*nqc pairs.
          call dgemm('T', 'N', ldp, ldq, n2, 1.0_dp, pag, n2, qar, n2, 0.0_dp, resp, ldp)
          call dgemm('T', 'N', ldp, ldq, n2, 1.0_dp, pbg, n2, qbr, n2, 1.0_dp, resp, ldp)
          call dgemm('T', 'N', ldp, ldq, n2, 1.0_dp, pcg, n2, qcr, n2, 1.0_dp, resp, ldp)

          ! resm((kq,m), (kp,n)) = (A+B+C)[Q,P][m,n]: the same tile with the
          ! roles of the two ranges exchanged.
          call dgemm('T', 'N', ldq, ldp, n2, 1.0_dp, qag, n2, par, n2, 0.0_dp, resm, ldq)
          call dgemm('T', 'N', ldq, ldp, n2, 1.0_dp, qbg, n2, pbr, n2, 1.0_dp, resm, ldq)
          call dgemm('T', 'N', ldq, ldp, n2, 1.0_dp, qcg, n2, pcr, n2, 1.0_dp, resm, ldq)
          call ptock(5)

          do kp = 0, npc - 1
            do kq = 0, nqc - 1
              do j = 0, n - 1
                do m = 0, n - 1
                  cmat(m, j) = &
                      resp(int(kp + npc*m, i8) + int(ldp, i8)*int(kq + nqc*j, i8)) &
                    - resm(int(kq + nqc*m, i8) + int(ldq, i8)*int(kp + npc*j, i8))
                end do
              end do
              call ptock(6)
              call finish_pair(l + (ip0 + kp)*nrun + iq0 + kq, pc0 + kp, qs + kq)
            end do
          end do

          iq0 = iq0 + nqc
        end do
        ip0 = ip0 + npc
      end do

      l = l + nprow * nrun
    end do

    if (prof_on) then
      tot = sum(prof_t)
      write(0, '(a)') '--- casscf_hess_bmat sections ---'
      write(0, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'nbf=', n, ' ncore=', nc, &
          ' nact=', na, ' npar=', np, ' npc=', npmax, ' nqc=', nqmax
      do ks = 1, nsec
        write(0, '(2x,a11,f10.1,a,f7.2,a)') sec_name(ks), prof_t(ks)*1.0e3_dp, &
            ' ms  ', 100.0_dp*prof_t(ks)/max(tot, 1.0e-30_dp), ' %'
      end do
      write(0, '(2x,a11,f10.1,a)') 'TOTAL', tot*1.0e3_dp, ' ms'
    end if

    deallocate(ymat, cmat, tmat, w1, w2, zmat, dm)
    deallocate(pag, pbg, pcg, par, pbr, pcr)
    deallocate(qag, qbg, qcg, qar, qbr, qcr, resp, resm)

  contains

    !> The maximal p-range x q-range rectangle of the pair list starting at
    !> `l0`: `npr` consecutive p, each carrying the same ascending run of `nr`
    !> values of q from `qq0`, laid out p-outermost.  Every entry of every row
    !> is verified, so a list that is not in that order simply yields a smaller
    !> rectangle -- at worst 1x1, which is the per-pair form.
    subroutine scan_block(l0, pp, qq0, nr, npr)
      integer, intent(in) :: l0
      integer, intent(out) :: pp, qq0, nr, npr
      integer :: jj, kk, base
      logical :: ok
      pp = int(pairs(2*l0))
      qq0 = int(pairs(2*l0 + 1))
      nr = 1
      do while (l0 + nr < np)
        if (int(pairs(2*(l0 + nr))) /= pp) exit
        if (int(pairs(2*(l0 + nr) + 1)) /= qq0 + nr) exit
        nr = nr + 1
      end do
      npr = 1
      do while (l0 + (npr + 1)*nr <= np)
        base = l0 + npr*nr
        ok = .true.
        do jj = 0, nr - 1
          kk = base + jj
          if (int(pairs(2*kk)) /= pp + npr .or. &
              int(pairs(2*kk + 1)) /= qq0 + jj) then
            ok = .false.
            exit
          end if
        end do
        if (.not. ok) exit
        npr = npr + 1
      end do
    end subroutine scan_block

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
      call ptock(7)

      ! ---- Z[m,n] = (D t - t D)[n,m] + 2 C1[m,n]
      call dgemm('N', 'N', n, n, n, 1.0_dp, dm, n, tmat, n, 0.0_dp, w1, n)
      call dgemm('N', 'N', n, n, n, 1.0_dp, tmat, n, dm, n, 0.0_dp, w2, n)
      call ptock(8)
      do mm = 0, n - 1
        do jj = 0, n - 1
          zmat(mm, jj) = (w1(jj, mm) - w2(jj, mm)) + 2.0_dp * cmat(mm, jj)
        end do
      end do
      call ptock(9)

      ! ---- B[k,l] = Z[p_k,q_k] - Z[q_k,p_k]
      do kk2 = 0, np - 1
        ii = int(pairs(2*kk2))
        jj = int(pairs(2*kk2 + 1))
        bmat(int(kk2, i8)*int(np, i8) + int(lp, i8)) = zmat(ii, jj) - zmat(jj, ii)
      end do
      call ptock(10)

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
      call ptock(11)
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
      call ptock(12)
    end subroutine finish_pair

    !> out((c,d), kk + nq*x) = src[x, qs+kk, c, d] for kk = 0..nq-1: the
    !> A-term operand over a contiguous range.  Each (x, q) slab is nbf^2
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
    !> contiguous range.  Batching turns the per-pair nbf-long run into an
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
    !> contiguous range.
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
