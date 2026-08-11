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
!> *smaller* with nbf.  On one thread this kernel is GEMM, and only GEMM.
!> (An earlier comment here claimed 80% of the routine was in its gathers.
!> That was true of the per-pair gathers this file replaced; it has not been
!> true since they were batched, and it sent a later reader looking in the
!> wrong place.)
!>
!> Read that table for what to OPTIMISE and it is right; read it for what to
!> THREAD and it is a trap.  A section worth 0.5 % on one thread is worth 13 %
!> on thirty-two if it does not thread, and that is exactly what the tail did.
!> The thread-count table further down is the one to use for scaling work.
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
!> ---------------------------------------------------------------------------
!> WHERE THE THREADS GO.  Batching the tiles made the routine fast; it did not
!> by itself make it scale, and the reason is worth writing down because it is
!> not the one a reader expects.  With the tiles batched and everything handed
!> to a threaded BLAS, the routine was flat above 8 cores and had these
!> sections at nbf=58 on 32 threads (168 ms total):
!>
!>     run_gemm   101.0 ms  60 %     the six tile GEMMs
!>     p_gather    20.2 ms  12 %
!>     ymat_gemm   18.8 ms  11 %
!>     pair tail   22.5 ms  13 %     seven small serial sections
!>     q_gather     5.3 ms   3 %
!>
!> Three separate causes, fixed separately and each measured on its own build.
!> Kernel wall time at nbf=58 ncore=2 nact=6 npar=412 (and at nbf=76 npar=556),
!> shipped ILP64 OpenMP OpenBLAS 0.3.15, OPENBLAS_CORETYPE=SkylakeX:
!>
!>     threads          1     8    16    32    48    64   128
!>     nbf=58
!>       batched-only 1028   252   184   164   170   187   165  ms
!>       + ymat split 1035   235   172   156   156   178   152
!>       + pair tail  1016   218   151   131   134   163   136
!>       + gemm split 1025   175   111    87    82   111   114
!>     nbf=76
!>       batched-only 3740   852   652   595   532   668   635  ms
!>       + gemm split 3714   645   393   316   288   368   327
!>
!>   1. ymat_gemm was M = N = nbf with K = nbf^3.  No BLAS threads a shape
!>      that narrow -- 34.0 GFlop/s on 1 thread and 77.4 on 32, measured on
!>      the shipped library at exactly this shape.  Its K is now split by hand.
!>      See the comment on that call.
!>
!>   2. The per-pair tail ran serially while everything around it threaded, so
!>      it went from 0.5 % of the routine at 1 thread to 13 % at 32.  Its
!>      pairs are independent and write disjoint output, so it is now an
!>      OpenMP loop with per-thread scratch.  Nothing is reduced and no
!>      arithmetic order changes.
!>
!>   3. The tile GEMMs themselves.  OpenBLAS 0.3.15 stops scaling on this box
!>      around 32 threads and it is NOT this kernel's shapes that do it -- a
!>      4000x4000x4000 square GEMM measures 1365 GFlop/s at 32 threads, 1444
!>      at 76 and 1605 at 128, for a machine whose AVX-512 peak is several
!>      times that.  Cutting the larger free dimension across threads by hand
!>      and letting each chunk call a self-serialised BLAS beats it by 1.2x to
!>      2.7x depending on the tile orientation and the thread count.  See
!>      `gemm_tn_split`.
!>
!> The nesting question is decided per call site, not globally.  OpenMP around
!> a LARGE BLAS call in this module was tried once and measured at 0.67x, and
!> that result stands: what is wrapped here is either a call BLAS cannot
!> thread (ymat, the pair tail's nbf^3 GEMMs) or a call BLAS threads badly
!> (the tile GEMMs, where the wrapping replaces BLAS's threading rather than
!> nesting inside it).  There is no second OpenMP runtime to worry about --
!> liboqp and libopenblaso64 share libgomp -- and OpenBLAS serialises itself
!> inside a parallel region, so no region here oversubscribes.
!>
!> WHERE IT SATURATES NOW, and what is next.  On this 2 x 38-core Xeon 8368
!> (152 logical CPUs are 76 physical cores plus SMT), both shapes peak at 48
!> threads and lose ~25 % by 64.  Above 48 nothing in this file improves: the
!> tile GEMMs stop gaining, the gathers get slower, and the pair tail's fork
!> and barrier costs grow against a fixed 412 or 556 pairs.  The three limits
!> in front of a further gain, in order:
!>
!>   * p_gather, 20 % at nbf=76 and nearly flat in the thread count.  It moves
!>     ~1.5 GB in strided nbf-long runs at about 21 GB/s, an order below what
!>     the sockets can do.  It is the largest single non-GEMM item left.
!>   * DGEMM throughput.  Even split by hand the tile GEMMs reach ~1.2-1.6
!>     TFlop/s.  A newer OpenBLAS (0.3.15 predates Ice Lake-SP; its cpu probe
!>     does not even recognise this part without OPENBLAS_CORETYPE) or a
!>     different BLAS is the lever here, not this file.
!>   * Problem size.  412 pairs over 76 cores is 5 pairs a core.  The tile
!>     rectangle is also capped by p_budget: at nbf=76 npc is 25, so M is
!>     1900 rather than the 4180 the pair list would allow.
!>
!> ---------------------------------------------------------------------------
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
#ifdef _OPENMP
  use omp_lib, only: omp_get_max_threads
#endif
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

  !> Threading threshold for the per-pair tail, in doubles of arithmetic
  !> (npc*nqc*nbf^3, which is what its two nbf^3 GEMMs cost).  Same reasoning
  !> as `omp_min_words` and the same order of magnitude; below it the tail is
  !> already under a millisecond and the fork costs more than it saves.
  integer(i8), parameter :: omp_min_pairwork = 262144_i8

  !> Below this many multiply-adds a tile GEMM is handed to BLAS whole rather
  !> than split by `gemm_tn_split`; there is nothing to spread.
  integer(i8), parameter :: split_min_madds = 16777216_i8
  !> Smallest free-dimension chunk `gemm_tn_split` will cut.  The sweep that
  !> set it is in that routine's comment: at M=3248 the optimum is one chunk
  !> per thread down to ~25 rows, and 12 rows is already 30% off the peak.
  integer, parameter :: split_min_rows = 16

  public :: casscf_hess_bmat

  ! -------------------------------------------------------------------------
  ! Section profiler, off unless OQP_HESS_PROF=1 is set in the environment.
  ! This routine is nearly all GEMM; that was not visible from the outside and
  ! cost a round of tuning aimed at the wrong three loops, so the split stays
  ! available.  When off, each probe is one predicted branch on a saved
  ! logical.
  !
  ! Sections 1..6 are WALL clock on the master thread and add up to TOTAL.
  ! Sections 7..13 sit inside the threaded pair tail, so they are accumulated
  ! per thread and reported as THREAD-seconds -- a breakdown of the work in
  ! section 6 (`pair_tail`), not an addition to it.  Their own probe is
  ! `ptockp`, whose mark is threadprivate; mixing the two would have made the
  ! tail's numbers meaningless the moment it was parallelised.
  ! -------------------------------------------------------------------------
  integer, parameter :: nsec = 13
  integer, parameter :: nwall = 6
  character(len=11), parameter :: sec_name(nsec) = [ character(len=11) :: &
      'setup', 'ymat_gemm', 'p_gather', 'q_gather', 'run_gemm', 'pair_tail', &
      'cmat_diff', 'fp_slab1_t', 'fp_z_gemm', 'fp_z_asm', 'fp_bextract', &
      'fp_gder', 'fp_fder' ]
  real(dp), save :: prof_t(nsec) = 0.0_dp
  integer(i8), save :: prof_rate = 1_i8, prof_mark = 0_i8
  logical, save :: prof_on = .false., prof_ready = .false.
  integer(i8), save :: prof_markp = 0_i8
  !$omp threadprivate(prof_markp)

contains

  !> Charge the wall time since the previous probe to section `k` (1..nwall).
  subroutine ptock(k)
    integer, intent(in) :: k
    integer(i8) :: now
    if (.not. prof_on) return
    call system_clock(now)
    prof_t(k) = prof_t(k) + real(now - prof_mark, dp) / real(prof_rate, dp)
    prof_mark = now
  end subroutine ptock

  !> Start this thread's clock for the per-pair sections.
  subroutine ptickp()
    if (.not. prof_on) return
    call system_clock(prof_markp)
  end subroutine ptickp

  !> Charge this thread's elapsed time to section `k` (> nwall).  Called from
  !> inside the threaded pair loop, so the accumulator update is atomic and
  !> the result is thread-seconds.
  subroutine ptockp(k)
    integer, intent(in) :: k
    integer(i8) :: now
    real(dp) :: dt
    if (.not. prof_on) return
    call system_clock(now)
    dt = real(now - prof_markp, dp) / real(prof_rate, dp)
    !$omp atomic update
    prof_t(k) = prof_t(k) + dt
    prof_markp = now
  end subroutine ptockp

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
    real(dp), allocatable :: ymat(:,:), ypart(:,:,:), dm(:,:)
    ! Operand buffers.  p?g / p?r are the ERI and 2-RDM gathered over a chunk
    ! of consecutive p, q?g / q?r the same over a range of q; the a/b/c suffix
    ! is which of the three interior slots the batched index occupies.
    real(dp), allocatable :: pag(:,:), pbg(:,:), pcg(:,:)
    real(dp), allocatable :: par(:,:), pbr(:,:), pcr(:,:)
    real(dp), allocatable :: qag(:,:), qbg(:,:), qcg(:,:)
    real(dp), allocatable :: qar(:,:), qbr(:,:), qcr(:,:)
    real(dp), allocatable :: resp(:), resm(:)
    integer :: nrun, nprow, maxrun, maxprow, nqmax, npmax, nsplit
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
#ifdef _OPENMP
    nsplit = omp_get_max_threads()
#else
    nsplit = 1
#endif
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

    allocate(ymat(0:n-1, 0:n-1), ypart(0:n-1, 0:n-1, 0:n-1))
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
    ! with leading dimension nbf^3, so this needs no repacking -- but as a
    ! single GEMM it is M = N = nbf with K = nbf^3, and no BLAS threads a shape
    ! that narrow: the free dimensions are all it has to cut, and there are
    ! only nbf of each.  Measured on the shipped ILP64 OpenBLAS at nbf=58
    ! (M=N=58, K=195112): 34.0 GFlop/s on 1 thread, 77.4 on 32 -- 2.3x for 32x
    ! the cores, which left it as 11% of the routine once everything around it
    ! had been fixed.
    !
    ! So the K is split here instead, into `n` chunks of nbf^2 rows each, one
    ! GEMM per chunk into its own nbf x nbf partial, and the partials summed.
    ! Each chunk is M = N = nbf with K = nbf^2 -- still narrow, but there are
    ! now nbf of them to hand out, and BLAS serialises itself inside the
    ! parallel region so nothing is oversubscribed.
    !
    ! The chunk count is nbf, NOT the thread count.  That is deliberate: the
    ! summation order over partials is then a property of the problem and not
    ! of OMP_NUM_THREADS, so the answer is bit-identical from 1 thread to 128.
    ! It is not bit-identical to the single un-split GEMM -- a different K
    ! blocking never is -- and it does not have to be: ymat reaches only B,
    ! whose standard is 1e-15 relative.  f_der and g_der are built from `eri`
    ! and `h1e` directly and never see this array, which is why the bit-exact
    ! pin on g_der is untouched by anything in this file's GEMM path.
    ! ymat(a,n) holds Y[n,a].
    !$omp parallel do default(shared) private(i) schedule(static) &
    !$omp   if (int(n, i8)**4 >= omp_min_words)
    do i = 0, n - 1
      call dgemm('T', 'N', n, n, n2, 1.0_dp, &
                 eri(int(i, i8)*int(n2, i8)), int(n3), &
                 rdm2(int(i, i8)*int(n2, i8)), int(n3), &
                 0.0_dp, ypart(0, 0, i), n)
    end do
    !$omp parallel do default(shared) private(i, j, ks) schedule(static) &
    !$omp   if (int(n, i8)**3 >= omp_min_words)
    do j = 0, n - 1
      do i = 0, n - 1
        ymat(i, j) = ypart(i, j, 0)
        do ks = 1, n - 1
          ymat(i, j) = ymat(i, j) + ypart(i, j, ks)
        end do
      end do
    end do
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
          call gemm_tn_split(ldp, ldq, n2, pag, n2, qar, n2, 0.0_dp, resp, ldp)
          call gemm_tn_split(ldp, ldq, n2, pbg, n2, qbr, n2, 1.0_dp, resp, ldp)
          call gemm_tn_split(ldp, ldq, n2, pcg, n2, qcr, n2, 1.0_dp, resp, ldp)

          ! resm((kq,m), (kp,n)) = (A+B+C)[Q,P][m,n]: the same tile with the
          ! roles of the two ranges exchanged.
          call gemm_tn_split(ldq, ldp, n2, qag, n2, par, n2, 0.0_dp, resm, ldq)
          call gemm_tn_split(ldq, ldp, n2, qbg, n2, pbr, n2, 1.0_dp, resm, ldq)
          call gemm_tn_split(ldq, ldp, n2, qcg, n2, pcr, n2, 1.0_dp, resm, ldq)
          call ptock(5)

          ! Everything left for this tile is per-pair and independent: pair
          ! (kp,kq) reads its own window of resp/resm and writes column
          ! l+(ip0+kp)*nrun+iq0+kq of B and slice l+... of f_der and g_der.
          ! No two pairs touch the same output element and nothing is reduced,
          ! so this loop parallelises without changing a single arithmetic
          ! order -- the bit-exact g_der pin survives it by construction.
          !
          ! The scratch is what has to move.  cmat/tmat/w1/w2/zmat were host
          ! -associated nbf x nbf buffers shared by every pair; they are now
          ! allocated inside the parallel region, which makes them one set per
          ! thread, and passed to `finish_pair` explicitly.
          !
          ! The two nbf^3 GEMMs in the tail run one-thread-per-pair rather
          ! than BLAS-threaded-per-GEMM.  That is the right nesting for them:
          ! at nbf=58 each is 2*58^3 = 0.4 MFlop, far too small for BLAS to
          ! spread, and the OpenMP OpenBLAS serialises itself inside a
          ! parallel region anyway, so there is no second runtime and no
          ! oversubscription.  (Wrapping OpenMP around a *large* BLAS call in
          ! this kernel is the opposite trade and was measured at 0.67x.)
          !$omp parallel default(shared) private(kp, kq, j, m) &
          !$omp   if (int(npc, i8)*int(nqc, i8)*int(n, i8)**3 >= omp_min_pairwork)
          block
            real(dp), allocatable :: cmat(:,:), tmat(:,:)
            real(dp), allocatable :: w1(:,:), w2(:,:), zmat(:,:)
            allocate(cmat(0:n-1, 0:n-1), tmat(0:n-1, 0:n-1), zmat(0:n-1, 0:n-1))
            allocate(w1(0:n-1, 0:n-1), w2(0:n-1, 0:n-1))
            call ptickp()
            !$omp do collapse(2) schedule(static)
            do kp = 0, npc - 1
              do kq = 0, nqc - 1
                do j = 0, n - 1
                  do m = 0, n - 1
                    cmat(m, j) = &
                        resp(int(kp + npc*m, i8) + int(ldp, i8)*int(kq + nqc*j, i8)) &
                      - resm(int(kq + nqc*m, i8) + int(ldq, i8)*int(kp + npc*j, i8))
                  end do
                end do
                call ptockp(7)
                call finish_pair(l + (ip0 + kp)*nrun + iq0 + kq, pc0 + kp, &
                                 qs + kq, cmat, tmat, w1, w2, zmat)
              end do
            end do
            !$omp end do
            deallocate(cmat, tmat, w1, w2, zmat)
          end block
          !$omp end parallel
          call ptock(6)

          iq0 = iq0 + nqc
        end do
        ip0 = ip0 + npc
      end do

      l = l + nprow * nrun
    end do

    if (prof_on) then
      tot = sum(prof_t(1:nwall))
      write(0, '(a)') '--- casscf_hess_bmat sections ---'
      ! nthr is what omp_get_max_threads() reported AT ENTRY, not what the
      ! environment asked for.  Those two have come apart before: a BLAS
      ! thread save/restore elsewhere in the library can move the process
      ! -wide OpenMP count, and without this field a whole-job thread scan
      ! that is silently pinned looks like a kernel that does not scale.
      write(0, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'nbf=', n, ' ncore=', nc, &
          ' nact=', na, ' npar=', np, ' npc=', npmax, ' nqc=', nqmax, &
          ' nthr=', nsplit
      do ks = 1, nwall
        write(0, '(2x,a11,f10.1,a,f7.2,a)') sec_name(ks), prof_t(ks)*1.0e3_dp, &
            ' ms  ', 100.0_dp*prof_t(ks)/max(tot, 1.0e-30_dp), ' %'
      end do
      write(0, '(2x,a11,f10.1,a)') 'TOTAL', tot*1.0e3_dp, ' ms'
      write(0, '(a)') '  (pair_tail, in thread-seconds:)'
      do ks = nwall + 1, nsec
        write(0, '(4x,a11,f10.1,a,f7.2,a)') sec_name(ks), prof_t(ks)*1.0e3_dp, &
            ' ms  ', 100.0_dp*prof_t(ks) &
                     / max(sum(prof_t(nwall+1:nsec)), 1.0e-30_dp), ' %'
      end do
    end if

    deallocate(ymat, ypart, dm)
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

    !> C := op(A)^T . B + beta C, cut across threads by hand instead of by
    !> BLAS: one chunk of the larger free dimension per thread, each chunk one
    !> BLAS call that then serialises itself.
    !>
    !> This exists because the shipped OpenBLAS stops scaling around 32
    !> threads on this box and it is not the shape's fault -- a 4000^3 square
    !> GEMM does the same.  Measured on the tile shapes at nbf=58, chunking by
    !> the number of threads against letting OpenBLAS thread the whole call:
    !>
    !>     M=3248 N=348  K=3364     32 thr  1007 ->  1247 GFlop/s   1.24x
    !>     (p-major, resp)         128 thr   861 ->  1618            1.88x
    !>     M=348  N=3248 K=3364     32 thr   615 ->  1657            2.69x
    !>     (q-major, resm)         128 thr   673 ->  1336            1.98x
    !>
    !> Which free dimension is cut matters as much as cutting one at all: on
    !> the q-major shape, splitting M=348 instead of N=3248 gives 320 GFlop/s
    !> at 32 threads, half of what BLAS manages unaided.  So the LARGER of the
    !> two is always the one that is cut, and never K -- a K split would be a
    !> reduction, and this is not.
    !>
    !> Chunk count.  One per thread is the optimum wherever the chunks stay
    !> above ~25 rows (nb=256 on M=3248, i.e. 12 rows, is already 30% off), so
    !> it is min(threads, M/split_min_rows).
    !>
    !> DETERMINISM.  Splitting a FREE dimension does not reorder any sum: each
    !> C(i,j) is still one BLAS call's accumulation over the whole k.  That is
    !> not merely an argument -- it was checked, on this shape, on this
    !> library: every chunk count from 1 to 128, at every OMP_NUM_THREADS from
    !> 1 to 128, gives a BITWISE identical C.  The result is therefore
    !> independent of the thread count, which the BLAS-threaded form it
    !> replaces was NOT: OpenBLAS's own k-blocking varies with its thread
    !> count and moved this product by ~2e-13 absolute between 1 and 32
    !> threads.  This is a strict improvement in reproducibility, not a
    !> concession.
    subroutine gemm_tn_split(mm, nn, kk, a, lda, b, ldb, beta, c, ldc)
      integer, intent(in) :: mm, nn, kk, lda, ldb, ldc
      ! Assumed size, so the chunk offsets below are element references and
      ! the columns after them stay contiguous by sequence association.
      real(dp), intent(in) :: a(0:*), b(0:*)
      real(dp), intent(in) :: beta
      real(dp), intent(inout) :: c(0:*)
      integer :: ic, nchunk, i0, ilen

      if (int(mm, i8)*int(nn, i8)*int(kk, i8) < split_min_madds &
          .or. nsplit <= 1) then
        call dgemm('T', 'N', mm, nn, kk, 1.0_dp, a, lda, b, ldb, beta, c, ldc)
        return
      end if

      if (mm >= nn) then
        nchunk = max(1, min(nsplit, mm / split_min_rows))
        !$omp parallel do default(shared) private(ic, i0, ilen) &
        !$omp   schedule(static) if (nchunk > 1)
        do ic = 0, nchunk - 1
          i0   = int(int(mm, i8)*int(ic, i8) / int(nchunk, i8))
          ilen = int(int(mm, i8)*int(ic + 1, i8) / int(nchunk, i8)) - i0
          if (ilen > 0) &
            call dgemm('T', 'N', ilen, nn, kk, 1.0_dp, &
                       a(int(i0, i8)*int(lda, i8)), lda, b, ldb, &
                       beta, c(int(i0, i8)), ldc)
        end do
      else
        nchunk = max(1, min(nsplit, nn / split_min_rows))
        !$omp parallel do default(shared) private(ic, i0, ilen) &
        !$omp   schedule(static) if (nchunk > 1)
        do ic = 0, nchunk - 1
          i0   = int(int(nn, i8)*int(ic, i8) / int(nchunk, i8))
          ilen = int(int(nn, i8)*int(ic + 1, i8) / int(nchunk, i8)) - i0
          if (ilen > 0) &
            call dgemm('T', 'N', mm, ilen, kk, 1.0_dp, a, lda, &
                       b(int(i0, i8)*int(ldb, i8)), ldb, &
                       beta, c(int(i0, i8)*int(ldc, i8)), ldc)
        end do
      end if
    end subroutine gemm_tn_split

    !> Everything downstream of C1 for one pair: the slot-1 slabs, the one-index
    !> derivative of h, Z, the Hessian column and the folded active integrals.
    !> `cmat` holds C1 on entry, however it was formed.  The five nbf x nbf
    !> scratch buffers are arguments rather than host-associated because the
    !> caller runs this loop across threads and each thread needs its own.
    subroutine finish_pair(lp, pp, qq, cmat, tmat, w1, w2, zmat)
      integer, intent(in) :: lp, pp, qq
      ! Explicit shape, not assumed shape: these go straight on to DGEMM,
      ! which has no explicit interface, so an assumed-shape dummy would make
      ! gfortran pass a contiguous temporary copy instead of the buffer.
      real(dp), intent(inout) :: cmat(0:n-1, 0:n-1), tmat(0:n-1, 0:n-1)
      real(dp), intent(inout) :: w1(0:n-1, 0:n-1), w2(0:n-1, 0:n-1)
      real(dp), intent(inout) :: zmat(0:n-1, 0:n-1)
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
      call ptockp(8)

      ! ---- Z[m,n] = (D t - t D)[n,m] + 2 C1[m,n]
      call dgemm('N', 'N', n, n, n, 1.0_dp, dm, n, tmat, n, 0.0_dp, w1, n)
      call dgemm('N', 'N', n, n, n, 1.0_dp, tmat, n, dm, n, 0.0_dp, w2, n)
      call ptockp(9)
      do mm = 0, n - 1
        do jj = 0, n - 1
          zmat(mm, jj) = (w1(jj, mm) - w2(jj, mm)) + 2.0_dp * cmat(mm, jj)
        end do
      end do
      call ptockp(10)

      ! ---- B[k,l] = Z[p_k,q_k] - Z[q_k,p_k]
      do kk2 = 0, np - 1
        ii = int(pairs(2*kk2))
        jj = int(pairs(2*kk2 + 1))
        bmat(int(kk2, i8)*int(np, i8) + int(lp, i8)) = zmat(ii, jj) - zmat(jj, ii)
      end do
      call ptockp(11)

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
      call ptockp(12)
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
      call ptockp(13)
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
