!> @file cc_lib.F90
!>
!> @brief Closed-shell (RHF-reference) CCSD and perturbative triples, CCSD(T).
!>
!> Spin-adapted formulation -- one closed-shell amplitude set instead of the
!> four spin blocks of a spin-orbital code, which cuts both the operation count
!> and the amplitude storage by roughly an order of magnitude at identical
!> energies.
!>
!>   CCSD amplitude equations: Hirata, Podeszwa, Tobita, Bartlett,
!>       J. Chem. Phys. 120, 2581 (2004).
!>   (T) correction: Raghavachari, Trucks, Pople, Head-Gordon,
!>       Chem. Phys. Lett. 157, 479 (1989), in the closed-shell
!>       a >= b >= c virtual-triple form.
!>
!> Performance design (follows current practice for distributed CCSD(T), e.g.
!> Kallay and co-workers, J. Chem. Theory Comput. 16, 366 (2020); Datta and
!> Gordon, J. Chem. Theory Comput. 17, 4799 (2021)):
!>
!>   * Every O(N^6)/O(N^7) contraction is cast as a single large DGEMM.  The
!>     index orders of the working arrays below are chosen so that the pair of
!>     contracted indices is always the fastest-running pair of one operand and
!>     the slowest of the other -- no transposes inside the hot loops.
!>   * The particle-particle ladder (the O(o^2 v^4) CCSD bottleneck) is blocked
!>     over the last virtual index.  Blocking serves three purposes at once: the
!>     t1-dressed ladder integrals are built one block at a time, so the second
!>     v^4 array never exists, and the blocks are the unit of both MPI and
!>     OpenMP work.  The block size is capped so there are several blocks per
!>     rank and per thread -- the memory target alone would hand back a single
!>     block for moderate v and leave every worker but one idle.
!>   * The (T) correction is parallelised over a >= b >= c virtual triples.
!>     That is ~v^3/6 independent tasks -- far more than the o^3/6 of an
!>     occupied-triple decomposition -- so the dynamic schedule stays balanced
!>     out to large rank counts.  MPI takes a round-robin stripe of triples and
!>     OpenMP consumes each rank's stripe with a dynamic schedule.
!>   * BLAS is pinned to a single thread for the whole solver.  Every hot region
!>     is already threaded here, and letting a threaded BLAS open a second pool
!>     underneath makes the idle pool spin-wait through the other's region --
!>     measured ~2x slower on four cores, and worse than serial at the largest
!>     thread count tested.  One pool, restored to the caller's setting on exit.
!>
!> Validated against PySCF RCCSD/CCSD(T) to better than 1e-10 Hartree.
module cc_lib

  use precision, only: dp
  use parallel, only: par_env_t
  use blas_thread, only: blas_thread_count, blas_thread_set
  use memory_info, only: oqp_available_memory_gb, OQP_MEMORY_SAFETY_FRACTION
  use io_constants, only: iw
  use, intrinsic :: iso_c_binding, only: c_int64_t

  implicit none

  private
  public :: cc_ccsd_t_energy
  public :: cc_options_t

  !> Ladder blocks handed to each MPI rank.  More than one so a rank that
  !> stalls does not hold up the reduction; small enough that the per-block
  !> DGEMM stays large.
  integer, parameter :: BLOCKS_PER_RANK = 4

  !> Runtime controls for the coupled-cluster solver.
  type :: cc_options_t
    integer  :: maxit     = 50        !< max CCSD iterations
    real(dp) :: conv      = 1.0e-7_dp !< convergence on the amplitude RMS
    integer  :: ndiis     = 8         !< DIIS subspace size (0 disables DIIS)
    logical  :: do_triples = .true.   !< evaluate the (T) correction
    integer  :: verbose   = 1         !< 0 silent, 1 iteration table
    integer  :: iw        = 6         !< output unit
  end type cc_options_t

contains

!###############################################################################

!> @brief Driver: CCSD iterations followed by the optional (T) correction.
!>
!> All two-electron integrals are supplied in the MO basis in chemist notation
!> and in the index order documented for each dummy argument.  @p eo and @p ev
!> are the occupied and virtual orbital energies of the (canonical) RHF
!> reference; the Fock matrix is therefore diagonal and its off-diagonal
!> occupied-virtual block, which would otherwise enter the singles equation and
!> the disconnected triples, vanishes.
subroutine cc_ccsd_t_energy(no, nv, eo, ev, oooo, ooov, oovv, ovov, ovvv, vvvv, &
                            pe, opts, e_ccsd, e_t, converged, &
                            time_ccsd, time_triples, bvv, nchol)

  integer, intent(in) :: no            !< correlated occupied orbitals
  integer, intent(in) :: nv            !< virtual orbitals
  real(dp), intent(in) :: eo(no), ev(nv)
  real(dp), intent(in) :: oooo(no,no,no,no)  !< (ij|kl)
  real(dp), intent(in) :: ooov(no,no,no,nv)  !< (ij|ka)
  real(dp), intent(in) :: oovv(no,no,nv,nv)  !< (ij|ab)
  real(dp), intent(in) :: ovov(no,nv,no,nv)  !< (ia|jb)
  real(dp), intent(in) :: ovvv(no,nv,nv,nv)  !< (ia|bc)
  real(dp), intent(in) :: vvvv(nv,nv,nv,nv)  !< (ab|cd)
  type(par_env_t), intent(inout) :: pe
  type(cc_options_t), intent(in) :: opts
  real(dp), intent(out) :: e_ccsd
  real(dp), intent(out) :: e_t
  logical, intent(out) :: converged
  !> Wall-clock seconds spent in the CCSD iterations and in the (T) correction.
  real(dp), optional, intent(out) :: time_ccsd, time_triples
  !> Cholesky vectors over the virtual pair block.  Supplying them makes the
  !> ladder rebuild its integrals instead of reading @p vvvv, so the caller may
  !> pass a zero-sized vvvv and never pay for the v^4 array.
  real(dp), optional, intent(in) :: bvv(:,:)
  integer, optional, intent(in) :: nchol

  real(dp), allocatable :: t1(:,:), t2(:,:,:,:)
  real(dp) :: t0, t1w
  integer(c_int64_t) :: nblas_save

  e_ccsd = 0.0_dp
  e_t = 0.0_dp
  converged = .false.
  if (present(time_ccsd)) time_ccsd = 0.0_dp
  if (present(time_triples)) time_triples = 0.0_dp

  if (no <= 0 .or. nv <= 0) return

  allocate(t1(no,nv), t2(no,no,nv,nv))

  ! Both the CCSD iterations and the (T) correction drive their own OpenMP
  ! regions and call BLAS from inside them.  Leaving BLAS threaded would run two
  ! pools against each other -- whichever is idle spin-waits through the other's
  ! region and burns the cores (measured ~2x slowdown on four cores).  Pin BLAS
  ! to one thread for the whole solver and restore the caller's setting after.
  nblas_save = blas_thread_count()
  if (nblas_save > 0) call blas_thread_set(1_c_int64_t)

  call cc_wall_time(t0)
  call ccsd_iterate(no, nv, eo, ev, oooo, ooov, oovv, ovov, ovvv, vvvv, &
                    pe, opts, t1, t2, e_ccsd, converged, bvv, nchol)
  call cc_wall_time(t1w)
  if (present(time_ccsd)) time_ccsd = t1w - t0

  ! Only if the amplitudes are actually converged.  The driver aborts on a
  ! failed iteration anyway, and (T) is the dominant O(N^7) cost, so running
  ! it first would spend the largest part of the job producing a number that
  ! is then thrown away.
  if (opts%do_triples .and. converged) then
    call cc_wall_time(t0)
    call triples_correction(no, nv, eo, ev, ooov, ovov, ovvv, t1, t2, pe, opts, e_t)
    call cc_wall_time(t1w)
    if (present(time_triples)) time_triples = t1w - t0
  end if

  if (nblas_save > 0) call blas_thread_set(nblas_save)

  deallocate(t1, t2)

end subroutine cc_ccsd_t_energy

!###############################################################################

!> @brief OpenMP threads available to this region (1 without OpenMP).
integer function nthreads_now() result(n)
!$ use omp_lib, only: omp_get_max_threads
  n = 1
!$ n = omp_get_max_threads()
end function nthreads_now

!###############################################################################

!> @brief Wall-clock seconds from the system clock.
subroutine cc_wall_time(t)
  real(dp), intent(out) :: t
  integer(8) :: c, r
  call system_clock(c, r)
  t = real(c, dp) / real(r, dp)
end subroutine cc_wall_time

!###############################################################################

!> @brief Spin-adapted closed-shell CCSD amplitude iterations with DIIS.
subroutine ccsd_iterate(no, nv, eo, ev, oooo, ooov, oovv, ovov, ovvv, vvvv, &
                        pe, opts, t1, t2, ecc, converged, bvv, nchol)

  integer, intent(in) :: no, nv
  real(dp), intent(in) :: eo(no), ev(nv)
  real(dp), intent(in) :: oooo(no,no,no,no), ooov(no,no,no,nv), oovv(no,no,nv,nv)
  real(dp), intent(in) :: ovov(no,nv,no,nv), ovvv(no,nv,nv,nv), vvvv(nv,nv,nv,nv)
  type(par_env_t), intent(inout) :: pe
  type(cc_options_t), intent(in) :: opts
  real(dp), intent(out) :: t1(no,nv), t2(no,no,nv,nv)
  real(dp), intent(out) :: ecc
  logical, intent(out) :: converged
  !> Cholesky vectors over the virtual pair block; forwarded to the ladder,
  !> which rebuilds its integrals from them instead of reading vvvv.
  real(dp), intent(in), optional :: bvv(:,:)
  integer, intent(in), optional :: nchol

  ! --- persistent working arrays -------------------------------------------
  real(dp), allocatable :: tau(:,:,:,:), t1n(:,:), t2n(:,:,:,:), half(:,:,:,:)
  real(dp), allocatable :: Lovov(:,:,:,:)   ! L(k,c,l,d) = 2(kc|ld) - (kd|lc)
  real(dp), allocatable :: ovov_s(:,:,:,:)  ! S(l,d,k,c) = (lc|kd)
  real(dp), allocatable :: Wklcd(:,:,:,:)   ! (k,l,c,d) = (kc|ld)
  real(dp), allocatable :: Lovvv(:,:,:,:)   ! 2(kd|ac) - (kc|ad), order (k,d,a,c)
  real(dp), allocatable :: Looov(:,:,:,:)   ! 2(ki|lc) - (li|kc), order (k,i,l,c)
  real(dp), allocatable :: Foo(:,:), Fvv(:,:), Fov(:,:), Loo(:,:), Lvv(:,:)
  real(dp), allocatable :: Woooo(:,:,:,:)
  real(dp), allocatable :: WA(:,:,:,:), WB(:,:,:,:)   ! (i,a,k,c) ring intermediates
  real(dp), allocatable :: M1(:,:,:,:), M2(:,:,:,:)
  real(dp), allocatable :: t2q(:,:,:,:), t2q2(:,:,:,:)
  real(dp), allocatable :: rbuf(:,:,:,:), rbuf2(:,:,:,:)
  real(dp), allocatable :: tmp2a(:,:,:,:), tmp2b(:,:,:,:)
  real(dp), allocatable :: taux(:,:,:,:), tauy(:,:,:,:)
  real(dp), allocatable :: dia(:,:), dijab(:,:,:,:)
  ! DIIS
  real(dp), allocatable :: dt1(:,:,:), dt2(:,:,:,:,:), et1(:,:,:), et2(:,:,:,:,:)

  integer :: i, j, k, l, a, b, c, d, it, nvec
  integer :: no2, nv2, nov
  real(dp) :: rms, eold

  no2 = no*no
  nv2 = nv*nv
  nov = no*nv

  allocate(tau(no,no,nv,nv), t1n(no,nv), t2n(no,no,nv,nv), half(no,no,nv,nv))
  allocate(Lovov(no,nv,no,nv), ovov_s(no,nv,no,nv), Wklcd(no,no,nv,nv))
  allocate(Lovvv(no,nv,nv,nv), Looov(no,no,no,nv))
  allocate(Foo(no,no), Fvv(nv,nv), Fov(no,nv), Loo(no,no), Lvv(nv,nv))
  allocate(Woooo(no,no,no,no))
  allocate(WA(no,nv,no,nv), WB(no,nv,no,nv), M1(no,nv,no,nv), M2(no,nv,no,nv))
  allocate(t2q(no,nv,no,nv), t2q2(no,nv,no,nv))
  allocate(rbuf(no,nv,no,nv), rbuf2(no,nv,no,nv))
  allocate(tmp2a(nv,nv,no,nv), tmp2b(nv,no,no,no))
  allocate(taux(nv,no,nv,no), tauy(no,no,nv,nv))
  allocate(dia(no,nv), dijab(no,no,nv,nv))

  ! --- reusable integral permutations --------------------------------------
  !$omp parallel do collapse(3) private(i,j,k,l,a,b,c,d) schedule(static)
  do d = 1, nv
    do l = 1, no
      do c = 1, nv
        do k = 1, no
          Lovov(k,c,l,d) = 2.0_dp*ovov(k,c,l,d) - ovov(k,d,l,c)
          ovov_s(k,c,l,d) = ovov(k,d,l,c)
          Wklcd(k,l,c,d) = ovov(k,c,l,d)
        end do
      end do
    end do
  end do
  !$omp end parallel do

  !$omp parallel do collapse(3) private(k,d,a,c) schedule(static)
  do c = 1, nv
    do a = 1, nv
      do d = 1, nv
        do k = 1, no
          Lovvv(k,d,a,c) = 2.0_dp*ovvv(k,d,a,c) - ovvv(k,c,a,d)
        end do
      end do
    end do
  end do
  !$omp end parallel do

  !$omp parallel do collapse(3) private(k,i,l,c) schedule(static)
  do c = 1, nv
    do l = 1, no
      do i = 1, no
        do k = 1, no
          Looov(k,i,l,c) = 2.0_dp*ooov(k,i,l,c) - ooov(l,i,k,c)
        end do
      end do
    end do
  end do
  !$omp end parallel do

  ! --- denominators ---------------------------------------------------------
  do a = 1, nv
    do i = 1, no
      dia(i,a) = eo(i) - ev(a)
    end do
  end do
  !$omp parallel do collapse(3) private(i,j,a,b) schedule(static)
  do b = 1, nv
    do a = 1, nv
      do j = 1, no
        do i = 1, no
          dijab(i,j,a,b) = eo(i) + eo(j) - ev(a) - ev(b)
        end do
      end do
    end do
  end do
  !$omp end parallel do

  ! --- MP2 start ------------------------------------------------------------
  t1 = 0.0_dp
  !$omp parallel do collapse(3) private(i,j,a,b) schedule(static)
  do b = 1, nv
    do a = 1, nv
      do j = 1, no
        do i = 1, no
          t2(i,j,a,b) = ovov(i,a,j,b) / dijab(i,j,a,b)
        end do
      end do
    end do
  end do
  !$omp end parallel do

  if (opts%ndiis > 0) then
    allocate(dt1(no,nv,opts%ndiis), dt2(no,no,nv,nv,opts%ndiis))
    allocate(et1(no,nv,opts%ndiis), et2(no,no,nv,nv,opts%ndiis))
  end if
  nvec = 0

  ecc = cc_energy(no, nv, ovov, t1, t2)
  if (opts%verbose > 0 .and. pe%rank == 0) then
    write(opts%iw,'(/4x,"CCSD iterations")')
    write(opts%iw,'(4x,"iter",9x,"E(corr)",11x,"dE",12x,"RMS")')
    write(opts%iw,'(4x,i4,f18.10,a)') 0, ecc, '        (MP2 guess)'
  end if

  converged = .false.

  do it = 1, opts%maxit
    eold = ecc

    ! tau_ij^ab = t2 + t1 t1
    !$omp parallel do collapse(3) private(i,j,a,b) schedule(static)
    do b = 1, nv
      do a = 1, nv
        do j = 1, no
          do i = 1, no
            tau(i,j,a,b) = t2(i,j,a,b) + t1(i,a)*t1(j,b)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! ================= one-particle intermediates =========================
    ! Foo(k,i) = sum_lcd L(k,c,l,d) tau(i,l,c,d)
    !$omp parallel do collapse(3) private(i,l,c,d) schedule(static)
    do i = 1, no
      do d = 1, nv
        do l = 1, no
          do c = 1, nv
            taux(c,l,d,i) = tau(i,l,c,d)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    call dgemm('n','n', no, no, nv*no*nv, 1.0_dp, Lovov, no, taux, nv*no*nv, &
               0.0_dp, Foo, no)

    ! Fvv(a,c) = - sum_kld L(k,c,l,d) tau(k,l,a,d)
    !   contract (k,l,d) : reuse tauy(k,l,d,a) and Lovov re-viewed as (c,(k,l,d))
    !$omp parallel do collapse(3) private(k,l,d,a) schedule(static)
    do a = 1, nv
      do d = 1, nv
        do l = 1, no
          do k = 1, no
            tauy(k,l,d,a) = tau(k,l,a,d)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    block
      real(dp), allocatable :: Lc(:,:), Fvvt(:,:)
      allocate(Lc(nv, no*no*nv), Fvvt(nv,nv))
      !$omp parallel do collapse(3) private(k,c,l,d) schedule(static)
      do d = 1, nv
        do l = 1, no
          do c = 1, nv
            do k = 1, no
              Lc(c, k + (l-1)*no + (d-1)*no*no) = Lovov(k,c,l,d)
            end do
          end do
        end do
      end do
      !$omp end parallel do
      ! DGEMM yields the (c,a) ordering; Fvv is stored as (a,c).
      call dgemm('n','n', nv, nv, no*no*nv, -1.0_dp, Lc, nv, tauy, no*no*nv, &
                 0.0_dp, Fvvt, nv)
      Fvv = transpose(Fvvt)
      deallocate(Lc, Fvvt)
    end block

    ! Fov(k,c) = sum_ld L(k,c,l,d) t1(l,d)
    call dgemv('n', nov, nov, 1.0_dp, Lovov, nov, t1, 1, 0.0_dp, Fov, 1)

    ! Loo = Foo + sum_lc [2(ki|lc) - (li|kc)] t1(l,c)
    Loo = Foo
    call dgemv('n', no2, nov, 1.0_dp, Looov, no2, t1, 1, 1.0_dp, Loo, 1)

    ! Lvv = Fvv + sum_kd [2(kd|ac) - (kc|ad)] t1(k,d)
    Lvv = Fvv
    call dgemv('t', nov, nv2, 1.0_dp, Lovvv, nov, t1, 1, 1.0_dp, Lvv, 1)

    ! ================= Woooo ==============================================
    ! W(k,l,i,j) = (ki|lj) + sum_c (ki|lc) t1(j,c) + sum_c (lj|kc) t1(i,c)
    !            + sum_cd (kc|ld) tau(i,j,c,d)
    call dgemm('n','t', no2, no2, nv2, 1.0_dp, Wklcd, no2, tau, no2, 0.0_dp, Woooo, no2)
    !$omp parallel do collapse(3) private(i,j,k,l,c) schedule(static)
    do j = 1, no
      do i = 1, no
        do l = 1, no
          do k = 1, no
            Woooo(k,l,i,j) = Woooo(k,l,i,j) + oooo(k,i,l,j)
            do c = 1, nv
              Woooo(k,l,i,j) = Woooo(k,l,i,j) + ooov(k,i,l,c)*t1(j,c) &
                                              + ooov(l,j,k,c)*t1(i,c)
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! ================= ring intermediates WA, WB ==========================
    ! WA(i,a,k,c) == Wvoov[a,k,i,c],  WB(i,a,k,c) == Wvovo[a,k,c,i]
    !$omp parallel do collapse(3) private(i,a,l,d) schedule(static)
    do d = 1, nv
      do l = 1, no
        do a = 1, nv
          do i = 1, no
            M1(i,a,l,d) = t2(i,l,a,d) - 0.5_dp*t2(i,l,d,a) - t1(i,d)*t1(l,a)
            M2(i,a,l,d) = -0.5_dp*t2(i,l,a,d)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    ! WA(i,a,k,c) = sum_ld M1(i,a,l,d) (ld|kc)  +  sum_ld M2(i,a,l,d) (lc|kd)
    call dgemm('n','n', nov, nov, nov, 1.0_dp, M1, nov, ovov, nov, 0.0_dp, WA, nov)
    call dgemm('n','n', nov, nov, nov, 1.0_dp, M2, nov, ovov_s, nov, 1.0_dp, WA, nov)

    !$omp parallel do collapse(3) private(i,a,l,d) schedule(static)
    do d = 1, nv
      do l = 1, no
        do a = 1, nv
          do i = 1, no
            M1(i,a,l,d) = -(0.5_dp*t2(i,l,d,a) + t1(i,d)*t1(l,a))
          end do
        end do
      end do
    end do
    !$omp end parallel do
    call dgemm('n','n', nov, nov, nov, 1.0_dp, M1, nov, ovov_s, nov, 0.0_dp, WB, nov)

    !$omp parallel do collapse(3) private(i,a,k,c,d,l) schedule(static)
    do c = 1, nv
      do k = 1, no
        do a = 1, nv
          do i = 1, no
            WA(i,a,k,c) = WA(i,a,k,c) + ovov(k,c,i,a)
            WB(i,a,k,c) = WB(i,a,k,c) + oovv(k,i,a,c)
            do d = 1, nv
              WA(i,a,k,c) = WA(i,a,k,c) + ovvv(k,c,a,d)*t1(i,d)
              WB(i,a,k,c) = WB(i,a,k,c) + ovvv(k,d,a,c)*t1(i,d)
            end do
            do l = 1, no
              WA(i,a,k,c) = WA(i,a,k,c) - ooov(l,i,k,c)*t1(l,a)
              WB(i,a,k,c) = WB(i,a,k,c) - ooov(k,i,l,c)*t1(l,a)
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! ================= T1 =================================================
    call ccsd_t1_equation(no, nv, ooov, oovv, ovov, ovvv, t1, t2, &
                          Foo, Fvv, Fov, t1n)

    ! ================= T2 =================================================
    ! Terms split into two groups.  `t2n` collects the contributions that are
    ! already invariant under (ij,ab) -> (ji,ba); `half` collects those that the
    ! working equations pair with their own transpose.  One symmetrisation at
    ! the end then covers the whole `half` group at once.
    t2n = 0.0_dp
    half = 0.0_dp

    ! (a) t1-dressed <ab|ic> and <ak|ij> driver terms
    !$omp parallel do collapse(3) private(a,b,i,c,k) schedule(static)
    do c = 1, nv
      do i = 1, no
        do b = 1, nv
          do a = 1, nv
            tmp2a(a,b,i,c) = ovvv(i,a,b,c)
            do k = 1, no
              tmp2a(a,b,i,c) = tmp2a(a,b,i,c) - oovv(k,i,b,c)*t1(k,a)
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do collapse(3) private(a,k,i,j,c) schedule(static)
    do j = 1, no
      do i = 1, no
        do k = 1, no
          do a = 1, nv
            tmp2b(a,k,i,j) = ooov(j,k,i,a)
            do c = 1, nv
              tmp2b(a,k,i,j) = tmp2b(a,k,i,j) + ovov(k,c,i,a)*t1(j,c)
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(3) private(i,j,a,b,c,k) schedule(static)
    do b = 1, nv
      do a = 1, nv
        do j = 1, no
          do i = 1, no
            do c = 1, nv
              half(i,j,a,b) = half(i,j,a,b) + tmp2a(a,b,i,c)*t1(j,c)
            end do
            do k = 1, no
              half(i,j,a,b) = half(i,j,a,b) - tmp2b(a,k,i,j)*t1(k,b)
            end do
            ! (ia|jb) is invariant under the pair exchange -- added once.
            t2n(i,j,a,b) = t2n(i,j,a,b) + ovov(i,a,j,b)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! (b) Woooo ladder :  sum_kl W(k,l,i,j) tau(k,l,a,b)
    call dgemm('t','n', no2, nv2, no2, 1.0_dp, Woooo, no2, tau, no2, 1.0_dp, t2n, no2)

    ! (c) particle-particle ladder, blocked over d and distributed over MPI
    call ladder_contraction(no, nv, vvvv, ovvv, t1, tau, pe, t2n, bvv, nchol)

    ! (d) Lvv / Loo dressing.  The Lvv contraction runs over the third index of
    ! t2, so it is issued as one DGEMM per b rather than a single flat call.
    do b = 1, nv
      call dgemm('n','t', no2, nv, nv, 1.0_dp, t2(1,1,1,b), no2, Lvv, nv, &
                 1.0_dp, half(1,1,1,b), no2)
    end do
    call dgemm('t','n', no, no*nv2, no, -1.0_dp, Loo, no, t2, no, 1.0_dp, half, no)

    ! (e) ring terms
    !$omp parallel do collapse(3) private(k,c,j,b) schedule(static)
    do b = 1, nv
      do j = 1, no
        do c = 1, nv
          do k = 1, no
            t2q(k,c,j,b)  = t2(k,j,c,b)
            t2q2(k,c,j,b) = t2(k,j,b,c)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(3) private(i,a,k,c) schedule(static)
    do c = 1, nv
      do k = 1, no
        do a = 1, nv
          do i = 1, no
            M1(i,a,k,c) = 2.0_dp*WA(i,a,k,c) - WB(i,a,k,c)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    call dgemm('n','n', nov, nov, nov, 1.0_dp, M1, nov, t2q, nov, 0.0_dp, rbuf, nov)
    call dgemm('n','n', nov, nov, nov, -1.0_dp, WA, nov, t2q2, nov, 1.0_dp, rbuf, nov)
    call dgemm('n','n', nov, nov, nov, -1.0_dp, WB, nov, t2q2, nov, 0.0_dp, rbuf2, nov)

    !$omp parallel do collapse(3) private(i,j,a,b) schedule(static)
    do b = 1, nv
      do a = 1, nv
        do j = 1, no
          do i = 1, no
            half(i,j,a,b) = half(i,j,a,b) + rbuf(i,a,j,b) + rbuf2(i,b,j,a)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! (f) fold in the transpose partner of the `half` group, then divide
    !$omp parallel do collapse(3) private(i,j,a,b) schedule(static)
    do b = 1, nv
      do a = 1, nv
        do j = 1, no
          do i = 1, no
            t2n(i,j,a,b) = (t2n(i,j,a,b) + half(i,j,a,b) + half(j,i,b,a)) &
                           / dijab(i,j,a,b)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do collapse(1) private(i,a) schedule(static)
    do a = 1, nv
      do i = 1, no
        t1n(i,a) = t1n(i,a) / dia(i,a)
      end do
    end do
    !$omp end parallel do

    ! --- convergence + DIIS ------------------------------------------------
    rms = sqrt((sum((t1n-t1)**2) + sum((t2n-t2)**2)) / real(no*nv + no2*nv2, dp))

    if (opts%ndiis > 0) then
      call diis_extrapolate(no, nv, opts%ndiis, nvec, t1, t2, t1n, t2n, &
                            dt1, dt2, et1, et2)
    else
      t1 = t1n
      t2 = t2n
    end if

    ecc = cc_energy(no, nv, ovov, t1, t2)

    if (opts%verbose > 0 .and. pe%rank == 0) then
      write(opts%iw,'(4x,i4,f18.10,2es16.5)') it, ecc, ecc-eold, rms
      flush(opts%iw)
    end if

    if (rms < opts%conv .and. abs(ecc-eold) < opts%conv) then
      converged = .true.
      exit
    end if
  end do

  if (opts%verbose > 0 .and. pe%rank == 0) then
    if (converged) then
      write(opts%iw,'(4x,"CCSD converged.")')
    else
      write(opts%iw,'(4x,"*** WARNING: CCSD did NOT converge in ",i0," iterations")') opts%maxit
    end if
  end if

  deallocate(tau, t1n, t2n, half, Lovov, ovov_s, Wklcd, Lovvv, Looov)
  deallocate(Foo, Fvv, Fov, Loo, Lvv, Woooo, WA, WB, M1, M2)
  deallocate(t2q, t2q2, rbuf, rbuf2, tmp2a, tmp2b, taux, tauy, dia, dijab)
  if (allocated(dt1)) deallocate(dt1, dt2, et1, et2)

end subroutine ccsd_iterate

!###############################################################################

!> @brief Singles residual (before division by the denominator).
subroutine ccsd_t1_equation(no, nv, ooov, oovv, ovov, ovvv, t1, t2, &
                            Foo, Fvv, Fov, t1n)
  integer, intent(in) :: no, nv
  real(dp), intent(in) :: ooov(no,no,no,nv), oovv(no,no,nv,nv)
  real(dp), intent(in) :: ovov(no,nv,no,nv), ovvv(no,nv,nv,nv)
  real(dp), intent(in) :: t1(no,nv), t2(no,no,nv,nv)
  real(dp), intent(in) :: Foo(no,no), Fvv(nv,nv), Fov(no,nv)
  real(dp), intent(out) :: t1n(no,nv)

  integer :: i, a, k, c, d, l
  real(dp) :: s

  !$omp parallel do collapse(2) private(i,a,k,c,d,l,s) schedule(static)
  do a = 1, nv
    do i = 1, no
      s = 0.0_dp
      do c = 1, nv
        s = s + Fvv(a,c)*t1(i,c)
      end do
      do k = 1, no
        s = s - Foo(k,i)*t1(k,a)
      end do
      do c = 1, nv
        do k = 1, no
          s = s + Fov(k,c)*(2.0_dp*t2(k,i,c,a) - t2(i,k,c,a)) &
                + Fov(k,c)*t1(i,c)*t1(k,a)                    &
                + 2.0_dp*ovov(k,c,i,a)*t1(k,c) - oovv(k,i,a,c)*t1(k,c)
        end do
      end do
      do d = 1, nv
        do c = 1, nv
          do k = 1, no
            s = s + (2.0_dp*ovvv(k,d,a,c) - ovvv(k,c,a,d)) &
                    * (t2(i,k,c,d) + t1(k,d)*t1(i,c))
          end do
        end do
      end do
      do c = 1, nv
        do l = 1, no
          do k = 1, no
            s = s - (2.0_dp*ooov(k,i,l,c) - ooov(l,i,k,c)) &
                    * (t2(k,l,a,c) + t1(l,c)*t1(k,a))
          end do
        end do
      end do
      t1n(i,a) = s
    end do
  end do
  !$omp end parallel do

end subroutine ccsd_t1_equation

!###############################################################################

!> @brief Particle-particle ladder  t2new(ij,ab) += sum_cd Wvvvv(a,b,c,d) tau(ij,cd).
!>
!> The t1-dressed ladder integral
!>     W(a,b,c,d) = (ac|bd) - sum_k (kd|ac) t1(k,b) - sum_k (kc|bd) t1(k,a)
!> is rebuilt one d-block at a time, so the second v^4 array is never
!> materialised.  Blocks are handed out round-robin to MPI ranks and the partial
!> results are reduced once at the end -- one collective per iteration.
subroutine ladder_contraction(no, nv, vvvv, ovvv, t1, tau, pe, t2n, bvv, nchol)

  integer, intent(in) :: no, nv
  real(dp), intent(in) :: vvvv(nv,nv,nv,nv), ovvv(no,nv,nv,nv)
  real(dp), intent(in) :: t1(no,nv), tau(no,no,nv,nv)
  type(par_env_t), intent(inout) :: pe
  real(dp), intent(inout) :: t2n(no,no,nv,nv)
  !> Cholesky vectors over the virtual-virtual block, (nv*nv, nchol), with
  !> (ac|bd) = sum_J bvv(ac,J) bvv(bd,J).  When present the ladder integrals
  !> are rebuilt from these per block and @p vvvv is never referenced, which is
  !> what lets the caller skip allocating the v^4 array at all.
  real(dp), intent(in), optional :: bvv(:,:)
  integer, intent(in), optional :: nchol

  real(dp), allocatable :: acc(:,:,:,:)
  integer :: nblk, bsize
  integer :: no2, nv2
  logical :: use_chol
  integer :: nchol_
  integer :: nthr_use, max_thr
  real(dp) :: target_dp, avail_gb, share_dp, per_thread

  no2 = no*no
  nv2 = nv*nv

  use_chol = present(bvv) .and. present(nchol)
  nchol_ = 0
  if (use_chol) nchol_ = nchol

  ! Target ~64 MB per dressed-integral block, at least one and at most all
  ! virtuals.  OQP_CC_LADDER_BLOCK overrides it (tuning and regression tests).
  !
  ! That target is per thread, so cap it by what this machine can give as
  ! well.  The cap only ever shrinks the tuned size: with room the behaviour
  ! is unchanged, and without it the job blocks smaller instead of dying in
  ! the allocator.
  target_dp = 8.0e6_dp
  avail_gb = oqp_available_memory_gb()
  if (avail_gb > 0.0_dp) then
    share_dp = OQP_MEMORY_SAFETY_FRACTION*avail_gb*1.073741824e9_dp/8.0_dp/3.0_dp
    share_dp = share_dp/real(max(1,nthreads_now()), dp)
    if (share_dp < target_dp) target_dp = max(real(nv2,dp)*real(nv,dp), share_dp)
  end if
  bsize = max(1, min(nv, int(target_dp / real(max(1,nv2*nv),dp))))

  ! The memory target alone says nothing about work granularity: for moderate
  ! nv it yields a SINGLE block, which would leave every rank but one idle.
  ! Under MPI, cap the block so there are several blocks per rank -- enough
  ! chunks to spread and to absorb jitter.  Serial runs keep the large blocks,
  ! where fewer, bigger DGEMMs are the faster choice.
  if (pe%size > 1) then
    bsize = max(1, min(bsize, nv / max(1, int(pe%size)*BLOCKS_PER_RANK)))
  end if

  block
    character(len=32) :: sval
    integer :: ln, bs_env
    call get_environment_variable("OQP_CC_LADDER_BLOCK", sval, ln)
    if (ln > 0) then
      read(sval,*,iostat=ln) bs_env
      if (ln == 0 .and. bs_env > 0) bsize = min(nv, bs_env)
    end if
  end block
  ! Enough blocks that the OpenMP schedule below has something to balance, even
  ! when memory alone would have asked for one big block.
  if (nthreads_now() > 1) then
    bsize = max(1, min(bsize, nv / max(1, nthreads_now()*BLOCKS_PER_RANK)))
  end if
  nblk = (nv + bsize - 1) / bsize

  allocate(acc(no,no,nv,nv), source=0.0_dp)

  ! How many threads this region can afford.  Wblk is genuine private working
  ! storage -- each thread builds a different block at the same time, so there
  ! is nothing to share -- and it shrinks with bsize, so the width is capped
  ! only when even the smallest useful block will not fit across the pool.
  nthr_use = max(1, nthreads_now())
  if (avail_gb > 0.0_dp) then
    per_thread = real(nv2,dp)*real(nv,dp)*real(bsize,dp)
    max_thr = int(OQP_MEMORY_SAFETY_FRACTION*avail_gb*1.073741824e9_dp/8.0_dp &
                  /3.0_dp/max(1.0_dp, per_thread))
    if (max_thr < nthr_use) then
      nthr_use = max(1, max_thr)
      write(iw,'(2X,A,I0,A,I0,A)') 'CCSD: ladder limited to ', nthr_use, &
          ' of ', nthreads_now(), ' threads by available memory'
    end if
  end if

  ! Parallelise over blocks rather than inside them.  Threading the block build
  ! and then calling a threaded DGEMM alternates two thread pools, and the idle
  ! one spin-waits through the other's region -- measured ~2x slowdown at four
  ! cores.  One pool owning whole blocks, with BLAS pinned to a single thread
  ! (restored by the caller), keeps every core on useful work.
  !
  ! The blocked index is b, an OUTPUT index, not the contracted d.  That is
  ! what keeps the per-thread cost proportional to the block: each thread
  ! writes a disjoint slice acc(:,:,:,b0:b1), so there is no private
  ! accumulator and no reduction at the end.  Blocking d instead made every
  ! thread accumulate the whole no^2 nv^2 output -- 500 MB per thread at
  ! nbf = 300, unshrinkable, and the largest thing in a wide run.
  !$omp parallel default(shared) num_threads(nthr_use)
  block
    ! Declared inside the construct, so each thread gets its own instance --
    ! an allocatable in a private() clause would arrive unallocated instead.
    real(dp), allocatable :: Wblk(:,:,:,:), Vblk(:,:), Bsub(:,:)
    integer :: a, b, c, d, bb, k, b0, b1, nbb, iblk, jj
    real(dp) :: s
    ! W((a,bb),(c,d)) for the b-block in hand.
    allocate(Wblk(nv,bsize,nv,nv))
    ! Cholesky path: V(ac,(bb,d)) for this block and the bvv rows it is built
    ! from.  Unlike a d-block, the (b,d) rows of a b-block are strided -- b
    ! runs inside each d -- so they are gathered rather than passed directly.
    if (use_chol) allocate(Vblk(nv2, bsize*nv), Bsub(bsize*nv, nchol_))

    !$omp do schedule(dynamic,1)
    do iblk = 1, nblk
      if (pe%size > 1) then
        if (mod(iblk-1, int(pe%size)) /= int(pe%rank)) cycle
      end if
      b0 = (iblk-1)*bsize + 1
      b1 = min(nv, b0 + bsize - 1)
      nbb = b1 - b0 + 1

      if (use_chol) then
        do jj = 1, nchol_
          do d = 1, nv
            do bb = 1, nbb
              Bsub(bb + (d-1)*nbb, jj) = bvv((b0+bb-1) + (d-1)*nv, jj)
            end do
          end do
        end do
        call dgemm('n','t', nv2, nbb*nv, nchol_, 1.0_dp, &
                   bvv, nv2, Bsub, bsize*nv, 0.0_dp, Vblk, nv2)
      end if

      do d = 1, nv
        do c = 1, nv
          do bb = 1, nbb
            b = b0 + bb - 1
            do a = 1, nv
              if (use_chol) then
                s = Vblk((c-1)*nv+a, bb + (d-1)*nbb)
              else
                s = vvvv(a,c,b,d)
              end if
              do k = 1, no
                s = s - ovvv(k,d,a,c)*t1(k,b) - ovvv(k,c,b,d)*t1(k,a)
              end do
              Wblk(a,bb,c,d) = s
            end do
          end do
        end do
      end do

      ! acc(ij,(a,b)) += sum_(c,d) tau(ij,(c,d)) * W((a,bb),(c,d)), written
      ! straight into this block's slice, which no other thread touches.
      !
      ! Wblk is dimensioned on bsize, so its (c,d) columns are nv*bsize apart
      ! whatever the current block holds.  The last block has nbb < bsize, and
      ! passing nv*nbb here walks the wrong stride -- silently, and differently
      ! for every block size, which is how it first showed up as an energy that
      ! depended on the thread count.
      call dgemm('n','t', no2, nv*nbb, nv2, 1.0_dp, tau, no2, &
                 Wblk, nv*bsize, 1.0_dp, acc(1,1,1,b0), no2)
    end do
    !$omp end do

    deallocate(Wblk)
    if (allocated(Vblk)) deallocate(Vblk, Bsub)
  end block
  !$omp end parallel

  if (pe%size > 1) call pe%allreduce(acc, size(acc))

  t2n = t2n + acc

  deallocate(acc)

end subroutine ladder_contraction

!###############################################################################

!> @brief Closed-shell CCSD correlation energy.
!>   E = sum_ijab [2(ia|jb) - (ib|ja)] (t2_ij^ab + t1_i^a t1_j^b)
pure function cc_energy(no, nv, ovov, t1, t2) result(e)
  integer, intent(in) :: no, nv
  real(dp), intent(in) :: ovov(no,nv,no,nv), t1(no,nv), t2(no,no,nv,nv)
  real(dp) :: e
  integer :: i, j, a, b
  e = 0.0_dp
  do b = 1, nv
    do a = 1, nv
      do j = 1, no
        do i = 1, no
          e = e + (2.0_dp*ovov(i,a,j,b) - ovov(i,b,j,a)) &
                  * (t2(i,j,a,b) + t1(i,a)*t1(j,b))
        end do
      end do
    end do
  end do
end function cc_energy

!###############################################################################

!> @brief DIIS extrapolation on the concatenated (t1,t2) amplitude vector.
subroutine diis_extrapolate(no, nv, ndiis, nvec, t1, t2, t1n, t2n, dt1, dt2, et1, et2)

  integer, intent(in) :: no, nv, ndiis
  integer, intent(inout) :: nvec
  real(dp), intent(inout) :: t1(no,nv), t2(no,no,nv,nv)
  real(dp), intent(in) :: t1n(no,nv), t2n(no,no,nv,nv)
  real(dp), intent(inout) :: dt1(no,nv,ndiis), dt2(no,no,nv,nv,ndiis)
  real(dp), intent(inout) :: et1(no,nv,ndiis), et2(no,no,nv,nv,ndiis)

  real(dp), allocatable :: bmat(:,:), rhs(:)
  integer, allocatable :: ipiv(:)
  integer :: i, j, n, info, islot

  ! shift the subspace when full (drop the oldest vector)
  if (nvec < ndiis) then
    nvec = nvec + 1
    islot = nvec
  else
    do i = 1, ndiis-1
      dt1(:,:,i) = dt1(:,:,i+1)
      dt2(:,:,:,:,i) = dt2(:,:,:,:,i+1)
      et1(:,:,i) = et1(:,:,i+1)
      et2(:,:,:,:,i) = et2(:,:,:,:,i+1)
    end do
    islot = ndiis
  end if

  dt1(:,:,islot) = t1n
  dt2(:,:,:,:,islot) = t2n
  et1(:,:,islot) = t1n - t1
  et2(:,:,:,:,islot) = t2n - t2

  n = nvec
  if (n < 2) then
    t1 = t1n
    t2 = t2n
    return
  end if

  allocate(bmat(n+1,n+1), rhs(n+1), ipiv(n+1))
  bmat = 0.0_dp
  do i = 1, n
    do j = 1, n
      bmat(i,j) = sum(et1(:,:,i)*et1(:,:,j)) + sum(et2(:,:,:,:,i)*et2(:,:,:,:,j))
    end do
  end do
  bmat(1:n, n+1) = -1.0_dp
  bmat(n+1, 1:n) = -1.0_dp
  bmat(n+1, n+1) = 0.0_dp
  rhs = 0.0_dp
  rhs(n+1) = -1.0_dp

  call dgesv(n+1, 1, bmat, n+1, ipiv, rhs, n+1, info)

  if (info /= 0) then
    t1 = t1n
    t2 = t2n
  else
    t1 = 0.0_dp
    t2 = 0.0_dp
    do i = 1, n
      t1 = t1 + rhs(i)*dt1(:,:,i)
      t2 = t2 + rhs(i)*dt2(:,:,:,:,i)
    end do
  end if

  deallocate(bmat, rhs, ipiv)

end subroutine diis_extrapolate

!###############################################################################

!> @brief Perturbative triples (T), closed-shell, over a >= b >= c.
!>
!> For each virtual triple the six index permutations of the connected W and
!> disconnected V arrays are formed with two DGEMMs each and contracted in the
!> 6x6 permutation pattern.  Work is striped round-robin over MPI ranks and the
!> stripe is consumed by a dynamic OpenMP schedule; BLAS runs single-threaded
!> inside the parallel region so the per-triple DGEMMs do not oversubscribe.
subroutine triples_correction(no, nv, eo, ev, ooov, ovov, ovvv, t1, t2, pe, opts, e_t)

!$ use omp_lib

  integer, intent(in) :: no, nv
  real(dp), intent(in) :: eo(no), ev(nv)
  real(dp), intent(in) :: ooov(no,no,no,nv), ovov(no,nv,no,nv), ovvv(no,nv,nv,nv)
  real(dp), intent(in) :: t1(no,nv), t2(no,no,nv,nv)
  type(par_env_t), intent(inout) :: pe
  type(cc_options_t), intent(in) :: opts
  real(dp), intent(out) :: e_t

  ! reordered drivers
  real(dp), allocatable :: Wvvov(:,:,:,:)  ! (i,f,a,b) = (ia|bf)
  real(dp), allocatable :: Wvooo(:,:,:,:)  ! (i,j,m,a) = (ia|jm)
  real(dp), allocatable :: t2c(:,:,:,:)    ! (f,j,k,c) = t2(k,j,c,f)
  real(dp), allocatable :: Wovov(:,:,:,:)  ! (i,j,a,b) = (ia|jb)
  real(dp), allocatable :: eijk(:,:,:)

  integer :: i, j, k, a, b, c, f, m, no3
  integer :: ntask, itask
  integer, allocatable :: ta(:), tb(:), tc(:)
  real(dp) :: et_local

  ! permutation bookkeeping: 1=abc 2=acb 3=bac 4=bca 5=cab 6=cba
  integer, parameter :: wid(6,6) = reshape( &
      [ 1,2,3,4,5,6, &
        2,1,5,6,3,4, &
        3,4,1,2,6,5, &
        4,3,6,5,1,2, &
        5,6,2,1,4,3, &
        6,5,4,3,2,1 ], [6,6] )

  e_t = 0.0_dp
  no3 = no*no*no

  allocate(Wvvov(no,nv,nv,nv), Wvooo(no,no,no,nv), t2c(nv,no,no,nv))
  allocate(Wovov(no,no,nv,nv), eijk(no,no,no))

  !$omp parallel do collapse(3) private(i,f,a,b) schedule(static)
  do b = 1, nv
    do a = 1, nv
      do f = 1, nv
        do i = 1, no
          Wvvov(i,f,a,b) = ovvv(i,a,b,f)
        end do
      end do
    end do
  end do
  !$omp end parallel do

  !$omp parallel do collapse(3) private(i,j,m,a) schedule(static)
  do a = 1, nv
    do m = 1, no
      do j = 1, no
        do i = 1, no
          Wvooo(i,j,m,a) = ooov(j,m,i,a)
        end do
      end do
    end do
  end do
  !$omp end parallel do

  !$omp parallel do collapse(3) private(f,j,k,c) schedule(static)
  do c = 1, nv
    do k = 1, no
      do j = 1, no
        do f = 1, nv
          t2c(f,j,k,c) = t2(k,j,c,f)
        end do
      end do
    end do
  end do
  !$omp end parallel do

  !$omp parallel do collapse(3) private(i,j,a,b) schedule(static)
  do b = 1, nv
    do a = 1, nv
      do j = 1, no
        do i = 1, no
          Wovov(i,j,a,b) = ovov(i,a,j,b)
        end do
      end do
    end do
  end do
  !$omp end parallel do

  do k = 1, no
    do j = 1, no
      do i = 1, no
        eijk(i,j,k) = eo(i) + eo(j) + eo(k)
      end do
    end do
  end do

  ! ---- flat task list over a >= b >= c ------------------------------------
  ntask = 0
  do a = 1, nv
    do b = 1, a
      do c = 1, b
        ntask = ntask + 1
      end do
    end do
  end do
  allocate(ta(ntask), tb(ntask), tc(ntask))
  itask = 0
  do a = 1, nv
    do b = 1, a
      do c = 1, b
        itask = itask + 1
        ta(itask) = a
        tb(itask) = b
        tc(itask) = c
      end do
    end do
  end do

  if (opts%verbose > 0 .and. pe%rank == 0) then
    write(opts%iw,'(/4x,"(T) correction over ",i0," virtual triples (a>=b>=c)")') ntask
    flush(opts%iw)
  end if

  et_local = 0.0_dp

  !$omp parallel default(shared) reduction(+:et_local)
  block
    real(dp), allocatable :: w(:,:,:,:), v(:,:,:,:), z(:,:,:), d3(:,:,:)
    integer :: p, q, t, ia, ib, ic, pa(6), pb(6), pc(6)
    real(dp) :: fac, acc

    allocate(w(no,no,no,6), v(no,no,no,6), z(no,no,no), d3(no,no,no))

    !$omp do schedule(dynamic,1)
    do t = 1, ntask
      if (pe%size > 1) then
        if (mod(t-1, int(pe%size)) /= int(pe%rank)) cycle
      end if
      ia = ta(t); ib = tb(t); ic = tc(t)

      pa = [ia, ia, ib, ib, ic, ic]
      pb = [ib, ic, ia, ic, ia, ib]
      pc = [ic, ib, ic, ia, ib, ia]

      fac = 1.0_dp
      if (ia == ic) then
        fac = 6.0_dp
      else if (ia == ib .or. ib == ic) then
        fac = 2.0_dp
      end if

      do k = 1, no
        do j = 1, no
          do i = 1, no
            d3(i,j,k) = (eijk(i,j,k) - ev(ia) - ev(ib) - ev(ic)) * fac
          end do
        end do
      end do

      ! W and V for the six permutations
      do p = 1, 6
        call dgemm('n','n', no, no*no, nv, 1.0_dp, Wvvov(1,1,pa(p),pb(p)), no, &
                   t2c(1,1,1,pc(p)), nv, 0.0_dp, w(1,1,1,p), no)
        call dgemm('n','n', no*no, no, no, -1.0_dp, Wvooo(1,1,1,pa(p)), no*no, &
                   t2(1,1,pb(p),pc(p)), no, 1.0_dp, w(1,1,1,p), no*no)
        do k = 1, no
          do j = 1, no
            do i = 1, no
              v(i,j,k,p) = Wovov(i,j,pa(p),pb(p)) * t1(k,pc(p))
            end do
          end do
        end do
      end do

      ! contract each Z against the six W permutations
      do p = 1, 6
        do k = 1, no
          do j = 1, no
            do i = 1, no
              z(i,j,k) = ( 4.0_dp*(w(i,j,k,p) + 0.5_dp*v(i,j,k,p))       &
                         + (w(j,k,i,p) + 0.5_dp*v(j,k,i,p))              &
                         + (w(k,i,j,p) + 0.5_dp*v(k,i,j,p))              &
                         - 2.0_dp*(w(k,j,i,p) + 0.5_dp*v(k,j,i,p))       &
                         - 2.0_dp*(w(i,k,j,p) + 0.5_dp*v(i,k,j,p))       &
                         - 2.0_dp*(w(j,i,k,p) + 0.5_dp*v(j,i,k,p)) )     &
                         / d3(i,j,k)
            end do
          end do
        end do

        acc = 0.0_dp
        do k = 1, no
          do j = 1, no
            do i = 1, no
              q = wid(1,p); acc = acc + w(i,j,k,q)*z(i,j,k)
              q = wid(2,p); acc = acc + w(i,k,j,q)*z(i,j,k)
              q = wid(3,p); acc = acc + w(j,i,k,q)*z(i,j,k)
              q = wid(4,p); acc = acc + w(j,k,i,q)*z(i,j,k)
              q = wid(5,p); acc = acc + w(k,i,j,q)*z(i,j,k)
              q = wid(6,p); acc = acc + w(k,j,i,q)*z(i,j,k)
            end do
          end do
        end do
        et_local = et_local + acc
      end do
    end do
    !$omp end do

    deallocate(w, v, z, d3)
  end block
  !$omp end parallel

  e_t = 2.0_dp * et_local
  if (pe%size > 1) call pe%allreduce(e_t, 1)

  deallocate(Wvvov, Wvooo, t2c, Wovov, eijk, ta, tb, tc)

end subroutine triples_correction

!###############################################################################

end module cc_lib
