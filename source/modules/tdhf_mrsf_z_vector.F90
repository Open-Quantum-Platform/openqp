module tdhf_mrsf_z_vector_mod

  use precision, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use zvector_common, only: sanitize_zvector_preconditioner
  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_z_vector_mod"
  real(kind=dp), parameter :: GMRES_DENOMINATOR_FLOOR = 1.0d-14
  real(kind=dp), parameter :: MRSF_ZVEC_DENOMINATOR_FLOOR = 1.0d-14

  ! Module-level work arrays for GMRES to avoid repeated allocation
  real(kind=8), allocatable :: gmres_wrk1(:,:), gmres_wrk2(:,:), gmres_wrk3(:,:)
  real(kind=8), allocatable, target :: gmres_pa(:,:,:)
  real(kind=8), allocatable :: gmres_ab1_mo_a(:,:), gmres_ab1_mo_b(:,:)
  logical :: gmres_work_allocated = .false.
  integer :: gmres_nbf = 0
  integer :: gmres_nocca = 0
  integer :: gmres_noccb = 0

  ! ----------------------------------------------------------------------------
  ! Optional per-iteration profiling of the MRSF z-vector (CPHF/CPKS) solve.
  ! Default OFF; enabled by setting env OQP_MRSF_ZV_TIMERS to a non-empty value.
  ! Wall-clock seconds are accumulated per section and reported at the end of
  ! each tdhf_mrsf_z_vector call (reset at the start of every call).
  ! ----------------------------------------------------------------------------
  logical, save :: zv_tmr_on   = .false.
  logical, save :: zv_tmr_init = .false.
  real(kind=dp), save :: zv_t_rhs   = 0.0_dp  !< RHS assembly
  real(kind=dp), save :: zv_t_int2  = 0.0_dp  !< per-iter 2e digestion (int2_driver%run)
  real(kind=dp), save :: zv_t_xc    = 0.0_dp  !< per-iter XC kernel (utddft_fxc)
  real(kind=dp), save :: zv_t_trans = 0.0_dp  !< per-iter AO/MO transforms + operator apply
  real(kind=dp), save :: zv_t_cg    = 0.0_dp  !< per-iter linear-solver vector algebra
  real(kind=dp), save :: zv_t_back  = 0.0_dp  !< relaxed-density / W back-projection
  integer, save       :: zv_n_iter  = 0       !< iterations performed by the chosen solver

  ! ----------------------------------------------------------------------------
  ! Performance controls for the MRSF z-vector solve (read once from env).
  !   OQP_MRSF_ZV_WARMSTART  : seed from the previous step's solution. DEFAULT ON
  !                            (=0/n/f to disable). Cannot change a converged
  !                            result -- only the iteration count.
  !   OQP_MRSF_ZV_PROG       : progressive (iteration-dependent) screening.
  !                            DEFAULT ON (=0/n/f to disable). Perturbs the
  !                            gradient by <~1e-8 (small systems) to ~7e-6 (large),
  !                            within the gradient gate but NOT bit-identical.
  !   OQP_MRSF_ZV_CONV       : override the convergence tol (default 1e-10 kept).
  !   OQP_MRSF_ZV_CUTOFF     : static loose 2e cutoff (off; superseded by _PROG).
  ! ----------------------------------------------------------------------------
  logical, save :: zv_cfg_init   = .false.
  logical, save :: zv_warm_on    = .false.
  real(kind=dp), save :: zv_conv_user   = -1.0_dp  !< <0 => not set
  real(kind=dp), save :: zv_cutoff_user = -1.0_dp  !< <0 => not set

  ! Progressive (iteration-dependent) integral screening for the CG solve.
  ! Loose cutoff while the residual is large (the search direction is inexact
  ! anyway), tightened toward the tight floor and PINNED exact once the residual
  ! is small, so the converged z-vector / gradient matches the all-tight run.
  ! tau(k) = clamp(prog_k * ||r_{k-1}||, tight_floor, prog_cap); tau=tight once
  ! ||r||^2 < prog_pin. Same coupling idea as feat/progressive-screening-scf.
  logical, save :: zv_prog_on  = .false.
  real(kind=dp), save :: zv_prog_k   = 1.0e-2_dp  !< tau = k*residual_norm
  real(kind=dp), save :: zv_prog_cap = 1.0e-6_dp  !< loosest tau (upper clamp)
  real(kind=dp), save :: zv_prog_pin = 1.0e-6_dp  !< pin tau=tight once error(=||r||^2) < this

  ! Cold-start initial guess: Jacobi x0 = M^-1 rhs instead of x0 = 0. The cold CG
  ! already spends one sigma build on A*x0 (=0 when x0=0), so this is free and
  ! puts the solve ~1 iteration ahead. DEFAULT ON; disable with OQP_MRSF_ZV_DIAGGUESS=0.
  ! Accuracy-safe (linear solve converges to the same A^-1 rhs) and protected by
  ! the CG safeguard (falls back to zero if it doesn't reduce the residual).
  logical, save :: zv_diag_guess = .true.

  ! Coarser DFT grid for the z-vector XC kernel than the SCF grid (the response
  ! tolerates a coarser grid). Shrinks the dominant per-iteration utddft_fxc cost.
  ! The grid is selected by the pruned-grid NAME (SG0 < SG1 < SG2(default) < SG3);
  ! the lever swaps to a coarser named grid for the z-vector only. DEFAULT OFF
  ! (env OQP_MRSF_ZV_COARSEGRID=1); grid choice via OQP_MRSF_ZV_GRID (default SG1).
  logical, save :: zv_coarse_on = .false.
  character(len=8), save :: zv_grid_name = 'SG1'

  ! Warm-start store (module-level; persists across the geometry/MD steps that
  ! share one process). The converged z-vector varies smoothly along a path, so
  ! the previous step's solution is a strong initial guess for the LINEAR CPHF
  ! solve. Correctness is unconditional: the iterative solver still converges to
  ! the same residual tolerance, so only the iteration count changes -- the
  ! guess can never bias the gradient. Keyed by target state; reset when the
  ! problem dimension (basis/occupation) changes.
  real(kind=dp), allocatable, save :: zv_warm(:,:)
  logical,       allocatable, save :: zv_warm_has(:)
  integer, save :: zv_warm_lzdim = 0
  ! Old MO coefficients from the step that filled the cache, for MO-basis
  ! projection of the warm guess (handles MO rotation between steps so warm-start
  ! also helps large optimizer steps, not just small MD steps).
  real(kind=dp), allocatable, save :: zv_warm_moa(:,:), zv_warm_mob(:,:)
  logical, save :: zv_warm_have_mo = .false.
  integer, save :: zv_warm_nbf = 0

contains

  !> Wall-clock seconds (monotonic), for the optional z-vector profiler.
  function zv_wtime() result(t)
    real(kind=dp) :: t
    integer(kind=8) :: c, r
    call system_clock(c, r)
    if (r > 0_8) then
      t = real(c, dp) / real(r, dp)
    else
      t = 0.0_dp
    end if
  end function zv_wtime

  !> Read the OQP_MRSF_ZV_TIMERS opt-in once and reset the accumulators.
  subroutine zv_timers_begin()
    character(len=8) :: e_
    if (.not. zv_tmr_init) then
      call get_environment_variable('OQP_MRSF_ZV_TIMERS', e_)
      zv_tmr_on = len_trim(e_) > 0
      zv_tmr_init = .true.
    end if
    zv_t_rhs = 0.0_dp; zv_t_int2 = 0.0_dp; zv_t_xc = 0.0_dp
    zv_t_trans = 0.0_dp; zv_t_cg = 0.0_dp; zv_t_back = 0.0_dp
    zv_n_iter = 0
  end subroutine zv_timers_begin

  !> Emit the accumulated per-section breakdown (no-op unless timers are on).
  subroutine zv_timers_report(log_unit, solver_name)
    integer, intent(in) :: log_unit
    character(len=*), intent(in) :: solver_name
    real(kind=dp) :: tot, sigma
    if (.not. zv_tmr_on) return
    sigma = zv_t_int2 + zv_t_xc + zv_t_trans
    tot = zv_t_rhs + sigma + zv_t_cg + zv_t_back
    write(log_unit,'(/1x,"==== MRSF Z-VECTOR PROFILE (",a,") ====")') trim(solver_name)
    write(log_unit,'(1x,"Iterations                 : ",i8)') zv_n_iter
    write(log_unit,'(1x,"RHS build            (s)   : ",f12.4)') zv_t_rhs
    write(log_unit,'(1x,"Per-iter 2e digestion(s)   : ",f12.4)') zv_t_int2
    write(log_unit,'(1x,"Per-iter XC kernel   (s)   : ",f12.4)') zv_t_xc
    write(log_unit,'(1x,"Per-iter transforms  (s)   : ",f12.4)') zv_t_trans
    write(log_unit,'(1x,"  -> sigma/Fock build(s)   : ",f12.4)') sigma
    write(log_unit,'(1x,"Per-iter CG algebra  (s)   : ",f12.4)') zv_t_cg
    write(log_unit,'(1x,"Back-projection      (s)   : ",f12.4)') zv_t_back
    write(log_unit,'(1x,"Sum of sections      (s)   : ",f12.4)') tot
    if (zv_n_iter > 0) &
      write(log_unit,'(1x,"Avg sigma / iteration(s)   : ",f12.4)') sigma/real(zv_n_iter,dp)
    write(log_unit,'(1x,"=========================================")')
    call flush(log_unit)
  end subroutine zv_timers_report

  !> Read the z-vector opt-in environment variables once.
  subroutine zv_read_config()
    character(len=32) :: e_
    integer :: ios
    if (zv_cfg_init) return
    ! Warm-start is ON by default; disable with OQP_MRSF_ZV_WARMSTART=0/n/f/off.
    call get_environment_variable('OQP_MRSF_ZV_WARMSTART', e_)
    zv_warm_on = .true.
    if (len_trim(e_) > 0) zv_warm_on = .not. (e_(1:1)=='0' .or. e_(1:1)=='n' .or. &
        e_(1:1)=='N' .or. e_(1:1)=='f' .or. e_(1:1)=='F')
    call get_environment_variable('OQP_MRSF_ZV_CONV', e_)
    if (len_trim(e_) > 0) then
      read(e_,*,iostat=ios) zv_conv_user
      if (ios /= 0) zv_conv_user = -1.0_dp
    end if
    call get_environment_variable('OQP_MRSF_ZV_CUTOFF', e_)
    if (len_trim(e_) > 0) then
      read(e_,*,iostat=ios) zv_cutoff_user
      if (ios /= 0) zv_cutoff_user = -1.0_dp
    end if
    ! Progressive screening is ON by default; disable with OQP_MRSF_ZV_PROG=0/n/f.
    call get_environment_variable('OQP_MRSF_ZV_PROG', e_)
    zv_prog_on = .true.
    if (len_trim(e_) > 0) zv_prog_on = .not. (e_(1:1)=='0' .or. e_(1:1)=='n' .or. &
        e_(1:1)=='N' .or. e_(1:1)=='f' .or. e_(1:1)=='F')
    call get_environment_variable('OQP_MRSF_ZV_PROG_K', e_)
    if (len_trim(e_) > 0) read(e_,*,iostat=ios) zv_prog_k
    call get_environment_variable('OQP_MRSF_ZV_PROG_CAP', e_)
    if (len_trim(e_) > 0) read(e_,*,iostat=ios) zv_prog_cap
    call get_environment_variable('OQP_MRSF_ZV_PROG_PIN', e_)
    if (len_trim(e_) > 0) read(e_,*,iostat=ios) zv_prog_pin
    ! Jacobi cold-start guess: default ON; disable with OQP_MRSF_ZV_DIAGGUESS=0/n/f.
    call get_environment_variable('OQP_MRSF_ZV_DIAGGUESS', e_)
    if (len_trim(e_) > 0) zv_diag_guess = .not. (e_(1:1)=='0' .or. e_(1:1)=='n' .or. &
        e_(1:1)=='N' .or. e_(1:1)=='f' .or. e_(1:1)=='F')
    ! Coarser response grid (default OFF; new feature, validate before defaulting).
    call get_environment_variable('OQP_MRSF_ZV_COARSEGRID', e_)
    zv_coarse_on = len_trim(e_) > 0 .and. (e_(1:1)=='1' .or. e_(1:1)=='y' .or. &
        e_(1:1)=='Y' .or. e_(1:1)=='t' .or. e_(1:1)=='T')
    call get_environment_variable('OQP_MRSF_ZV_GRID', e_)
    if (len_trim(e_) > 0) zv_grid_name = trim(adjustl(e_))
    zv_cfg_init = .true.
  end subroutine zv_read_config

  !> Progressive screening threshold for a CG step given the current squared
  !> residual `error` and the tight floor. Returns the tight floor once pinned.
  function zv_prog_tau(error, tight) result(tau)
    real(kind=dp), intent(in) :: error, tight
    real(kind=dp) :: tau, resid
    if (.not. ieee_is_finite(error) .or. error < zv_prog_pin) then
      tau = tight
      return
    end if
    resid = sqrt(max(error, 0.0_dp))
    tau = zv_prog_k * resid
    if (tau < tight) tau = tight
    if (tau > zv_prog_cap) tau = zv_prog_cap
  end function zv_prog_tau

  !> Store a converged xk (and the current MOs, for later MO-basis projection)
  !> into the warm-start store for `state`. Reallocates if the problem dimension
  !> changed (e.g. a different basis/molecule).
  subroutine zv_store_guess(xk, lzdim, state, nstate, mo_a, mo_b, nbf)
    real(kind=dp), intent(in) :: xk(:)
    integer, intent(in) :: lzdim, state, nstate, nbf
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:)
    integer :: ncol
    if (.not. zv_warm_on) return
    if (any(.not. ieee_is_finite(xk))) return
    ncol = max(nstate, state)
    if (zv_warm_lzdim /= lzdim .or. .not. allocated(zv_warm) .or. &
        .not. allocated(zv_warm_has)) then
      if (allocated(zv_warm))     deallocate(zv_warm)
      if (allocated(zv_warm_has)) deallocate(zv_warm_has)
      allocate(zv_warm(lzdim, ncol), source=0.0_dp)
      allocate(zv_warm_has(ncol), source=.false.)
      zv_warm_lzdim = lzdim
    else if (size(zv_warm_has) < state) then
      return
    end if
    zv_warm(:,state) = xk
    zv_warm_has(state) = .true.
    ! Cache the current MOs (shared across states this step) for MO projection.
    if (zv_warm_nbf /= nbf .or. .not. allocated(zv_warm_moa)) then
      if (allocated(zv_warm_moa)) deallocate(zv_warm_moa)
      if (allocated(zv_warm_mob)) deallocate(zv_warm_mob)
      allocate(zv_warm_moa(nbf,nbf), zv_warm_mob(nbf,nbf))
      zv_warm_nbf = nbf
    end if
    zv_warm_moa(:,:) = mo_a(1:nbf,1:nbf)
    zv_warm_mob(:,:) = mo_b(1:nbf,1:nbf)
    zv_warm_have_mo = .true.
  end subroutine zv_store_guess

  !> Inverse of sfrogen: gather the z-vector amplitudes from MO-basis matrices
  !> ava (alpha) / avb (beta) at the same index positions sfrogen scatters to.
  !> The doc-virt block (which sfrogen writes to both spins) is averaged.
  subroutine zv_sfrogen_gather(ava, avb, xk, nocca, noccb)
    real(kind=dp), intent(in)  :: ava(:,:), avb(:,:)
    real(kind=dp), intent(out) :: xk(:)
    integer, intent(in) :: nocca, noccb
    integer :: ij, i, j, k, nbf
    nbf = ubound(ava,1)
    ij = 0
    do i = noccb+1, nocca        ! doc-socc  (beta)
      do j = 1, noccb
        ij = ij+1; xk(ij) = avb(j,i)
      end do
    end do
    do k = nocca+1, nbf          ! doc-virt  (sfrogen wrote both spins)
      do j = 1, noccb
        ij = ij+1; xk(ij) = 0.5_dp*(ava(j,k)+avb(j,k))
      end do
    end do
    do k = nocca+1, nbf          ! socc-virt (alpha)
      do i = noccb+1, nocca
        ij = ij+1; xk(ij) = ava(i,k)
      end do
    end do
  end subroutine zv_sfrogen_gather

  ! Initialize GMRES work arrays
  subroutine init_gmres_work(nbf, nocca, noccb)
    use messages, only: show_message, with_abort
    implicit none
    integer, intent(in) :: nbf, nocca, noccb
    integer :: nvira, nvirb, ok
    
    if (nbf <= 0 .or. nocca < 0 .or. noccb < 0) then
      call show_message('Invalid GMRES work dimensions before allocation', with_abort)
    end if

    nvira = nbf - nocca
    nvirb = nbf - noccb
    if (nvira <= 0 .or. nvirb <= 0) then
      call show_message('Invalid GMRES occupied/virtual dimensions before allocation', with_abort)
    end if

    if (gmres_work_allocated) then
      ! Check if dimensions match and every reusable array is still allocated.
      if (gmres_nbf == nbf .and. gmres_nocca == nocca .and. gmres_noccb == noccb) then
        if (gmres_work_arrays_allocated()) then
          return  ! Arrays already allocated with correct size
        end if
        call cleanup_gmres_work()
      else
        ! Deallocate old arrays before reallocating
        call cleanup_gmres_work()
      end if
    end if
    
    allocate(gmres_wrk1(nbf,nbf), &
             gmres_wrk2(nbf,nbf), &
             gmres_wrk3(nbf,nbf), &
             gmres_pa(nbf,nbf,2), &
             gmres_ab1_mo_a(nocca,nvira), &
             gmres_ab1_mo_b(noccb,nvirb), &
             stat=ok)
    
    if (ok /= 0) then
      call show_message('Cannot allocate GMRES work arrays', with_abort)
    end if
    
    gmres_work_allocated = .true.
    gmres_nbf = nbf
    gmres_nocca = nocca
    gmres_noccb = noccb
    
  end subroutine init_gmres_work
  
  ! Return true only when every reusable GMRES work array is allocated.
  logical function gmres_work_arrays_allocated()
    implicit none

    gmres_work_arrays_allocated = allocated(gmres_wrk1) .and. &
                                  allocated(gmres_wrk2) .and. &
                                  allocated(gmres_wrk3) .and. &
                                  allocated(gmres_pa) .and. &
                                  allocated(gmres_ab1_mo_a) .and. &
                                  allocated(gmres_ab1_mo_b)
  end function gmres_work_arrays_allocated

  ! Cleanup GMRES work arrays
  subroutine cleanup_gmres_work()
    implicit none
    
    if (allocated(gmres_wrk1)) deallocate(gmres_wrk1)
    if (allocated(gmres_wrk2)) deallocate(gmres_wrk2)
    if (allocated(gmres_wrk3)) deallocate(gmres_wrk3)
    if (allocated(gmres_pa)) deallocate(gmres_pa)
    if (allocated(gmres_ab1_mo_a)) deallocate(gmres_ab1_mo_a)
    if (allocated(gmres_ab1_mo_b)) deallocate(gmres_ab1_mo_b)
    
    gmres_work_allocated = .false.
    gmres_nbf = 0
    gmres_nocca = 0
    gmres_noccb = 0
    
  end subroutine cleanup_gmres_work

  ! GMRES solver for the z-vector equation
  subroutine gmres_solve(apply_operator, apply_precond, b, x, n, restart, max_iter, tol, &
                         infos, basis, molGrid, int2_driver, &
                         nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                         fa, fb, scale_exch, dft, error_out, iter_out, iw)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use mod_dft_molgrid, only: dft_grid_t
    
    implicit none
    
    interface
      subroutine apply_operator(x_in, x_out, infos, basis, molGrid, int2_driver, &
                               nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                               fa, fb, scale_exch, dft)
        use precision, only: dp
        use types, only: information
        use basis_tools, only: basis_set
        use int2_compute, only: int2_compute_t
        use mod_dft_molgrid, only: dft_grid_t
        real(kind=dp), intent(in) :: x_in(:)
        real(kind=dp), intent(out) :: x_out(:)
        type(information), intent(inout) :: infos
        type(basis_set), pointer :: basis
        type(dft_grid_t), intent(inout) :: molGrid
        type(int2_compute_t), intent(inout) :: int2_driver
        integer, intent(in) :: nocca, noccb, nbf
        real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
        real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
        logical, intent(in) :: dft
      end subroutine
      
      subroutine apply_precond(x_in, x_out)
        use precision, only: dp
        real(kind=dp), intent(in) :: x_in(:)
        real(kind=dp), intent(out) :: x_out(:)
      end subroutine
    end interface
    
    real(kind=dp), intent(in) :: b(:)
    real(kind=dp), intent(inout) :: x(:)
    integer, intent(in) :: n, restart, max_iter, iw
    real(kind=dp), intent(in) :: tol
    type(information), intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(dft_grid_t), intent(inout) :: molGrid
    type(int2_compute_t), intent(inout) :: int2_driver
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
    real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
    logical, intent(in) :: dft
    real(kind=dp), intent(out) :: error_out
    integer, intent(out) :: iter_out
    
    ! Local variables
    real(kind=dp), allocatable :: V(:,:)     ! Krylov basis
    real(kind=dp), allocatable :: H(:,:)     ! Hessenberg matrix
    real(kind=dp), allocatable :: c(:), s(:) ! Givens rotation coefficients
    real(kind=dp), allocatable :: g(:)       ! RHS for least squares
    real(kind=dp), allocatable :: y(:)       ! Solution of least squares
    real(kind=dp), allocatable :: r(:)       ! Residual
    real(kind=dp), allocatable :: w(:)       ! Work vector
    real(kind=dp), allocatable :: Ax(:)      ! A*x
    
    real(kind=dp) :: beta, h_ij, temp, error, error_initial, true_residual
    real(kind=dp) :: max_abs_overlap
    real(kind=dp), parameter :: gmres_reorth_threshold = 1.0d-10
    integer :: i, j, k, iter, m, restart_count, inner_iter
    logical :: converged, unstable, happy_breakdown
    
    ! Initialize GMRES work arrays ONCE at the beginning
    if (n <= 0 .or. restart <= 0) then
      write(iw,'(" GMRES: invalid dimensions provided (n/restart)")')
      error_out = huge(1.0_dp)
      iter_out = 0
      return
    end if

    if (size(b) /= n .or. size(x) /= n) then
      write(iw,'(" GMRES: vector size does not match problem size")')
      error_out = huge(1.0_dp)
      iter_out = 0
      return
    end if

    if (size(mo_a,1) /= nbf .or. size(mo_b,1) /= nbf) then
      write(iw,'(" GMRES: invalid basis-size arguments")')
      error_out = huge(1.0_dp)
      iter_out = 0
      return
    end if

    call init_gmres_work(nbf, nocca, noccb)
    
    ! Allocate workspace
    m = min(restart, n)
    allocate(V(n, m+1))
    allocate(H(m+1, m))
    allocate(c(m))
    allocate(s(m))
    allocate(g(m+1))
    allocate(y(m))
    allocate(r(n))
    allocate(w(n))
    allocate(Ax(n))
    
    iter_out = 0
    restart_count = 0
    converged = .false.
    unstable = .false.
    happy_breakdown = .false.
    
    write(iw,'(/," GMRES Solver Parameters:")')
    write(iw,'("   Problem size        : ", I8)') n
    write(iw,'("   Restart dimension   : ", I4)') m
    write(iw,'("   Max iterations      : ", I4)') max_iter
    write(iw,'("   Convergence tol     : ", 1p,e10.3)') tol
    write(iw,'(/," Iteration   Inner  Residual Norm   Reduction")')
    write(iw,'(" ---------   -----  -------------   ---------")')
    call flush(iw)
    
    ! Compute initial residual
    call apply_operator(x, Ax, infos, basis, molGrid, int2_driver, &
                       nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                       fa, fb, scale_exch, dft)
    if (any(.not. ieee_is_finite(Ax)) .or. any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(b))) then
      write(iw,'(" GMRES: initial residual has non-finite input")')
      call flush(iw)
      error_out = huge(1.0_dp)
      error_initial = huge(1.0_dp)
      iter_out = 0
      unstable = .true.
    else
    r = b - Ax
    error_initial = sqrt(dot_product(r, r))
    if (error_initial == 0.0_dp) then
      write(iw,'(" GMRES: initial residual is exactly zero")')
      call flush(iw)
      error_out = error_initial
      iter_out = 0
      converged = .true.
    else
    if (.not. ieee_is_finite(error_initial)) then
      write(iw,'(" GMRES: initial residual is non-finite; aborting iterative safety")')
      call flush(iw)
      error_out = huge(1.0_dp)
      iter_out = 0
      unstable = .true.
    else if (error_initial < GMRES_DENOMINATOR_FLOOR) then
      write(iw,'(" GMRES: initial residual is already below denominator floor")')
      call flush(iw)
      error_out = error_initial
      iter_out = 0
      converged = .true.
      end if
    end if
    end if
    
    if (unstable .or. converged) then
      ! Safety fallback or exact-zero residual handled above; skip GMRES iterations
    else
      do iter = 1, max_iter
      
      restart_count = restart_count + 1
      
      ! Compute initial residual r = b - A*x
      call apply_operator(x, Ax, infos, basis, molGrid, int2_driver, &
                         nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                         fa, fb, scale_exch, dft)
      if (any(.not. ieee_is_finite(Ax)) .or. any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(b))) then
        write(iw,'(" GMRES: non-finite values during restart residual setup")')
        unstable = .true.
        exit
      end if
      r = b - Ax
      if (any(.not. ieee_is_finite(r))) then
        write(iw,'(" GMRES: non-finite residual during restart")')
        unstable = .true.
        exit
      end if
      true_residual = sqrt(dot_product(r, r))
      if (.not. ieee_is_finite(true_residual)) then
        write(iw,'(" GMRES: non-finite true residual during restart")')
        unstable = .true.
        exit
      end if
      if (true_residual < tol) then
        error_out = true_residual
        error = true_residual
        converged = .true.
        exit
      end if
      
      ! Apply preconditioner to residual
      call apply_precond(r, V(:,1))
      if (any(.not. ieee_is_finite(V(:,1))) .or. any(.not. ieee_is_finite(r))) then
        write(iw,'(" GMRES: non-finite preconditioned residual")')
        unstable = .true.
        exit
      end if
      
      beta = sqrt(dot_product(V(:,1), V(:,1)))
      if (any(.not. ieee_is_finite(V(:,1)))) then
        write(iw,'(" GMRES: non-finite basis vector prior to scaling")')
        unstable = .true.
        exit
      end if
      if (.not. ieee_is_finite(beta) .or. beta < GMRES_DENOMINATOR_FLOOR) then
        write(iw,'(" GMRES: degenerate residual norm at restart ", I3)') restart_count
        unstable = .true.
        exit
      end if
      
      ! Report the true residual; beta is only the preconditioned Arnoldi seed norm.
      error = true_residual
      if (iter == 1) then
        write(iw,'(I6,8x,"  0",2x,1p,F13.8,1x,F13.8)') &
              restart_count, error, error/error_initial
      end if
      
      V(:,1) = V(:,1) / beta
      g = 0.0_dp
      g(1) = beta
      
      ! Reset H matrix for this restart
      H = 0.0_dp
      
      ! Arnoldi process
      inner_iter = 0
      happy_breakdown = .false.
      do j = 1, m
        inner_iter = j
        
        ! Apply operator to V_j
        call apply_operator(V(:,j), w, infos, basis, molGrid, int2_driver, &
                           nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                           fa, fb, scale_exch, dft)
        if (any(.not. ieee_is_finite(w)) .or. any(.not. ieee_is_finite(V(:,j)))) then
          write(iw,'(" GMRES: non-finite basis/operator values at inner step", I3)') j
          unstable = .true.
          exit
        end if
        
        ! Apply preconditioner
        call apply_precond(w, V(:,j+1))
        if (any(.not. ieee_is_finite(V(:,j+1))) .or. any(.not. ieee_is_finite(w))) then
          write(iw,'(" GMRES: non-finite preconditioned vector at inner step", I3)') j
          unstable = .true.
          exit
        end if
        
        ! Modified Gram-Schmidt orthogonalization
        do i = 1, j
          H(i,j) = dot_product(V(:,j+1), V(:,i))
          if (.not. ieee_is_finite(H(i,j))) then
            unstable = .true.
            exit
          end if
          V(:,j+1) = V(:,j+1) - H(i,j) * V(:,i)
        end do
        if (unstable) then
          write(iw,'(" GMRES: non-finite H entry during orthogonalization")')
          exit
        end if
        
        ! Reorthogonalize only when residual overlap drift is measurable.
        max_abs_overlap = 0.0_dp
        do i = 1, j
          temp = dot_product(V(:,j+1), V(:,i))
          if (.not. ieee_is_finite(temp)) then
            unstable = .true.
            exit
          end if
          max_abs_overlap = max(max_abs_overlap, abs(temp))
        end do
        if (unstable) then
          write(iw,'(" GMRES: non-finite overlap during re-orthogonalization check")')
          exit
        end if
        if (max_abs_overlap > gmres_reorth_threshold) then
          ! GMRES MGS reorthogonalization pass
          do i = 1, j
            temp = dot_product(V(:,j+1), V(:,i))
            H(i,j) = H(i,j) + temp
            if (.not. ieee_is_finite(temp) .or. .not. ieee_is_finite(H(i,j))) then
              unstable = .true.
              exit
            end if
            V(:,j+1) = V(:,j+1) - temp * V(:,i)
          end do
          if (unstable) then
            write(iw,'(" GMRES: non-finite H entry during re-orthogonalization")')
            exit
          end if
        end if
        
        H(j+1,j) = sqrt(dot_product(V(:,j+1), V(:,j+1)))
        if (any(.not. ieee_is_finite(V(:,j+1)))) then
          unstable = .true.
          exit
        end if
        
        ! Check for breakdown
        if (.not. ieee_is_finite(H(j+1,j))) then
          write(iw,'(" GMRES: non-finite Arnoldi norm at iteration ", I3)') j
          inner_iter = j
          unstable = .true.
          exit
        else if (abs(H(j+1,j)) < GMRES_DENOMINATOR_FLOOR) then
          write(iw,'(" GMRES: happy breakdown at iteration ", I3, "; recomputing true residual")') j
          inner_iter = j
          H(j+1,j) = 0.0_dp
          happy_breakdown = .true.
        else
          V(:,j+1) = V(:,j+1) / H(j+1,j)
        end if
        
        ! Apply previous Givens rotations
        do i = 1, j-1
          temp = c(i) * H(i,j) + s(i) * H(i+1,j)
          H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j)
          H(i,j) = temp
        end do
        
        if (.not. ieee_is_finite(H(j,j)) .or. .not. ieee_is_finite(H(j+1,j))) then
          unstable = .true.
          exit
        end if
        ! Compute new Givens rotation
        call givens_rotation(H(j,j), H(j+1,j), c(j), s(j))
        
        ! Apply new Givens rotation
        H(j,j) = c(j) * H(j,j) + s(j) * H(j+1,j)
        H(j+1,j) = 0.0_dp
        
        temp = c(j) * g(j) + s(j) * g(j+1)
        g(j+1) = -s(j) * g(j) + c(j) * g(j+1)
        g(j) = temp
        
        ! Check convergence
        error = abs(g(j+1))
        iter_out = iter_out + 1
        
        ! Print progress every 5 inner iterations or at convergence
        if (mod(j, 5) == 0 .or. error < tol .or. j == m) then
          write(iw,'(I6,8x,I3,2x,1p,F13.8,1x,F13.8)') &
                restart_count, j, error, error/error_initial
          call flush(iw)
        end if
        
        if (error < tol) then
          converged = .true.
          inner_iter = j
          exit
        end if

        if (happy_breakdown) then
          inner_iter = j
          exit
        end if
        
        if (iter_out >= max_iter) then
          inner_iter = j
          exit
        end if
      end do
      
      if (inner_iter > 0 .and. .not. unstable) then
        ! Solve upper triangular system for y
        call back_substitution(H(1:inner_iter,1:inner_iter), g(1:inner_iter), &
                               y(1:inner_iter), inner_iter, unstable)
      end if
      
      if (unstable) then
        error_out = huge(1.0_dp)
        exit
      end if

      ! Update solution: x = x + V*y
      do i = 1, inner_iter
        x = x + y(i) * V(:,i)
      end do
      if (any(.not. ieee_is_finite(x))) then
        write(iw,'(" GMRES: non-finite solution update")')
        unstable = .true.
        error_out = huge(1.0_dp)
        exit
      end if

      call recompute_gmres_true_residual(true_residual, unstable)
      if (unstable) then
        error_out = huge(1.0_dp)
        exit
      end if
      error_out = true_residual
      converged = true_residual < tol
      error = true_residual
      
      if (converged .or. iter_out >= max_iter) exit
      
      ! Print restart information
      if (.not. converged .and. inner_iter == m) then
        write(iw,'(" GMRES: Restarting (restart #", I3, ")")') restart_count
        call flush(iw)
      end if
      
    end do
    end if
    
    ! Final status
    write(iw,'(" ---------   -----  -------------   ---------")')
    if (converged) then
      write(iw,'(" GMRES converged in ", I4, " iterations (", I3, " restarts)")') &
            iter_out, restart_count-1
      write(iw,'(" Final residual norm: ", 1p,e13.6)') error_out
    else
      if (unstable) then
        write(iw,'(" GMRES terminated due to numerical instability")')
      else
        write(iw,'(" GMRES did not converge within ", I4, " iterations")') max_iter
      end if
      write(iw,'(" Final residual norm: ", 1p,e13.6)') error_out
    end if
    if (ieee_is_finite(error_initial) .and. abs(error_initial) >= GMRES_DENOMINATOR_FLOOR) then
      write(iw,'(" Relative reduction : ", 1p,e13.6)') error_out/error_initial
    else
      write(iw,'(" Relative reduction : not available")')
    end if
    call flush(iw)
    
    ! Clean up local arrays
    deallocate(V, H, c, s, g, y, r, w, Ax)
    
    ! NOTE: Do NOT clean up GMRES work arrays here - they will be cleaned in main routine
    
  contains
    
    subroutine recompute_gmres_true_residual(true_residual, unstable)
      real(kind=dp), intent(out) :: true_residual
      logical, intent(out) :: unstable

      call apply_operator(x, Ax, infos, basis, molGrid, int2_driver, &
                         nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                         fa, fb, scale_exch, dft)
      if (any(.not. ieee_is_finite(Ax)) .or. any(.not. ieee_is_finite(x)) .or. &
          any(.not. ieee_is_finite(b))) then
        write(iw,'(" GMRES: non-finite values during true residual recomputation")')
        true_residual = huge(1.0_dp)
        unstable = .true.
        return
      end if

      r = b - Ax
      if (any(.not. ieee_is_finite(r))) then
        write(iw,'(" GMRES: non-finite true residual vector")')
        true_residual = huge(1.0_dp)
        unstable = .true.
        return
      end if

      true_residual = sqrt(dot_product(r, r))
      if (.not. ieee_is_finite(true_residual)) then
        write(iw,'(" GMRES: non-finite true residual norm")')
        true_residual = huge(1.0_dp)
        unstable = .true.
        return
      end if

      unstable = .false.
    end subroutine recompute_gmres_true_residual

    subroutine givens_rotation(a, b, c, s)
      real(kind=dp), intent(in) :: a, b
      real(kind=dp), intent(out) :: c, s
      real(kind=dp) :: r, scale
      
      if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b) .or. &
          (abs(a) < GMRES_DENOMINATOR_FLOOR .and. abs(b) < GMRES_DENOMINATOR_FLOOR)) then
        c = 1.0_dp
        s = 0.0_dp
      else
        scale = max(abs(a), abs(b))
        r = scale * sqrt((a/scale)**2 + (b/scale)**2)
        c = a / r
        s = b / r
      end if
    end subroutine givens_rotation
    
    subroutine back_substitution(A, b, x, n, unstable)
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: A(n,n), b(n)
      real(kind=dp), intent(out) :: x(n)
      logical, intent(out) :: unstable
      integer :: i, j
      real(kind=dp) :: rhs

      if (n <= 0) then
        unstable = .true.
        return
      end if

      if (.not. ieee_is_finite(A(n,n)) .or. abs(A(n,n)) < GMRES_DENOMINATOR_FLOOR) then
        unstable = .true.
        return
      end if
      if (.not. ieee_is_finite(b(n))) then
        unstable = .true.
        return
      end if

      x(n) = b(n) / A(n,n)
      if (.not. ieee_is_finite(x(n))) then
        unstable = .true.
        return
      end if
      do i = n-1, 1, -1
        if (.not. ieee_is_finite(A(i,i)) .or. abs(A(i,i)) < GMRES_DENOMINATOR_FLOOR) then
          unstable = .true.
          return
        end if

        if (.not. ieee_is_finite(b(i))) then
          unstable = .true.
          return
        end if
        rhs = b(i)
        do j = i+1, n
          if (.not. ieee_is_finite(A(i,j))) then
            unstable = .true.
            return
          end if
          rhs = rhs - A(i,j) * x(j)
        end do
        if (.not. ieee_is_finite(rhs)) then
          unstable = .true.
          return
        end if
        x(i) = rhs / A(i,i)
        if (.not. ieee_is_finite(x(i))) then
          unstable = .true.
          return
        end if
      end do
      unstable = .false.
    end subroutine back_substitution

    
  end subroutine gmres_solve

  ! Apply the z-vector operator (A*x) - Modified to use module-level arrays
  subroutine apply_z_operator(x_in, x_out, infos, basis, molGrid, int2_driver, &
                              nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                              fa, fb, scale_exch, dft)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_tdgrd_data_t
    use tdhf_sf_lib, only: sfrogen, sfrolhs
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t
    use tdhf_lib, only: mntoia
    
    implicit none
    
    real(kind=dp), intent(in) :: x_in(:)
    real(kind=dp), intent(out) :: x_out(:)
    type(information), intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(dft_grid_t), intent(inout) :: molGrid
    type(int2_compute_t), intent(inout) :: int2_driver
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
    real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
    logical, intent(in) :: dft
    
    ! Local variables
    real(kind=dp), pointer :: ab1(:,:,:)
    type(int2_tdgrd_data_t), allocatable, target :: int2_data
    integer :: nvira, nvirb
    real(kind=dp) :: t0_op

    nvira = nbf - nocca
    nvirb = nbf - noccb

    if (any(.not. ieee_is_finite(x_in))) then
      x_out = ieee_value(0.0_dp, ieee_quiet_nan)
      write(*,'(" MRSF z-vector operator rejected non-finite input")')
      return
    end if
    
    ! Ensure work arrays are initialized and have correct dimensions
    if (.not. gmres_work_allocated .or. &
        gmres_nbf /= nbf .or. &
        gmres_nocca /= nocca .or. &
        gmres_noccb /= noccb) then
      call init_gmres_work(nbf, nocca, noccb)
    end if
    
    ! Clear work arrays
    gmres_wrk1 = 0.0_dp
    gmres_wrk2 = 0.0_dp
    gmres_wrk3 = 0.0_dp
    gmres_pa = 0.0_dp
    gmres_ab1_mo_a = 0.0_dp
    gmres_ab1_mo_b = 0.0_dp
    
    ! Generate density matrices from x_in
    if (zv_tmr_on) then
      t0_op = zv_wtime()
      zv_n_iter = zv_n_iter + 1
    end if
    call sfrogen(gmres_wrk1, gmres_wrk2, x_in, nocca, noccb)

    ! Transform to AO basis
    call orthogonal_transform('t', nbf, mo_a, gmres_wrk1, gmres_pa(:,:,1), gmres_wrk3)
    call orthogonal_transform('t', nbf, mo_b, gmres_wrk2, gmres_pa(:,:,2), gmres_wrk3)
    if (zv_tmr_on) zv_t_trans = zv_t_trans + (zv_wtime() - t0_op)

    ! Initialize ERI calculation with proper allocation
    allocate(int2_data)
    int2_data = int2_tdgrd_data_t( &
        d2 = gmres_pa, &
        int_apb = .true., &
        int_amb = .false., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)

    if (zv_tmr_on) t0_op = zv_wtime()
    call int2_driver%run(int2_data, &
          cam=dft.and.infos%dft%cam_flag, &
          alpha=infos%dft%cam_alpha, &
          beta=infos%dft%cam_beta,&
          mu=infos%dft%cam_mu)
    if (zv_tmr_on) zv_t_int2 = zv_t_int2 + (zv_wtime() - t0_op)
    ab1 => int2_data%apb(:,:,:,1)

    call symmetrize_matrix(gmres_pa(:,:,1), nbf)
    call symmetrize_matrix(gmres_pa(:,:,2), nbf)

    if (dft) then
      if (zv_tmr_on) t0_op = zv_wtime()
      call utddft_fxc( &
          basis = basis, &
          molGrid = molGrid, &
          isVecs = .true., &
          wfa = mo_a, &
          wfb = mo_b, &
          fxa = ab1(:,:,1:1), &
          fxb = ab1(:,:,2:2), &
          dxa = gmres_pa(:,:,1:1), &
          dxb = gmres_pa(:,:,2:2), &
          nmtx = 1, &
          threshold = 1.0d-15, &
          infos = infos)
      if (zv_tmr_on) zv_t_xc = zv_t_xc + (zv_wtime() - t0_op)
    end if

    if (any(.not. ieee_is_finite(ab1))) then
      x_out = ieee_value(0.0_dp, ieee_quiet_nan)
      write(*,'(" MRSF z-vector operator rejected non-finite response")')
      call int2_data%clean()
      deallocate(int2_data)
      return
    end if
    
    ! Transform to MO basis - Fixed to use correct mo_b for beta
    if (zv_tmr_on) t0_op = zv_wtime()
    call mntoia(ab1(:,:,1), gmres_ab1_mo_a, mo_a, mo_a, nocca, nocca)
    call mntoia(ab1(:,:,2), gmres_ab1_mo_b, mo_b, mo_b, noccb, noccb)

    ! Apply the operator
    call sfrolhs(x_out, x_in, mo_energy_a, fa, fb, gmres_ab1_mo_a, gmres_ab1_mo_b, &
                 nocca, noccb)
    if (zv_tmr_on) zv_t_trans = zv_t_trans + (zv_wtime() - t0_op)

    call int2_data%clean()
    deallocate(int2_data)
    
  end subroutine apply_z_operator

  ! Apply preconditioner (simple diagonal preconditioner)
  subroutine apply_z_precond(x_in, x_out, xminv)
    use precision, only: dp
    implicit none
    real(kind=dp), intent(in) :: x_in(:), xminv(:)
    real(kind=dp), intent(out) :: x_out(:)

    if (any(.not. ieee_is_finite(x_in)) .or. any(.not. ieee_is_finite(xminv))) then
      x_out = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if

    x_out = xminv * x_in

  end subroutine apply_z_precond

  subroutine tdhf_mrsf_z_vector_C(c_handle) bind(C, name="tdhf_mrsf_z_vector")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_z_vector(inf)
  end subroutine tdhf_mrsf_z_vector_C

  subroutine tdhf_mrsf_z_vector(infos)
    use precision, only: dp
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use util, only: measure_time

    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t
    use tdhf_lib, only: int2_tdgrd_data_t
    use tdhf_mrsf_lib, only: int2_mrsf_data_t
    use tdhf_lib, only: iatogen, mntoia
    use tdhf_sf_lib, only: sfrorhs, &
      sfromcal, sfrogen, sfrolhs, pcgrbpini, &
      pcgb, sfropcal, sfdmat
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform, &
            orthogonal_transform_sym
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: pack_matrix, unpack_matrix

    use tdhf_mrsf_lib, only: &
      mrinivec, mrsfcbc, mrsfxvec, mrsfsp, mrsfrowcal, &
      mrsfqrorhs, mrsfqropcal, mrsfqrowcal
    use oqp_linalg
    use printing, only: print_module_info
    use minres_mod, only: minres_t, MINRES_OK, MINRES_CONVERGED


    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_mrsf_z_vector"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: ok

    real(kind=dp), allocatable :: ab1_mo_a(:,:)
    real(kind=dp), allocatable :: ab1_mo_b(:,:)
    real(kind=dp), allocatable :: xm(:)
    real(kind=dp), pointer :: ab1(:,:,:)
    real(kind=dp), allocatable :: fa(:,:), fb(:,:)
    real(kind=dp), pointer :: bvec(:,:,:)
    real(kind=dp), pointer :: wmo(:,:)
    real(kind=dp), allocatable :: bvec_mo_d(:,:)
    real(kind=dp), allocatable, target :: &
      fmrst1(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)

    integer :: nocca, nvira, noccb, nvirb
    integer :: nbf, nbf_tri
    integer :: iter, gmres_iter
    type(minres_t) :: mr
    integer :: minres_iter
    integer, target :: minres_dummy
    real(kind=dp) :: cnvtol, scale_exch, scale_exch2
    logical :: roref = .false.
    integer :: mrst

    type(int2_compute_t) :: int2_driver
    type(int2_mrsf_data_t), allocatable, target :: int2_data_st
    type(int2_td_data_t), allocatable, target :: int2_data_q
    class(int2_td_data_t), allocatable, target :: int2_data
    type(dft_grid_t) :: molGrid

  ! scr data
    real(kind=dp), allocatable, target :: wrk1(:,:), wrk2(:,:), wrk3(:,:)
    real(kind=dp), pointer :: wrk1t(:)

  ! SF-TD Gradient data
    real(kind=dp), allocatable :: &
      rhs(:), lhs(:), xminv(:), xk(:), pk(:), errv(:), &
      hxa(:,:), hxb(:,:), tij(:,:), ppija(:,:), ppijb(:,:), tab(:,:)
    real(kind=dp), allocatable, target :: pa(:,:,:)
    integer :: nsocc, lzdim, xvec_dim

  ! General data
    real(kind=dp) :: alpha, error, pap
    real(kind=dp) :: t0_zv
    real(kind=dp) :: zv_rc_save
    logical :: zv_rc_loosen
    character :: zv_gname_save(16)
    integer :: ig
    character(len=10) :: solver_name

    logical :: dft, mrsf_zvector_breakdown
    integer :: scf_type, mol_mult, target_state

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      fock_a(:), mo_a(:,:), mo_energy_a(:), &
      fock_b(:), mo_b(:,:), &
      td_p(:,:), td_t(:,:), ta(:), tb(:), td_abxc(:,:), &
      td_mrsf_den(:,:,:), bvec_mo(:,:), wao(:), mrsf_energies(:)
    character(len=*), parameter :: tags_alloc(4) = (/ character(len=80) :: &
      OQP_WAO, OQP_td_mrsf_density, OQP_td_p, OQP_td_abxc /)
    character(len=*), parameter :: tags_required(8) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_E_MO_A, OQP_VEC_MO_A, OQP_FOCK_B, OQP_VEC_MO_B, OQP_td_bvec_mo, OQP_td_t, &
      OQP_td_energies /)

    mol_mult = infos%mol_prop%mult
    if (mol_mult/=3) call show_message(&
            'MRSF-TDDFT are available for ROHF/UHF ref.&
            &with ONLY triplet multiplicity(mult=3)', with_abort)

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.

    dft = infos%control%hamilton == 20
    mrsf_zvector_breakdown = .false.

  ! Optional per-iteration profiler (env OQP_MRSF_ZV_TIMERS; default off)
    call zv_timers_begin()
  ! Optional performance opt-ins (env-gated; default off)
    call zv_read_config()

  ! Files open
  ! 3. LOG: Write: Main output file
    open (unit=iw, file=infos%log_filename, position="append")
  !
    call print_module_info('MRSF_TDHF_Z_Vector','Solving Z-Vector for MRSF-TDDFT')

  ! Readings

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    ! Build the z-vector's DFT grid -- optionally coarser than the SCF grid
    ! (env OQP_MRSF_ZV_COARSEGRID). The grid params are restored immediately; the
    ! built molGrid carries the coarse grid through all of the z-vector's XC calls.
    if (dft) then
      zv_gname_save = infos%dft%grid_pruned_name
      if (zv_coarse_on) then
        infos%dft%grid_pruned = .true.
        infos%dft%grid_pruned_name = ' '
        do ig = 1, min(len_trim(zv_grid_name), 16)
          infos%dft%grid_pruned_name(ig) = zv_grid_name(ig:ig)
        end do
        write(iw,'(" MRSF z-vector coarse response grid: ",a)') trim(zv_grid_name)
      end if
      call dft_initialize(infos, basis, molGrid)
      if (zv_coarse_on) infos%dft%grid_pruned_name = zv_gname_save
    end if

  ! Parameter it should be inputed later
    mrst = infos%tddft%mult
    cnvtol = infos%tddft%zvconv
  ! Opt-in: override the z-vector convergence tolerance (env OQP_MRSF_ZV_CONV).
  ! Default 1e-10 is on the SQUARED residual and is typically over-converged for
  ! gradients; right-sizing it cuts iterations. Default unset = input zvconv.
    if (zv_conv_user > 0.0_dp) cnvtol = zv_conv_user

    nocca = infos%mol_prop%nelec_A
    nvira = nbf-noccA
    noccb = infos%mol_prop%nelec_B
    nvirb = nbf-noccb
    nsocc = nocca-noccb
    lzdim = noccb*(nsocc+nvira)+nsocc*nvira

    if(mrst==1 .or. mrst==3) then
      xvec_dim = nocca*nvirb
      allocate(&
    ! for Z-vector
        fmrst1(1,7,nbf,nbf), &
        bvec_mo_d(xvec_dim,1), &
        hxa(nbf,nocca), &
        hxb(nbf,nbf), &
    ! for gradient
        tij(nocca,nocca), &
        tab(nvirb,nvirb), &
        stat=ok, &
        source=0.0_dp)
    else if(mrst==5) then
      xvec_dim = noccb*nvira
      allocate(&
    ! for Z-vector
        hxa(nbf,nbf), &
        hxb(nbf,noccb), &
  ! for gradient
        tij(noccb,noccb), &
        tab(nvira,nvira), &
        stat=ok, &
        source=0.0_dp)
    endif
    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    allocate(&
  ! for Z-vector
      xminv(lzdim), &
      rhs(lzdim), &
      lhs(lzdim), &
      xm(lzdim), &
      xk(lzdim), &
      pk(lzdim), &
      errv(lzdim), &
   ! For gradient
      pa(nbf,nbf,2), &
      ppija(nocca,nocca), &
      ppijb(noccb,noccb), &
   ! Allocate TDDFT variables
      fa(nbf,nbf), &           ! Temporary matrix for diagonalization
      fb(nbf,nbf), &           ! Temporary matrix for diagonalization
      ab1_mo_a(nocca,nvira), &
      ab1_mo_b(noccb,nvirb), &
!   For scratch
      wrk1(nbf,nbf), &
      wrk2(nbf,nbf), &
      wrk3(nbf,nbf), &
      stat=ok, &
      source=0.0_dp)

    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    call infos%dat%remove_records(tags_alloc)

    call infos%dat%reserve_data(OQP_WAO, TA_TYPE_REAL64, nbf_tri, comment=OQP_WAO_comment)
    call infos%dat%reserve_data(OQP_td_mrsf_density, TA_TYPE_REAL64, nbf*nbf*7, (/7, nbf, nbf /), comment=OQP_td_mrsf_density)
    call infos%dat%reserve_data(OQP_td_p, TA_TYPE_REAL64, nbf_tri*2, (/ nbf_tri, 2 /), comment=OQP_td_p)
    call infos%dat%reserve_data(OQP_td_abxc, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_td_abxc)

    call data_has_tags(infos%dat, tags_alloc, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
    call tagarray_get_data(infos%dat, OQP_td_mrsf_density, td_mrsf_den)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_td_t, td_t)
    call tagarray_get_data(infos%dat, OQP_td_energies, mrsf_energies)

    ta => td_t(:,1)
    tb => td_t(:,2)

    target_state = min(infos%tddft%target_state, infos%tddft%nstate)
    if (target_state /=infos%tddft%target_state) then
      write(*,'(/1x,66("-")&
               &/1x,"WARNING: Target state has been changed to the max available nstates"/&
               &/1x,66("-")/)')
    end if

    ! Determine solver name for output (0=CG, 1=GMRES legacy, 2=MINRES, 3=AUTO)
    select case (infos%tddft%z_solver)
    case (3)
      solver_name = "AUTO"
    case (2)
      solver_name = "MINRES"
    case (1)
      solver_name = "GMRES"
    case default
      solver_name = "CG"
    end select

    ! Save unrelaxed density matrices and the `b=A*x` vector for target state
    if (mrst==1 .or. mrst==3 ) then
      call mrsfxvec(infos, bvec_mo(:,target_state), bvec_mo_d(:,1))
      call sfdmat(bvec_mo_d(:,1), td_abxc, mo_a, ta, tb, nocca, noccb)
    else if (mrst==5 ) then
      call sfdmat(bvec_mo(:,target_state), td_abxc, mo_a, tb, ta, noccb, nocca)
    end if

    bvec(1:nbf,1:nbf,1:1) => td_abxc

  ! Initialize ERI calculations
  ! Opt-in: loosen the 2e integral cutoff for the WHOLE z-vector build (RHS,
  ! per-iteration sigma, and relaxed-density tail) via env OQP_MRSF_ZV_CUTOFF.
  ! Gradients are more cutoff-sensitive than energies, so the safe value is
  ! tighter than the response's; find it empirically. Restored before return so
  ! the next step's SCF/response is unaffected. Default unset = exact (unchanged).
    zv_rc_save = infos%control%int2e_cutoff
    ! Progressive screening (below) keeps init at the TIGHT cutoff so the pair
    ! list is a full superset; it ramps the run-time threshold per iteration.
    ! The static loosen and progressive are therefore mutually exclusive.
    zv_rc_loosen = zv_cutoff_user > 0.0_dp .and. .not. zv_prog_on
    if (zv_rc_loosen) infos%control%int2e_cutoff = max(zv_rc_save, zv_cutoff_user)
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    if (zv_prog_on) write(iw,'(" MRSF z-vector progressive screening ON: ", &
        &"tau=clamp(",1p,e8.1,"*||r||, ",e8.1," , ",e8.1,"), pin@||r||^2<",e8.1)') &
        zv_prog_k, zv_rc_save, zv_prog_cap, zv_prog_pin

    write(*,'(/1x,71("-")&
             &/19x,"MRSF-DFT ENERGY GRADIENT CALCULATION"&
             &/1x,71("-")/)')

    write(iw,fmt='(5x,a/&
                  &5x,16("-")/&
                  &5x,a,x,i0,x,f17.10,x,"Hartree"/&
                  &5x,a,x,i0/&
                  &5x,a,x,e10.4/&
                  &5x,a,x,i0/&
                  &5x,a,x,a)') &
        'Z-vector options' &
      , 'Target state       is', target_state, infos%mol_energy%energy+mrsf_energies(target_state) &
      , 'Multiplicity       is', infos%tddft%mult &
      , 'Convergence        is', infos%tddft%zvconv &
      , 'Maximum iterations is', infos%control%maxit_zv &
      , 'Solver method      is', trim(solver_name)
    call flush(iw)

    ! ======================================================================
    ! Step 1: assemble the z-vector right-hand side and Fock/density pieces.
    ! ======================================================================
    if (zv_tmr_on) t0_zv = zv_wtime()
    call build_mrsf_zvector_rhs()
    if (zv_tmr_on) zv_t_rhs = zv_t_rhs + (zv_wtime() - t0_zv)

    write(*,'(/3x,25("-")&
             &/6x,"START Z-VECTOR LOOP (",A,")"&
             &/3x,25("-")/)') trim(solver_name)
    call flush(iw)

    call sfromcal(xm, xminv, mo_energy_a, fa, fb, nocca, noccb)
    call sanitize_zvector_preconditioner(xm, xminv, iw, MRSF_ZVEC_DENOMINATOR_FLOOR, "MRSF")

    ! ======================================================================
    ! Step 2: solve the z-vector linear system.
    !   0 = CG (default)   1 = GMRES (legacy)   2 = MINRES   3 = AUTO (CG->MINRES->GMRES)
    ! ======================================================================
    select case (infos%tddft%z_solver)
    case (2)
      call run_mrsf_minres_zvector()
    case (1)
      call run_mrsf_gmres_zvector()
    case (3)
      call run_mrsf_zvector_auto()
    case default
      call run_mrsf_cg_zvector()
    end select

    ! Progressive screening ramps the cutoff during the CG loop; restore the
    ! tight floor so the relaxed-density/W back-projection (which sets the
    ! gradient) is computed exactly.
    if (zv_prog_on) call int2_driver%set_cutoff(zv_rc_save)

! -----------------------------------------------
    if (mrsf_zvector_breakdown) then
       infos%mol_energy%Z_Vector_converged=.false.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector breakdown"&
             &/3x,24("-")/)')
       call flush(iw)
       call int2_data%clean()
       call int2_driver%clean()
       if (zv_rc_loosen) infos%control%int2e_cutoff = zv_rc_save
       if (dft) call dftclean(infos)
       call measure_time(print_total=1, log_unit=iw)
       call cleanup_gmres_work()
       close(iw)
       return
    end if

    if (error>cnvtol) then
       infos%mol_energy%Z_Vector_converged=.false.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector not converged"&
             &/3x,24("-")/)')
       write(iw,'(" MRSF z-vector solver reached the maximum iterations; solver = ",A)') trim(solver_name)
       write(iw,'(" final residual = ",1p,e13.6)') error
    else
       infos%mol_energy%Z_Vector_converged=.true.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector converged"&
             &/3x,24("-")/)')
       ! Cache the converged solution for warm-starting the next step (no-op
       ! unless OQP_MRSF_ZV_WARMSTART is set).
       call zv_store_guess(xk, lzdim, target_state, int(infos%tddft%nstate), &
                           mo_a, mo_b, nbf)
    endif

    call flush(iw)

    ! ======================================================================
    ! Step 3: build the relaxed density (td_p) and energy-weighted density (wao).
    ! ======================================================================
    if (zv_tmr_on) t0_zv = zv_wtime()
    call build_mrsf_relaxed_density_and_w()
    if (zv_tmr_on) zv_t_back = zv_t_back + (zv_wtime() - t0_zv)

    call zv_timers_report(iw, solver_name)

    call int2_driver%clean()
    if (zv_rc_loosen) infos%control%int2e_cutoff = zv_rc_save

    if (dft) call dftclean(infos)

    call measure_time(print_total=1, log_unit=iw)
    ! Clean up GMRES work arrays
    call cleanup_gmres_work()

    close(iw)

  contains

    ! GMRES z-vector solve.  Reports breakdown via mrsf_zvector_breakdown so the
    ! auto driver can detect failure uniformly across solvers.
    subroutine run_mrsf_gmres_zvector()
      call zv_warm_seed()
      call gmres_solve( &
          apply_operator = apply_z_operator, &
          apply_precond = lambda_precond, &
          b = rhs, &
          x = xk, &
          n = lzdim, &
          restart = min(int(infos%tddft%gmres_dim), lzdim), &
          max_iter = int(infos%control%maxit_zv), &
          tol = cnvtol, &
          infos = infos, basis = basis, molGrid = molGrid, &
          int2_driver = int2_driver, &
          nocca = nocca, noccb = noccb, nbf = nbf, &
          mo_a = mo_a, mo_b = mo_b, mo_energy_a = mo_energy_a, &
          fa = fa, fb = fb, scale_exch = scale_exch, dft = dft, &
          error_out = error, iter_out = gmres_iter, iw = iw)
      if (.not. ieee_is_finite(error) .or. error > cnvtol) then
        if (.not. ieee_is_finite(error)) mrsf_zvector_breakdown = .true.
      end if
      write(iw,'(/," Final Summary:")')
      write(iw,'(" GMRES total iterations: ", I4)') gmres_iter
      write(iw,'(" Final error norm      : ", 1p,e13.6)') error
      write(iw,'(" Convergence criterion : ", 1p,e13.6)') cnvtol
      call flush(iw)
      call cleanup_gmres_work()
    end subroutine run_mrsf_gmres_zvector

    ! MINRES z-vector solve: symmetric short-recurrence solver (Paige-Saunders)
    ! that stays stable when (A+B) turns indefinite, at CG-like cost.  Uses the
    ! same apply_z_operator / apply_z_precond as CG and GMRES.
    subroutine run_mrsf_minres_zvector()
      ! NOTE: this MINRES wrapper solves A x = rhs starting from the zero vector
      ! (mr%init seeds from rhs), so warm-start has no effect here; seeding is a
      ! no-op kept for uniformity. Warm-start benefits CG (default) and GMRES.
      call zv_warm_seed()
      call mr%init(b=rhs, update=minres_apply_op, precond=minres_apply_pc, &
                   dat=minres_dummy, tol=cnvtol)
      minres_iter = 0
      if (mr%errcode == MINRES_OK) then
        do iter = 1, infos%control%maxit_zv
          call mr%step()
          minres_iter = iter
          if (mr%errcode /= MINRES_OK) exit
        end do
      end if
      if (mr%errcode == MINRES_CONVERGED .or. mr%errcode == MINRES_OK) then
        xk = mr%x
        error = mr%error
      else
        write(iw,'(" MRSF MINRES Z-Vector breakdown (errcode=",I0,")")') int(mr%errcode)
        mrsf_zvector_breakdown = .true.
        error = huge(1.0_dp)
      end if
      call mr%clean()
      write(iw,'(/," Final Summary:")')
      write(iw,'(" MINRES total iterations: ", I4)') minres_iter
      write(iw,'(" Final error norm       : ", 1p,e13.6)') error
      write(iw,'(" Convergence criterion  : ", 1p,e13.6)') cnvtol
      call flush(iw)
    end subroutine run_mrsf_minres_zvector

    ! AUTO driver: try CG (cheap, needs SPD), fall back to MINRES (CG-cost,
    ! robust on indefinite operators), then GMRES (general, priciest).  rhs and
    ! the preconditioner xminv are solver-independent and reused; only xk and
    ! the breakdown flag reset between attempts.  solver_name is updated to the
    ! solver that actually converged so the summary reflects reality.
    subroutine run_mrsf_zvector_auto()
      call run_mrsf_cg_zvector()
      if (.not. mrsf_zvector_breakdown .and. ieee_is_finite(error) &
          .and. error <= cnvtol) then
        solver_name = "AUTO(CG)"
        return
      end if
      write(iw,'(/," [AUTO] CG did not converge (breakdown/maxit); ", &
                 &"falling back to MINRES")')
      call flush(iw)
      mrsf_zvector_breakdown = .false.
      call run_mrsf_minres_zvector()
      if (.not. mrsf_zvector_breakdown .and. ieee_is_finite(error) &
          .and. error <= cnvtol) then
        solver_name = "AUTO(MINRES)"
        return
      end if
      write(iw,'(/," [AUTO] MINRES did not converge; falling back to GMRES")')
      call flush(iw)
      mrsf_zvector_breakdown = .false.
      call run_mrsf_gmres_zvector()
      if (.not. mrsf_zvector_breakdown .and. ieee_is_finite(error) &
          .and. error <= cnvtol) then
        solver_name = "AUTO(GMRES)"
      else
        solver_name = "AUTO(failed)"
      end if
    end subroutine run_mrsf_zvector_auto

    ! Preconditioned CG z-vector solve (default path).  All state is
    ! reached by host association, matching the inline version exactly.
    subroutine run_mrsf_cg_zvector()
      real(kind=dp) :: t0, rhs_norm2
      logical :: warm_used

      ! ============================================
      ! ORIGINAL CONJUGATE GRADIENT SOLVER
      ! ============================================

      ! Initial guess: zero (cold) or the previous step's solution (warm-start).
      ! With a warm xk the initial residual machinery below stays correct -- it
      ! builds A*xk and r0 = rhs - A*xk exactly as for the cold start.
      call zv_warm_seed(warm_used)

      if (zv_tmr_on) t0 = zv_wtime()
      call sfrogen(wrk1, wrk2, xk, nocca, noccb)
      ! Alpha
      call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
      ! Beta
      call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)
      if (zv_tmr_on) zv_t_trans = zv_t_trans + (zv_wtime() - t0)

      !****** INITIAL (A+B)*xk : residual r0 = rhs - A*xk **********************
      ! IMPORTANT: this MUST apply the SAME operator the CG loop iterates with
      ! (int2_tdgrd_data_t, int_amb=.false., beta channel = apb(:,:,2)), or the
      ! recurrence residual `errv` drifts from the true residual and CG converges
      ! to the wrong xk. For the cold start (xk=0) the density and hence lhs are
      ! zero regardless, so this is bit-identical to the previous code; it only
      ! matters once a warm-start guess (xk/=0) is used.
      call int2_data%clean()
      deallocate(int2_data)
      int2_data = int2_tdgrd_data_t( &
          d2 = pa, &
          int_apb = .true., &
          int_amb = .false., &
          tamm_dancoff = .false., &
          scale_exchange = scale_exch)

      if (zv_tmr_on) t0 = zv_wtime()
      call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)
      if (zv_tmr_on) zv_t_int2 = zv_t_int2 + (zv_wtime() - t0)
      ab1 => int2_data%apb(:,:,:,1)

      if (zv_tmr_on) t0 = zv_wtime()
      call symmetrize_matrix(pa(:,:,1), nbf)
      call symmetrize_matrix(pa(:,:,2), nbf)
      call utddft_fxc( &
          basis = basis, &
          molGrid = molGrid, &
          isVecs = .true., &
          wfa = mo_a, &
          wfb = mo_b, &
          fxa = ab1(:,:,1:1), &
          fxb = ab1(:,:,2:2), &
          dxa = pa(:,:,1:1), &
          dxb = pa(:,:,2:2), &
          nmtx = 1, &
          threshold = 1.0d-15, &
          infos = infos)
      if (zv_tmr_on) zv_t_xc = zv_t_xc + (zv_wtime() - t0)

      if (zv_tmr_on) t0 = zv_wtime()
      !   ALPHA: AO(M,N) -> MO(IA+) ... LPTMOA
      call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

      call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

      call sfrolhs(lhs, xk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, &
                   nocca, noccb)
      if (zv_tmr_on) zv_t_trans = zv_t_trans + (zv_wtime() - t0)

      call pcgrbpini(errv, pk, error, rhs, xminv, lhs)
      ! Warm-start safeguard: if the seeded guess did not reduce the residual
      ! below the cold-start value ||rhs||^2 (e.g. large MO rotation between
      ! steps, as in degenerate systems), discard it and restart from zero. The
      ! cold residual is exact (A*0 = 0 => r0 = rhs), so the fallback needs no
      ! extra Fock build -- it just reuses rhs.
      if (warm_used) then
        rhs_norm2 = dot_product(rhs, rhs)
        if (.not. ieee_is_finite(error) .or. error >= rhs_norm2) then
          write(iw,'(" MRSF z-vector warm-start: guess rejected (r0^2 ",1p,e10.3, &
                   &" >= cold ",e10.3,"); restarting from zero")') error, rhs_norm2
          call flush(iw)
          xk = 0.0_dp
          lhs = 0.0_dp
          call pcgrbpini(errv, pk, error, rhs, xminv, lhs)
        end if
      end if
      if (.not. ieee_is_finite(error) .or. any(.not. ieee_is_finite(errv)) .or. &
          any(.not. ieee_is_finite(pk)) .or. any(.not. ieee_is_finite(lhs))) then
        write(iw,'(" MRSF CG Z-Vector breakdown: non-finite initial PCG state")')
        mrsf_zvector_breakdown = .true.
        error = huge(1.0_dp)
      end if

      write(iw,'(" Initial error =",3x,1p,e10.3,1x,"/",1p,e10.3)') error, cnvtol
      call flush(iw)

      ! -----------------------------------------------

      do iter = 1, infos%control%maxit_zv
        zv_n_iter = iter

        if (zv_tmr_on) t0 = zv_wtime()
        call sfrogen(wrk1, wrk2, pk, nocca, noccb)
        !     Alpha
        call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
        !     Beta
        call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)
        if (zv_tmr_on) zv_t_trans = zv_t_trans + (zv_wtime() - t0)

        ! Progressive screening: loosen the cutoff for this sigma build based on
        ! the residual carried in from the previous step (the initial residual on
        ! iter 1). Pinned to the tight floor once the residual is small, so the
        ! converged solution is exact.
        if (zv_prog_on) call int2_driver%set_cutoff(zv_prog_tau(error, zv_rc_save))

        !     (A+B)*PK
        call int2_data%clean()
        deallocate(int2_data)
        int2_data = int2_tdgrd_data_t( &
            d2 = pa, &
            int_apb = .true., &
            int_amb = .false., &
            tamm_dancoff = .false., &
            scale_exchange = scale_exch)

        if (zv_tmr_on) t0 = zv_wtime()
        call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)
        if (zv_tmr_on) zv_t_int2 = zv_t_int2 + (zv_wtime() - t0)
        ab1 => int2_data%apb(:,:,:,1)

        !ab1 = ab1/2
        call symmetrize_matrix(pa(:,:,1), nbf)
        call symmetrize_matrix(pa(:,:,2), nbf)
        if (zv_tmr_on) t0 = zv_wtime()
        call utddft_fxc( &
            basis = basis, &
            molGrid = molGrid, &
            isVecs = .true., &
            wfa = mo_a, &
            wfb = mo_b, &
            fxa = ab1(:,:,1:1), &
            fxb = ab1(:,:,2:2), &
            dxa = pa(:,:,1:1), &
            dxb = pa(:,:,2:2), &
            nmtx = 1, &
            threshold = 1.0d-15, &
            infos = infos)
        if (zv_tmr_on) zv_t_xc = zv_t_xc + (zv_wtime() - t0)

        if (zv_tmr_on) t0 = zv_wtime()
        !     ALPHA: AO(M,N) -> MO(IA+) ... LPTMOA
        call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

        call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

        call sfrolhs(lhs, pk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, &
                     nocca, noccb)
        if (zv_tmr_on) zv_t_trans = zv_t_trans + (zv_wtime() - t0)

        if (any(.not. ieee_is_finite(lhs)) .or. any(.not. ieee_is_finite(pk))) then
          write(iw,'(" MRSF CG Z-Vector breakdown: non-finite lhs/search direction")')
          mrsf_zvector_breakdown = .true.
          error = huge(1.0_dp)
          exit
        end if

        if (zv_tmr_on) t0 = zv_wtime()
        pap = dot_product(pk, lhs)
        if (.not. ieee_is_finite(pap) .or. abs(pap) < MRSF_ZVEC_DENOMINATOR_FLOOR) then
          write(iw,'(" MRSF CG Z-Vector breakdown: unsafe p^T A p denominator")')
          mrsf_zvector_breakdown = .true.
          error = huge(1.0_dp)
          exit
        end if

        alpha = 1.0_dp / pap
        if (.not. ieee_is_finite(alpha)) then
          write(iw,'(" MRSF CG Z-Vector breakdown: non-finite alpha")')
          mrsf_zvector_breakdown = .true.
          error = huge(1.0_dp)
          exit
        end if

        xk = xk + pk * alpha
        errv = errv - alpha*lhs
        if (any(.not. ieee_is_finite(xk)) .or. any(.not. ieee_is_finite(errv))) then
          write(iw,'(" MRSF CG Z-Vector breakdown: non-finite solution/residual update")')
          mrsf_zvector_breakdown = .true.
          error = huge(1.0_dp)
          exit
        end if

        error = dot_product(errv, errv)
        if (.not. ieee_is_finite(error)) then
          write(iw,'(" MRSF CG Z-Vector breakdown: non-finite residual norm")')
          mrsf_zvector_breakdown = .true.
          error = huge(1.0_dp)
          exit
        end if
        write(iw,'(" Iter#",I2," Error =",&
              &3x,1p,e10.3,1x,"/",1p,e10.3)') &
                iter, error, cnvtol
        call flush(iw)

        if (error<cnvtol) then
          if (zv_tmr_on) zv_t_cg = zv_t_cg + (zv_wtime() - t0)
          exit
        end if

        call pcgb(pk, errv, xminv)
        if (any(.not. ieee_is_finite(pk))) then
          write(iw,'(" MRSF CG Z-Vector breakdown: non-finite search direction after preconditioner")')
          mrsf_zvector_breakdown = .true.
          error = huge(1.0_dp)
          exit
        end if
        if (zv_tmr_on) zv_t_cg = zv_t_cg + (zv_wtime() - t0)

      end do

      ! Always leave int2_driver at the tight cutoff on exit. Progressive
      ! screening ramps it during the loop; if CG exits UNCONVERGED (e.g. the
      ! AUTO solver then falls back to MINRES/GMRES, which reuse int2_driver via
      ! apply_z_operator and do not ramp), the fallback must see the tight
      ! operator -- otherwise it would solve/check convergence against the
      ! screened one. (On convergence the last step is already pinned tight.)
      if (zv_prog_on) call int2_driver%set_cutoff(zv_rc_save)
          end subroutine run_mrsf_cg_zvector

    ! Lambda wrapper for preconditioner
    subroutine lambda_precond(x_in, x_out)
      real(kind=dp), intent(in) :: x_in(:)
      real(kind=dp), intent(out) :: x_out(:)
      call apply_z_precond(x_in, x_out, xminv)
    end subroutine lambda_precond

    ! Seed xk for the solver: zero (cold), the raw previous solution, or — when
    ! the MOs have rotated between steps — the previous solution PROJECTED into
    ! the current MO basis via the geometry-stable AO z-density:
    !   xk_old -> av_old(MO,old) -> Pz = C_old av_old C_old^T (AO)
    !          -> C_new^T S Pz S C_new (MO,new) -> gather -> xk.
    ! Same-geometry round-trip is exact (C^T S C = I). The CG safeguard still
    ! rejects the guess if it doesn't beat the cold residual, so this can only
    ! help (large steps) or be ignored — never corrupt the gradient.
    subroutine zv_warm_seed(used)
      use mathlib, only: unpack_matrix
      logical, intent(out), optional :: used
      real(kind=dp), contiguous, pointer :: smptr(:)
      real(kind=dp), allocatable :: smat(:,:), ava(:,:), avb(:,:), pz(:,:), &
                                    tmp(:,:), avn1(:,:), avn2(:,:)
      integer :: st, ok
      logical :: have, got_warm
      if (present(used)) used = .false.
      xk = 0.0_dp
      got_warm = .false.
      st = target_state

      ! ---- try a warm-start guess from the previous step (short-circuit guards) ----
      have = zv_warm_on .and. allocated(zv_warm) .and. allocated(zv_warm_has) &
             .and. zv_warm_lzdim == lzdim
      if (have) have = (st >= 1 .and. st <= size(zv_warm_has))
      if (have) have = zv_warm_has(st)
      if (have) have = .not. any(.not. ieee_is_finite(zv_warm(:,st)))
      if (have) then
        if (zv_warm_have_mo .and. zv_warm_nbf == nbf) then
          ! MO-projected seed: xk_old -> Pz(AO,old) -> C_new^T S Pz S C_new -> gather
          allocate(smat(nbf,nbf), ava(nbf,nbf), avb(nbf,nbf), pz(nbf,nbf), &
                   tmp(nbf,nbf), avn1(nbf,nbf), avn2(nbf,nbf), stat=ok)
          if (ok == 0) then
            call tagarray_get_data(infos%dat, OQP_SM, smptr)
            call unpack_matrix(smptr, smat, nbf, 'U')
            call sfrogen(ava, avb, zv_warm(:,st), nocca, noccb)
            call orthogonal_transform('t', nbf, zv_warm_moa, ava, pz, wrk3)
            call dgemm('n','n',nbf,nbf,nbf,1.0_dp,smat,nbf,pz,nbf,0.0_dp,tmp,nbf)
            call dgemm('n','n',nbf,nbf,nbf,1.0_dp,tmp,nbf,smat,nbf,0.0_dp,pz,nbf)
            call orthogonal_transform('n', nbf, mo_a, pz, avn1, wrk3)
            call orthogonal_transform('t', nbf, zv_warm_mob, avb, pz, wrk3)
            call dgemm('n','n',nbf,nbf,nbf,1.0_dp,smat,nbf,pz,nbf,0.0_dp,tmp,nbf)
            call dgemm('n','n',nbf,nbf,nbf,1.0_dp,tmp,nbf,smat,nbf,0.0_dp,pz,nbf)
            call orthogonal_transform('n', nbf, mo_b, pz, avn2, wrk3)
            call zv_sfrogen_gather(avn1, avn2, xk, nocca, noccb)
            if (any(.not. ieee_is_finite(xk))) xk = zv_warm(:,st)
            deallocate(smat, ava, avb, pz, tmp, avn1, avn2)
            got_warm = .true.
            write(iw,'(" MRSF z-vector warm-start: MO-projected seed (state ",i0,")")') st
          end if
        else
          xk = zv_warm(:,st); got_warm = .true.
          write(iw,'(" MRSF z-vector warm-start: raw seed (state ",i0,")")') st
        end if
      end if

      if (got_warm) then
        if (present(used)) used = .true.
        call flush(iw)
        return
      end if

      ! ---- cold start: Jacobi (diagonal) guess x0 = M^-1 rhs, else zero ----
      ! Free (the cold initial Fock build runs regardless) and ~1 iteration ahead.
      ! Marked "used" so the CG safeguard validates it and falls back to zero if
      ! it does not reduce the residual.
      if (zv_diag_guess) then
        xk = xminv * rhs
        if (any(.not. ieee_is_finite(xk))) then
          xk = 0.0_dp
        else
          if (present(used)) used = .true.
          write(iw,'(" MRSF z-vector: Jacobi (diagonal) initial guess")')
          call flush(iw)
        end if
      end if
    end subroutine zv_warm_seed

    ! minres_matvec(y, x, dat) wrappers: y is the output, x the input, dat is
    ! unused (the solver context is reached by host association).
    subroutine minres_apply_op(y, x, dat)
      use iso_c_binding, only: c_ptr
      real(kind=dp) :: x(:)
      real(kind=dp) :: y(:)
      type(c_ptr) :: dat
      call apply_z_operator(x, y, infos, basis, molGrid, int2_driver, &
                            nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                            fa, fb, scale_exch, dft)
    end subroutine minres_apply_op

    subroutine minres_apply_pc(y, x, dat)
      use iso_c_binding, only: c_ptr
      real(kind=dp) :: x(:)
      real(kind=dp) :: y(:)
      type(c_ptr) :: dat
      call apply_z_precond(x, y, xminv)
    end subroutine minres_apply_pc


    ! Build the relaxed (z-vector) density (-> td_p) and energy-weighted
    ! density W (-> wao) from the converged z-vector xk.  Host association
    ! preserves behavior versus the previous inline tail.
    subroutine build_mrsf_relaxed_density_and_w()
      if (mrst==1 .or. mrst==3) then

        call sfropcal(wrk1, wrk2, tij, tab, xk, nocca, noccb)

      else if (mrst==5) then

        call mrsfqropcal(wrk1, wrk2, tab, tij, xk, nocca, noccb)

      end if

   !  Update density for alpha
      call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)

   !  Update density for beta
      call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)
      call int2_data%clean()
      deallocate(int2_data)
      int2_data = int2_tdgrd_data_t( &
          d2 = pa, &
          int_apb = .true., &
          int_amb = .false., &
          tamm_dancoff = .false., &
          scale_exchange = scale_exch)

      call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)
      ab1 => int2_data%apb(:,:,:,1)

      call symmetrize_matrix(pa(:,:,1), nbf)
      call symmetrize_matrix(pa(:,:,2), nbf)
      call pack_matrix(pa(:,:,1), td_p(:,1))
      call pack_matrix(pa(:,:,2), td_p(:,2))

      td_p = 0.5_dp*td_p

      call utddft_fxc( &
          basis = basis, &
          molGrid = molGrid, &
          isVecs = .true., &
          wfa = mo_a, &
          wfb = mo_b, &
          fxa = ab1(:,:,1:1), &
          fxb = ab1(:,:,2:2), &
          dxa = pa(:,:,1:1), &
          dxb = pa(:,:,2:2), &
          nmtx = 1, &
          threshold = 1.0d-15, &
          infos = infos)

  !   ALPHA AO(M,N) -> MO(I-,J-) ... LPPIJA
      call dgemm('n', 'n', nbf, nocca, nbf,  &
                 1.0_dp, ab1(:,:,1), nbf,  &
                         mo_a, nbf,  &
                 0.0_dp, wrk2, nbf)
      call dgemm('t', 'n', nocca, nocca, nbf,  &
                 1.0_dp, mo_a,  nbf,  &
                         wrk2,  nbf,  &
                 0.0_dp, ppija, nocca)
  !   BETA: AO(M,N) -> MO(I-,J-) ... LPPIJB
      call dgemm('n', 'n', nbf, noccb, nbf,  &
                 1.0_dp, ab1(:,:,2), nbf,  &
                         mo_b, nbf,  &
                 0.0_dp, wrk2, nbf)
      call dgemm('t', 'n', noccb, noccb, nbf,  &
                 1.0_dp, mo_b,  nbf,  &
                         wrk2,  nbf,  &
                 0.0_dp, ppijb, noccb)

  !   Calculate W (in MO basis)
      wmo => wrk3
      wmo = 0
      if (mrst==1 .or. mrst==3) then

        call mrsfrowcal(wmo, mo_energy_a, fa, fb, xk, &
                        hxa, hxb, ppija, ppijb, &
                        nocca, noccb)

      else if (mrst==5) then

        call mrsfqrowcal(wmo, mo_energy_a, fa, fb, xk, &
                         hxa, hxb, ppija, ppijb, &
                         nocca, noccb)

      end if

      call orthogonal_transform('t', nbf, mo_a, wmo, wrk2, wrk1)
      call symmetrize_matrix(wrk2, nbf)
      call pack_matrix(wrk2, wao)
      wao = wao*0.5_dp
  !   ROHF, half one more time:
      wao = wao*0.5_dp
    end subroutine build_mrsf_relaxed_density_and_w


    ! Assemble the z-vector right-hand side (rhs) and the ROHF Fock/density
    ! intermediates it needs.  Host association preserves behavior versus the
    ! previous inline RHS construction.
    subroutine build_mrsf_zvector_rhs()
    ! Prepare for ROHF
      ! Fock matrices A and B
      if( roref )then
          wrk1t(1:nbf*nbf) => wrk1
    !   Alapha
        call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, wrk1)
        call unpack_matrix(wrk1t, fa)

    !   Beta
        call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, wrk1)
        call unpack_matrix(wrk1t, fb)
      end if

    ! Make density like part
      call unpack_matrix(ta, pa(:,:,1))
      call unpack_matrix(tb, pa(:,:,2))

    ! Initialize ERI calculations
      scale_exch = 1.0_dp
      scale_exch2 = 1.0_dp
      if (dft) then
         scale_exch = infos%dft%HFscale    !> Reference HF exchange
         scale_exch2 = infos%tddft%HFscale !> Response HF exchange
      end if

      if (mrst==1 .or. mrst==3 ) then

        int2_data_st = int2_mrsf_data_t( &
            d3 = fmrst1, &
            tamm_dancoff = .true., &
            scale_exchange = scale_exch2, &
            scale_coulomb = scale_exch2)

      else if( mrst==5  )then

        int2_data_q = int2_td_data_t( &
            d2=bvec, &
            int_apb = .false., &
            int_amb = .false., &
            tamm_dancoff = .true., &
            scale_exchange = scale_exch2)

      end if

      int2_data = int2_tdgrd_data_t( &
          d2 = pa, &
          int_apb = .true., &
          int_amb = .false., &
          tamm_dancoff = .false., &
          scale_exchange = scale_exch)

      call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)
      ab1 => int2_data%apb(:,:,:,1)

      pa = pa*2
      call utddft_fxc( &
          basis = basis, &
          molGrid = molGrid, &
          isVecs = .true., &
          wfa = mo_a, &
          wfb = mo_b, &
          fxa = ab1(:,:,1:1), &
          fxb = ab1(:,:,2:2), &
          dxa = pa(:,:,1:1), &
          dxb = pa(:,:,2:2), &
          nmtx = 1, &
          threshold = 1.0d-15, &
          infos = infos)

  !   ALPHA: AO(M,N) -> MO(IA+)
      call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

      call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

      if (mrst==1 .or. mrst==3) then

        call iatogen(bvec_mo(:,target_state), wrk1, nocca, noccb)
        call mrsfcbc(infos, mo_a, mo_a, wrk1, fmrst1(1,:,:,:))

        fmrst1(1,7,:,:) = td_abxc

        td_mrsf_den(1:7,:,:) = fmrst1(1,1:7,:,:)

      ! Initialize ERI calculations
        call int2_driver%run(int2_data_st, &
              cam = dft.and.infos%dft%cam_flag, &
              alpha = infos%tddft%cam_alpha, &
              alpha_coulomb = infos%tddft%cam_alpha, &
              beta = infos%tddft%cam_beta,&
              beta_coulomb = infos%tddft%cam_beta, &
              mu = infos%tddft%cam_mu)
        fmrst2 => int2_data_st%f3(:,:,:,:,1)! ado2v, ado1v, adco1, adco2, ao21v, aco12, agdlr

      ! Scaling factor if triplet
        if (mrst==3) fmrst2(:,1:6,:,:) = -1.0_dp*fmrst2(:,1:6,:,:)

        ! Spin pair coupling
        if (infos%tddft%spc_coco /= infos%tddft%hfscale) &
           fmrst2(:,6,:,:) = fmrst2(:,6,:,:) * infos%tddft%spc_coco / infos%tddft%hfscale
        if (infos%tddft%spc_ovov /= infos%tddft%hfscale) &
           fmrst2(:,5,:,:) = fmrst2(:,5,:,:) * infos%tddft%spc_ovov / infos%tddft%hfscale
        if (infos%tddft%spc_coov /= infos%tddft%hfscale) &
           fmrst2(:,1:4,:,:) = fmrst2(:,1:4,:,:) * infos%tddft%spc_coov / infos%tddft%hfscale

        call orthogonal_transform('n', nbf, mo_a, fmrst2(1,7,:,:), wrk2, wrk1)

        call mrsfxvec(infos, bvec_mo(:,target_state), bvec_mo_d(:,1))

        call iatogen(bvec_mo_d(:,1), wrk3, nocca, noccb)

        call dgemm('n', 't', nbf, nocca, nbf, &
                   2.0_dp, wrk2, nbf, &
                           wrk3, nbf, &
                   0.0_dp, hxa, nbf)
        call dgemm('t', 'n', nbf, nbf, nocca, &
                   2.0_dp, wrk2, nbf, &
                           wrk3, nbf, &
                   0.0_dp, hxb, nbf)

     ! spin pair ov-ov, co-co, co-ov coupling
        call mrsfsp(hxa, hxb, mo_a, mo_a, wrk3, fmrst2(1,:,:,:), nocca, noccb)

     !  Unrelaxed difference density matries T_ij and T_ab
     !  Ta(i+,j+):= -X(i+,a-)*X(j+,a-) for singlet and triplet
        call dgemm('n', 't', nocca, nocca, nvirb, &
                  -1.0_dp, bvec_mo_d, nocca, &
                           bvec_mo_d, nocca, &
                   0.0_dp, tij, nocca)

     !  Tb(a-,b-):= X(i+,a-)*X(i+,b-) for singlet and triplet
        call dgemm('t', 'n', nvirb, nvirb, nocca, &
                   1.0_dp, bvec_mo_d, nocca, &
                           bvec_mo_d, nocca, &
                   0.0_dp, tab, nvirb)

        call sfrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, &
                     Tij, Tab, Fa, Fb, nocca, noccb)

      else if(mrst==5) then

     !  Initialize ERI calculations
        call int2_driver%run(int2_data_q, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%tddft%cam_alpha, &
              beta=infos%tddft%cam_beta,&
              mu=infos%tddft%cam_mu)

        call orthogonal_transform('n', nbf, mo_a, int2_data_q%amb(:,:,1,1), wrk2, wrk1)

        call iatogen(bvec_mo(:,target_state),wrk3,noccb,nocca)

        call dgemm('t', 'n', nbf, nbf, noccb, &
                   2.0_dp, wrk2, nbf, &
                           wrk3, nbf, &
                   0.0_dp, hxa, nbf)
        call dgemm('n', 't', nbf, noccb, nbf, &
                   2.0_dp, wrk2, nbf, &
                           wrk3, nbf, &
                   0.0_dp, hxb, nbf)

     !  Unrelaxed difference density matries T_ij and T_ab
     !  Ta(i+,j+):= -X(i+,a-)*X(j+,a-) for singlet and triplet
        call dgemm('n', 't', noccb, noccb, nvira, &
                  -1.0_dp, bvec_mo(:,target_state), noccb, &
                           bvec_mo(:,target_state), noccb, &
                   0.0_dp, tij, noccb)

     !  Tb(a-,b-):= X(i+,a-)*X(i+,b-) for singlet and triplet
        call dgemm('t', 'n', nvira, nvira, noccb, &
                   1.0_dp, bvec_mo(:,target_state), noccb, &
                           bvec_mo(:,target_state), noccb, &
                   0.0_dp, tab, nvira)

        call mrsfqrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, &
                        tab, tij, fa, fb, nocca, noccb)
      end if
    end subroutine build_mrsf_zvector_rhs

  end subroutine tdhf_mrsf_z_vector

end module tdhf_mrsf_z_vector_mod
