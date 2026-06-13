#define UNUSED_DUMMY(x) if (.false.) then ; if (size(shape(x))<0) continue ; end if

module int2_compute
  use precision, only: dp

  use int2e_libint, only: libint_t, libint2_active
  use int2e_rys, only: int2_rys_data_t
  use basis_tools, only: basis_set
  use atomic_structure_m, only: atomic_structure
  use int2_pairs, only: &
    int2_cutoffs_t, &
    int2_pair_storage
  use messages, only: show_message, WITH_ABORT
  use parallel, only: par_env_t

  implicit none

!###############################################################################

  character(len=*), parameter :: module_name = "int2_compute"

!###############################################################################

  integer, parameter :: ERR_CAM_PARAM = 1

!###############################################################################

  private
  public int2_compute_t
  public int2_storage_t
  public int2_compute_data_t
  public int2_fock_data_t
  public int2_rhf_data_t
  public int2_urohf_data_t
  public ints_exchange
  public petite_quartet_weight
  public load_petite_shell_map

!###############################################################################

  type eri_data_t
    logical :: attenuated_ints = .false.
    logical :: rys_only = .false.
    integer :: ids(4)
    integer :: flips(4)
    integer :: am(4)
    integer :: nbf(4)
    real(kind=dp) :: mu2 = 1.0d99
    ! Petite-list orbit weight (q4); 1 unless symmetry reduction is active.
    real(kind=dp) :: weight = 1.0d0
    ! Weighted screening thresholds: needed by the non-abelian full-group
    ! tier (orbit members have unequal magnitudes); the abelian tier keeps
    ! the C1-identical unweighted cutoff for exact cancellation.
    logical :: weighted_cutoff = .false.
    real(kind=dp), pointer :: pints(:,:,:,:)
    type(libint_t), allocatable :: erieval(:)
    type(int2_rys_data_t), allocatable :: gdat
    real(kind=dp), allocatable :: ints(:)
  end type eri_data_t

!###############################################################################

  type :: int2_storage_t
    integer :: ncur = 0
    integer :: buf_size = 0
    integer :: thread_id = 1
    integer(2), allocatable :: ids(:,:)
    real(kind=dp), allocatable, dimension(:) :: ints
  contains
    procedure, pass :: init => int2_storage_init
    procedure, pass :: clean => int2_storage_clean
  end type

!###############################################################################

  type, abstract :: int2_compute_data_t
    logical :: multipass = .false.
    integer :: num_passes = 1
    integer :: cur_pass = 1
    real(kind=dp) :: scale_coulomb = 1.0d0
    real(kind=dp) :: scale_exchange = 1.0d0
    type(par_env_t) :: pe
  contains
!    procedure, pass :: storeints => int2_compute_data_t_storeints
    procedure(int2_compute_data_parallel_start), deferred, pass :: parallel_start
    procedure(int2_compute_data_parallel_stop), deferred, pass :: parallel_stop
    procedure(int2_compute_data_update), deferred, pass :: update
    procedure(int2_compute_data_clean), deferred, pass :: clean
    procedure, pass :: screen_ij => int2_compute_data_t_screen_ij
    procedure, pass :: screen_ijkl => int2_compute_data_t_screen_ijkl
    generic :: screen => screen_ij, screen_ijkl
  end type

!###############################################################################

  type, abstract, extends(int2_compute_data_t) :: int2_fock_data_t
    integer :: nshells = 0
    integer :: fockdim = 0
    integer :: nthreads = 1
    integer :: nfocks = 1
    real(kind=dp), allocatable :: f(:,:,:)
    real(kind=dp), allocatable :: dsh(:,:)
    real(kind=dp) :: max_den = 1.0d0
    real(kind=dp), pointer :: d(:,:) => null()
    !> Memory mode for the per-thread Fock accumulator.  Default (.false.) keeps
    !> one full Fock copy PER THREAD (fast, no contention) and reduces at the end
    !> -- O(fockdim*nthreads) memory.  When .true. a SINGLE shared Fock is used
    !> with atomic updates -- O(fockdim) memory, ~Nthreads x smaller, for very
    !> large systems.  Auto-enabled when the replicated buffers would exceed
    !> OQP_FOCK_MEM_MB (default 4096); forced by OQP_FOCK_ATOMIC=1/0.
    logical :: atomic_fock = .false.
  contains
    procedure :: parallel_start => int2_fock_data_t_parallel_start
    procedure :: parallel_stop => int2_fock_data_t_parallel_stop
    procedure :: clean => int2_fock_data_t_clean
    procedure :: init_screen => int2_fock_data_t_init_screen
    procedure :: screen_ij => int2_fock_data_t_screen_ij
    procedure :: screen_ijkl => int2_fock_data_t_screen_ijkl
    procedure :: int2_fock_data_t_parallel_start
    procedure :: int2_fock_data_t_parallel_stop
  end type

  type, extends(int2_fock_data_t) :: int2_rhf_data_t
  contains
    procedure :: parallel_start => int2_rhf_data_t_parallel_start
    procedure :: update => int2_rhf_data_t_update
  end type

  type, extends(int2_fock_data_t) :: int2_urohf_data_t
  contains
    procedure :: parallel_start => int2_urohf_data_t_parallel_start
    procedure :: update => int2_urohf_data_t_update
  end type

!###############################################################################

  type :: int2_compute_t

    type(basis_set), pointer :: basis
    type(atomic_structure), pointer :: atoms
    real(kind=dp), allocatable :: schwarz_ints_regular(:,:)
    real(kind=dp), allocatable :: schwarz_ints_attenuated(:,:)
    real(kind=dp), contiguous, pointer :: schwarz_ints(:,:) => null()

    logical :: schwarz = .true.
    integer :: buf_size = 50000

    integer :: skipped = 0

    type(int2_cutoffs_t) :: cutoffs
    type(int2_pair_storage) :: ppairs

    logical :: attenuated = .false.
    !> Force the native Rys path for all L>2 quartets, bypassing libint even
    !> when it is compiled in.  Consumers whose validated reference data was
    !> produced with the Rys kernels (e.g. the NMR magnetic response) set this.
    logical :: rys_only = .false.
    real(kind=dp) :: mu = 1.0d99

    ! Symmetry petite list (loaded from tagarray when pyoqp enables
    ! use_integral_symmetry). The shell map is 1-based, stored flat with
    ! shell index fastest: map(shell, op) = sym_shell_map((op-1)*nshell+shell).
    logical :: petite = .false.
    logical :: sym_full = .false.
    integer :: sym_nops = 0
    integer(8), contiguous, pointer :: sym_shell_map(:) => null()

    type(par_env_t) :: pe
  contains

    private
!    procedure, pass :: storeints => int2_compute_data_t_storeints
    procedure, public, pass :: init => int2_compute_t_init
    procedure, public, pass :: enable_petite => int2_compute_t_enable_petite
    procedure, public, pass :: set_screening => int2_compute_t_set_screening
    procedure, public, pass :: clean => int2_compute_t_clean
    procedure, public, pass :: run => int2_run
    procedure, public, pass :: run_generic => int2_twoei
    procedure, public, pass :: run_cam => int2_run_cam
  end type

!###############################################################################

  abstract interface

    subroutine int2_process(this)
      import :: int2_compute_t
      implicit none
      class(int2_compute_t), intent(inout) :: this
    end subroutine

    subroutine int2_compute_data_parallel_start(this, basis, nthreads)
      import :: int2_compute_data_t, basis_set
      implicit none
      class(int2_compute_data_t), target, intent(inout) :: this
      type(basis_set), intent(in) :: basis
      integer, intent(in) :: nthreads
    end subroutine

    subroutine int2_compute_data_parallel_stop(this)
      import :: int2_compute_data_t
      implicit none
      class(int2_compute_data_t), intent(inout) :: this
    end subroutine

    subroutine int2_compute_data_clean(this)
      import :: int2_compute_data_t
      implicit none
      class(int2_compute_data_t), intent(inout) :: this
    end subroutine

    subroutine int2_compute_data_update(this, buf)
      import :: int2_compute_data_t, int2_storage_t, dp
      implicit none
      class(int2_compute_data_t), intent(inout) :: this
      type(int2_storage_t), intent(inout) :: buf
    end subroutine

  end interface

!###############################################################################

contains

!###############################################################################

  subroutine int2_compute_t_clean(this)
    implicit none
    class(int2_compute_t), intent(inout) :: this
    call this%ppairs%clean()
    this%basis => null()
    this%atoms => null()
    if (allocated(this%schwarz_ints_regular)) deallocate(this%schwarz_ints_regular)
    if (allocated(this%schwarz_ints_attenuated)) deallocate(this%schwarz_ints_attenuated)
    this%skipped = 0
  end subroutine

!###############################################################################

  subroutine int2_compute_t_init(this, basis, infos)
    use types, only: information
    use oqp_tagarray_driver

    use int2e_libint, only: &
      libint_static_init

    implicit none

    character(len=*), parameter :: subroutine_name = "int2_compute_t_init"

    class(int2_compute_t), intent(inout) :: this
    type(information), target, intent(inout) :: infos
    type(basis_set), target :: basis
    real(kind=dp) :: cutoff
    real(kind=dp), parameter :: ec1 = 1.0d-02, ec2 = 1.0d-04
    real(kind=dp), parameter :: cx1 = 25.0d+00
    integer(4) :: status

    this%basis => basis
    this%atoms => infos%atoms
    allocate(this%schwarz_ints_regular(basis%nshell,basis%nshell), source=0d0)
    this%skipped = 0

    call basis%init_shell_centers()

    cutoff = infos%control%int2e_cutoff
    call this%cutoffs%set(cutoff, ec1*cutoff, ec2*cutoff, cx1*log(10.0d0))

    call this%pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)
    if (libint2_active) call libint_static_init

    call this%ppairs%alloc(basis, this%cutoffs)
    call this%ppairs%compute(basis, this%cutoffs)

    ! Note: the petite-list reduction is NOT loaded here. It is only valid
    ! for totally symmetric densities (SCF Fock), so the SCF caller opts in
    ! explicitly via enable_petite(); response/CPHF builders must not.
    this%petite = .false.
    this%sym_nops = 0
    this%sym_shell_map => null()

  end subroutine int2_compute_t_init

!###############################################################################

!> @brief Opt into the symmetry petite-list reduction (skeleton Fock).
!> @detail Loads the shell map written by pyoqp when use_integral_symmetry
!>   is enabled. Only valid when the contracted density is totally
!>   symmetric and the caller symmetrizes the resulting skeleton matrix.
  subroutine int2_compute_t_enable_petite(this, infos)
    use types, only: information
    use oqp_tagarray_driver
    use tagarray, only: TA_OK

    implicit none

    class(int2_compute_t), intent(inout) :: this
    type(information), target, intent(inout) :: infos

    integer(8), contiguous, pointer :: petite_flag(:)
    integer(4) :: status

    call tagarray_get_data(infos%dat, OQP_sym_petite, petite_flag, status=status)
    if (status /= TA_OK) return
    if (petite_flag(1) == 0) return

    call tagarray_get_data(infos%dat, OQP_sym_shell_map, this%sym_shell_map, status=status)
    if (status /= TA_OK) then
      this%sym_shell_map => null()
      return
    end if
    if (mod(size(this%sym_shell_map), this%basis%nshell) /= 0) then
      ! Stale map from a different basis: ignore (fail-safe to C1).
      this%sym_shell_map => null()
      return
    end if

    this%sym_nops = int(size(this%sym_shell_map)/this%basis%nshell)
    this%petite = this%sym_nops > 1

    block
      real(kind=dp), contiguous, pointer :: blocks(:)
      call tagarray_get_data(infos%dat, OQP_sym_op_blocks, blocks, status=status)
      this%sym_full = status == TA_OK
    end block

  end subroutine int2_compute_t_enable_petite

!###############################################################################

!> @brief Load the petite-list shell map written by pyoqp, if enabled.
!> @detail nops = 0 when the reduction is disabled or the map is stale.
  subroutine load_petite_shell_map(infos, nshell, map, nops)
    use types, only: information
    use oqp_tagarray_driver
    use tagarray, only: TA_OK

    implicit none

    type(information), target, intent(inout) :: infos
    integer, intent(in) :: nshell
    integer(8), contiguous, pointer, intent(out) :: map(:)
    integer, intent(out) :: nops

    integer(8), contiguous, pointer :: petite_flag(:)
    integer(4) :: status

    map => null()
    nops = 0

    call tagarray_get_data(infos%dat, OQP_sym_petite, petite_flag, status=status)
    if (status /= TA_OK) return
    if (petite_flag(1) == 0) return

    call tagarray_get_data(infos%dat, OQP_sym_shell_map, map, status=status)
    if (status /= TA_OK) then
      map => null()
      return
    end if
    if (mod(size(map), nshell) /= 0) then
      map => null()
      return
    end if

    nops = int(size(map)/nshell)
    if (nops < 2) then
      map => null()
      nops = 0
    end if

  end subroutine load_petite_shell_map

!###############################################################################

!> @brief Petite-list weight of a canonical shell quartet.
!> @detail Returns 0 if (i,j,k,l) is not the lexicographically largest
!>   member of its orbit under the abelian group (the caller skips it),
!>   otherwise the orbit size |G|/|stabilizer|. The quartet must be in
!>   canonical order (i>=j, (i,j)>=(k,l), k>=l), which the int2_twoei
!>   loop guarantees.
  integer function petite_quartet_weight(map, nops, nshell, i, j, k, l) result(q4)
    implicit none
    integer(8), intent(in) :: map(:)
    integer, intent(in) :: nops, nshell, i, j, k, l

    integer :: iop, nstab, t, base
    integer :: pi, pj, pk, pl, a1, a2, b1, b2

    nstab = 0
    do iop = 1, nops
      base = (iop-1)*nshell
      pi = int(map(base + i))
      pj = int(map(base + j))
      pk = int(map(base + k))
      pl = int(map(base + l))

      a1 = max(pi, pj); a2 = min(pi, pj)
      b1 = max(pk, pl); b2 = min(pk, pl)
      if (a1 < b1 .or. (a1 == b1 .and. a2 < b2)) then
        t = a1; a1 = b1; b1 = t
        t = a2; a2 = b2; b2 = t
      end if

      if (a1 /= i) then
        if (a1 > i) then
          q4 = 0
          return
        end if
      else if (a2 /= j) then
        if (a2 > j) then
          q4 = 0
          return
        end if
      else if (b1 /= k) then
        if (b1 > k) then
          q4 = 0
          return
        end if
      else if (b2 /= l) then
        if (b2 > l) then
          q4 = 0
          return
        end if
      else
        nstab = nstab + 1
      end if
    end do

    q4 = nops / nstab

  end function petite_quartet_weight

!###############################################################################

  subroutine int2_compute_t_set_screening(this)
    implicit none
    class(int2_compute_t), intent(inout) :: this
    call ints_exchange(this%basis, this%schwarz_ints_regular, rys_only=this%rys_only)
  end subroutine int2_compute_t_set_screening

!###############################################################################

  subroutine int2_run(this, int2_consumer, stat, cam, alpha, beta, mu, &
                      alpha_coulomb, beta_coulomb)

    implicit none

    class(int2_compute_t), intent(inout) :: this
    class(int2_compute_data_t), intent(inout) :: int2_consumer
    logical, optional, intent(in) :: cam
    integer, optional, intent(out) :: stat
    real(kind=dp), optional, intent(in) :: alpha, beta, mu, &
                                           alpha_coulomb, beta_coulomb

    logical :: do_cam = .false.
    if (present(cam)) do_cam = cam

    if (present(stat)) stat = 0

    if (do_cam) then
      if (present(alpha).and.present(beta).and.present(mu)) then
        if (present(alpha_coulomb).and.present(beta_coulomb)) then
          call this%run_cam(int2_consumer, alpha, beta, mu, &
                            alpha_coulomb, beta_coulomb)
        else
          call this%run_cam(int2_consumer, alpha, beta, mu)
        end if
      else if (present(stat)) then
          stat = ERR_CAM_PARAM
      else
          call show_message("No CAM parameters given", WITH_ABORT)
      end if
    else
      call this%run_generic(int2_consumer)
    end if

  end subroutine int2_run

!###############################################################################

  subroutine int2_run_cam(this, int2_consumer, alpha, beta, mu, &
                          alpha_coulomb, beta_coulomb)

    implicit none

    class(int2_compute_t), intent(inout) :: this
    class(int2_compute_data_t), intent(inout) :: int2_consumer
    real(kind=dp), intent(in) :: alpha, beta, mu
    real(kind=dp), intent(in), optional :: alpha_coulomb, beta_coulomb

    logical :: attenuated_save
    real(kind=dp) :: mu_save

    attenuated_save = this%attenuated
    mu_save = this%mu

    int2_consumer%multipass = .true.
    int2_consumer%num_passes = 2

    int2_consumer%cur_pass = 1
    ! Regular Coulomb and exchange
    this%attenuated = .false.
    int2_consumer%scale_coulomb = 1.0d0
    if (present(alpha_coulomb)) &
       int2_consumer%scale_coulomb = alpha_coulomb
    int2_consumer%scale_exchange = alpha
    call this%run_generic(int2_consumer)

    int2_consumer%cur_pass = 2

    ! Short-range exchange:
    this%attenuated = .true.
    this%mu = mu
    int2_consumer%scale_coulomb = 0.0d0
    if (present(beta_coulomb)) &
       int2_consumer%scale_coulomb = beta_coulomb
    int2_consumer%scale_exchange = beta
    call this%run_generic(int2_consumer)

    int2_consumer%multipass = .false.
    int2_consumer%num_passes = 1
    int2_consumer%cur_pass = 1

    this%attenuated = attenuated_save
    this%mu = mu_save

  end subroutine int2_run_cam

!###############################################################################
!###############################################################################

  subroutine int2_twoei(this, int2_consumer)

    use int2e_libint, ONLY: libint2_init_eri, libint2_cleanup_eri
    use, intrinsic :: iso_c_binding, only: C_NULL_PTR, C_INT, c_int64_t
    use types, only: information
    use constants, only: NUM_CART_BF
    use blas_thread, only: blas_thread_count, blas_thread_set
!$  use omp_lib

    implicit none

    class(int2_compute_t), target, intent(inout) :: this
    class(int2_compute_data_t), intent(inout) :: int2_consumer

    integer :: i, j, k, l, ij_pair, npairs
    integer, allocatable :: pair_i(:), pair_j(:)
    integer :: lmax
    integer :: nint
    integer :: q4
    real(kind=dp) :: test

    integer :: nshell

    integer :: nthreads, ithread, thr_nshq, jork, nschwz
    real(kind=dp) :: tim0, tim1, tim2, tim3, tim4
    character(len=*), parameter :: &
        dbgfmt1 = '(/2x,&
                    &"Thread | Number of |",19X,"Timing",&
                    &/1x," number | quartets  |  Integrals |   F update  ",&
                    &"|   Schwartz  |    Total    |")', &
        dbgfmt2 = '(i5,4x,"|",i10," |",4(f9.2," s | "))'
    logical, parameter :: oflag = .false.
!    logical, parameter :: oflag = .true.
    logical :: omp
    logical :: zero_shq
    type(int2_storage_t) :: int2_storage
    type(eri_data_t), allocatable :: eri_data
    integer :: ok
    integer(c_int64_t) :: nBlasThreads

    nshell = this%basis%nshell
    npairs = nshell*(nshell+1)/2
    allocate(pair_i(npairs), pair_j(npairs))
    call int2_build_shell_pair_map(nshell, pair_i, pair_j)

    ! Hardwire BLAS to a single thread for the duration of the OpenMP 2e build.
    ! The Fock build is OpenMP-parallel and calls no BLAS itself, but a threaded
    ! BLAS keeps an idle worker pool that spins and oversubscribes the cores
    ! against the integral threads (measured ~1.5x slowdown at 28 threads). This
    ! uses the BLAS library's own runtime setter (openblas_set_num_threads /
    ! MKL_Set_Num_Threads / BLIS) so it CANNOT be overridden by a stray
    ! OPENBLAS_NUM_THREADS/MKL_NUM_THREADS in the environment.  The previous
    ! thread count is restored on exit, so diagonalisation and other BLAS-heavy
    ! phases outside this routine keep their full threading.  No-op (-1) for
    ! reference BLAS / Apple Accelerate where no setter is exported.
    nBlasThreads = -1
    nBlasThreads = blas_thread_count()
    if (nBlasThreads > 0) call blas_thread_set(1_c_int64_t)


    ! preparations for screening
    if (this%schwarz) then
      if (this%attenuated) then
        if (.not.allocated(this%schwarz_ints_attenuated)) then
          allocate(this%schwarz_ints_attenuated, mold=this%schwarz_ints_regular)
          call ints_exchange(this%basis, this%schwarz_ints_attenuated, this%mu**2, &
                             rys_only=this%rys_only)
        end if
        this%schwarz_ints => this%schwarz_ints_attenuated
      else
        this%schwarz_ints => this%schwarz_ints_regular
      end if
    end if

    omp = .false.
!$  omp = .true.

    nschwz = 0
    lmax = maxval(this%basis%am)
    if (lmax < 0 .or. lmax > 6) &
            call show_message("Basis set agular momentum exceeds max. supported", WITH_ABORT)

    if (oflag) write(*,dbgfmt1)

!$omp parallel &
!$omp   private( &
!$omp   i, j, k, l, ij_pair, jork, q4, &
!$omp   tim0, tim1, tim2, tim3, tim4, ithread,   &
!$omp   test, &
!$omp   int2_storage, &
!$omp   eri_data, &
!$omp   nthreads, &
!$omp   zero_shq) &
!$omp   shared(int2_consumer) &
!$omp   reduction(+:nschwz, nint, thr_nshq)

    nint = 0
    ithread  = 0
    nthreads = 1
!$  nthreads = omp_get_num_threads()
!$  ithread  = omp_get_thread_num()

!$  if (oflag) then
!$    thr_nshq = 0
!$    tim1     = 0.0d0
!$    tim2     = 0.0d0
!$    tim3     = 0.0d0
!$    tim4     = 0.0d0
!$    tim0     = omp_get_wtime()
!$  end if

!$omp master
    call int2_consumer%parallel_start(this%basis, nthreads)
!$omp end master
!$omp barrier

    allocate(eri_data)
    allocate(eri_data%ints(NUM_CART_BF(lmax)**4), source=0.0d0)
    allocate(eri_data%gdat)

    eri_data%attenuated_ints = this%attenuated
    eri_data%rys_only = this%rys_only
    eri_data%mu2 = this%mu**2

    if (libint2_active) then
        allocate(eri_data%erieval(this%basis%mxcontr**4))
        call libint2_init_eri(eri_data%erieval, int(4, C_INT), C_NULL_PTR)
    end if
    call eri_data%gdat%init(lmax, this%cutoffs, ok)

    call int2_storage%init(this%buf_size)
    int2_storage%thread_id = ithread + 1

!$omp barrier

! Distribute shell pairs through one flat dynamic workshare per int2_twoei
! parallel region.  The previous structure opened a fresh dynamic workshare
! for every (i,j) shell pair and needed a barrier after every shell pair to
! keep OpenMP runtimes from overlapping workshares with different k-loop
! bounds.  A flat shell-pair workshare keeps scheduling granular without any
! synchronization point inside the shell-pair loop.
!$omp do schedule(dynamic,1)
    do ij_pair = 1, npairs
        i = pair_i(ij_pair)
        j = pair_j(ij_pair)
        if (this%pe%size>1) then
          if (mod(ij_pair, this%pe%size) /= this%pe%rank) cycle
        end if

        if (this%schwarz) then
          test = int2_consumer%screen_ij(this%schwarz_ints, i, j)
          ! With the petite list active the surviving representative carries
          ! up to |G| weight, so the pair-level skip must be conservative by
          ! the same factor (orbit members of non-abelian operations live in
          ! different shell pairs with different bounds).
          if (this%petite .and. this%sym_full) test = test*real(this%sym_nops, dp)
          if (test < this%cutoffs%integral_cutoff) then
            nschwz = nschwz + i*(i-1)/2+j
            cycle
          end if
        end if

        do k = 1, i

          jork = k
          if (i==k) jork = j

          do l = 1, jork

!$          if (oflag) tim1 = omp_get_wtime()
            ! Petite list: keep only the orbit representative, weighted by
            ! the orbit size; the skeleton Fock is symmetrized afterwards.
            if (this%petite) then
              q4 = petite_quartet_weight(this%sym_shell_map, this%sym_nops, &
                                         nshell, i, j, k, l)
              if (q4 == 0) cycle
              eri_data%weight = real(q4, dp)
              eri_data%weighted_cutoff = this%sym_full
            else
              eri_data%weight = 1.0d0
              eri_data%weighted_cutoff = .false.
            end if

            if (this%schwarz) then
              test = int2_consumer%screen_ijkl(this%schwarz_ints, i, j, k, l)
              ! Screen the weighted contribution: non-abelian orbit members
              ! have unequal element magnitudes, so the unweighted threshold
              ! would leak a systematic cutoff-level error into the skeleton.
              if (eri_data%weighted_cutoff) test = test*eri_data%weight
              if (test < this%cutoffs%integral_cutoff) then
                nschwz = nschwz+1
                cycle
              end if
            end if

            thr_nshq = thr_nshq + 1
!$          if (oflag) tim2 = tim2 + omp_get_wtime() - tim1

            eri_data%ids(:) = [i,j,k,l]
            call shellquartet(this%basis, this%ppairs, this%cutoffs, eri_data, zero_shq)
!$          if (oflag) tim3 = tim3 + omp_get_wtime() - tim1

            if (zero_shq) cycle

            call int2_compute_data_t_storeints(int2_consumer, &
                    this%basis, eri_data, int2_storage, this%cutoffs%integral_cutoff, nint)
!$          if (oflag) tim4 = tim4 + omp_get_wtime() - tim1

          end do
        end do
    end do
!$omp end do

    call int2_consumer%update(int2_storage)

!$  if (oflag) tim0 = omp_get_wtime() - tim0

!  Debug timing output
    if (omp.and.oflag) then
!$omp do ordered
      do i = 0, nthreads-1
!$omp ordered
        write(*,dbgfmt2) ithread,thr_nshq,tim3-tim2,tim4-tim3,tim2,tim0
!$omp end ordered
      end do
!$omp end do
    end if

    if (libint2_active) then
        call libint2_cleanup_eri(eri_data%erieval)
        deallocate(eri_data%erieval)
    end if
    call eri_data%gdat%clean()
    call int2_storage%clean()

    deallocate(eri_data)
!$omp end parallel
    ! restore BLAS threading for diagonalisation / other BLAS-heavy phases
    if (nBlasThreads > 0) call blas_thread_set(nBlasThreads)
    deallocate(pair_i, pair_j)
    call int2_consumer%pe%init(this%pe%comm, this%pe%use_mpi)
    call int2_consumer%parallel_stop()

    this%skipped = nschwz

  contains

    subroutine int2_build_shell_pair_map(nshell, shell_pair_i, shell_pair_j)
      !> Build the flat shell-pair work list, ordered by DESCENDING estimated cost
      !> so the dynamic OpenMP schedule hands out the expensive quartets first and
      !> the cheap ones fill the tail -- minimising the end-of-region load-imbalance
      !> barrier wait ("adaptive dynamic dispatch").  Cost ~ nbf_i*nbf_j*i captures
      !> both the per-quartet kernel size (angular momentum) and the k,l iteration
      !> count (~i), and -- when available -- the Schwarz magnitude of the bra pair
      !> (its SPARSITY: diffuse/negligible pairs are mostly screened and do little
      !> work, so they sink to the tail).  A counting sort over log2(cost) classes
      !> keeps this O(n); reordering only changes work distribution, never results.
      implicit none
      integer, intent(in) :: nshell
      integer, intent(out) :: shell_pair_i(:), shell_pair_j(:)
      integer :: i, j, ij_pair, c, nbfi, nbfj, ami, amj
      integer, parameter :: NCLASS = 64, OFFS = 42
      integer :: cnt(0:NCLASS), off(0:NCLASS), cls
      real(kind=dp) :: cost, sw
      logical :: use_sw

      use_sw = this%schwarz .and. allocated(this%schwarz_ints_regular)

      ! Pass 1: count pairs per cost-class (class 0 = most expensive).
      cnt = 0
      do i = nshell, 1, -1
        ami = this%basis%am(i); nbfi = (ami+1)*(ami+2)/2
        do j = 1, i
          amj = this%basis%am(j); nbfj = (amj+1)*(amj+2)/2
          sw = 1.0d0
          if (use_sw) sw = max(this%schwarz_ints_regular(i,j), 1.0d-30)
          cost = sw * real(nbfi*nbfj,dp) * real(i,dp)
          cls = NCLASS - max(0, min(NCLASS, int(log(cost)/log(2.0d0)) + OFFS))
          cnt(cls) = cnt(cls) + 1
        end do
      end do

      ! Prefix offsets (ascending class index = descending cost).
      off(0) = 0
      do c = 1, NCLASS
        off(c) = off(c-1) + cnt(c-1)
      end do

      ! Pass 2: scatter pairs into cost-sorted positions (stable within a class,
      ! preserving the i-descending build order as the secondary key).
      do i = nshell, 1, -1
        ami = this%basis%am(i); nbfi = (ami+1)*(ami+2)/2
        do j = 1, i
          amj = this%basis%am(j); nbfj = (amj+1)*(amj+2)/2
          sw = 1.0d0
          if (use_sw) sw = max(this%schwarz_ints_regular(i,j), 1.0d-30)
          cost = sw * real(nbfi*nbfj,dp) * real(i,dp)
          cls = NCLASS - max(0, min(NCLASS, int(log(cost)/log(2.0d0)) + OFFS))
          off(cls) = off(cls) + 1
          ij_pair = off(cls)
          shell_pair_i(ij_pair) = i
          shell_pair_j(ij_pair) = j
        end do
      end do
    end subroutine int2_build_shell_pair_map

  end subroutine int2_twoei

!###############################################################################
!###############################################################################

  function int2_compute_data_t_screen_ij(this, xints, i, j) result(res)
    implicit none
    class(int2_compute_data_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j
    res = 1
    return
    UNUSED_DUMMY( this  )
    UNUSED_DUMMY( xints )
    UNUSED_DUMMY( i )
    UNUSED_DUMMY( j )
  end function int2_compute_data_t_screen_ij

!###############################################################################

  function int2_compute_data_t_screen_ijkl(this, xints, i, j, k, l) result(res)
    implicit none
    class(int2_compute_data_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j, k, l
    res = 1
    return
    UNUSED_DUMMY( this  )
    UNUSED_DUMMY( xints )
    UNUSED_DUMMY( i )
    UNUSED_DUMMY( j )
    UNUSED_DUMMY( k )
    UNUSED_DUMMY( l )
  end function int2_compute_data_t_screen_ijkl

!###############################################################################
!###############################################################################

  function int2_fock_data_t_screen_ij(this, xints, i, j) result(res)
    implicit none
    class(int2_fock_data_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j
    integer :: ij
    res = xints(i,j) * this%max_den
  end function int2_fock_data_t_screen_ij

!###############################################################################

  function int2_fock_data_t_screen_ijkl(this, xints, i, j, k, l) result(res)
    implicit none
    class(int2_fock_data_t), intent(in) :: this
    real(kind=dp), contiguous, intent(in) :: xints(:,:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j, k, l
    integer :: ij, kl
    res = xints(i,j)*xints(k,l)
    if (allocated(this%dsh)) then
      res = res * max(4*this%dsh(i,j), 4*this%dsh(k,l), this%dsh(j,l), this%dsh(j,k), this%dsh(i,l), this%dsh(i,k))
    end if
  end function int2_fock_data_t_screen_ijkl

!###############################################################################
!> @brief Computes the shell density matrix from the given density matrix and basis set.
!>
!> @details This subroutine calculates the shell density matrix `shell_density`
!>          by determining the maximum density values for each shell pair from
!>          the provided density matrix `density_matrix` and the basis set information.
!>
!> @param[out] shell_density  The computed shell density matrix.
!> @param[in]  density_matrix The input density matrix in AO basis.
!> @param[in]  basis          The basis set information including shell data.
!>
  subroutine shlden(shell_density, density_matrix, basis)
    use basis_tools, only: basis_set
    implicit none
    type(basis_set), intent(in) :: basis             !< Basis set information
    real(kind=dp), intent(out) :: shell_density(:,:) !< Density matrix compressed to shells
    real(kind=dp), intent(in) :: density_matrix(:,:) !< Density matrix in AO basis

    integer :: shell_i, shell_j, i, j, ij
    integer :: mini, maxi, minj, maxj, fock_index
    real(kind=dp) :: dmax

    ! Initialize shell density to zero
    shell_density = 0

    ! Loop over each Fock matrix
    do fock_index = 1, ubound(density_matrix,2)

      ! Loop over shells
      do shell_i = 1, basis%nshell
        mini = basis%ao_offset(shell_i)
        maxi = mini + basis%naos(shell_i) - 1

        ! Loop over shells again for symmetry
        do shell_j = 1, shell_i
          minj = basis%ao_offset(shell_j)
          maxj = minj + basis%naos(shell_j) - 1
          dmax = 0.0d0

          ! Loop over basis functions within shells
          do i = mini, maxi
            if (shell_i == shell_j) maxj = i
            do j = minj, maxj
              ij = i*(i-1)/2 + j
              dmax = max(abs(density_matrix(ij,fock_index)), dmax)
            end do
          end do

          shell_density(shell_i, shell_j) = max(dmax, shell_density(shell_i, shell_j))
        end do
      end do
    end do

    ! Symmetrize
    do shell_i = 1, basis%nshell
      do shell_j = shell_i+1, basis%nshell
        shell_density(shell_i, shell_j) = shell_density(shell_j, shell_i)
      end do
    end do
  end subroutine shlden

!###############################################################################

  subroutine shellquartet(basis, ppairs, cutoffs, eri_data, zero_shq)
    use io_constants, only: iw
    use int2e_rotaxis, only: genr22, genr22_pure, genr22_reduce_pure
    use int2e_libint, only: libint_compute_eri, libint_print_eri
    use int2e_rys, only: int2_rys_compute, int2_rys_reduce_pure, rys_print_eri
    use iso_c_binding, only: c_f_pointer
    use constants, only: HARMONIC_ACTIVE
    use int2_pure_generated, only: int2_project_pure_block
    implicit none
    type(basis_set), intent(in) :: basis
    type(int2_pair_storage), intent(in) :: ppairs
    type(int2_cutoffs_t), intent(in) :: cutoffs
    type(eri_data_t), intent(inout), target :: eri_data
    logical, intent(out) :: zero_shq
    logical :: rotspd, libint, rys
    integer :: nbf(4)
    integer :: s_, fp_, orig_, am_s(4), pure_s(4), nbf_s(4), nbf_out_s(4)
    logical, parameter :: dbg_output = .false.
!    logical, parameter :: dbg_output = .true.
    integer :: max_am
    logical :: err

    zero_shq = .false.

    eri_data%am = basis%am(eri_data%ids)
    max_am = maxval(eri_data%am)
    eri_data%nbf = (eri_data%am+1)*(eri_data%am+2)/2

    rotspd = max_am <= 2
    libint = .not.rotspd.and.libint2_active.and..not.eri_data%attenuated_ints &
             .and..not.eri_data%rys_only
    rys = .not.rotspd.and..not.libint

    if (rotspd) then

      if (HARMONIC_ACTIVE .and. &
          any(basis%harmonic(eri_data%ids) == 1 .and. eri_data%am >= 2)) then
        ! direct pure-spherical output: lab rotation and c2s fused per index
        if (eri_data%attenuated_ints) then
          call genr22_pure(basis, ppairs, eri_data%ints, eri_data%ids, eri_data%flips, &
                           cutoffs, eri_data%nbf, eri_data%mu2)
        else
          call genr22_pure(basis, ppairs, eri_data%ints, eri_data%ids, eri_data%flips, &
                           cutoffs, eri_data%nbf)
        end if
      else
        if (eri_data%attenuated_ints) then
          ! erf-attenuated integrals
          call genr22(basis, ppairs, eri_data%ints, eri_data%ids, eri_data%flips, cutoffs, eri_data%mu2)
        else
          ! regular integrals
          call genr22(basis, ppairs, eri_data%ints, eri_data%ids, eri_data%flips, cutoffs)
        end if

        eri_data%nbf = eri_data%nbf(eri_data%flips)
        call genr22_reduce_pure(basis, eri_data%ids, eri_data%flips, eri_data%ints, eri_data%nbf)
      end if
      eri_data%pints(1:eri_data%nbf(4), 1:eri_data%nbf(3), 1:eri_data%nbf(2), 1:eri_data%nbf(1)) => eri_data%ints

    else if (libint) then

      call libint_compute_eri(basis, ppairs, cutoffs, eri_data%ids, 0, eri_data%erieval, eri_data%flips, zero_shq)
      if (zero_shq) return
      nbf = eri_data%nbf(eri_data%flips)
      eri_data%nbf = nbf
      call c_f_pointer(eri_data%erieval(1)%targets(1), eri_data%pints, shape=nbf([4,3,2,1]))
      call normalize_ints(nbf, eri_data%am(eri_data%flips), eri_data%pints)

    else if (rys) then

      call eri_data%gdat%set_ids(basis, eri_data%ids)
      if (eri_data%attenuated_ints) then
        ! erf-attenuated integrals
        call int2_rys_compute(eri_data%ints, eri_data%gdat, ppairs, zero_shq, &
                              mu2=eri_data%mu2, basis=basis, direct_pure=.true.)
      else
        ! regular integrals
        call int2_rys_compute(eri_data%ints, eri_data%gdat, ppairs, zero_shq, &
                              basis=basis, direct_pure=.true.)
      end if

      if (zero_shq) return
      nbf = eri_data%gdat%nbf
      eri_data%flips = eri_data%gdat%flips
      if (eri_data%gdat%direct_pure) then
        eri_data%nbf = nbf
        eri_data%pints(1:nbf(4), 1:nbf(3), 1:nbf(2), 1:nbf(1)) => eri_data%ints
      else
        eri_data%pints(1:nbf(4), 1:nbf(3), 1:nbf(2), 1:nbf(1)) => eri_data%ints
        call normalize_ints(nbf, eri_data%gdat%am, eri_data%pints)
        call int2_rys_reduce_pure(basis, eri_data%gdat, eri_data%ints, nbf)
        eri_data%nbf = nbf
        eri_data%pints(1:nbf(4), 1:nbf(3), 1:nbf(2), 1:nbf(1)) => eri_data%ints
      end if

    end if

    ! Libint returns a normalized Cartesian target buffer. Stage it into the
    ! owned ERI buffer before reducing harmonic-flagged shells to pure output.
    if (HARMONIC_ACTIVE .and. libint) then
      do s_ = 1, 4
        fp_ = 5 - s_                       ! pints storage dim s_ <-> flipped position fp_
        orig_ = eri_data%flips(fp_)        ! original shell index at that position
        am_s(s_)   = eri_data%am(orig_)
        pure_s(s_) = basis%harmonic(eri_data%ids(orig_))
        nbf_s(s_)  = eri_data%nbf(fp_)
      end do
      if (any(pure_s == 1 .and. am_s >= 2)) then
        eri_data%ints(1:product(nbf_s)) = reshape(eri_data%pints, [product(nbf_s)])
        call int2_project_pure_block(eri_data%ints, am_s, pure_s, nbf_s, nbf_out_s)
        do s_ = 1, 4
          eri_data%nbf(5 - s_) = nbf_out_s(s_)
        end do
        eri_data%pints(1:eri_data%nbf(4), 1:eri_data%nbf(3), &
                       1:eri_data%nbf(2), 1:eri_data%nbf(1)) => eri_data%ints
      end if
    end if

    if (dbg_output) then
      block
          integer :: len1, len2, len3, len4
          write (iw, '("shells", 4i5, " :", 4i5, " :", 4i5, " :", 4i5)') &
                  eri_data%ids, &
                  eri_data%am, &
                  eri_data%am(eri_data%flips), &
                  eri_data%flips
      end block
      if (libint) call libint_print_eri(basis, eri_data%ids, 0, eri_data%erieval, eri_data%flips)
      if (rys) call rys_print_eri(eri_data%gdat, eri_data%pints)

    end if

  end subroutine shellquartet

!###############################################################################

  subroutine normalize_ints(nbf, am, ints)
    use constants, only: shells_pnrm2
    implicit none
    integer, intent(in) :: nbf(4), am(4)
    real(kind=dp), contiguous, intent(inout) :: ints(:,:,:,:)

    integer :: na, nb, nc, nd
    real(kind=dp), pointer:: pnorma(:), pnormb(:), pnormc(:), pnormd(:)

    pnorma => shells_pnrm2(:,am(1))
    pnormb => shells_pnrm2(:,am(2))
    pnormc => shells_pnrm2(:,am(3))
    pnormd => shells_pnrm2(:,am(4))

    do concurrent (na = 1:nbf(1), nb = 1:nbf(2), nc = 1:nbf(3), nd = 1:nbf(4))
      ints(nd,nc,nb,na) = ints(nd,nc,nb,na) &
                        * pnorma(na) * pnormb(nb) &
                        * pnormc(nc) * pnormd(nd)
    end do

  end subroutine

!###############################################################################

  subroutine int2_storage_init(this, buf_size)

    implicit none

    class(int2_storage_t), intent(inout) :: this
    integer, intent(in) :: buf_size

    if(allocated(this%ids)) call this%clean()

    this%buf_size = buf_size
    this%ncur = 0

    allocate( this%ids  (4,this%buf_size) )
    allocate( this%ints(this%buf_size) )

  end subroutine

!###############################################################################

  subroutine int2_storage_clean(this)

    implicit none

    class(int2_storage_t), intent(inout) :: this

    deallocate(this%ids)
    deallocate(this%ints)

    this%buf_size = 0
    this%ncur = 0

  end subroutine

!###############################################################################
!###############################################################################

  subroutine int2_fock_data_t_parallel_start(this, basis, nthreads)
    implicit none
    class(int2_fock_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads
    integer :: nsh, ncopy, si, sj, npair_sh, nsig
    character(len=32) :: sval
    integer :: ln
    real(kind=dp) :: repl_mb, cap_mb, frac_sig, spthr, sptol

    this%nthreads = nthreads

    nsh = basis%nshell

    if (allocated(this%dsh)) then
        if (size(this%dsh) /= nsh*nsh) deallocate(this%dsh)
    end if

    if (.not.allocated(this%dsh)) then
        allocate(this%dsh(nsh,nsh), source=0.0d0)
    else
        this%dsh = 0
    end if

!   Form the shell density first -- it drives both screening and the Fock-memory
!   mode decision below.
    call this%init_screen(basis)

!   Decide Fock accumulator memory mode.  Replicated (one copy/thread) is fastest
!   but costs fockdim*nfocks*nthreads*8 bytes; for very large systems that blows
!   up.  A single shared atomic Fock uses ~Nthreads x less memory, but per-integral
!   atomics contend -- cheaply only when few quartets survive, i.e. when the DENSITY
!   IS SPARSE.  So switch to the low-memory atomic buffer only when (a) the
!   replicated buffers would exceed OQP_FOCK_MEM_MB (default 4096) AND (b) the
!   shell density is sparse (significant-pair fraction < OQP_FOCK_SPARSITY,
!   default 0.5).  Dense density keeps the fast replicated path.  OQP_FOCK_ATOMIC
!   forces the choice.
    cap_mb = 4096.0d0
    call get_environment_variable("OQP_FOCK_MEM_MB", sval, ln)
    if (ln > 0) read(sval,*,iostat=ln) cap_mb
    spthr = 0.5d0
    call get_environment_variable("OQP_FOCK_SPARSITY", sval, ln)
    if (ln > 0) read(sval,*,iostat=ln) spthr

    repl_mb = real(this%fockdim,dp)*this%nfocks*nthreads*8.0d0/1.048576d6

    ! density sparsity = fraction of shell pairs carrying significant density
    sptol = 1.0d-4
    npair_sh = nsh*(nsh+1)/2
    nsig = 0
    do si = 1, nsh
      do sj = 1, si
        if (this%dsh(si,sj) > sptol*this%max_den) nsig = nsig + 1
      end do
    end do
    frac_sig = real(nsig,dp) / real(max(1,npair_sh),dp)

    this%atomic_fock = (nthreads > 1) .and. (repl_mb > cap_mb) .and. (frac_sig < spthr)
    call get_environment_variable("OQP_FOCK_ATOMIC", sval, ln)
    if (ln > 0) then
      this%atomic_fock = (sval(1:1)=='1' .or. sval(1:1)=='y' .or. sval(1:1)=='Y' &
                          .or. sval(1:1)=='t' .or. sval(1:1)=='T')
    end if
    ncopy = nthreads
    if (this%atomic_fock) ncopy = 1

    if (this%cur_pass == 1 .and. this%atomic_fock) then
      write(*,'(2x,a,f9.1,a,f9.1,a,i0,a,f5.2,a)') &
        "Fock accumulator: shared+atomic (low-memory) using ", &
        real(this%fockdim,dp)*this%nfocks*8.0d0/1.048576d6, " MB vs ", &
        repl_mb, " MB replicated (", nthreads, " threads), density frac_sig=", &
        frac_sig, ""
    end if

    if (this%cur_pass == 1) then
      if (allocated(this%f)) then
          if ( any((shape(this%f) - [this%fockdim, this%nfocks, ncopy])/=0) ) then
              deallocate(this%f)
          end if
      end if
      if (.not.allocated(this%f)) then
          allocate(this%f(this%fockdim, this%nfocks, ncopy), source=0.0d0)
      else
          this%f = 0
      end if
    end if

  end subroutine

!###############################################################################

  subroutine int2_fock_data_t_init_screen(this, basis)
    implicit none
    class(int2_fock_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis

!   Form shell density
    call shlden(this%dsh, this%d, basis)
    this%max_den = maxval(abs(this%dsh))

  end subroutine

!###############################################################################

  subroutine int2_rhf_data_t_parallel_start(this, basis, nthreads)
    implicit none
    class(int2_rhf_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads

    this%fockdim = basis%nbf*(basis%nbf+1) / 2
    this%nfocks = ubound(this%d, size(shape(this%d)))

    call this%int2_fock_data_t_parallel_start(basis, nthreads)

  end subroutine

!###############################################################################

  subroutine int2_urohf_data_t_parallel_start(this, basis, nthreads)

    implicit none
    class(int2_urohf_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: nthreads


    this%fockdim = basis%nbf*(basis%nbf+1) / 2
    this%nfocks = 2

    call this%int2_fock_data_t_parallel_start(basis, nthreads)

  end subroutine

!###############################################################################

  subroutine int2_fock_data_t_parallel_stop(this)

    implicit none
    class(int2_fock_data_t), intent(inout) :: this

    call this%pe%barrier()
    if (this%cur_pass /= this%num_passes) return
    ! atomic_fock already accumulated into the single shared copy (dim3==1);
    ! only the replicated mode needs the cross-thread reduction.
    if (this%nthreads /= 1 .and. .not.this%atomic_fock) then
      this%f(:,:,lbound(this%f,3)) = sum(this%f, dim=size(shape(this%f)))
    end if

    call this%pe%allreduce(this%f(:,:,1), &
                         size(this%f(:,:,1)))
    call this%pe%barrier()
    this%nthreads = 1
  end subroutine

!###############################################################################

  subroutine int2_fock_data_t_clean(this)
    implicit none
    class(int2_fock_data_t), intent(inout) :: this
    deallocate(this%f)
    deallocate(this%dsh)
    nullify(this%d)
  end subroutine

!###############################################################################

  subroutine int2_rhf_data_t_update(this, buf)
    implicit none
    class(int2_rhf_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: ii, jj, kk, ll, ij, ik, il, jk, jl, kl, n, ii2, jj2, kk2
    real(kind=dp) :: xval1, xval4, val, val1, val4
    real(kind=dp) :: aij, akl, aik, ajl, ail, ajk
    integer :: ifock, mythread

    xval1 = this%scale_exchange
    xval4 = 4 * this%scale_coulomb
    mythread = buf%thread_id
    if (this%atomic_fock) mythread = 1

    do ifock = 1, this%nfocks
      do n = 1, buf%ncur
        ii = buf%ids(1,n)
        jj = buf%ids(2,n)
        kk = buf%ids(3,n)
        ll = buf%ids(4,n)
        val = buf%ints(n)

        ii2 = ii*(ii-1)/2
        jj2 = jj*(jj-1)/2
        kk2 = kk*(kk-1)/2

        ij = ii2+jj
        ik = ii2+kk
        il = ii2+ll
        jk = jj2+kk
        jl = jj2+ll
        kl = kk2+ll
        if (jj<kk) jk = kk2 + jj
        if (jj<ll) jl = ll*(ll-1)/2 + jj

        val1 = val*xval1
        val4 = val*xval4

        if (this%atomic_fock) then
          ! single shared Fock: atomic accumulation (low-memory mode).
          ! Contributions are precomputed into locals so the atomic statement's
          ! RHS references no component of `this` (gfortran atomic requirement).
          aij = val4*this%d(kl,ifock); akl = val4*this%d(ij,ifock)
          aik = -val1*this%d(jl,ifock); ajl = -val1*this%d(ik,ifock)
          ail = -val1*this%d(jk,ifock); ajk = -val1*this%d(il,ifock)
          !$omp atomic update
          this%f(ij,ifock,1) = this%f(ij,ifock,1) + aij
          !$omp atomic update
          this%f(kl,ifock,1) = this%f(kl,ifock,1) + akl
          !$omp atomic update
          this%f(ik,ifock,1) = this%f(ik,ifock,1) + aik
          !$omp atomic update
          this%f(jl,ifock,1) = this%f(jl,ifock,1) + ajl
          !$omp atomic update
          this%f(il,ifock,1) = this%f(il,ifock,1) + ail
          !$omp atomic update
          this%f(jk,ifock,1) = this%f(jk,ifock,1) + ajk
        else
          this%f(ij,ifock,mythread) = this%f(ij,ifock,mythread) + val4*this%d(kl,ifock)
          this%f(kl,ifock,mythread) = this%f(kl,ifock,mythread) + val4*this%d(ij,ifock)
          this%f(ik,ifock,mythread) = this%f(ik,ifock,mythread) - val1*this%d(jl,ifock)
          this%f(jl,ifock,mythread) = this%f(jl,ifock,mythread) - val1*this%d(ik,ifock)
          this%f(il,ifock,mythread) = this%f(il,ifock,mythread) - val1*this%d(jk,ifock)
          this%f(jk,ifock,mythread) = this%f(jk,ifock,mythread) - val1*this%d(il,ifock)
        end if
      end do
    end do

    buf%ncur = 0

  end subroutine int2_rhf_data_t_update

!###############################################################################

  subroutine int2_urohf_data_t_update(this, buf)
    implicit none
    class(int2_urohf_data_t), intent(inout) :: this
    type(int2_storage_t), intent(inout) :: buf
    integer :: ii, jj, kk, ll, ij, ik, il, jk, jl, kl, n, ii2, jj2, kk2
    real(kind=dp) :: xval2, xval4, val, val1, val4, cij, ckl
    real(kind=dp) :: a1ik, a1jl, a1il, a1jk, a2ik, a2jl, a2il, a2jk
    integer :: mythread

    xval2 = 2 * this%scale_exchange
    xval4 = 4 * this%scale_coulomb

    mythread = buf%thread_id
    if (this%atomic_fock) mythread = 1

      do n = 1, buf%ncur
        ii = buf%ids(1,n)
        jj = buf%ids(2,n)
        kk = buf%ids(3,n)
        ll = buf%ids(4,n)
        val = buf%ints(n)

      ii2 = ii*(ii-1)/2
      jj2 = jj*(jj-1)/2
      kk2 = kk*(kk-1)/2

      ij = ii2+jj
      ik = ii2+kk
      il = ii2+ll
      jk = jj2+kk
      jl = jj2+ll
      kl = kk2+ll
      if (jj<kk) jk = kk2 + jj
      if (jj<ll) jl = ll*(ll-1)/2 + jj

      val1 = val*xval2
      val4 = val*xval4

      cij = val4*sum(this%d(ij,1:2))
      ckl = val4*sum(this%d(kl,1:2))

      if (this%atomic_fock) then
        ! locals so atomic RHS references no component of `this`
        a1ik = -val1*this%d(jl,1); a1jl = -val1*this%d(ik,1)
        a1il = -val1*this%d(jk,1); a1jk = -val1*this%d(il,1)
        a2ik = -val1*this%d(jl,2); a2jl = -val1*this%d(ik,2)
        a2il = -val1*this%d(jk,2); a2jk = -val1*this%d(il,2)
        !$omp atomic update
        this%f(ij,1,1) = this%f(ij,1,1) + ckl
        !$omp atomic update
        this%f(kl,1,1) = this%f(kl,1,1) + cij
        !$omp atomic update
        this%f(ik,1,1) = this%f(ik,1,1) + a1ik
        !$omp atomic update
        this%f(jl,1,1) = this%f(jl,1,1) + a1jl
        !$omp atomic update
        this%f(il,1,1) = this%f(il,1,1) + a1il
        !$omp atomic update
        this%f(jk,1,1) = this%f(jk,1,1) + a1jk
        !$omp atomic update
        this%f(ij,2,1) = this%f(ij,2,1) + ckl
        !$omp atomic update
        this%f(kl,2,1) = this%f(kl,2,1) + cij
        !$omp atomic update
        this%f(ik,2,1) = this%f(ik,2,1) + a2ik
        !$omp atomic update
        this%f(jl,2,1) = this%f(jl,2,1) + a2jl
        !$omp atomic update
        this%f(il,2,1) = this%f(il,2,1) + a2il
        !$omp atomic update
        this%f(jk,2,1) = this%f(jk,2,1) + a2jk
      else
      this%f(ij,1,mythread) = this%f(ij,1,mythread) + ckl
      this%f(kl,1,mythread) = this%f(kl,1,mythread) + cij
      this%f(ik,1,mythread) = this%f(ik,1,mythread) - val1*this%d(jl,1)
      this%f(jl,1,mythread) = this%f(jl,1,mythread) - val1*this%d(ik,1)
      this%f(il,1,mythread) = this%f(il,1,mythread) - val1*this%d(jk,1)
      this%f(jk,1,mythread) = this%f(jk,1,mythread) - val1*this%d(il,1)

      this%f(ij,2,mythread) = this%f(ij,2,mythread) + ckl
      this%f(kl,2,mythread) = this%f(kl,2,mythread) + cij
      this%f(ik,2,mythread) = this%f(ik,2,mythread) - val1*this%d(jl,2)
      this%f(jl,2,mythread) = this%f(jl,2,mythread) - val1*this%d(ik,2)
      this%f(il,2,mythread) = this%f(il,2,mythread) - val1*this%d(jk,2)
      this%f(jk,2,mythread) = this%f(jk,2,mythread) - val1*this%d(il,2)
      end if
    end do

    buf%ncur = 0

  end subroutine int2_urohf_data_t_update

!###############################################################################

  subroutine ints_exchange(basis, schwarz_ints, mu2, rys_only)
    use int2e_rotaxis, only: genr22, genr22_pure, genr22_reduce_pure
    use int2e_libint, only: libint2_init_eri, libint2_cleanup_eri
    use int2e_libint, only: libint_compute_eri, libint_print_eri
    use int2e_libint, only: libint_t, libint2_active
    use int2e_rys, only: int2_rys_compute, int2_rys_reduce_pure
    use types, only: information
    use constants, only: HARMONIC_ACTIVE, NUM_CART_BF
    use int2_pure_generated, only: int2_project_pure_block
    use, intrinsic :: iso_c_binding, only: C_NULL_PTR, C_INT,  c_f_pointer

    implicit none

    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(inout) :: schwarz_ints(:,:)
    real(kind=dp), optional, intent(in) :: mu2
    logical, optional, intent(in) :: rys_only

    real(kind=dp), parameter :: &
      ic_exchng  = 1.0d-15, &
      ei1_exchng = 1.0d-17, &
      ei2_exchng = 1.0d-17, &
      cux_exchng = 50.0

    integer :: flips(4)
    integer :: shell_ids(4)
    integer :: nbf(4)
    integer :: lmax
    integer :: ish, jsh
    integer :: i, j, jmax
    integer :: am(4), max_am
    integer :: s_, fp_, orig_, am_s(4), pure_s(4), nbf_s(4), nbf_out_s(4)
    integer :: ok
    logical :: rotspd, libint, zero_shq, rys
    logical :: attenuated, rys_only_
    real(kind=dp) :: vmax
    real(kind=dp), allocatable, target :: ints(:)
    real(kind=dp), pointer :: pints(:,:,:,:)
    type(libint_t), allocatable :: erieval(:)
    type(int2_rys_data_t) :: gdat
    type(int2_cutoffs_t) :: cutoffs
    type(int2_pair_storage) :: ppairs

    attenuated = present(mu2)
    rys_only_ = .false.
    if (present(rys_only)) rys_only_ = rys_only

    lmax = maxval(basis%am)
    if (lmax < 0 .or. lmax > 6) call show_message("Basis set agular momentum exceeds max. supported", WITH_ABORT)

!   Set very tight cutoff
    call cutoffs%set(ic_exchng, ei1_exchng, ei2_exchng, cux_exchng)

    if (libint2_active) then
        allocate(erieval(basis%mxcontr**4))
        call libint2_init_eri(erieval, int(4, C_INT), C_NULL_PTR)
    end if
    allocate(ints(NUM_CART_BF(lmax)**4), source=0.0d0)
    call gdat%init(lmax, cutoffs, ok)

    call ppairs%alloc(basis, cutoffs)
    call ppairs%compute(basis, cutoffs)

    do ish = 1, basis%nshell
      do jsh = 1, ish
        shell_ids = [ish, jsh, ish, jsh]

        am = basis%am(shell_ids)
        max_am = maxval(am)

        rotspd = max_am <= 2
        libint = .not.rotspd.and.libint2_active.and..not.attenuated.and..not.rys_only_
        rys = .not.rotspd.and..not.libint
        if (rotspd) then
          if (HARMONIC_ACTIVE .and. any(basis%harmonic(shell_ids) == 1 .and. am >= 2)) then
            ! direct pure-spherical output: lab rotation and c2s fused per index
            if (attenuated) then
              call genr22_pure(basis, ppairs, ints, shell_ids, flips, cutoffs, nbf, mu2)
            else
              call genr22_pure(basis, ppairs, ints, shell_ids, flips, cutoffs, nbf)
            end if
          else
            if (attenuated) then
              call genr22(basis, ppairs, ints, shell_ids, flips, cutoffs, mu2)
            else
              call genr22(basis, ppairs, ints, shell_ids, flips, cutoffs)
            end if
            nbf = NUM_CART_BF(am(flips))
            call genr22_reduce_pure(basis, shell_ids, flips, ints, nbf)
          end if
          vmax = maxval(abs(ints(1:product(nbf))))
        else if (libint) then
          call libint_compute_eri(basis, ppairs, cutoffs, shell_ids, 0, erieval, flips, zero_shq)
          ! flips is intent(out) of libint_compute_eri; only valid afterwards
          nbf = NUM_CART_BF(am(flips))
          if (zero_shq) then
            vmax = 0.0_dp
          else
            call c_f_pointer(erieval(1)%targets(1), pints, shape=nbf([4,3,2,1]))
            call normalize_ints(nbf, am(flips), pints)
            if (HARMONIC_ACTIVE) then
              do s_ = 1, 4
                fp_ = 5 - s_
                orig_ = flips(fp_)
                am_s(s_)   = am(orig_)
                pure_s(s_) = basis%harmonic(shell_ids(orig_))
                nbf_s(s_)  = nbf(fp_)
              end do
              if (any(pure_s == 1 .and. am_s >= 2)) then
                ints(1:product(nbf_s)) = reshape(pints, [product(nbf_s)])
                call int2_project_pure_block(ints, am_s, pure_s, nbf_s, nbf_out_s)
                do s_ = 1, 4
                  nbf(5 - s_) = nbf_out_s(s_)
                end do
                vmax = maxval(abs(ints(1:product(nbf))))
              else
                vmax = maxval(abs(pints))
              end if
            else
              vmax = maxval(abs(pints))
            end if
          end if
        else if (rys) then
          call gdat%set_ids(basis, shell_ids)
          if (attenuated) then
            call int2_rys_compute(ints, gdat, ppairs, zero_shq, &
                                  mu2=mu2, basis=basis, direct_pure=.true.)
          else
            call int2_rys_compute(ints, gdat, ppairs, zero_shq, &
                                  basis=basis, direct_pure=.true.)
          end if
          if (zero_shq) then
            vmax = 0.0_dp
          else
            nbf = gdat%nbf
            if (gdat%direct_pure) then
              pints(1:nbf(4), 1:nbf(3), 1:nbf(2), 1:nbf(1)) => ints
            else
              pints(1:nbf(4), 1:nbf(3), 1:nbf(2), 1:nbf(1)) => ints
              call normalize_ints(nbf, gdat%am, pints)
              call int2_rys_reduce_pure(basis, gdat, ints, nbf)
            end if
            vmax = maxval(abs(ints(1:product(nbf))))
          end if
        end if
        schwarz_ints(ish, jsh) = sqrt(vmax)
        schwarz_ints(jsh, ish) = sqrt(vmax)
      end do
    end do

    if (libint2_active) then
        call libint2_cleanup_eri(erieval)
        deallocate(erieval)
    end if
    call gdat%clean()
  end subroutine ints_exchange

!###############################################################################

  subroutine int2_compute_data_t_storeints(consumer, basis, eri_data, &
                  buf, cutoff, nint)
    use precision, only: dp
    use basis_tools, only: basis_set

    implicit none

    class(int2_compute_data_t), intent(inout) :: consumer
    type(int2_storage_t), intent(inout) :: buf
    type(eri_data_t), target, intent(inout) :: eri_data

    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: cutoff
    integer, intent(inout) :: nint

    logical :: iandj, kandl, same
    integer :: itmp, nij, nkl, maxj, maxl, &
      i, i1, ii, loci, &
      j, j1, jj, locj, &
      k, k1, kk, lock, &
      l, l1, ll, locl
    integer :: ids(4), flips(4), nbf(4)
    real(kind=dp) :: val

    ids = eri_data%ids(eri_data%flips)

    nbf = eri_data%nbf

    same = ids(1)==ids(3) .and. ids(2)==ids(4)
    iandj = ids(1)==ids(2)
    kandl = ids(3)==ids(4)

    loci = basis%ao_offset(ids(1))-1
    locj = basis%ao_offset(ids(2))-1
    lock = basis%ao_offset(ids(3))-1
    locl = basis%ao_offset(ids(4))-1

    nij = 0
    maxj = nbf(2)
    do i = 1, nbf(1)
      if (iandj) maxj = i


jc:   do j = 1, maxj
        nij = nij+1

        nkl = nij

        maxl = nbf(4)
        do k = 1, nbf(3)
          if (kandl) maxl = k

          if (same) then ! account for non-unique permutations
             itmp = min(maxl,nkl)
             if (itmp==0) cycle jc
             maxl = itmp
             nkl = nkl-itmp
          end if

          do l = 1, maxl

            ! Element cutoff on the weighted magnitude: the stored value
            ! carries the orbit weight, so the effective threshold matches
            ! the C1 path even when non-abelian orbit members have unequal
            ! element magnitudes.
            val = eri_data%pints(l,k,j,i)
            if (eri_data%weighted_cutoff) then
              if (abs(val)*eri_data%weight < cutoff) cycle
            else
              ! Abelian tier: C1-identical element set (exact cancellation).
              if (abs(val) < cutoff) cycle
            end if
            val = val*eri_data%weight
            nint = nint + 1

            i1 = i+loci
            j1 = j+locj
            k1 = k+lock
            l1 = l+locl

            if (i1<j1) then ! sort <ij|
              j1 = i+loci
              i1 = j+locj
            end if

            if (k1<l1) then ! sort |kl>
              l1 = k+lock
              k1 = l+locl
            end if

            ii = i1
            jj = j1
            kk = k1
            ll = l1

            if (ii<kk) then ! sort <ij|kl>
              ii = k1
              jj = l1
              kk = i1
              ll = j1
            else if (ii==kk .and. jj<ll) then ! sort <ij|il>
              ii = i1
              jj = l1
              kk = k1
              ll = j1
            end if

!           Account for identical permutations.
            if (ii==jj) val = val*0.5d0
            if (kk==ll) val = val*0.5d0
            if (ii==kk .and. jj==ll) val = val*0.5d0

            buf%ncur = buf%ncur + 1
            buf%ids(1,buf%ncur) = int(ii, 2)
            buf%ids(2,buf%ncur) = int(jj, 2)
            buf%ids(3,buf%ncur) = int(kk, 2)
            buf%ids(4,buf%ncur) = int(ll, 2)
            buf%ints(buf%ncur) = val

            if (buf%ncur==buf%buf_size) call consumer%update(buf)
          end do
        end do
      end do jc
    end do
  end subroutine int2_compute_data_t_storeints

!###############################################################################

end module int2_compute
