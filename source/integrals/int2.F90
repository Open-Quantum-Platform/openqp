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

!###############################################################################

  type eri_data_t
    logical :: attenuated_ints = .false.
    integer :: ids(4)
    integer :: flips(4)
    integer :: am(4)
    integer :: nbf(4)
    real(kind=dp) :: mu2 = 1.0d99
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
    real(kind=dp) :: mu = 1.0d99

    type(par_env_t) :: pe
  contains

    private
!    procedure, pass :: storeints => int2_compute_data_t_storeints
    procedure, public, pass :: init => int2_compute_t_init
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

  end subroutine int2_compute_t_init

!###############################################################################

  subroutine int2_compute_t_set_screening(this)
    implicit none
    class(int2_compute_t), intent(inout) :: this
    call ints_exchange(this%basis, this%schwarz_ints_regular)
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
    use, intrinsic :: iso_c_binding, only: C_NULL_PTR, C_INT
    use types, only: information
    use constants, only: NUM_CART_BF
!$  use omp_lib

    implicit none

    class(int2_compute_t), target, intent(inout) :: this
    class(int2_compute_data_t), intent(inout) :: int2_consumer

    integer :: i, j, k, l, ij
    integer :: lmax
    integer :: nint
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

    nshell = this%basis%nshell


    ! preparations for screening
    if (this%schwarz) then
      if (this%attenuated) then
        if (.not.allocated(this%schwarz_ints_attenuated)) then
          allocate(this%schwarz_ints_attenuated, mold=this%schwarz_ints_regular)
          call ints_exchange(this%basis, this%schwarz_ints_attenuated, this%mu**2)
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
!$omp   i, j, k, l, ij, jork, &
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
    eri_data%mu2 = this%mu**2

    if (libint2_active) then
        allocate(eri_data%erieval(this%basis%mxcontr**4))
        call libint2_init_eri(eri_data%erieval, int(4, C_INT), C_NULL_PTR)
    end if
    call eri_data%gdat%init(lmax, this%cutoffs, ok)

    call int2_storage%init(this%buf_size)
    int2_storage%thread_id = ithread + 1

!$omp barrier
    if (this%pe%size>1) then
      ij= 0
    end if
    do i = this%basis%nshell, 1, -1
      do j = 1, i
        if (this%pe%size>1) then
          ij = ij +1
          if (mod(ij, this%pe%size) /= this%pe%rank) cycle
        end if

        if (this%schwarz) then
          test = int2_consumer%screen_ij(this%schwarz_ints, i, j)
          if (test < this%cutoffs%integral_cutoff) then
            nschwz = nschwz + i*(i-1)/2+j
            cycle
          end if
        end if

!$omp do schedule(dynamic,2)
        do k = 1, i

          jork = k
          if (i==k) jork = j

          do l = 1, jork

!$          if (oflag) tim1 = omp_get_wtime()
            if (this%schwarz) then
              test = int2_consumer%screen_ijkl(this%schwarz_ints, i, j, k, l)
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
!$omp end do nowait
      end do
    end do

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
    call int2_consumer%pe%init(this%pe%comm, this%pe%use_mpi)
    call int2_consumer%parallel_stop()

    this%skipped = nschwz

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
    use int2e_rotaxis, only: genr22
    use int2e_libint, only: libint_compute_eri, libint_print_eri
    use int2e_rys, only: int2_rys_compute, rys_print_eri
    use iso_c_binding, only: c_f_pointer
    implicit none
    type(basis_set), intent(in) :: basis
    type(int2_pair_storage), intent(in) :: ppairs
    type(int2_cutoffs_t), intent(in) :: cutoffs
    type(eri_data_t), intent(inout), target :: eri_data
    logical, intent(out) :: zero_shq
    logical :: rotspd, libint, rys
    integer :: nbf(4)
    logical, parameter :: dbg_output = .false.
!    logical, parameter :: dbg_output = .true.
    integer :: max_am
    logical :: err

    zero_shq = .false.

    eri_data%am = basis%am(eri_data%ids)
    max_am = maxval(eri_data%am)
    eri_data%nbf = (eri_data%am+1)*(eri_data%am+2)/2

    rotspd = max_am <= 2
    libint = .not.rotspd.and.libint2_active.and..not.eri_data%attenuated_ints
    rys = .not.rotspd.and..not.libint

    if (rotspd) then

      if (eri_data%attenuated_ints) then
        ! erf-attenuated integrals
        call genr22(basis, ppairs, eri_data%ints, eri_data%ids, eri_data%flips, cutoffs, eri_data%mu2)
      else
        ! regular integrals
        call genr22(basis, ppairs, eri_data%ints, eri_data%ids, eri_data%flips, cutoffs)
      end if

      eri_data%nbf = eri_data%nbf(eri_data%flips)
      eri_data%pints(1:eri_data%nbf(4), 1:eri_data%nbf(3), 1:eri_data%nbf(2), 1:eri_data%nbf(1)) => eri_data%ints

    else if (libint) then

      call libint_compute_eri(basis, ppairs, cutoffs, eri_data%ids, 0, eri_data%erieval, eri_data%flips, zero_shq)
      if (zero_shq) return
      nbf = eri_data%nbf(eri_data%flips)
      eri_data%nbf = nbf
      call c_f_pointer(eri_data%erieval(1)%targets(1), eri_data%pints, shape=nbf([4,3,2,1]))
      call normalize_ints(nbf, eri_data%gdat%am, eri_data%pints)

    else if (rys) then

      call eri_data%gdat%set_ids(basis, eri_data%ids)
      if (eri_data%attenuated_ints) then
        ! erf-attenuated integrals
        call int2_rys_compute(eri_data%ints, eri_data%gdat, ppairs, zero_shq, eri_data%mu2)
      else
        ! regular integrals
        call int2_rys_compute(eri_data%ints, eri_data%gdat, ppairs, zero_shq)
      end if

      if (zero_shq) return
      nbf = eri_data%gdat%nbf
      eri_data%nbf = nbf
      eri_data%flips = eri_data%gdat%flips
      eri_data%pints(1:nbf(4), 1:nbf(3), 1:nbf(2), 1:nbf(1)) => eri_data%ints
      call normalize_ints(nbf, eri_data%gdat%am, eri_data%pints)

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
    integer :: nsh

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

    if (this%cur_pass == 1) then
      if (allocated(this%f)) then
          if ( any((shape(this%f) - [this%fockdim, this%nfocks, nthreads])/=0) ) then
              deallocate(this%f)
          end if
      end if
      if (.not.allocated(this%f)) then
          allocate(this%f(this%fockdim, this%nfocks, nthreads), source=0.0d0)
      else
          this%f = 0
      end if
    end if

    call this%init_screen(basis)

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
    if (this%nthreads /= 1) then
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
    integer :: ifock, mythread

    xval1 = this%scale_exchange
    xval4 = 4 * this%scale_coulomb
    mythread = buf%thread_id

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

        this%f(ij,ifock,mythread) = this%f(ij,ifock,mythread) + val4*this%d(kl,ifock)
        this%f(kl,ifock,mythread) = this%f(kl,ifock,mythread) + val4*this%d(ij,ifock)
        this%f(ik,ifock,mythread) = this%f(ik,ifock,mythread) - val1*this%d(jl,ifock)
        this%f(jl,ifock,mythread) = this%f(jl,ifock,mythread) - val1*this%d(ik,ifock)
        this%f(il,ifock,mythread) = this%f(il,ifock,mythread) - val1*this%d(jk,ifock)
        this%f(jk,ifock,mythread) = this%f(jk,ifock,mythread) - val1*this%d(il,ifock)
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
    integer :: mythread

    xval2 = 2 * this%scale_exchange
    xval4 = 4 * this%scale_coulomb

    mythread = buf%thread_id

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
    end do

    buf%ncur = 0

  end subroutine int2_urohf_data_t_update

!###############################################################################

  subroutine ints_exchange(basis, schwarz_ints, mu2)
    use int2e_rotaxis, only: genr22
    use int2e_libint, only: libint2_init_eri, libint2_cleanup_eri
    use int2e_libint, only: libint_compute_eri, libint_print_eri
    use int2e_libint, only: libint_t, libint2_active
    use int2e_rys, only: int2_rys_compute
    use types, only: information
    use constants, only: NUM_CART_BF
    use, intrinsic :: iso_c_binding, only: C_NULL_PTR, C_INT,  c_f_pointer

    implicit none

    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(inout) :: schwarz_ints(:,:)
    real(kind=dp), optional, intent(in) :: mu2

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
    integer :: ok
    logical :: rotspd, libint, zero_shq, rys
    logical :: attenuated
    real(kind=dp) :: vmax
    real(kind=dp), allocatable, target :: ints(:)
    real(kind=dp), pointer :: pints(:,:,:,:)
    type(libint_t), allocatable :: erieval(:)
    type(int2_rys_data_t) :: gdat
    type(int2_cutoffs_t) :: cutoffs
    type(int2_pair_storage) :: ppairs

    attenuated = present(mu2)

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
        libint = .not.rotspd.and.libint2_active.and..not.attenuated
        rys = .not.rotspd.and..not.libint
        if (rotspd) then
          if (attenuated) then
            call genr22(basis, ppairs, ints, shell_ids, flips, cutoffs, mu2)
          else
            call genr22(basis, ppairs, ints, shell_ids, flips, cutoffs)
          end if
          nbf = am(flips)
          nbf = (nbf+1)*(nbf+2)/2
          vmax = maxval(abs(ints(1:product(nbf))))
        else if (libint) then
          nbf = am(flips)
          call libint_compute_eri(basis, ppairs, cutoffs, shell_ids, 0, erieval, flips, zero_shq)
          call c_f_pointer(erieval(1)%targets(1), pints, shape=nbf([4,3,2,1]))
          vmax = maxval(abs(pints))
        else if (rys) then
          call gdat%set_ids(basis, shell_ids)
          if (attenuated) then
            call int2_rys_compute(ints, gdat, ppairs, zero_shq, mu2)
          else
            call int2_rys_compute(ints, gdat, ppairs, zero_shq)
          end if
          nbf = gdat%nbf
          pints(1:nbf(4), 1:nbf(3), 1:nbf(2), 1:nbf(1)) => ints
          call normalize_ints(nbf, gdat%am, pints)
          vmax = maxval(abs(pints))
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

            val = eri_data%pints(l,k,j,i)

            if (abs(val)<cutoff) cycle
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
