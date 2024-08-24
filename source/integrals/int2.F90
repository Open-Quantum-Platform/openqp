#define UNUSED_DUMMY(x) if (.false.) then ; if (size(shape(x))<0) continue ; end if

module int2_compute
  use precision, only: dp
  use openqp_config, only: basis_max_contraction

  use int2e_libint, only: libint_t, libint2_active
  use int2e_rys, only: int2_rys_data_t
  use basis_tools, only: basis_set
  use atomic_structure_m, only: atomic_structure
  use int2_pairs, only: &
    int2_cutoffs_t, &
    int2_pair_storage
  use messages, only: show_message, WITH_ABORT

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
  public schwdn
  public schwd2n

!###############################################################################

!  The memory for temporary integral savings
!  nangm = 2L+1 for L angular momentum
  integer, parameter :: nangm = 28, maxg = nangm**4

  type eri_data_t
    integer :: len_int
    logical :: attenuated_ints = .false.
    real(kind=dp) :: mu2 = 1.0d99
    real(kind=dp), allocatable :: ghondo(:)
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
    real(kind=dp), allocatable :: dsh(:)
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
    real(kind=dp), contiguous, pointer :: xints(:)

    logical :: schwarz = .true.
    integer :: buf_size = 50000

    integer :: skipped = 0

    integer :: mx_int_dim = 1

    type(int2_cutoffs_t) :: cutoffs
    type(int2_pair_storage) :: ppairs

    logical :: attenuated = .false.
    real(kind=dp) :: mu = 1.0d99

  contains

    private
!    procedure, pass :: storeints => int2_compute_data_t_storeints
    procedure, public, pass :: init => int2_compute_t_init
    procedure, public, pass :: set_screening => int2_compute_t_set_screening
    procedure, public, pass :: clean => int2_compute_t_clean
    procedure, public, pass :: run => int2_run
    procedure, public, pass :: run_generic => int2_twoei
    procedure, public, pass :: run_cam => int2_run_cam

    procedure, pass :: exchng => int2_compute_t_exchng
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
    this%xints => null()
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
    integer :: lmax, mxlen
    integer(4) :: status

    this%basis => basis
    this%atoms => infos%atoms
    call tagarray_get_data(infos%dat, OQP_XINTS, this%xints, status)
    call check_status(status, module_name, subroutine_name, OQP_XINTS)
    this%skipped = 0

    call basis%init_shell_centers()

    lmax = maxval(basis%ktype) - 1
    lmax = min(lmax, 4)
    this%mx_int_dim = max(4, (lmax*lmax+3*lmax+2)/2)

    cutoff = infos%control%int2e_cutoff
    call this%cutoffs%set(cutoff, ec1*cutoff, ec2*cutoff, cx1*log(10.0d0))

    if (libint2_active) call libint_static_init

    call this%ppairs%alloc(basis, this%cutoffs)
    call this%ppairs%compute(basis, this%cutoffs)

  end subroutine int2_compute_t_init

!###############################################################################

  subroutine int2_compute_t_set_screening(this)
    implicit none
    class(int2_compute_t), intent(inout) :: this
    call this%exchng()
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
!$  use omp_lib

    implicit none

    class(int2_compute_t), intent(inout) :: this
    class(int2_compute_data_t), intent(inout) :: int2_consumer

    integer :: i, j, k, l
    integer :: lmax, nangm, ngth(4)
    integer :: nint
    real(kind=dp) :: test
    integer, parameter :: lmax_vals(0:*) = [1, 4, 6, 10, 15, 21, 28]

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

    omp = .false.
!$  omp = .true.

    nschwz = 0
    lmax = maxval(this%basis%ktype) - 1
    if (lmax < 0 .or. lmax > 6) call abort
    nangm = lmax_vals(lmax)
    ngth(4) = 1
    ngth(3) = ngth(4)*nangm
    ngth(2) = ngth(3)*nangm
    ngth(1) = ngth(2)*nangm

    if (oflag) write(*,dbgfmt1)

!$omp parallel &
!$omp   private( &
!$omp   i, j, k, l, jork, &
!$omp   tim0, tim1, tim2, tim3, tim4, ithread,   &
!$omp   test, &
!$omp   int2_storage, &
!$omp   eri_data, &
!$omp   nthreads, &
!$omp   zero_shq) &
!$omp   shared(int2_consumer) &
!$omp   reduction(+:nschwz, nint, thr_nshq)

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
    allocate(eri_data%ghondo(maxg), source=0.0d0)
    allocate(eri_data%ints(maxg), source=0.0d0)
    allocate(eri_data%gdat)
    eri_data%len_int = this%mx_int_dim

    eri_data%attenuated_ints = this%attenuated
    eri_data%mu2 = this%mu**2

    if (libint2_active) then
        allocate(eri_data%erieval(basis_max_contraction**4))
        call libint2_init_eri(eri_data%erieval, int(4, C_INT), C_NULL_PTR)
    end if
    call eri_data%gdat%init(lmax, this%cutoffs, ok)

    call int2_storage%init(this%buf_size)
    int2_storage%thread_id = ithread + 1

!$omp barrier

    do i = this%basis%nshell, 1, -1
      do j = 1, i

        if (this%schwarz) then
          test = int2_consumer%screen_ij(this%xints, i, j)
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
              test = int2_consumer%screen_ijkl(this%xints, i, j, k, l)
              if (test < this%cutoffs%integral_cutoff) then
                nschwz = nschwz+1
                cycle
              end if
            end if
            thr_nshq = thr_nshq + 1
!$          if (oflag) tim2 = tim2 + omp_get_wtime() - tim1

            call shellquartet(this%basis, this%ppairs, this%cutoffs, eri_data, i, j, k, l, zero_shq)
!$          if (oflag) tim3 = tim3 + omp_get_wtime() - tim1

            if (zero_shq) cycle

            call int2_compute_data_t_storeints(int2_consumer, &
                    this%basis, eri_data, int2_storage, i, j, k, l, ngth, this%cutoffs%integral_cutoff, nint)
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

    call int2_consumer%parallel_stop()

    this%skipped = nschwz

  end subroutine int2_twoei

!###############################################################################

  real(kind=dp) function schwdn(dsh, ish, jsh, ksh, lsh)

    implicit none

    real(kind=dp), intent(in) :: dsh(*)
    integer, intent(in) :: ish, jsh, ksh, lsh

    integer :: ij, ik, ill, kl, jk, jl

!      ----- FIND MAXIMUM DENSITY CONTRIBUTION TO THIS SHELL SET -----
!      -DSH- IS THE DENSITY MATRIX ALREADY COMPRESSED TO SHELLS

    ij = ish*(ish-1)/2+jsh
    ik = ish*(ish-1)/2+ksh
    ill = ish*(ish-1)/2+lsh
    kl = ksh*(ksh-1)/2+lsh
    jk = jsh*(jsh-1)/2+ksh
    jl = jsh*(jsh-1)/2+lsh
    if (jsh < ksh) jk = ksh*(ksh-1)/2+jsh
    if (jsh < lsh) jl = lsh*(lsh-1)/2+jsh
    schwdn = max(4*dsh(ij), 4*dsh(kl), &
                 dsh(jl), dsh(jk), dsh(ill), dsh(ik))
  end function schwdn

!###############################################################################

  real(kind=dp) function schwd2n(dsh, i, j, k, l)

    implicit none

    real(kind=dp), intent(in) :: dsh(:,:)
    integer, intent(in) :: i, j, k, l

    schwd2n = max(4*dsh(i,j), 4*dsh(k,l), dsh(j,l), dsh(j,k), dsh(i,l), dsh(i,k))
  end function schwd2n

!###############################################################################
!###############################################################################

  function int2_compute_data_t_screen_ij(this, xints, i, j) result(res)
    implicit none
    class(int2_compute_data_t), intent(in) :: this
    real(kind=dp), intent(in) :: xints(:)
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
    real(kind=dp), intent(in) :: xints(:)
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
    real(kind=dp), intent(in) :: xints(:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j
    integer :: ij
    res = 1
    ij = (i*i-i)/2 + j
    res = res * xints(ij) * this%max_den
  end function int2_fock_data_t_screen_ij

!###############################################################################

  function int2_fock_data_t_screen_ijkl(this, xints, i, j, k, l) result(res)
    implicit none
    class(int2_fock_data_t), intent(in) :: this
    real(kind=dp), intent(in) :: xints(:)
    real(kind=dp) :: res
    integer, intent(in) :: i, j, k, l
    integer :: ij, kl
    res = 1
    ij = (i*i-i)/2 + j
    kl = (k*k-k)/2 + l
    res = res * xints(ij)*xints(kl)
    if (allocated(this%dsh)) res = res * schwdn(this%dsh, i, j, k, l)
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
    real(kind=dp), intent(out) :: shell_density(:)   !< Shell density matrix
    real(kind=dp), intent(in) :: density_matrix(:,:) !< Density matrix in AO basis

    integer :: ij_index, num_electrons, shell_i, shell_j, i, j, ij
    integer :: mini, maxi, minj, maxj, fock_index
    real(kind=dp) :: dmax

    ! Initialize shell density to zero
    shell_density = 0

    ! Loop over each Fock matrix
    do fock_index = 1, ubound(density_matrix,2)
      ij_index = 0
      num_electrons = 0

      ! Loop over shells
      do shell_i = 1, basis%nshell
        mini = basis%kloc(shell_i)
        maxi = mini+basis%kmax(shell_i) - basis%kmin(shell_i)

        ! Loop over shells again for symmetry
        do shell_j = 1, shell_i
          minj = basis%kloc(shell_j)
          maxj = minj + basis%kmax(shell_j) - basis%kmin(shell_j)
          ij_index = ij_index+1
          dmax = 0.0d0

          ! Loop over basis functions within shells
          do i = mini, maxi
            if (shell_i == shell_j) maxj = i
            do j = minj, maxj
              ij = (i - num_electrons)*(i-num_electrons - 1)/2 + j - num_electrons
              dmax = max(abs(density_matrix(ij,fock_index)), dmax)
            end do
          end do

          shell_density(ij_index) = max(dmax, shell_density(ij_index))
        end do
      end do
    end do
  end subroutine shlden

!###############################################################################

  subroutine shellquartet(basis, ppairs, cutoffs, eri_data, ish, jsh, ksh, lsh, zero_shq)
    use io_constants, only: iw
    use int2e_rotaxis, only: genr22, idrotax
    use int2e_libint, only: libint_compute_eri, libint_print_eri, libint_to_ghondo
    use int2e_rys, only: int2_rys_compute, rys_print_eri, rys_to_ghondo
    implicit none
    type(basis_set), intent(in) :: basis
    type(int2_pair_storage), intent(in) :: ppairs
    type(int2_cutoffs_t), intent(in) :: cutoffs
    type(eri_data_t), intent(inout), target :: eri_data
    integer, intent(in) :: ish, jsh, ksh, lsh
    logical, intent(out) :: zero_shq
    integer :: &
      i, ii, ij, ijk, ijkl, &
      j, k, l, &
      maxi, maxj, maxk, maxl, &
      mini, minj, mink, minl
    logical :: rotspd, libint, rys
    logical :: iandj, kandl, same
    integer :: ihondo, jhondo, khondo, lhondo, &
               ijkrotax, ijn, ijrotax, jmax, ijklrotax, irotax, kln, lmax
    integer :: flips(4)
    integer :: shell_ids(4)
    integer :: nbf(4)
    logical, parameter :: dbg_output = .false.
!    logical, parameter :: dbg_output = .true.
    real(kind=dp), pointer :: ghondo_p(:,:,:,:)
    real(kind=dp), pointer :: ints_p(:,:,:,:)
    integer :: max_am
    real(kind=dp) :: grotspd(1296)

    zero_shq = .false.

    ghondo_p(1:eri_data%len_int, 1:eri_data%len_int, 1:eri_data%len_int, 1:eri_data%len_int) => eri_data%ghondo

!  Rotated axis code

    if (dbg_output) eri_data%ghondo = 0.0d0

    shell_ids = [ish, jsh, ksh, lsh]
    iandj = ish == jsh
    kandl = ksh == lsh
    same = ish == ksh .and. jsh == lsh

    mini = basis%kmin(ish)
    maxi = basis%kmax(ish)
    minj = basis%kmin(jsh)
    maxj = basis%kmax(jsh)
    mink = basis%kmin(ksh)
    maxk = basis%kmax(ksh)
    minl = basis%kmin(lsh)
    maxl = basis%kmax(lsh)

    max_am = maxval(basis%ktype([ish,jsh,ksh,lsh])) - 1

    rotspd = max_am <= 2
    libint = .not.rotspd.and.libint2_active
    rys = .not.rotspd.and..not.libint

    if (rotspd) then

      if (eri_data%attenuated_ints) then
        ! erf-attenuated integrals
        call genr22(basis, ppairs, grotspd, shell_ids, flips, eri_data%mu2)
      else
        ! regular integrals
        call genr22(basis, ppairs, grotspd, shell_ids, flips)
      end if


!  Save to output array with hondo indexing

      ijn = 0
      jmax = maxj
      do i = mini, maxi
        ihondo = (i-mini)+1
        irotax = idrotax(flips(1), i)+1
        if (iandj) jmax = i
        jc: do j = minj, jmax
          jhondo = (j-minj)+1
          ijrotax = idrotax(flips(2), j)+irotax
          ijn = ijn+1
          lmax = maxl
          kln = 0
          do k = mink, maxk
            khondo = (k-mink)+1
            ijkrotax = idrotax(flips(3), k)+ijrotax
            if (kandl) lmax = k
            do l = minl, lmax
              kln = kln+1
              if (same .and. kln > ijn) cycle jc
              lhondo = (l-minl)+1
              ijklrotax = idrotax(flips(4), l)+ijkrotax
              ghondo_p(lhondo,khondo,jhondo,ihondo) = grotspd(ijklrotax)
            end do
          end do
        end do jc
      end do

      if (dbg_output) then
        block
            integer :: len1, len2, len3, len4
            write (iw, '("shells", 4i5, " :", 4i5, " :", 4i5, " :", 4i5)') &
                    ish, jsh, ksh, lsh, &
                    basis%ktype(shell_ids)-1, &
                    basis%ktype(shell_ids(flips))-1, &
                    flips
            len1 = 1
            len2 = eri_data%len_int
            len3 = len2*len2
            len4 = len2*len2*len2
            ii = 1
            do i = mini, maxi
              ij = ii
              do j = minj, maxj
                ijk = ij
                do k = mink, maxk
                  ijkl = ijk
                  do l = minl, maxl
                    write (iw, '("elem (",2i3," |",2i3,") = ",es30.15)') &
                            i-mini+1, j-minj+1, &
                            k-mink+1, l-minl+1, &
                            eri_data%ghondo(ijkl)
                    ijkl = ijkl+len1
                  end do
                  ijk = ijk+len2
                end do
                ij = ij+len3
              end do
              ii = ii+len4
            end do
        end block
      end if

    else if (libint) then

      call libint_compute_eri(basis, ppairs, cutoffs, shell_ids, 0, eri_data%erieval, flips, zero_shq)
      if (zero_shq) return
      call libint_to_ghondo(basis, shell_ids, 0, eri_data%erieval, flips, ghondo_p)

      if (dbg_output) then
        write (*, '("shells", 4i5, " :", 4i5, " :", 4i5, " :", 4i5)') &
                ish, jsh, ksh, lsh, &
                basis%ktype(shell_ids)-1, &
                basis%ktype(shell_ids(flips))-1, &
                flips
        call libint_print_eri(basis, shell_ids, 0, eri_data%erieval, flips)
      end if

    else if (rys) then

      call eri_data%gdat%set_ids(basis, shell_ids)
      if (eri_data%attenuated_ints) then
        ! erf-attenuated integrals
        call int2_rys_compute(eri_data%ints, eri_data%gdat, ppairs, zero_shq, eri_data%mu2)
      else
        ! regular integrals
        call int2_rys_compute(eri_data%ints, eri_data%gdat, ppairs, zero_shq)
      end if
      if (zero_shq) return
      nbf = eri_data%gdat%nbf
      ints_p(1:nbf(4), 1:nbf(3), 1:nbf(2), 1:nbf(1)) => eri_data%ints
      call rys_to_ghondo(eri_data%gdat, ints_p, ghondo_p)

      if (dbg_output) then
        write (*, '("dbg: shells", 4i5, " :", 4i5, " :", 4i5, " :", 4i5)') &
                ish, jsh, ksh, lsh, &
                basis%ktype(shell_ids)-1, &
                basis%ktype(shell_ids(eri_data%gdat%flips))-1, &
                eri_data%gdat%flips
        call rys_print_eri(eri_data%gdat, ints_p)
      end if
    end if

  end subroutine shellquartet

!###############################################################################

  subroutine compare_ints(basis, ish, jsh, ksh, lsh, lint, g1, g2)
    use io_constants, only: iw
    implicit none
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: ish, jsh, ksh, lsh, lint
    real(kind=dp), intent(in) :: g1(:), g2(:)
    real(kind=dp), parameter :: tol = 1.0d-6
    integer :: ii, ij, ijk, ijkl
    integer :: i, j, k, l
    integer :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer :: jmax, lmax
    integer :: ijn, kln
    logical :: iandj, kandl, same
    integer :: len1, len2, len3, len4

    len1 = 1
    len2 = lint
    len3 = lint*lint
    len4 = lint*lint*lint

    mini = basis%kmin(ish)
    maxi = basis%kmax(ish)
    minj = basis%kmin(jsh)
    maxj = basis%kmax(jsh)
    mink = basis%kmin(ksh)
    maxk = basis%kmax(ksh)
    minl = basis%kmin(lsh)
    maxl = basis%kmax(lsh)

    iandj = ish == jsh
    kandl = ksh == lsh
    same = ish == ksh .and. jsh == lsh

    ijn = 0
    ii = 1
    do i = mini, maxi
      ij = ii
      jmax = maxj
      if (iandj) jmax = i
      do j = minj, jmax
        ijk = ij
        ijn = ijn+1
        kln = 0
        do k = mink, maxk
          ijkl = ijk
          lmax = maxl
          if (kandl) lmax = k
          do l = minl, lmax
            kln = kln+1
            if (same .and. kln > ijn) cycle
            if (abs(g1(ijkl)-g2(ijkl))>tol) then
              write (iw, '("shell <", 2i4, " |",2i4,">" x "elem (",2i3," |",2i3,") = ",2es30.15)') &
                    ish, jsh, ksh, lsh, &
                    i-mini+1, j-minj+1, &
                    k-mink+1, l-minl+1, &
                    g1(ijkl), g2(ijkl)
            end if
            ijkl = ijkl+len1
          end do
          ijk = ijk+len2
        end do
        ij = ij+len3
      end do
      ii = ii+len4
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
        if (size(this%dsh) /= nsh*(nsh+1)/2) deallocate(this%dsh)
    end if

    if (.not.allocated(this%dsh)) then
        allocate(this%dsh(nsh*(nsh+1)/2), source=0.0d0)
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
    if (this%cur_pass /= this%num_passes) return
    if (this%nthreads == 1) return
    this%f(:,:,lbound(this%f,2)) = sum(this%f, dim=size(shape(this%f)))
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

  subroutine int2_compute_t_exchng(this)
    use int2e_rotaxis, only: &
      genr22, idrotax
    use int2e_libint, ONLY: libint2_init_eri, libint2_cleanup_eri
    use int2e_libint, only: libint_compute_eri, libint_print_eri, libint_to_ghondo
    use int2e_libint, only: libint_t, libint2_active
    use int2e_rys, only: int2_rys_compute
    use types, only: information
    use, intrinsic :: iso_c_binding, only: C_NULL_PTR, C_INT

    implicit none

    class(int2_compute_t), intent(inout) :: this

    logical iandj

    integer :: nangm, ngth(4)

    real(kind=dp) :: grotspd(1296)

    real(kind=dp), parameter :: &
      ic_exchng  = 1.0d-15, &
      ei1_exchng = 1.0d-17, &
      ei2_exchng = 1.0d-17, &
      cux_exchng = 50.0

    integer, parameter :: lmax_vals(0:*) = [1, 4, 6, 10, 15, 21, 28]
    integer :: flips(4)
    integer :: shell_ids(4)
    integer :: nbf(4)
    integer :: lmax
    integer :: ijij, ish, jsh
    integer :: i, j, ijn, jmax, nn
    integer :: mini, maxi, minj
    integer :: max_am
    integer :: ok
    logical :: rotspd, libint, zero_shq, rys
    real(kind=dp) :: icsv, ei1sv, ei2sv, cuxsv
    real(kind=dp) :: val, vmax
    real(kind=dp), allocatable, target :: ghondo(:)
    real(kind=dp), pointer :: ghondo_p(:,:,:,:)
    type(libint_t), allocatable :: erieval(:)
    type(int2_rys_data_t) :: gdat

    lmax = maxval(this%basis%ktype) - 1
    if (lmax < 0 .or. lmax > 6) then
      print *, 'lmax=', lmax
      flush (6)
      call abort
    end if
    nangm = lmax_vals(lmax)
    ngth(4) = 1
    ngth(3) = ngth(4)*nangm
    ngth(2) = ngth(3)*nangm
    ngth(1) = ngth(2)*nangm


!   Set very tight cutoff
    call this%cutoffs%get(icsv, ei1sv, ei2sv, cuxsv)
    call this%cutoffs%set(ic_exchng, ei1_exchng, ei2_exchng, cux_exchng)

    if (libint2_active) then
        allocate(erieval(basis_max_contraction**4))
        call libint2_init_eri(erieval, int(4, C_INT), C_NULL_PTR)
    end if
    allocate(ghondo(maxg), source=0.0d0)
    call gdat%init(lmax, this%cutoffs, ok)

    ijij = 0
    do ish = 1, this%basis%nshell
      do jsh = 1, ish
        ijij = ijij+1
        shell_ids = [ish, jsh, ish, jsh]

        max_am = maxval(this%basis%ktype(shell_ids)) - 1

        rotspd = max_am <= 2
        libint = .not.rotspd.and.libint2_active
        rys = .not.rotspd.and..not.libint
        if (rotspd) then
          call genr22(this%basis, this%ppairs, grotspd, shell_ids, flips)

          vmax = 0.0d0
          mini = this%basis%kmin(ish)
          minj = this%basis%kmin(jsh)
          maxi = this%basis%kmax(ish)
          jmax = this%basis%kmax(jsh)
          iandj = ish == jsh
          ijn = 0
          do i = mini, maxi
            if (iandj) jmax = i
            do j = minj, jmax
              nn = sum(idrotax(flips(1:4:2), i))+ &
                   sum(idrotax(flips(2:4:2), j))+1
              val = grotspd(nn)
              if (val > vmax) vmax = val
            end do
          end do
        else if (libint) then
          ghondo_p(1:this%mx_int_dim, 1:this%mx_int_dim, 1:this%mx_int_dim, 1:this%mx_int_dim) => ghondo
          call libint_compute_eri(this%basis, this%ppairs, this%cutoffs, shell_ids, 0, erieval, flips, zero_shq)
          call libint_to_ghondo(this%basis, shell_ids, 0, erieval, flips, ghondo_p)
          vmax = maxval(abs(ghondo_p))
        else if (rys) then
          call gdat%set_ids(this%basis, shell_ids)
          call int2_rys_compute(ghondo, gdat, this%ppairs, zero_shq)
          nbf = gdat%nbf
          ghondo_p(1:nbf(4), 1:nbf(3), 1:nbf(2), 1:nbf(1)) => ghondo
          vmax = maxval(abs(ghondo_p))
        end if
        this%xints(ijij) = sqrt(vmax)
      end do
    end do

    if (libint2_active) then
        call libint2_cleanup_eri(erieval)
        deallocate(erieval)
    end if
    call gdat%clean()
    call this%cutoffs%set(icsv, ei1sv, ei2sv, cuxsv)
  end subroutine int2_compute_t_exchng

!###############################################################################

  subroutine int2_compute_data_t_storeints(consumer, basis, eri_data, &
                  buf, ish, jsh, ksh, lsh, strides, cutoff, nint)
    use precision, only: dp
    use basis_tools, only: basis_set

    implicit none

    class(int2_compute_data_t), intent(inout) :: consumer
    type(int2_storage_t), intent(inout) :: buf
    type(eri_data_t), intent(inout) :: eri_data

    type(basis_set), intent(in) :: basis
    real(kind=dp), intent(in) :: cutoff
    integer, intent(in) :: &
      ish, jsh, ksh, lsh, strides(4)
    integer, intent(inout) :: nint

    logical :: iandj, kandl, same
    integer :: itmp, nij, nkl, maxj2, maxl2, &
      i, i1, i2, ii, maxi, mini, loci, i_index, &
      j, j1, j2, jj, maxj, minj, locj, ij_index, &
      k, k1, k2, kk, maxk, mink, lock, ijk_index, &
      l, l1, l2, ll, maxl, minl, locl, ijkl_index
    integer :: istride, jstride, kstride, lstride
    real(kind=dp) :: val

    istride = strides(1)
    jstride = strides(2)
    kstride = strides(3)
    lstride = strides(4)

    same = ish==ksh .and. jsh==lsh
    iandj = ish==jsh
    kandl = ksh==lsh

    mini = basis%kmin(ish)
    minj = basis%kmin(jsh)
    mink = basis%kmin(ksh)
    minl = basis%kmin(lsh)
    maxi = basis%kmax(ish)
    maxj = basis%kmax(jsh)
    maxk = basis%kmax(ksh)
    maxl = basis%kmax(lsh)
    loci = basis%kloc(ish)-mini
    locj = basis%kloc(jsh)-minj
    lock = basis%kloc(ksh)-mink
    locl = basis%kloc(lsh)-minl

    nij = 0
    maxj2 = maxj
    i_index = 1
    do i = mini, maxi
      if (iandj) maxj2 = i

      i1 = i+loci

      ij_index = i_index
      i_index = i_index+istride

jc:   do j = minj, maxj2
        nij = nij+1
        maxl2 = maxl

        j1 = j+locj
        i2 = i1
        j2 = j1
        if (i1<j1) then ! sort <ij|
          i2 = j1
          j2 = i1
        end if

        ijk_index = ij_index
        ij_index = ij_index+jstride

        nkl = nij

        do k = mink, maxk
          if (kandl) maxl2 = k

          k1 = k+lock

          if (same) then ! account for non-unique permutations
             itmp = min(maxl2-minl+1,nkl)
             if (itmp==0) cycle jc
             maxl2 = minl+itmp-1
             nkl = nkl-itmp
          end if

          ijkl_index = ijk_index
          ijk_index = ijk_index+kstride

          do l = minl, maxl2

            val = eri_data%ghondo( ijkl_index )
            ijkl_index = ijkl_index+lstride

            if (abs(val)<cutoff) cycle
            nint = nint + 1

            l1 = l+locl
            k2 = k1
            l2 = l1

            if (k2<l2) then ! sort |kl>
              k2 = l1
              l2 = k1
            end if

            ii = i2
            jj = j2
            kk = k2
            ll = l2

            if (ii<kk) then ! sort <ij|kl>
              ii = k2
              jj = l2
              kk = i2
              ll = j2
            else if (ii==kk .and. jj<ll) then ! sort <ij|il>
              jj = l2
              ll = j2
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
