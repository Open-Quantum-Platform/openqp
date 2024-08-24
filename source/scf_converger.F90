module scf_converger

  use precision, only: dp
  implicit none

  private
  public scf_conv
  public scf_conv_result

  integer, parameter, public :: conv_none = 1
  integer, parameter, public :: conv_cdiis = 2
  integer, parameter, public :: conv_ediis = 3
  integer, parameter, public :: conv_adiis = 4

  integer, parameter, public :: conv_dtype_fock    = 1
  integer, parameter, public :: conv_dtype_err     = 2
  integer, parameter, public :: conv_dtype_density = 3
  integer, parameter, public :: conv_dtype_energy  = 4

  integer, parameter, public :: conv_state_not_initialized = 0
  integer, parameter, public :: conv_state_initialized = 1

  integer, parameter, public :: conv_name_maxlen = 32


!> @brief Storage of Fock matrix convergence history for SCF convergers
!> @detail Is stores up to a `num_slots` of Fock, density, DIIS error matrices and SCF energies
!>  in cyclic buffers.
!>  Data can be accessed via `put` and `get` subroutines using the datatype id: conv_dtype_*
!>  Before putting new data, procedure `next_slot` should be called
  type :: converger_data
    integer :: slot = 0  !< current slot
    integer :: num_saved  !< amount of occupied slots
    integer :: num_slots !< total number of slots
    integer :: num_focks !< dimension of Fock and density matrices, = 1 for RHF, 2 for UHF
    integer :: ldim !< total size of the Fock matrix
    real(kind=dp), allocatable :: focks(:,:,:) !< array of Fock matrices
    real(kind=dp), allocatable :: errs(:,:,:) !< array of DIIS error matrices
    real(kind=dp), allocatable :: densities(:,:,:) !< array of density matrices
    real(kind=dp), allocatable :: energies(:) !< array of SCF energies
  contains
    procedure, private, pass :: init => converger_data_init
    procedure, private, pass :: clean => converger_data_clean
    procedure, private, pass :: next_slot => conv_data_next_slot
    procedure, private, pass :: discard_last => conv_data_discard
    procedure, private, pass :: conv_data_save_vec
    procedure, private, pass :: conv_data_save_scalar
    procedure, private, pass :: get => conv_data_get_val
    procedure, private, pass :: get_e => conv_data_get_energy

    generic, private :: put => conv_data_save_vec, conv_data_save_scalar
  end type

!> @brief Base type for SCF converger results
!> @detail It is used by main SCF convergence driver `scf_conv`
!>  The extending type should provide the following interfaces:
!>    init  : preliminary initialization of internal subconverger data
!>    clean : destructor
!>    setup : setting-up of the sub-converger equations, taking into the account the new data
!>            added in main SCF driver
!>    run   : solving the equations, returns `scf_conv_result` datatype. Multiple subsequent calls to
!>            this procedure without re-running `setup` should not change internal state and have to give same results
  type :: scf_conv_result
    integer :: ierr = 5
    real(kind=dp) :: error = 1.0e99_dp
    type(converger_data), pointer :: dat => null()
    character(len=conv_name_maxlen) :: active_converger_name = ''
  contains
    procedure, pass :: get_fock => conv_result_dummy_get_fock
    procedure, pass :: get_density => conv_result_dummy_get_density
  end type

!> @brief SCF converger results for interpolation methods like DIIS
!> @details The updated Fock/density will be computed as
!>   \f$ F_{n+1} = \sum_{i=1}^{n}{F_i c_i} } \f$
  type, extends(scf_conv_result) :: scf_conv_interp_result
    real(kind=dp), allocatable :: coeffs(:) !< interpolation coefficients \f$ c_i \f$
    integer :: pstart
    integer :: pend
  contains
    procedure, pass :: get_fock    => conv_result_interp_get_fock
    procedure, pass :: get_density => conv_result_interp_get_density
    procedure, private, pass :: interpolate => conv_result_interpolate
  end type

!> @brief Base type for real SCF convergers (subconvergers)
!> @detail It is used by main SCF convergence driver `scf_conv`
!>  The extending type should provide the following interfaces:
!>    init  : preliminary initialization of internal subconverger data
!>    clean : destructor
!>    setup : setting-up of the sub-converger equations, taking into the account the new data
!>            added in main SCF driver
!>    run   : solving the equations, returns `scf_conv_result` datatype. Multiple subsequent calls to
!>            this procedure without re-running `setup` should not change internal state and have to give same results
  type, abstract :: subconverger
    integer :: last_setup = 1024  !< Number of SCF iterations `setup` has been called
    integer :: iter  = 0 !< Number of iterations passed
    character(len=conv_name_maxlen) :: conv_name = '' !< Converger name
    type(converger_data), pointer :: dat !< Pointer to Fock, density, and DIIS error matrices data
  contains
    private
    procedure, pass :: subconverger_init
    procedure, pass :: subconverger_clean
    procedure, public, pass :: init => subconverger_init
    procedure, public, pass :: clean => subconverger_clean
    procedure(subconverger_run), pass, deferred :: run
    procedure(subconverger_setup), pass, deferred :: setup
  end type

!> @brief Container type for subconvergers
!>   Used to create allocatable array of allocatable types
  type :: subconverger_
    class(subconverger), allocatable :: s
  end type

!> @brief Main driver for converging SCF problems
!> @detail Manages and runs different real convergers depending on the current
!>  state of SCF optimization
  type :: scf_conv
    integer :: step  !< Current step
    real(kind=dp), pointer :: overlap(:,:) => null()  !< Pointer to the full-format overlap matrix (S)
    real(kind=dp), pointer :: overlap_sqrt(:,:) => null()  !< Pointer to the full-format S^(1/2) matrix
    type(converger_data) :: dat !< Storage of Fock, density, and DIIS error matrices for previous steps
    type(subconverger_), allocatable :: sconv(:) !< Array of converger methods
    real(kind=dp), allocatable :: thresholds(:) !< Thresholds to initiate convergers from `sconv`
                                                !< if `current_error < thresholds(i)`, then `sconv(i)` is used
    integer :: iter_space_size = 10  !< Default size of subconverger problem size
    integer :: verbose = 0  !< Verbosity
    integer :: state = 0 !< Current state of converger:
                         !< 0 - not initialized
                         !< 1 - initialized, normal
    real(kind=dp) :: current_error = 1.0e99_dp !< Maximum absolute value of current DIIS error matrix
  contains
    procedure, pass :: init => scf_conv_init
    procedure, pass :: clean => scf_conv_clean
    procedure, pass :: add_data => scf_conv_add_data
    procedure, private, pass :: select_method => scf_conv_select
    procedure, pass :: run => scf_conv_run
    procedure, private, pass :: compute_error => scf_conv_compute_error
  end type

!> @brief Dummy (steepest descent) converger
!> @detail Used for convenience and does nothing
  type, extends(subconverger) :: noconv_converger
  contains
    procedure, pass :: init => noconv_init
    procedure, pass :: run => noconv_run
    procedure, pass :: setup => noconv_setup
  end type

!> @brief Commutator DIIS (C-DIIS) converger
!> @detail Solves constraint minimization problem:
!>   \f$ min \{ Ax, \sum_i{x_i} = 1 \} \f$
!>   where \f$ A_{ij} = Tr([F_i D_i S_i],\ [F_j D_j S_j]) \f$
  type, extends(subconverger) :: cdiis_converger
    integer :: maxdiis
    real(kind=dp), allocatable :: a(:,:) !< C-DIIS A-matrix
    integer :: verbose = 0  !< Verbosity parameter
    integer :: old_dim = 0  !< Actual dimension of A matrix from previous step
  contains
    procedure, pass :: init => cdiis_init
    procedure, pass :: clean => cdiis_clean
    procedure, pass :: run => cdiis_run
    procedure, pass :: setup => cdiis_setup
  end type

!> @brief Datatype to pass optimization parameters to NLOpt
!>   for solving E/A-DIIS-like constraint minimization problems:
!>   \f$ min \{ x^T A x + bx,\ x_i \geq 0,\ \sum_i{x_i} = 1 \} \f$
  type :: ediis_opt_data
    real(kind=8), pointer :: b(:) => null()
    real(kind=8), pointer :: A(:,:) => null()
    procedure(eadiis_f), pointer, nopass :: fun => null()
  end type

!> @brief Energy DIIS (E-DIIS) converger
!> @detail Optimizes the following function:
!>  \f$ E_\mathrm{E-DIIS} = \sum_{i} c_i E_i - 0.5 \sum_{i,j} c_i c_j Tr( (D_i-D_j) (F_i-F_j) ) \f$
!> under the constraints:
!>  \f$ \sum_i {c_i} = 1, c_i \geq 0 \f$
  type, extends(cdiis_converger) :: ediis_converger
    real(kind=dp), allocatable :: b(:) !< energy history
    real(kind=dp), allocatable :: xlog(:,:) !< interpolation coefficients history
    type(ediis_opt_data) :: t
    procedure(eadiis_f), pointer, nopass :: fun => null()
  contains
    procedure, pass :: init => ediis_init
    procedure, pass :: clean => ediis_clean
    procedure, pass :: setup => ediis_setup
    procedure, pass :: run => ediis_run
  end type

!> @brief A-DIIS subconverger type
!> @detail Optimizes the following function:
!>  \f$ E_\mathrm{A-DIIS} = \sum_{i} c_i Tr((D_i-D_n)F_n) + 2 \sum_{i,j} c_i c_j Tr( (D_i-D_n) (F_j-F_n) ) \f$
!> under the constraints:
!>  \f$ \sum_i {c_i} = 1, c_i \geq 0 \f$
  type, extends(ediis_converger) :: adiis_converger
  contains
    procedure, pass :: init => adiis_init
    procedure, pass :: setup => adiis_setup
  end type

  abstract interface
    subroutine subconverger_run(self, res)
      import
      class(subconverger), target :: self
      class(scf_conv_result), allocatable, intent(out) :: res
    end subroutine

    subroutine subconverger_setup(self)
      import
      class(subconverger) :: self
    end subroutine

    subroutine eadiis_f(val, n, t, grad, need_gradient, d)
      import
      implicit none
      real(kind=8) :: val, t(*), grad(*)
      integer(kind=4), intent(in) :: n, need_gradient
      type(ediis_opt_data) :: d
    end subroutine
  end interface

!==============================================================================

contains

!==============================================================================

!> @brief Form the new Fock matrix
  subroutine conv_result_dummy_get_fock(self, matrix, istat)
    use precision, only: dp
    implicit none
    class(scf_conv_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)
    istat = 0

  end subroutine conv_result_dummy_get_fock

!==============================================================================

!> @brief Form the new density matrix
  subroutine conv_result_dummy_get_density(self, matrix, istat)
    use precision, only: dp
    implicit none
    class(scf_conv_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)
    istat = 0

  end subroutine conv_result_dummy_get_density

!==============================================================================

!> @brief Form the interpolated Fock matrix
!> @detail F_n = \sum_{n-maxdiis+1}^{n} F_i * C_i
  subroutine conv_result_interp_get_fock(self, matrix, istat)
    use precision, only: dp
    implicit none
    class(scf_conv_interp_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    if (self%ierr == 0) then
      call self%interpolate(matrix, conv_dtype_fock, istat)
    end if

  end subroutine conv_result_interp_get_fock

!==============================================================================

!> @brief Form the interpolated density matrix
!> @detail D_n = \sum_{n-maxdiis+1}^{n} D_i * C_i
  subroutine conv_result_interp_get_density(self, matrix, istat)
    use precision, only: dp
    implicit none
    class(scf_conv_interp_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    if (self%ierr == 0) then
      call self%interpolate(matrix, conv_dtype_density, istat)
    end if

  end subroutine conv_result_interp_get_density

!==============================================================================

!> @brief Form the interpolated matrix
!> @detail F_n = \sum_{n-maxdiis+1}^{n} F_i * C_i
  subroutine conv_result_interpolate(self, matrix, datatype, istat)
    use precision, only: dp

    implicit none

    class(scf_conv_interp_result), intent(in) :: self
    integer, intent(in) :: datatype
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)
    integer :: i, ifock

    if (self%ierr/=0) then
        istat = self%ierr
        return
    end if

    matrix = 0
    do ifock = 1, self%dat%num_focks
      do i = self%pstart, self%pend
        matrix(:,ifock) = matrix(:,ifock) &
                          + self%coeffs(i)*self%dat%get(i,datatype,ifock)
      end do
    end do
    istat = 0
  end subroutine conv_result_interpolate

!==============================================================================
!==============================================================================

!> @brief Finalize scf_conv datatype
  subroutine scf_conv_clean(self)
    class(scf_conv) :: self
    self%verbose = 0
    self%step = 0
    self%state = conv_state_not_initialized

    call self%dat%clean()

    if (allocated(self%thresholds)) deallocate(self%thresholds)
    if (allocated(self%sconv)) deallocate(self%sconv)
    nullify(self%overlap)
    nullify(self%overlap_sqrt)
  end subroutine

!==============================================================================

!> @brief Initializes the SCF converger driver
!> @param[in] ldim  Number of orbitals
!> @param[in] maxvec  Size of SCF converger linear space, e.g. number of DIIS vectors
!> @param[in] subconvergers  Array of SCF subconverger codes
!> @param[in] thresholds  Thresholds to initate subconverger
!>                         subconvergers[i] runs when current error is less than thresholds[i]
!> @param[in] overlap  overlap matrix (S) in full format
!> @param[in] overlap_sqrt  S^(1/2) matrix in full format
!> @param[in] num_focks  1 if R/ROHF, 2 if UHF
!> @param[in] verbose  sets verbosity
  subroutine scf_conv_init(self, ldim, maxvec, subconvergers, thresholds, &
                  overlap, overlap_sqrt, num_focks, verbose)
    class(scf_conv) :: self
    integer, intent(in) :: ldim
    integer, optional, intent(in) :: maxvec
    integer, optional, intent(in) :: subconvergers(:)
    real(kind=dp), optional, intent(in) :: thresholds(:)
    real(kind=dp), optional, target, intent(in) :: overlap(:,:), overlap_sqrt(:,:)
    integer, optional, intent(in) :: num_focks
    integer, optional, intent(in) :: verbose
    integer :: nfocks
    integer :: istat
    integer :: i

    if (self%state/=conv_state_not_initialized) call self%clean()

    if (present(thresholds)) then
      allocate(self%thresholds(0:ubound(thresholds,1)))
      self%thresholds(0:) = [thresholds, 0.0d0]
    else
      allocate(self%thresholds(0:1))
      self%thresholds(0:) = [1.0d0, 0.0d0]
    end if

    self%overlap => null()
    if (present(overlap)) then
        self%overlap => overlap
    end if

    self%overlap_sqrt => null()
    if (present(overlap_sqrt)) then
        self%overlap_sqrt => overlap_sqrt
    end if

    nfocks = 1
    if (present(num_focks)) nfocks = num_focks

    self%verbose = 0
    if (present(verbose)) then
        self%verbose = verbose
    end if

    self%iter_space_size = 15
    if (present(maxvec)) then
      self%iter_space_size = maxvec
    end if

    self%step = 0

    if (present(subconvergers)) then
        allocate(self%sconv(0:ubound(subconvergers,1)))
        allocate(noconv_converger :: self%sconv(0)%s)
        call self%sconv(0)%s%init(self)
        do i = 1, ubound(subconvergers,1)
          select case (subconvergers(i))
          case (conv_none)
            allocate(noconv_converger :: self%sconv(i)%s)
          case (conv_cdiis)
            allocate(cdiis_converger :: self%sconv(i)%s)
          case (conv_ediis)
            allocate(ediis_converger :: self%sconv(i)%s)
          case (conv_adiis)
            allocate(adiis_converger :: self%sconv(i)%s)
          end select
          call self%sconv(i)%s%init(self)
        end do
    end if

    call self%dat%init(ldim, nfocks, self%iter_space_size, istat)

    if (istat == 0) self%state = conv_state_initialized
  end subroutine

!==============================================================================

!> @brief Store the data from the new SCF iteration
!> @param[in] f  Fock matrix/matrices
!> @param[in] dens  Density matrix/matrices
!> @param[in] e  SCF energy
  subroutine  scf_conv_add_data(self, f, dens, e)
    use precision, only: dp

    implicit none

    class(scf_conv), intent(inout) :: self
!    real(kind=dp), intent(in) :: err(:,:)
    real(kind=dp), intent(in) :: f(:,:)
    real(kind=dp), intent(in) :: dens(:,:)
    real(kind=dp), optional, intent(in) :: e
    integer :: i

    call self%dat%next_slot()

!   Save the current Fock and density matrices
    call self%dat%put(f,    conv_dtype_fock)
    call self%dat%put(dens, conv_dtype_density)
    if (present(e)) call self%dat%put(e, conv_dtype_energy)

    self%current_error = self%compute_error()

    do i = lbound(self%sconv, 1), ubound(self%sconv, 1)
      self%sconv(i)%s%last_setup = self%sconv(i)%s%last_setup + 1
    end do

  end subroutine scf_conv_add_data

!==============================================================================

!> @brief Computes the new guess to the SCF wavefunction
!> @param[out] conv_result  results of the calculation
  subroutine scf_conv_run(self, conv_result)

    implicit none

    class(scf_conv), target :: self
    class(subconverger), pointer :: conv
    class(scf_conv_result), allocatable :: conv_result

    class(scf_conv_result), allocatable :: tmp_result

!    integer :: ierr

    conv => self%select_method(self%current_error)

    if (self%state == conv_state_not_initialized) then
        conv_result = scf_conv_result( &
                error = self%current_error &
            )
        return
    end if

!    if (conv%iter == 0) then
!      write (iw, "(10x, '* * *   Initiating ',A,' Converger   * * *')") &
!            trim(conv%conv_name)
!    end if

    conv%iter = conv%iter + 1
    self%step = self%step + 1

!!    call self%adiis_print()
!!   Update DIIS equation matrix
    call conv%setup()
!
!   Nothing left to do on the first iteration, exit
    if (self%step == 1) then
        conv_result = scf_conv_result( &
                ierr = 0, &
                active_converger_name = 'SD', &
                error = self%current_error &
            )
        return
    end if
!
!   Solve the set of DIIS linear equations
    call conv%run(tmp_result)

    tmp_result%error = self%current_error
!
!    if (self%verbose > 1) call self%print_solution

    call move_alloc(from=tmp_result, to=conv_result)

  end subroutine

!==============================================================================

!> @brief Select subconverger basing on the current DIIS error value
!> @param[in] error DIIS error value
!> @return  pointer to selected converger
  function scf_conv_select(self, error) result(conv)
    implicit none
    class(scf_conv), target :: self
    real(kind=dp), intent(in) :: error
    class(subconverger), pointer :: conv
    integer :: i, nconv

    nconv = ubound(self%thresholds, 1)
    do i = 0, nconv
      if (error > self%thresholds(i)) exit
    end do

    conv => self%sconv(min(i,nconv))%s

    ! Continue using the 'SD' converger if
    ! already initiated
    if (i == 0 .and. self%step > 0) then
        conv => self%sconv(1)%s
    end if

    ! Use SD by default if no other converger selected
    if (.not.associated(conv)) then
        conv => self%sconv(0)%s
    end if

  end function

!==============================================================================

!> @brief Calculate the DIIS error matrix: \f$ \mathrm{Err} = FDS - SDF \f$
!> @details This routine is general for RHF, ROHF, and UHF.
!>    Since each of `F`, `D`, `S` are symmetric, this means calculate \f$ FDS \f$,
!>    and then subtract the transpose from that result.
!>    Before entry, `F`, `D` and `S` must be expanded to square storage.
!> @return  DIIS error value
  function scf_conv_compute_error(self) result(diis_error)
    use mathlib, only: antisymmetrize_matrix
    use mathlib, only: unpack_matrix, pack_matrix
    use precision, only: dp
    use oqp_linalg

    implicit none

    class(scf_conv) :: self
    real(kind=dp) :: diis_error

!   all are (nbf, nbf) square matrices
    real(kind=dp), pointer, dimension(:) :: f, d, err
    real(kind=dp), allocatable :: fock_full(:,:), dens_full(:,:), err_full(:,:), wrk(:,:)
    integer :: nbf, ifock, nfocks

    nfocks = self%dat%num_focks

    nbf = self%dat%ldim

    allocate(fock_full(nbf,nbf))
    allocate(dens_full(nbf,nbf))
    allocate(err_full(nbf,nbf))
    allocate(wrk(nbf,nbf))

    diis_error = 0

    do ifock = 1, nfocks
        f => self%dat%get(-1,conv_dtype_fock,ifock)
        d => self%dat%get(-1,conv_dtype_density,ifock)
        call unpack_matrix(f,fock_full,nbf,'u')
        call unpack_matrix(d,dens_full,nbf,'u')
        ! F*D
        call dsymm('l','u',nbf,nbf,1.0d0, &
                fock_full,nbf,dens_full,nbf,0.0d0,wrk,nbf)
        ! (F*D)*S
        call dgemm('n','n',nbf,nbf,nbf,1.0d0, &
                wrk,nbf,self%overlap,nbf,0.0d0,err_full,nbf)
        ! F*D*S - S*D*F
        call antisymmetrize_matrix(err_full, nbf)

!       MV: This step is not really necessary
!       Put error matrix into consistent orthonormal basis
!       Pulay uses S**-1/2, but here we use Q, Q obeys Q-dagger*S*Q=I
!       E-orth = Q-dagger * E * Q, FCKA is used as a scratch `nbf` vector.
        call dgemm('t','n',nbf,nbf,nbf,1.0d0, &
                self%overlap_sqrt,nbf,err_full,nbf,0.0d0,wrk,nbf)
        call dgemm('n','n',nbf,nbf,nbf,1.0d0, &
                wrk,nbf,self%overlap_sqrt,nbf,0.0d0,err_full,nbf)

        err => self%dat%get(-1,conv_dtype_err,ifock)
        call pack_matrix(err_full, nbf, err, 'u')

!       Compute DIIS error (infinity norm of error matrix)
        diis_error = diis_error + maxval(abs(err))
    end do

    deallocate(fock_full)
    deallocate(dens_full)
    deallocate(err_full)
    deallocate(wrk)
  end function scf_conv_compute_error

!==============================================================================
!==============================================================================

!> @brief Allocate memory and initialize converger_data type
!> @param[in] ldim  number of orbitals
!> @param[in] nfocks  1 if R/ROHF, 2 if UHF
!> @param[in] nslots  max. number of stored datasets from previous SCF iterations
!> @param[out] istat  success status, nonzero if error occured
  subroutine converger_data_init(self, ldim, nfocks, nslots, istat)
    class(converger_data) :: self
    integer, intent(in) :: ldim, nfocks, nslots
    integer, intent(out) :: istat
    integer :: nbf_tri
    if (allocated(self%focks)) call self%clean()
    nbf_tri = (ldim+1)*ldim/2
    self%ldim = ldim
    self%num_focks = nfocks
    self%num_slots = nslots
    self%num_saved = 0
    self%slot = 0
    allocate( &
        self%focks(nbf_tri,nfocks,nslots), &
        self%densities(nbf_tri,nfocks,nslots), &
        self%energies(nslots), &
        self%errs(nbf_tri,nfocks,nslots), &
        stat=istat)
  end subroutine

!==============================================================================

!> @brief converger_data type finalization
  subroutine converger_data_clean(self)
    class(converger_data) :: self
    self%ldim = 0
    self%num_focks = 0
    self%num_slots = 0
    self%num_saved = 0
    if (allocated(self%focks)    ) deallocate(self%focks)
    if (allocated(self%densities)) deallocate(self%densities)
    if (allocated(self%energies) ) deallocate(self%energies)
    if (allocated(self%errs)     ) deallocate(self%errs)
  end subroutine

!==============================================================================

!> @brief Switch to the next slot in cyclic buffer
  subroutine conv_data_next_slot(self)
    class(converger_data) :: self
    self%slot = mod(self%slot, self%num_slots)+1
    self%num_saved = min(self%num_saved+1, self%num_slots)
  end subroutine

!==============================================================================

!> @brief Discard latest data
  subroutine conv_data_discard(self)
    class(converger_data) :: self
    self%slot = mod(self%slot-2, self%num_slots)+1
    self%num_saved = min(self%num_saved-1, 1)
  end subroutine

!==============================================================================

!> @brief Save a vector/matrix to the converger_data
!> @param[in] val data to be saved
!> @param[in] dtype the type of the data (fock/density/error)
  subroutine conv_data_save_vec(self, val, dtype)
    class(converger_data) :: self
    real(kind=dp) :: val(:,:)
    integer, intent(in) :: dtype
    integer :: im, nmatrices
    nmatrices = ubound(val,2)
    do im = 1, nmatrices
      select case(dtype)
      case (conv_dtype_fock)
          self%focks(:,im,self%slot) = val(:,im)
      case (conv_dtype_err)
          self%errs(:,im,self%slot) = val(:,im)
      case (conv_dtype_density)
          self%densities(:,im,self%slot) = val(:,im)
      case default
      end select
    end do
  end subroutine

!==============================================================================

!> @brief Save a double precision scalar to the converger_data
!> @param[in] val data to be saved
!> @param[in] dtype the type of the data
  subroutine conv_data_save_scalar(self, val, dtype)
    class(converger_data) :: self
    real(kind=dp) :: val
    integer, intent(in) :: dtype
    select case(dtype)
    case (conv_dtype_energy)
        self%energies(self%slot) = val
    case default
    end select
  end subroutine

!==============================================================================

!> @brief Get a vector/matrix from the converger_data
!> @param[in] n Data slot ID. The allowed values are:
!>              1 (oldest) to num_saved (latest)
!>              -1 is a shortcut to latest
!> @param[in] dtype the type of the data (fock/density/error)
!> @param[in] matrix_id 1 - alpha, 2 - beta Fock/density/error matrix
!> @result pointer to the data
  function conv_data_get_val(self, n, dtype, matrix_id) result(res)
    class(converger_data), target :: self
    real(kind=dp), pointer :: res(:)
    integer, intent(in) :: n, dtype, matrix_id
    integer :: slot, num_saved

    if (n == -1) then
!     Pick last one
      slot = self%slot
    else
!     Pick n-th starting from the oldest
      num_saved = self%num_saved
      slot = modulo(self%slot-num_saved+n-1, self%num_slots)+1
    end if
    select case(dtype)
    case (conv_dtype_fock)
      res => self%focks(:,matrix_id,slot)
    case (conv_dtype_err)
      res => self%errs(:,matrix_id,slot)
    case (conv_dtype_density)
      res => self%densities(:,matrix_id,slot)
    case default
    end select
  end function

!==============================================================================

!> @brief Get a double precision scalar from the converger_data
!> @param[in] n Data slot ID. The allowed values are:
!>              1 (oldest) to num_saved (latest)
!>              -1 is a shortcut to latest
!> @param[in] dtype the type of the data (fock/density/error)
!> @result stored value
  function conv_data_get_energy(self, n, dtype) result(res)
    class(converger_data), target :: self
    real(kind=dp) :: res
    integer, intent(in) :: n, dtype
    integer :: slot, num_saved

    if (n == -1) then
!     Pick last one
      slot = self%slot
    else
!     Pick n-th starting from the oldest
      num_saved = self%num_saved
      slot = modulo(self%slot-num_saved+n-1, self%num_slots)+1
    end if
    select case(dtype)
    case (conv_dtype_energy)
      res = self%energies(slot)
    case default
    end select
  end function

!==============================================================================

!> @brief Initialize subconverger
!> @detail This subroutine takes SCF converger driver as argument.
!>   It should be initialized and include all the required parameters.
!> @param[in] params  current SCF converger driver
  subroutine subconverger_init(self, params)
    implicit none
    class(subconverger) :: self
    type(scf_conv), target :: params
    self%iter = 0
    self%last_setup = 1024
    self%dat => params%dat
  end subroutine

!==============================================================================

!> @brief Finalize subconverger
  subroutine subconverger_clean(self)
    class(subconverger) :: self
    self%iter = 0
    self%last_setup = 1024
    self%dat => null()
  end subroutine

!==============================================================================

!> @brief Initialize SD subconverger
!> @param[in] params current SCF converger driver
  subroutine noconv_init(self, params)
    implicit none
    class(noconv_converger) :: self
    type(scf_conv), target :: params

    if (self%iter>0) call self%clean()

    call self%subconverger_init(params)
    self%conv_name = 'SD'

  end subroutine

!> @brief Computes the new guess to the SCF wavefunction
!> @param[out] res  results of the calculation
  subroutine noconv_run(self, res)
    class(noconv_converger), target :: self
    class(scf_conv_result), allocatable, intent(out) :: res
    res = scf_conv_result(ierr = 0, active_converger_name = 'SD')
  end subroutine

!==============================================================================

!> @brief Prepare subconverger to run
  subroutine noconv_setup(self)
    class(noconv_converger) :: self
  end subroutine

!==============================================================================

!> @brief Initialize C-DIIS subconverger
!> @detail This subroutine takes SCF converger driver as argument.
!>   It should be initialized and include all the required parameters.
!> @param[in] params current SCF converger driver
  subroutine cdiis_init(self, params)
    implicit none
    class(cdiis_converger) :: self
    type(scf_conv), target :: params

    if (self%iter>0) call self%clean()

    call self%subconverger_init(params)

    self%conv_name = 'C-DIIS'
    self%verbose = params%verbose
    self%maxdiis = params%iter_space_size
    self%dat => params%dat

    allocate(self%a(self%maxdiis,self%maxdiis), source=0.0d0)
  end subroutine

!==============================================================================

!> @brief Finalize C-DIIS subconverger
  subroutine cdiis_clean(self)
    class(cdiis_converger) :: self
    call self%subconverger_clean()

    self%verbose = 0
    if (allocated(self%a)) deallocate(self%a)
  end subroutine

!==============================================================================

!> @brief Computes the new guess to the SCF wavefunction using C-DIIS
!> @param[out] res  results of the calculation
  subroutine cdiis_run(self, res)
    use mathlib, only: solve_linear_equations
    use io_constants, only: iw

    implicit none

    class(cdiis_converger), target :: self
    class(scf_conv_result), allocatable, intent(out) :: res

    integer :: i, na, cur, info
    real(kind=dp) :: a_loc(self%maxdiis+1, self%maxdiis+1)
    real(kind=dp) :: x_loc(self%maxdiis+1)
    real(kind=dp), allocatable :: x(:)

    allocate(scf_conv_interp_result :: res)
    res%dat => self%dat
    res%active_converger_name = self%conv_name
    select type (res)
    class is (scf_conv_interp_result)
        res%pstart = 1
        res%pend = self%dat%num_saved
    end select

    res%ierr = 3  ! need to set up DIIS equations first
    if (.not.self%last_setup == 0) return

    !na = min(self%iter, self%maxdiis)
    na = self%dat%num_saved

!   Solve the set of DIIS linear equations
    a_loc(:self%maxdiis,:self%maxdiis) = self%a
    a_loc(:,na+1) = -1
    do i = na, 1, -1
!     Helper index, needed for dimension reduction in case of instability
      cur = na-i+1

      x_loc = 0.0d0
      x_loc(na+1) = -1.0d0
      info= 0
      call solve_linear_equations(a_loc(cur, cur), x_loc(cur:), i+1, 1, self%maxdiis+1, info)
      if (info <= 0) exit
      write (iw, *) 'Reducing DIIS Equation size by 1 for numerical stability'
    end do

    if (info<0) then ! illegal value of DSYSV argument
        res%ierr = 2
    else if (info > 0) then
        res%ierr = 1 ! singular DIIS matrix
    else
        res%ierr = 0 ! normal exit
        x = x_loc(1:self%maxdiis)
        select type (res)
        class is (scf_conv_interp_result)
            call move_alloc(from=x, to=res%coeffs)
        end select
    end if

  end subroutine

!==============================================================================

!> @brief Prepare C-DIIS subconverger to run
  subroutine cdiis_setup(self)
    use precision, only: dp

    implicit none

    class(cdiis_converger) :: self
    integer :: i, j, na, maxdiis, ifock, nfocks
    real(kind=dp) :: factor

    maxdiis = self%maxdiis
    !na = min(self%iter, maxdiis)
    na = self%dat%num_saved
    nfocks = self%dat%num_focks

!   Factor to account RHF/UHF cases
    factor = 1.0d0/nfocks

    if (self%last_setup>1) then
!     DIIS matrix is rather old, generate it from scratch
      self%old_dim = na
      self%a = 0
      do ifock = 1, nfocks
        do i = 1, na
        do j = 1, i
          self%a(j,i) = self%a(j,i) + factor*dot_product(&
                            self%dat%get(j,conv_dtype_err,ifock), &
                            self%dat%get(i,conv_dtype_err,ifock))
        end do
        end do
      end do

    else if (self%last_setup == 1) then
!     DIIS matrix is old by 1 iteration, just update it

!     If the current number of iterations exceeds the dimension of A matrix:
!     discard the data of oldest iteration by shifting the bottom-rigth square
!     to the top-left corner
      if (self%old_dim >= maxdiis) then
        self%a(1:maxdiis-1,1:maxdiis-1) = self%a(2:maxdiis,2:maxdiis)
      end if

      self%old_dim = na
!     Compute new elements (`na`-th column)
      self%a(:,na) = 0
      do ifock = 1, nfocks
        do i = 1, na
          self%a(i,na) = self%a(i,na) + factor*dot_product(&
                            self%dat%get(na,conv_dtype_err,ifock), &
                            self%dat%get(i,conv_dtype_err,ifock))
        end do
      end do
    else if (self%last_setup == 0) then
!       DIIS matrix is already prepared nothing to do here
        continue
    else
!       Should not get here
        continue
    end if

!   Set last column to 1 for the Lagrange multiplier part
!    self%a(:na+1,na+1) = -1

!   DIIS matrix is already prepared nothing to do here
    self%last_setup = 0

  end subroutine

!==============================================================================
!==============================================================================

!> @brief Initialize E-DIIS subconverger
!> @detail This subroutine takes SCF converger driver as argument.
!>   It should be initialized and include all the required parameters.
!> @param[in] params current SCF converger driver
  subroutine ediis_init(self, params)
    class(ediis_converger) :: self
    type(scf_conv), target :: params

    call self%cdiis_converger%init(params)
    self%conv_name = 'E-DIIS'
    self%fun => ediis_fun
    allocate(self%b(self%maxdiis))
    allocate(self%xlog(self%maxdiis,self%maxdiis), source=0.0d0)
  end subroutine

!> @brief Finalize E-DIIS subconverger
  subroutine ediis_clean(self)
    class(ediis_converger) :: self
    call self%cdiis_converger%clean()

    if (allocated(self%b)) deallocate(self%b)
    if (allocated(self%xlog)) deallocate(self%xlog)
  end subroutine

!> @brief Prepare E-DIIS subconverger to run
  subroutine ediis_setup(self)
    use precision, only: dp

    implicit none

    class(ediis_converger) :: self
    integer :: i, j, na, maxdiis, ifock, nfocks
    real(kind=dp) :: factor

    maxdiis = self%maxdiis
    !na = min(self%iter, maxdiis)
    na = self%dat%num_saved
    nfocks = self%dat%num_focks

!   Factor to account RHF/UHF cases
    factor = 1.0d0/nfocks

    if (self%last_setup>1) then
!     DIIS matrix is rather old, generate it from scratch
      self%old_dim = na
      self%a = 0
      do ifock = 1, nfocks
        do i = 1, na
        do j = 1, i-1
          self%a(j,i) = self%a(j,i) + factor*dot_product(&
                    self%dat%get(j,conv_dtype_density,ifock) &
                  - self%dat%get(i,conv_dtype_density,ifock), &
                    self%dat%get(j,conv_dtype_fock,   ifock) &
                  - self%dat%get(i,conv_dtype_fock,   ifock))
          self%a(i,j) = self%a(j,i)
        end do
        end do
      end do

    else if (self%last_setup == 1) then
!     DIIS matrix is old by 1 iteration, just update it

!     If the current number of iterations exceeds the dimension of A matrix:
!     discard the data of oldest iteration by shifting the bottom-rigth square
!     to the top-left corner
      if (self%old_dim >= maxdiis) then
        self%a(1:maxdiis-1,1:maxdiis-1) = self%a(2:maxdiis,2:maxdiis)
      end if

      self%old_dim = na
!     Compute new elements (`na`-th column)
      self%a(:,na) = 0
      do ifock = 1, nfocks
        do j = 1, na
          self%a(j,na) = self%a(j,na) + factor*dot_product(&
                    self%dat%get(j, conv_dtype_density,ifock) &
                  - self%dat%get(na,conv_dtype_density,ifock), &
                    self%dat%get(j, conv_dtype_fock,   ifock) &
                  - self%dat%get(na,conv_dtype_fock,   ifock))
        end do
      end do
      self%a(na,1:na-1) = self%a(1:na-1,na)

    else if (self%last_setup == 0) then
!       DIIS matrix is already prepared nothing to do here
        continue
    else
!       Should not get here
        continue
    end if

    do i = 1, na
      self%b(i) = self%dat%get_e(i,conv_dtype_energy)
    end do

!   DIIS matrix is already prepared nothing to do here
    self%last_setup = 0

  end subroutine

!==============================================================================

!> @brief Computes the new guess to the SCF wavefunction using E-DIIS
!> @param[out] res  results of the calculation
  subroutine ediis_run(self, res)
    use io_constants, only: iw
    use nlopt

    implicit none

    class(ediis_converger), target :: self
    class(scf_conv_result), allocatable, intent(out) :: res

    real(kind=dp), parameter :: tol = 1.0d-5
    real(kind=dp), parameter :: constrtol = 1.0d-8
    real(kind=dp) :: minf, minf_min
    real(kind=dp), allocatable :: x(:), xmin(:)
    integer(kind=4) :: ires
    integer :: na
    integer :: i, j
    logical :: is_a_repeat
    integer :: opt_global, opt_lbfgs
    type(ediis_opt_data) :: t

    allocate(scf_conv_interp_result :: res)
    res%dat => self%dat
    res%active_converger_name = self%conv_name
    select type (res)
    class is (scf_conv_interp_result)
        res%pstart = 1
        res%pend = self%dat%num_saved
    end select

    res%ierr = 3  ! need to set up DIIS equations first
    if (.not.self%last_setup == 0) return

    !na = min(self%iter, self%maxdiis)
    na = self%dat%num_saved

    if (self%iter > self%maxdiis) then
      do i = 1, self%maxdiis-1
        self%xlog(:,i) = cshift(self%xlog(:,i+1), 1)
        self%xlog(i+1:,i) = 0
      end do
    end if

    allocate(x(na), xmin(na))

    opt_global = 0
    minf_min = huge(1.0d0)

!   Initialize E-DIIS equation parameters for NLOpt
    t = ediis_opt_data(A=self%a(1:na,1:na), b=self%b(1:na), fun=self%fun)
!    call nlo_optimize(ires, opt_lbfgs, x(:na), minf)
!    call nlo_destroy(opt_lbfgs)

!   Set up the Improved Stochastic Ranking Evolution Strategy
!   It will run coarse global optimization, which will be further refined via L-BFGS
    call nlo_create(opt_global, NLOPT_GN_ISRES, na)
!   Max. number of calls to the objective function
    call nlo_set_maxeval(ires, opt_global, 100*(na+1));
!   Relative convergence tolerance for arguments
    call nlo_set_xtol_rel(ires, opt_global, 1.0d-2)
!   Absolute convergence tolerance for function value
    call nlo_set_ftol_abs(ires, opt_global, 1.0d-4)
!   Relative convergence tolerance for function value
    call nlo_set_ftol_rel(ires, opt_global, 1.0d-4)

!   Sum of coeffs equal to 1
    call nlo_add_equality_constraint(ires, opt_global, eadiis_constraints, 0, constrtol)
!   0 <= c_i <= 1
    call nlo_set_lower_bounds1(ires, opt_global, 0.0d0)
    call nlo_set_upper_bounds1(ires, opt_global, 1.0d0)
!   Objective function
    call nlo_set_min_objective(ires, opt_global, eadiis_fun, t)

    x = 0.0d0
    x(:na) = 1.0d0/na
    call nlo_optimize(ires, opt_global, x(:na), minf)
    call nlo_destroy(opt_global)

!   Refine the results of global optimization via L-BFGS
    opt_lbfgs = 0
    call nlo_create(opt_lbfgs, NLOPT_LD_LBFGS, na)
    call nlo_set_xtol_rel(ires, opt_lbfgs, tol)
    call nlo_set_ftol_abs(ires, opt_lbfgs, tol*tol)
    call nlo_set_ftol_rel(ires, opt_lbfgs, tol*tol)

!   Here, the modified E-DIIS equations are used, because L-BFGS does not support
!   equality constraints
!   They utilize the following variable substitution:
!   c_i = t_i^2/\sum_i{t_i^2}
    call nlo_set_min_objective(ires, opt_lbfgs, eadiis_objective, t)

    call nlo_optimize(ires, opt_lbfgs, x(:na), minf)

!   Because we used modified E-DIIS equations, we need to compute
!   coefficients `c` from `t`:
    x = x**2/sum(x**2)

!   Get prediction of the new SCF energy
    call eadiis_fun(minf, int(na, 4), x, x, int(0,4), t)

    if (ires < 0) then
      if (self%verbose > 2) then
        write(iw, '(10X,"*** nlopt0 failed:",I4," ***")') ires
      end if
    elseif (minf < minf_min) then
        is_a_repeat = &
          any( [( norm2(self%xlog(:na,j)-x(:na))<1.0d-4, &
                    j = 1, min(self%iter, self%maxdiis) )] )&
          .or.any(1.0d0-x(1:na-1) < 1.0d-4)
        if (.not.is_a_repeat) then
            minf_min = minf
            xmin = 0
            xmin = x
            if (self%verbose > 2) then
              write(iw,'(A,*(F15.6))') 'nlopt0: improving x at ', xmin(:)
              write(iw,'(A,*(F15.6))') 'nlopt0: improved val = ', minf
            end if
        elseif (self%verbose > 2) then
            write(iw,'(A,*(F15.6))') 'nlopt0: found rep at ', x(:)
            write(iw,'(A,*(F15.6))') 'nlopt0: rep val = ', minf
        end if
    else
        if (self%verbose > 2) then
          write(iw,'(A,*(F15.6))') 'nlopt0: found min at ', x(:)
          write(iw,'(A,*(F15.6))') 'nlopt0: min val = ', minf
        end if
    end if

!   If no solution found, try the alternative:
!   Start from the trivial guess [x(1:n-1)=0, x(n) = 1]
!   then run two L-BFGS iterations and average with
!   [x(1:i-i), x(i+1:n) = 0, x(i) = 1] vector and run few L-BFGS steps again
!   for all [ i = n-1, 1 ] and then [i = 1, n]
    if ( minf_min > 1d99) then
      call nlo_set_maxeval(ires, opt_lbfgs, 2);
      call nlo_set_xtol_rel(ires, opt_lbfgs, 0.1d0)

      xmin = 0
      xmin(na) = 1.0d0
      do i = na-1, 1, -1
        x = 0
        x(i) = 1
        xmin = (xmin + x) / (1+sum(x))
        call nlo_optimize(ires, opt_lbfgs, xmin, minf)
        if (ires < 0) exit
        xmin = xmin**2/sum(xmin**2)
        if (ires < 0) then
          if (self%verbose > 2) then
            write(iw, '(10X,"*** nlopt2 failed:",I4," ***")') ires
          end if
          exit
        end if
      end do
      if (ires >= 0) then
        do i = 1, na
          x = 0
          x(i) = 1
          xmin = (xmin + x) / (1+sum(x))
          call nlo_optimize(ires, opt_lbfgs, xmin, minf)
          if (ires < 0) then
            if (self%verbose > 2) then
              write(iw, '(10X,"*** nlopt2 failed:",I4," ***")') ires
            end if
            exit
          end if
          xmin = xmin**2/sum(xmin**2)
        end do
      end if
!     If still no success, the default is minimum energy + small contribution from others:
      if (ires < 0) then
        xmin = 1
        xmin(minloc(self%b(1:na-1))) = 10
        xmin = xmin/sum(xmin)
        if (self%verbose > 2) then
          write(iw, *) 'nlopt2: unoptimal default'
        end if
      end if
      call eadiis_fun(minf_min, int(na, 4), xmin, xmin, int(0,4), t)

    end if

    minf = minf_min

    self%xlog(:na,min(self%iter, self%maxdiis)) = xmin(1:na)

    call nlo_destroy(opt_lbfgs)

    res%ierr = 0
    select type (res)
    class is (scf_conv_interp_result)
        call move_alloc(from=xmin, to=res%coeffs)
    end select

  end subroutine

!==============================================================================

!> @brief Modified E/A-DIIS objective function wrapper, which allows to use unconstrained optimization
!> @details The following variable substitution is used:
!>  \f$ c_i = t_i^2 / \sum_i{t_i^2} \f$
!> @note This is standard interface to work with NLOpt library
!> @param[out] val function value
!> @param[in] n dimension of the problem
!> @param[in] t vector of arguments
!> @param[out] grad vector of function gradient
!> @param[in] need_gradient flag to turn on computing gradient, 0 - gradient not computed
!> @param[in] d datatype storing function parameters
  subroutine eadiis_objective(val, n, t, grad, need_gradient, d)
    implicit none
    integer(kind=4), intent(in) :: n, need_gradient
    real(kind=8) :: val, t(n), grad(n)
    type(ediis_opt_data) :: d
    real(kind=8) :: x(n), tnorm, jac(n,n)
    integer :: i, j
    tnorm = 1/sum(t(1:n)**2)
    x = (t(1:n)**2)*tnorm

    call d%fun(val, n, x, grad, need_gradient, d)

    jac = 0
    if (need_gradient/=0) then
      do i = 1, n
        jac(i,i) = 1
        do j = 1, n
          jac(j,i) = 2*tnorm*t(j)*(jac(j,i)-x(i))
        end do
      end do
      grad(1:n) = matmul(jac, grad(1:n))
    end if
  end subroutine

!==============================================================================

!> @brief Non-modified E/A-DIIS objective function wrapper
!> @note This is standard interface to work with NLOpt library
!> @param[out] val function value
!> @param[in] n dimension of the problem
!> @param[in] t vector of arguments
!> @param[out] grad vector of function gradient
!> @param[in] need_gradient flag to turn on computing gradient, 0 - gradient not computed
!> @param[in] d datatype storing function parameters
  subroutine eadiis_fun(val, n, x, grad, need_gradient, d)
    implicit none
    integer(kind=4), intent(in) :: n, need_gradient
    real(kind=8) :: val, x(n), grad(n)
    type(ediis_opt_data) :: d

    call d%fun(val, n, x, grad, need_gradient, d)
  end subroutine

!==============================================================================

!> @brief E-DIIS objective function calculation
!> @note This is standard interface to work with NLOpt library
!> @param[out] val function value
!> @param[in] n dimension of the problem
!> @param[in] t vector of arguments
!> @param[out] grad vector of function gradient
!> @param[in] need_gradient flag to turn on computing gradient, 0 - gradient not computed
!> @param[in] d datatype storing function parameters
  subroutine ediis_fun(val, n, x, grad, need_gradient, d)
    implicit none
    real(kind=8) :: val, x(*), grad(*)
    integer(kind=4), intent(in) :: n, need_gradient
    type(ediis_opt_data) :: d

    if (need_gradient/=0) then
       grad(1:n) = d%b(1:n) - matmul(d%a(1:n,1:n),x(1:n))
    end if
    val = dot_product(x(1:n),d%b(1:n)) &
            - 0.5*dot_product(x(1:n), matmul(d%a(1:n,1:n),x(1:n)))
  end subroutine

!==============================================================================

!> @brief E/A-DIIS constraints
!> @note This is standard interface to work with NLOpt library
!> @param[out] val function value
!> @param[in] n dimension of the problem
!> @param[in] t vector of arguments
!> @param[out] grad vector of function gradient
!> @param[in] need_gradient flag to turn on computing gradient, 0 - gradient not computed
!> @param[in] d datatype storing function parameters
  subroutine eadiis_constraints(val, n, x, grad, need_gradient, d)
    implicit none
    integer(kind=4) :: need_gradient
    integer(kind=4) :: n
    real(kind=8) :: val, x(n), grad(n)
    class(ediis_converger) :: d
    if (need_gradient.ne.0) then
      grad = 1
    end if
    val = sum(x) - 1
  end subroutine

!==============================================================================
!==============================================================================

!> @brief Initialize A-DIIS subconverger
!> @detail This subroutine takes SCF converger driver as argument.
!>   It should be initialized and include all the required parameters.
!> @param[in] params current SCF converger driver
  subroutine adiis_init(self, params)
    class(adiis_converger) :: self
    type(scf_conv), target :: params

    call self%cdiis_converger%init(params)
    self%conv_name = 'A-DIIS'
    self%fun => adiis_fun
    allocate(self%b(self%maxdiis))
    allocate(self%xlog(self%maxdiis,self%maxdiis), source=0.0d0)
  end subroutine

!> @brief Prepare A-DIIS subconverger to run
  subroutine adiis_setup(self)
    use precision, only: dp
!    use io_constants, only: iw

    implicit none

    class(adiis_converger) :: self
    integer :: i, j, na, maxdiis, ifock, nfocks
    real(kind=dp) :: factor

    maxdiis = self%maxdiis
    !na = min(self%iter, maxdiis)
    na = self%dat%num_saved
    nfocks = self%dat%num_focks

!   Factor to account RHF/UHF cases
    factor = 1.0d0/nfocks

    self%old_dim = na
    self%a = 0
    self%b = 0
    do ifock = 1, nfocks
      do i = 1, na
      do j = 1, na
        self%a(j,i) = self%a(j,i) + 0.5d0*factor*dot_product(&
            self%dat%get(j,conv_dtype_density,ifock) - self%dat%get(-1,conv_dtype_density,ifock), &
            self%dat%get(i,conv_dtype_fock,   ifock) - self%dat%get(-1,conv_dtype_fock,   ifock))
      end do
      end do
      do i = 1, na
        self%b(i) = self%b(i) + factor*dot_product(&
            self%dat%get(i,conv_dtype_density,ifock) &
          - self%dat%get(-1,conv_dtype_density,ifock), &
            self%dat%get(-1,conv_dtype_fock,ifock))
      end do
    end do


!   DIIS matrix is already prepared nothing to do here
    self%last_setup = 0

!    call diis_print_equation(self)

  end subroutine

!==============================================================================

!> @brief A-DIIS objective function calculation
!> @note This is standard interface to work with NLOpt library
!> @param[out] val function value
!> @param[in] n dimension of the problem
!> @param[in] t vector of arguments
!> @param[out] grad vector of function gradient
!> @param[in] need_gradient flag to turn on computing gradient, 0 - gradient not computed
!> @param[in] d datatype storing function parameters
  subroutine adiis_fun(val, n, x, grad, need_gradient, d)
    implicit none
    type(ediis_opt_data) :: d
    real(kind=8) :: val, x(*), grad(*)
    integer(kind=4), intent(in) :: n, need_gradient

    if (need_gradient/=0) then
       grad(1:n) = 2.0*d%b(1:n) + matmul(d%a(1:n,1:n),x(1:n)) &
                                + matmul(x(1:n),d%a(1:n,1:n))
    end if
    val = d%b(n) &!d%get_e(-1,conv_dtype_energy) &
          + 2.0*dot_product(x(1:n),d%b(1:n)) &
          + dot_product(x(1:n), matmul(d%a(1:n,1:n),x(1:n)))
  end subroutine

!==============================================================================
!==============================================================================

!> @brief Debug printing of the DIIS equation data
  subroutine diis_print_equation(self)
    use printing, only: print_square
    use io_constants, only: iw
    implicit none
    class(cdiis_converger) :: self
    integer :: num_saved
    !num_saved = min(self%iter, self%maxdiis)
    num_saved = self%dat%num_saved
    write(iw,'("(dbg) --------------------------------------------------")')
    write(iw,'("(dbg) DIIS iteration / max.dim. : ",G0," / ",G0)') self%iter, self%maxdiis
    write(iw,'("(dbg)",I4," Fock sets stored")') num_saved
    write(iw,'("(dbg) --------------------------------------------------")')
    write(iw,'("(dbg) Current DIIS equation matrix:")')
    call print_square(self%a, num_saved, num_saved, ubound(self%a,1), tag='(dbg)')
    write(iw,'("(dbg) --------------------------------------------------")')
    select type(self)
    class is (ediis_converger)
    write(iw,'("(dbg) Current DIIS equation vector:")')
    write(iw,'(*(ES15.7,","))') self%b(1:num_saved)
    end select
  end subroutine diis_print_equation

end module scf_converger
