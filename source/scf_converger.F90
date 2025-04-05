module scf_converger

  use precision, only: dp
  use io_constants, only: iw
  use mathlib, only: traceprod_sym_packed, orb_to_dens
  use messages, only: show_message, with_abort
  use mathlib, only: antisymmetrize_matrix, unpack_matrix, pack_matrix
  implicit none

  private
  public :: scf_conv
  public :: scf_conv_result
  public :: soscf_converger ! Used to provide parametres for soscf

  !> Constants for converger types
  integer, parameter, public :: conv_none  = 1
  integer, parameter, public :: conv_cdiis = 2
  integer, parameter, public :: conv_ediis = 3
  integer, parameter, public :: conv_adiis = 4
  integer, parameter, public :: conv_soscf = 5

  !> Converger state constants
  integer, parameter, public :: conv_state_not_initialized = 0
  integer, parameter, public :: conv_state_initialized     = 1

  !> Maximum length for converger names
  integer, parameter, public :: conv_name_maxlen = 32

  !> @brief Type to encapsulate SCF iteration data
  !> @detail Stores Fock, density, DIIS error matrices, and SCF energy for a single SCF iteration.
  !>         Designed for efficient memory use and data locality in a ring buffer.
  type :: scf_data_t
    real(kind=dp), allocatable :: focks(:,:)      !< Fock matrices (nbf_tri, num_focks)
    real(kind=dp), allocatable :: densities(:,:)  !< Density matrices (nbf_tri, num_focks)
    real(kind=dp), allocatable :: errs(:,:)       !< DIIS error matrices (nbf_tri, num_focks)
    real(kind=dp), allocatable :: mo_a(:,:)       !< MO coefficients alpha (nbf, nbf)
    real(kind=dp), allocatable :: mo_b(:,:)       !< MO coefficients beta (nbf, nbf)
    real(kind=dp), allocatable :: mo_e_a(:)       !< MO energies alpha (nbf)
    real(kind=dp), allocatable :: mo_e_b(:)       !< MO energies beta (nbf)
    real(kind=dp), allocatable :: occ_a(:)        !< Alpha orbital occupations (nbf)
    real(kind=dp), allocatable :: occ_b(:)        !< Beta orbital occupations (nbf)
    real(kind=dp)              :: energy = 0.0_dp !< SCF energy
    integer                    :: nocc_a = 0      !< Number of occupied alpha orbitals
    integer                    :: nocc_b = 0      !< Number of occupied beta orbitals
  contains
    procedure, private, pass :: init  => scf_data_init
    procedure, private, pass :: clean => scf_data_clean
  end type scf_data_t

  !> @brief Storage of SCF iteration history using a ring buffer of scf_data_t
  !> @detail Manages a cyclic buffer of SCF data up to num_slots. Data is stored in a single
  !>         array of scf_data_t for better encapsulation and cache efficiency.
  type :: converger_data
    integer                   :: slot = 0         !< Current slot (1-based index)
    integer                   :: num_saved = 0    !< Number of occupied slots
    integer                   :: num_slots = 0    !< Total number of slots in buffer
    integer                   :: num_focks = 0    !< Number of Fock matrices per iteration (1 for RHF/ROHF, 2 for UHF)
    integer                   :: ldim = 0         !< Size of square matrices (nbf)
    type(scf_data_t), allocatable :: buffer(:)    !< Ring buffer of SCF data
  contains
    procedure, private, pass :: init         => conv_data_init
    procedure, private, pass :: clean        => conv_data_clean
    procedure, private, pass :: next_slot    => conv_data_next_slot
    procedure, private, pass :: discard_last => conv_data_discard
    procedure, private, pass :: put          => conv_data_put
    procedure, private, pass :: get_fock     => conv_data_get_fock
    procedure, private, pass :: get_mo_a     => conv_data_get_mo_a
    procedure, private, pass :: get_mo_b     => conv_data_get_mo_b
    procedure, private, pass :: get_mo_e_a   => conv_data_get_mo_e_a
    procedure, private, pass :: get_mo_e_b   => conv_data_get_mo_e_b
    procedure, private, pass :: get_occ_a    => conv_data_get_occ_a
    procedure, private, pass :: get_occ_b    => conv_data_get_occ_b
    procedure, private, pass :: get_density  => conv_data_get_density
    procedure, private, pass :: get_err      => conv_data_get_err
    procedure, private, pass :: get_energy   => conv_data_get_energy
    procedure, private, pass :: get_slot     => conv_data_get_slot
    procedure, private, pass :: get_nocc_a   => conv_data_get_nocc_a
    procedure, private, pass :: get_nocc_b   => conv_data_get_nocc_b
  end type converger_data

  !> @brief Base type for SCF converger results
  !> @detail Used by the main SCF convergence driver `scf_conv` to return results.
  !>  The extending type should provide the following interfaces:
  !>    init  : preliminary initialization of internal subconverger data
  !>    clean : destructor
  !>    setup : setting-up of the sub-converger equations, taking into the account the new data
  !>            added in main SCF driver
  !>    run   : solving the equations, returns `scf_conv_result` datatype. Multiple subsequent calls to
  !>            this procedure without re-running `setup` should not change internal state and have to give same results
  type :: scf_conv_result
    integer :: ierr = 5                            !< Error status (0 = success, nonzero = failure)
    real(kind=dp) :: error = 1.0e99_dp             !< Current DIIS error magnitude
    type(converger_data), pointer :: dat => null() !< Pointer to SCF data storage
    character(len=conv_name_maxlen) :: active_converger_name = '' !< Name of active converger
  contains
    procedure, pass :: get_error   => conv_result_get_error
    procedure, pass :: get_fock    => conv_result_dummy_get_fock
    procedure, pass :: get_density => conv_result_dummy_get_density
    procedure, pass :: get_mo_a    => conv_result_dummy_get_mo_a
    procedure, pass :: get_mo_b    => conv_result_dummy_get_mo_b
    procedure, pass :: get_mo_e_a  => conv_result_dummy_get_mo_e_a
    procedure, pass :: get_mo_e_b  => conv_result_dummy_get_mo_e_b
  end type scf_conv_result

  !> @brief SCF converger results for interpolation methods like DIIS
  !> @details The updated Fock/density is computed as F_n = \sum_{i} F_i * c_i
  type, extends(scf_conv_result) :: scf_conv_interp_result
    real(kind=dp), allocatable :: coeffs(:) !< Interpolation coefficients c_i
    integer :: pstart                       !< Start index of interpolation range
    integer :: pend                         !< End index of interpolation range
  contains
    procedure, pass :: get_fock    => conv_result_interp_get_fock
    procedure, pass :: get_density => conv_result_interp_get_density
    procedure, private, pass :: interpolate => conv_result_interpolate
  end type scf_conv_interp_result

  !> @brief SCF converger results for SOSCF method
  type, extends(scf_conv_result) :: scf_conv_soscf_result
  contains
    procedure, pass :: get_fock    => conv_result_soscf_get_fock
    procedure, pass :: get_mo_a    => conv_result_soscf_get_mo_a
    procedure, pass :: get_mo_b    => conv_result_soscf_get_mo_b
    procedure, pass :: get_mo_e_a  => conv_result_soscf_get_mo_e_a
    procedure, pass :: get_mo_e_b  => conv_result_soscf_get_mo_e_b
  end type scf_conv_soscf_result

  !> @brief Base type for real SCF convergers (subconvergers)
  !> @detail Used by main SCF convergence driver `scf_conv`.
  !>  The extending type should provide the following interfaces:
  !>    init  : preliminary initialization of internal subconverger data
  !>    clean : destructor
  !>    setup : setting-up of the sub-converger equations, taking into the account the new data
  !>            added in main SCF driver
  !>    run   : solving the equations, returns `scf_conv_result` datatype. Multiple subsequent calls to
  !>            this procedure without re-running `setup` should not change internal state and have to give same results
  type, abstract :: subconverger
    integer :: last_setup = 1024                      !< Number of SCF iterations since last setup
    integer :: iter = 0                               !< Number of iterations passed
    character(len=conv_name_maxlen) :: conv_name = '' !< Converger name
    type(converger_data), pointer :: dat => null()    !< Pointer to SCF data
  contains
    private
    procedure, pass :: subconverger_init
    procedure, pass :: subconverger_clean
    procedure, public, pass :: init  => subconverger_init
    procedure, public, pass :: clean => subconverger_clean
    procedure(subconverger_run), pass, deferred :: run
    procedure(subconverger_setup), pass, deferred :: setup
  end type subconverger

  !> @brief Container type for subconvergers
  !> @detail Used to create allocatable array of allocatable types
  type :: subconverger_
    class(subconverger), allocatable :: s
  end type subconverger_

  !> @brief Main driver for converging SCF problems
  !> @detail Manages and runs different real convergers depending on the current
  !>         state of SCF optimization.
  type :: scf_conv
    integer :: step = 0                           !< Current SCF iteration step
    real(kind=dp), pointer :: overlap(:,:) => null()      !< Pointer to full-format overlap matrix (S)
    real(kind=dp), pointer :: overlap_sqrt(:,:) => null() !< Pointer to full-format S^(1/2) matrix
    type(converger_data) :: dat                   !< Storage of SCF iteration history
    type(subconverger_), allocatable :: sconv(:)  !< Array of subconverger methods
    real(kind=dp), allocatable :: thresholds(:)   !< Thresholds to initiate subconvergers
    integer :: iter_space_size = 10               !< Default size of subconverger problem space
    integer :: verbose = 0                        !< Verbosity level
    integer :: state = 0                          !< Current state (0 = not initialized, 1 = initialized)
    real(kind=dp) :: current_error = 1.0e99_dp    !< Maximum absolute value of current DIIS error
  contains
    procedure, pass :: init     => scf_conv_init
    procedure, pass :: clean    => scf_conv_clean
    procedure, pass :: add_data => scf_conv_add_data
    procedure, private, pass :: select_method => scf_conv_select
    procedure, pass :: run      => scf_conv_run
    procedure, private, pass :: compute_error => scf_conv_compute_error
  end type scf_conv

  !> @brief Dummy (steepest descent) converger
  !> @detail Used for convenience, does nothing but return the latest data.
  type, extends(subconverger) :: noconv_converger
  contains
    procedure, pass :: init  => noconv_init
    procedure, pass :: run   => noconv_run
    procedure, pass :: setup => noconv_setup
  end type noconv_converger

  !> @brief Commutator DIIS (C-DIIS) converger
  !> @detail Solves constraint minimization problem: min { Ax, \sum_i x_i = 1 }
  !>         where A_{ij} = Tr([F_i D_i S_i], [F_j D_j S_j])
  type, extends(subconverger) :: cdiis_converger
    integer :: maxdiis                            !< Maximum number of DIIS vectors
    real(kind=dp), allocatable :: a(:,:)          !< C-DIIS A-matrix
    integer :: verbose = 0                        !< Verbosity parameter
    integer :: old_dim = 0                        !< Dimension of A matrix from previous step
  contains
    procedure, pass :: init  => cdiis_init
    procedure, pass :: clean => cdiis_clean
    procedure, pass :: run   => cdiis_run
    procedure, pass :: setup => cdiis_setup
  end type cdiis_converger

  !> @brief Datatype to pass optimization parameters to NLOpt for E/A-DIIS
  !> @detail Used for solving constraint minimization problems in E/A-DIIS:
  !>         \f$ min \{ x^T A x + bx,\ x_i \geq 0,\ \sum_i{x_i} = 1 \} \f$
  type :: ediis_opt_data
    real(kind=8), pointer :: b(:) => null()       !< Energy or linear term coefficients
    real(kind=8), pointer :: A(:,:) => null()     !< Quadratic term matrix
    procedure(eadiis_f), pointer, nopass :: fun => null() !< Objective function pointer
  end type ediis_opt_data

  !> @brief Energy DIIS (E-DIIS) converger
  !> @detail Optimizes the following function:
  !>         \f$ E_\mathrm{E-DIIS} = \sum_{i} c_i E_i - 0.5 \sum_{i,j} c_i c_j Tr( (D_i-D_j) (F_i-F_j) ) \f$
  !>         under the constraints:
  !>         \f$ \sum_i {c_i} = 1, c_i \geq 0 \f$
  type, extends(cdiis_converger) :: ediis_converger
    real(kind=dp), allocatable :: b(:)            !< Energy history
    real(kind=dp), allocatable :: xlog(:,:)       !< Interpolation coefficients history
    type(ediis_opt_data) :: t                     !< Optimization data for NLOpt
    procedure(eadiis_f), pointer, nopass :: fun => null() !< Pointer to objective function
  contains
    procedure, pass :: init  => ediis_init
    procedure, pass :: clean => ediis_clean
    procedure, pass :: setup => ediis_setup
    procedure, pass :: run   => ediis_run
  end type ediis_converger

  !> @brief A-DIIS subconverger type
  !> @detail Optimizes the following function:
  !>         \f$ E_\mathrm{A-DIIS} = \sum_{i} c_i Tr((D_i-D_n)F_n) + 2 \sum_{i,j} c_i c_j Tr( (D_i-D_n) (F_j-F_n) ) \f$
  !>         under the constraints:
  !>         \f$ \sum_i {c_i} = 1, c_i \geq 0 \f$
  type, extends(ediis_converger) :: adiis_converger
  contains
    procedure, pass :: init  => adiis_init
    procedure, pass :: setup => adiis_setup
  end type adiis_converger

  abstract interface
    subroutine subconverger_run(self, res)
      import
      class(subconverger), target, intent(inout) :: self
      class(scf_conv_result), allocatable, intent(out) :: res
    end subroutine

    subroutine subconverger_setup(self)
      import
      class(subconverger), intent(inout) :: self
    end subroutine

    subroutine eadiis_f(val, n, t, grad, need_gradient, d)
      import
      implicit none
      real(kind=8) :: val, t(*), grad(*)
      integer(kind=4), intent(in) :: n, need_gradient
      type(ediis_opt_data), intent(in) :: d
    end subroutine
  end interface

  !> @brief SOSCF subconverger type
  !> @detail Implementation is based on
  !          T. H. Fischer, J. Almlof, Journal of Physical Chemistry, 96(24), 9768-9774 (1992)
  !          F. Neese, Chemical Physics Letters, 325(1-3), 93-98 (2000)
  type, extends(subconverger) :: soscf_converger

    integer :: verbose = 0            !< Verbosity parameter
    integer :: nfocks = 0             !< Number of Focks (RHF = 1, UHF = 2)
    integer :: nbf = 0                !< Number of basis functions
    integer :: nbf_tri = 0            !< Number of basis functions in triangular format
    integer :: scf_type = 0           !< SCF type: 1:RHF, 2:UHF, 3:ROHF
    integer :: nocc_a = 0             !< Number of occupied alpha orbitals
    integer :: nocc_b = 0             !< Number of occupied beta orbitals
    integer :: nvec = 0               !< Size of the gradient vector

    real(kind=dp), pointer :: overlap(:, :) => null()
    real(kind=dp), pointer :: overlap_invsqrt(:, :) => null()

    ! SOSCF parameters from scf_driver:
    integer :: soscf_start = -1
    integer :: soscf_freq = 2
    logical :: soscf_diis_alternate = .False.
    integer :: max_iter = 10                  !< Maximum micro-iterations
    integer :: min_iter = 1                   !< Minimum micro-iterations
    real(kind=dp) :: soscf_conv = 1.0e-3_dp
    real(kind=dp) :: hess_thresh = 1.0e-10_dp !< Orbital Hessian threshold
    real(kind=dp) :: grad_thresh = 1.0e-3_dp  !< Gradient threshold
    real(kind=dp) :: level_shift = 0.2_dp     !< Level shifting parameter
    logical :: use_lineq = .false.            !< Use linear equations (vs BFGS)

    ! L-BFGS history
    real(kind=dp), allocatable :: s_history(:,:)  !< Step history (nvec, m_max)
    real(kind=dp), allocatable :: rho_history(:)  !< Curvature reciprocals (m_max)
    real(kind=dp), allocatable :: y_history(:,:)  !< Gradient difference history (nvec, m_max)
    real(kind=dp), allocatable :: grad_prev(:)    !< Previous gradient (nvec)
    real(kind=dp), allocatable :: x_prev(:)       !< Previous rotation parameters (nvec)
    real(kind=dp), allocatable :: h_inv(:)        !< Initial inverse Hessian diagonal (nvec)
    real(kind=dp), allocatable :: work_1(:,:)     !< Work matrix (nbf, nbf)
    real(kind=dp), allocatable :: work_2(:,:)     !< Work matrix (nbf, nbf)
    real(kind=dp), allocatable :: mo_a(:,:)       !< MOs matrix (nbf, nbf)
    real(kind=dp), allocatable :: mo_b(:,:)       !< MOs matrix (nbf, nbf)
    real(kind=dp), allocatable :: dens_a(:)        !< MOs matrix (nbf_tri)
    real(kind=dp), allocatable :: dens_b(:)        !< MOs matrix (nbf_tri)
    integer :: m_max = 0                          !< Maximum number of stored history steps
    integer :: m_history = 0                      !< Number of stored history steps
    logical :: first_macro = .true.               !< Flag for first macro-iteration
  contains
    procedure, pass :: init  => soscf_init
    procedure, pass :: clean => soscf_clean
    procedure, pass :: setup => soscf_setup
    procedure, pass :: run   => soscf_run
  end type soscf_converger

contains

!==============================================================================
! scf_data_t Methods
!==============================================================================

  !> @brief Initialize an scf_data_t instance
  !> @param[inout] self The scf_data_t object to initialize.
  !> @param[in] ldim Size of triangular matrices (nbf*(nbf+1)/2)
  !> @param[in] nfocks Number of Fock matrices (1 for RHF/ROHF, 2 for UHF)
  !> @param[out] istat Status code (0 for success, nonzero for allocation failure)
  subroutine scf_data_init(self, ldim, nfocks, istat)
    class(scf_data_t), intent(inout) :: self
    integer, intent(in) :: ldim, nfocks
    integer, intent(out) :: istat
    integer :: nbf_tri

    istat = 0
    nbf_tri = ldim * (ldim + 1) / 2
    if (allocated(self%focks)) call self%clean()
    allocate(self%focks(nbf_tri, nfocks), &
             self%densities(nbf_tri, nfocks), &
             self%errs(nbf_tri, nfocks), &
             stat=istat)
    ! Note: MO coefficients mo_a, mo_b and energies mo_e_a, mo_e_b
    ! Note: are allocated on demand in conv_data_put

    self%focks     = 0.0_dp
    self%densities = 0.0_dp
    self%errs      = 0.0_dp
    self%energy    = 0.0_dp
    self%nocc_a     = 0
    self%nocc_b     = 0
  end subroutine scf_data_init

  !> @brief Clean up an scf_data_t instance
  subroutine scf_data_clean(self)
    class(scf_data_t), intent(inout) :: self

    if (allocated(self%focks))     deallocate(self%focks)
    if (allocated(self%densities)) deallocate(self%densities)
    if (allocated(self%errs))      deallocate(self%errs)
    self%energy = 0.0_dp
    self%nocc_a  = 0
    self%nocc_b  = 0
  end subroutine scf_data_clean

!==============================================================================
! converger_data Methods
!==============================================================================

  !> @brief Set up a ring buffer to store SCF iteration history.
  !> @param[inout] self The converger_data object to initialize.
  !> @param[in] ldim Number of basis functions (nbf)
  !> @param[in] nfocks Number of Fock matrices per iteration (1 for RHF/ROHF, 2 for UHF)
  !> @param[in] nslots Maximum number of SCF iterations to store
  !> @param[out] istat Success status (0 for success, nonzero for error)
  subroutine conv_data_init(self, ldim, nfocks, nslots, istat)
    class(converger_data), intent(inout) :: self
    integer, intent(in) :: ldim, nfocks, nslots
    integer, intent(out) :: istat
    integer :: i

    istat = 0
    if (allocated(self%buffer)) call self%clean()
    self%ldim      = ldim
    self%num_focks = nfocks
    self%num_slots = nslots
    self%num_saved = 0
    self%slot      = 0
    allocate(self%buffer(nslots), stat=istat)
    do i = 1, nslots
      call self%buffer(i)%init(ldim, nfocks, istat)
    end do
  end subroutine conv_data_init

  !> @brief Finalize converger_data and deallocate memory
  subroutine conv_data_clean(self)
    class(converger_data), intent(inout) :: self
    integer :: i

    if (allocated(self%buffer)) then
      do i = 1, self%num_slots
        call self%buffer(i)%clean()
      end do
      deallocate(self%buffer)
    end if
    self%ldim      = 0
    self%num_focks = 0
    self%num_slots = 0
    self%num_saved = 0
    self%slot      = 0
  end subroutine conv_data_clean

  !> @brief Advance to the next slot in the ring buffer
  subroutine conv_data_next_slot(self)
    class(converger_data), intent(inout) :: self

    self%slot = mod(self%slot, self%num_slots) + 1
    self%num_saved = min(self%num_saved + 1, self%num_slots)
  end subroutine conv_data_next_slot

  !> @brief Discard the most recent data entry
  subroutine conv_data_discard(self)
    class(converger_data), intent(inout) :: self

    self%slot = mod(self%slot - 2, self%num_slots) + 1
    self%num_saved = min(self%num_saved - 1, 1)
  end subroutine conv_data_discard

  !> @brief Store SCF data for the current iteration
  !> @param[in] fock Fock matrices (optional)
  !> @param[in] dens Density matrices (optional)
  !> @param[in] energy SCF energy (optional)
  !> @param[in] mo_a Alpha MO coefficients (optional)
  !> @param[in] mo_b Beta MO coefficients (optional)
  !> @param[in] mo_e_a Alpha MO energies (optional)
  !> @param[in] mo_e_b Beta MO energies (optional)
  !> @param[in] nocc_a Number of occupied alpha orbitals (optional)
  !> @param[in] nocc_b Number of occupied beta orbitals (optional)
  subroutine conv_data_put(self, fock, dens, energy, mo_a, mo_b, &
                           mo_e_a, mo_e_b, occ_a, occ_b, &
                           nocc_a, nocc_b)
    class(converger_data), intent(inout) :: self
    real(kind=dp), intent(in), optional :: fock(:,:)
    real(kind=dp), intent(in), optional :: dens(:,:)
    real(kind=dp), intent(in), optional :: energy
    real(kind=dp), intent(in), optional :: mo_a(:,:), mo_b(:,:)
    real(kind=dp), intent(in), optional :: mo_e_a(:), mo_e_b(:)
    real(kind=dp), intent(in), optional :: occ_a(:), occ_b(:)
    integer, intent(in), optional :: nocc_a, nocc_b

    integer :: slot, nbf

    nbf = self%ldim
    slot = self%slot
    if (present(fock)) self%buffer(slot)%focks = fock
    if (present(dens)) self%buffer(slot)%densities = dens
    if (present(energy)) self%buffer(slot)%energy = energy
    if (present(nocc_a)) self%buffer(slot)%nocc_a = nocc_a
    if (present(nocc_b)) self%buffer(slot)%nocc_b = nocc_b
    if (present(mo_a)) then
      if (.not. allocated(self%buffer(slot)%mo_a)) then
        allocate(self%buffer(slot)%mo_a(nbf, nbf))
      end if
      self%buffer(slot)%mo_a = mo_a
    end if
    if (present(mo_b)) then
      if (.not. allocated(self%buffer(slot)%mo_b)) then
        allocate(self%buffer(slot)%mo_b(nbf, nbf))
      end if
      self%buffer(slot)%mo_b = mo_b
    end if
    if (present(mo_e_a)) then
      if (.not. allocated(self%buffer(slot)%mo_e_a)) then
        allocate(self%buffer(slot)%mo_e_a(nbf))
      end if
      self%buffer(slot)%mo_e_a = mo_e_a
    end if
    if (present(mo_e_b)) then
      if (.not. allocated(self%buffer(slot)%mo_e_b)) then
        allocate(self%buffer(slot)%mo_e_b(nbf))
      end if
      self%buffer(slot)%mo_e_b = mo_e_b
    end if
    if (present(occ_a)) then
      if (.not. allocated(self%buffer(slot)%occ_a)) then
        allocate(self%buffer(slot)%occ_a(size(occ_a)))
      end if
      self%buffer(slot)%occ_a = occ_a
    end if
    if (present(occ_b)) then
      if (.not. allocated(self%buffer(slot)%occ_b)) then
        allocate(self%buffer(slot)%occ_b(size(occ_b)))
      end if
      self%buffer(slot)%occ_b = occ_b
    end if
  end subroutine conv_data_put

  !> @brief Get Fock matrix for a specific iteration and matrix ID
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @param[in] matrix_id Fock matrix index (1 = alpha, 2 = beta)
  !> @return Pointer to the Fock matrix
  function conv_data_get_fock(self, n, matrix_id) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n, matrix_id
    real(kind=dp), pointer :: res(:)
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%focks(:, matrix_id)
  end function conv_data_get_fock

  !> @brief Get MO coefficients for alpha orbitals for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return Pointer to the alpha MO coefficients
  function conv_data_get_mo_a(self, n) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n
    real(kind=dp), pointer :: res(:,:)
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%mo_a
  end function conv_data_get_mo_a

  !> @brief Get MO coefficients for beta orbitals for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return Pointer to the beta MO coefficients
  function conv_data_get_mo_b(self, n) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n
    real(kind=dp), pointer :: res(:,:)
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%mo_b
  end function conv_data_get_mo_b

  !> @brief Get alpha MO energies for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return Pointer to the alpha MO energies
  function conv_data_get_mo_e_a(self, n) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n
    real(kind=dp), pointer :: res(:)
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%mo_e_a
  end function conv_data_get_mo_e_a

  !> @brief Get beta MO energies for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return Pointer to the beta MO energies
  function conv_data_get_mo_e_b(self, n) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n
    real(kind=dp), pointer :: res(:)
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%mo_e_b
  end function conv_data_get_mo_e_b

  !> @brief Get density matrix for a specific iteration and matrix ID
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @param[in] matrix_id Density matrix index (1 = alpha, 2 = beta)
  !> @return Pointer to the density matrix
  function conv_data_get_density(self, n, matrix_id) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n, matrix_id
    real(kind=dp), pointer :: res(:)
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%densities(:, matrix_id)
  end function conv_data_get_density

  !> @brief Get DIIS error matrix for a specific iteration and matrix ID
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @param[in] matrix_id Error matrix index (1 = alpha, 2 = beta)
  !> @return Pointer to the error matrix
  function conv_data_get_err(self, n, matrix_id) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n, matrix_id
    real(kind=dp), pointer :: res(:)
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%errs(:, matrix_id)
  end function conv_data_get_err

  !> @brief Get SCF energy for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return SCF energy value
  function conv_data_get_energy(self, n) result(res)
    class(converger_data), intent(in) :: self
    integer, intent(in) :: n
    real(kind=dp) :: res
    integer :: slot

    slot = self%get_slot(n)
    res = self%buffer(slot)%energy
  end function conv_data_get_energy

  !> @brief Compute slot index from iteration number
  !> @param[in] n Iteration number (1 = oldest, -1 = latest)
  !> @return Slot index in the ring buffer
  function conv_data_get_slot(self, n) result(slot)
    class(converger_data), intent(in) :: self
    integer, intent(in) :: n
    integer :: slot, num_saved
    num_saved = self%num_saved
    if (n == -1) then
      slot = self%slot
    else
      slot = modulo(self%slot - num_saved + n - 1, self%num_slots) + 1
    end if
  end function

  !> @brief Get alpha orbital occupations for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return Alpha orbital occupations
  function conv_data_get_occ_a(self, n) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n
    real(kind=dp), pointer :: res(:)
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%occ_a
  end function conv_data_get_occ_a

  !> @brief Get beta orbital occupations for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return Beta orbital occupations
  function conv_data_get_occ_b(self, n) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n
    real(kind=dp), pointer :: res(:)
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%occ_b
  end function conv_data_get_occ_b

  !> @brief Get number of occupied alpha orbitals for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return Number of occupied alpha orbitals
  function conv_data_get_nocc_a(self, n) result(res)
    class(converger_data), intent(in) :: self
    integer, intent(in) :: n
    integer :: res, slot

    slot = self%get_slot(n)
    res = self%buffer(slot)%nocc_a
  end function conv_data_get_nocc_a

  !> @brief Get number of occupied beta orbitals for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return Number of occupied beta orbitals
  function conv_data_get_nocc_b(self, n) result(res)
    class(converger_data), intent(in) :: self
    integer, intent(in) :: n
    integer :: res, slot

    slot = self%get_slot(n)
    res = self%buffer(slot)%nocc_b
  end function conv_data_get_nocc_b

!==============================================================================
! scf_conv_result Methods
!==============================================================================

  !> @brief Form the new Fock matrix
  !> @detail Placeholder for derived types to override. Does nothing by default.
  subroutine conv_result_dummy_get_fock(self, matrix, istat)
    class(scf_conv_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    istat = 0
  end subroutine conv_result_dummy_get_fock

  !> @brief Form the new density matrix
  !> @detail Placeholder for derived types to override. Does nothing by default.
  subroutine conv_result_dummy_get_density(self, matrix, istat)
    class(scf_conv_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    istat = 0
  end subroutine conv_result_dummy_get_density

  subroutine conv_result_dummy_get_mo_a(self, matrix, istat)
    class(scf_conv_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    istat = 0
  end subroutine conv_result_dummy_get_mo_a

  subroutine conv_result_dummy_get_mo_b(self, matrix, istat)
    class(scf_conv_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    istat = 0
  end subroutine conv_result_dummy_get_mo_b

  subroutine conv_result_dummy_get_mo_e_a(self, vector, istat)
    class(scf_conv_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: vector(:)

    istat = 0
  end subroutine conv_result_dummy_get_mo_e_a

  subroutine conv_result_dummy_get_mo_e_b(self, vector, istat)
    class(scf_conv_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: vector(:)

    istat = 0
  end subroutine conv_result_dummy_get_mo_e_b

  !> @brief Form the interpolated Fock matrix
  !> @detail F_n = \sum_{i=start}^{end} F_i * c_i
  subroutine conv_result_interp_get_fock(self, matrix, istat)
    class(scf_conv_interp_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    if (self%ierr == 0) then
      call self%interpolate(matrix, 'fock', istat)
    else
      istat = self%ierr
    end if
  end subroutine conv_result_interp_get_fock

  !> @brief Get error value from result
  !> @return Current error value
  !> @brief Get the current convergence error
  !> @return Current error value for the active converger
  function conv_result_get_error(self) result(err)
    class(scf_conv_result), intent(in) :: self
    real(kind=dp) :: err

    err = self%error
  end function conv_result_get_error

!==============================================================================
! DIIS Result Methods
!==============================================================================

  !> @brief Form the interpolated density matrix
  !> @detail D_n = \sum_{i=start}^{end} D_i * c_i
  subroutine conv_result_interp_get_density(self, matrix, istat)
    class(scf_conv_interp_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    if (self%ierr == 0) then
      call self%interpolate(matrix, 'density', istat)
    else
      istat = self%ierr
    end if
  end subroutine conv_result_interp_get_density

  !> @brief Form the interpolated matrix (Fock or density)
  !> @param[inout] matrix Output matrix (Fock or density)
  !> @param[in] datatype 'fock' or 'density' to specify which matrix to interpolate
  !> @param[out] istat Status code (0 = success)
  subroutine conv_result_interpolate(self, matrix, datatype, istat)
    class(scf_conv_interp_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)
    character(len=*), intent(in) :: datatype
    integer :: i, ifock

    if (self%ierr /= 0) then
      istat = self%ierr
      return
    end if

    matrix = 0.0_dp
    do ifock = 1, self%dat%num_focks
      do i = self%pstart, self%pend
        select case (datatype)
        case ('fock')
          matrix(:, ifock) = matrix(:, ifock) + self%coeffs(i) * self%dat%get_fock(i, ifock)
        case ('density')
          matrix(:, ifock) = matrix(:, ifock) + self%coeffs(i) * self%dat%get_density(i, ifock)
        end select
      end do
    end do

    istat = 0
  end subroutine conv_result_interpolate

!==============================================================================
! SOSCF Result Methods
!==============================================================================

  !> @brief Get Fock matrix from SOSCF result
  subroutine conv_result_soscf_get_fock(self, matrix, istat)
    class(scf_conv_soscf_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    integer :: i

    if (self%ierr /= 0) then
      istat = self%ierr
      return
    end if

    do i = 1, self%dat%num_focks
      matrix(:,i) = self%dat%get_fock(-1, i)
    end do

  end subroutine conv_result_soscf_get_fock

  !> @brief Get alpha MO coefficients from SOSCF result
  subroutine conv_result_soscf_get_mo_a(self, matrix, istat)
    class(scf_conv_soscf_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    if (self%ierr /= 0) then
      istat = self%ierr
      return
    end if

    matrix = self%dat%buffer(self%dat%slot)%mo_a
    istat = 0
  end subroutine conv_result_soscf_get_mo_a

  !> @brief Get beta MO coefficients from SOSCF result
  subroutine conv_result_soscf_get_mo_b(self, matrix, istat)
    class(scf_conv_soscf_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: matrix(:,:)

    if (self%ierr /= 0) then
      istat = self%ierr
      return
    end if

    matrix = self%dat%buffer(self%dat%slot)%mo_b
    istat = 0
  end subroutine conv_result_soscf_get_mo_b

  !> @brief Get alpha orbital energies from SOSCF result
  subroutine conv_result_soscf_get_mo_e_a(self, vector, istat)
    class(scf_conv_soscf_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: vector(:)

    if (self%ierr /= 0) then
      istat = self%ierr
      return
    end if

    vector = self%dat%buffer(self%dat%slot)%mo_e_a
    istat = 0
  end subroutine conv_result_soscf_get_mo_e_a

  !> @brief Get beta orbital energies from SOSCF result
  subroutine conv_result_soscf_get_mo_e_b(self, vector, istat)
    class(scf_conv_soscf_result), intent(in) :: self
    integer, intent(out) :: istat
    real(kind=dp), intent(inout) :: vector(:)

    if (self%ierr /= 0) then
      istat = self%ierr
      return
    end if

    vector = self%dat%buffer(self%dat%slot)%mo_e_b
    istat = 0
  end subroutine conv_result_soscf_get_mo_e_b

!==============================================================================
! scf_conv Methods
!==============================================================================

  !> @brief Finalize scf_conv datatype
  subroutine scf_conv_clean(self)
    class(scf_conv), intent(inout) :: self

    self%verbose = 0
    self%step = 0
    self%state = conv_state_not_initialized
    call self%dat%clean()
    if (allocated(self%thresholds)) deallocate(self%thresholds)
    if (allocated(self%sconv)) deallocate(self%sconv)
    nullify(self%overlap)
    nullify(self%overlap_sqrt)
  end subroutine scf_conv_clean

  !> @brief Initializes the SCF converger driver
  !> @param[in] ldim Number of orbitals
  !> @param[in] maxvec Size of SCF converger linear space (e.g., number of DIIS vectors)
  !> @param[in] subconvergers Array of SCF subconverger codes
  !> @param[in] thresholds Thresholds to initiate subconvergers, where
  !>                       subconvergers[i] runs when current error is less than thresholds[i]
  !> @param[in] overlap Overlap matrix (S) in full format
  !> @param[in] overlap_sqrt S^(1/2) matrix in full format
  !> @param[in] num_focks 1 if R/ROHF, 2 if UHF
  !> @param[in] verbose Verbosity level
  subroutine scf_conv_init(self, ldim, maxvec, subconvergers, thresholds, &
                          overlap, overlap_sqrt, num_focks, verbose)
    class(scf_conv), intent(inout) :: self
    integer, intent(in) :: ldim
    integer, optional, intent(in) :: maxvec
    integer, optional, intent(in) :: subconvergers(:)
    real(kind=dp), optional, intent(in) :: thresholds(:)
    real(kind=dp), optional, target, intent(in) :: overlap(:,:), overlap_sqrt(:,:)
    integer, optional, intent(in) :: num_focks
    integer, optional, intent(in) :: verbose
    integer :: nfocks, istat, i

    if (self%state /= conv_state_not_initialized) call self%clean()

    if (present(thresholds)) then
      allocate(self%thresholds(0:ubound(thresholds,1)))
      self%thresholds(0:) = [thresholds, 0.0_dp]
    else
      allocate(self%thresholds(0:1))
      self%thresholds(0:) = [1.0_dp, 0.0_dp]
    end if

    self%overlap => null()
    if (present(overlap)) self%overlap => overlap

    self%overlap_sqrt => null()
    if (present(overlap_sqrt)) self%overlap_sqrt => overlap_sqrt

    nfocks = 1
    if (present(num_focks)) nfocks = num_focks

    self%verbose = 0
    if (present(verbose)) self%verbose = verbose

    self%iter_space_size = 15
    if (present(maxvec)) self%iter_space_size = maxvec

    self%step = 0

    call self%dat%init(ldim, nfocks, self%iter_space_size, istat)
    if (istat /= 0) then
      self%state = conv_state_not_initialized
      return
    end if

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
        case (conv_soscf)
          allocate(soscf_converger :: self%sconv(i)%s)
        end select
        call self%sconv(i)%s%init(self)
      end do
    end if

    self%state = conv_state_initialized

  end subroutine scf_conv_init

  !> @brief Store data from the new SCF iteration
  !> @param[in] f Fock matrix/matrices
  !> @param[in] dens Density matrix/matrices
  !> @param[in] e SCF energy (optional)
  subroutine scf_conv_add_data(self, f, dens, e, mo_a, mo_b, mo_e_a, mo_e_b, &
                               occ_a, occ_b, nocc_a, nocc_b)
    class(scf_conv), intent(inout) :: self
    real(kind=dp), intent(in), optional :: f(:,:)    ! Fock matrices
    real(kind=dp), intent(in), optional :: dens(:,:) ! Density matrices
    real(kind=dp), intent(in), optional :: e         ! SCF energy
    real(kind=dp), intent(in), optional :: mo_a(:,:) ! Alpha MO coefficients
    real(kind=dp), intent(in), optional :: mo_b(:,:) ! Beta MO coefficients
    real(kind=dp), intent(in), optional :: mo_e_a(:) ! Alpha MO energies
    real(kind=dp), intent(in), optional :: mo_e_b(:) ! Beta MO energies
    real(kind=dp), intent(in), optional :: occ_a(:)  ! Alpha orbital occupations
    real(kind=dp), intent(in), optional :: occ_b(:)  ! Beta orbital occupations
    integer, intent(in), optional :: nocc_a          ! Number of occupied alpha orbitals
    integer, intent(in), optional :: nocc_b          ! Number of occupied beta orbitals
    integer :: i

    call self%dat%next_slot()

    ! Save the current Fock and density matrices
    if (present(f)) call self%dat%put(fock=f)
    if (present(dens)) call self%dat%put(dens=dens)
    if (present(e)) call self%dat%put(energy=e)
    if (present(mo_a)) call self%dat%put(mo_a=mo_a)
    if (present(mo_b)) call self%dat%put(mo_b=mo_b)
    if (present(mo_e_a)) call self%dat%put(mo_e_a=mo_e_a)
    if (present(mo_e_b)) call self%dat%put(mo_e_b=mo_e_b)
    if (present(occ_a)) call self%dat%put(occ_a=occ_a)
    if (present(occ_b)) call self%dat%put(occ_b=occ_b)
    if (present(nocc_a)) call self%dat%put(nocc_a=nocc_a)
    if (present(nocc_b)) call self%dat%put(nocc_b=nocc_b)

    ! Compute the current error
    self%current_error = self%compute_error()

    ! Update subconverger states
    do i = lbound(self%sconv, 1), ubound(self%sconv, 1)
      self%sconv(i)%s%last_setup = self%sconv(i)%s%last_setup + 1
    end do
  end subroutine scf_conv_add_data

  !> @brief Computes the new guess to the SCF wavefunction
  !> @param[out] conv_result Results of the calculation
  subroutine scf_conv_run(self, conv_result)
    class(scf_conv), target, intent(inout) :: self
    class(scf_conv_result), allocatable, intent(out) :: conv_result
    class(subconverger), pointer :: conv
    class(scf_conv_result), allocatable :: tmp_result

    conv => self%select_method(self%current_error)

    if (self%state == conv_state_not_initialized) then
      conv_result = scf_conv_result(error=self%current_error)
      return
    end if

    conv%iter = conv%iter + 1
    self%step = self%step + 1

    call conv%setup()
    ! Nothing left to do on the first iteration, exit
    if (self%step == 1) then
      conv_result = scf_conv_result( &
                      ierr=0, &
                      active_converger_name='SD', &
                      error=self%current_error)
      return
    end if

    ! Solve the set of DIIS linear equations
    call conv%run(tmp_result)
    call move_alloc(from=tmp_result, to=conv_result)
  end subroutine scf_conv_run

  !> @brief Select subconverger based on the current DIIS error value
  !> @param[in] error DIIS error value
  !> @return Pointer to selected converger
  function scf_conv_select(self, error) result(conv)
    class(scf_conv), target, intent(in) :: self
    real(kind=dp), intent(in) :: error
    class(subconverger), pointer :: conv
    integer :: i, nconv
    logical :: use_soscf_now

    ! Initial selection based on error thresholds
    nconv = ubound(self%thresholds, 1)
    do i = 0, nconv
      if (error > self%thresholds(i)) exit
    end do

    ! Initially select converger based on error thresholds
    conv => self%sconv(min(i, nconv))%s

    ! Check if any of the convergers is SOSCF
    use_soscf_now = .false.
    do i = 0, nconv
      select type (sc => self%sconv(i)%s)
      type is (soscf_converger)
        ! SOSCF exists in the set of convergers
        ! Check if we should use it based on current iteration
        if (sc%soscf_start >= 0 .and. self%step >= sc%soscf_start) then
          use_soscf_now = .true.
        end if
        exit
      end select
    end do

    ! Handle alternating between DIIS and SOSCF if configured: soscf_diis_alternate
    if (use_soscf_now) then
      select type (selected_conv => conv)
      type is (soscf_converger)
        ! Current selection is SOSCF - check if we should alternate
        if (selected_conv%soscf_diis_alternate) then
          ! Apply SOSCF based on frequency, otherwise use DIIS
          if (mod(self%step - selected_conv%soscf_start, selected_conv%soscf_freq) == 0) then
            ! This iteration aligns with SOSCF frequency - keep SOSCF
            conv => selected_conv
          else
            ! Use DIIS instead - find the first non-SOSCF converger
            do i = 0, nconv
              select type (sc => self%sconv(i)%s)
              class default
                conv => sc
                exit
              end select
            end do
          end if
        end if
      class default
        ! Current selection is not SOSCF - check if we should alternate
        ! Find a SOSCF converger
        do i = 0, nconv
          select type (sc => self%sconv(i)%s)
          type is (soscf_converger)
            if (sc%soscf_diis_alternate) then
              ! Apply SOSCF based on frequency, otherwise keep non-SOSCF
              if (mod(self%step - sc%soscf_start, sc%soscf_freq) == 0) then
                ! This iteration aligns with SOSCF frequency - use SOSCF
                conv => sc
              end if ! Otherwise, retain the non-SOSCF converger (no action needed)
              exit
            end if
          end select
        end do
      end select
    end if

    ! Continue using the 'SD' converger if
    ! already initiated
    if (i == 0 .and. self%step > 0) then
      conv => self%sconv(1)%s
    end if

    ! Use SD by default if no other converger selected
    if (.not. associated(conv)) conv => self%sconv(0)%s
  end function scf_conv_select

  !> @brief Calculate the DIIS error matrix: \f$ \mathrm{Err} = FDS - SDF \f$
  !> @details This routine is general for RHF, ROHF, and UHF.
  !>          Since each of `F`, `D`, `S` are symmetric, this means calculate \f$ FDS \f$,
  !>          and then subtract the transpose from that result.
  !>          Before entry, `F`, `D` and `S` must be expanded to square storage.
  !> @return DIIS error value (infinity norm across all matrices)
  function scf_conv_compute_error(self) result(diis_error)
    use mathlib, only: antisymmetrize_matrix, unpack_matrix, pack_matrix
    use oqp_linalg
    class(scf_conv), target, intent(inout) :: self
    real(kind=dp) :: diis_error
    real(kind=dp), pointer :: f(:), d(:), err(:)
!   all are (nbf, nbf) square matrices
    real(kind=dp), allocatable :: fock_full(:,:), dens_full(:,:), err_full(:,:), wrk(:,:)
    integer :: nbf, ifock, nfocks, slot

    nfocks = self%dat%num_focks
    nbf = self%dat%ldim
    slot = self%dat%slot

    allocate(fock_full(nbf, nbf), &
             dens_full(nbf, nbf), &
             err_full(nbf, nbf), &
             wrk(nbf, nbf))

    diis_error = 0.0_dp
    do ifock = 1, nfocks
      f => self%dat%get_fock(-1, ifock)
      d => self%dat%get_density(-1, ifock)
      call unpack_matrix(f, fock_full, nbf, 'u')
      call unpack_matrix(d, dens_full, nbf, 'u')
      ! F*D
      call dsymm('l', 'u', nbf, nbf, 1.0_dp, fock_full, nbf, dens_full, nbf, 0.0_dp, wrk, nbf)
      ! (F*D)*S
      call dgemm('n', 'n', nbf, nbf, nbf, 1.0_dp, wrk, nbf, self%overlap, nbf, 0.0_dp, err_full, nbf)
      ! F*D*S - S*D*F
      call antisymmetrize_matrix(err_full, nbf)

      ! MV: This step is not really necessary
      ! Put error matrix into consistent orthonormal basis
      ! Pulay uses S**-1/2, but here we use Q, Q obeys Q-dagger*S*Q=I
      ! E-orth = Q-dagger * E * Q, FCKA is used as a scratch `nbf` vector.
      call dgemm('t', 'n', nbf, nbf, nbf, 1.0_dp, self%overlap_sqrt, nbf, err_full, nbf, 0.0_dp, wrk, nbf)
      call dgemm('n', 'n', nbf, nbf, nbf, 1.0_dp, wrk, nbf, self%overlap_sqrt, nbf, 0.0_dp, err_full, nbf)
      err => self%dat%buffer(slot)%errs(:, ifock)
      call pack_matrix(err_full, nbf, err, 'u')
      ! Compute DIIS error (infinity norm of error matrix)
      diis_error = diis_error + maxval(abs(err))
    end do

    deallocate(fock_full, dens_full, err_full, wrk)
  end function scf_conv_compute_error

!==============================================================================
! Subconverger Methods
!==============================================================================

  !> @brief Initialize subconverger
  !> @detail This subroutine takes SCF converger driver as argument.
  !>         It should be initialized and include all the required parameters.
  !> @param[in] params Current SCF converger driver
  subroutine subconverger_init(self, params)
    class(subconverger), intent(inout) :: self
    type(scf_conv), target, intent(in) :: params

    self%iter = 0
    self%last_setup = 1024
    self%dat => params%dat
  end subroutine subconverger_init

!> @brief Finalize subconverger
  subroutine subconverger_clean(self)
    class(subconverger), intent(inout) :: self

    self%iter = 0
    self%last_setup = 1024
    self%dat => null()
  end subroutine subconverger_clean

!==============================================================================
! noconv_converger Methods
!==============================================================================

!> @brief Initialize SD subconverger
!> @param[in] params current SCF converger driver
  subroutine noconv_init(self, params)
    class(noconv_converger), intent(inout) :: self
    type(scf_conv), target, intent(in) :: params

    if (self%iter > 0) call self%clean()
    call self%subconverger_init(params)
    self%conv_name = 'SD'
  end subroutine noconv_init

!> @brief Computes the new guess to the SCF wavefunction
!> @param[out] res results of the calculation
  subroutine noconv_run(self, res)
    class(noconv_converger), target, intent(inout) :: self
    class(scf_conv_result), allocatable, intent(out) :: res
    real(kind=dp) :: diis_error
    integer :: ifock

    allocate(scf_conv_result :: res)
    res = scf_conv_result(ierr=0, active_converger_name='SD', dat=self%dat)

    ! Compute error
    diis_error = 0.0_dp
    do ifock = 1, self%dat%num_focks
      diis_error = max(diis_error, maxval(abs(self%dat%get_err(-1, ifock))))
    end do
    res%error = diis_error
  end subroutine noconv_run

!> @brief Prepare subconverger to run
  subroutine noconv_setup(self)
    class(noconv_converger), intent(inout) :: self

    self%last_setup = 0
  end subroutine noconv_setup

!==============================================================================
! cdiis_converger Methods
!==============================================================================

!> @brief Initialize C-DIIS subconverger
!> @detail This subroutine takes SCF converger driver as argument.
!>         It should be initialized and include all the required parameters.
!> @param[in] params Current SCF converger driver
  subroutine cdiis_init(self, params)
    class(cdiis_converger), intent(inout) :: self
    type(scf_conv), target, intent(in) :: params

    if (self%iter > 0) call self%clean()
    call self%subconverger_init(params)
    self%conv_name = 'C-DIIS'
    self%verbose = params%verbose
    self%maxdiis = params%iter_space_size
    allocate(self%a(self%maxdiis, self%maxdiis), source=0.0_dp)
  end subroutine cdiis_init

!> @brief Finalize C-DIIS subconverger
  subroutine cdiis_clean(self)
    class(cdiis_converger), intent(inout) :: self

    call self%subconverger_clean()
    self%verbose = 0
    if (allocated(self%a)) deallocate(self%a)
  end subroutine cdiis_clean

  !> @brief Computes the new guess to the SCF wavefunction using C-DIIS
  !> @param[out] res Results of the calculation
  subroutine cdiis_run(self, res)
    use mathlib, only: solve_linear_equations
    use io_constants, only: iw
    class(cdiis_converger), target, intent(inout) :: self
    class(scf_conv_result), allocatable, intent(out) :: res
      integer :: i, na, cur, info, ifock
    real(kind=dp) :: a_loc(self%maxdiis+1, self%maxdiis+1)
    real(kind=dp) :: x_loc(self%maxdiis+1)
    real(kind=dp), allocatable :: x(:)
    real(kind=dp), pointer :: err(:)
    real(kind=dp) :: diis_error

    allocate(scf_conv_interp_result :: res)
    res%dat => self%dat
    res%active_converger_name = self%conv_name
    select type (res)
    class is (scf_conv_interp_result)
      res%pstart = 1
      res%pend = self%dat%num_saved
    end select

    res%ierr = 3 ! Need to set up DIIS equations first
    if (self%last_setup /= 0) return

    na = self%dat%num_saved
    ! Solve the set of DIIS linear equations
    a_loc(:self%maxdiis, :self%maxdiis) = self%a
    a_loc(:, na+1) = -1.0_dp
    do i = na, 1, -1
      ! Helper index, needed for dimension reduction in case of instability
      cur = na - i + 1
      x_loc = 0.0_dp
      x_loc(na+1) = -1.0_dp
      info = 0
      call solve_linear_equations(a_loc(cur:, cur:), x_loc(cur:), i+1, 1, self%maxdiis+1, info)
      if (info <= 0) exit
      write(iw, *) 'Reducing DIIS Equation size by 1 for numerical stability'
    end do

    if (info < 0) then
      res%ierr = 2 ! Illegal value in DSYSV
    else if (info > 0) then
      res%ierr = 1 ! Singular DIIS matrix
    else
      res%ierr = 0 ! normal exit
      x = x_loc(1:self%maxdiis)
      select type (res)
      class is (scf_conv_interp_result)
        call move_alloc(from=x, to=res%coeffs)
      end select
    end if

    ! Compute DIIS error from the latest iteration
    diis_error = 0.0_dp
    do ifock = 1, self%dat%num_focks
      err => self%dat%get_err(-1, ifock)
      diis_error = max(diis_error, maxval(abs(err)))
    end do
    res%error = diis_error
  end subroutine cdiis_run

!> @brief Prepare C-DIIS subconverger to run
  subroutine cdiis_setup(self)
    class(cdiis_converger), intent(inout) :: self
    integer :: i, j, na, maxdiis, ifock, nfocks
    real(kind=dp) :: factor

    maxdiis = self%maxdiis
    na = self%dat%num_saved
    nfocks = self%dat%num_focks
    ! Factor to account RHF/UHF cases
    factor = 1.0_dp / nfocks

    if (self%last_setup > 1) then
      ! DIIS matrix is rather old, generate it from scratch
      self%old_dim = na
      self%a = 0.0_dp
      do ifock = 1, nfocks
        do i = 1, na
          do j = 1, i
            self%a(j,i) = self%a(j,i) + factor * dot_product( &
                            self%dat%get_err(j, ifock), &
                            self%dat%get_err(i, ifock))
          end do
        end do
      end do
    else if (self%last_setup == 1) then
      ! DIIS matrix is old by 1 iteration, just update it

      ! If the current number of iterations exceeds the dimension of A matrix:
      ! discard the data of oldest iteration by shifting the bottom-rigth square
      ! to the top-left corner
      if (self%old_dim >= maxdiis) then
        self%a(1:maxdiis-1, 1:maxdiis-1) = self%a(2:maxdiis, 2:maxdiis)
      end if
      self%old_dim = na
!     Compute new elements (`na`-th column)
      self%a(:, na) = 0.0_dp
      do ifock = 1, nfocks
        do i = 1, na
          self%a(i, na) = self%a(i, na) + factor * dot_product( &
                            self%dat%get_err(na, ifock), &
                            self%dat%get_err(i, ifock))
        end do
      end do
    end if

    ! DIIS matrix is already prepared nothing to do here
    self%last_setup = 0
  end subroutine cdiis_setup

!==============================================================================
! ediis_converger Methods
!==============================================================================

!> @brief Initialize E-DIIS subconverger
!> @detail This subroutine takes SCF converger driver as argument.
!>   It should be initialized and include all the required parameters.
!> @param[in] params current SCF converger driver
  subroutine ediis_init(self, params)
    class(ediis_converger), intent(inout) :: self
    type(scf_conv), target, intent(in) :: params

    call self%cdiis_converger%init(params)
    self%conv_name = 'E-DIIS'
    self%fun => ediis_fun
    allocate(self%b(self%maxdiis))
    allocate(self%xlog(self%maxdiis, self%maxdiis), source=0.0_dp)
  end subroutine ediis_init

!> @brief Finalize E-DIIS subconverger
  subroutine ediis_clean(self)
    class(ediis_converger), intent(inout) :: self

    call self%cdiis_converger%clean()
    if (allocated(self%b)) deallocate(self%b)
    if (allocated(self%xlog)) deallocate(self%xlog)
  end subroutine ediis_clean

!> @brief Prepare E-DIIS subconverger to run
  subroutine ediis_setup(self)
    class(ediis_converger), intent(inout) :: self
    integer :: i, j, na, maxdiis, ifock, nfocks
    real(kind=dp) :: factor

    maxdiis = self%maxdiis
    na = self%dat%num_saved
    nfocks = self%dat%num_focks
    ! Factor to account RHF/UHF cases
    factor = 1.0_dp / nfocks

    if (self%last_setup > 1) then
      ! DIIS matrix is rather old, generate it from scratch
      self%old_dim = na
      self%a = 0.0_dp
      do ifock = 1, nfocks
        do i = 1, na
          do j = 1, i-1
            self%a(j,i) = self%a(j,i) + factor * dot_product( &
                            self%dat%get_density(j, ifock) - self%dat%get_density(i, ifock), &
                            self%dat%get_fock(j, ifock) - self%dat%get_fock(i, ifock))
            self%a(i,j) = self%a(j,i)
          end do
        end do
      end do
    else if (self%last_setup == 1) then
      ! DIIS matrix is old by 1 iteration, just update it

      ! If the current number of iterations exceeds the dimension of A matrix:
      ! discard the data of oldest iteration by shifting the bottom-rigth square
      ! to the top-left corner
      if (self%old_dim >= maxdiis) then
        self%a(1:maxdiis-1, 1:maxdiis-1) = self%a(2:maxdiis, 2:maxdiis)
      end if
      self%old_dim = na
!     Compute new elements (`na`-th column)
      self%a(:, na) = 0.0_dp
      do ifock = 1, nfocks
        do j = 1, na
          self%a(j, na) = self%a(j, na) + factor * dot_product( &
                            self%dat%get_density(j, ifock) - self%dat%get_density(na, ifock), &
                            self%dat%get_fock(j, ifock) - self%dat%get_fock(na, ifock))
        end do
      end do
      self%a(na, 1:na-1) = self%a(1:na-1, na)
    end if

    do i = 1, na
      self%b(i) = self%dat%get_energy(i)
    end do

    ! DIIS matrix is already prepared nothing to do here
    self%last_setup = 0
  end subroutine ediis_setup

  !> @brief Computes the new guess to the SCF wavefunction using E-DIIS
  !> @param[out] res  results of the calculation
  subroutine ediis_run(self, res)
    use io_constants, only: iw
    use nlopt
    class(ediis_converger), target, intent(inout) :: self
    class(scf_conv_result), allocatable, intent(out) :: res
    real(kind=dp), parameter :: tol = 1.0e-5_dp
    real(kind=dp), parameter :: constrtol = 1.0e-8_dp
    real(kind=dp) :: minf, minf_min
    real(kind=dp), allocatable :: x(:), xmin(:)
    real(kind=dp), pointer :: err(:)
    real(kind=dp) :: diis_error
    integer(kind=4) :: ires
    integer :: na, i, j, ifock
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

    res%ierr = 3 ! Need to set up DIIS equations first
    if (self%last_setup /= 0) return

    na = self%dat%num_saved
    if (self%iter > self%maxdiis) then
      do i = 1, self%maxdiis-1
        self%xlog(:, i) = cshift(self%xlog(:, i+1), 1)
        self%xlog(i+1:, i) = 0.0_dp
      end do
    end if

    allocate(x(na), xmin(na))
    opt_global = 0
    minf_min = huge(1.0_dp)

    ! Initialize E-DIIS equation parameters for NLOpt
    t = ediis_opt_data(A=self%a(1:na, 1:na), b=self%b(1:na), fun=self%fun)
    ! Set up the Improved Stochastic Ranking Evolution Strategy
    ! It will run coarse global optimization, which will be further refined via L-BFGS
    call nlo_create(opt_global, NLOPT_GN_ISRES, na)
    ! Max. number of calls to the objective function
    call nlo_set_maxeval(ires, opt_global, 100*(na+1))
    ! Relative convergence tolerance for arguments
    call nlo_set_xtol_rel(ires, opt_global, 1.0e-2_dp)
    ! Absolute convergence tolerance for function value
    call nlo_set_ftol_abs(ires, opt_global, 1.0e-4_dp)
    ! Relative convergence tolerance for function value
    call nlo_set_ftol_rel(ires, opt_global, 1.0e-4_dp)
    ! Sum of coeffs equal to 1
    call nlo_add_equality_constraint(ires, opt_global, eadiis_constraints, 0, constrtol)
    ! 0 <= c_i <= 1
    call nlo_set_lower_bounds1(ires, opt_global, 0.0_dp)
    call nlo_set_upper_bounds1(ires, opt_global, 1.0_dp)
    ! Objective function
    call nlo_set_min_objective(ires, opt_global, eadiis_fun, t)

    x = 1.0_dp / na
    call nlo_optimize(ires, opt_global, x(:na), minf)
    call nlo_destroy(opt_global)

    ! Refine the results of global optimization via L-BFGS
    opt_lbfgs = 0
    call nlo_create(opt_lbfgs, NLOPT_LD_LBFGS, na)
    call nlo_set_xtol_rel(ires, opt_lbfgs, tol)
    call nlo_set_ftol_abs(ires, opt_lbfgs, tol*tol)
    call nlo_set_ftol_rel(ires, opt_lbfgs, tol*tol)
    ! Here, the modified E-DIIS equations are used, because L-BFGS does not support
    ! equality constraints
    ! They utilize the following variable substitution:
    ! c_i = t_i^2/\sum_i{t_i^2}
    call nlo_set_min_objective(ires, opt_lbfgs, eadiis_objective, t)
    call nlo_optimize(ires, opt_lbfgs, x(:na), minf)
    ! Because we used modified E-DIIS equations, we need to compute
    ! coefficients `c` from `t`:
    x = x**2 / sum(x**2)
    ! Get prediction of the new SCF energy
    call eadiis_fun(minf, int(na, 4), x, x, int(0,4), t)

    if (ires < 0) then
      if (self%verbose > 2) write(iw, '(10X,"*** nlopt0 failed:",I4," ***")') ires
    elseif (minf < minf_min) then
      is_a_repeat = any([(norm2(self%xlog(:na, j) - x(:na)) < 1.0e-4_dp, &
                          j = 1, min(self%iter, self%maxdiis))]) .or. &
                    any(1.0_dp - x(1:na-1) < 1.0e-4_dp)
      if (.not. is_a_repeat) then
        minf_min = minf
        xmin = x
        if (self%verbose > 2) then
          write(iw, '(A,*(F15.6))') 'nlopt0: improving x at ', xmin(:)
          write(iw, '(A,*(F15.6))') 'nlopt0: improved val = ', minf
        end if
      elseif (self%verbose > 2) then
        write(iw, '(A,*(F15.6))') 'nlopt0: found rep at ', x(:)
        write(iw, '(A,*(F15.6))') 'nlopt0: rep val = ', minf
      end if
    else
      if (self%verbose > 2) then
        write(iw, '(A,*(F15.6))') 'nlopt0: found min at ', x(:)
        write(iw, '(A,*(F15.6))') 'nlopt0: min val = ', minf
      end if
    end if

    ! If no solution found, try the alternative:
    ! Start from the trivial guess [x(1:n-1)=0, x(n) = 1]
    ! then run two L-BFGS iterations and average with
    ! [x(1:i-i), x(i+1:n) = 0, x(i) = 1] vector and run few L-BFGS steps again
    ! for all [ i = n-1, 1 ] and then [i = 1, n]
    if (minf_min > 1.0e99_dp) then
      call nlo_set_maxeval(ires, opt_lbfgs, 2)
      call nlo_set_xtol_rel(ires, opt_lbfgs, 0.1_dp)
      xmin = 0.0_dp
      xmin(na) = 1.0_dp
      do i = na-1, 1, -1
        x = 0.0_dp
        x(i) = 1.0_dp
        xmin = (xmin + x) / (1 + sum(x))
        call nlo_optimize(ires, opt_lbfgs, xmin, minf)
        if (ires < 0) then
          if (self%verbose > 2) then
            write(iw, '(10X,"*** nlopt2 failed:",I4," ***")') ires
          end if
          exit
        end if
        xmin = xmin**2 / sum(xmin**2)
      end do
      if (ires >= 0) then
        do i = 1, na
          x = 0.0_dp
          x(i) = 1.0_dp
          xmin = (xmin + x) / (1 + sum(x))
          call nlo_optimize(ires, opt_lbfgs, xmin, minf)
          if (ires < 0) then
            if (self%verbose > 2) then
              write(iw, '(10X,"*** nlopt2 failed:",I4," ***")') ires
            end if
            exit
          end if
          xmin = xmin**2 / sum(xmin**2)
        end do
      end if
      ! If still no success, the default is minimum energy + small contribution from others:
      if (ires < 0) then
        xmin = 1.0_dp
        xmin(minloc(self%b(1:na-1))) = 10.0_dp
        xmin = xmin / sum(xmin)
        if (self%verbose > 2) write(iw, *) 'nlopt2: unoptimal default'
      end if
      call eadiis_fun(minf_min, int(na, 4), xmin, xmin, int(0,4), t)
    end if

    minf = minf_min
    self%xlog(:na, min(self%iter, self%maxdiis)) = xmin(1:na)
    call nlo_destroy(opt_lbfgs)

    res%ierr = 0
    select type (res)
    class is (scf_conv_interp_result)
      call move_alloc(from=xmin, to=res%coeffs)
    end select

    ! compute diis error from the latest iteration
    diis_error = 0.0_dp
    do ifock = 1, self%dat%num_focks
      err => self%dat%get_err(-1, ifock)
      diis_error = max(diis_error, maxval(abs(err)))
    end do
    res%error = diis_error
  end subroutine ediis_run

!==============================================================================
! adiis_converger Methods
!==============================================================================

  !> @brief Initialize A-DIIS subconverger
  !> @detail This subroutine takes SCF converger driver as argument.
  !>   It should be initialized and include all the required parameters.
  !> @param[in] params current SCF converger driver
  subroutine adiis_init(self, params)
    class(adiis_converger), intent(inout) :: self
    type(scf_conv), target, intent(in) :: params

    call self%cdiis_converger%init(params)
    self%conv_name = 'A-DIIS'
    self%fun => adiis_fun
    allocate(self%b(self%maxdiis))
    allocate(self%xlog(self%maxdiis, self%maxdiis), source=0.0_dp)
  end subroutine adiis_init

  !> @brief Prepare A-DIIS subconverger to run
  subroutine adiis_setup(self)
    class(adiis_converger), intent(inout) :: self
    integer :: i, j, na, maxdiis, ifock, nfocks
    real(kind=dp) :: factor

    maxdiis = self%maxdiis
    na = self%dat%num_saved
    nfocks = self%dat%num_focks
    ! Factor to account RHF/UHF cases
    factor = 1.0_dp / nfocks

    self%old_dim = na
    self%a = 0.0_dp
    self%b = 0.0_dp
    do ifock = 1, nfocks
      do i = 1, na
        do j = 1, na
          self%a(j,i) = self%a(j,i) + 0.5_dp * factor * dot_product( &
                          self%dat%get_density(j, ifock) - self%dat%get_density(-1, ifock), &
                          self%dat%get_fock(j, ifock) - self%dat%get_fock(-1, ifock))
        end do
        self%b(i) = self%b(i) + factor * dot_product( &
                        self%dat%get_density(i, ifock) - self%dat%get_density(-1, ifock), &
                        self%dat%get_fock(-1, ifock))
      end do
    end do

    ! DIIS matrix is already prepared nothing to do here
    self%last_setup = 0
  end subroutine adiis_setup

!==============================================================================
! Optimization Helper Routines
!==============================================================================

  !> @brief Modified E/A-DIIS objective function wrapper, which allows to use unconstrained optimization
  !> @details The following variable substitution is used:
  !>          \f$ c_i = t_i^2 / \sum_i{t_i^2} \f$
  !> @note This is standard interface to work with NLOpt library
  !> @param[out] val Function value
  !> @param[in] n Dimension of the problem
  !> @param[in] t Vector of arguments
  !> @param[out] grad Vector of function gradient
  !> @param[in] need_gradient Flag to turn on computing gradient, 0 - gradient not computed
  !> @param[in] d Datatype storing function parameters
  subroutine eadiis_objective(val, n, t, grad, need_gradient, d)
    integer(kind=4), intent(in) :: n, need_gradient
    real(kind=8) :: val, t(n), grad(n)
    type(ediis_opt_data), intent(in) :: d
    real(kind=8) :: x(n), tnorm, jac(n,n)
    integer :: i, j

    tnorm = 1.0_dp / sum(t(1:n)**2)
    x = (t(1:n)**2) * tnorm
    call d%fun(val, n, x, grad, need_gradient, d)

    if (need_gradient /= 0) then
      jac = 0.0_dp
      do i = 1, n
        jac(i,i) = 1.0_dp
        do j = 1, n
          jac(j,i) = 2.0_dp * tnorm * t(j) * (jac(j,i) - x(i))
        end do
      end do
      grad(1:n) = matmul(jac, grad(1:n))
    end if
  end subroutine eadiis_objective

  !> @brief Non-modified E/A-DIIS objective function wrapper
  !> @note This is standard interface to work with NLOpt library
  !> @param[out] val Function value
  !> @param[in] n Dimension of the problem
  !> @param[in] t Vector of arguments
  !> @param[out] grad Vector of function gradient
  !> @param[in] need_gradient Flag to turn on computing gradient, 0 - gradient not computed
  !> @param[in] d Datatype storing function parameters
  subroutine eadiis_fun(val, n, x, grad, need_gradient, d)
    integer(kind=4), intent(in) :: n, need_gradient
    real(kind=8) :: val, x(n), grad(n)
    type(ediis_opt_data), intent(in) :: d

    call d%fun(val, n, x, grad, need_gradient, d)
  end subroutine eadiis_fun

  !> @brief E-DIIS objective function calculation
  !> @note This is standard interface to work with NLOpt library
  !> @param[out] val Function value
  !> @param[in] n Dimension of the problem
  !> @param[in] t Vector of arguments
  !> @param[out] grad Vector of function gradient
  !> @param[in] need_gradient Flag to turn on computing gradient, 0 - gradient not computed
  !> @param[in] d Datatype storing function parameters
  subroutine ediis_fun(val, n, x, grad, need_gradient, d)
    integer(kind=4), intent(in) :: n, need_gradient
    real(kind=8) :: val, x(*), grad(*)
    type(ediis_opt_data), intent(in) :: d

    if (need_gradient /= 0) then
      grad(1:n) = d%b(1:n) - matmul(d%A(1:n, 1:n), x(1:n))
    end if
    val = dot_product(x(1:n), d%b(1:n)) - &
          0.5_dp * dot_product(x(1:n), matmul(d%A(1:n, 1:n), x(1:n)))
  end subroutine ediis_fun

  !> @brief E/A-DIIS constraints
  !> @note This is standard interface to work with NLOpt library
  !> @param[out] val Function value
  !> @param[in] n Dimension of the problem
  !> @param[in] t Vector of arguments
  !> @param[out] grad Vector of function gradient
  !> @param[in] need_gradient Flag to turn on computing gradient, 0 - gradient not computed
  !> @param[in] d Datatype storing function parameters
  subroutine eadiis_constraints(val, n, x, grad, need_gradient, d)
    integer(kind=4), intent(in) :: n, need_gradient
    real(kind=8) :: val, x(n), grad(n)
    class(ediis_converger), intent(in) :: d

    if (need_gradient /= 0) grad = 1.0_dp
    val = sum(x) - 1.0_dp
  end subroutine eadiis_constraints

  !> @brief A-DIIS objective function calculation
  !> @note This is standard interface to work with NLOpt library
  !> @param[out] val Function value
  !> @param[in] n Dimension of the problem
  !> @param[in] t Vector of arguments
  !> @param[out] grad Vector of function gradient
  !> @param[in] need_gradient Flag to turn on computing gradient, 0 - gradient not computed
  !> @param[in] d Datatype storing function parameters
  subroutine adiis_fun(val, n, x, grad, need_gradient, d)
    integer(kind=4), intent(in) :: n, need_gradient
    real(kind=8) :: val, x(*), grad(*)
    type(ediis_opt_data), intent(in) :: d

    if (need_gradient /= 0) then
      grad(1:n) = 2.0_dp * d%b(1:n) + matmul(d%A(1:n, 1:n), x(1:n)) + &
                  matmul(x(1:n), d%A(1:n, 1:n))
    end if
    val = d%b(n) + 2.0_dp * dot_product(x(1:n), d%b(1:n)) + &
          dot_product(x(1:n), matmul(d%A(1:n, 1:n), x(1:n)))
  end subroutine adiis_fun

!==============================================================================
! Debug Printing (optional)
!==============================================================================

!> @brief Debug printing of the DIIS equation data
  subroutine diis_print_equation(self)
    use printing, only: print_square
    use io_constants, only: iw
    class(cdiis_converger), intent(in) :: self
    integer :: num_saved

    num_saved = self%dat%num_saved
    write(iw, '("(dbg) --------------------------------------------------")')
    write(iw, '("(dbg) DIIS iteration / max.dim. : ",G0," / ",G0)') self%iter, self%maxdiis
    write(iw, '("(dbg)",I4," Fock sets stored")') num_saved
    write(iw, '("(dbg) --------------------------------------------------")')
    write(iw, '("(dbg) Current DIIS equation matrix:")')
    call print_square(self%a, num_saved, num_saved, ubound(self%a,1), tag='(dbg)')
    write(iw, '("(dbg) --------------------------------------------------")')
    select type(self)
    class is (ediis_converger)
      write(iw, '("(dbg) Current DIIS equation vector:")')
      write(iw, '(*(ES15.7,","))') self%b(1:num_saved)
    end select
  end subroutine diis_print_equation

!==============================================================================
! soscf_converger Methods
!==============================================================================

  !> @brief Initialize the SOSCF converger
  !> @param[inout] self The SOSCF converger instance
  !> @param[in] params SCF convergence parameters from the driver
  subroutine soscf_init(self, params)
    class(soscf_converger), intent(inout) :: self
    type(scf_conv), target, intent(in) :: params

    integer :: nvec, m_max, istat, nocc_a, nocc_b

    ! Initialize base class
    call self%subconverger_init(params)
    self%conv_name = 'SOSCF'
    self%nfocks    = params%dat%num_focks
    self%verbose   = params%verbose
    self%nbf       = params%dat%ldim
    self%nbf_tri   = self%nbf * (self%nbf + 1) / 2
    self%m_max     = params%dat%num_slots
    self%m_history = 0
    self%dat       => params%dat
    self%overlap   => params%overlap
    self%overlap_invsqrt => params%overlap_sqrt
    self%first_macro = .true.

    ! Allocate working matrices
    if (.not. allocated(self%rho_history)) &
      allocate(self%rho_history(self%m_max), stat=istat, source=0.0_dp)
    if (.not. allocated(self%work_1)) &
      allocate(self%work_1(self%nbf, self%nbf), stat=istat, source=0.0_dp)
    if (.not. allocated(self%work_2)) &
      allocate(self%work_2(self%nbf, self%nbf), stat=istat, source=0.0_dp)

  end subroutine soscf_init

  !> @brief Clean up SOSCF converger
  !> @param[inout] self The SOSCF converger instance
  subroutine soscf_clean(self)
    class(soscf_converger), intent(inout) :: self
    if (allocated(self%s_history)) deallocate(self%s_history)
    if (allocated(self%rho_history)) deallocate(self%rho_history)
    if (allocated(self%y_history)) deallocate(self%y_history)
    if (allocated(self%grad_prev))    deallocate(self%grad_prev)
    if (allocated(self%x_prev))    deallocate(self%x_prev)
    if (allocated(self%h_inv))     deallocate(self%h_inv)
    if (allocated(self%work_1)) deallocate(self%work_1)
    if (allocated(self%work_2)) deallocate(self%work_2)
    if (allocated(self%mo_a)) deallocate(self%mo_a)
    if (allocated(self%mo_b)) deallocate(self%mo_b)
    if (allocated(self%dens_a)) deallocate(self%dens_a)
    if (allocated(self%dens_b)) deallocate(self%dens_b)
    call self%subconverger_clean()
  end subroutine soscf_clean

  !> @brief Setup SOSCF for the current iteration
  !> @param[inout] self The SOSCF converger instance
subroutine soscf_setup(self)

  class(soscf_converger), intent(inout) :: self

  ! Data verification on each setup call
  self%last_setup = 0
end subroutine soscf_setup

  !> @brief Run the SOSCF convergence step
  !> @param[inout] self The SOSCF converger instance
  !> @param[out] res Convergence result
subroutine soscf_run(self, res)
  use mathlib, only: unpack_matrix
  use oqp_linalg

  class(soscf_converger), target, intent(inout) :: self
  class(scf_conv_result), allocatable, intent(out) :: res

  real(kind=dp) :: soscf_error

  ! Allocate result as SOSCF-specific type
  allocate(scf_conv_soscf_result :: res)
  res%dat => self%dat
  res%active_converger_name = self%conv_name

  res%ierr = 0

  ! Main micro-iteration loop
! write(iw,*) "check max_iter=",self%max_iter

  soscf_error = 1.0e-2_dp  ! Dummy
  res%error = soscf_error

end subroutine soscf_run

  !> @brief Compute the orbital gradient for SOSCF
  !> @param[inout] self The SOSCF converger instance
  !> @param[out] gradient The computed gradient vector
subroutine soscf_compute_gradient(self, gradient)
  use mathlib, only: unpack_matrix, orthogonal_transform_sym

  class(soscf_converger), intent(inout) :: self
  real(kind=dp), intent(out) :: gradient(:)

end subroutine soscf_compute_gradient

  !> @brief Update the inverse Hessian and compute displacement using BFGS
  !> @param[inout] self The SOSCF converger instance
  !> @param[in] gradient Current orbital gradient
  !> @param[out] displacement New displacement vector
subroutine soscf_update_hessian(self, gradient, displacement)

  class(soscf_converger), intent(inout) :: self
  real(kind=dp), intent(in) :: gradient(:)
  real(kind=dp), intent(out) :: displacement(:)


end subroutine soscf_update_hessian

end module scf_converger
