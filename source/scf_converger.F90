!===============================================================================
! MODULE: scf_converger
!===============================================================================
!
! DESCRIPTION:
!   The scf_converger module implements a framework for managing convergence
!   in Self-Consistent Field (SCF) calculations. It provides a main driver
!   type `scf_conv` to coordinate multiple convergence methods, including
!   variants of Direct Inversion in the Iterative Subspace (DIIS) and
!   Second-Order SCF (SOSCF), storing iteration data in a ring buffer structure.
!
! MEMBERS:
!   - conv_none  [INTEGER]: Constant (1) for no convergence method.
!   - conv_cdiis [INTEGER]: Constant (2) for Commutator DIIS.
!   - conv_ediis [INTEGER]: Constant (3) for Energy DIIS.
!   - conv_adiis [INTEGER]: Constant (4) for Approximate DIIS.
!   - conv_soscf [INTEGER]: Constant (5) for Second-Order SCF.
!   - conv_state_not_initialized [INTEGER]: Constant (0) for uninitialized state.
!   - conv_state_initialized     [INTEGER]: Constant (1) for initialized state.
!   - conv_name_maxlen           [INTEGER]: Constant (32) for maximum length of converger names.
!
! DEPENDENCIES:
!   - precision: Provides `dp` for double precision real numbers.
!   - io_constants: Provides `iw` for output unit.
!   - messages: Provides `show_message` and `with_abort` for error handling.
!   - mathlib: Provides functions for matrix operations:
!     `traceprod_sym_packed`: Calculates trace of product of two packed symmetric matrices.
!     `orb_to_dens`: Converts orbitals to density matrix.
!     `antisymmetrize_matrix`: Makes matrix antisymmetric (A = A - A^T).
!     `unpack_matrix`: Converts packed triangular to full square matrix.
!     `pack_matrix`: Converts full square to packed triangular matrix.
!   - `scf_addons`: Provides Pseudo-Fractional Occupation Number (pFON) functionality.
!
! DATA MANAGEMENT:
!   The module uses a ring buffer system implemented in the `converger_data` type
!   for efficiently storing SCF iteration history. This designimproves memory efficiency
!   and data locality compared to traditional vector growth approaches.
!   The buffer system supports both forward and backward indices for intuitive access
!   (1 = oldest, -1 = latest).
!
! NOTES:
!   - The module defines public types: `scf_conv`, `scf_conv_result`, and `soscf_converger`.
!   - All other entities are private unless explicitly made public.
!   - The code supports handling of Fock matrices for RHF/ROHF (1 matrix) and UHF (2 matrices).
!   - ROHF Fock matrix has already gone through the Guest-Saunders transformation.
!   - The module supports integration with the pFON (pseudo-Fractional Occupation Number)
!     methodology.
!
! HISTORY:
!   - [pre-2022] Initial Development - Vladimir Mironov
!     Established the core framework for SCF convergence management,
!     including the main driver `scf_conv` and initial DIIS methods.
!   - [2025] Ring Buffer - Konstantin Komarov
!     Implemented the ring buffer structure via the `converger_data` type,
!     improving memory efficiency for storing SCF iteration history.
!   - [2025] SOSCF Converger - Konstantin Komarov
!     Added the Second-Order SCF (SOSCF) method through the `soscf_converger` type,
!     enabling convergence using orbital optimization techniques.
!
!===============================================================================

!===============================================================================
! TYPE: scf_conv - MAIN DRIVER TYPE
!===============================================================================
!
! DESCRIPTION:
!   The `scf_conv` type is the main driver for SCF convergence. It manages a
!   ring buffer of iteration data via `converger_data`, maintains an array of
!   subconvergers, and selects a convergence method based on the current error
!   and predefined thresholds.
!
! MEMBERS:
!   step          [INTEGER]: Current SCF iteration step.
!   overlap       [REAL(dp), POINTER]: Pointer to the full-format overlap matrix (S).
!   overlap_sqrt  [REAL(dp), POINTER]: Pointer to the full-format S^(1/2) matrix.
!   dat           [TYPE(converger_data)]: Stores SCF iteration history in a ring buffer.
!   sconv         [TYPE(subconverger_), ALLOCATABLE]: Array of subconverger instances.
!   thresholds    [REAL(dp), ALLOCATABLE]: Error thresholds for selecting subconvergers.
!   iter_space_size [INTEGER]: Size of the subconverger problem space, default is 10.
!   verbose       [INTEGER]: Verbosity level, default is 0.
!   state         [INTEGER]: Current state (0 = not initialized, 1 = initialized), default is 0.
!   current_error [REAL(dp)]: Maximum absolute value of the current DIIS error, default is 1.0e99_dp.
!
! METHODS:
!   init          - Initializes the driver with iteration parameters and subconvergers.
!   clean         - Deallocates resources and resets the driver.
!   add_data      - Adds data from a new SCF iteration to the ring buffer.
!   select_method - Selects the active subconverger based on the current error (private).
!   run           - Executes the selected subconverger and returns a result.
!   compute_error - Computes the current DIIS error (private).
!
! USAGE EXAMPLE IS BELOW.
!
!===============================================================================

!===============================================================================
! TYPE: converger_data - DATA MANAGEMENT TYPE
!===============================================================================
!
! DESCRIPTION:
!   The `converger_data` type implements a ring buffer for efficient storage
!   of SCF iteration history. It manages data locality and memory usage by
!   maintaining a fixed-size circular buffer rather than continuously growing
!   arrays, improving performance for calculations with many iterations.
!
! MEMBERS:
!   slot          [INTEGER]: Current slot in the ring buffer (1-based index).
!   num_saved     [INTEGER]: Number of occupied slots in the buffer.
!   num_slots     [INTEGER]: Total capacity of the buffer.
!   num_focks     [INTEGER]: Number of Fock matrices per iteration (1 for RHF/ROHF, 2 for UHF).
!   ldim          [INTEGER]: Size of square matrices (nbf).
!   nelec_a       [INTEGER]: Number of alpha electrons.
!   nelec_b       [INTEGER]: Number of beta electrons.
!   buffer        [TYPE(scf_data_t), ALLOCATABLE]: Array of SCF data containers.
!
! METHODS:
!   init         - Initializes the ring buffer with specified dimensions.
!   clean        - Deallocates buffer resources.
!   next_slot    - Advances to the next slot in the buffer (private).
!   discard_last - Removes the most recent iteration data (private).
!   put          - Stores new SCF data in the current slot (private).
!   get_fock     - Retrieves Fock matrix for specific iteration (private).
!   get_mo_a     - Retrieves alpha MO coefficients (private).
!   get_mo_b     - Retrieves beta MO coefficients (private).
!   get_mo_e_a   - Retrieves alpha MO energies (private).
!   get_mo_e_b   - Retrieves beta MO energies (private).
!   get_pfon     - Retrieves pFON object pointer (private).
!   get_density  - Retrieves density matrix (private).
!   get_err      - Retrieves error matrix (private).
!   get_energy   - Retrieves SCF energy (private).
!   get_slot     - Computes slot index from iteration number (private).
!
! USAGE NOTES:
!   - Data access supports both positive indices (1 = oldest) and
!     negative indices (-1 = latest) for intuitive programming.
!   - Each slot contains complete SCF iteration data including Fock matrices,
!     density matrices, MO coefficients, energies, and error matrices.
!   - Integration with pFON allows storage of fractional occupation data
!     when dealing with challenging electronic structures.
!
!===============================================================================

!===============================================================================
! TYPE: subconverger - ABSTRACT BASE TYPE
!===============================================================================
!
! DESCRIPTION:
!   The `subconverger` type is an abstract base type defining the interface
!   for all convergence methods used by `scf_conv`. Each subtype must implement
!   the deferred procedures `run` and `setup`.
!
! MEMBERS:
!   last_setup [INTEGER]: Number of SCF iterations since last setup, default is 1024.
!   iter       [INTEGER]: Number of iterations passed, default is 0.
!   conv_name  [CHARACTER(len=32)]: Name of the converger, default is an empty string.
!   dat        [TYPE(converger_data), POINTER]: Pointer to SCF data.
!
! METHODS:
!   init  - Initializes the subconverger (implemented as `subconverger_init`).
!   clean - Deallocates resources (implemented as `subconverger_clean`).
!   run   - Executes the convergence method, returns a `scf_conv_result` (deferred).
!   setup - Prepares the subconverger with new data (deferred).
!
!===============================================================================

!===============================================================================
! TYPE: noconv_converger - DUMMY CONVERGER
!===============================================================================
!
! DESCRIPTION:
!   The `noconv_converger` type extends `subconverger` and implements a dummy
!   convergence method that returns the latest iteration data without modification.
!
! METHODS:
!   init  - Initializes with the parent `scf_conv` instance.
!   run   - Returns the latest Fock matrix data in a `scf_conv_result`.
!   setup - Does nothing (empty implementation).
!
!===============================================================================

!===============================================================================
! TYPE: cdiis_converger - COMMUTATOR DIIS
!===============================================================================
!
! DESCRIPTION:
!   The `cdiis_converger` type extends `subconverger` and implements the
!   Commutator DIIS (C-DIIS) method, solving a constrained minimization problem
!   min { Ax, \sum_i x_i = 1 } where A_{ij} = Tr([F_i D_i S_i], [F_j D_j S_j]).
!
! MEMBERS:
!   maxdiis [INTEGER]: Maximum number of DIIS vectors.
!   a       [REAL(dp), ALLOCATABLE]: C-DIIS A-matrix for error overlaps.
!   verbose [INTEGER]: Verbosity level, default is 0.
!   old_dim [INTEGER]: Dimension of A-matrix from the previous step, default is 0.
!
! METHODS:
!   init  - Initializes the C-DIIS instance.
!   clean - Deallocates the A-matrix.
!   run   - Solves the DIIS equations and returns a `scf_conv_interp_result`.
!   setup - Updates the A-matrix with new error data.
!
!===============================================================================

!===============================================================================
! TYPE: ediis_converger - ENERGY DIIS
!===============================================================================
!
! DESCRIPTION:
!   The `ediis_converger` type extends `cdiis_converger` and implements the
!   Energy DIIS (E-DIIS) method, optimizing the function
!   E_E-DIIS = \sum_i c_i E_i - 0.5 \sum_{i,j} c_i c_j Tr((D_i-D_j)(F_i-F_j))
!   with constraints \sum_i c_i = 1, c_i \geq 0.
!
! MEMBERS:
!   b    [REAL(dp), ALLOCATABLE]: Energy history from previous iterations.
!   xlog [REAL(dp), ALLOCATABLE]: History of interpolation coefficients.
!   t    [TYPE(ediis_opt_data)]: Optimization data structure for NLOpt.
!   fun  [PROCEDURE(eadiis_f), POINTER]: Pointer to the objective function.
!
! METHODS:
!   init  - Initializes the E-DIIS instance.
!   clean - Deallocates additional arrays (b, xlog).
!   setup - Prepares the optimization problem with energy and matrix data.
!   run   - Solves the optimization problem and returns a `scf_conv_interp_result`.
!
!===============================================================================

!===============================================================================
! TYPE: adiis_converger - APPROXIMATE DIIS
!===============================================================================
!
! DESCRIPTION:
!   The `adiis_converger` type extends `ediis_converger` and implements the
!   Approximate DIIS (A-DIIS) method, optimizing the function
!   E_A-DIIS = \sum_i c_i Tr((D_i-D_n)F_n) + 2 \sum_{i,j} c_i c_j Tr((D_i-D_n)(F_j-F_n))
!   with constraints \sum_i c_i = 1, c_i \geq 0.
!
! METHODS:
!   init  - Initializes the A-DIIS instance (inherits from `ediis_converger`).
!   setup - Prepares the modified optimization problem specific to A-DIIS.
!
! NOTES:
!   - Inherits all members and most methods from `ediis_converger`.
!   - Only `setup` is overridden to adjust the objective function.
!
!===============================================================================

!===============================================================================
! TYPE: soscf_converger - SECOND-ORDER SCF
!===============================================================================
!
! DESCRIPTION:
!   The `soscf_converger` type extends `subconverger` and implements the
!   Second-Order SCF (SOSCF) method, using orbital gradients and an L-BFGS
!   approximation of the Hessian for convergence.
!
! MEMBERS:
!   verbose         [INTEGER]: Verbosity level, default is 0.
!   nfocks          [INTEGER]: Number of Fock matrices (1 for RHF, 2 for UHF), default is 0.
!   nbf             [INTEGER]: Number of basis functions, default is 0.
!   nbf_tri         [INTEGER]: Triangular size of basis functions (nbf*(nbf+1)/2), default is 0.
!   scf_type        [INTEGER]: SCF type (1=RHF, 2=UHF, 3=ROHF), default is 0.
!   nocc_a          [INTEGER]: Number of occupied alpha orbitals, default is 0.
!   nocc_b          [INTEGER]: Number of occupied beta orbitals, default is 0.
!   nvec            [INTEGER]: Size of the gradient vector, default is 0.
!   overlap         [REAL(dp), POINTER]: Pointer to overlap matrix, null by default.
!   overlap_invsqrt [REAL(dp), POINTER]: Pointer to inverse square root of overlap, null by default.
!   max_iter        [INTEGER]: Maximum micro-iterations, default is 10.
!   min_iter        [INTEGER]: Minimum micro-iterations, default is 1.
!   hess_thresh     [REAL(dp)]: Orbital Hessian threshold, default is 1.0e-10_dp.
!   grad_thresh     [REAL(dp)]: Gradient threshold, default is 1.0e-5_dp.
!   level_shift     [REAL(dp)]: Level shifting parameter, default is 0.2_dp.
!   s_history       [REAL(dp), ALLOCATABLE]: L-BFGS step history (nvec, m_max).
!   rho_history     [REAL(dp), ALLOCATABLE]: L-BFGS curvature reciprocals (m_max).
!   y_history       [REAL(dp), ALLOCATABLE]: L-BFGS gradient difference history (nvec, m_max).
!   grad            [REAL(dp), ALLOCATABLE]: Current gradient (nvec).
!   step            [REAL(dp), ALLOCATABLE]: Current step (nvec).
!   grad_prev       [REAL(dp), ALLOCATABLE]: Previous gradient (nvec).
!   x_prev          [REAL(dp), ALLOCATABLE]: Previous rotation parameters (nvec).
!   h_inv           [REAL(dp), ALLOCATABLE]: Initial inverse Hessian diagonal (nvec).
!   work_1          [REAL(dp), ALLOCATABLE]: Working matrix (nbf, nbf).
!   work_2          [REAL(dp), ALLOCATABLE]: Working matrix (nbf, nbf).
!   work_3          [REAL(dp), ALLOCATABLE]: Working matrix (nbf, nbf).
!   mo_a            [REAL(dp), ALLOCATABLE]: Alpha MO coefficients (nbf, nbf).
!   mo_b            [REAL(dp), ALLOCATABLE]: Beta MO coefficients (nbf, nbf).
!   dens_a          [REAL(dp), ALLOCATABLE]: Alpha density matrix (nbf_tri).
!   dens_b          [REAL(dp), ALLOCATABLE]: Beta density matrix (nbf_tri).
!   m_max           [INTEGER]: Maximum number of stored history steps, default is 0.
!   m_history       [INTEGER]: Current number of stored history steps, default is 0.
!   first_macro     [LOGICAL]: Flag for first macro-iteration, default is .true.
!
! METHODS:
!   init           - Initializes the SOSCF instance.
!   clean          - Deallocates arrays and resets pointers.
!   setup          - Prepares data for SOSCF iterations.
!   run            - Executes SOSCF micro-iterations and returns a `scf_conv_soscf_result`.
!   init_hess_inv  - Computes the initial inverse Hessian diagonal.
!   calc_orb_grad  - Calculates the orbital gradient.
!   lbfgs_step     - Computes the orbital rotation step using L-BFGS.
!   line_search    - Performs a line search to optimize the step size.
!   rotate_orbs    - Applies the orbital rotation.
!
!===============================================================================

!===============================================================================
! TYPE: scf_conv_result - BASE RESULT TYPE
!===============================================================================
!
! DESCRIPTION:
!   The `scf_conv_result` type is a base type for returning results from
!   convergence methods executed by `scf_conv`. It provides default
!   implementations that do nothing, overridden by subtypes.
!
! MEMBERS:
!   ierr           [INTEGER]: Error status (0 = success, nonzero = failure), default is 5.
!   error          [REAL(dp)]: Current DIIS error magnitude, default is 1.0e99_dp.
!   dat            [TYPE(converger_data), POINTER]: Pointer to SCF data, null by default.
!   active_converger_name [CHARACTER(len=32)]: Name of the active converger, default is empty.
!
! METHODS:
!   get_error    - Returns the current error value.
!   get_fock     - Placeholder returning no data (overridden by subtypes).
!   get_density  - Placeholder returning no data (overridden by subtypes).
!   get_mo_a     - Placeholder returning no data (overridden by subtypes).
!   get_mo_b     - Placeholder returning no data (overridden by subtypes).
!   get_mo_e_a   - Placeholder returning no data (overridden by subtypes).
!   get_mo_e_b   - Placeholder returning no data (overridden by subtypes).
!
!===============================================================================

!===============================================================================
! TYPE: scf_conv_interp_result - INTERPOLATION RESULT TYPE
!===============================================================================
!
! DESCRIPTION:
!   The `scf_conv_interp_result` type extends `scf_conv_result` and provides
!   results for interpolation-based methods (e.g., DIIS), computing updated
!   Fock or density matrices as F_n = \sum_i F_i * c_i or D_n = \sum_i D_i * c_i.
!
! MEMBERS:
!   coeffs [REAL(dp), ALLOCATABLE]: Interpolation coefficients c_i.
!   pstart [INTEGER]: Start index of the interpolation range.
!   pend   [INTEGER]: End index of the interpolation range.
!
! METHODS:
!   get_fock    - Returns the interpolated Fock matrix.
!   get_density - Returns the interpolated density matrix.
!   interpolate - Computes the interpolated matrix (Fock or density, private).
!
!===============================================================================

!===============================================================================
! TYPE: scf_conv_soscf_result - SOSCF RESULT TYPE
!===============================================================================
!
! DESCRIPTION:
!   The `scf_conv_soscf_result` type extends `scf_conv_result` and provides
!   results specific to the SOSCF method, returning the latest MO coefficients,
!   and MO energies.
!
! METHODS:
!   get_mo_a    - Returns the updated alpha MO coefficients.
!   get_mo_b    - Returns the updated beta MO coefficients.
!   get_mo_e_a  - Returns the updated alpha MO energies.
!   get_mo_e_b  - Returns the updated beta MO energies.
!
!===============================================================================

!===============================================================================
! EXAMPLE USAGE
!===============================================================================
!
!   use scf_converger
!   type(scf_conv) :: conv
!   class(scf_conv_result), allocatable :: res
!   ! Initialize converger
!   call conv%init(ldim=nbf,
!                  maxvec=15, &
!                  subconvergers=[conv_cdiis, conv_ediis], &
!                  thresholds=[2.0_dp, 1.0_dp], &
!                  overlap=s_mat,
!                  overlap_sqrt=s_sqrt, &
!                  num_focks=1,
!                  verbose=1)
!   ! In SCF loop:
!   do iter = 1, max_iter
!     ! Build Fock matrix...
!
!     ! Add data to converger
!     call conv%add_data(f=fock, dens=density, e=energy)
!
!     ! Run convergence step and get new Fock matrix
!     call conv%run(conv_res)
!     call conv_res%get_fock(matrix=fock, istat=stat)
!
!     ! Check convergence
!     diis_error = conv_res%get_error()
!     if (diis_error < thresh) exit
!   end do
!   ! Clean converger
!   call conv%clean()
!===============================================================================
! EXAMPLE USAGE
!===============================================================================
module scf_converger

  use precision, only: dp
  use io_constants, only: iw
  use mathlib, only: traceprod_sym_packed, orb_to_dens
  use messages, only: show_message, with_abort
  use mathlib, only: antisymmetrize_matrix, unpack_matrix, pack_matrix
  use scf_addons, only: pfon_t

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
    real(kind=dp)              :: energy = 0.0_dp !< SCF energy
    type(pfon_t), pointer      :: pfon_obj => null() !< Pseudo-Fractional Occupation Number (pFON) type
    real(kind=dp), pointer     :: occ_a(:) => null() !< Alpha orbital occupations (nbf)
    real(kind=dp), pointer     :: occ_b(:) => null() !< Beta orbital occupations (nbf)
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
    integer                   :: nelec_a = 0      !< Number of occupied alpha orbitals
    integer                   :: nelec_b = 0      !< Number of occupied beta orbitals
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
    procedure, private, pass :: get_pfon     => conv_data_get_pfon
    procedure, private, pass :: get_density  => conv_data_get_density
    procedure, private, pass :: get_err      => conv_data_get_err
    procedure, private, pass :: get_energy   => conv_data_get_energy
    procedure, private, pass :: get_slot     => conv_data_get_slot
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
    procedure, pass :: get_rms_grad => conv_result_dummy_get_rms_g
    procedure, pass :: get_rms_dp  => conv_result_dummy_get_rms_dp
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
    real(kind=dp) :: rms_grad = 1
    real(kind=dp) :: rms_dp = 1
  contains
    procedure, pass :: get_mo_a    => conv_result_soscf_get_mo_a
    procedure, pass :: get_mo_b    => conv_result_soscf_get_mo_b
    procedure, pass :: get_mo_e_a  => conv_result_soscf_get_mo_e_a
    procedure, pass :: get_mo_e_b  => conv_result_soscf_get_mo_e_b
    procedure, pass :: get_rms_grad => conv_result_soscf_get_rms_g
    procedure, pass :: get_rms_dp  => conv_result_soscf_get_rms_dp
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
    integer :: scf_type = 0
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

!    ! SOSCF parameters from scf_driver:
    real(kind=dp) :: hess_thresh = 1.0e-10_dp !< Orbital Hessian threshold
    real(kind=dp) :: grad_thresh = 1.0e-3_dp  !< Gradient threshold
    real(kind=dp) :: level_shift = 0.0_dp     !< Level shifting parameter
    real(kind=dp) :: rms_grad_prev = 0.0_dp   !< previous Gradient norm
    integer :: soscf_reset_mod = 1            !< Set the SOSCF Hessian reset mode.
    logical :: use_lineq = .false.            !< Use linear equations (vs BFGS)

    ! L-BFGS history
    real(kind=dp), allocatable :: s_history(:,:)  !< Step history (nvec, m_max)
    real(kind=dp), allocatable :: rho_history(:)  !< Curvature reciprocals (m_max)
    real(kind=dp), allocatable :: y_history(:,:)  !< Gradient difference history (nvec, m_max)
    real(kind=dp), allocatable :: upd_history(:,:)! h_inv*dgrad history (nvec, m_max)
    real(kind=dp), allocatable :: grad(:)         !< Gradient (nvec)
    real(kind=dp), allocatable :: step(:)         !< Step (nvec)
    real(kind=dp), allocatable :: grad_prev(:)    !< Previous gradient (nvec)
    real(kind=dp), allocatable :: x_prev(:)       !< Previous rotation parameters (nvec)
    real(kind=dp), allocatable :: h_inv(:)        !< Initial inverse Hessian diagonal (nvec)
    real(kind=dp), allocatable :: work_1(:,:)     !< Work matrix (nbf, nbf)
    real(kind=dp), allocatable :: work_2(:,:)     !< Work matrix (nbf, nbf)
    real(kind=dp), allocatable :: work_3(:,:)     !< Work matrix (nbf, nbf)
    real(kind=dp), allocatable :: mo_a(:,:)       !< MOs matrix (nbf, nbf)
    real(kind=dp), allocatable :: mo_b(:,:)       !< MOs matrix (nbf, nbf)
    real(kind=dp), allocatable :: dens_a(:)       !< MOs matrix (nbf_tri)
    real(kind=dp), allocatable :: dens_b(:)       !< MOs matrix (nbf_tri)
    real(kind=dp), allocatable :: dens_a_old(:)       !< MOs matrix (nbf_tri)
    real(kind=dp), allocatable :: dens_b_old(:)       !< MOs matrix (nbf_tri)
    integer :: m_max = 0                          !< Maximum number of stored history steps
    integer :: m_history = 0                      !< Number of stored history steps
    logical :: first_macro = .true.               !< Flag for first macro-iteration
  contains
    procedure, pass :: init  => soscf_init
    procedure, pass :: clean => soscf_clean
    procedure, pass :: setup => soscf_setup
    procedure, pass :: run   => soscf_run
    procedure, private, pass :: init_hess_inv => init_hess_inv
    procedure, private, pass :: calc_orb_grad => calc_orb_grad
    procedure, private, pass :: bfgs_step => bfgs_step
    procedure, private, pass :: rotate_orbs => rotate_orbs
    procedure, private, pass :: line_search => line_search
    procedure, private, pass :: rms_density => rms_density
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

    self%densities = 0.0_dp
    self%focks     = 0.0_dp
    self%errs      = 0.0_dp
    self%energy    = 0.0_dp
  end subroutine scf_data_init

  !> @brief Clean up an scf_data_t instance
  subroutine scf_data_clean(self)
    class(scf_data_t), intent(inout) :: self

    if (allocated(self%focks))     deallocate(self%focks)
    if (allocated(self%densities)) deallocate(self%densities)
    if (allocated(self%errs))      deallocate(self%errs)
    self%occ_a => null()
    self%occ_b => null()
    self%energy = 0.0_dp
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
  subroutine conv_data_init(self, ldim, nfocks, nslots, istat, nelec_a, nelec_b)
    class(converger_data), intent(inout) :: self
    integer, intent(in) :: ldim, nfocks, nslots
    integer, intent(out) :: istat
    integer, intent(in), optional :: nelec_a, nelec_b
    integer :: i

    istat = 0
    if (allocated(self%buffer)) call self%clean()
    self%ldim      = ldim
    self%num_focks = nfocks
    self%num_slots = nslots
    self%num_saved = 0
    self%slot      = 0
    if (present(nelec_a)) self%nelec_a = nelec_a
    if (present(nelec_b)) self%nelec_b = nelec_b
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
  subroutine conv_data_put(self, fock, dens, energy, mo_a, mo_b, &
                           mo_e_a, mo_e_b, pfon_obj)
    class(converger_data), intent(inout) :: self
    real(kind=dp), intent(in), optional :: fock(:,:)
    real(kind=dp), intent(in), optional :: dens(:,:)
    real(kind=dp), intent(in), optional :: energy
    real(kind=dp), intent(in), optional :: mo_a(:,:), mo_b(:,:)
    real(kind=dp), intent(in), optional :: mo_e_a(:), mo_e_b(:)
    type(pfon_t), optional, target, intent(in) :: pfon_obj

    integer :: slot, nbf

    nbf = self%ldim
    slot = self%slot
    if (present(fock)) self%buffer(slot)%focks = fock
    if (present(dens)) self%buffer(slot)%densities = dens
    if (present(energy)) self%buffer(slot)%energy = energy
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
    if (present(pfon_obj)) then
      self%buffer(slot)%pfon_obj => pfon_obj
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

  !> @brief Get pfon object for a specific iteration
  !> @param[in] n Slot ID (1 = oldest, -1 = latest)
  !> @return Pointer to the pfon object
  function conv_data_get_pfon(self, n) result(res)
    class(converger_data), intent(in), target :: self
    integer, intent(in) :: n
    type(pfon_t), pointer :: res
    integer :: slot

    slot = self%get_slot(n)
    res => self%buffer(slot)%pfon_obj
  end function conv_data_get_pfon

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

  function conv_result_dummy_get_rms_g(self) result(istat)
    class(scf_conv_result), intent(in) :: self
    real(kind=dp) :: istat

    istat = 0
  end function conv_result_dummy_get_rms_g

  function conv_result_dummy_get_rms_dp(self) result(istat)
    class(scf_conv_result), intent(in) :: self
    real(kind=dp) :: istat

    istat = 0
  end function conv_result_dummy_get_rms_dp

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

  function conv_result_soscf_get_rms_g(self) result(rms)
    class(scf_conv_soscf_result), intent(in) :: self
    real(kind=dp) :: rms

    rms = self%rms_grad
  end function conv_result_soscf_get_rms_g

  function conv_result_soscf_get_rms_dp(self) result(rms)
    class(scf_conv_soscf_result), intent(in) :: self
    real(kind=dp) :: rms

    rms = self%rms_dp
  end function conv_result_soscf_get_rms_dp


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
  subroutine scf_conv_init(self, ldim, nelec_a, nelec_b, maxvec, subconvergers, thresholds, &
                          overlap, overlap_sqrt, num_focks, scf_type ,verbose)
    class(scf_conv), intent(inout) :: self
    integer, intent(in) :: ldim
    integer, optional, intent(in) :: nelec_a
    integer, optional, intent(in) :: nelec_b
    integer, optional, intent(in) :: maxvec
    integer, optional, intent(in) :: subconvergers(:)
    real(kind=dp), optional, intent(in) :: thresholds(:)
    real(kind=dp), optional, target, intent(in) :: overlap(:,:), overlap_sqrt(:,:)
    integer, optional, intent(in) :: num_focks
    integer, optional, intent(in) :: verbose
    integer :: nfocks, istat, i
    integer, optional, intent(in) :: scf_type

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

    if (present(scf_type)) self%scf_type = scf_type

    self%verbose = 0
    if (present(verbose)) self%verbose = verbose

    self%iter_space_size = 15
    if (present(maxvec)) self%iter_space_size = maxvec

    self%step = 0

    if (present(nelec_a) .and. present(nelec_b)) then
      call self%dat%init(ldim, nfocks, self%iter_space_size, istat, nelec_a, nelec_b)
    else
      call self%dat%init(ldim, nfocks, self%iter_space_size, istat)
    end if
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
                               pfon)
    class(scf_conv), intent(inout) :: self
    real(kind=dp), intent(in), optional :: f(:,:)    ! Fock matrices
    real(kind=dp), intent(in), optional :: dens(:,:) ! Density matrices
    real(kind=dp), intent(in), optional :: e         ! SCF energy
    real(kind=dp), intent(in), optional :: mo_a(:,:) ! Alpha MO coefficients
    real(kind=dp), intent(in), optional :: mo_b(:,:) ! Beta MO coefficients
    real(kind=dp), intent(in), optional :: mo_e_a(:) ! Alpha MO energies
    real(kind=dp), intent(in), optional :: mo_e_b(:) ! Beta MO energies
    type(pfon_t), pointer, intent(in), optional :: pfon ! Pseudo-Fractional Occupation Number (pFON) object
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
    if (present(pfon)) call self%dat%put(pfon_obj=pfon)

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

  !> @brief Select subconverger basing on the current DIIS error value
  !> @param[in] error DIIS error value
  !> @return  pointer to selected converger
  function scf_conv_select(self, error) result(conv)
    implicit none
    class(scf_conv), target :: self
    real(kind=dp), intent(in) :: error
    class(subconverger), pointer :: conv
    integer :: i, nconv
    logical :: use_soscf_now

    nconv = ubound(self%thresholds, 1)
    do i = 0, nconv
      if (error > self%thresholds(i)) exit
    end do

    conv => self%sconv(min(i,nconv))%s
    use_soscf_now = .false.
    do i = 0, nconv
      select type (sc => self%sconv(i)%s)
      type is (soscf_converger)
        ! SOSCF exists in the set of convergers
        ! Check if we should use it based on current iteration
        use_soscf_now = .true.
        exit
      end select
    end do
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

    integer :: nvec, m_max, istat, nvir_a, nvir_b

    ! Initialize base class
    call self%subconverger_init(params)
    self%conv_name = 'SOSCF'
    self%nfocks    = params%dat%num_focks
    self%verbose   = params%verbose
    self%nbf       = params%dat%ldim
    self%nocc_a    = params%dat%nelec_a
    self%nocc_b    = params%dat%nelec_b
    self%nbf_tri   = self%nbf * (self%nbf + 1) / 2
    self%m_max     = params%dat%num_slots
    self%m_history = 0
    self%dat       => params%dat
    self%overlap   => params%overlap
    self%overlap_invsqrt => params%overlap_sqrt
    self%first_macro = .true.
    self%scf_type = params%scf_type

    nvir_a = self%nbf - self%nocc_a
    nvir_b = self%nbf - self%nocc_b

    ! Calculate gradient vector size (nvec)
    select case (self%scf_type)
    case (1)  ! RHF
      self%nvec = self%nocc_a * nvir_a
    case (2)  ! UHF
      self%nvec = self%nocc_a * nvir_a + self%nocc_b * nvir_b
    case (3)  ! ROHF
      self%nvec = self%nocc_b * nvir_b + (self%nocc_a - self%nocc_b) * nvir_a
    end select

    ! Allocate working matrices
    istat = 0
    if (.not. allocated(self%rho_history)) &
      allocate(self%rho_history(self%m_max), stat=istat, source=0.0_dp)
    if (.not. allocated(self%work_1)) &
      allocate(self%work_1(self%nbf, self%nbf), stat=istat, source=0.0_dp)
    if (.not. allocated(self%work_2)) &
      allocate(self%work_2(self%nbf, self%nbf), stat=istat, source=0.0_dp)
    if (.not. allocated(self%work_3)) &
      allocate(self%work_3(self%nbf, self%nbf), stat=istat, source=0.0_dp)
    ! Allocate L-BFGS history arrays
    if (.not. allocated(self%s_history)) &
      allocate(self%s_history(self%nvec, self%m_max), stat=istat, source=0.0_dp)
    if (.not. allocated(self%grad)) &
      allocate(self%grad(self%nvec), stat=istat, source=0.0_dp)
    if (.not. allocated(self%step)) &
      allocate(self%step(self%nvec), stat=istat, source=0.0_dp)
    if (.not. allocated(self%grad_prev)) &
      allocate(self%grad_prev(self%nvec), stat=istat, source=0.0_dp)
    if (.not. allocated(self%y_history)) &
      allocate(self%y_history(self%nvec, self%m_max), stat=istat, source=0.0_dp)
    if (.not. allocated(self%upd_history)) &
      allocate(self%upd_history(self%nvec, self%m_max), stat=istat, source=0.0_dp)
    if (.not. allocated(self%x_prev)) &
      allocate(self%x_prev(self%nvec), stat=istat, source=0.0_dp)
    if (.not. allocated(self%h_inv)) &
      allocate(self%h_inv(self%nvec), stat=istat, source=0.0_dp)

    ! Allocate SOSCF arrays
    if (.not. allocated(self%mo_a)) &
      allocate(self%mo_a(self%nbf, self%nbf), stat=istat, source=0.0_dp)
    if (.not. allocated(self%dens_a)) &
      allocate(self%dens_a(self%nbf_tri), stat=istat, source=0.0_dp)
    if (.not. allocated(self%dens_a_old)) &
      allocate(self%dens_a_old(self%nbf_tri), stat=istat, source=0.0_dp)
    if (self%scf_type > 1) then
      if (.not. allocated(self%mo_b)) &
        allocate(self%mo_b(self%nbf, self%nbf), stat=istat, source=0.0_dp)
      if (.not. allocated(self%dens_b)) &
        allocate(self%dens_b(self%nbf_tri), stat=istat, source=0.0_dp)
      if (.not. allocated(self%dens_b_old)) &
        allocate(self%dens_b_old(self%nbf_tri), stat=istat, source=0.0_dp)
    end if

    if (istat /= 0) then
      write(iw, '(A)') 'ERROR: Failed to allocate arrays in soscf_init'
      stop
    end if

  end subroutine soscf_init

  !> @brief Clean up SOSCF converger
  !> @param[inout] self The SOSCF converger instance
  subroutine soscf_clean(self)
    class(soscf_converger), intent(inout) :: self
    if (allocated(self%rho_history)) deallocate(self%rho_history)
    if (allocated(self%work_1)) deallocate(self%work_1)
    if (allocated(self%work_2)) deallocate(self%work_2)
    if (allocated(self%work_3)) deallocate(self%work_3)
    if (allocated(self%s_history)) deallocate(self%s_history)
    if (allocated(self%y_history)) deallocate(self%y_history)
    if (allocated(self%upd_history)) deallocate(self%upd_history)
    if (allocated(self%grad)) deallocate(self%grad)
    if (allocated(self%grad_prev)) deallocate(self%grad_prev)
    if (allocated(self%x_prev)) deallocate(self%x_prev)
    if (allocated(self%h_inv)) deallocate(self%h_inv)
    if (allocated(self%mo_a)) deallocate(self%mo_a)
    if (allocated(self%dens_a)) deallocate(self%dens_a)
    if (allocated(self%mo_b)) deallocate(self%mo_b)
    if (allocated(self%dens_b)) deallocate(self%dens_b)
    if (allocated(self%dens_a_old)) deallocate(self%dens_a_old)
    if (allocated(self%dens_b_old)) deallocate(self%dens_b_old)
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

    ! Local variables
    type(pfon_t), pointer :: pfon => null()
    real(kind=dp), pointer :: occ_a(:) => null()
    real(kind=dp), pointer :: occ_b(:) => null()
    real(kind=dp), pointer :: fock_ao_a(:), fock_ao_b(:)
    real(kind=dp), pointer :: mo_e_a(:), mo_e_b(:)
    real(kind=dp) :: grad_norm, alpha, sy, rms_dp
    real(kind=dp) :: grad_norm_ratio
    integer :: iter, istat
    integer :: i

    ! Allocate result object
    allocate(scf_conv_soscf_result :: res)
    res%dat => self%dat
    res%active_converger_name = self%conv_name

    if (self%last_setup /= 0) then
      if (self%verbose > 0) write(iw, '(A)') 'SOSCF: Setup not called, returning'
      return
    end if

    res%ierr = 0
    istat = 0


    ! --- Step 1: Extract current data from converger_data ---
    mo_e_a => self%dat%get_mo_e_a(-1)
    fock_ao_a => self%dat%get_fock(-1, 1)
    self%mo_a = self%dat%get_mo_a(-1)
    self%dens_a = self%dat%get_density(-1, 1)
    if (self%scf_type > 1) then
      mo_e_b => self%dat%get_mo_e_b(-1)
      fock_ao_b => self%dat%get_fock(-1, 2)
      self%mo_b = self%dat%get_mo_b(-1)
      self%dens_b = self%dat%get_density(-1, 2)
    end if
    if (self%scf_type == 1) then
      call self%rms_density( d_new_a = self%dens_a,       &
                             d_old_a = self%dens_a_old,   &
                             rms_dp  = rms_dp )
      self%dens_a_old = self%dens_a
    else
      call self%rms_density(self%dens_a, self%dens_b,            &
                              self%dens_a_old, self%dens_b_old, rms_dp)
      self%dens_a_old = self%dens_a
      self%dens_b_old = self%dens_b
    end if

    pfon => self%dat%get_pfon(-1)
    if (associated(pfon)) then
      ! Get occupation arrays directly from pfon
      occ_a => pfon%occ_a
      if (self%scf_type > 1) then
        occ_b => pfon%occ_b
      end if
    else
      ! Allocate local occupation arrays
      allocate(occ_a(self%nbf), source=0.0_dp)
      if (self%scf_type > 1) then
        allocate(occ_b(self%nbf), source=0.0_dp)
      end if

      ! Set initial occupation numbers based on SCF type
      select case (self%scf_type)
      case (1)
        occ_a(1:self%nocc_a) = 2.0_dp
      case (2)
        occ_a(1:self%nocc_a) = 1.0_dp
        occ_b(1:self%nocc_b) = 1.0_dp
      case (3)
        occ_a(1:self%nocc_b) = 2.0_dp              ! Closed shells
        occ_a(self%nocc_b+1:self%nocc_a) = 1.0_dp  ! Open shells
        occ_b(1:self%nocc_b) = 2.0_dp              ! Closed shells only
      end select
    end if

    ! --- Step 2: Initialize LBFGS history ---
    if (self%first_macro) then
      self%s_history = 0.0_dp
      self%y_history = 0.0_dp
      self%upd_history = 0.0_dp
      self%grad_prev = 0.0_dp
      self%grad = 0.0_dp
      self%step = 0.0_dp
      self%x_prev = 0.0_dp
      self%m_history = 0
      if (self%m_history == 0) &
        call self%init_hess_inv(mo_e_a, mo_e_b)
      self%first_macro = .false.
      if (self%verbose>1) then
        write(iw, '(A,E20.10)') 'DEBUG: soscf_run: Input: grad_thresh=', self%grad_thresh
        write(iw, '(A,E20.10)') 'DEBUG: soscf_run: Input: level_shift=', self%level_shift
      end if
    end if

    call self%calc_orb_grad(self%grad, fock_ao_a, fock_ao_b, self%mo_a, self%mo_b)
    grad_norm = sqrt(dot_product(self%grad, self%grad)/self%nvec)
    if (self%m_history == 0 )  self%rms_grad_prev = grad_norm
    if (self%verbose>1) &
      write(iw, '(A,I3,A,E20.10)') 'DEBUG soscf_run:1:calc_orb_grad: iter=', 0, ' grad_norm=', grad_norm

    if (self%verbose>1) &
      write(iw, '(A,E20.10)') 'DEBUG soscf_run: init: grad_thresh=', self%grad_thresh
    ! Check convergence
    if ( grad_norm < self%grad_thresh) then
      if (self%verbose>1) &
        write(iw, '(A,I3,A,E20.10)') 'DEBUG soscf_run: loop exit: iter=',iter, ' grad_norm=', grad_norm
    end if

    ! Compute step using J. Phys. Chem. 1992, 96, 9768-9774
    call self%bfgs_step(self%step)
    if (self%scf_type == 1) then
       call self%rotate_orbs(self%step, self%nocc_a, self%nocc_a, self%mo_a)
      if (associated(pfon)) then
        call compute_mo_energies(self, fock_ao_a, self%mo_a, mo_e_a, self%work_1, self%work_2)
        call pfon%compute_occupations(mo_e_a)
        call pfon%build_density(self%dens_a, self%mo_a, self%work_1, self%work_2)
      else
        call orb_to_dens(self%dens_a, self%mo_a, occ_a, self%nocc_a, self%nbf, self%nbf)
      end if
    elseif(self%scf_type == 2) then
      call self%rotate_orbs(self%step, self%nocc_a, self%nocc_a, self%mo_a)
      call self%rotate_orbs(self%step(self%nocc_a * (self%nbf - self%nocc_a) +1 : self%nvec)&
              , self%nocc_b, self%nocc_b, self%mo_b)
      if (associated(pfon)) then
        call compute_mo_energies(self, fock_ao_a, self%mo_a, mo_e_a, self%work_1, self%work_2)
        call compute_mo_energies(self, fock_ao_b, self%mo_b, mo_e_b, self%work_1, self%work_2)
        call pfon%compute_occupations(mo_e_a, mo_e_b)
        call pfon%build_density(self%dens_a, self%mo_a, self%work_1, self%work_2, self%dens_b, self%mo_b)
      else
        call orb_to_dens(self%dens_a, self%mo_a, occ_a, self%nocc_a, self%nbf, self%nbf)
        call orb_to_dens(self%dens_b, self%mo_b, occ_b, self%nocc_b, self%nbf, self%nbf)
      end if
    elseif(self%scf_type == 3) then
      call self%rotate_orbs(self%step, self%nocc_a, self%nocc_b, self%mo_a)
      self%mo_b(1:self%nbf, 1:self%nbf) = self%mo_a(1:self%nbf, 1:self%nbf)
      if (associated(pfon)) then
        call compute_mo_energies(self, fock_ao_a, self%mo_a, mo_e_a, self%work_1, self%work_2)
        call compute_mo_energies(self, fock_ao_b, self%mo_b, mo_e_b, self%work_1, self%work_2)
        call pfon%compute_occupations(mo_e_a, mo_e_b)
        call pfon%build_density(self%dens_a, self%mo_a, self%work_1, self%work_2, self%dens_b, self%mo_b)
      else
        call orb_to_dens(self%dens_a, self%mo_a, occ_a, self%nocc_a, self%nbf, self%nbf)
        call orb_to_dens(self%dens_b, self%mo_b, occ_b, self%nocc_b, self%nbf, self%nbf)
      end if
    end if
    self%m_history = self%m_history + 1
    self%grad_prev = self%grad

    ! --- Step 4: Update result and converger_data ---
    res%ierr = 0
    res%error = grad_norm
    select type (res)
    class is (scf_conv_soscf_result)
      res%rms_grad = res%error
      res%rms_dp = rms_dp
    end select

    ! Update MO coefficients and compute MO energies
    self%dat%buffer(self%dat%slot)%mo_a = self%mo_a
    if (self%scf_type == 2) self%dat%buffer(self%dat%slot)%mo_b = self%mo_b
    call compute_mo_energies(self, fock_ao_a, self%mo_a, &
                             self%dat%buffer(self%dat%slot)%mo_e_a, self%work_1, self%work_2)
    if (self%scf_type == 2) then
      call compute_mo_energies(self, fock_ao_b, self%mo_b, &
                               self%dat%buffer(self%dat%slot)%mo_e_b, self%work_1, self%work_2)
    end if

    if (self%soscf_reset_mod /= 0 .and. mod(self%m_history, self%soscf_reset_mod) == 0) then
      if (self%rms_grad_prev > 1.0d-12) then
        grad_norm_ratio = grad_norm / self%rms_grad_prev
        self%rms_grad_prev = grad_norm
        if (grad_norm_ratio > 0.95_dp) then
          write(iw, '(8X, "Resetting Hessian, gradient norm ratio = ", F8.5)') grad_norm_ratio
          self%m_history = 0
          self%first_macro = .true.
        end if
      end if
    end if


  end subroutine soscf_run

  subroutine rms_density(self, d_new_a, d_new_b,            &
                                d_old_a, d_old_b, rms_dp, max_dp)
     use, intrinsic :: iso_fortran_env, only : dp => real64
     implicit none
     class(soscf_converger), intent(in)       :: self
     real(dp), intent(in)                     :: d_new_a(:), d_old_a(:)
     real(dp), intent(in),  optional          :: d_new_b(:), d_old_b(:)
     real(dp), intent(out)                    :: rms_dp
     real(dp), intent(out), optional          :: max_dp

     integer  :: ntri, k
     real(dp) :: sum_sq, sum_sq_b
     real(dp) :: max_loc, max_loc_b
     real(dp) :: diff

     ntri = self%nbf_tri

     sum_sq  = 0.0_dp
     max_loc = 0.0_dp
  !$omp   parallel do default(shared) private(k,diff)             &
  !$omp&  reduction(+:sum_sq) reduction(max:max_loc)
     do k = 1, ntri
        diff    = d_new_a(k) - d_old_a(k)
        sum_sq  = sum_sq + diff*diff
        max_loc = max(max_loc, abs(diff))
     end do
  !$omp   end parallel do

     rms_dp = sqrt( sum_sq / real(ntri, dp) )
     if (present(max_dp)) max_dp = max_loc

     if (present(d_new_b) .and. present(d_old_b)) then
        sum_sq_b  = 0.0_dp
        max_loc_b = 0.0_dp
  !$omp   parallel do default(shared) private(k,diff)             &
  !$omp&  reduction(+:sum_sq_b) reduction(max:max_loc_b)
        do k = 1, ntri
           diff       = d_new_b(k) - d_old_b(k)
           sum_sq_b   = sum_sq_b + diff*diff
           max_loc_b  = max(max_loc_b, abs(diff))
        end do
  !$omp   end parallel do

        rms_dp = 0.5_dp * ( rms_dp + sqrt( sum_sq_b / real(ntri, dp) ) )
        if (present(max_dp)) max_dp = max(max_dp, max_loc_b)
     end if
  end subroutine rms_density

  !> @brief Computes the initial diagonal inverse Hessian for SOSCF
  !> @details Approximates the inverse Hessian diagonal using orbital energy differences,
  !>          adjusted for SCF type:
  !>          - RHF: h_inv = 0.25 / (ε_a - ε_i) for closed-shell.
  !>          - UHF: h_inv combines alpha and beta energy differences.
  !>          - ROHF: h_inv varies by orbital region (closed-virtual, open-virtual).
  !>          Includes level-shifting for stability.
  !> @param[out] self%h_inv Diagonal inverse Hessian (size depends on scf_type)
  !> @param[in] mo_e_a Alpha orbital energies (size: nbf)
  !> @param[in] mo_e_b Beta orbital energies (size: nbf, ignored for RHF/ROHF)
  subroutine init_hess_inv(self, mo_e_a, mo_e_b)
    implicit none
    class(soscf_converger) :: self
    real(kind=dp), pointer, intent(in) :: mo_e_a(:)
    real(kind=dp), pointer, intent(in) :: mo_e_b(:)

    real(kind=dp) :: diff, scale
    integer :: i, a, k, istart

    associate (nbf       => self%nbf, &
               nocc_a    => self%nocc_a, &
               nocc_b    => self%nocc_b, &
               scf_type  => self%scf_type, &
               lvl_shift => self%level_shift, &
               thresh    => self%hess_thresh)
      select case (scf_type)
      case (1) ! RHF: Closed-shell system
        ! Single set of orbitals, nocc_a = nocc_b, 4 * (ε_a - ε_i) scaling
        k = 0
        do i = 1, nocc_a
          if (i <= nocc_b ) then
            istart = nocc_b+1
            scale = 0.25_dp  ! Closed-virtual: like RHF
          else
            istart = nocc_a +1
            scale = 0.25_dp   ! Open-virtual: like UHF (singly occupied)
          end if
          do a= istart, nbf
            k = k +1
            diff = mo_e_a(a) - mo_e_a(i)
            if (abs(diff) < thresh) then
               diff = sign(thresh + lvl_shift, diff)
            end if
            self%h_inv(k) = scale /diff
          end do
        end do

      case (2) ! UHF: Unrestricted, separate alpha and beta orbitals
        ! Alpha occupied -> alpha virtual, beta occupied -> beta virtual
        k = 0
        ! Alpha rotations
        do i = 1, nocc_a
          do a= nocc_a+1, nbf
            k = k +1
            diff = mo_e_a(a) - mo_e_a(i)
            if (abs(diff) < thresh) then
               diff = sign(thresh + lvl_shift, diff)
            end if
            self%h_inv(k) = 0.5_dp /diff
          end do
        end do
        ! Beta rotations
        do i = 1, nocc_b
          do a= nocc_b+1, nbf
            k = k +1
            diff = mo_e_a(a) - mo_e_a(i)
            if (abs(diff) < thresh) then
               diff = sign(thresh + lvl_shift, diff)
            end if
            self%h_inv(k) = 0.5_dp /diff
          end do
        end do


      case (3) ! ROHF: Restricted open-shell
        ! Regions: closed (j <= nocc_b)
        !          open (nocc_b < j <= nocc_a)
        !          virtual (j > nocc_a)
        k = 0
        do i = 1, nocc_a
          if (i <= nocc_b ) then
            istart = nocc_b+1
            scale = 0.25_dp  ! Closed-virtual: like RHF
          else
            istart = nocc_a +1
            scale = 0.50_dp   ! Open-virtual: like UHF (singly occupied)
          end if
          do a= istart, nbf
            k = k +1
            if (i <= nocc_b) then
              diff = mo_e_a(a) - mo_e_a(i)
              if (abs(diff) < thresh) then
                 diff = sign(thresh + lvl_shift, diff)
              end if
            end if
            self%h_inv(k) = scale /diff
          end do
        end do
      end select
    end associate
  end subroutine init_hess_inv

  !> @brief Computes orbital gradient.
  !> @details Calculates the gradient for occupied-virtual orbital rotations:
  !>          - RHF: g(i,a) = 4 * F(i,a), single gradient.
  !>          - UHF: g_a(i,a) = 2 * F_a(i,a), g_b(i,a) = 2 * F_b(i,a), separate gradients.
  !>          - ROHF: g(c,v) = 2 * F_b(c,v) for doubly occ -> singly occ,
  !>                  g(o,v) = 2 * (F_a(o,v) + F_b(o,v)) doubly occ ->virt,
  !>                  g(o,v) = 2 * F_a(o,v) for singly occ -> virt
  !> @param[out] grad_a Alpha gradient for RHF, ROHF, and UHF
  !> @param[out] grad_b Beta gradient for UHF (optional, size: nocc_b * (nbf - nocc_b))
  !> @param[in] fock_a Packed alpha/single Fock matrix (size: nbf*(nbf+1)/2)
  !> @param[in] fock_b Packed beta Fock matrix for UHF (optional, size: nbf*(nbf+1)/2)
  subroutine calc_orb_grad(self, grad, fock_a, fock_b, mo_a, mo_b)
    implicit none
    class(soscf_converger)             :: self
    real(kind=dp), intent(out)         :: grad(:)
    real(kind=dp), pointer, intent(in) :: fock_a(:)
    real(kind=dp), pointer, intent(in) :: fock_b(:)
    real(kind=dp), intent(in)          :: mo_a(:,:)
    real(kind=dp), intent(in)          :: mo_b(:,:)

    integer :: i, a, k, nvir_a, nvir_b, istart

    if (.not. ASSOCIATED(fock_a)) &
      call show_message('Failed to use Fock array in calc_orb_grad', &
                        with_abort)

    associate (nbf    => self%nbf, &
               nocc_a => self%nocc_a, &
               nocc_b => self%nocc_b)

      nvir_a = nbf - nocc_a
      nvir_b = nbf - nocc_b
      grad = 0.0_dp

      select case (self%scf_type)
      case (1)  ! RHF
        ! Unpack AO-basis Fock matrix
        self%work_1 = 0.0_dp
        call unpack_matrix(fock_a, self%work_1)
        ! Convert AO Fock to MO Fock.
        ! self%work_2 = F_ao * C_virt
        call dgemm('N', 'N', nbf, nvir_a, nbf, &
                   1.0_dp, self%work_1, nbf, &
                           mo_a(:, nocc_a+1:), nbf, &
                   0.0_dp, self%work_2, nbf)
        ! self%work_1 = C_occ^T * self%work_2
        call dgemm('T', 'N', nocc_a, nvir_a, nbf, &
                   1.0_dp, mo_a(:, 1:nocc_a), nbf, &
                           self%work_2, nbf, &
                   0.0_dp, self%work_1(1:nocc_a, 1:nvir_a), nocc_a)
        k = 0
        do i = 1, nocc_a
          if (i <= nocc_b ) then
            istart = nbf - nocc_b
          else
            istart = nbf - nocc_a
          end if
          do a= 1, istart
            k = k +1
            grad(k) = 4.0_dp * self%work_1(i, a)
          end do
        end do

      case (2)  ! UHF
        if (.not. associated(fock_b)) &
          call show_message('Failed to use Fock arrays in calc_orb_grad', &
                            with_abort)
        ! Alpha
        call unpack_matrix(fock_a, self%work_1)
        call dgemm('N', 'N', nbf, nvir_a, nbf, &
                   1.0_dp, self%work_1, nbf, &
                           mo_a(:, nocc_a+1:), nbf, &
                   0.0_dp, self%work_2, nbf)
        call dgemm('T', 'N', nocc_a, nvir_a, nbf, &
                   1.0_dp, mo_a(:, 1:nocc_a), nbf, &
                           self%work_2, nbf, &
                   0.0_dp, self%work_1(1:nocc_a, 1:nvir_a), nocc_a)
        k = 0
        do i = 1, nocc_a
          do a= 1, nvir_a
            k = k +1
            grad(k) = 2.0_dp * self%work_1(i, a)
          end do
        end do

        ! Beta
        self%work_1 = 0.0_dp
        call unpack_matrix(fock_b, self%work_1)
        call dgemm('N', 'N', nbf, nvir_b, nbf, &
                   1.0_dp, self%work_1, nbf, &
                           mo_b(:, nocc_b+1:), nbf, &
                   0.0_dp, self%work_2, nbf)
        call dgemm('T', 'N', nocc_b, nvir_b, nbf, &
                   1.0_dp, mo_b(:, 1:nocc_b), nbf, &
                           self%work_2, nbf, &
                   0.0_dp, self%work_1(1:nocc_b, 1:nvir_b), nocc_b)
        do i = 1, nocc_b
          do a= 1, nvir_b
            k = k +1
            grad(k) = 2.0_dp * self%work_1(i, a)
          end do
        end do

      case (3)  ! ROHF
        self%work_1 = 0.0_dp
        self%work_3 = 0.0_dp
        call unpack_matrix(fock_a, self%work_1)
        call unpack_matrix(fock_b, self%work_3)
        call dgemm('N', 'N', nbf, nvir_a, nbf, &
                   1.0_dp, self%work_1, nbf, &
                           mo_a(:, nocc_a+1:), nbf, &
                   0.0_dp, self%work_2, nbf)
        call dgemm('T', 'N', nocc_a, nvir_a, nbf, &
                   1.0_dp, mo_a(:, 1:nocc_a), nbf, &
                           self%work_2, nbf, &
                   0.0_dp, self%work_1(1:nocc_a, 1:nvir_a), nocc_a)
        self%work_2 = 0.0_dp
        call dgemm('N', 'N', nbf, nvir_b, nbf, &
                   1.0_dp, self%work_3, nbf, &
                           mo_a(:, nocc_b+1:), nbf, &
                   0.0_dp, self%work_2, nbf)
        call dgemm('T', 'N', nocc_b, nvir_b, nbf, &
                   1.0_dp, mo_a(:, 1:nocc_b), nbf, &
                           self%work_2, nbf, &
                   0.0_dp, self%work_3(1:nocc_b, 1:nvir_b), nocc_b)
        k = 0
        do i = 1, nocc_a
          if (i <= nocc_b ) then
            istart = nbf - nocc_b
          else
            istart = nbf - nocc_a
          end if
          do a= 1, istart
            k = k +1
            if (i <= nocc_b) then
              if (a <= (nocc_a - nocc_b)) then
                grad(k) = 2.0_dp * self%work_3(i,a)
              else
                grad(k) = 2.0_dp * (self%work_3(i,a) + self%work_1(i, a - nocc_a + nocc_b))
              end if
            else
              grad(k) = 2.0_dp * self%work_1(i, a)
            end if
          end do
        end do
      end select
    end associate
  end subroutine calc_orb_grad

  subroutine bfgs_step(self, step)
    implicit none
    class(soscf_converger) :: self
    real(kind=dp), intent(out) :: step(:)

    real(kind=dp) :: displn(self%nvec)
    real(kind=dp) :: dgrad(self%nvec)
    real(kind=dp) :: updti(self%nvec)
    real(kind=dp) :: alpha, beta, t, t1, t2, t3, t4, scale, norm_disp
    real(kind=dp) :: s1, s2, s3, s4, s5, s6
    integer :: i, j

    ! Initialize displacement and preconditioned gradient difference
    if (self%m_history < 1) then
      scale = 1.0_dp
      step = -scale * self%h_inv * self%grad
      where (isnan(step)) step = 0.0_dp
    else
      displn = self%h_inv * self%grad
      dgrad = self%grad-self%grad_prev
      updti = self%h_inv * dgrad
      if (self%m_history > 1) then
        do i = 1, self%m_history-1
          s1 = dot_product(self%s_history(:,i), self%y_history(:,i))
          s2 = dot_product(self%y_history(:,i), self%upd_history(:,i))
          s3 = dot_product(self%s_history(:,i), self%grad)
          s4 = dot_product(self%upd_history(:,i), self%grad)
          s5 = dot_product(self%s_history(:,i), dgrad)
          s6 = dot_product(self%upd_history(:,i), dgrad)

          s1 = 1.0d0 / s1
          s2 = 1.0d0 / s2
          t = 1.0d0 + s1 / s2
          t2 = s1 * s3
          t4 = s1 * s5
          t1 = t * t2 - s1 * s4
          t3 = t * t4 - s1 * s6
          displn = displn + t1 * self%s_history(:,i) - t2 * self%upd_history(:,i)
          updti  = updti  + t3 * self%s_history(:,i) - t4 * self%upd_history(:,i)
        end do
      end if

      ! Final correction using current dgrad and updti
      s1 = dot_product(step, dgrad)
      s2 = dot_product(dgrad, updti)
      s3 = dot_product(step, self%grad)
      s4 = dot_product(updti, self%grad)

      s1 = 1.0_dp / s1
      s2 = 1.0_dp /s2
      t = 1.0_dp + s1 / s2
      t2 = s1 * s3
      t1 = t * t2 - s1 * s4
      displn = displn + t1 * step - t2 * updti
      self%s_history(:, self%m_history) = self%step
      self%y_history(:, self%m_history) = dgrad
      self%upd_history(:, self%m_history) = updti
      step = -displn
    end if
    norm_disp = sqrt(dot_product(step, step)/self%nvec)
    if (norm_disp > 0.1)  then
       step = step*0.1/norm_disp
    end if
  end subroutine bfgs_step

  subroutine rotate_orbs(self, step, nocc_a, nocc_b, mo)
    implicit none

    class(soscf_converger) :: self
    real(kind=dp), intent(in)        :: step(:)
    integer, intent(in)              :: nocc_a, nocc_b
    real(kind=dp), intent(inout)     :: mo(:,:)

    integer            :: nbf, i, idx
    logical :: second_term

    nbf = self%nbf
    self%work_1 = 0
    self%work_2 = 0
    idx = 0
    second_term = .true.
    if (self%scf_type == 3) then! ROHF
      second_term = .false.
    end if
    call exp_scaling(self%work_1, step, idx, nocc_a, nocc_b, nbf, second_term)
    call orthonormalize(self%work_1, nbf)

    call dgemm('N','N', nbf, nbf, nbf, 1.0_dp, mo, nbf, self%work_1, nbf, 0.0_dp, self%work_2, nbf)
    mo = self%work_2

  contains

    subroutine exp_scaling(G, step, idx, nocc_a, nocc_b, nbf, second_term)
      real(kind=dp), intent(out)     :: G(:,:)
      real(kind=dp), intent(in)      :: step(:)
      integer, intent(inout)         :: idx
      integer, intent(in)            :: nocc_a, nocc_b, nbf

      real(kind=dp), allocatable     :: K(:,:), K2(:,:)
      integer                        :: occ, virt, istart, i
      logical, intent(inout)         :: second_term

      allocate(K(nbf,nbf), source=0.0_dp)
      allocate(K2(nbf,nbf), source=0.0_dp)

      do occ = 1, nocc_a
        istart = merge(nocc_b+1, nocc_a+1, occ <= nocc_b)
        do virt = istart, nbf
          idx = idx + 1
          K(virt, occ) =  step(idx)
          K(occ, virt) = -step(idx)
        end do
      end do

      G = 0.0_dp
      do i = 1, nbf
        G(i,i) = 1.0_dp
      end do
      G = G + K
      if (second_term) then
        call dgemm('N','T', nbf, nbf, nbf, 1.0_dp, K, nbf, K, nbf, 0.0_dp, K2, nbf)
        G = G + 0.5_dp * K2
      end if

      deallocate(K, K2)
    end subroutine exp_scaling

    subroutine orthonormalize(G, nbf)
      real(kind=dp), intent(inout) :: G(:,:)
      integer, intent(in)          :: nbf
      integer                      :: i, j
      real(kind=dp)                :: norm, dot

      do i = 1, nbf
        norm = sqrt(dot_product(G(:,i), G(:,i)))
        call dscal(nbf, 1.0_dp / norm, G(1,i), 1)
        if (i == nbf) cycle
        do j = i+1, nbf
          dot = dot_product(G(:,i), G(:,j))
          call daxpy(nbf, -dot, G(1,i), 1, G(1,j), 1)
        end do
      end do
    end subroutine orthonormalize

  end subroutine rotate_orbs

  !> @brief Computes MO energies from the Fock matrix and updated MO coefficients.
  !> @details Projects the Fock matrix onto the MO basis: F_mo = C^T * F * C,
  !>          taking diagonal elements as approximate MO energies for monitoring.
  !> @param[in] self The SOSCF converger object (for nbf).
  !> @param[in] fock Packed Fock matrix (nbf_tri).
  !> @param[in] mo_coeffs MO coefficients (nbf, nbf).
  !> @param[out] mo_energies MO energies (nbf).
  subroutine compute_mo_energies(self, fock, mo_coeffs, mo_energies, work_1, work_2)
    use mathlib, only: unpack_matrix
    implicit none
    class(soscf_converger), intent(in) :: self
    real(kind=dp), intent(in) :: fock(:)
    real(kind=dp), intent(in) :: mo_coeffs(:, :)
    real(kind=dp), intent(out) :: mo_energies(:)
    real(kind=dp), intent(inout) :: work_1(:, :)
    real(kind=dp), intent(inout) :: work_2(:, :)

    real(kind=dp), allocatable :: f_full(:, :), temp(:, :), f_mo(:, :)
    integer :: nbf, i

    nbf = self%nbf

    call unpack_matrix(fock, work_1)
    call dgemm('T', 'N', nbf, nbf, nbf, &
               1.0_dp, mo_coeffs, nbf, &
                       work_1, nbf, &
               0.0_dp, work_2, nbf)
    call dgemm('N', 'N', nbf, nbf, nbf, &
               1.0_dp, work_2, nbf, &
                       mo_coeffs, nbf, &
               0.0_dp, work_1, nbf)

    do i = 1, nbf
      mo_energies(i) = work_1(i, i)
    end do

  end subroutine compute_mo_energies

  !==============================================================================
  ! Line Search Subroutine for SOSCF
  !==============================================================================
  !> @brief Finds the optimal step size for orbital updates in the SOSCF method
  !>        to minimize the system energy along the LBFGS search direction.
  !>
  !> @detail This subroutine implements an adaptive line search algorithm that
  !>         combines backtracking with polynomial interpolation (quadratic or cubic)
  !>         to efficiently find a step size satisfying the Armijo and Wolfe conditions.
  !>
  !> The algorithm operates as follows:
  !> 1. Initialization:
  !>    - Computes the initial energy `e0` and directional derivative `dphi0` along the search direction.
  !>    - Sets an initial step size `alpha_try` with bounds to prevent excessively small or large values.
  !> 2. Search Loop:
  !>    - Applies an orbital rotation using the current step size `alpha_try * step`.
  !>    - Updates density matrices and calculates the new energy `e_new`.
  !>    - Evaluates the Armijo condition (sufficient decrease) and Wolfe condition (curvature).
  !>    - If the Armijo condition is satisfied, the search terminates.
  !>    - Otherwise, adjusts the step size using polynomial interpolation (quadratic on the first iteration,
  !>      cubic on subsequent iterations) or simple backtracking.
  !> 3. Additional Checks:
  !>    - Tracks the best step size (`alpha_best`) corresponding to the lowest energy encountered.
  !>    - Terminates if the maximum number of iterations or a minimum step size threshold is reached.
  !> 4. Exit Conditions:
  !>    - Satisfaction of the Armijo condition.
  !>    - Negligible energy change between iterations (optional, enforced after a minimum number of iterations).
  !>    - Reaching the maximum number of iterations or a minimum step size.
  !>
  !> @section Supported Interpolation Options:
  !>    - Backtracking: Simple step size reduction (interp_type = 0).
  !>    - Quadratic: Quadratic polynomial interpolation (interp_type = 1).
  !>    - Cubic: Cubic polynomial interpolation (default, interp_type = 2).
  !>
  !> @param[inout] self        SOSCF converger instance containing data and parameters.
  !> @param[inout] alpha       Initial step size guess (in), optimal step size (out).
  !> @param[in]    grad        Orbital gradient in the orthogonal basis.
  !> @param[in]    step        LBFGS search direction.
  !> @param[in]    fock_a      Packed alpha Fock matrix.
  !> @param[in]    fock_b      Packed beta Fock matrix (UHF only).
  !> @param[in]    occ_a       Alpha orbital occupations.
  !> @param[in]    occ_b       Beta orbital occupations (UHF only).
  !> @param[inout] mo_a        Alpha molecular orbital coefficients (updated during trials).
  !> @param[inout] mo_b        Beta molecular orbital coefficients (updated during trials).
  !> @param[in]    interp_type Optional: Interpolation method (0=backtracking, 1=quadratic, 2=cubic; default=2).
  !==============================================================================
  subroutine line_search(self, alpha, grad, step, &
                         fock_a, fock_b, occ_a, occ_b, &
                         mo_a, mo_b, pfon, interp_type)
    class(soscf_converger) :: self
    real(kind=dp), intent(inout) :: alpha
    real(kind=dp), intent(in) :: grad(:)
    real(kind=dp), intent(in) :: step(:)
    real(kind=dp), pointer, intent(in) :: fock_a(:)
    real(kind=dp), pointer, intent(in) :: fock_b(:)
    real(kind=dp), intent(inout) :: mo_a(:,:)
    real(kind=dp), intent(inout) :: mo_b(:,:)
    real(kind=dp), pointer, intent(in) :: occ_a(:)
    real(kind=dp), pointer, intent(in) :: occ_b(:)
    real(kind=dp), pointer :: mo_e_a(:)
    real(kind=dp), pointer :: mo_e_b(:)
    type(pfon_t), pointer :: pfon
    integer, optional, intent(in) :: interp_type

    ! Local variables
    real(kind=dp) :: e0, e_new, dphi0, alpha_try, alpha_prev, e_prev
    real(kind=dp) :: alpha_best, e_best ! Best alpha value found
    real(kind=dp) :: dphi_new           ! New directional derivative for Wolfe condition
    real(kind=dp) :: a, b, c, disc      ! For polynomial interpolation

    ! Line search parameters
    real(kind=dp), parameter :: c1 = 1.0e-4_dp ! Armijo condition parameter
    real(kind=dp), parameter :: c2 = 0.9_dp    ! Wolfe condition parameter
    real(kind=dp), parameter :: min_alpha = 1.0e-6_dp ! Minimum step size
    real(kind=dp), parameter :: max_alpha = 2.0_dp    ! Maximum step size
    real(kind=dp), parameter :: threshold = 1.0e-10_dp ! Small number threshold
    real(kind=dp), parameter :: energy_tol = 1.0e-10_dp ! Energy change threshold
    integer, parameter :: max_iter = 25 ! Maximum line search iterations
    integer, parameter :: min_iter_for_energy_check = 16 ! Minimum number of iterations

    integer :: i, j, interp, istat
    logical :: found_better, wolfe_satisfied, armijo_satisfied

    ! Set interpolation method
    interp = 2  ! Default to cubic interpolation
    if (present(interp_type)) interp = min(max(interp_type, 0), 2)

    associate (nbf => self%nbf, &
               nocc_a => self%nocc_a, &
               nocc_b => self%nocc_b, &
               mo_e_a => self%dat%buffer(self%dat%slot)%mo_e_a, &
               mo_e_b => self%dat%buffer(self%dat%slot)%mo_e_b)

      ! Calculate initial energy
      e0 = traceprod_sym_packed(self%dens_a, fock_a, nbf)
      if (self%scf_type > 1) &  ! UHF
        e0 = e0 + traceprod_sym_packed(self%dens_b, fock_b, nbf)

      ! Calculate directional derivative (should be negative for descent)
      dphi0 = dot_product(grad, step)

      ! Initialize search parameters
      alpha_try = min(max(alpha, 0.1_dp), 1.0_dp)  ! Initial step size
      alpha_prev = 0.0_dp  ! Previous step size
      e_prev = e0  ! Previous energy

      ! Initialize best solution tracking
      alpha_best = alpha_try
      e_best = e0
      found_better = .false.

      if (self%verbose>1) then
        write(iw, '(A,E20.10)') 'DEBUG line_search: init: dot_product(grad, step)=', dphi0
        write(iw, '(A,E20.10)') 'DEBUG line_search: init: e0=', e0
        write(iw, '(A,E20.10)') 'DEBUG line_search: init: alpha_try=', alpha_try
      end if

      ! Line search loop
      do i = 1, max_iter

        if (self%verbose>1) then
          write(iw, '(A,I3)') 'DEBUG line_search: loop: first 3x3 MO_a: Iter=', i
          do j = 1, min(3, nbf)
            write(iw, '(3E20.10)') self%mo_a(j, 1:min(3, nbf))
          end do
        end if

        ! Apply orbital rotation with current step size
        if (self%scf_type == 1) then
          call self%rotate_orbs(step * alpha_try, self%nocc_a, self%nocc_a, mo_a)

          if (associated(pfon)) then
            call compute_mo_energies(self, fock_a, mo_a, mo_e_a, self%work_1, self%work_2)
            call pfon%compute_occupations(mo_e_a)
            call pfon%build_density(self%dens_a, mo_a, self%work_1, self%work_2)
          else
            call orb_to_dens(self%dens_a, mo_a, occ_a, nocc_a, nbf, nbf)
          end if
        else
          call self%rotate_orbs(step * alpha_try, self%nocc_a, self%nocc_a, mo_a)
          call self%rotate_orbs(step * alpha_try, self%nocc_b, self%nocc_b, mo_b)
          if (associated(pfon)) then
            call compute_mo_energies(self, fock_a, mo_a, mo_e_a, self%work_1, self%work_2)
            call compute_mo_energies(self, fock_b, mo_b, mo_e_b, self%work_1, self%work_2)
            call pfon%compute_occupations(mo_e_a, mo_e_b)
            call pfon%build_density(self%dens_a, mo_a, self%work_1, self%work_2, self%dens_b, mo_b)
          else
            call orb_to_dens(self%dens_a, mo_a, occ_a, nocc_a, nbf, nbf)
            call orb_to_dens(self%dens_b, mo_b, occ_b, nocc_b, nbf, nbf)
          end if
        end if

        if (self%verbose>1) then
          call unpack_matrix(fock_a, self%work_1)
          call dgemm('T', 'N', nbf, nbf, nbf, 1.0_dp, mo_a, nbf, self%work_1, nbf, 0.0_dp, self%work_2, nbf)
          call dgemm('N', 'N', nbf, nbf, nbf, 1.0_dp, self%work_2, nbf, mo_a, nbf, 0.0_dp, self%work_1, nbf)
          write(iw, '(A,E20.10)') 'DEBUG: Max |F_mo(occ,virt)| after rotation = ', maxval(abs(self%work_1(1:nocc_a, nocc_a+1:)))

          call unpack_matrix(fock_b, self%work_1)
          call dgemm('T', 'N', nbf, nbf, nbf, 1.0_dp, mo_b, nbf, self%work_1, nbf, 0.0_dp, self%work_2, nbf)
          call dgemm('N', 'N', nbf, nbf, nbf, 1.0_dp, self%work_2, nbf, mo_b, nbf, 0.0_dp, self%work_1, nbf)
          write(iw, '(A,E20.10)') 'DEBUG: Max |F_mo(occ,virt)| after rotation = ', maxval(abs(self%work_1(1:nocc_a, nocc_a+1:)))

          ! Debug check of density matrix trace
          call unpack_matrix(self%dens_a, self%work_1)
          call dgemm('N', 'N', nbf, nbf, nbf, &
                     1.0_dp, self%work_1, nbf, &
                             self%overlap, nbf, &
                     0.0_dp, self%work_2, nbf)
          write(iw, '(A,E20.10)') &
            'DEBUG orb_to_dens Alpha: Tr(D*S) after update:', &
            sum([(self%work_2(i,i), i=1,nbf)])

          call unpack_matrix(self%dens_b, self%work_1)
          call dgemm('N', 'N', nbf, nbf, nbf, &
                     1.0_dp, self%work_1, nbf, &
                             self%overlap, nbf, &
                     0.0_dp, self%work_2, nbf)
          write(iw, '(A,E20.10)') &
            'DEBUG orb_to_dens Beta: Tr(D*S) after update:', &
            sum([(self%work_2(i,i), i=1,nbf)])
        end if

        ! Calculate new energy
        e_new = traceprod_sym_packed(self%dens_a, fock_a, nbf)
        if (self%scf_type > 1) &  ! UHF
          e_new = e_new + traceprod_sym_packed(self%dens_b, fock_b, nbf)

        if (self%verbose>1) then
          write(iw, '(A,I3,A,E20.10)') 'DEBUG line_search: loop: iter=', i, ' e_new=', e_new
          write(iw, '(A,I3,A,E20.10)') 'DEBUG line_search: loop: iter=', i, ' condition=', &
            e_new - (e0 + c1 * alpha_try * dphi0)
        end if

        ! Calculate gradient at new point for Wolfe condition
        call self%calc_orb_grad(self%grad, fock_a, fock_b, mo_a, mo_b)
        dphi_new = dot_product(self%grad, step)
        if (self%verbose>1) then
          write(iw, '(A,I3,A,E20.10)') 'DEBUG line_search: loop: iter=', i, &
            ' grad_norm=', sqrt(dot_product(self%grad,self%grad))
          write(iw, '(A,I3,A,E20.10)') 'DEBUG line_search: loop: iter=', i, &
            ' dphi_new=', dphi_new
        end if

        ! Check Armijo condition (sufficient decrease)
        armijo_satisfied = (e_new <= e0 + c1 * alpha_try * dphi0)

        ! Check Wolfe condition (curvature)
        wolfe_satisfied = (abs(dphi_new) <= c2 * abs(dphi0))

        if (self%verbose>1) then
          write(iw, '(A,I3,A,L1,A,L1)') 'DEBUG line_search: loop: iter=', i, &
                                        ' Armijo=', armijo_satisfied, &
                                        ' Wolfe=', wolfe_satisfied
        end if

        ! Track best solution found
        if (e_new < e_best) then
          e_best = e_new
          alpha_best = alpha_try
          found_better = .true.
        end if

        ! Check convergence criteria
        if (armijo_satisfied) then
          alpha = alpha_try
          if (self%verbose>1) then
            write(iw, '(A,E20.10)') 'DEBUG line_search: exit: alpha=', alpha
            write(iw, '(A,I3,A)') 'DEBUG line_search: loop: iter=', i, ' converged - Armijo satisfied'
          end if
          exit
        end if

        ! Additional exit condition based on small energy change
        if (i > min_iter_for_energy_check .and. abs(e_new - e_prev) < energy_tol) then
          if (self%verbose > 1) then
            write(iw, '(A,I3,A,E20.10)') &
              'DEBUG line_search: Energy change below threshold at iter=', i, &
              ', delta E=', abs(e_new - e_prev)
          end if
          alpha = alpha_try  ! Accept the current step
          exit
        end if

        ! Interpolation for next step size
        if (i == 1 .and. interp >= 1) then
          ! Quadratic interpolation: phi(alpha) = a*alpha^2 + dphi0*alpha + e0
          a = (e_new - e0 - dphi0 * alpha_try) / (alpha_try**2)
          if (a > 0.0_dp) then  ! Ensure minimum exists
            alpha_try = -dphi0 / (2.0_dp * a)
            ! Safeguard to avoid too small or large steps
            alpha_try = min(max(alpha_try, 0.2_dp * alpha_try), 0.8_dp * alpha_try)
          else
            alpha_try = alpha_try * 0.5_dp  ! Fallback to backtracking
          end if
        else if (i > 1 .and. interp == 2) then
          ! Cubic interpolation: phi(alpha) = a*alpha^3 + b*alpha^2 + dphi0*alpha + e0
          associate (a1 => alpha_prev, &
                     a2 => alpha_try, &
                     p1 => e_prev, &
                     p2 => e_new)
            if (abs(a2 - a1) > threshold) then
              ! Compute cubic coefficients
              block
                real(kind=dp) :: den, r1, r2
                den = (a1 - a2) * a1 * a2
                if (abs(den) > threshold) then
                  a = ((p2 - e0 - dphi0*a2)*a1**2 - (p1 - e0 - dphi0*a1)*a2**2) / (den**2)
                  b = ((p1 - e0 - dphi0*a1)*a2 - (p2 - e0 - dphi0*a2)*a1) / den
                else
                  a = 0.0_dp
                  b = 0.0_dp
                end if

                ! Solve derivative: 3a*alpha^2 + 2b*alpha + dphi0 = 0
                disc = b**2 - 3.0_dp * a * dphi0
                if (abs(a) > threshold .and. disc > -threshold) then
                  ! Ensure discriminant is treated as non-negative within tolerance
                  disc = max(disc, 0.0_dp)
                  r1 = (-b + sqrt(disc)) / (3.0_dp * a)  ! First root
                  r2 = (-b - sqrt(disc)) / (3.0_dp * a)  ! Second root

                  ! Select the root that positive and within bounds
                  if (r1 > 0.0_dp .and. r1 <= a2) then
                    alpha_try = r1
                  else if (r2 > 0.0_dp .and. r2 <= a2) then
                    alpha_try = r2
                  else
                    alpha_try = a2 * 0.5_dp  ! Fallback if no valid root
                  end if
                else if (abs(b) > threshold) then
                  ! Fall back to quadratic: 2b*alpha + dphi0 = 0
                  alpha_try = -dphi0 / (2.0_dp * b)
                  if (alpha_try <= 0.0_dp .or. alpha_try > a2) alpha_try = a2 * 0.5_dp
                else
                  ! Flat or linear: backtrack
                  alpha_try = a2 * 0.5_dp
                end if
              end block
              alpha_try = min(max(alpha_try, min_alpha), max_alpha)
            else
              alpha_try = a2 * 0.5_dp
            end if
          end associate
        else
          alpha_try = alpha_try * 0.5_dp  ! Simple backtracking
        end if

        ! Enforce min/max step size bounds
        alpha_try = min(max(alpha_try, min_alpha), max_alpha)

        if (self%verbose>1) then
          write(iw, '(A,I3,A,E20.10)') &
            'DEBUG line_search: loop: iter=', i, &
            ' alpha_try=', alpha_try
          write(iw, '(A,I3,A,E20.10)') &
            'DEBUG line_search: loop: iter=', i, &
            ' E_prev=', e_prev
          write(iw, '(A,I3,A,E20.10)') &
            'DEBUG line_search: loop: iter=', i, &
            ' E_new=', e_new
        end if

        ! Check for minimum step size
        if (alpha_try < min_alpha) then
          if (self%verbose>1) then
            write(iw, '(A,I3,A)') 'DEBUG line_search: loop: iter=', i, ' alpha < min_alpha'
          end if
          if (found_better) then
            alpha = alpha_best
            if (self%verbose>1) then
              write(iw, '(A,E20.10,A,E20.10)') 'DEBUG line_search: Using best alpha=', &
                alpha_best, ' with energy=', e_best
            end if
          else
            alpha = alpha_prev
          end if
          exit
        end if

        ! Store current values for next iteration
        alpha_prev = alpha_try
        e_prev = e_new

      end do  ! End of line search loop

      ! Final handling for max iterations exit
      if (i > max_iter) then
        if (self%verbose>1) &
          write(iw, '(A,I3,A)') 'DEBUG line_search: Maximum iterations (', max_iter, ') reached'
        if (found_better) then
          alpha = alpha_best
          if (self%verbose>1) &
            write(iw, '(A,E20.10,A,E20.10)') 'DEBUG line_search: Using best alpha=', &
              alpha_best, ' with energy=', e_best
        else
          alpha = alpha_prev  ! Use last attempt if no better solution found
        end if
      end if

    end associate
  end subroutine line_search

end module scf_converger
