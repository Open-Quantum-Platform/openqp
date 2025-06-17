!===============================================================================
! MODULE: scf_addons
!===============================================================================
!
! DESCRIPTION:
!   The scf_addons module provides specialized functionality to enhance SCF
!   convergence for challenging electronic systems. It implements three major
!   techniques:
!       - pseudo-Fractional Occupation Numbers (pFON),
!       - Maximum Overlap Method (MOM), and
!       - level shifting.
!   These methods help with systems exhibiting near-degeneracies, state flipping,
!   or convergence difficulties.
!
! MEMBERS:
!   - pfon_t [TYPE]: Encapsulates pFON functionality for managing fractional
!                    occupations based on temperature-dependent Fermi-Dirac
!                    distributions.
!
! DEPENDENCIES:
!   - precision: Provides `dp` for double precision real numbers.
!   - io_constants: Provides `iw` for output unit.
!   - mathlib: For matrix operations including pack_matrix, unpack_matrix.
!   - messages: For error handling.
!
! PUBLIC INTERFACES:
!   - pfon_t: Type for managing pseudo-Fractional Occupation Numbers.
!   - apply_mom: Implements Maximum Overlap Method for orbital tracking.
!   - level_shift_fock: Applies level shifting to the Fock matrix.
!
! NOTES:
!   - The module's functionality is designed to be used within SCF iterations.
!   - Methods work with RHF, UHF, and ROHF wavefunctions.
!   - The pFON implementation follows the approach described in:
!     https://doi.org/10.1063/1.478177
!
! HISTORY:
!   - [2025] Initial Module Creation - Konstantin Komarov
!     Established this module by extracting and refactoring auxiliary SCF
!     functionality from the main `scf`` module for better code organization
!     and maintainability. Implemented the `pfon_t`` type as a proper object
!     to encapsulate the pFON functionality.
!   - [January 2025] pFON Implementation - Alireza Lashkaripour
!     Developed the pseudo-Fractional Occupation Number (pFON) method
!     functionality that was later integrated into this module.
!   - [2023-2025] Advanced Convergence Methods - Konstantin Komarov
!     Implemented the Maximum Overlap Method (MOM) and level shifting
!     techniques to improve convergence for challenging electronic systems.
!
!===============================================================================

!===============================================================================
! TYPE: pfon_t - PSEUDO-FRACTIONAL OCCUPATION NUMBERS
!===============================================================================
!
! DESCRIPTION:
!   The `pfon_t` type encapsulates functionality for managing fractional
!   occupation numbers in SCF calculations using a temperature-dependent
!   Fermi-Dirac distribution. This technique smooths convergence for systems
!   with near-degeneracies by allowing partial orbital occupations.
!
! MEMBERS:
!   active         [LOGICAL]: Whether pFON is currently enabled.
!   temp           [REAL(dp)]: Current temperature for Fermi-Dirac distribution.
!   beta           [REAL(dp)]: Inverse temperature parameter (1/(kB*T)).
!   last_cooled_temp [REAL(dp)]: Last temperature at which cooling occurred.
!   cooling_rate   [REAL(dp)]: Rate of temperature decrease per iteration.
!   nsmear         [INTEGER]: Number of orbitals to smear around the Fermi level.
!   occ_a          [REAL(dp), POINTER]: Alpha orbital occupations array.
!   occ_b          [REAL(dp), POINTER]: Beta orbital occupations array.
!   scf_type       [INTEGER]: SCF calculation type (1=RHF, 2=UHF, 3=ROHF).
!   nelec          [INTEGER]: Total number of electrons.
!   nelec_a        [INTEGER]: Number of alpha electrons.
!   nelec_b        [INTEGER]: Number of beta electrons.
!   nbf            [INTEGER]: Number of basis functions.
!
! METHODS:
!   init                 - Initializes pFON parameters based on control settings.
!   adjust_temperature   - Dynamically adjusts temperature during iterations.
!   compute_occupations  - Wrapper for helper `pfon_occupations` function.
!                          Calculates fractional occupations from orbital energies.
!   build_density        - Wrapper for helper `build_pfon_density` function.
!                          Constructs density matrices using fractional occupations.
!
! HELPER FUNCTIONS:
!   pfon_occupations     - Standalone function that computes fractional occupations
!   build_pfon_density   - Constructs density matrices from MO coefficients and
!                          fractional occupations.
!
! ALGORITHM:
!   1. Start with high temperature (typically 2000K) to allow significant
!      fractional occupation and smooth energy surface
!   2. Gradually decrease temperature during iterations (cooling_rate parameter)
!   3. Compute Fermi level as average of HOMO and LUMO energies
!   4. Calculate occupations using Fermi-Dirac distribution:
!      n_i = 2/(1+exp((ε_i-εF)/kT)) for RHF
!      n_i = 1/(1+exp((ε_i-εF)/kT)) for UHF/ROHF
!   5. Normalize occupations to preserve total electron count
!   6. Use occupations to build weighted density matrices
!   7. Set temperature to 1K for final iteration to obtain integer occupations
!
! USAGE NOTES:
!   - Temperature gradually decreases during SCF iterations to facilitate convergence
!   - Final iteration typically uses T=0K to obtain integer occupations
!   - Works with all SCF types (RHF, UHF, ROHF) with appropriate occupation patterns
!   - Particularly effective for systems with small HOMO-LUMO gaps or
!     near-degenerate orbital energies
!
!===============================================================================

!===============================================================================
! TYPE: pfon_t - PSEUDO-FRACTIONAL OCCUPATION NUMBERS
!===============================================================================
!
! DESCRIPTION:
!   The `pfon_t` type encapsulates functionality for managing fractional
!   occupation numbers in SCF calculations using a temperature-dependent
!   Fermi-Dirac distribution. This technique smooths convergence for systems
!   with near-degeneracies by allowing partial orbital occupations.
!
! MEMBERS:
!   active         [LOGICAL]: Whether pFON is currently enabled.
!   temp           [REAL(dp)]: Current temperature for Fermi-Dirac distribution.
!   beta           [REAL(dp)]: Inverse temperature parameter (1/(kB*T)).
!   last_cooled_temp [REAL(dp)]: Last temperature at which cooling occurred.
!   cooling_rate   [REAL(dp)]: Rate of temperature decrease per iteration.
!   nsmear         [INTEGER]: Number of orbitals to smear around the Fermi level.
!   occ_a          [REAL(dp), POINTER]: Alpha orbital occupations array.
!   occ_b          [REAL(dp), POINTER]: Beta orbital occupations array.
!   scf_type       [INTEGER]: SCF calculation type (1=RHF, 2=UHF, 3=ROHF).
!   nelec          [INTEGER]: Total number of electrons.
!   nelec_a        [INTEGER]: Number of alpha electrons.
!   nelec_b        [INTEGER]: Number of beta electrons.
!   nbf            [INTEGER]: Number of basis functions.
!
! METHODS:
!   init                 - Initializes pFON parameters based on control settings.
!   adjust_temperature   - Dynamically adjusts temperature during iterations.
!   compute_occupations  - Calculates fractional occupations from orbital energies.
!   build_density        - Constructs density matrices using fractional occupations.
!
! NOTES:
!   - Works with RHF, UHF, ROHF with appropriate occupation patterns.
!   - Works with Second-Order SCF convergence method.
!
!===============================================================================

!===============================================================================
! SUBROUTINE: apply_mom - MAXIMUM OVERLAP METHOD
!===============================================================================
!
! DESCRIPTION:
!   Implements the Maximum Overlap Method (MOM) to maintain consistent orbital
!   ordering between SCF iterations. This helps prevent oscillations and state
!   flipping during convergence, especially for open-shell systems or cases
!   with near-degeneracies.
!
! PARAMETERS:
!   infos         [TYPE(information)]: System information.
!   v_prev        [REAL(dp)]: Previous iteration's MO coefficients.
!   e_prev        [REAL(dp)]: Previous iteration's orbital energies.
!   v_curr        [REAL(dp)]: Current iteration's MO coefficients (reordered on output).
!   e_curr        [REAL(dp)]: Current iteration's orbital energies (reordered on output).
!   s_ao          [REAL(dp)]: Overlap matrix in AO basis.
!   n_occ         [INTEGER]: Number of occupied orbitals.
!   spin_label    [CHARACTER(*)]: Identifier for spin channel ("Alpha" or "Beta").
!   work          [REAL(dp)]: Work array for intermediate calculations.
!   s_mo          [REAL(dp)]: Work array for MO overlap matrix.
!
! ALGORITHM:
!   1. Computes overlap between previous and current MOs: S_MO = V_prev^T * S * V_curr
!   2. For each orbital (occupied+1), finds maximum overlap match
!   3. Reorders current orbitals to maximize consistency with previous iteration
!   4. Ensures proper tracking of HOMO/LUMO and other important orbitals
!
! HELPER ROUTINES:
!   reorder_orbitals - Internal subroutine that performs the actual orbital swapping
!
!===============================================================================

!===============================================================================
! SUBROUTINE: level_shift_fock - VIRTUAL ORBITAL SHIFTING
!===============================================================================
!
! DESCRIPTION:
!   Applies level shifting to the virtual orbitals in the Fock matrix to increase
!   the HOMO-LUMO gap and improve SCF convergence. This technique is particularly
!   useful for systems with small HOMO-LUMO gaps or near-degeneracies.
!
! PARAMETERS:
!   fock_ao       [REAL(dp)]: Fock matrix in AO basis (triangular format).
!   mo_coefs      [REAL(dp)]: MO coefficients.
!   smat_full     [REAL(dp)]: Full overlap matrix.
!   nocc          [INTEGER]: Number of occupied orbitals.
!   nbf           [INTEGER]: Number of basis functions.
!   vshift        [REAL(dp)]: Level shift parameter value.
!   work1, work2  [REAL(dp)]: Work arrays for intermediate calculations.
!
! ALGORITHM:
!   1. Transforms Fock from AO to MO basis: F_MO = C^T * F_AO * C
!   2. Adds shift to diagonal elements corresponding to virtual orbitals
!   3. Transforms modified Fock back to AO basis for use in SCF
!
! USAGE NOTES:
!   - Typically applied in early iterations and gradually reduced
!   - Often combined with DIIS for optimal convergence
!
! NOTES:
!   - Works with RHF and UHF calculations. The ROHF case is handled through
!   the `form_rohf_fock` function in the `scf` module.
!
!===============================================================================
module scf_addons
  use precision, only: dp

  character(len=*), parameter :: module_name = "scf_addons"

  private

  public :: pfon_t
  public :: apply_mom
  public :: level_shift_fock

  !> @brief Type to encapsulate pFON (pseudo-Fractional Occupation Number) functionality
  !> @detail Provides methods for managing fractional occupation numbers in SCF calculations,
  !>         including temperature control, occupation computation, and density building.
  type :: pfon_t
    private
    logical :: active = .false.                   ! Whether pFON is enabled
    real(kind=dp), public :: temp                 ! Current temperature
    real(kind=dp), public :: beta                 ! Inverse temperature (1/(kB * temp))
    real(kind=dp) :: last_cooled_temp = 0.0_dp    ! Last temperature at which cooling occurred
    real(kind=dp) :: cooling_rate = 50.0_dp       ! Temperature cooling rate
    integer :: nsmear = 0                         ! Number of smearing steps
    real(kind=dp), pointer, public :: occ_a(:) => null()  ! Alpha occupations
    real(kind=dp), pointer, public :: occ_b(:) => null()  ! Beta occupations
    integer :: scf_type = 1  ! SCF calculation type (1=RHF, 2=UHF, 3=ROHF)
    integer :: nelec = 0     ! Total number of electrons
    integer :: nelec_a = 0   ! Number of alpha electrons
    integer :: nelec_b = 0   ! Number of beta electrons
    integer :: nbf = 0       ! Number of basis functions
  contains
    procedure :: init => pfon_init
    procedure :: adjust_temperature => pfon_adjust_temperature
    procedure :: compute_occupations => pfon_compute_occupations
    procedure :: build_density => pfon_build_density
  end type pfon_t

contains

  !> @brief Applies the Maximum Overlap Method (MOM) to reorder orbitals.
  !> @detail Reorders the current iteration’s orbitals to maximize overlap
  !>         with the previous iteration’s orbitals,
  !>         ensuring consistent electronic state tracking during SCF convergence
  !>         (useful for avoiding state flipping).
  !> @param[in] infos System information.
  !> @param[in] v_prev Previous iteration’s MO coefficients.
  !> @param[in] e_prev Previous iteration’s orbital energies.
  !> @param[inout] v_curr Current iteration’s MO coefficients (reordered on output).
  !> @param[inout] e_curr Current iteration’s orbital energies (reordered on output).
  !> @param[in] s_ao Overlap matrix in AO basis.
  !> @param[in] n_occ Number of occupied orbitals.
  !> @param[in] spin_label Identifier for spin channel ("Alpha" or "Beta").
  !> @param[inout] work Work array for intermediate calculations (nbf x nbf).
  !> @param[inout] s_mo Work array for MO overlap matrix (nbf x nbf).
  subroutine apply_mom(infos, v_prev, e_prev, v_curr, e_curr, s_ao, n_occ, &
                       spin_label, work, s_mo)
    use precision, only: dp
    use io_constants, only: iw
    use types, only: information

    implicit none

    ! Input/output parameters
    type(information), intent(in) :: infos
    real(kind=dp), intent(in),    dimension(:,:) :: v_prev
    real(kind=dp), intent(in),    dimension(:)   :: e_prev
    real(kind=dp), intent(inout), dimension(:,:) :: v_curr
    real(kind=dp), intent(inout), dimension(:)   :: e_curr
    real(kind=dp), intent(in),    dimension(:,:) :: s_ao
    integer,       intent(in)                    :: n_occ
    character(*),  intent(in)                    :: spin_label
    real(kind=dp), intent(inout), dimension(:,:) :: work
    real(kind=dp), intent(inout), dimension(:,:) :: s_mo

    ! Local variables
    integer :: i, j, k, ip1, nbf
    integer :: max_idx
    real(kind=dp) :: max_overlap, overlap
    logical, allocatable :: reordered(:)

    nbf = size(v_curr, 1)

    if (infos%control%verbose>=1) then
      if (infos%control%rstctmo) then
        write(IW, fmt='(/,"Applying Reodering for ",A," spin channel")') trim(spin_label)
      else
        write(IW, fmt='(/,"Applying MOM for ",A," spin channel")') trim(spin_label)
      end if
    end if

    ! Allocate reordered flag array
    allocate(reordered(nbf), source=.false.)

    ! Calculate overlap between previous and current MOs: s_mo = v_prev^T * s_ao * v_curr
    call dgemm('t', 'n', nbf, nbf, nbf, 1.0_dp, v_prev, nbf, s_ao, nbf, 0.0_dp, work, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, 1.0_dp, work, nbf, v_curr, nbf, 0.0_dp, s_mo, nbf)

    ! Normalize columns to ensure proper comparison
    do i = 1, nbf
      s_mo(:,i) = s_mo(:,i) / max(norm2(s_mo(:,i)), 1.0e-10_dp)
    end do

    ! First, identify the best match for each orbital from the previous iteration
    ! Focus particularly on occupied orbitals and the HOMO-LUMO region
    ! Print information about important orbitals (HOMO, LUMO)
    if (infos%control%verbose>1) then
      write(IW,fmt='(1X,"MOM reordering for ",A," orbitals:")') trim(spin_label)
      write(IW,fmt='(1X,"Old Index → New Index   | Overlap |  Status")')
      write(IW,fmt='(1X,"--------------------------------------------")')
    end if

    ! First pass: check which orbitals need reordering
    do i = 1, nbf
      max_overlap = 0.0_dp
      max_idx = i  ! Default to no change

      ! Find the orbital with maximum overlap
      do j = 1, nbf
        if (.not. reordered(j)) then
          overlap = abs(s_mo(i,j))
          if (overlap > max_overlap) then
            max_overlap = overlap
            max_idx = j
          end if
        end if
      end do

      ! Mark the orbital as reordered and print info for occupied orbitals
      reordered(max_idx) = .true.
      if (infos%control%verbose>1) then
        ! Print info for important orbitals or those being reordered
        if (((i <= n_occ+1) .or. (i /= max_idx)).and. infos%control%verbose>=1) then
          write(IW, fmt='(3X,I3,5X,"→",5X,I3,5X,"| ",F7.5," |")', advance='no') &
            i, max_idx, max_overlap

          ! Add label for HOMO/LUMO
          if (i == n_occ)   write(IW, fmt='(1X,"HOMO")', advance='no')
          if (i == n_occ+1) write(IW, fmt='(1X,"LUMO")', advance='no')

          ! Add status message
          if (i /= max_idx .and. max_overlap < 0.9_dp) then
            write(IW, fmt='(1X,"Reordered (warning: low overlap)")')
          else if (i /= max_idx) then
            write(IW, fmt='(1X,"Reordered")')
          else if (max_overlap < 0.9_dp) then
            write(IW, fmt='(1X,"Unchanged (warning: low overlap)")')
          else
            write(IW, fmt='(1X,"Unchanged")')
          end if
        end if
      end if
    end do

    ! Check if all orbitals were successfully assigned
    if (.not. all(reordered)) then
      write(IW, fmt='(/,"WARNING: Some orbitals could not be properly reordered!")')
      write(IW, fmt='("This may indicate a significant change in electronic structure.")')
    end if

    ! Apply the reordering
    call reorder_orbitals(v_curr, e_curr, s_mo, nbf, &
                          start_mo=1, &
                          end_mo=n_occ+1)

    deallocate(reordered)
  end subroutine apply_mom

  !> @brief Reorders orbitals based on overlap with the previous iteration.
  !> @detail Internal helper routine for 'apply_mom' that swaps orbital
  !>         coefficients and energies to maximize overlap,
  !>         focusing on a specified range of molecular orbitals.
  !> @param[inout] v MO coefficients (reordered on output).
  !> @param[inout] e Orbital energies (reordered on output).
  !> @param[in] smo Overlap matrix between previous and current MOs.
  !> @param[in] nbf Number of basis functions.
  !> @param[in] start_mo First MO to reorder.
  !> @param[in] end_mo Last MO to reorder.
  subroutine reorder_orbitals(v, e, smo, nbf, start_mo, end_mo)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(inout) :: v(nbf,*)
    real(kind=dp), intent(inout) :: e(*)
    real(kind=dp), intent(in) :: smo(nbf,*)
    integer, intent(in) :: nbf, start_mo, end_mo

    integer :: i, j, k, ip1
    integer, allocatable :: reorder_idx(:)
    real(kind=dp) :: smax, tmp_e

    ! Allocate array for reordering indices
    allocate(reorder_idx(nbf), source=0)

    ! Determine the reordering indices based on maximum overlap
    do i = 1, nbf
      smax = 0.0_dp
      reorder_idx(i) = 0

      ! Find maximum overlap
      do j = 1, nbf
        ! Skip already assigned orbitals
        if (any(reorder_idx(1:i-1) == j)) cycle

        if (abs(smo(i,j)) > smax) then
          smax = abs(smo(i,j))
          reorder_idx(i) = j
        end if
      end do

      ! Ensure sign consistency
      if (smo(i, reorder_idx(i)) < 0.0_dp) then
        v(:, reorder_idx(i)) = -v(:, reorder_idx(i))
      end if
    end do

    ! Apply reordering for the specified range
    do i = start_mo, end_mo
      j = reorder_idx(i)

      ! Swap orbital coefficients
      call dswap(nbf, v(1,i), 1, v(1,j), 1)

      ! Swap orbital energies
      tmp_e = e(i)
      e(i) = e(j)
      e(j) = tmp_e

      ! Update reordering indices for remaining swaps
      ip1 = i + 1
      do k = ip1, end_mo
        if (reorder_idx(k) == i) reorder_idx(k) = j
      end do
    end do

    deallocate(reorder_idx)
  end subroutine reorder_orbitals

  !> @brief Computes fractional occupation numbers using
  !>        the pseudo-Fractional Occupation Number (pFON) method.
  !> @detail Implements the pFON method to assign fractional occupations
  !>         via a Fermi-Dirac distribution, smoothing near-degenerate states.
  !>         Reference: https://doi.org/10.1063/1.478177
  !> @author Alireza Lashkaripour, January 2025
  !> @param[in] mo_energy Orbital energies.
  !> @param[in] nbf Number of basis functions.
  !> @param[in] nelec Total number of electrons.
  !> @param[inout] occ Occupation numbers (updated on output).
  !> @param[in] beta_pfon Inverse temperature parameter (1/(kB * T)).
  !> @param[in] scf_type SCF type (1=RHF, 2=UHF, 3=ROHF).
  !> @param[in] nsmear Number of orbitals to smear around the Fermi level.
  !> @param[in] is_beta Flag indicating beta spin calculation (optional).
  !> @param[in] nelec_a Number of alpha electrons (for UHF/ROHF).
  !> @param[in] nelec_b Number of beta electrons (for UHF/ROHF).
  subroutine pfon_occupations(mo_energy, nbf, nelec, occ, beta_pfon, &
                              scf_type, nsmear, is_beta, nelec_a, nelec_b)
    use precision, only: dp
    implicit none

    integer, intent(in) :: nbf
    integer, intent(in) :: nelec, nsmear
    real(kind=dp), intent(in) :: beta_pfon
    real(kind=dp), intent(in) :: mo_energy(nbf)
    real(kind=dp), intent(inout) :: occ(nbf)
    integer, intent(in) :: scf_type ! 1,2,3 RHF,UHF,ROHF
    logical, intent(in), optional :: is_beta
    integer, intent(in), optional :: nelec_a, nelec_b
    real(kind=dp) :: eF, sum_occ
    integer :: i, i_homo, i_lumo, i_low, i_high
    real(kind=dp) :: tmp
    logical :: is_beta_calc
    integer :: n_electrons, n_double, n_single

    is_beta_calc = .false.
    if (present(is_beta)) is_beta_calc = is_beta

    select case (scf_type)
    case(1) ! RHF
      i_homo = max(1, nelec/2)
      n_electrons = nelec

    case(2) ! UHF
      if (.not. present(nelec_a) .or. .not. present(nelec_b)) then
        stop 'UHF requires nelec_a and nelec_b'
      end if
      ! UHF: completely independent alpha and beta
      if (is_beta_calc) then
        i_homo = max(1, nelec_b)
        n_electrons = nelec_b
      else
        i_homo = max(1, nelec_a)
        n_electrons = nelec_a
      end if

    case(3) ! ROHF
      if (.not. present(nelec_a) .or. .not. present(nelec_b)) then
        stop 'ROHF requires nelec_a and nelec_b'
      end if
      ! ROHF: same spatial orbitals, different occupations
      n_double = nelec_b
      n_single = nelec_a - nelec_b
      if (is_beta_calc) then
        i_homo = n_double
        n_electrons = nelec_b
      else
        i_homo = n_double + n_single
        n_electrons = nelec_a
      end if
    end select

    i_lumo = i_homo + 1
    if (i_lumo > nbf) i_lumo = nbf

    ! Calculate Fermi level
    eF = 0.5_dp * (mo_energy(i_homo) + mo_energy(i_lumo))

    if (nsmear <= 0) then
      do i = 1, nbf
        tmp = beta_pfon * (mo_energy(i) - eF)
        if (scf_type == 1) then  ! RHF
          occ(i) = 2.0_dp / (1.0_dp + exp(tmp))
        else  ! UHF or ROHF
          occ(i) = 1.0_dp / (1.0_dp + exp(tmp))
        end if
      end do
    else
      i_low = max(1, i_homo - nsmear)
      i_high = min(nbf, i_lumo + nsmear)

      ! Special handling for ROHF
      if (scf_type == 3) then
        if (is_beta_calc) then
          do i = 1, n_double
            occ(i) = 1.0_dp
          end do
          do i = n_double + 1, nbf
            occ(i) = 0.0_dp
          end do
        else
          do i = 1, n_double
            occ(i) = 1.0_dp
          end do
          do i = n_double + 1, n_double + n_single
            occ(i) = 1.0_dp
          end do
          do i = n_double + n_single + 1, nbf
            occ(i) = 0.0_dp
          end do
        end if

        ! Apply smearing only around the Fermi level
        do i = i_low, i_high
          tmp = beta_pfon * (mo_energy(i) - eF)
          occ(i) = occ(i) / (1.0_dp + exp(tmp))
        end do
      else
        ! RHF/UHF handling
        do i = 1, i_low - 1
          if (scf_type == 1) then
            occ(i) = 2.0_dp
          else
            occ(i) = 1.0_dp
          end if
        end do

        do i = i_high + 1, nbf
          occ(i) = 0.0_dp
        end do

        do i = i_low, i_high
          tmp = beta_pfon * (mo_energy(i) - eF)
          if (scf_type == 1) then
            occ(i) = 2.0_dp / (1.0_dp + exp(tmp))
          else
            occ(i) = 1.0_dp / (1.0_dp + exp(tmp))
          end if
        end do
      end if
    end if

    ! Normalize occupations
    sum_occ = sum(occ(1:nbf))
    if (sum_occ < 1.0e-14_dp) then
      sum_occ = 1.0_dp
    end if
    occ(1:nbf) = occ(1:nbf) * (real(n_electrons,dp) / sum_occ)

  end subroutine pfon_occupations

  !> @brief Builds density matrices using fractional occupation numbers for the pFON method.
  !> @detail Constructs density matrices from molecular orbital coefficients
  !>         and fractional occupations.
  !> @param[inout] pdmat Density matrices (triangular format, updated on output).
  !> @param[in] mo_a Alpha MO coefficients.
  !> @param[in] mo_b Beta MO coefficients (UHF only).
  !> @param[in] occ_a Alpha occupation numbers.
  !> @param[in] occ_b Beta occupation numbers (UHF/ROHF).
  !> @param[in] scf_type SCF type (1=RHF, 2=UHF, 3=ROHF).
  !> @param[in] nbf Number of basis functions.
  !> @param[in] nelec_a Number of alpha electrons.
  !> @param[in] nelec_b Number of beta electrons.
  !> @param[inout] dtmp Work array for density matrix construction.
  !> @param[inout] work Additional work array.
  subroutine build_pfon_density(pdmat_a, mo_a, occ_a, scf_type, nbf, dtmp, work, &
                                pdmat_b, mo_b, occ_b, nelec_a, nelec_b)
    use precision, only: dp
    use mathlib, only: pack_matrix
    implicit none

    real(kind=dp), intent(inout) :: pdmat_a(:)
    real(kind=dp), intent(in) :: mo_a(:,:)
    real(kind=dp), intent(in) :: occ_a(:)
    integer, intent(in) :: nbf, scf_type
    real(kind=dp), intent(inout) :: dtmp(:,:), work(:,:)
    real(kind=dp), intent(inout), optional :: pdmat_b(:)
    real(kind=dp), intent(in), optional :: mo_b(:,:)
    real(kind=dp), intent(in), optional :: occ_b(:)
    integer, intent(in), optional :: nelec_a, nelec_b

    integer :: i, mu, nu
    integer :: n_double, n_single
    real(kind=dp) :: occ_factor


    select case(scf_type)
    case(1)  ! RHF
      ! Scale MO coefficients by square root of occupation numbers
      do i = 1, nbf
        if (occ_a(i) > 1.0e-14_dp) then
            call dger(nbf, nbf, occ_a(i), mo_a(:,i), 1, mo_a(:,i), 1, dtmp, nbf)
        end if
      end do
      pdmat_a = 0.0_dp
      call pack_matrix(dtmp, pdmat_a)

    case(2)  ! UHF
      do i = 1, nbf
        if (occ_a(i) > 1.0e-14_dp) then
          call dger(nbf, nbf, occ_a(i), mo_a(:,i), 1, mo_a(:,i), 1, dtmp, nbf)
        end if
      end do
      pdmat_a = 0.0_dp
      call pack_matrix(dtmp, pdmat_a)

      dtmp(:,:) = 0.0_dp
      do i = 1, nbf
        if (occ_b(i) > 1.0e-14_dp) then
          call dger(nbf, nbf, occ_b(i), mo_b(:,i), 1, mo_b(:,i), 1, dtmp, nbf)
        end if
      end do
      pdmat_b = 0.0_dp
      call pack_matrix(dtmp, pdmat_b)

    case(3)  ! ROHF
      n_double = nelec_b
      n_single = nelec_a - nelec_b

      dtmp(:,:) = 0.0_dp
      do i = 1, nbf
        if (occ_a(i) > 1.0e-14_dp) then
          if (i <= n_double) then
            occ_factor = occ_a(i)
          else if (i <= n_double + n_single) then
            occ_factor = 1.0_dp
          else
            occ_factor = occ_a(i)  ! Virtual orbitals
          end if

          ! dtmp += occ_factor * mo_a(:,i) * mo_a(:,i)^T
!         call dger(nbf, nbf, occ_factor, mo_a(:,i), 1, mo_a(:,i), 1, dtmp, nbf)
          do mu = 1, nbf
            do nu = 1, nbf
              dtmp(mu,nu) = dtmp(mu,nu) + occ_factor * mo_a(mu,i)*mo_a(nu,i)
            end do
          end do
        end if
      end do
      pdmat_a = 0.0_dp
      call pack_matrix(dtmp, pdmat_a)

      dtmp(:,:) = 0.0_dp
      do i = 1, nbf
        if (occ_b(i) > 1.0e-14_dp) then
          if (i <= n_double) then
            occ_factor = occ_b(i)
          else
            occ_factor = 0.0_dp
          end if

          ! dtmp += occ_factor * mo_a(:,i) * mo_a(:,i)^T
!         call dger(nbf, nbf, occ_factor, mo_a(:,i), 1, mo_a(:,i), 1, dtmp, nbf)
          do mu = 1, nbf
            do nu = 1, nbf
              dtmp(mu,nu) = dtmp(mu,nu) + occ_factor * mo_a(mu,i)*mo_a(nu,i)
            end do
          end do
        end if
      end do
      pdmat_b = 0.0_dp
      call pack_matrix(dtmp, pdmat_b)
    end select

  end subroutine build_pfon_density

  !> @brief Applies level shifting to the Fock matrix for improved SCF convergence.
  !> @detail Modifies the diagonal elements of the Fock matrix in the MO basis
  !>         for virtual orbitals by adding a shift parameter,
  !>         then transforms the result back to the AO basis.
  !> @param[inout] fock_ao Fock matrix in AO basis (triangular format, updated on output).
  !> @param[in] mo_coefs MO coefficients.
  !> @param[in] smat_full Full overlap matrix.
  !> @param[in] nocc Number of occupied orbitals.
  !> @param[in] nbf Number of basis functions.
  !> @param[in] vshift Level shift parameter value.
  subroutine level_shift_fock(fock_ao, mo_coefs, smat_full, nocc, nbf, vshift, &
                              work1, work2)
    use precision, only: dp
    use mathlib, only: orthogonal_transform_sym, &
                       orthogonal_transform2, &
                       unpack_matrix, &
                       pack_matrix

    implicit none

    integer, intent(in) :: nocc, nbf
    real(kind=dp), intent(inout) :: fock_ao(:)
    real(kind=dp), intent(in) :: mo_coefs(:,:)
    real(kind=dp), intent(in) :: smat_full(:,:)
    real(kind=dp), intent(in) :: vshift
    real(kind=dp), intent(inout) :: work1(:,:)
    real(kind=dp), intent(inout) :: work2(:,:)

    ! Local variables
    real(kind=dp), allocatable :: fock_mo_full(:,:), fock_mo(:), work_matrix(:,:)
    integer :: i, nbf_tri

    nbf_tri = nbf*(nbf+1)/2

    work1 = 0.0_dp
    work2 = 0.0_dp

    ! Allocate work arrays
    allocate(fock_mo_full(nbf, nbf), &
             fock_mo(nbf_tri), &
             work_matrix(nbf, nbf), &
             source=0.0_dp)

    ! Transform Fock from AO to MO basis: F_MO = C^T * F_AO * C
    call orthogonal_transform_sym(nbf, nbf, fock_ao, mo_coefs, nbf, fock_mo)

    ! Unpack triangular matrices to full format
    call unpack_matrix(fock_mo, fock_mo_full)

    ! Apply level shift to virtual orbitals in F_MO
    do i = nocc + 1, nbf
      fock_mo_full(i, i) = fock_mo_full(i, i) + vshift
    end do

    ! Back-transform ROHF Fock matrix to AO basis
    call dsymm('l', 'u', nbf, nbf, &
               1.0_dp, smat_full, nbf, &
                       mo_coefs, nbf, &
               0.0_dp, work1, nbf)
    call orthogonal_transform2('t', nbf, nbf, work1, nbf, fock_mo_full, nbf, &
                               work_matrix, nbf, work2)

    ! Pack the result back to triangular form
    call pack_matrix(work_matrix, fock_ao)

    deallocate(fock_mo_full, fock_mo, work_matrix)
  end subroutine level_shift_fock

  !> @brief Initialize pFON parameters
  !> @detail Sets up temperature, inverse temperature (beta),
  !>         and other pFON parameters based on input controls.
  !> @param[in] control Control structure containing pFON settings
  !> @param[in] nbf Number of basis functions
  !> @param[in] nelec Total number of electrons
  !> @param[in] nelec_a Number of alpha electrons
  !> @param[in] nelec_b Number of beta electrons
  !> @param[in] scf_type SCF type (1=RHF, 2=UHF, 3=ROHF)
  !> @param[inout] occ_a Pointer to alpha occupations array
  !> @param[inout] occ_b Pointer to beta occupations array (only for UHF/ROHF)
  subroutine pfon_init(this, control, nbf, nelec, nelec_a, nelec_b, scf_type, occ_a, occ_b)
    use types, only: control_parameters
    use constants, only: kB_HaK

    implicit none

    class(pfon_t), intent(inout) :: this
    type(control_parameters), intent(in) :: control
    integer, intent(in) :: nbf, nelec, nelec_a, nelec_b, scf_type
    real(dp), target, intent(inout) :: occ_a(:)
    real(dp), target, optional, intent(inout) :: occ_b(:)

    this%active = control%pfon
    if (.not. this%active) return

    this%nbf = nbf
    this%nelec = nelec
    this%nelec_a = nelec_a
    this%nelec_b = nelec_b
    this%scf_type = scf_type

    ! Set temperature parameters
    this%temp = control%pfon_start_temp
    if (this%temp <= 0.0_dp) this%temp = 2000.0_dp  ! Default temperature
    this%beta = 1.0_dp / (kB_HaK * this%temp)
    this%cooling_rate = control%pfon_cooling_rate
    if (this%cooling_rate <= 0.0_dp) this%cooling_rate = 50.0_dp

    ! Set number of orbitals to smear
    this%nsmear = int(control%pfon_nsmear)

    ! Set pointers to occupation arrays
    this%occ_a => occ_a
    if (present(occ_b)) this%occ_b => occ_b

  end subroutine pfon_init

  !> @brief Adjust pFON temperature based on convergence status
  !> @detail Dynamically modifies the temperature and beta parameters during SCF
  !>         iterations, reducing temperature as convergence improves.
  !> @param[in] iter Current SCF iteration
  !> @param[in] maxit Maximum number of SCF iterations
  !> @param[in] diis_error Current DIIS error
  !> @param[in] conv Convergence threshold
  subroutine pfon_adjust_temperature(this, iter, maxit, diis_error, conv)
    use constants, only: kB_HaK
    class(pfon_t), intent(inout) :: this
    integer, intent(in) :: iter, maxit
    real(dp), intent(in) :: diis_error, conv

    if (.not. this%active) return

    if (iter == maxit) then
      ! Final iteration: set temperature to zero for pure integer occupations
      this%temp = 0.0_dp
    else if (abs(diis_error) < 10.0_dp * conv) then
      ! Near convergence: set to minimum temperature (1K)
      if (this%temp > 1.0_dp) then
        this%last_cooled_temp = this%temp
      end if
      this%temp = 1.0_dp
    else
      ! Not converged yet: continue cooling temperature
      if (this%temp == 1.0_dp .and. this%last_cooled_temp > 1.0_dp) then
        this%temp = this%last_cooled_temp
      end if
      this%temp = this%temp - this%cooling_rate
      if (this%temp < 1.0_dp) then
        this%temp = 1.0_dp
      end if
      this%last_cooled_temp = this%temp
    end if

    ! Calculate beta = 1/(kB*T) for Fermi-Dirac distribution
    if (this%temp > 1.0e-12_dp) then
      this%beta = 1.0_dp / (kB_HaK * this%temp)
    else
      this%beta = 1.0e20_dp  ! Zero temperature
    end if
  end subroutine pfon_adjust_temperature

  !> @brief Compute fractional occupation numbers using pFON method
  !> @detail Uses current orbital energies to calculate occupations
  !>         via a Fermi-Dirac distribution.
  !> @param[in] mo_energy_a Alpha orbital energies
  !> @param[in] mo_energy_b Beta orbital energies (only for UHF)
  subroutine pfon_compute_occupations(this, mo_energy_a, mo_energy_b)
    class(pfon_t), intent(inout) :: this
    real(dp), intent(in) :: mo_energy_a(:)
    real(dp), intent(in), optional :: mo_energy_b(:)

    if (.not. this%active) return

    ! Calculate alpha occupations
    call pfon_occupations(mo_energy_a, this%nbf, this%nelec, this%occ_a, &
                          this%beta, this%scf_type, this%nsmear, &
                          is_beta=.false., nelec_a=this%nelec_a, nelec_b=this%nelec_b)

    ! Calculate beta occupations if needed
    if (this%scf_type > 1 .and. associated(this%occ_b)) then
      if (this%scf_type == 2 .and. present(mo_energy_b)) then
        ! UHF case - use separate beta orbital energies
        call pfon_occupations(mo_energy_b, this%nbf, this%nelec, this%occ_b, &
                            this%beta, this%scf_type, this%nsmear, &
                            is_beta=.true., nelec_a=this%nelec_a, nelec_b=this%nelec_b)
      else
        ! ROHF case - use same orbital energies for alpha and beta
        call pfon_occupations(mo_energy_a, this%nbf, this%nelec, this%occ_b, &
                            this%beta, this%scf_type, this%nsmear, &
                            is_beta=.true., nelec_a=this%nelec_a, nelec_b=this%nelec_b)
      end if
    end if
  end subroutine pfon_compute_occupations

  !> @brief Build density matrices using fractional occupation numbers
  !> @detail Constructs density matrices for the current SCF iteration
  !>         using fractional occupations and MO coefficients.
  !> @param[inout] this pFON type instance.
  !> @param[out] pdmat Density matrices (triangular format).
  !> @param[in] mo_a Alpha MO coefficients.
  !> @param[inout] work1 Work array 1.
  !> @param[inout] work2 Work array 2.
  !> @param[in] mo_b Beta MO coefficients (optional for UHF).
  subroutine pfon_build_density(this, pdmat_a, mo_a, work1, work2, pdmat_b, mo_b)
    class(pfon_t), intent(inout) :: this
    real(kind=dp), intent(out) :: pdmat_a(:)
    real(kind=dp), intent(in) :: mo_a(:,:)
    real(kind=dp), intent(inout) :: work1(:,:)
    real(kind=dp), intent(inout) :: work2(:,:)
    real(kind=dp), intent(out), optional :: pdmat_b(:)
    real(kind=dp), intent(in), optional :: mo_b(:,:)

    if (.not. this%active) return

    ! Nullify work arrays
    work1 = 0.0_dp
    work2 = 0.0_dp

    ! Call the existing build_pfon_density function with appropriate parameters
    select case (this%scf_type)
    case (1) ! RHF
      call build_pfon_density(pdmat_a, mo_a, this%occ_a, this%scf_type, this%nbf, &
                              work1, work2)
    case (2) ! UHF
      call build_pfon_density(pdmat_a, mo_a, this%occ_a, this%scf_type, this%nbf, &
                              work1, work2, pdmat_b, mo_b, this%occ_b)
    case (3) ! ROHF
      call build_pfon_density(pdmat_a, mo_a, this%occ_a, this%scf_type, this%nbf, &
                              work1, work2, pdmat_b, mo_b, this%occ_b, this%nelec_a, this%nelec_b)
    end select
  end subroutine pfon_build_density

end module scf_addons
