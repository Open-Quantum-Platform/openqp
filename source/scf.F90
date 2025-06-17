module scf
  use precision, only: dp

  character(len=*), parameter :: module_name = "scf"

  public :: scf_driver
  public :: fock_jk

contains

  !==============================================================================
  ! Main SCF Driver Subroutine
  !==============================================================================
  !> @brief Performs self-consistent field (SCF) calculations for Hartree-Fock (HF)
  !>        and Density Functional Theory (DFT) methods.
  !> @detail This subroutine implements SCF iterations for
  !>         Restricted (RHF), Unrestricted (UHF), and Restricted Open-Shell (ROHF)
  !>
  !> @section Supported SCF Options for RHF/UHF/ROHF:
  !>          - MOM: Maximum Overlap Method for orbital consistency between SCF iterations.
  !>          - pFON: Pseudo-Fractional Occupation Number for near-degenerate states.
  !>          - Vshift: Level shifting for virtual orbitals.
  !>
  !> @section accelerators Supported SCF Convergence Accelerators:
  !>          - A-DIIS: Augmented Direct Inversion in the Iterative Subspace.
  !>          - E-DIIS: Energy-based DIIS.
  !>          - C-DIIS: Commutator-based DIIS.
  !>          - V-DIIS: Variable DIIS with dynamic switching.
  !>          - SOSCF:  Second-Order SCF convergence method.
  !>
  !> @author Vladimir Mironov - Original author, developed initial SCF module
  !>                            and DIIS drivers (pre-2022).
  !> @author Konstantin Komarov - Added Guest-Saunders ROHF Fock transformation,
  !>                              MOM, Vshift, SOSCF, wrapped pFON into `pfon_t` type,
  !>                              optimized memory usage, added documentation and comments,
  !>                              and cleaned up the code (2023-2025).
  !> @author Mohsen Mazaherifar - Added MPI parallelization support (2023-2024).
  !> @author Alireza Lashkaripour - Implemented pFON functionality (January 2025).
  !>
  !> @date Initial version: pre-2022; Major updates: 2023-2025.
  !>
  !> @param[in]     basis    Basis set information.
  !> @param[inout]  infos    System information and calculation parameters
  !>                         (updated with converged energy and wavefunction).
  !> @param[in]     molGrid  Molecular grid for DFT calculations.
  subroutine scf_driver(basis, infos, molGrid)
    USE precision, only: dp
    use oqp_tagarray_driver
    use constants, only: kB_HaK
    use types, only: information
    use int2_compute, only: int2_compute_t, int2_fock_data_t, &
                            int2_rhf_data_t, int2_urohf_data_t
    use dft, only: dftexcor
    use mod_dft_molgrid, only: dft_grid_t
    use messages, only: show_message, WITH_ABORT
    use guess, only: get_ab_initio_density, get_ab_initio_orbital
    use util, only: measure_time, e_charge_repulsion
    use printing, only: print_mo_range
    use mathlib, only: traceprod_sym_packed, matrix_invsqrt
    use mathlib, only: unpack_matrix
    use io_constants, only: IW
    use basis_tools, only: basis_set
    use scf_converger, only: scf_conv_result, scf_conv, &
                             conv_cdiis, conv_ediis, conv_soscf
    use scf_addons, only: pfon_t, apply_mom, level_shift_fock

    implicit none

    character(len=*), parameter :: subroutine_name = "scf_driver"

    !==============================================================================
    ! Input/Output Arguments
    !==============================================================================
    type(basis_set), intent(in) :: basis              ! Basis set information
    type(information), target, intent(inout) :: infos ! System information & parameters
    type(dft_grid_t), intent(in) :: molGrid           ! Molecular grid for DFT

    !================================================o=o===========================
    ! Matrix Dimensions and Basic Parameters
    !==============================================================================
    integer :: nbf      ! Number of basis functions
    integer :: nbf_tri  ! Size of triangular matrices (nbf*(nbf+1)/2)
    integer :: nbf2     ! Square matrix size (nbf*nbf)
    integer :: nfocks   ! Number of Fock matrices (1 for RHF, 2 for UHF/ROHF)
    integer :: nschwz   ! Number of skipped integrals (integral screening)
    integer :: ok       ! Status flag for memory allocation

    !==============================================================================
    ! SCF Type Parameters
    !==============================================================================
    integer :: scf_type                 ! Type of SCF calculation
    integer, parameter :: scf_rhf  = 1  ! Restricted Hartree-Fock
    integer, parameter :: scf_uhf  = 2  ! Unrestricted Hartree-Fock
    integer, parameter :: scf_rohf = 3  ! Restricted Open-Shell Hartree-Fock
    character(16) :: scf_name = ""      ! Name of the SCF method (RHF/UHF/ROHF)
    logical :: is_dft                   ! True if using DFT, false for HF
    real(kind=dp) :: scalefactor        ! Scaling factor for HF exchange

    !==============================================================================
    ! Electron Counting Parameters
    !==============================================================================
    integer :: nelec    ! Total number of electrons
    integer :: nelec_a  ! Number of alpha electrons
    integer :: nelec_b  ! Number of beta electrons

    !==============================================================================
    ! Iteration Control Parameters
    !==============================================================================
    integer :: i, ii, iter  ! Loop counters and current iteration number
    integer :: maxit        ! Maximum number of SCF iterations

    !==============================================================================
    ! Energy Components
    !==============================================================================
    real(kind=dp) :: ehf      ! Electronic energy (HF part)
    real(kind=dp) :: ehf1     ! One-electron energy
    real(kind=dp) :: nenergy  ! Nuclear repulsion energy
    real(kind=dp) :: etot     ! Total SCF energy
    real(kind=dp) :: e_old    ! Energy from previous iteration
    real(kind=dp) :: psinrm   ! Wavefunction normalization
    real(kind=dp) :: vne      ! Nucleus-electron potential energy
    real(kind=dp) :: vnn      ! Nucleus-nucleus potential energy
    real(kind=dp) :: vee      ! Electron-electron potential energy
    real(kind=dp) :: vtot     ! Total potential energy
    real(kind=dp) :: virial   ! Virial ratio (V/T)
    real(kind=dp) :: tkin     ! Kinetic energy
    real(kind=dp) :: eexc     ! Exchange-correlation energy for DFT
    real(kind=dp) :: totele   ! Total electron density for DFT
    real(kind=dp) :: totkin   ! Total kinetic energy for DFT

    !==============================================================================
    ! DIIS Convergence Acceleration Parameters
    !==============================================================================
    integer :: diis_nfocks       ! Number of Fock matrices for DIIS
    integer :: soscf_nfocks      ! Number of Fock matrices for SOSCF
    integer :: diis_reset        ! Frequency of DIIS reset
    integer :: maxdiis           ! Maximum number of DIIS vectors
    real(kind=dp) :: diis_error  ! DIIS error matrix norm
    character(len=6), dimension(5) :: diis_name          ! Names of DIIS methods
    real(kind=dp), parameter :: ethr_cdiis_big = 2.0_dp  ! DIIS error threshold for C-DIIS
    real(kind=dp), parameter :: ethr_ediis = 1.0_dp      ! DIIS error threshold for E-DIIS
    logical :: diis_reset_condition                      ! Flag for DIIS reset condition

    !==============================================================================
    ! SOSCF Convergence Acceleration Parameters
    !==============================================================================
    logical :: use_soscf            ! Flag to use SOSCF method
    real(kind=dp) :: rms_grad       ! RMS of gradient
    real(kind=dp) :: rms_dp         ! RMS of density different

    !==============================================================================
     ! Virtual Orbital Shift Parameters (for ROHF)
    !==============================================================================
    real(kind=dp) :: vshift        ! Virtual orbital energy shift
    real(kind=dp) :: H_U_gap       ! HOMO-LUMO gap
    real(kind=dp) :: H_U_gap_crit  ! HOMO-LUMO gap critical value
    logical :: vshift_last_iter    ! Flag for last iteration with vshift

    !==============================================================================
    ! MOM (Maximum Overlap Method) Parameters
    !==============================================================================
    logical :: do_mom            ! Flag to enable MOM method
    logical :: initial_mom_iter  ! Flag for first MOM iteration
    logical :: mom_active        ! Flag indicating MOM is currently active
    real(kind=dp), allocatable :: mo_a_prev(:,:)  ! Previous alpha MO coefficients
    real(kind=dp), allocatable :: mo_e_a_prev(:)  ! Previous alpha orbital energies
    real(kind=dp), allocatable :: mo_b_prev(:,:)  ! Previous beta MO coefficients (UHF only)
    real(kind=dp), allocatable :: mo_e_b_prev(:)  ! Previous beta orbital energies (UHF only)

    !==============================================================================
    ! pFON (pseudo-Fractional Occupation Number) Parameters
    !==============================================================================
    logical :: do_pfon        ! Flag to use pFON method
    logical :: do_pfon_final  ! Flag to trigger extra iteration at 1K
    type(pfon_t), pointer :: pfon  ! pFON handler object
    real(kind=dp), allocatable, target :: occ_a(:), occ_b(:) ! Orbital occupations for alpha/beta

    !==============================================================================
    ! Matrices and Vectors for SCF Calculation
    !==============================================================================
    real(kind=dp), allocatable, target :: smat_full(:,:)  ! Full overlap matrix
    real(kind=dp), allocatable, target :: pdmat(:,:)  ! Density matrices in triangular format
    real(kind=dp), allocatable, target :: pfock(:,:)  ! Fock matrices in triangular format
    real(kind=dp), allocatable, target :: rohf_bak(:,:)  ! Backup for ROHF Fock
    real(kind=dp), allocatable, target :: dold(:,:)  ! Old density for incremental builds
    real(kind=dp), allocatable, target :: fold(:,:)  ! Old Fock for incremental builds
    real(kind=dp), allocatable :: pfxc(:,:)  ! DFT exchange-correlation matrix
    real(kind=dp), allocatable :: qmat(:,:)  ! Orthogonalization matrix
    real(kind=dp), allocatable, target :: work1(:,:)  ! Work matrix 1
    real(kind=dp), allocatable, target :: work2(:,:)  ! Work matrix 2
    !==============================================================================
    ! Matrices and Vectors for SCF Calculation
    !==============================================================================
    logical :: do_rstctmo
    real(kind=dp), allocatable, dimension(:,:) :: mo_a_for_rstctmo, mo_b_for_rstctmo
    real(kind=dp), allocatable :: mo_energy_a_for_rstctmo(:)

    !==============================================================================
    ! Tag Arrays for Accessing Data
    !==============================================================================
    real(kind=dp), contiguous, pointer :: smat(:), hcore(:), tmat(:), &
                                          fock_a(:), fock_b(:), &
                                          dmat_a(:), dmat_b(:), &
                                          mo_energy_b(:), mo_energy_a(:), &
                                          mo_a(:,:), mo_b(:,:)
    character(len=*), parameter :: tags_general(3) = &
      (/ character(len=80) :: OQP_SM, OQP_TM, OQP_Hcore /)
    character(len=*), parameter :: tags_alpha(4) = &
      (/ character(len=80) :: OQP_FOCK_A, OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_beta(4) = &
      (/ character(len=80) :: OQP_FOCK_B, OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B /)

    !==============================================================================
    ! SCF Convergence Accelerator Objects
    !==============================================================================
    type(scf_conv) :: conv                           ! SCF convergence driver
    class(scf_conv_result), allocatable :: conv_res  ! SCF convergence result
    integer :: stat                                  ! Status flag for DIIS/SOSCF

    !==============================================================================
    ! Integral Evaluation Objects
    !==============================================================================
    type(int2_compute_t) :: int2_driver                ! Two-electron integral driver
    class(int2_fock_data_t), allocatable :: int2_data  ! Two-electron integral data

    !==============================================================================
    ! Extract Calculation Parameters from Input
    !==============================================================================
    ! Set SCF type (RHF, UHF, or ROHF) and
    ! configure parameters based on SCF type
    select case (infos%control%scftype)
    case (1)
      scf_type = scf_rhf
      scf_name = "RHF"
      nfocks = 1
      diis_nfocks = 1
      soscf_nfocks = 1
    case (2)
      scf_type = scf_uhf
      scf_name = "UHF"
      nfocks = 2
      diis_nfocks = 2
      soscf_nfocks = 2
    case (3)
      scf_type = scf_rohf
      scf_name = "ROHF"
      nfocks = 2
      diis_nfocks = 1
      soscf_nfocks = 2
    end select

    ! Get electron counts
    nelec = infos%mol_prop%nelec
    nelec_a = infos%mol_prop%nelec_a
    nelec_b = infos%mol_prop%nelec_b

    ! Get matrix dimensions
    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2
    nbf2 = nbf*nbf

    ! Get iteration parameters
    maxit = infos%control%maxit

    ! Determine calculation type (HF or DFT)
    is_dft = infos%control%hamilton >= 20

    ! Set HF exchange scaling factor for DFT
    if (is_dft) then
      scalefactor = infos%dft%HFscale
    else
      scalefactor = 1.0_dp
    end if

    !==============================================================================
    ! Retrieve Tag Arrays and Allocate Memory
    !==============================================================================
    ! Get general tag arrays
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_Hcore, hcore)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call tagarray_get_data(infos%dat, OQP_TM, tmat)

    ! Get alpha-spin tag arrays
    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    ! Get beta-spin tag arrays if needed
    if (nfocks > 1) then
      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
      call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    end if

    ! Allocate main work arrays
    ok = 0
    allocate(smat_full(nbf, nbf), &
             pdmat(nbf_tri, nfocks), &
             pfock(nbf_tri, nfocks), &
             rohf_bak(nbf_tri, nfocks), &
             qmat(nbf, nbf), &
             work1(nbf,nbf), &
             work2(nbf,nbf), &
             stat=ok, &
             source=0.0_dp)
    if(ok/=0) call show_message('Cannot allocate memory for SCF', WITH_ABORT)

    ! Allocate incremental Fock building arrays if enabled
    if (infos%control%scf_incremental /= 0) then
      allocate(dold(nbf_tri, nfocks), &
               fold(nbf_tri, nfocks), &
               stat=ok, &
               source=0.0_dp)
      if(ok/=0) call show_message('Cannot allocate memory for SCF: 2',WITH_ABORT)
    end if

    ! Allocate DFT arrays if needed
    if (is_dft) then
      allocate(pfxc(nbf_tri, nfocks), &
               stat=ok, &
               source=0.0_dp)
      if(ok/=0) call show_message('Cannot allocate memory for temporary vectors',WITH_ABORT)
    end if

    !==============================================================================
    ! Initialize pFON Parameters
    !==============================================================================
    do_pfon = .false.
    do_pfon = infos%control%pfon

    if (do_pfon) then
      ! Flag to trigger extra iteration at 1K
      do_pfon_final = .false.

      ! Allocate and initialize occupation arrays
      allocate(occ_a(nbf), source=0.0_dp, stat=ok)

      if (nfocks > 1) then  ! For UHF and ROHF
        allocate(occ_b(nbf), source=0.0_dp, stat=ok)
      end if

      if(ok/=0) call show_message('Cannot allocate memory for occupation arrays',WITH_ABORT)

      ! Set initial occupation numbers based on SCF type
      select case (scf_type)
      case (scf_rhf)
        occ_a(1:nelec/2) = 2.0_dp
      case (scf_uhf)
        occ_a(1:nelec_a) = 1.0_dp
        occ_b(1:nelec_b) = 1.0_dp
      case (scf_rohf)
        occ_a(1:nelec_b) = 2.0_dp          ! Closed shells
        occ_a(nelec_b+1:nelec_a) = 1.0_dp  ! Open shells
        occ_b(1:nelec_b) = 2.0_dp          ! Closed shells only
      end select

      ! Allocate pFON object
      allocate(pfon)

      ! Initialize pFON object
      if (nfocks > 1) then
        call pfon%init(infos%control, nbf, nelec, nelec_a, nelec_b, scf_type, occ_a, occ_b)
      else
        call pfon%init(infos%control, nbf, nelec, nelec_a, nelec_b, scf_type, occ_a)
      end if
    end if

    !==============================================================================
    ! Initialize MOM parameters
    !==============================================================================
    do_mom = infos%control%mom

    if (do_mom) then
      initial_mom_iter = .true.
      mom_active = .false.

      ! Allocate storage for previous iteration's orbitals
      allocate(mo_a_prev(nbf,nbf), source=0.0_dp)
      allocate(mo_e_a_prev(nbf), source=0.0_dp)

      ! For UHF, we need separate storage for beta orbitals
      if (scf_type == scf_uhf) then
        allocate(mo_b_prev(nbf,nbf), source=0.0_dp)
        allocate(mo_e_b_prev(nbf), source=0.0_dp)
      end if
    end if
    !==============================================================================
    ! Initialize XAS parameters
    !==============================================================================
    do_rstctmo = infos%control%rstctmo
    if(do_mom .and. do_rstctmo) call show_message('* Error: Use either MOM or RSTCTMO',WITH_ABORT)
    if (do_rstctmo) then
      allocate(mo_a_for_rstctmo(nbf,nbf), &
                mo_b_for_rstctmo(nbf,nbf), &
                mo_energy_a_for_rstctmo(nbf), &
                source=0.0_dp)
    end if
    !==============================================================================
    ! Initialize Vshift Parameters (currently only works for ROHF)
    !==============================================================================
    vshift = infos%control%vshift
    vshift_last_iter = .false.
    H_U_gap_crit = 0.02_dp   ! Small gap threshold (in a.u.)

    !==============================================================================
    ! Initialize SCF Calculation
    !==============================================================================
    call measure_time(print_total=1, log_unit=IW)

    ! Prepare orthogonalization matrix (S^-1/2)
    call matrix_invsqrt(smat, qmat, nbf)

    ! Compute Nuclear-Nuclear energy
    nenergy = e_charge_repulsion(infos%atoms%xyz, infos%atoms%zn - infos%basis%ecp_zn_num)

    ! During guess, the Hcore, Q nd Overlap matrices were formed.
    ! Using these, the initial orbitals (VEC) and density (Dmat) were subsequently computed.
    ! Now we are going to calculate ERI(electron repulsion integrals) to form a new FOCK
    ! matrix.

    ! Initialize ERI calculations and screening
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    call flush(IW)

    ! Initialize density matrices for integral evaluation
    select case (scf_type)
    case (scf_rhf)
      pdmat(:,1) = dmat_a
      allocate(int2_rhf_data_t :: int2_data)
      int2_data = int2_rhf_data_t(nfocks=1, &
                                  d=pdmat, &
                                  scale_exchange=scalefactor)
    case (scf_uhf, scf_rohf)
      pdmat(:,1) = dmat_a
      pdmat(:,2) = dmat_b
      allocate(int2_urohf_data_t :: int2_data)
      int2_data = int2_urohf_data_t(nfocks=2, &
                                    d=pdmat, &
                                    scale_exchange=scalefactor)
    end select

    ! Convert overlap matrix to full format for DIIS/SOSCF
    call unpack_matrix(smat, smat_full, nbf, 'U')

    !==============================================================================
    ! Configure SCF Convergence Accelerator (DIIS/SOSCF)
    !==============================================================================
    ! Configuration is determined by soscf_type, diis_type, and vshift:
    !
    ! soscf_type = 0 (Pure DIIS):
    !   ├── diis_type = 5 (V-DIIS):
    !   │   ├── vshift unset (0.0): Sets vshift = 0.1
    !   │   └── Uses: [C-DIIS, E-DIIS, C-DIIS]
    !   │       Thresholds: [2.0, 1.0, vdiis_cdiis_switch]
    !   ├── vshift set (non-zero):
    !   │   └── Uses: [C-DIIS, E-DIIS, C-DIIS]
    !   │       Thresholds: [2.0, 1.0, vshift_cdiis_switch]
    !   └── Otherwise:
    !       └── Uses: diis_type method (1=C-DIIS, 2=E-DIIS, 3=A-DIIS)
    !           Threshold: diis_method_threshold
    !
    ! soscf_type = 1 (Pure SOSCF):
    !   └── Uses: SOSCF
    !       Starts: From first iteration
    !
    ! soscf_type = 2 (DIIS + SOSCF):
    !   ├── diis_type = 5 (V-DIIS):
    !   │   ├── vshift unset (0.0): Sets vshift = 0.1
    !   │   └── Uses: [C-DIIS, E-DIIS, C-DIIS, SOSCF]
    !   │       Thresholds: [2.0, 1.0, vdiis_cdiis_switch, soscf_conv]
    !   └── Otherwise:
    !       └── Uses: [diis_type method, SOSCF]
    !           Thresholds: [diis_method_threshold, soscf_conv]
    !
    ! Additional Note:
    ! - MOM activates if mom = .true. and DIIS error < mom_switch, handled outside this block.
    !==============================================================================

    ! SOSCF options
    use_soscf = .false.

    ! DIIS options
    maxdiis = infos%control%maxdiis
    diis_error = 2.0_dp
    diis_name = [character(len=6) :: "none", "c-DIIS", "e-DIIS", "a-DIIS", "v-DIIS"]
    diis_reset = infos%control%diis_reset_mod

    ! Initialize SCF Convergence Accelerator
    select case (infos%control%soscf_type)
    case (0) ! Pure DIIS Accelerators
      ! Set up DIIS convergence accelerator
      if (infos%control%diis_type == 5) then
        ! v-DIIS setup: combination of E-DIIS and C-DIIS with vshift
        call conv%init(ldim=nbf, &
                       maxvec=maxdiis, &
                       subconvergers=[conv_cdiis, &
                                      conv_ediis, &
                                      conv_cdiis], &
                       thresholds   =[ethr_cdiis_big, &
                                      ethr_ediis, &
                                      infos%control%vdiis_cdiis_switch], &
                       overlap=smat_full, &
                       overlap_sqrt=qmat, &
                       num_focks=diis_nfocks, &
                       verbose=infos%control%verbose)
        if (infos%control%vshift == 0.0_dp) then
          infos%control%vshift = 0.1_dp
          vshift = 0.1_dp
          write(IW, '(X,A)') 'Setting Vshift = 0.1 a.u., since VDIIS is chosen without Vshift value.'
        end if
      elseif (infos%control%vshift /= 0.0_dp) then
        ! Custom vshift setup with E-DIIS and C-DIIS
        call conv%init(ldim=nbf, &
                       maxvec=maxdiis, &
                       subconvergers=[conv_cdiis, &
                                      conv_ediis, &
                                      conv_cdiis], &
                       thresholds   =[ethr_cdiis_big, &
                                      ethr_ediis, &
                                      infos%control%vshift_cdiis_switch], &
                       overlap=smat_full, &
                       overlap_sqrt=qmat, &
                       num_focks=diis_nfocks, &
                       verbose=infos%control%verbose)
      else
        ! Standard DIIS setup from input
        call conv%init(ldim=nbf, &
                       maxvec=maxdiis, &
                       subconvergers=[infos%control%diis_type], &
                       thresholds   =[infos%control%diis_method_threshold], &
                       overlap=smat_full, &
                       overlap_sqrt=qmat, &
                       num_focks=diis_nfocks, &
                       verbose=infos%control%verbose)
      end if

    case (1) ! Pure SOSCF Accelerator
      use_soscf = .true.
      ! Pure SOSCF strategy
      call conv%init(ldim=nbf, nelec_a=nelec_a, nelec_b=nelec_b, &
                     maxvec=infos%control%maxit, &
                     subconvergers=[conv_soscf], &
                     thresholds   =[huge(1.0_dp)], &  ! SOSCF runs from first iteration
                     overlap=smat_full, &
                     overlap_sqrt=qmat, &
                     num_focks=soscf_nfocks, &
                     scf_type=infos%control%scftype, &
                     verbose=infos%control%verbose)
      ! Configure the SOSCF converger with SOSCF input parameters
      call set_soscf_parametres(infos, conv)
    case (2) ! DIIS+SOSCF Accelerator
      use_soscf = .true.
      if (infos%control%diis_type == 5) then
        ! V-DIIS + SOSCF hybrid strategy
        call conv%init(ldim=nbf, nelec_a=nelec_a, nelec_b=nelec_b, &
                       maxvec=maxdiis, &
                       subconvergers=[conv_cdiis, &
                                      conv_ediis, &
                                      conv_cdiis, &
                                      conv_soscf], &
                       thresholds   =[ethr_cdiis_big, &
                                      ethr_ediis, &
                                      infos%control%vdiis_cdiis_switch], &
!                                      infos%control%soscf_conv], &
                       overlap=smat_full, &
                       overlap_sqrt=qmat, &
                       num_focks=diis_nfocks, &
                       verbose=infos%control%verbose)
        if (infos%control%vshift == 0.0_dp) then
          infos%control%vshift = 0.1_dp
          vshift = 0.1_dp
          write(IW, '(X,A)') 'Setting Vshift = 0.1 a.u., since VDIIS is chosen without Vshift value.'
        end if
      else
        ! Regular DIIS + SOSCF hybrid strategy
        call conv%init(ldim=nbf, nelec_a=nelec_a, nelec_b=nelec_b, &
                       maxvec=maxdiis, &
                       subconvergers=[infos%control%diis_type, &
                                      conv_soscf], &
                       thresholds   =[infos%control%diis_method_threshold], &
!                                      infos%control%soscf_conv], &
                       overlap=smat_full, &
                       overlap_sqrt=qmat, &
                       num_focks=diis_nfocks, &
                       verbose=infos%control%verbose)
      end if

      ! Configure the SOSCF converger with SOSCF input parameters
      call set_soscf_parametres(infos, conv)
    end select


    ! Initialize DFT exchange-correlation energy
    eexc = 0.0_dp
    e_old = 0.0_dp

    !==============================================================================
    ! Print SCF Options
    !==============================================================================
    write(IW,'(/5X,"SCF options"/ &
               &5X,18("-")/ &
               &5X,"SCF type = ",A,5x,"MaxIT = ",I5/, &
               &5X,"MaxDIIS = ",I5,17x,"Conv = ",F14.10/, &
               &5X,"DIIS Type = ",A/, &
               &5X,"vDIIS_cDIIS_Switch = ",F8.5,3x,"vDIIS_vshift_Switch = ",F8.5/, &
               &5X,"DIIS Reset Mod = ",I5,10x,"DIIS Reset Conv = ",F12.8/, &
               &5X,"VShift = ",F8.5,15X,"VShift_cDIIS_Switch = ",F8.5)') &
               & scf_name, infos%control%maxit, &
               & infos%control%maxdiis, infos%control%conv, &
               & diis_name(infos%control%diis_type), &
               & infos%control%vdiis_cdiis_switch, infos%control%vdiis_vshift_switch, &
               & infos%control%diis_reset_mod, infos%control%diis_reset_conv, &
               & infos%control%vshift, infos%control%vshift_cdiis_switch
    write(IW,'(5X,"MOM = ",L5,21X,"MOM_Switch = ",F8.5)') &
               & infos%control%mom, infos%control%mom_switch
    write(IW,'(5X,"pFON = ",L5,20X,"pFON Start Temp. = ",F9.2,/, &
               5X, "pFON Cooling Rate = ", F9.2,2X," pFON Num. Smearing = ",F8.5)') &
               infos%control%pfon, infos%control%pfon_start_temp, &
               infos%control%pfon_cooling_rate, infos%control%pfon_nsmear
    if (use_soscf) then
      write(IW,'(/5X,"SOSCF options"/ &
                 &5X,18("-"))')
      write(IW,'(5X,"SOSCF enabled = ",L5,16X,"SOSCF type = ",I1)') &
             & use_soscf, infos%control%soscf_type
    end if

    ! Initial message for SCF iterations
    if (infos%control%pfon) then
      write(IW,fmt="&
            &(/3x,'Direct SCF iterations begin.'/, &
            &  3x,113('='),/ &
            &  4x,'Iter',9x,'Energy',12x,'Delta E',9x,'Int Skip',5x,'DIIS Error',5x,'Shift',5x,'Method',5x,'pFON'/ &
            &  3x,113('='))")
    elseif (infos%control%soscf_type == 1) then
      write(IW,fmt="&
            &(/3x,'Direct SCF iterations begin.'/, &
            &  3x,107('='),/ &
            &  4x,'Iter',9x,'Energy',12x,'Delta E',9x,'Int Skip',5x,'Grad. RMS',6x,'Den. RMS',7x,'Shift',5x,'Method'/ &
            &  3x,107('='))")
    else
      write(IW,fmt="&
            &(/3x,'Direct SCF iterations begin.'/, &
            &  3x,93('='),/ &
            &  4x,'Iter',9x,'Energy',12x,'Delta E',9x,'Int Skip',5x,'DIIS Error',5x,'Shift',5x,'Method'/ &
            &  3x,93('='))")
    end if
    call flush(IW)

    !==============================================================================
    ! Begin Main SCF Iteration Loop
    !==============================================================================
    do iter = 1, maxit
      if (do_rstctmo) then
        mo_energy_a_for_rstctmo = mo_energy_a
        mo_a_for_rstctmo = mo_a
      end if
      !----------------------------------------------------------------------------
      ! Update pFON Temperature (if enabled)
      !----------------------------------------------------------------------------
      if (do_pfon) then
        if (do_pfon_final) then
          ! Extra iteration: set temperature to 1K
          pfon%temp = 1.0_dp
          pfon%beta = 1.0_dp / (kB_HaK * pfon%temp)
          write(IW, "(10x, 'Extra SCF iteration with Temp = 1K')")
        else
          ! Normal temperature adjustment
          call pfon%adjust_temperature(iter, maxit, diis_error, infos%control%conv)
        end if
      end if

      !----------------------------------------------------------------------------
      ! Initialize Fock Matrices for Current Iteration
      !----------------------------------------------------------------------------
      pfock = 0.0_dp

      !----------------------------------------------------------------------------
      ! Compute Difference Density Matrix for Incremental Fock Build
      !----------------------------------------------------------------------------
      ! This is the main advantage of direct SCF - allows better ERI screening
      if (infos%control%scf_incremental /= 0) then
        pdmat = pdmat - dold
      end if

      !----------------------------------------------------------------------------
      ! Compute Two-Electron Contribution to Fock Matrix
      !----------------------------------------------------------------------------
      call int2_driver%run(int2_data, &
                           cam=is_dft.and.infos%dft%cam_flag, &
                           alpha=infos%dft%cam_alpha, &
                           beta=infos%dft%cam_beta,&
                           mu=infos%dft%cam_mu)
      nschwz = int2_driver%skipped

      !----------------------------------------------------------------------------
      ! Recover Full Fock and Density from Difference Matrices (if incremental)
      !----------------------------------------------------------------------------
      if (infos%control%scf_incremental /= 0) then
        pdmat = pdmat + dold
        int2_data%f(:,:,1) = int2_data%f(:,:,1) + fold
        fold = int2_data%f(:,:,1)
        dold = pdmat
      end if

      !----------------------------------------------------------------------------
      ! Scale Two-Electron Terms in Fock Matrix
      !----------------------------------------------------------------------------
      ! Scale off-diagonal elements by 0.5, diagonal by 1.0
      pfock(:,:) = 0.5_dp * int2_data%f(:,:,1)
      ii = 0
      do i = 1, nbf
        ii = ii + i
        pfock(ii,1:nfocks) = 2.0_dp * pfock(ii,1:nfocks)
      end do

      !----------------------------------------------------------------------------
      ! Add One-Electron Core Hamiltonian to Fock Matrix
      !----------------------------------------------------------------------------
      do i = 1, nfocks
        pfock(:,i) = pfock(:,i) + hcore
      end do

      !----------------------------------------------------------------------------
      ! Compute HF Energy Components
      !----------------------------------------------------------------------------
      ehf = 0.0_dp
      ehf1 = 0.0_dp

      ! Compute one and two-electron energies
      do i = 1, nfocks
        ehf1 = ehf1 + traceprod_sym_packed(pdmat(:,i), hcore, nbf)
        ehf = ehf + traceprod_sym_packed(pdmat(:,i), pfock(:,i), nbf)
      end do

      ! Total HF energy = 0.5*(E1 + E2) (to avoid double-counting)
      ehf = 0.5_dp * (ehf + ehf1)
      etot = ehf + nenergy

      !----------------------------------------------------------------------------
      ! Compute DFT Exchange-Correlation Contribution (if DFT)
      !----------------------------------------------------------------------------
      if (is_dft) then
        if (scf_type == scf_rhf) then
          call dftexcor(basis, molgrid, 1, pfxc, pfxc, mo_a, mo_a, &
                        nbf, nbf_tri, eexc, totele, totkin, infos)
        else if (scf_type == scf_uhf) then
          call dftexcor(basis, molgrid, 2, pfxc(:,1), pfxc(:,2), mo_a, mo_b, &
                        nbf, nbf_tri, eexc, totele, totkin, infos)
        else if (scf_type == scf_rohf) then
          ! ROHF does not have MO_B. So we copy MO_A to MO_B.
          mo_b = mo_a
          call dftexcor(basis, molgrid, 2, pfxc(:,1), pfxc(:,2), mo_a, mo_b, &
                        nbf, nbf_tri, eexc, totele, totkin, infos)
        end if

        ! Add XC contribution to Fock and total energy
        pfock = pfock + pfxc
        etot = etot + eexc
      end if

      !----------------------------------------------------------------------------
      ! Form Special ROHF Fock Matrix and Apply Vshift (if ROHF calculation)
      !----------------------------------------------------------------------------
      if (scf_type == scf_rohf) then
        ! Store the original alpha Fock matrix before ROHF transformation
        ! This is needed to preserve it for energy evaluation and printing
        rohf_bak(:,1) = pfock(:,1)
        rohf_bak(:,2) = pfock(:,2)

        ! Turn off level shifting for the final iteration if requested
        if (vshift_last_iter) vshift = 0.0_dp

        ! Apply the Guest-Saunders ROHF Fock transformation
        ! This creates a modified Fock matrix with proper coupling between
        ! closed-shell, open-shell, and virtual orbital spaces
        call form_rohf_fock(pfock(:,1),pfock(:,2), mo_a, smat_full, &
                            nelec_a, nelec_b, nbf, vshift, work1, work2)

        ! Combine alpha and beta densities for ROHF
        ! This is needed because ROHF uses a single set of MOs for both spins,
        ! so we need the total density for the next iteration
        pdmat(:,1) = pdmat(:,1) + pdmat(:,2)
      end if

      !----------------------------------------------------------------------------
      ! Apply Vshift for RHF/UHF (if enabled)
      !----------------------------------------------------------------------------
      if (vshift > 0.0_dp .and. scf_type /= scf_rohf) then
        ! Turn off level shifting for the final iteration if requested
        if (vshift_last_iter) vshift = 0.0_dp

        ! Apply level shifting based on SCF type
        select case (scf_type)
        case (scf_rhf)
          ! RHF: One Fock matrix with doubly occupied orbitals
          call level_shift_fock(pfock(:,1), mo_a, smat_full, nelec/2, nbf, vshift, &
                                work1, work2)

        case (scf_uhf)
          ! UHF: Two Fock matrices with separate alpha and beta occupations
          call level_shift_fock(pfock(:,1), mo_a, smat_full, nelec_a, nbf, vshift, &
                                work1, work2)
          if (nelec_b > 0) then
            call level_shift_fock(pfock(:,2), mo_b, smat_full, nelec_b, nbf, vshift, &
                                  work1, work2)
          end if
        end select
      end if

      !----------------------------------------------------------------------------
      ! Pass Fock and Density to Convergence Accelerator
      !----------------------------------------------------------------------------
      if (use_soscf) then
        ! SOSCF: Pass MO coefficients and energies for orbital rotation
        ! Note: Fock matrix is fixed in SOSCF;
        !       rebuilt in next iteration,
        !       Fock is not retrieved here.
        if (do_pfon) then
          select case (scf_type)
          case (scf_rhf)
            call conv%add_data( &
                     f=pfock(:,1:diis_nfocks), &
                     dens=pdmat(:,1:diis_nfocks), &
                     e=etot, &
                     mo_a=mo_a, &
                     mo_e_a=mo_energy_a, &
                     pfon=pfon)
          case (scf_rohf)
            call conv%add_data( &
                     f=pfock(:,1:soscf_nfocks), &
                     dens=pdmat(:,1:soscf_nfocks), &
                     e=etot, &
                     mo_a=mo_a, &
                     mo_b=mo_b, &
                     mo_e_a=mo_energy_a, &
                     mo_e_b=mo_energy_b, &
                     pfon=pfon)
          case (scf_uhf)
            call conv%add_data( &
                     f=pfock(:,1:diis_nfocks), &
                     dens=pdmat(:,1:diis_nfocks), &
                     e=etot, &
                     mo_a=mo_a, &
                     mo_b=mo_b, &
                     mo_e_a=mo_energy_a, &
                     mo_e_b=mo_energy_b, &
                     pfon=pfon)
          end select
        else
          select case (scf_type)
          case (scf_rhf)
            call conv%add_data( &
                     f=pfock(:,1:diis_nfocks), &
                     dens=pdmat(:,1:diis_nfocks), &
                     e=etot, &
                     mo_a=mo_a, &
                     mo_e_a=mo_energy_a)
          case (scf_rohf)
            call conv%add_data( &
                     f=rohf_bak(:,1:soscf_nfocks), &
                     dens=pdmat(:,1:soscf_nfocks), &
                     e=etot, &
                     mo_a=mo_a, &
                     mo_b=mo_b, &
                     mo_e_a=mo_energy_a, &
                     mo_e_b=mo_energy_b)
          case (scf_uhf)
            call conv%add_data( &
                     f=pfock(:,1:diis_nfocks), &
                     dens=pdmat(:,1:diis_nfocks), &
                     e=etot, &
                     mo_a=mo_a, &
                     mo_b=mo_b, &
                     mo_e_a=mo_energy_a, &
                     mo_e_b=mo_energy_b)
          end select
        end if
      else
        ! DIIS: Only pass Fock and density matrices
        call conv%add_data( &
                 f=pfock(:,1:diis_nfocks), &
                 dens=pdmat(:,1:diis_nfocks), &
                 e=etot)
      end if

      !----------------------------------------------------------------------------
      ! Run Convergence Accelerator (DIIS/SOSCF)
      !----------------------------------------------------------------------------
      call conv%run(conv_res)
      diis_error = conv_res%get_error()
      if (use_soscf) then
        rms_grad = conv_res%get_rms_grad()
        rms_dp = conv_res%get_rms_dp()
      end if


      !----------------------------------------------------------------------------
      ! Print Current Energy
      !----------------------------------------------------------------------------
      ! Print iteration information
      if (infos%control%pfon) then
        write(IW,fmt="(4x,i4.1,2x,f17.10,1x,f17.10,1x,i13,1x,f14.8,5x,f5.3,5x,a,5x,a,f9.2)") &
              iter, etot, etot - e_old, nschwz, diis_error, vshift, &
              trim(conv_res%active_converger_name), "Temp:", pfon%temp
        write(IW,fmt="(100x,a,f9.2)") "Beta:", pfon%beta
      elseif (infos%control%soscf_type == 1) then
        write(IW,'(4x,i4.1,2x,f17.10,1x,f17.10,1x,i13,1x,f14.8,1x,f14.8,5x,f5.3,5x,a)') &
              iter, etot, etot - e_old, nschwz, rms_grad, rms_dp, vshift, &
              trim(conv_res%active_converger_name)
      else
        write(IW,'(4x,i4.1,2x,f17.10,1x,f17.10,1x,i13,1x,f14.8,5x,f5.3,5x,a)') &
              iter, etot, etot - e_old, nschwz, diis_error, vshift, &
              trim(conv_res%active_converger_name)
      end if
      call flush(IW)

      !----------------------------------------------------------------------------
      ! Update VDIIS Parameters (if using VDIIS)
      !----------------------------------------------------------------------------
      if ((infos%control%diis_type == 5) .and. &
          (diis_error < infos%control%vdiis_vshift_switch)) then
        vshift = 0.0_dp
      else if ((infos%control%diis_type == 5) .and. &
               (diis_error >= infos%control%vdiis_vshift_switch)) then
        vshift = infos%control%vshift
      end if

      e_old = etot

      !----------------------------------------------------------------------------
      ! Check for SCF Convergence
      !----------------------------------------------------------------------------
      if ((abs(diis_error) < infos%control%conv) .and. (vshift == 0.0_dp)) then
        ! Fully converged - exit loop
        if (do_pfon) then
          if (pfon%temp > 1.0_dp + 1.0e-6_dp) then
            do_pfon_final = .true.
          else
            exit
          end if
        else
          exit
        end if
      elseif ((abs(diis_error) < infos%control%conv) .and. (vshift /= 0.0_dp)) then
        ! Converged but need one more iteration with vshift=0
        write(IW,"(3x,64('-')/10x,'Performing a last SCF with zero VSHIFT.')")
        vshift_last_iter = .true.
      elseif (vshift_last_iter) then
        ! Only for ROHF case the final iteration with vshift=0 complete - exit loop
        call get_ab_initio_orbital(pfock(:,1), mo_a, mo_energy_a, qmat)
        exit
      end if

      !----------------------------------------------------------------------------
      ! Reset DIIS
      !----------------------------------------------------------------------------
      ! Vshift=0.0 and slow cases
      ! Check if DIIS reset is needed
      diis_reset_condition = (((iter/diis_reset) >= 1) .and. &
                              (modulo(iter,diis_reset) == 0) .and. &
                              (diis_error > infos%control%diis_reset_conv) .and. &
                              (infos%control%vshift == 0.0_dp) .and. &
                              (infos%control%soscf_type /= 1))

      if (diis_reset_condition) then
        ! Resetting DIIS for difficult cases
        write(IW,"(3x,64('-')/10x,'Resetting DIIS.')")
        call conv_res%get_fock(matrix=pfock(:,1:diis_nfocks), istat=stat)
        call conv%init(ldim=nbf, &
                       maxvec=maxdiis, &
                       subconvergers=[conv_cdiis], &
                       thresholds   =[ethr_cdiis_big], &
                       overlap=smat_full, &
                       overlap_sqrt=qmat, &
                       num_focks=diis_nfocks, &
                       verbose=infos%control%verbose)
        ! After resetting DIIS, we need to skip SD
        call conv%add_data(f=pfock(:,1:diis_nfocks), &
                           dens=pdmat(:,1:diis_nfocks), &
                           e=Etot)
        call conv%run(conv_res)
      end if

      !----------------------------------------------------------------------------
      ! Update Fock or Orbitals and Eigenvalues Based on Active Converger
      !----------------------------------------------------------------------------
      if (use_soscf .and. trim(conv_res%active_converger_name) == 'SOSCF') then
        ! SOSCF: Retrieve updated MOs and energies directly
        ! Note: Fock matrix is fixed in SOSCF;
        !       rebuilt in next iteration,
        !       Fock is not retrieved here.
        if (int2_driver%pe%rank == 0) then
          ! Retrieve updated Alpha Orbitals
          call conv_res%get_mo_a(mo_a, istat=stat)
          ! Retrieve updated Energies of Alpha Orbitals
          call conv_res%get_mo_e_a(mo_energy_a, istat=stat)
          if (scf_type == scf_uhf .and. nelec_b /= 0) then
            ! Retrieve updated Beta Orbitals and its Energies
            call conv_res%get_mo_b(mo_b, stat)
            call conv_res%get_mo_e_b(mo_energy_b, stat)
          elseif (scf_type == scf_rohf) then
            mo_b = mo_a
            mo_energy_b = mo_energy_a
          end if
          if (stat /= 0) then
            call show_message('Error retrieving SOSCF results', WITH_ABORT)
          end if
        end if
      else
        ! DIIS: Retrieve updated Fock directly
        ! Form the interpolated the Fock/Density matrix
        call conv_res%get_fock(matrix=pfock(:,1:diis_nfocks), istat=stat)
        if (int2_driver%pe%rank == 0) then
           ! Compute New Alpha Orbitals
           call get_ab_initio_orbital(pfock(:,1), mo_a, mo_energy_a, qmat)
           if (scf_type == scf_uhf .and. nelec_b /= 0) then
              ! Only UHF has beta orbitals.
              call get_ab_initio_orbital(pfock(:,2), mo_b, mo_energy_b, qmat)
           end if
        end if
      end if

      ! Broadcast updated orbitals to all processes
      call int2_driver%pe%bcast(mo_a, size(mo_a))
      call int2_driver%pe%bcast(mo_energy_a, size(mo_energy_a))
      if (scf_type == scf_uhf .and. nelec_b /= 0) then
        call int2_driver%pe%bcast(mo_b, size(mo_b))
        call int2_driver%pe%bcast(mo_energy_b, size(mo_energy_b))
      end if

      !----------------------------------------------------------------------------
      ! Calculate pFON Occupations (if enabled)
      !----------------------------------------------------------------------------
      if (do_pfon) then
        if (scf_type == scf_uhf) then
          call pfon%compute_occupations(mo_energy_a, mo_energy_b)
        else
          call pfon%compute_occupations(mo_energy_a)
        end if
      end if

      !----------------------------------------------------------------------------
      ! Apply MOM (Maximum Overlap Method) if enabled
      !----------------------------------------------------------------------------
      if (do_mom) then
        ! Check if MOM should be activated based on convergence
        if (diis_error < infos%control%mom_switch) then
          mom_active = .true.
        end if

        ! Apply MOM if active and not the first iteration
        if (mom_active .and. .not. initial_mom_iter) then
          ! Apply MOM for alpha spin channel
          call apply_mom(infos, mo_a_prev, mo_e_a_prev, &
                         mo_a, mo_energy_a, smat_full, nelec_a, &
                         "Alpha", work1, work2)

          ! Apply MOM for beta spin channel (UHF only)
          if (scf_type == scf_uhf .and. nelec_b > 0) then
            call apply_mom(infos, mo_b_prev, mo_e_b_prev, &
                           mo_b, mo_energy_b, smat_full, nelec_b, &
                           "Beta", work1, work2)
          end if
        end if

        ! Store current orbitals for next iteration
        mo_a_prev = mo_a
        mo_e_a_prev = mo_energy_a

        if (scf_type == scf_uhf) then
           mo_b_prev = mo_b
           mo_e_b_prev = mo_energy_b
        end if

        initial_mom_iter = .false.
      end if
      !----------------------------------------------------------------------------
      ! Apply XAS if enabled
      !----------------------------------------------------------------------------
      if (do_rstctmo) then
        call apply_mom(infos, mo_a_for_rstctmo, mo_energy_a_for_rstctmo, &
                       mo_a, mo_energy_a, smat_full, nelec_a, "Alpha", work1, work2)
      end if
      !----------------------------------------------------------------------------
      ! Build New Density Matrix from Updated Orbitals
      !----------------------------------------------------------------------------
      if (int2_driver%pe%rank == 0) then
        if (do_pfon) then
          ! pFON density build with fractional occupations
          select case(scf_type)
          case (scf_rhf)
            call pfon%build_density(pdmat(:,1), mo_a, work1, work2)
          case (scf_uhf)
            call pfon%build_density(pdmat(:,1), mo_a, work1, work2, pdmat(:,2), mo_b)
          case (scf_rohf)
            call pfon%build_density(pdmat(:,1), mo_a, work1, work2, pdmat(:,2), mo_b)
          end select
        else
          ! Standard density build with integer occupations
          call get_ab_initio_density(pdmat(:,1),mo_a,pdmat(:,2),mo_b,infos,basis)
        end if
      end if
      call int2_driver%pe%bcast(pdmat, size(pdmat))

      !----------------------------------------------------------------------------
      ! Check HOMO-LUMO Gap for Convergence Prediction
      !----------------------------------------------------------------------------
      if ((iter > 10)) then
        ! Calculate HOMO-LUMO gap
        select case (scf_type)
        case (scf_rhf)
          H_U_gap = mo_energy_a(nelec/2) - mo_energy_a(nelec/2-1)
        case (scf_uhf)
          H_U_gap = mo_energy_a(nelec_a+1) - mo_energy_a(nelec_a)
          if (nelec_b > 0) then
            H_U_gap = min(H_U_gap, &
                          mo_energy_b(nelec_b+1) - mo_energy_b(nelec_b))
          end if
        case (scf_rohf)
          H_U_gap = mo_energy_a(nelec_a+1) - mo_energy_a(nelec_a)
        end select

        if (H_U_gap < H_U_gap_crit .and. vshift>0) then
          ! Small gap detected, enable level shifting
          write(IW,"(3x,64('-')/10x,&
                   &'Small HOMO-LUMO gap detected (',&
                   &F8.5,' au). Apply level vshift=',F6.4)") H_U_gap, vshift
        end if
      end if

    ! End of Main SCF Iteration Loop
    end do

    !----------------------------------------------------------------------------
    ! Clean Convergence Accelerator (DIIS/SOSCF)
    !----------------------------------------------------------------------------
    call conv%clean()

    !==============================================================================
    ! Post-SCF Processing and Final Output
    !==============================================================================

    !----------------------------------------------------------------------------
    ! Report SCF Convergence Status
    !----------------------------------------------------------------------------
    if (iter > maxit) then
      write(IW,"(3x,64('-')/10x,'SCF is not converged ....')")
      infos%mol_energy%SCF_converged = .false.
    else
      write(IW,"(3x,64('-')/10x,'SCF convergence achieved ....')")
      infos%mol_energy%SCF_converged = .true.
    end if

    write(IW,"(/' Final ',A,' energy is',F20.10,' after',I4,' iterations'/)") trim(scf_name), etot, iter

    !----------------------------------------------------------------------------
    ! Print DFT-Specific Information (if DFT)
    !----------------------------------------------------------------------------
    if (is_dft) then
      write(IW,*)
      write(IW,"(' DFT: XC energy              = ',F20.10)") eexc
      write(IW,"(' DFT: total electron density = ',F20.10)") totele
      write(IW,"(' DFT: number of electrons    = ',I9,/)") nelec
    end if

    !----------------------------------------------------------------------------
    ! Broadcast Final MOs and Energies to All Processes
    !----------------------------------------------------------------------------
    call int2_driver%pe%bcast(pdmat, size(pdmat))
    call int2_driver%pe%bcast(mo_a, size(mo_a))
    call int2_driver%pe%bcast(mo_energy_a, size(mo_energy_a))
    if (scf_type == scf_uhf .and. nelec_b /= 0) then
        call int2_driver%pe%bcast(mo_b, size(mo_b))
        call int2_driver%pe%bcast(mo_energy_b, size(mo_energy_b))
    end if

    !----------------------------------------------------------------------------
    ! Save Final Fock and Density Matrices
    !----------------------------------------------------------------------------
    select case (scf_type)
    case (scf_rhf)
      fock_a = pfock(:,1)
      dmat_a = pdmat(:,1)
    case (scf_uhf)
      fock_a = pfock(:,1)
      fock_b = pfock(:,2)
      dmat_a = pdmat(:,1)
      dmat_b = pdmat(:,2)
    case (scf_rohf)
      fock_a = rohf_bak(:,1)
      call mo_to_ao(fock_b, pfock(:,2), smat_full, mo_a, nbf, nbf, work1, work2)
      dmat_a = pdmat(:,1) - pdmat(:,2)
      dmat_b = pdmat(:,2)
      mo_b = mo_a
      mo_energy_b = mo_energy_a
    end select

    !----------------------------------------------------------------------------
    ! Print Molecular Orbitals
    !----------------------------------------------------------------------------
    call print_mo_range(basis, infos, mostart=1, moend=nbf)

    !----------------------------------------------------------------------------
    ! Calculate Final Energy Components
    !----------------------------------------------------------------------------
    psinrm = 0.0_dp
    tkin = 0.0_dp
    do i = 1, diis_nfocks
      psinrm = psinrm + traceprod_sym_packed(pdmat(:,i),smat,nbf)/nelec
      tkin = tkin + traceprod_sym_packed(pdmat(:,i),tmat,nbf)
    end do

    ! Calculate energy components
    vne = ehf1 - tkin
    vee = etot - ehf1 - nenergy
    vnn = nenergy
    vtot = vne + vnn + vee
    virial = - vtot/tkin

    !----------------------------------------------------------------------------
    ! Print Final Energy Components
    !----------------------------------------------------------------------------
    call print_scf_energy(psinrm, ehf1, nenergy, etot, vee, vne, vnn, vtot, tkin, virial)

    !----------------------------------------------------------------------------
    ! Save Results to infos Structure
    !----------------------------------------------------------------------------
    infos%mol_energy%energy = etot
    infos%mol_energy%psinrm = psinrm
    infos%mol_energy%ehf1 = ehf1
    infos%mol_energy%vee = vee
    infos%mol_energy%nenergy = nenergy
    infos%mol_energy%vne = vne
    infos%mol_energy%vnn = vnn
    infos%mol_energy%vtot = vtot
    infos%mol_energy%tkin = tkin
    infos%mol_energy%virial = virial
    infos%mol_energy%energy = etot

    !----------------------------------------------------------------------------
    ! Clean Up Resources
    !----------------------------------------------------------------------------
    call int2_driver%clean()
    call measure_time(print_total=1, log_unit=IW)

    !----------------------------------------------------------------------------
    ! Clean up in the finalization section
    !----------------------------------------------------------------------------
    if (do_mom) then
       if (allocated(mo_a_prev))   deallocate(mo_a_prev)
       if (allocated(mo_e_a_prev)) deallocate(mo_e_a_prev)
       if (allocated(mo_b_prev))   deallocate(mo_b_prev)
       if (allocated(mo_e_b_prev)) deallocate(mo_e_b_prev)
    end if
  end subroutine scf_driver

  !> @brief Prints the final energy components of the SCF calculation.
  !> @detail Outputs a detailed breakdown of energy terms, including one-electron,
  !>         two-electron, nuclear repulsion, and total energies, as well as potential
  !>         (electron-electron, nucleus-electron, nucleus-nucleus, total) and
  !>         kinetic contributions, and the virial ratio.
  !> @param[in] psinrm Wavefunction normalization factor.
  !> @param[in] ehf1 One-electron energy.
  !> @param[in] enuclear Nuclear repulsion energy.
  !> @param[in] etot Total SCF energy.
  !> @param[in] vee Electron-electron potential energy.
  !> @param[in] vne Nucleus-electron potential energy.
  !> @param[in] vnn Nucleus-nucleus potential energy.
  !> @param[in] vtot Total potential energy.
  !> @param[in] tkin Kinetic energy.
  !> @param[in] virial Virial ratio (V/T).
  subroutine print_scf_energy(psinrm, ehf1, enuclear, etot, vee, vne, vnn, vtot, tkin, virial)
     use precision, only: dp
     use io_constants, only: iw
     implicit none
     real(kind=dp) :: psinrm, ehf1, enuclear, etot, vee, vne, vnn, vtot, tkin, virial

     write(IW,"(/10X,17('=')/10X,'Energy components'/10X,17('=')/)")
     write(IW,"('         Wavefunction normalization =',F19.10)") psinrm
     write(IW,*)
     write(IW,"('                One electron energy =',F19.10)") ehf1
     write(IW,"('                Two electron energy =',F19.10)") vee
     write(IW,"('           Nuclear repulsion energy =',F19.10)") enuclear
     write(IW,"(38X,18('-'))")
     write(IW,"('                       TOTAL energy =',F19.10)") etot
     write(IW,*)
     write(IW,"(' Electron-electron potential energy =',F19.10)") vee
     write(IW,"('  Nucleus-electron potential energy =',F19.10)") vne
     write(IW,"('   Nucleus-nucleus potential energy =',F19.10)") vnn
     write(IW,"(38X,18('-'))")
     write(IW,"('             TOTAL potential energy =',F19.10)") vtot
     write(IW,"('               TOTAL kinetic energy =',F19.10)") tkin
     write(IW,"('                 Virial ratio (V/T) =',F19.10)") virial
     write(IW,*)
  end subroutine print_scf_energy

  !> @In this implementation, we don’t need these parameters— they were added,
  !> @but they appear to be unnecessary right now.
  !> @brief Configures parameters for the Second-Order SCF (SOSCF) convergence accelerator.
  !> @detail Sets SOSCF-specific parameters.
  !> @author Konstantin Komarov, 2023
  !> @param[in] infos System information and control parameters.
  !> @param[inout] conv SCF convergence driver object.
  subroutine set_soscf_parametres(infos, conv)
    use types, only: information
    use scf_converger, only: scf_conv, soscf_converger

    type(information), target, intent(inout) :: infos
    type(scf_conv) :: conv

    integer :: i

    ! Through accessing the SOSCF converger set its parameters:
    do i = lbound(conv%sconv, 1), ubound(conv%sconv, 1)
      select type (sc => conv%sconv(i)%s)
        type is (soscf_converger)
          sc%level_shift = infos%control%soscf_lvl_shift
          sc%soscf_reset_mod = infos%control%soscf_reset_mod
      end select
    end do

  end subroutine set_soscf_parametres

  !> @brief Forms the ROHF Fock matrix in the MO basis using the Guest-Saunders method.
  !> @detail Transforms alpha and beta Fock matrices from the AO basis to the MO basis,
  !>         constructs the ROHF Fock matrix following the Guest-Saunders approach,
  !>         and optionally applies a level shift to virtual orbitals.
  !>         Reference: M. F. Guest, V. Saunders. Mol. Phys. 28, 819 (1974).
  !> @author Konstantin Komarov, 2023
  !> @param[inout] fock_a_ao Alpha Fock matrix in AO basis (triangular format).
  !> @param[inout] fock_b_ao Beta Fock matrix in AO basis (triangular format).
  !> @param[in] mo_a Alpha MO coefficients.
  !> @param[in] smat_full Full overlap matrix.
  !> @param[in] nocca Number of occupied alpha orbitals.
  !> @param[in] noccb Number of occupied beta orbitals.
  !> @param[in] nbf Number of basis functions.
  !> @param[in] vshift Level shift parameter for virtual orbitals.
  !> @param[inout] work1 Work array 1 (nbf x nbf).
  !> @param[inout] work2 Work array 2 (nbf x nbf).
  subroutine form_rohf_fock(fock_a_ao, fock_b_ao, &
                            mo_a, smat_full, &
                            nocca, noccb, nbf, vshift, &
                            work1, work2)
    use precision, only: dp
    use mathlib, only: orthogonal_transform_sym, &
                       orthogonal_transform2, &
                       unpack_matrix, &
                       pack_matrix

    implicit none

    real(kind=dp), intent(inout), dimension(:) :: fock_a_ao
    real(kind=dp), intent(inout), dimension(:) :: fock_b_ao
    real(kind=dp), intent(in), dimension(:,:) :: mo_a
    real(kind=dp), intent(in), dimension(:,:) :: smat_full
    real(kind=dp), intent(inout), dimension(:,:) :: work1
    real(kind=dp), intent(inout), dimension(:,:) :: work2
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: vshift

    real(kind=dp), allocatable, dimension(:) :: fock_mo
    real(kind=dp), allocatable, dimension(:,:) :: &
          work_matrix, fock, fock_a, fock_b
    real(kind=dp) :: acc, aoo, avv, bcc, boo, bvv
    integer :: i, nbf_tri

    acc = 0.5_dp; aoo = 0.5_dp; avv = 0.5_dp
    bcc = 0.5_dp; boo = 0.5_dp; bvv = 0.5_dp
    nbf_tri = nbf * (nbf + 1) / 2

    ! Allocate full matrices
    allocate(work_matrix(nbf, nbf), &
             fock(nbf, nbf), &
             fock_mo(nbf_tri), &
             fock_a(nbf, nbf), &
             fock_b(nbf, nbf), &
             source=0.0_dp)

    ! Transform alpha and beta Fock matrices to MO basis
    call orthogonal_transform_sym(nbf, nbf, fock_a_ao, mo_a, nbf, fock_mo)
    fock_a_ao(:nbf_tri) = fock_mo(:nbf_tri)

    call orthogonal_transform_sym(nbf, nbf, fock_b_ao, mo_a, nbf, fock_mo)
    fock_b_ao(:nbf_tri) = fock_mo(:nbf_tri)

    ! Unpack triangular matrices to full matrices
    call unpack_matrix(fock_a_ao, fock_a)
    call unpack_matrix(fock_b_ao, fock_b)

    ! Construct ROHF Fock matrix in MO basis using Guest-Saunders method
    associate ( na => nocca &
              , nb => noccb &
      )
      fock(1:nb, 1:nb) = acc * fock_a(1:nb, 1:nb) &
                       + bcc * fock_b(1:nb, 1:nb)
      fock(nb+1:na, nb+1:na) = aoo * fock_a(nb+1:na, nb+1:na) &
                             + boo * fock_b(nb+1:na, nb+1:na)
      fock(na+1:nbf, na+1:nbf) = avv * fock_a(na+1:nbf, na+1:nbf) &
                               + bvv * fock_b(na+1:nbf, na+1:nbf)
      fock(1:nb, nb+1:na) = fock_b(1:nb, nb+1:na)
      fock(nb+1:na, 1:nb) = fock_b(nb+1:na, 1:nb)
      fock(1:nb, na+1:nbf) = 0.5_dp * (fock_a(1:nb, na+1:nbf) &
                                     + fock_b(1:nb, na+1:nbf))
      fock(na+1:nbf, 1:nb) = 0.5_dp * (fock_a(na+1:nbf, 1:nb) &
                                     + fock_b(na+1:nbf, 1:nb))
      fock(nb+1:na, na+1:nbf) = fock_a(nb+1:na, na+1:nbf)
      fock(na+1:nbf, nb+1:na) = fock_a(na+1:nbf, nb+1:na)

      ! Apply Vshift to the diagonal
      do i = nb+1, na
        fock(i,i) = fock(i,i) + vshift * 0.5_dp
      end do
      do i = na+1, nbf
        fock(i,i) = fock(i,i) + vshift
      end do
    end associate

    ! Back-transform ROHF Fock matrix to AO basis
    call dsymm('l', 'u', nbf, nbf, &
               1.0_dp, smat_full, nbf, &
                       mo_a, nbf, &
               0.0_dp, work1, nbf)
    call orthogonal_transform2('t', nbf, nbf, work1, nbf, fock, nbf, &
                               work_matrix, nbf, work2)
    call pack_matrix(work_matrix, fock_a_ao)

    deallocate(work_matrix, fock, fock_mo, fock_a, fock_b)

  end subroutine form_rohf_fock

  !> @brief Back-transforms a symmetric operator from the MO basis to the AO basis.
  !> @detail Computes the transformation Fao = S * V * Fmo * (S * V)^T,
  !>         where V are the MO coefficients and S is the overlap matrix,
  !>         typically used for converting the Fock matrix or similar
  !>         operators from MO to AO representation.
  !> @param[out] Fao Operator in AO basis (triangular format).
  !> @param[in] Fmo Operator in MO basis (triangular format).
  !> @param[in] smat_full Full overlap matrix in AO basis.
  !> @param[in] v MO coefficients.
  !> @param[in] nmo Number of molecular orbitals.
  !> @param[in] nbf Number of basis functions.
  !> @param[inout] sv Work array for S * V.
  !> @param[inout] work Work array for intermediate calculations.
  subroutine mo_to_ao(fao, fmo, smat_full, v, nmo, nbf, sv, work)
    use precision, only: dp
    use mathlib, only: pack_matrix, unpack_matrix
    use oqp_linalg

    implicit none

    real(kind=dp), intent(out) :: fao(:)
    real(kind=dp), intent(in) :: fmo(:)
    real(kind=dp), intent(in) :: smat_full(:,:)
    real(kind=dp), intent(in) :: v(*)
    real(kind=dp), intent(in) :: sv(*), work(*)
    integer, intent(in) :: nmo, nbf

    integer :: nbf2
    real(kind=dp), allocatable :: ftmp(:,:)

    allocate(ftmp(nbf,nbf))

    call unpack_matrix(fmo, ftmp)

    ! compute S*V
    call dsymm('l', 'u', nbf, nmo, &
               1.0_dp, smat_full, nbf, &
                       v,  nbf, &
               0.0_dp, sv, nbf)

    ! compute (S * V) * Fmo
    call dsymm('r', 'u', nbf, nmo, &
               1.0d0, ftmp, nbf, &
                      sv,  nbf, &
               0.0d0, work, nbf)

    ! compute ((S * V) * Fmo) * (S * V)^T
    call dgemm('n', 't', nbf, nbf, nmo, &
               1.0d0, work, nbf, &
                      sv, nbf, &
               0.0d0, ftmp, nbf)

    nbf2 = nbf*(nbf+1)/2
    call pack_matrix(ftmp, fao(:nbf2))

    deallocate(ftmp)

  end subroutine mo_to_ao

  !> @brief Computes the two-electron part (Coulomb and exchange) of the Fock matrix.
  !> @detail Forms the Coulomb (J) and exchange (K) contributions to the Fock matrix
  !>         using two-electron integrals,
  !>         with optional scaling of the exchange term for hybrid DFT methods.
  !> @param[in] basis Basis set information.
  !> @param[in] d Density matrices (triangular format).
  !> @param[inout] f Fock matrices to be updated (triangular format).
  !> @param[in] scalefactor Optional scaling factor for exchange (default = 1.0).
  !> @param[inout] infos System information.
  subroutine fock_jk(basis, d, f, scalefactor, infos)
    use precision, only: dp
    use io_constants, only: iw
    use util, only: measure_time
    use basis_tools, only: basis_set
    use types, only: information
    use int2_compute, only: int2_compute_t, int2_fock_data_t, &
                            int2_rhf_data_t, int2_urohf_data_t

    implicit none

    type(basis_set), intent(in) :: basis
    type(information), intent(inout) :: infos
    real(kind=dp), optional, intent(in) :: scalefactor

    integer :: i, ii, nf, nschwz
    real(kind=dp) :: scalef
    real(kind=dp), target, intent(in) :: d(:,:)
    real(kind=dp), intent(inout) :: f(:,:)

    type(int2_compute_t) :: int2_driver
    class(int2_fock_data_t), allocatable :: int2_data

    ! Initial Settings
    scalef = 1.0d0
    if (present(scalefactor)) scalef = scalefactor

    call measure_time(print_total=1, log_unit=IW)

    write(IW,"(/3x,'Form Two-Electron J and K Fock')")

    ! Initialize ERI calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()
    int2_data = int2_rhf_data_t(nfocks=1, d=d, scale_exchange=scalefactor)

    call flush(IW)

    ! Constructing two electron Fock matrix
    call int2_driver%run(int2_data)
    nschwz = int2_driver%skipped

    ! Scaling (everything except diagonal is halved)
    f =  0.5 * int2_data%f(:,:,1)

    do nf = 1, ubound(f,2)
      ii = 0
      do i = 1, basis%nbf
         ii = ii + i
         f(ii,nf) = 2*f(ii,nf)
      end do
    end do

    call int2_driver%clean()

  end subroutine fock_jk

end module scf
