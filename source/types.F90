!  17 Aug 12 - CHC - Initial file
module types
!  Definition of types

  use precision, only: dp
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t, c_double, c_char, c_bool, c_int
  use tagarray, only: container_t
  use functionals, only: functional_t
  use atomic_structure_m, only: atomic_structure
  use parallel, only: MPI_COMM_NULL
  use basis_tools, only: basis_set

  implicit none

  private

!     The information of a system
  type, public, bind(C) :: molecule
    integer(c_int64_t) :: natom      = 0 !< The number of atom
    integer(c_int64_t) :: charge     = 0 !< Molecular charge
    integer(c_int64_t) :: nelec      = 0 !< The number of electron
    integer(c_int64_t) :: nelec_A    = 0 !< The number of alpha electron
    integer(c_int64_t) :: nelec_B    = 0 !< The number of beta electron
    integer(c_int64_t) :: mult       = 0 !< Spin multiplicity
    integer(c_int64_t) :: nvelec     = 0 !< The number of valence electron
    integer(c_int64_t) :: nocc       = 0 !< The number of occupied orbitals
    !< nOCC = nelec/2 for RHF
    !< nOCC = nelec_A for ROHF/UHF with mult=3
    !< nOCC = nelec/2 for ROHF/UHF with mult=1
  end type molecule

  type, public, bind(C) :: dft_parameters
    character(kind=c_char) :: XC_functional_name(20)  !< Name of XC functional
    real(c_double) :: hfscale               = 1.0_dp  !< HF scale for global hybrids
    real(c_double) :: cam_alpha             = 0.0_dp  !< alpha coefficient
    real(c_double) :: cam_beta              = 0.0_dp  !< beta coefficient
    real(c_double) :: cam_mu                = 0.0_dp  !< mu coefficient
    real(c_double) :: MP2SS_Scale           = 0.0_dp  !< coefficient of non-local MP2 same-spin correlation
    real(c_double) :: MP2OS_Scale           = 0.0_dp  !< coefficient of non-local MP2 opposite-spin correlation
    logical(c_bool) :: cam_flag             = .false. !< switch Coulomb-Attenuating Method (CAM-) \frac{1}{r_{12}} = \frac{1-[cam_{alpha}+cam_{beta}*\erf(cam_{mu}*r_{12})]}{r_{12}} + \frac{cam_{alpha}+cam_{beta}*\erf(cam_{mu}*r_{12})}{r_{12}}
    logical(c_bool) :: dh_flag              = .false. !< logical flag of using DH-DFT functionals
    logical(c_bool) :: grid_pruned          = .false. !< true if pruned grid (e.g. sg1) is used
    logical(c_bool) :: grid_ao_pruned           = .true.  !< true if grid is pruned by AO distance
    real(c_double) :: grid_ao_threshold         = 0.0_dp  !< Prune grid AOs with threshold
    real(c_double) :: grid_ao_sparsity_ratio    = 0.9_dp  !< Prune grid AOs if sparsity exceeds 10%
    character(c_char) :: grid_pruned_name(16)   = ''  !< prune grid name
    integer(c_int64_t) :: grid_num_ang_grids    = 0      !< number of angular grids, >1 means pruned grid
    integer(c_int64_t) :: grid_rad_size         = 96     !< number of radial grid pts.
    integer(c_int64_t) :: grid_ang_size         = 302    !< number of angular grid pts.
    real(c_double)     :: grid_density_cutoff   = 0.0d0  !< grid DFT density cutoff
    integer(c_int64_t) :: dft_partfun           = 0  !< partition function type in grid-based DFT
                                                     !< -  0 (default) - SSF original polynomial
                                                     !< -  1           - Becke's 4th degree polynomial
                                                     !< -  2           - Modified SSF ( erf(x/(1-x**2)) )
                                                     !< -  3           - Modified SSF (smoothstep-2, 5th order)
                                                     !< -  4           - Modified SSF (smoothstep-3, 7rd order)
                                                     !< -  5           - Modified SSF (smoothstep-4, 9th order)
                                                     !< -  6           - Modified SSF (smoothstep-5,11th order)
                                                     !< Note: Becke's polynomial with 1 iteration is actually
                                                     !< a smoothstep-1 polynomial (3*x**2 - 2*x**3)
    integer(c_int64_t) :: rad_grid_type         = 0  !< type of the radial grid in DFT:
                                                     !< - 0 (default) - Euler-Maclaurin grid (Murray et al.)
                                                     !< - 1           - Log3 grid (Mura and Knowles)
                                                     !< - 2           - Treutler and Ahlrichs radial grid
                                                     !< - 3           - Becke's grid

    integer(c_int64_t) :: dft_bfc_algo          = 0  !< type of the Becke's fuzzy cell method
                                                     !< - 0 (default) - SSF-like algorithm
                                                     !< - 1           - Becke's algorithm
    logical(c_bool)    :: dft_wt_der            = .false.   !< .TRUE. if quadrature weights derivative
                                                            !<  contribution to the nuclear gradient is needed
                                                            !< Weight derivatives are not always required,
                                                            !< especially if the fine grid is used


  end type dft_parameters

  type, public, bind(C) :: energy_results
    real(c_double) :: energy         = 0.0_dp !< Total energy
    real(c_double) :: enuc           = 0.0_dp !< Nuclear repulsion energy
    real(c_double) :: psinrm         = 0.0_dp !< wavefunction normalization
    real(c_double) :: ehf1           = 0.0_dp !< one-electron energy
    real(c_double) :: vee            = 0.0_dp !< two-electron energy
    real(c_double) :: nenergy        = 0.0_dp !< nuclear repulsion energy
    real(c_double) :: etot           = 0.0_dp !< total energy
    real(c_double) :: vne            = 0.0_dp !< nucleus-electron potential energy
    real(c_double) :: vnn            = 0.0_dp !< nucleus-nucleus potential energy
    real(c_double) :: vtot           = 0.0_dp !< total potential energy
    real(c_double) :: tkin           = 0.0_dp !< total kinetic energy
    real(c_double) :: virial         = 0.0_dp !< virial ratio (v/t)
    real(c_double) :: excited_energy = 0.0_dp !< targeted excited state energy
    logical(c_bool) :: SCF_converged = .false. !< Convergence checking Flag for SCF
    logical(c_bool) :: Davidson_converged = .false. !< Convergence checking Flag for Davidson Iteration
    logical(c_bool) :: Z_Vector_converged = .false. !< Convergence checking Flag for Z-Vector Iteration
  end type energy_results

  type, public, bind(C) :: control_parameters
    integer(c_int64_t) :: hamilton = 10      !< The method of calculations: 10=HF, 20=DFT
    integer(c_int64_t) :: scftype  = 1       !< Refence wavefuction, 1= RHF 2= UHF 3= ROHF
    character(c_char)  :: runtype(20)  = ''  !<  Run type: energy, grad, etc.
    integer(c_int64_t) :: guess    = 1       !< used guess
    integer(c_int64_t) :: active_basis = 0    !< Choose data basis: 0 -> info%basis, 1 -> info%alt_basis
    integer(c_int64_t) :: maxit    = 3       !< The maximum number of iterations
    integer(c_int64_t) :: maxit_dav = 50     !< The maximum number of iterations in Davidson eigensolver
    integer(c_int64_t) :: maxit_zv = 50      !< The maximum number of CG iterations in Z-vector subroutines
    integer(c_int64_t) :: maxdiis  = 7               !< The maximum number of diis equations
    integer(c_int64_t) :: diis_reset_mod = 10        !< The maximum number of diis iteration before resetting
    real(c_double) :: diis_reset_conv = 0.005_dp     !< Convergency criteria of DIIS reset
    real(c_double) :: diis_method_threshold = 2.0_dp !< DIIS threshold for switching DIIS method
    integer(c_int64_t) :: diis_type = 5              !< 1: none, 2: cdiis, 3: ediis, 4: adiis, 5: vdiis
    real(c_double) :: vdiis_cdiis_switch = 0.3_dp    !< The threshold for selecting cdiis
    real(c_double) :: vdiis_vshift_switch = 0.003_dp !< The threshold for setting vshift = 0
    real(c_double) :: vshift_cdiis_switch = 0.3_dp   !< The threshold for selecting cdiis for vshift
    real(c_double) :: vshift = 0.0_dp                !< Virtual orbital shift for ROHF
    logical(c_bool) :: mom = .false.                 !< Maximum Overlap Method for SCF Convergency
    logical(c_bool) :: pfon = .false.                !< Pseudo-Fractional Occupation Number Method (pFON) for scf
    real(c_double) :: mom_switch = 0.003_dp          !< Turn on criteria of DIIS error
    real(c_double) :: pfon_start_temp = 2000.0_dp    !< Starting tempreature for pFON
    real(c_double) :: pfon_cooling_rate = 50.0_dp    !< Tempreature cooling rate for pFON
    real(c_double) :: pfon_nsmear = 5.0_dp           !< Num. of smearing orbitals for pFON if = 0, all 
    real(c_double) :: conv = 1e-6_dp                 !< Convergency criteria of SCF
    integer(c_int64_t) :: scf_incremental = 1        !< Enable/disable incremental Fock build
    real(c_double) :: int2e_cutoff = 5e-11_dp        !< 2e-integrals cutoff
    integer(c_int64_t) :: esp = 0                    !< (R)ESP charges, 0 - skip, 1 - ESP, 2 - RESP
    integer(c_int64_t) :: resp_target = 0            !< RESP charges target: 0 - zero, 1 - Mulliken
    real(c_double) :: resp_constr = 0.01             !< RESP charges constraint
    logical(c_bool) :: basis_set_issue = .false.     !< Basis set issue flag

    real(c_double) :: conf_print_threshold = 5.0d-02 !< The threshold for configuration printout
    logical(c_bool) :: rstctmo = .false.               !< Restrict new MO similar to previous MO. This is similar to MOM method
    ! SOSCF parameters
    integer(c_int64_t) :: soscf_type = 0       !< SOSCF type: 0=off, 1=SOSCF only, 2=SOSCF+DIIS
    real(c_double) :: soscf_lvl_shift = 0.0_dp !< Level shifting parameter for SOSCF
    integer(c_int64_t) :: soscf_reset_mod = 0  !< Reset the orbital Hessian. If it is zero, we don't reset by default.
    integer(c_int64_t) :: verbose = 1          !< Controls output verbosity: 0 for minimal, 1+ for detailed.
  end type control_parameters

  type, public, bind(c) :: tddft_parameters
    integer(c_int64_t) :: nstate = 1       !< Number of excited states
    integer(c_int64_t) :: target_state = 1 !< Target excited state for properties calculation, ground state == 0
    integer(c_int64_t) :: maxvec = 50      !< Max number of trial vectors
    integer(c_int64_t) :: mult = 1         !< MRSF multiplicity
    real(c_double) :: cnvtol = 1.0e-10_dp  !< convergence tolerance in the iterative TD-DFT step
    real(c_double) :: zvconv = 1.0e-10_dp  !< convergence tolerance in Z-vector equation
    logical(c_bool) :: debug_mode = .false.!< Debug print
    logical(c_bool) :: tda = .false.       !< switch for Tamm-Dancoff approximation
    integer(c_int64_t) :: tlf = 2          !< truncated Leibniz formula (TLF) approximation algorithm,
                                           !< 0   - zeroth-order (scales as O(n^2)) DO NOT WORK
                                           !< 1   - first-order (scales as O(n^3))
                                           !< 2   - second-order (scales as O(n^3))
    real(c_double) :: HFScale = 1.0_dp     !< HF scale for global hybrids in response calculations
    real(c_double) :: cam_alpha = 0.0_dp   !< alpha coefficient
    real(c_double) :: cam_beta = 0.0_dp    !< beta coefficient
    real(c_double) :: cam_mu = 0.0_dp      !< mu coefficient
    real(c_double) :: spc_coco = 0.0_dp    !< Spin-pair coupling parameter MRSF (C=closed, O=open, V=virtual MOs)
    real(c_double) :: spc_ovov = 0.0_dp    !< Spin-pair coupling parameter MRSF (C=closed, O=open, V=virtual MOs)
    real(c_double) :: spc_coov = 0.0_dp    !< Spin-pair coupling parameter MRSF (C=closed, O=open, V=virtual MOs)
    type(c_ptr) :: ixcore                  !< orbital index responsible for excitation (ixcore=1 means that it computes 
    integer(c_int64_t) :: ixcore_len = 1   !< length of ixcore
  end type tddft_parameters

  type, public, bind(c) :: mpi_communicator
    integer(c_int) :: comm = MPI_COMM_NULL       !< MPI communicator
    logical(c_bool) :: debug_mode = .false.
    logical(c_bool) :: usempi = .false.
  end type mpi_communicator

  type, public, bind(c) :: electron_shell
    integer(c_int) :: id = 0
    integer(c_int) :: element_id = -1
!    integer(c_int) :: num_expo = 0
    integer(c_int) :: ang_mom = 0
    integer(c_int) :: ecp_nam = 0
    type(c_ptr) :: num_expo
    type(c_ptr) :: expo
    type(c_ptr) :: coef
    type(c_ptr) :: ecp_am
    type(c_ptr) :: ecp_rex
    type(c_ptr) :: ecp_coord
    type(c_ptr) :: ecp_zn
  end type electron_shell

  type, public :: information
    type(molecule) :: mol_prop
    type(energy_results) :: mol_energy
    type(dft_parameters) :: dft
    type(control_parameters) :: control
    type(atomic_structure) :: atoms
    type(functional_t) :: functional
    type(tddft_parameters) :: tddft
    type(container_t) :: dat
    type(basis_set) :: basis
    type(basis_set) :: alt_basis
    character(len=:), allocatable :: log_filename
    type(mpi_communicator) :: mpiinfo
    type(electron_shell) :: elshell
  contains
    generic :: set_atoms => set_atoms_arr, set_atoms_atm
    procedure, pass :: set_atoms_arr => info_set_atoms_arr
    procedure, pass :: set_atoms_atm => info_set_atoms_atm
  end type information

contains

  function info_set_atoms_arr(this, natoms, x, y, z, q, mass) result(ok)
    class(information) :: this
    integer(c_int64_t) :: natoms
    real(c_double) :: x(*), y(*), z(*), q(*)
    real(c_double), optional :: mass(*)
    integer(c_int) :: ok

    integer :: i

    ok = this%atoms%init(natoms)
    if (ok/=0) return

    do i = 1, natoms
      this%atoms%xyz(1,i) = x(i)
      this%atoms%xyz(2,i) = y(i)
      this%atoms%xyz(3,i) = z(i)
      this%atoms%zn(i) = q(i)
      if (present(mass)) this%atoms%mass(i) = mass(i)
    end do
    this%mol_prop%natom = natoms
  end function

  function info_set_atoms_atm(this, atoms) result(ok)
    class(information) :: this
    type(atomic_structure) :: atoms
    integer(c_int) :: ok

    ok = 1
    this%atoms = atoms

  end function

end module types

