module guess_sap_mod

  implicit none

  character(len=*), parameter :: module_name = "guess_sap_mod"

contains

  subroutine guess_sap_C(c_handle) bind(C, name="guess_sap")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call guess_sap(inf)
  end subroutine guess_sap_C

  subroutine guess_sap(infos)
    use precision, only: dp
    use types, only: information
    use io_constants, only: IW
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use guess, only: get_ab_initio_density, get_ab_initio_orbital
    use qmat_cache, only: get_qmat_cached
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use printing, only: print_module_info
    use parallel, only: par_env_t
    use iso_c_binding, only: c_char, c_null_char
    use mod_dft_molgrid, only: dft_grid_t
    use dft, only: dft_initialize
    use sap_lut, only: sap_table_t
    use mod_dft_gridint_sap, only: sap_potential_matrix

    implicit none

    character(len=*), parameter :: subroutine_name = "guess_sap"

    type(information), target, intent(inout) :: infos
    integer :: nbf, nbf2, ok, i, maxz

    real(kind=dp), allocatable :: qmat(:,:), fock(:), vsap(:)
    type(basis_set), pointer :: basis
    type(dft_grid_t) :: molGrid
    character(kind=c_char) :: saved_xcname(20)
    type(sap_table_t) :: sap
    character(len=:), allocatable :: sap_file
    logical :: err
    integer, parameter :: root = 0
    type(par_env_t) :: pe

  ! tagarray
    real(kind=dp), contiguous, pointer :: &
      Tmat(:), Smat(:), &
      dmat_a(:), mo_a(:,:), mo_energy_a(:), &
      dmat_b(:), mo_b(:,:), mo_energy_b(:)
    character(len=*), parameter :: tags_alpha(3) = (/ character(len=80) :: &
      OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_beta(3) = (/ character(len=80) :: &
      OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B /)
    character(len=*), parameter :: tags_general(3) = (/ character(len=80) :: &
      OQP_SM, OQP_TM, OQP_hbasis_filename /)
    character(len=1,kind=c_char), contiguous, pointer :: sap_filename(:)

    open (unit=IW, file=infos%log_filename, position="append")
    call print_module_info('Guess_SAP', &
        'Initial guess using Superposition of Atomic Potentials')

    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

  ! Resolve the SAP data file path passed from Python (via the hbasis tag)
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_hbasis_filename, sap_filename)
    allocate(character(ubound(sap_filename,1)) :: sap_file)
    do i = 1, ubound(sap_filename,1)
      sap_file(i:i) = sap_filename(i)
    end do

  ! Load the radial SAP table
    call sap%load(sap_file, err)
    if (err) call show_message('Guess_SAP: cannot read SAP data file '//trim(sap_file), WITH_ABORT)

  ! Warn about elements beyond the tabulated range
    maxz = maxval(nint(infos%atoms%zn))
    if (maxz > sap%zmax) then
      write(IW, '(1x,a,i0,a,i0,a)') &
        'Guess_SAP warning: element Z=', maxz, &
        ' exceeds tabulated SAP range (Z<=', sap%zmax, &
        '); those atoms contribute only the core Hamiltonian'
    end if

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    allocate(qmat(nbf,nbf), fock(nbf2), vsap(nbf2), stat=ok)
    if (ok /= 0) call show_message('Guess_SAP: cannot allocate memory', WITH_ABORT)

    ! clean previous data
    call infos%dat%remove_records(tags_alpha)
    call infos%dat%remove_records(tags_beta)

    ! load general data (overlap and kinetic-energy matrices are built by int1e).
    ! SAP uses the kinetic matrix T (not Hcore): the screened potential
    ! V_SAP = -Z_eff(r)/r already includes the nuclear attraction, so the guess
    ! Fock is F = T + V_SAP. Using Hcore would double-count nuclear attraction.
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call tagarray_get_data(infos%dat, OQP_TM, tmat)

    ! allocate alpha
    call infos%dat%reserve_data(OQP_DM_A, TA_TYPE_REAL64, nbf2, comment=OQP_DM_A_comment)
    call infos%dat%reserve_data(OQP_E_MO_A, TA_TYPE_REAL64, nbf, comment=OQP_E_MO_A_comment)
    call infos%dat%reserve_data(OQP_VEC_MO_A, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_VEC_MO_A_comment)
    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    ! UHF/ROHF
    if (infos%control%scftype >= 2) then
      call infos%dat%reserve_data(OQP_DM_B, TA_TYPE_REAL64, nbf2, comment=OQP_DM_B_comment)
      call infos%dat%reserve_data(OQP_E_MO_B, TA_TYPE_REAL64, nbf, comment=OQP_E_MO_B_comment)
      call infos%dat%reserve_data(OQP_VEC_MO_B, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_VEC_MO_B_comment)
      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
      call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    end if

  ! Build the DFT molecular grid. The XC functional name is not yet set at
  ! guess time (it is initialized later for the SCF), so temporarily blank it
  ! and request grid setup without a functional (need_functional=.false.) so
  ! dft_set_options skips any LibXC functional setup.
    saved_xcname = infos%dft%xc_functional_name
    infos%dft%xc_functional_name = c_null_char
    call dft_initialize(infos, basis, molGrid, verbose=.false., need_functional=.false.)
    infos%dft%xc_functional_name = saved_xcname

  ! Integrate the SAP potential matrix on the grid
    vsap = 0.0_dp
    call sap_potential_matrix(basis, molGrid, vsap, nbf, sap, infos)

  ! Guess Fock = T + V_SAP  (V_SAP already contains the screened nuclear term)
    fock(1:nbf2) = tmat(1:nbf2) + vsap(1:nbf2)

  ! Solve F C = eps S C
    call get_qmat_cached(infos, smat, qmat, nbf)
    call get_ab_initio_orbital(fock, MO_A, MO_Energy_A, QMat)
    if (infos%control%scftype >= 2) MO_B = MO_A

  ! Density matrix
    if (infos%control%scftype == 1) then
      call get_ab_initio_density(Dmat_A, MO_A, infos=infos, basis=basis)
    else
      call get_ab_initio_density(Dmat_A, MO_A, Dmat_B, MO_B, infos, basis)
    end if

    call sap%clean()

    write (IW, '(/1x,a/)') '...... End Of Initial Orbital Guess ......'
    call measure_time(print_total=1, log_unit=iw)
    close(IW)

  end subroutine guess_sap

end module guess_sap_mod
