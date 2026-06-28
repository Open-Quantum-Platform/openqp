module hf_energy_mod

  implicit none

  character(len=*), parameter :: module_name = "hf_energy_mod"

  private

  public hf_energy

contains

  subroutine hf_energy_C(c_handle) bind(C, name="hf_energy")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call hf_energy(inf)
  end subroutine hf_energy_C

  subroutine hf_energy(infos)
    use io_constants, only: iw
    use basis_tools, only: basis_set
    use messages, only: show_message
    use scf, only: scf_driver
    use dft, only: dft_initialize, dftclean
    use types, only: information
    use oqp_tagarray_driver
    use strings, only: Cstring, fstring
    use mod_dft_molgrid, only: dft_grid_t
    use printing, only: print_module_info
    use iso_c_binding, only: c_char, c_null_char, c_bool

    implicit none

    character(len=*), parameter :: subroutine_name = "hf_energy_mod"

    type(information), target, intent(inout) :: infos

    integer :: nbf2, nbf, nsh2
    logical :: urohf, dft
    type(basis_set), pointer :: basis
    type(dft_grid_t) :: molGrid
    ! Progressive coarse->fine XC grid ramp (optional)
    type(dft_grid_t) :: coarseGrid
    logical :: ps_have_coarse, ps_on_env
    integer :: ps_grid_rad, ps_grid_ang, ps_el
    integer(8) :: ps_rad0, ps_ang0
    logical(c_bool) :: ps_pruned0
    character(kind=c_char) :: ps_xcname0(20)
    character(len=64) :: ps_ev

    urohf = infos%control%scftype == 2 .or. infos%control%scftype == 3
    dft = infos%control%hamilton == 20

!   3. LOG: Write: Main output file
    open (unit=iw, file=infos%log_filename, position="append")
!
    call print_module_info('HF_DFT_Energy','Computing HF/DFT SCF Energy')

!   Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

!   Allocate memory
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nsh2 = (basis%nshell**2+basis%nshell)/2

    ! clean data
    call infos%dat%remove_records((/ character(len=80) :: OQP_FOCK_A, OQP_FOCK_B /))

    call infos%dat%reserve_data(OQP_FOCK_A, TA_TYPE_REAL64, nbf2, comment=OQP_FOCK_A_comment)
    call check_status(infos%dat%get_status(), module_name, subroutine_name, OQP_FOCK_A)

    if (urohf) then
      call infos%dat%reserve_data(OQP_FOCK_B, TA_TYPE_REAL64, nbf2, comment=OQP_FOCK_B_comment)
      call check_status(infos%dat%get_status(), module_name, subroutine_name, OQP_FOCK_B)
    end if

!   Prepare dft grid
    if (dft) call dft_initialize(infos, basis, molGrid, verbose=.true.)

!   Optional progressive coarse grid for the SCF descent (gated by scf_pscreen +
!   pscreen_grid_rad/ang, env-overridable). Built here because dft_initialize needs
!   basis intent(inout); the functional name is blanked so libxc is NOT re-initialised
!   (need_functional=.false. plus an empty name). scf_driver picks coarse vs full grid
!   per iteration and pins to the full grid in the convergence tail.
    ps_grid_rad = int(infos%control%pscreen_grid_rad)
    ps_grid_ang = int(infos%control%pscreen_grid_ang)
    call get_environment_variable("OQP_PSCREEN_GRID_RAD", ps_ev, ps_el)
    if (ps_el > 0) read(ps_ev,*,iostat=ps_el) ps_grid_rad
    call get_environment_variable("OQP_PSCREEN_GRID_ANG", ps_ev, ps_el)
    if (ps_el > 0) read(ps_ev,*,iostat=ps_el) ps_grid_ang
    ps_on_env = infos%control%scf_pscreen /= 0
    call get_environment_variable("OQP_PSCREEN", ps_ev, ps_el)
    if (ps_el > 0) ps_on_env = (ps_ev(1:1)=='1' .or. ps_ev(1:1)=='t' .or. ps_ev(1:1)=='T' &
                                .or. ps_ev(1:1)=='y' .or. ps_ev(1:1)=='Y')
    ps_have_coarse = dft .and. ps_on_env .and. ps_grid_rad > 0 .and. ps_grid_ang > 0
    if (ps_have_coarse) then
      ps_rad0 = infos%dft%grid_rad_size
      ps_ang0 = infos%dft%grid_ang_size
      ps_pruned0 = infos%dft%grid_pruned
      ps_xcname0 = infos%dft%XC_functional_name
      infos%dft%grid_rad_size = ps_grid_rad
      infos%dft%grid_ang_size = ps_grid_ang
      infos%dft%grid_pruned   = .false.
      infos%dft%XC_functional_name = c_null_char
      call dft_initialize(infos, basis, coarseGrid, verbose=.false., need_functional=.false.)
      infos%dft%grid_rad_size = ps_rad0
      infos%dft%grid_ang_size = ps_ang0
      infos%dft%grid_pruned   = ps_pruned0
      infos%dft%XC_functional_name = ps_xcname0
      write(iw,'(5x,a,i0,a,i0)') &
        'Progressive XC coarse grid (descent): NRAD=', ps_grid_rad, '  NLEB=', ps_grid_ang
    end if

!   Run HF/DFT calculation
    if (ps_have_coarse) then
      call scf_driver(basis, infos, molGrid, coarseGrid)
    else
      call scf_driver(basis, infos, molGrid)
    end if

!   Cleanup
    if (dft) call dftclean(infos)

    close(iw)

  end subroutine hf_energy

end module hf_energy_mod
