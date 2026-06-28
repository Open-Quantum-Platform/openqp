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
    use dft, only: dft_initialize, dftclean, dft_setup_descent_grid
    use types, only: information
    use oqp_tagarray_driver
    use strings, only: Cstring, fstring
    use mod_dft_molgrid, only: dft_grid_t
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "hf_energy_mod"

    type(information), target, intent(inout) :: infos

    integer :: nbf2, nbf, nsh2
    logical :: urohf, dft
    type(basis_set), pointer :: basis
    type(dft_grid_t) :: molGrid
    ! Coarse->fine XC grid ramp: optional coarse "descent" grid (built by the
    ! unified policy in dft_setup_descent_grid; scf_driver ramps it to the
    ! production grid in the convergence tail).
    type(dft_grid_t) :: coarseGrid
    logical :: have_coarse

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

!   Coarse->fine XC grid ramp: build the optional coarse "descent" grid under the
!   single unified policy (default-on coarse-to-fine, or an opt-in progressive-
!   screening request). scf_driver picks coarse vs production per iteration and
!   pins to the production grid in the convergence tail.
    have_coarse = .false.
    if (dft) call dft_setup_descent_grid(infos, basis, molGrid, coarseGrid, have_coarse)

!   Run HF/DFT calculation
    if (have_coarse) then
      call scf_driver(basis, infos, molGrid, coarseGrid)
    else
      call scf_driver(basis, infos, molGrid)
    end if

!   Cleanup
    if (dft) call dftclean(infos)

    close(iw)

  end subroutine hf_energy

end module hf_energy_mod
