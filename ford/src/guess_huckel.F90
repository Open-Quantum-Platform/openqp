!> @brief Extended Huckel initial guess drivers
!
!> @details The standard (`guess_huckel`) and modified (`guess_modhuckel`)
!>          extended Huckel guesses share the same driver; they differ only
!>          in the off-diagonal Wolfsberg-Helmholz formula, selected by the
!>          `modified` flag.
module guess_huckel_mod

  implicit none

  character(len=*), parameter :: module_name = "guess_huckel_mod"

contains

  subroutine guess_huckel_C(c_handle) bind(C, name="guess_huckel")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call guess_huckel_driver(inf, modified=.false.)
  end subroutine guess_huckel_C

  subroutine guess_modhuckel_C(c_handle) bind(C, name="guess_modhuckel")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call guess_huckel_driver(inf, modified=.true.)
  end subroutine guess_modhuckel_C

  subroutine guess_huckel_driver(infos, modified)
    use precision, only: dp
    use types, only: information
    use io_constants, only: IW
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use guess, only: get_ab_initio_density
    use huckel, only: huckel_guess
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use printing, only: print_module_info
    use iso_c_binding, only: c_char
    use parallel, only: par_env_t

    implicit none

    character(len=*), parameter :: subroutine_name = "guess_huckel_driver"

    type(information), target, intent(inout) :: infos
    logical, intent(in) :: modified

    integer :: i, nbf, nbf2

    type(basis_set), pointer :: basis
    type(basis_set) :: huckel_basis
    character(len=:), allocatable :: basis_file
    logical :: err
    integer , parameter :: root = 0
    type(par_env_t) :: pe
  ! tagarray
    real(kind=dp), contiguous, pointer :: &
      Smat(:), &
      dmat_a(:), mo_a(:,:), mo_energy_a(:), &
      dmat_b(:), mo_b(:,:), mo_energy_b(:)
    character(len=*), parameter :: tags_alpha(3) = (/ character(len=80) :: &
      OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_beta(3) = (/ character(len=80) :: &
      OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B /)
    character(len=*), parameter :: tags_general(2) = (/ character(len=80) :: &
      OQP_SM, OQP_hbasis_filename /)
    character(len=1,kind=c_char), contiguous, pointer :: basis_filename(:)

  ! The Huckel (MINI) basis set file name is passed from Python via tagarray
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_hbasis_filename, basis_filename)
    allocate(character(ubound(basis_filename,1)) :: basis_file)
    do i = 1, ubound(basis_filename,1)
       basis_file(i:i) = basis_filename(i)
    end do

    open (unit=IW, file=infos%log_filename, position="append")

    if (modified) then
      call print_module_info('Guess_ModHuckel','Initial guess using modified Huckel theory')
    else
      call print_module_info('Guess_Huckel','Initial guess using Huckel theory')
    end if

  ! load the Huckel minimal basis set on the master process
    basis => infos%basis
    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    err = .false.
    if (pe%rank == root) then
      call huckel_basis%from_file(basis_file, infos%atoms, err)
    end if

  ! Checking error of basis set reading..
    infos%control%basis_set_issue = err
    call pe%bcast(infos%control%basis_set_issue, 1)
    if (infos%control%basis_set_issue) then
      call show_message('Failed to read the Huckel (MINI) basis set from '//basis_file, WITH_ABORT)
    end if

    basis%atoms => infos%atoms

    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    ! load general data
    call tagarray_get_data(infos%dat, OQP_SM, smat)

    ! allocate alpha
    call infos%dat%alloc_or_die(OQP_DM_A, (/ nbf2 /), dmat_a, description=OQP_DM_A_comment)
    call infos%dat%alloc_or_die(OQP_E_MO_A, (/ nbf /), mo_energy_a, description=OQP_E_MO_A_comment)
    call infos%dat%alloc_or_die(OQP_VEC_MO_A, (/ nbf, nbf /), mo_a, description=OQP_VEC_MO_A_comment)

  ! UHF/ROHF
    if (infos%control%scftype >= 2) then
      ! allocate beta
      call infos%dat%alloc_or_die(OQP_DM_B, (/ nbf2 /), dmat_b, description=OQP_DM_B_comment)
      call infos%dat%alloc_or_die(OQP_E_MO_B, (/ nbf /), mo_energy_b, description=OQP_E_MO_B_comment)
      call infos%dat%alloc_or_die(OQP_VEC_MO_B, (/ nbf, nbf /), mo_b, description=OQP_VEC_MO_B_comment)
    end if

  ! Calculate Huckel MOs and the corresponding density matrix.
  ! All heavy work runs on the master process only; the results are
  ! broadcast below. mo_energy_a receives the Huckel eigenvalues for
  ! the projected orbitals (approximate guess orbital energies).
    if (pe%rank == root) then
      call huckel_guess(Smat, MO_A, infos, basis, huckel_basis, &
                        modified=modified, mo_energy=mo_energy_a)

  !   For ROHF/UHF the beta orbitals start identical to alpha
      if (infos%control%scftype >= 2) then
        MO_B = MO_A
        mo_energy_b = mo_energy_a
      end if

      if (infos%control%scftype == 1) then
        call get_ab_initio_density(Dmat_A, MO_A, infos=infos, basis=basis)
      else
        call get_ab_initio_density(Dmat_A, MO_A, Dmat_B, MO_B, infos, basis)
      end if
    end if

    ! Broadcast MO, MO energy and density matrices to all processes
    call pe%bcast(MO_A, nbf*nbf)
    call pe%bcast(mo_energy_a, nbf)
    call pe%bcast(Dmat_A, nbf2)
    if (infos%control%scftype >= 2) then
      call pe%bcast(MO_B, nbf*nbf)
      call pe%bcast(mo_energy_b, nbf)
      call pe%bcast(Dmat_B, nbf2)
    end if
    call pe%barrier()

    write (iw, '(/x,a,/)') '...... End of initial orbital guess ......'
    call measure_time(print_total=1, log_unit=iw)
    close(iw)
  end subroutine guess_huckel_driver

end module guess_huckel_mod
