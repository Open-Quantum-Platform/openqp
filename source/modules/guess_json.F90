module guess_json_mod

  implicit none

  character(len=*), parameter :: module_name = "guess_json_mod"

contains

  subroutine guess_json_C(c_handle) bind(C, name="guess_json")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    print *, "HELLLO"
    call guess_json(inf)
  end subroutine guess_json_C

  subroutine guess_json(infos)
    use precision, only: dp
    use types, only: information
    use io_constants, only: IW
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use guess, only: get_ab_initio_density
    use mathlib, only: matrix_invsqrt
    use huckel, only: huckel_guess
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use strings, only: Cstring, fstring
    use printing, only: print_module_info
    use oqp_tagarray_driver
    use iso_c_binding, only: c_char
    use parallel, only: par_env_t

    implicit none

    character(len=*), parameter :: subroutine_name = "guess_json"

    type(information), target, intent(inout) :: infos
    integer :: i, nbf, nbf2

    type(basis_set), pointer :: basis
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
    character(len=*), parameter :: tags_general(1) = (/ character(len=80) :: &
      OQP_SM /)

  ! Files open
  ! 1. XYZ: Read : Geometric data, ATOMS
  ! 3. LOG: Read Write: Main output file
  !
    open (unit=IW, file=infos%log_filename, position="append")

    print_module_info("Loading JSON", "Using stored SCF guess")

  ! load basis set
    basis => infos%basis
    call pe%init(infos%mpiinfo%comm, infos%mpiinfo%usempi)

    basis%atoms => infos%atoms
  !  Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nbf2 =nbf*(nbf+1)/2

    ! load general data
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM, smat)

    ! load alpha data
    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

    ! allocate beta
    call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
    call infos%dat%reserve_data(OQP_DM_B, TA_TYPE_REAL64, nbf2, comment=OQP_DM_B_comment)
    call infos%dat%reserve_data(OQP_E_MO_B, TA_TYPE_REAL64, nbf, comment=OQP_E_MO_B_comment)
    call infos%dat%reserve_data(OQP_VEC_MO_B, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_VEC_MO_B_comment)

   ! load beta
    call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)


  !  For ROHF/UHF
    if (INFOS%control%scftype >= 2) MO_B = MO_A

! Calculate Density Matrix
    if (pe%rank == root) then
    ! RHF
      if (infos%control%scftype == 1) then
        call get_ab_initio_density(Dmat_A, MO_A, infos=infos, basis=basis)
    ! ROHF/UHF
      else
        call get_ab_initio_density(Dmat_A, MO_A, Dmat_B, MO_B, infos, basis)
      endif
    endif
    ! Broadcast MO and density matrices to all processes
    call pe%bcast(MO_A, nbf*nbf)
    if (infos%control%scftype >= 2) then
      call pe%bcast(MO_B, nbf*nbf)
    endif
    ! Broadcast the density matrices to all processes
    if (infos%control%scftype == 1) then
      call pe%bcast(Dmat_A, nbf2)
    else
      call pe%bcast(Dmat_A, nbf2)
      call pe%bcast(Dmat_B, nbf2)
    endif
    call pe%barrier()
    write (iw, '(/x,a,/)') '...... End of initial orbital guess ......'
    call measure_time(print_total=1, log_unit=iw)
    close(iw)
  end subroutine guess_json

end module guess_json_mod
