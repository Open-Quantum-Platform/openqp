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
    call guess_huckel(inf)
  end subroutine guess_huckel_C

  subroutine guess_huckel(infos)
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

    implicit none

    character(len=*), parameter :: subroutine_name = "guess_huckel"

    type(information), target, intent(inout) :: infos
    integer :: i, nbf, nbf2

    type(basis_set), pointer :: basis
    type(basis_set) :: huckel_basis
    character(len=:), allocatable :: basis_file
    logical :: err

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
  !
  ! Section of Tagarray for the log filename
  ! We are getting lot file name from Python via tagarray
  !
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, with_abort)
    call tagarray_get_data(infos%dat, OQP_hbasis_filename, basis_filename)
    allocate(character(ubound(basis_filename,1)) :: basis_file)
    do i = 1, ubound(basis_filename,1)
       basis_file(i:i) = basis_filename(i)
    end do
   !
  ! Files open
  ! 1. XYZ: Read : Geometric data, ATOMS
  ! 3. LOG: Read Write: Main output file
  !
    open (unit=IW, file=infos%log_filename, position="append")
  !
    call print_module_info('Guess_Huckel','Initial guess using Huckel theory')
  !
  ! Readings
  ! load basis set
    basis => infos%basis
    call huckel_basis%from_file(basis_file, infos%atoms, err)
  ! Checking error of basis set reading..
    infos%control%basis_set_issue = err
    basis%atoms => infos%atoms
  !  Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nbf2 =nbf*(nbf+1)/2

    ! clean data
    call infos%dat%remove_records(tags_alpha)
    call infos%dat%remove_records(tags_beta)

    ! load general data
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM, smat)

    ! allocate alpha
    call infos%dat%reserve_data(OQP_DM_A, TA_TYPE_REAL64, nbf2, comment=OQP_DM_A_comment)
    call infos%dat%reserve_data(OQP_E_MO_A, TA_TYPE_REAL64, nbf, comment=OQP_E_MO_A_comment)
    call infos%dat%reserve_data(OQP_VEC_MO_A, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_VEC_MO_A_comment)

    ! load alpha data
    call data_has_tags(infos%dat, tags_alpha, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)

  ! UHF/ROHF
    if (infos%control%scftype >= 2) then
      ! allocate beta
      call infos%dat%reserve_data(OQP_DM_B, TA_TYPE_REAL64, nbf2, comment=OQP_DM_B_comment)
      call infos%dat%reserve_data(OQP_E_MO_B, TA_TYPE_REAL64, nbf, comment=OQP_E_MO_B_comment)
      call infos%dat%reserve_data(OQP_VEC_MO_B, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_VEC_MO_B_comment)

      ! load beta
      call data_has_tags(infos%dat, tags_beta, module_name, subroutine_name, WITH_ABORT)
      call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
      call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
      call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    end if

  ! Calculate Huckel MO
    call huckel_guess(Smat, MO_A, infos, basis, huckel_basis)
  !  For ROHF/UHF
    if (INFOS%control%scftype >= 2) MO_B = MO_A
  !
  ! Calculate Density Matrix
  ! RHF
    if (INFOS%control%scftype == 1) then
      call get_ab_initio_density(Dmat_A, MO_A, infos=infos, basis=basis)
  ! ROHF/UHF
    else
      call get_ab_initio_density(Dmat_A, MO_A, Dmat_B, MO_B, infos, basis)
    end if
  !
    write (iw, '(/x,a,/)') '...... End of initial orbital guess ......'
    call measure_time(print_total=1, log_unit=iw)
    close(iw)
  end subroutine guess_huckel

end module guess_huckel_mod
