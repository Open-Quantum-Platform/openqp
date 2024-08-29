module guess_hcore_mod

  implicit none

  character(len=*), parameter :: module_name = "guess_hcore_mod"

contains

  subroutine guess_hcore_C(c_handle) bind(C, name="guess_hcore")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call guess_hcore(inf)
  end subroutine guess_hcore_C

  subroutine guess_hcore(infos)
    use precision, only: dp
    use types, only: information
    use io_constants, only: IW
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use guess, only: get_ab_initio_density, &
                     get_ab_initio_orbital
    use mathlib, only: matrix_invsqrt
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use strings, only: Cstring, fstring
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "guess_hcore"

    type(information), target, intent(inout) :: infos
  !
    integer :: nbf2, nbf, ok
  !
    real(kind=dp), allocatable :: qmat(:,:)
    type(basis_set), pointer :: basis
  !
  ! tagarray
    real(kind=dp), contiguous, pointer :: &
      Hcore(:), Smat(:), &
      dmat_a(:), mo_a(:,:), mo_energy_a(:), &
      dmat_b(:), mo_b(:,:), mo_energy_b(:)
    character(len=*), parameter :: tags_alpha(3) = (/ character(len=80) :: &
      OQP_DM_A, OQP_E_MO_A, OQP_VEC_MO_A /)
    character(len=*), parameter :: tags_beta(3) = (/ character(len=80) :: &
      OQP_DM_B, OQP_E_MO_B, OQP_VEC_MO_B /)
    character(len=*), parameter :: tags_general(2) = (/ character(len=80) :: &
      OQP_SM, OQP_Hcore /)

  ! Files open
  ! 1. XYZ: Read : Geometric data, ATOMS
  ! 3. LOG: Read Write: Main output file
  !
    open (unit=IW, file=infos%log_filename, position="append")
  !
  !
    call print_module_info('guess_Hcore','Initial Guess using H Matrix')

  ! load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

  !  Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2

    allocate (qmat(nbf,nbf), stat=ok)
    if (ok /= 0) call show_message('Cannot allocate memory', WITH_ABORT)

    ! clean data
    call infos%dat%remove_records(tags_alpha)
    call infos%dat%remove_records(tags_beta)

    ! load general data
    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_SM, smat)
    call tagarray_get_data(infos%dat, OQP_Hcore, hcore)

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

  !  End of Readings.................................
  !
    call matrix_invsqrt(smat, qmat, nbf)
  !
  !  Calculate Hcore MO
    call Get_ab_initio_orbital(Hcore, MO_A, MO_Energy_A, QMat)
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
    write (IW, 9090)
    call measure_time(print_total=1, log_unit=iw)
    close(IW)
  9090 format(/1x, '...... End Of Initial Orbital Guess ......'/)

  end subroutine guess_hcore

end module guess_hcore_mod
