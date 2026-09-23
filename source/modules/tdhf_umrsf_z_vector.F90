!> UMRSF-TDDFT Z-vector / response preparation.
!>
!> The Python gradient pipeline calls the TD z-vector entry point before the
!> gradient kernel.  Keep the expensive spin-coupled alpha/beta UMRSF response
!> solve here and cache the resulting response-gradient contribution for
!> tdhf_umrsf_gradient.
module tdhf_umrsf_z_vector_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_umrsf_z_vector_mod"

contains

  subroutine tdhf_umrsf_z_vector_C(c_handle) bind(C, name="tdhf_umrsf_z_vector")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    use precision, only: dp
    use io_constants, only: iw
    use messages, only: show_message, WITH_ABORT
    use oqp_tagarray_driver, only: OQP_umrsf_response_gradient, &
      OQP_umrsf_response_gradient_comment, TA_TYPE_REAL64, tagarray_get_data
    use tdhf_umrsf_gradient_mod, only: tdhf_umrsf_build_response_gradient
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    real(kind=dp), allocatable :: response_grad(:,:)
    real(kind=dp), contiguous, pointer :: cached_response(:,:)
    integer :: natom
    character(len=*), parameter :: tags_response(1) = (/ character(len=80) :: &
      OQP_umrsf_response_gradient /)

    inf => oqp_handle_get_info(c_handle)
    inf%mol_energy%Z_Vector_converged = .false.
    natom = ubound(inf%atoms%zn, 1)

    open(unit=iw, file=inf%log_filename, position="append")
    write(iw,'(/2x,a)') 'UMRSF Z-vector: solving coupled alpha/beta response'
    close(iw)

    call tdhf_umrsf_build_response_gradient(inf, response_grad)
    if (.not. allocated(response_grad)) then
      call show_message('UMRSF z-vector did not produce a response gradient.', WITH_ABORT)
      return
    end if

    call inf%dat%erase(tags_response)
    call tagarray_reserve_data(inf%dat, OQP_umrsf_response_gradient, TA_TYPE_REAL64, 3*natom, &
                               (/ 3, natom /), comment=OQP_umrsf_response_gradient_comment)
    call tagarray_get_data(inf%dat, OQP_umrsf_response_gradient, cached_response)
    cached_response = response_grad

    open(unit=iw, file=inf%log_filename, position="append")
    write(iw,'(2x,a)') 'UMRSF Z-vector: coupled alpha/beta response-gradient cache ready'
    close(iw)
    inf%mol_energy%Z_Vector_converged = .true.
  end subroutine tdhf_umrsf_z_vector_C

end module tdhf_umrsf_z_vector_mod
