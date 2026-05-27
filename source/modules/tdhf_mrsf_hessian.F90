module tdhf_mrsf_hessian_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_hessian_mod"

contains

!###############################################################################

  subroutine tdhf_mrsf_hessian_C(c_handle) bind(C, name="tdhf_mrsf_hessian")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_hessian(inf)
  end subroutine tdhf_mrsf_hessian_C

!###############################################################################

  subroutine tdhf_mrsf_hessian(infos)
    use types, only: information
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos

    ! Analytic MRSF-TDDFT Hessian kernel scaffold reached. The C ABI is present
    ! only for build/link integration. Runtime support must remain guarded until
    ! the MRSF gradient/Z-vector finite-difference baseline is validated across
    ! multiple molecules and roots; do not return placeholder zeros or fall back
    ! to a numerical Hessian from this native path.
    call show_message(&
      'Analytic MRSF-TDDFT Hessian kernel scaffold reached; implementation ' // &
      'awaits validated MRSF gradient/Z-vector finite-difference baseline.', &
      WITH_ABORT)
  end subroutine tdhf_mrsf_hessian

end module tdhf_mrsf_hessian_mod
