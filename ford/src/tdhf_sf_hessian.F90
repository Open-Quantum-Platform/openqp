module tdhf_sf_hessian_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_sf_hessian_mod"

contains

!###############################################################################

  subroutine tdhf_sf_hessian_C(c_handle) bind(C, name="tdhf_sf_hessian")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call tdhf_sf_hessian(inf)
  end subroutine tdhf_sf_hessian_C

!###############################################################################

  subroutine tdhf_sf_hessian(infos)
    use types, only: information
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos

    ! Analytic SF-TDDFT Hessian kernel scaffold reached. The C ABI is present
    ! for build/link integration, but the scientific kernel is deliberately
    ! guarded so `[hess] type=analytical` cannot return placeholder zeros or
    ! fall back to the numerical Hessian.
    call show_message(&
      'Analytic SF-TDDFT Hessian kernel scaffold reached; implementation is not available yet.', &
      WITH_ABORT)
  end subroutine tdhf_sf_hessian

end module tdhf_sf_hessian_mod
