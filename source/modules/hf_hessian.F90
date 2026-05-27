module hf_hessian_mod

  implicit none

  character(len=*), parameter :: module_name = "hf_hessian_mod"

contains

!###############################################################################

  subroutine hf_hessian_C(c_handle) bind(C, name="hf_hessian")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information

    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf

    inf => oqp_handle_get_info(c_handle)
    call hf_hessian(inf)
  end subroutine hf_hessian_C

!###############################################################################

  subroutine hf_hessian(infos)
    use types, only: information
    use messages, only: show_message, WITH_ABORT

    implicit none

    type(information), target, intent(inout) :: infos

    ! Analytic HF/DFT Hessian kernel scaffold reached.  The C ABI is now
    ! intentionally present for build/link integration, but the scientific
    ! kernel is still guarded so `[hess] type=analytical` cannot silently fall
    ! back to the numerical Hessian or return placeholder zeros.
    call show_message(&
      'Analytic HF/DFT Hessian kernel scaffold reached; implementation is not available yet.', &
      WITH_ABORT)
  end subroutine hf_hessian

end module hf_hessian_mod
