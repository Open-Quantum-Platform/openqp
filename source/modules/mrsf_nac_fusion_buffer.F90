!> Small dependency-neutral exchange buffer for the experimental fused
!> MRSF gradient/NAC adjoint solve.  The legacy gradient module cannot USE the
!> NAC driver directly because the driver consumes kernels from that module.
module mrsf_nac_fusion_buffer_mod

  use precision, only: dp

  implicit none
  private
  public :: mrsf_nac_fusion_set_rhs, mrsf_nac_fusion_get_rhs, &
            mrsf_nac_fusion_set_solution, mrsf_nac_fusion_take_solution

  real(kind=dp), allocatable, save :: gradient_rhs_buffer(:)
  real(kind=dp), allocatable, save :: gradient_solution_buffer(:)

contains

  subroutine mrsf_nac_fusion_set_rhs(rhs)
    real(kind=dp), intent(in) :: rhs(:)

    if (allocated(gradient_rhs_buffer)) deallocate(gradient_rhs_buffer)
    allocate(gradient_rhs_buffer(size(rhs)), source=rhs)
    if (allocated(gradient_solution_buffer)) &
      deallocate(gradient_solution_buffer)
  end subroutine mrsf_nac_fusion_set_rhs

  subroutine mrsf_nac_fusion_get_rhs(rhs)
    use messages, only: show_message, WITH_ABORT
    real(kind=dp), intent(out) :: rhs(:)

    if (.not. allocated(gradient_rhs_buffer)) &
      call show_message('Fused gradient/NAC RHS buffer is empty.', WITH_ABORT)
    if (size(rhs) /= size(gradient_rhs_buffer)) &
      call show_message('Fused gradient/NAC RHS buffer has the wrong size.', &
                        WITH_ABORT)
    rhs = gradient_rhs_buffer
  end subroutine mrsf_nac_fusion_get_rhs

  subroutine mrsf_nac_fusion_set_solution(solution)
    real(kind=dp), intent(in) :: solution(:)

    if (allocated(gradient_solution_buffer)) &
      deallocate(gradient_solution_buffer)
    allocate(gradient_solution_buffer(size(solution)), source=solution)
  end subroutine mrsf_nac_fusion_set_solution

  subroutine mrsf_nac_fusion_take_solution(solution)
    use messages, only: show_message, WITH_ABORT
    real(kind=dp), intent(out) :: solution(:)

    if (.not. allocated(gradient_solution_buffer)) &
      call show_message('Fused gradient/NAC solution buffer is empty.', &
                        WITH_ABORT)
    if (size(solution) /= size(gradient_solution_buffer)) &
      call show_message('Fused gradient/NAC solution buffer has the wrong size.', &
                        WITH_ABORT)
    solution = gradient_solution_buffer
    deallocate(gradient_solution_buffer)
    if (allocated(gradient_rhs_buffer)) deallocate(gradient_rhs_buffer)
  end subroutine mrsf_nac_fusion_take_solution

end module mrsf_nac_fusion_buffer_mod
