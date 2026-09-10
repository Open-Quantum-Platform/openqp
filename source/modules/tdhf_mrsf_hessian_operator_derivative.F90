module tdhf_mrsf_hessian_operator_derivative_mod

  use precision, only: dp
  use types, only: information
  use int2_compute, only: int2_compute_t
  use tdhf_mrsf_sigma_mod, only: prepare_mrsf_response_scaling
  use tdhf_mrsf_hessian_one_e_action_mod, only: &
    build_mrsf_one_e_operator_derivative_action
  use tdhf_mrsf_hessian_two_e_action_mod, only: &
    build_mrsf_two_e_operator_derivative_action

  implicit none

  private
  public :: build_mrsf_operator_derivative_action

contains

!###############################################################################

  subroutine build_mrsf_operator_derivative_action(infos,int2_driver,mo_a, &
      mo_b,dmo_a,dmo_b,fock_a_ao,fock_b_ao,dfock_a_ao,dfock_b_ao,x_packed, &
      dax,status)
    ! Total first nuclear derivative (dA/dR)X in the packed spin-adapted MRSF
    ! response space.  The one-electron Fock difference and all seven
    ! two-electron density channels are differentiated before the physical
    ! amplitude response is solved.

    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: int2_driver
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:),dmo_b(:,:,:)
    real(kind=dp), intent(in) :: fock_a_ao(:,:),fock_b_ao(:,:), &
      dfock_a_ao(:,:,:),dfock_b_ao(:,:,:),x_packed(:)
    real(kind=dp), intent(out) :: dax(:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: one_e(:,:),two_e(:,:)
    real(kind=dp) :: response_scale
    integer :: ncoord,packed,local_status

    packed=size(x_packed)
    ncoord=size(dmo_a,3)
    status=0
    dax=0.0_dp
    if(packed<=0 .or. ncoord<=0 .or. any(shape(dax)/=[packed,ncoord])) then
      status=-1
      return
    end if
    call prepare_mrsf_response_scaling(infos,response_scale,local_status)
    if(local_status/=0) then
      status=-2
      return
    end if
    allocate(one_e(packed,ncoord),two_e(packed,ncoord),source=0.0_dp)
    call build_mrsf_one_e_operator_derivative_action(infos,mo_a,mo_b,dmo_a, &
      dmo_b,fock_a_ao,fock_b_ao,dfock_a_ao,dfock_b_ao,x_packed,one_e, &
      local_status)
    if(local_status==0) call build_mrsf_two_e_operator_derivative_action( &
      infos,int2_driver,mo_a,mo_b,dmo_a,dmo_b,x_packed,response_scale, &
      infos%tddft%spc_coco,infos%tddft%spc_ovov,infos%tddft%spc_coov, &
      two_e,local_status)
    if(local_status==0) then
      dax=one_e+two_e
    else
      status=-3
    end if
    deallocate(one_e,two_e)
  end subroutine build_mrsf_operator_derivative_action

end module tdhf_mrsf_hessian_operator_derivative_mod
