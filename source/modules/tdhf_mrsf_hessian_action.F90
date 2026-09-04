module tdhf_mrsf_hessian_action_mod

  use precision, only: dp
  use types, only: information
  use tdhf_lib, only: iatogen
  use tdhf_mrsf_lib, only: mrsfcbc,mrsfmntoia
  use tdhf_mrsf_hessian_eri_derivative_mod, only: &
    build_mrsf_explicit_eri_derivative

  implicit none

  private
  public :: build_mrsf_explicit_eri_derivative_action

contains

!###############################################################################

  subroutine build_mrsf_explicit_eri_derivative_action(infos,mo_a,mo_b, &
      x_packed,response_scale,spc_coco,spc_ovov,spc_coov,dax,status)
    ! Build the derivative-ERI component of (dA/dR)X at fixed MOs and fixed
    ! physical MRSF amplitude.  This calls the same spin-adapted density and
    ! adjoint maps as the production sigma while the derivative integral
    ! contraction itself is supplied independently by fock_deriv.

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),x_packed(:)
    real(kind=dp), intent(in) :: response_scale,spc_coco,spc_ovov,spc_coov
    real(kind=dp), intent(out) :: dax(:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable, target :: density(:,:,:)
    real(kind=dp), allocatable :: derivative_fock(:,:,:,:),x_matrix(:,:)
    integer :: coordinate,nbf,nocca,noccb,ncoord,packed

    nbf=infos%basis%nbf
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    packed=nocca*(nbf-noccb)
    ncoord=3*size(infos%atoms%xyz,2)
    status=0
    dax=0.0_dp
    if(nocca-noccb/=2 .or. size(x_packed)/=packed .or. &
       any(shape(dax)/=[packed,ncoord]) .or. &
       any(shape(mo_a)/=[nbf,nbf]) .or. any(shape(mo_b)/=[nbf,nbf])) then
      status=-1
      return
    end if
    allocate(density(7,nbf,nbf),derivative_fock(7,nbf,nbf,ncoord), &
      x_matrix(nbf,nbf),source=0.0_dp)
    call iatogen(x_packed,x_matrix,nocca,noccb)
    call mrsfcbc(infos,mo_a,mo_b,x_matrix,density)
    call build_mrsf_explicit_eri_derivative(infos,density,response_scale, &
      spc_coco,spc_ovov,spc_coov,infos%tddft%mult,derivative_fock,status)
    if(status==0) then
      do coordinate=1,ncoord
        call mrsfmntoia(infos,derivative_fock(:,:,:,coordinate),dax, &
          mo_a,mo_b,coordinate)
      end do
    end if
    deallocate(density,derivative_fock,x_matrix)
  end subroutine build_mrsf_explicit_eri_derivative_action

end module tdhf_mrsf_hessian_action_mod
