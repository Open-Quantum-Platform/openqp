module tdhf_mrsf_hessian_two_e_action_mod

  use precision, only: dp
  use types, only: information
  use int2_compute, only: int2_compute_t
  use tdhf_lib, only: iatogen
  use tdhf_mrsf_lib, only: mrsfcbc,mrsfmntoia
  use tdhf_mrsf_density_contraction_mod, only: &
    contract_mrsf_seven_density_batch
  use tdhf_mrsf_hessian_eri_derivative_mod, only: &
    build_mrsf_explicit_eri_derivative
  use tdhf_mrsf_hessian_orbital_maps_mod, only: &
    build_mrsf_density_orbital_derivative, &
    build_mrsf_adjoint_orbital_derivative

  implicit none

  private
  public :: build_mrsf_two_e_operator_derivative_action

contains

!###############################################################################

  subroutine build_mrsf_two_e_operator_derivative_action(infos,int2_driver, &
      mo_a,mo_b,dmo_a,dmo_b,x_packed,response_scale,spc_coco,spc_ovov, &
      spc_coov,dax,status)
    ! Complete two-electron part of (dA/dR)X at fixed physical amplitude:
    !   (i) derivative ERIs acting on the unperturbed seven densities,
    !  (ii) ordinary ERIs acting on their orbital derivatives, and
    ! (iii) derivative of the adjoint MO backprojection.

    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: int2_driver
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:),dmo_b(:,:,:)
    real(kind=dp), intent(in) :: x_packed(:)
    real(kind=dp), intent(in) :: response_scale,spc_coco,spc_ovov,spc_coov
    real(kind=dp), intent(out) :: dax(:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable, target :: density(:,:,:),density_batch(:,:,:,:)
    real(kind=dp), allocatable :: density_derivative(:,:,:,:),base_batch(:,:,:,:), &
      base_fock_batch(:,:,:,:),response_fock_batch(:,:,:,:), &
      explicit_fock(:,:,:,:),total_fock(:,:,:,:),adjoint_derivative(:,:), &
      x_matrix(:,:)
    integer :: coordinate,nbf,nocca,noccb,ncoord,packed

    nbf=infos%basis%nbf; nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b; packed=nocca*(nbf-noccb)
    ncoord=size(dmo_a,3); status=0; dax=0.0_dp
    if(ncoord<=0 .or. nocca-noccb/=2 .or. size(x_packed)/=packed .or. &
       any(shape(dax)/=[packed,ncoord]) .or. &
       any(shape(dmo_a)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dmo_b)/=[nbf,nbf,ncoord])) then
      status=-1
      return
    end if
    ! The derivative-ERI primitive below is full-range.  A CAM/LC response
    ! requires separate short- and long-range derivative contractions and must
    ! not silently reuse it.
    if(infos%control%hamilton==20 .and. infos%dft%cam_flag) then
      status=-2
      return
    end if
    allocate(density(7,nbf,nbf),density_batch(ncoord,7,nbf,nbf), &
      density_derivative(7,nbf,nbf,ncoord),base_batch(1,7,nbf,nbf), &
      base_fock_batch(1,7,nbf,nbf),response_fock_batch(ncoord,7,nbf,nbf), &
      explicit_fock(7,nbf,nbf,ncoord),total_fock(7,nbf,nbf,ncoord), &
      adjoint_derivative(packed,ncoord),x_matrix(nbf,nbf),source=0.0_dp)
    call iatogen(x_packed,x_matrix,nocca,noccb)
    call mrsfcbc(infos,mo_a,mo_b,x_matrix,density)
    base_batch(1,:,:,:)=density
    call contract_mrsf_seven_density_batch(infos,int2_driver,base_batch, &
      response_scale,spc_coco,spc_ovov,spc_coov,base_fock_batch,status)
    if(status/=0) then
      call cleanup()
      return
    end if
    call build_mrsf_explicit_eri_derivative(infos,density,response_scale, &
      spc_coco,spc_ovov,spc_coov,infos%tddft%mult,explicit_fock,status)
    if(status/=0) then
      call cleanup()
      return
    end if
    call build_mrsf_density_orbital_derivative(infos,mo_a,mo_b,dmo_a,dmo_b, &
      x_packed,density_derivative,status)
    if(status/=0) then
      call cleanup()
      return
    end if
    do coordinate=1,ncoord
      density_batch(coordinate,:,:,:)=density_derivative(:,:,:,coordinate)
    end do
    call contract_mrsf_seven_density_batch(infos,int2_driver,density_batch, &
      response_scale,spc_coco,spc_ovov,spc_coov,response_fock_batch,status)
    if(status/=0) then
      call cleanup()
      return
    end if
    do coordinate=1,ncoord
      total_fock(:,:,:,coordinate)=explicit_fock(:,:,:,coordinate)+ &
        response_fock_batch(coordinate,:,:,:)
      call mrsfmntoia(infos,total_fock(:,:,:,coordinate),dax,mo_a,mo_b,coordinate)
    end do
    call build_mrsf_adjoint_orbital_derivative(infos,mo_a,mo_b,dmo_a,dmo_b, &
      base_fock_batch(1,:,:,:),adjoint_derivative,status)
    if(status==0) dax=dax+adjoint_derivative
    call cleanup()

  contains

    subroutine cleanup()
      if(allocated(density)) deallocate(density)
      if(allocated(density_batch)) deallocate(density_batch)
      if(allocated(density_derivative)) deallocate(density_derivative)
      if(allocated(base_batch)) deallocate(base_batch)
      if(allocated(base_fock_batch)) deallocate(base_fock_batch)
      if(allocated(response_fock_batch)) deallocate(response_fock_batch)
      if(allocated(explicit_fock)) deallocate(explicit_fock)
      if(allocated(total_fock)) deallocate(total_fock)
      if(allocated(adjoint_derivative)) deallocate(adjoint_derivative)
      if(allocated(x_matrix)) deallocate(x_matrix)
    end subroutine cleanup

  end subroutine build_mrsf_two_e_operator_derivative_action

end module tdhf_mrsf_hessian_two_e_action_mod
