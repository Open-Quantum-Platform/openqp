module tdhf_mrsf_hessian_orbital_maps_mod

  use precision, only: dp
  use types, only: information
  use tdhf_lib, only: iatogen
  use tdhf_mrsf_lib, only: mrsfcbc,mrsfmntoia

  implicit none

  private
  public :: build_mrsf_density_orbital_derivative
  public :: build_mrsf_adjoint_orbital_derivative

contains

!###############################################################################

  subroutine build_mrsf_density_orbital_derivative(infos,mo_a,mo_b,dmo_a, &
      dmo_b,x_packed,density_derivative,status)
    ! Exact directional derivative of the spin-adapted seven-density map with
    ! respect to the ROHF orbital connection.  mrsfcbc is bilinear in the MO
    ! coefficients, so its centered polarization identity is algebraically
    ! exact: quadratic terms in dC cancel without a nuclear displacement.

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:),dmo_b(:,:,:)
    real(kind=dp), intent(in) :: x_packed(:)
    real(kind=dp), intent(out) :: density_derivative(:,:,:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: x_matrix(:,:),plus_density(:,:,:), &
      minus_density(:,:,:),plus_a(:,:),minus_a(:,:),plus_b(:,:),minus_b(:,:)
    integer :: coordinate,nbf,nocca,noccb,ncoord,packed

    nbf=infos%basis%nbf; nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b; packed=nocca*(nbf-noccb)
    ncoord=size(dmo_a,3); status=0; density_derivative=0.0_dp
    if(ncoord<=0 .or. nocca-noccb/=2 .or. size(x_packed)/=packed .or. &
       any(shape(mo_a)/=[nbf,nbf]) .or. any(shape(mo_b)/=[nbf,nbf]) .or. &
       any(shape(dmo_a)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dmo_b)/=[nbf,nbf,ncoord]) .or. &
       any(shape(density_derivative)/=[7,nbf,nbf,ncoord])) then
      status=-1
      return
    end if
    allocate(x_matrix(nbf,nbf),plus_density(7,nbf,nbf), &
      minus_density(7,nbf,nbf),plus_a(nbf,nbf),minus_a(nbf,nbf), &
      plus_b(nbf,nbf),minus_b(nbf,nbf))
    call iatogen(x_packed,x_matrix,nocca,noccb)
    do coordinate=1,ncoord
      plus_a=mo_a+dmo_a(:,:,coordinate); minus_a=mo_a-dmo_a(:,:,coordinate)
      plus_b=mo_b+dmo_b(:,:,coordinate); minus_b=mo_b-dmo_b(:,:,coordinate)
      call mrsfcbc(infos,plus_a,plus_b,x_matrix,plus_density)
      call mrsfcbc(infos,minus_a,minus_b,x_matrix,minus_density)
      density_derivative(:,:,:,coordinate)=0.5_dp*(plus_density-minus_density)
    end do
    deallocate(x_matrix,plus_density,minus_density,plus_a,minus_a,plus_b,minus_b)
  end subroutine build_mrsf_density_orbital_derivative

!###############################################################################

  subroutine build_mrsf_adjoint_orbital_derivative(infos,mo_a,mo_b,dmo_a, &
      dmo_b,channel_fock,packed_derivative,status)
    ! Exact directional derivative of the adjoint seven-channel map
    ! mrsfmntoia with respect to the orbital connection, again evaluated by an
    ! algebraic polarization identity at fixed AO channel Fock matrices.

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:),dmo_b(:,:,:)
    real(kind=dp), intent(in) :: channel_fock(:,:,:)
    real(kind=dp), intent(out) :: packed_derivative(:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: plus_result(:,:),minus_result(:,:),plus_a(:,:), &
      minus_a(:,:),plus_b(:,:),minus_b(:,:)
    integer :: coordinate,nbf,nocca,noccb,ncoord,packed

    nbf=infos%basis%nbf; nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b; packed=nocca*(nbf-noccb)
    ncoord=size(dmo_a,3); status=0; packed_derivative=0.0_dp
    if(ncoord<=0 .or. nocca-noccb/=2 .or. &
       any(shape(channel_fock)/=[7,nbf,nbf]) .or. &
       any(shape(dmo_a)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dmo_b)/=[nbf,nbf,ncoord]) .or. &
       any(shape(packed_derivative)/=[packed,ncoord])) then
      status=-1
      return
    end if
    allocate(plus_result(packed,1),minus_result(packed,1), &
      plus_a(nbf,nbf),minus_a(nbf,nbf),plus_b(nbf,nbf),minus_b(nbf,nbf))
    do coordinate=1,ncoord
      plus_result=0.0_dp; minus_result=0.0_dp
      plus_a=mo_a+dmo_a(:,:,coordinate); minus_a=mo_a-dmo_a(:,:,coordinate)
      plus_b=mo_b+dmo_b(:,:,coordinate); minus_b=mo_b-dmo_b(:,:,coordinate)
      call mrsfmntoia(infos,channel_fock,plus_result,plus_a,plus_b,1)
      call mrsfmntoia(infos,channel_fock,minus_result,minus_a,minus_b,1)
      packed_derivative(:,coordinate)=0.5_dp*(plus_result(:,1)-minus_result(:,1))
    end do
    deallocate(plus_result,minus_result,plus_a,minus_a,plus_b,minus_b)
  end subroutine build_mrsf_adjoint_orbital_derivative

end module tdhf_mrsf_hessian_orbital_maps_mod
