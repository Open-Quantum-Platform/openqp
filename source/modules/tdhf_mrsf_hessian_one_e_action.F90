module tdhf_mrsf_hessian_one_e_action_mod

  use precision, only: dp
  use types, only: information
  use tdhf_lib, only: iatogen
  use tdhf_mrsf_lib, only: mrsfesum

  implicit none

  private
  public :: transform_mo_operator_derivative
  public :: build_mrsf_one_e_operator_derivative_action

contains

!###############################################################################

  subroutine transform_mo_operator_derivative(mo,dmo,operator_ao, &
      operator_ao_derivative,operator_mo_derivative,status)
    ! Total derivative of C^T F C in a moving atom-centred basis.  dC is the
    ! full orbital connection satisfying differentiated orthonormality; dF is
    ! the total AO derivative, including explicit integral and relaxed-density
    ! response terms supplied by the ROHF/ROKS nuclear-response builder.

    real(kind=dp), intent(in) :: mo(:,:),dmo(:,:,:),operator_ao(:,:), &
      operator_ao_derivative(:,:,:)
    real(kind=dp), intent(out) :: operator_mo_derivative(:,:,:)
    integer, intent(out) :: status

    integer :: coordinate,nbf,ncoord

    nbf=size(mo,1)
    ncoord=size(dmo,3)
    status=0
    operator_mo_derivative=0.0_dp
    if(nbf<=0 .or. ncoord<=0 .or. size(mo,2)/=nbf .or. &
       any(shape(dmo)/=[nbf,nbf,ncoord]) .or. &
       any(shape(operator_ao)/=[nbf,nbf]) .or. &
       any(shape(operator_ao_derivative)/=[nbf,nbf,ncoord]) .or. &
       any(shape(operator_mo_derivative)/=[nbf,nbf,ncoord])) then
      status=-1
      return
    end if
    do coordinate=1,ncoord
      operator_mo_derivative(:,:,coordinate)= &
        matmul(transpose(dmo(:,:,coordinate)),matmul(operator_ao,mo))+ &
        matmul(transpose(mo),matmul(operator_ao_derivative(:,:,coordinate),mo))+ &
        matmul(transpose(mo),matmul(operator_ao,dmo(:,:,coordinate)))
    end do
  end subroutine transform_mo_operator_derivative

!###############################################################################

  subroutine build_mrsf_one_e_operator_derivative_action(infos,mo_a,mo_b, &
      dmo_a,dmo_b,fock_a_ao,fock_b_ao,dfock_a_ao,dfock_b_ao,x_packed,dax, &
      status)
    ! One-electron contribution to (dA/dR)X in the packed MRSF response space.
    ! mrsfesum is linear in the alpha and beta MO-basis Fock matrices, so the
    ! exact directional derivative at fixed physical amplitude is obtained by
    ! applying it once to d(C_a^T F_a C_a) and d(C_b^T F_b C_b).
    !
    ! The electronic response remains the founding spin-adapted two-SOMO
    ! construction initiated by Hiroya Nakata; no configuration expansion is
    ! introduced at this derivative boundary.

    type(information), intent(in) :: infos
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:),dmo_b(:,:,:)
    real(kind=dp), intent(in) :: fock_a_ao(:,:),fock_b_ao(:,:), &
      dfock_a_ao(:,:,:),dfock_b_ao(:,:,:),x_packed(:)
    real(kind=dp), intent(out) :: dax(:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: x_matrix(:,:),dfa_mo(:,:,:),dfb_mo(:,:,:)
    integer :: coordinate,nbf,nocca,noccb,ncoord,packed,local_status

    nbf=infos%basis%nbf
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    packed=nocca*(nbf-noccb)
    ncoord=size(dmo_a,3)
    status=0
    dax=0.0_dp
    if(nbf<=0 .or. ncoord<=0 .or. nocca-noccb/=2 .or. &
       size(x_packed)/=packed .or. any(shape(dax)/=[packed,ncoord]) .or. &
       any(shape(mo_a)/=[nbf,nbf]) .or. any(shape(mo_b)/=[nbf,nbf]) .or. &
       any(shape(dmo_a)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dmo_b)/=[nbf,nbf,ncoord]) .or. &
       any(shape(fock_a_ao)/=[nbf,nbf]) .or. &
       any(shape(fock_b_ao)/=[nbf,nbf]) .or. &
       any(shape(dfock_a_ao)/=[nbf,nbf,ncoord]) .or. &
       any(shape(dfock_b_ao)/=[nbf,nbf,ncoord])) then
      status=-1
      return
    end if
    allocate(x_matrix(nbf,nbf),dfa_mo(nbf,nbf,ncoord), &
      dfb_mo(nbf,nbf,ncoord),source=0.0_dp)
    call transform_mo_operator_derivative(mo_a,dmo_a,fock_a_ao, &
      dfock_a_ao,dfa_mo,local_status)
    if(local_status==0) call transform_mo_operator_derivative(mo_b,dmo_b, &
      fock_b_ao,dfock_b_ao,dfb_mo,local_status)
    if(local_status/=0) then
      status=-2
      deallocate(x_matrix,dfa_mo,dfb_mo)
      return
    end if
    call iatogen(x_packed,x_matrix,nocca,noccb)
    do coordinate=1,ncoord
      call mrsfesum(infos,x_matrix,dfa_mo(:,:,coordinate), &
        dfb_mo(:,:,coordinate),dax,coordinate)
    end do
    deallocate(x_matrix,dfa_mo,dfb_mo)
  end subroutine build_mrsf_one_e_operator_derivative_action

end module tdhf_mrsf_hessian_one_e_action_mod
