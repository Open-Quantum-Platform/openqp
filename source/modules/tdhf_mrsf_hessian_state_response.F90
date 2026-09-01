module tdhf_mrsf_hessian_state_response_mod

  use precision, only: dp
  use types, only: information
  use int2_compute, only: int2_compute_t
  use tdhf_mrsf_hessian_fock_response_mod, only: &
    build_orbital_density_derivatives,build_rohf_total_fock_derivatives
  use tdhf_mrsf_hessian_operator_derivative_mod, only: &
    build_mrsf_operator_derivative_action
  use tdhf_mrsf_hessian_amplitude_mod, only: &
    solve_mrsf_tda_amplitude_derivatives
  use tdhf_mrsf_hessian_mo_response_mod, only: &
    canonicalize_rohf_common_response

  implicit none

  private
  public :: solve_mrsf_first_nuclear_response

contains

!###############################################################################

  subroutine solve_mrsf_first_nuclear_response(infos,int2_driver,mo_a,mo_b, &
      dmo_a,dmo_b,dmo_common,fock_a_ao,fock_b_ao,hcore_derivative,omega, &
      x_packed,dax,dx,domega,residual_max,status,dvxc_a,dvxc_b,dpa_out, &
      dpb_out,dfock_a_out,dfock_b_out)
    ! End-to-end first nuclear response of an isolated spin-adapted MRSF-TDA
    ! state.  ROHF/ROKS orbital response, total spin-Fock response, the seven
    ! density derivative action, and the projected amplitude equation meet at
    ! this boundary.  No numerical nuclear displacement is used here.

    type(information), target, intent(inout) :: infos
    type(int2_compute_t), intent(inout) :: int2_driver
    real(kind=dp), intent(in) :: mo_a(:,:),mo_b(:,:),dmo_a(:,:,:), &
      dmo_b(:,:,:)
    real(kind=dp), intent(inout) :: dmo_common(:,:,:)
    real(kind=dp), intent(in) :: fock_a_ao(:,:),fock_b_ao(:,:), &
      hcore_derivative(:,:,:),omega,x_packed(:)
    real(kind=dp), intent(out) :: dax(:,:),dx(:,:),domega(:),residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: dvxc_a(:,:,:),dvxc_b(:,:,:)
    real(kind=dp), intent(out), optional :: dpa_out(:,:,:),dpb_out(:,:,:), &
      dfock_a_out(:,:,:),dfock_b_out(:,:,:)

    real(kind=dp), allocatable :: pa(:,:),pb(:,:),dpa(:,:,:),dpb(:,:,:), &
      dfock_a(:,:,:),dfock_b(:,:,:),fa_mo(:,:),fb_mo(:,:)
    integer :: nbf,nocca,noccb,ncoord,packed,local_status
    logical :: dft

    nbf=infos%basis%nbf
    nocca=infos%mol_prop%nelec_a
    noccb=infos%mol_prop%nelec_b
    ncoord=size(dmo_a,3)
    packed=nocca*(nbf-noccb)
    dft=infos%control%hamilton==20
    status=0
    dax=0.0_dp
    dx=0.0_dp
    domega=0.0_dp
    residual_max=0.0_dp
    if(nbf<=0 .or. ncoord<=0 .or. nocca-noccb/=2 .or. &
       size(x_packed)/=packed .or. any(shape(dax)/=[packed,ncoord]) .or. &
       any(shape(dx)/=[packed,ncoord]) .or. size(domega)/=ncoord .or. &
       any(shape(dmo_common)/=[nbf,nbf,ncoord])) then
      status=-1
      return
    end if
    if(present(dpa_out)) then
      if(any(shape(dpa_out)/=[nbf,nbf,ncoord])) then
        status=-1
        return
      end if
    end if
    if(present(dpb_out)) then
      if(any(shape(dpb_out)/=[nbf,nbf,ncoord])) then
        status=-1
        return
      end if
    end if
    if(present(dfock_a_out)) then
      if(any(shape(dfock_a_out)/=[nbf,nbf,ncoord])) then
        status=-1
        return
      end if
    end if
    if(present(dfock_b_out)) then
      if(any(shape(dfock_b_out)/=[nbf,nbf,ncoord])) then
        status=-1
        return
      end if
    end if
    if(dft .and. (.not.present(dvxc_a) .or. .not.present(dvxc_b))) then
      status=-2
      return
    end if
    allocate(pa(nbf,nbf),pb(nbf,nbf),dpa(nbf,nbf,ncoord), &
      dpb(nbf,nbf,ncoord),dfock_a(nbf,nbf,ncoord), &
      dfock_b(nbf,nbf,ncoord),fa_mo(nbf,nbf),fb_mo(nbf,nbf),source=0.0_dp)
    pa=matmul(mo_a(:,1:nocca),transpose(mo_a(:,1:nocca)))
    pb=matmul(mo_b(:,1:noccb),transpose(mo_b(:,1:noccb)))
    call build_orbital_density_derivatives(mo_a,dmo_a,nocca,dpa,local_status)
    if(local_status==0) call build_orbital_density_derivatives( &
      mo_b,dmo_b,noccb,dpb,local_status)
    if(local_status/=0) status=-10+local_status
    if(status==0) then
      if(dft) then
        call build_rohf_total_fock_derivatives(infos,hcore_derivative,pa,pb, &
          dpa,dpb,dfock_a,dfock_b,local_status,dvxc_a,dvxc_b)
      else
        call build_rohf_total_fock_derivatives(infos,hcore_derivative,pa,pb, &
          dpa,dpb,dfock_a,dfock_b,local_status)
      end if
    end if
    if(status==0 .and. local_status/=0) status=-20+local_status
    if(status==0) call canonicalize_rohf_common_response(mo_a,dmo_common, &
      fock_a_ao,fock_b_ao,dfock_a,dfock_b,nocca,noccb,local_status)
    if(status==0 .and. local_status/=0) status=-25+local_status
    fa_mo=matmul(transpose(mo_a),matmul(fock_a_ao,mo_a))
    fb_mo=matmul(transpose(mo_b),matmul(fock_b_ao,mo_b))
    if(status==0) call build_mrsf_operator_derivative_action( &
      infos,int2_driver,mo_a,mo_b,dmo_common,dmo_common, &
      fock_a_ao,fock_b_ao, &
      dfock_a,dfock_b,x_packed,dax,local_status)
    if(status==0 .and. local_status/=0) status=-30+local_status
    if(status==0) call solve_mrsf_tda_amplitude_derivatives( &
      infos,int2_driver,mo_a,mo_b,fa_mo,fb_mo,omega,x_packed,dax,dx, &
      domega,residual_max,local_status)
    if(status==0 .and. local_status/=0) status=-40+local_status
    if(status==0) then
      if(present(dpa_out)) dpa_out=dpa
      if(present(dpb_out)) dpb_out=dpb
      if(present(dfock_a_out)) dfock_a_out=dfock_a
      if(present(dfock_b_out)) dfock_b_out=dfock_b
    end if
    deallocate(pa,pb,dpa,dpb,dfock_a,dfock_b,fa_mo,fb_mo)
  end subroutine solve_mrsf_first_nuclear_response

end module tdhf_mrsf_hessian_state_response_mod
