module tdhf_hessian_rows_mod
  use precision, only: dp
  implicit none
  private
  public :: build_tdhf_response_rows_hf
contains

  subroutine build_tdhf_response_rows_hf(infos, umat, eps_deriv, dground, &
      dprel, dwexc, duexc, dvexc, rows, rows_one, rows_two)
    use types, only: information
    use basis_tools, only: basis_set
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_VEC_MO_A, &
      OQP_E_MO_A, OQP_TD_P, OQP_TD_XPY, OQP_TD_XMY
    use mathlib, only: unpack_matrix, pack_matrix, symmetrize_matrix, orthogonal_transform
    use tdhf_lib, only: iatogen
    use grd1, only: grad_ee_overlap, grad_ee_kinetic, grad_en_hellman_feynman, &
      grad_en_pulay, grad_1e_ecp
    use messages, only: show_message, WITH_ABORT

    type(information),target,intent(inout)::infos
    real(dp),intent(in)::umat(:,:,:),eps_deriv(:,:),dground(:,:,:),dprel(:,:,:), &
      dwexc(:,:,:),duexc(:,:,:),dvexc(:,:,:)
    real(dp),intent(out)::rows(:,:),rows_one(:,:),rows_two(:,:)
    type(basis_set),pointer::basis
    real(dp),contiguous,pointer::dpk(:),c(:,:),eps(:),ppk(:,:),xpy(:,:),xmy(:,:)
    real(dp),allocatable,target::d0(:,:,:),p0(:,:,:),u0(:,:,:),v0(:,:,:)
    real(dp),allocatable::dplus(:,:,:),pplus(:,:,:),uplus(:,:,:),vplus(:,:,:)
    real(dp),allocatable::dminus(:,:,:),pminus(:,:,:),uminus(:,:,:),vminus(:,:,:)
    real(dp),allocatable::w0(:,:),dw0(:,:),dc(:,:),dens(:),gplus(:,:),gminus(:,:),gtmp(:,:),wrk(:,:)
    integer::i,k,nbf,nbf2,ncart,natom,nocc

    basis=>infos%basis; basis%atoms=>infos%atoms
    nbf=basis%nbf; nbf2=nbf*(nbf+1)/2; natom=size(basis%atoms%xyz,2); ncart=3*natom
    nocc=infos%mol_prop%nocc
    if(any(shape(rows)/=[ncart,ncart]) .or. any(shape(rows_one)/=[ncart,ncart]) .or. &
       any(shape(rows_two)/=[ncart,ncart])) &
      call show_message('TDHF response rows have wrong shape.',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_DM_A,dpk); call tagarray_get_data(infos%dat,OQP_VEC_MO_A,c)
    call tagarray_get_data(infos%dat,OQP_E_MO_A,eps); call tagarray_get_data(infos%dat,OQP_TD_P,ppk)
    call tagarray_get_data(infos%dat,OQP_TD_XPY,xpy); call tagarray_get_data(infos%dat,OQP_TD_XMY,xmy)
    allocate(d0(nbf,nbf,1),p0(nbf,nbf,1),u0(nbf,nbf,1),v0(nbf,nbf,1), &
             wrk(nbf,nbf),w0(nbf,nbf))
    call unpack_matrix(dpk,d0(:,:,1)); call unpack_matrix(ppk(:,1),p0(:,:,1))
    call iatogen(xpy(:,infos%tddft%target_state),wrk,nocc,nocc); call symmetrize_matrix(wrk,nbf); wrk=0.5_dp*wrk
    call orthogonal_transform('t',nbf,c,wrk,u0(:,:,1),w0)
    call iatogen(xmy(:,infos%tddft%target_state),wrk,nocc,nocc)
    call orthogonal_transform('t',nbf,c,wrk,v0(:,:,1),w0)
    allocate(dplus(nbf,nbf,1),pplus(nbf,nbf,1),uplus(nbf,nbf,1),vplus(nbf,nbf,1), &
             dminus(nbf,nbf,1),pminus(nbf,nbf,1),uminus(nbf,nbf,1),vminus(nbf,nbf,1))
    allocate(dw0(nbf,nbf),dc(nbf,nbf),dens(nbf2),gplus(3,natom),gminus(3,natom),gtmp(3,natom))
    w0=0.0_dp
    do i=1,nocc; w0=w0+2.0_dp*eps(i)*spread(c(:,i),2,nbf)*spread(c(:,i),1,nbf); end do
    rows=0.0_dp; rows_one=0.0_dp; rows_two=0.0_dp
    do k=1,ncart
      call pack_matrix(-2.0_dp*dwexc(:,:,k),dens)
      gplus=0.0_dp; call grad_ee_overlap(basis,dens,gplus)
      call pack_matrix(2.0_dp*dprel(:,:,k),dens)
      call grad_en_hellman_feynman(basis,basis%atoms%xyz,basis%atoms%zn-basis%ecp_zn_num,dens,gplus)
      call grad_ee_kinetic(basis,dens,gplus)
      call grad_en_pulay(basis,basis%atoms%xyz,basis%atoms%zn-basis%ecp_zn_num,dens,gplus)
      call grad_1e_ecp(infos,basis,basis%atoms%xyz,dens,gplus)
      rows_one(:,k)=reshape(gplus,[ncart])
      dplus(:,:,1)=d0(:,:,1); dminus(:,:,1)=d0(:,:,1)
      pplus(:,:,1)=p0(:,:,1)+dprel(:,:,k); pminus(:,:,1)=p0(:,:,1)-dprel(:,:,k)
      uplus(:,:,1)=u0(:,:,1)+duexc(:,:,k); uminus(:,:,1)=u0(:,:,1)-duexc(:,:,k)
      vplus(:,:,1)=v0(:,:,1)+dvexc(:,:,k); vminus(:,:,1)=v0(:,:,1)-dvexc(:,:,k)
      call two_e_gradient(infos,basis,dplus,pplus,uplus,vplus,gminus)
      call two_e_gradient(infos,basis,dminus,pminus,uminus,vminus,gtmp)
      rows_two(:,k)=reshape(0.5_dp*(gminus-gtmp),[ncart])
      ! The complete RHF Hessian supplies the D*D^R ground-state response,
      ! while the excitation part additionally contains P*D^R.  Isolate only
      ! that mixed term by changing the sign of D^R around zero; the quadratic
      ! D^R*D^R term is even and cancels, as do the unchanged X/Y terms.
      dplus(:,:,1)=dground(:,:,k); dminus(:,:,1)=-dground(:,:,k)
      pplus(:,:,1)=p0(:,:,1); pminus(:,:,1)=p0(:,:,1)
      uplus(:,:,1)=u0(:,:,1); uminus(:,:,1)=u0(:,:,1)
      vplus(:,:,1)=v0(:,:,1); vminus(:,:,1)=v0(:,:,1)
      call two_e_gradient(infos,basis,dplus,pplus,uplus,vplus,gminus)
      call two_e_gradient(infos,basis,dminus,pminus,uminus,vminus,gtmp)
      rows_two(:,k)=rows_two(:,k)+reshape(0.5_dp*(gminus-gtmp),[ncart])
      rows(:,k)=rows_one(:,k)+rows_two(:,k)
    end do
    deallocate(d0,p0,u0,v0,dplus,pplus,uplus,vplus,dminus,pminus,uminus,vminus, &
      w0,dw0,dc,dens,gplus,gminus,gtmp,wrk)
  contains
    subroutine two_e_gradient(inf,bas,d,p,u,v,g)
      use grd2, only: grd2_driver
      use tdhf_gradient_mod, only: grd2_tdhf_compute_data_t
      type(information),target,intent(inout)::inf; type(basis_set),intent(in)::bas
      real(dp),target,intent(inout)::d(:,:,:),p(:,:,:),u(:,:,:),v(:,:,:)
      real(dp),intent(out)::g(:,:)
      type(grd2_tdhf_compute_data_t)::gc
      real(dp)::scale_exch
      scale_exch=1.0_dp
      if(inf%control%hamilton>=20) scale_exch=inf%dft%hfscale
      g=0.0_dp; gc=grd2_tdhf_compute_data_t(d2=d,p2=p,xpy2=u,xmy2=v,hfscale=scale_exch,nbf=bas%nbf)
      call gc%init(); call gc%build_cart(bas); call grd2_driver(inf,bas,g,gc); call gc%clean()
    end subroutine
  end subroutine
end module tdhf_hessian_rows_mod
