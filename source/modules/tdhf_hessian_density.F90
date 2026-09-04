module tdhf_hessian_density_mod

  use precision, only: dp
  implicit none
  private
  public :: build_tdhf_relaxed_density_derivatives

contains

  subroutine build_tdhf_relaxed_density_derivatives(infos, umat, eps_deriv, &
      omega, domega, u0, v0, du, dv, z0, dz, dp_relaxed, dw_ao, du_ao, dv_ao)
    use types, only: information
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_VEC_MO_A, OQP_E_MO_A
    use tdhf_hessian_z_rhs_mod, only: differentiated_channel
    use tdhf_hessian_gxc_derivative_mod, only: build_gxc_and_derivative
    use messages, only: show_message, WITH_ABORT

    type(information), target, intent(inout) :: infos
    real(kind=dp), intent(in) :: umat(:,:,:), eps_deriv(:,:), omega, domega(:)
    real(kind=dp), intent(in) :: u0(:), v0(:), du(:,:), dv(:,:), z0(:), dz(:,:)
    real(kind=dp), intent(out) :: dp_relaxed(:,:,:), dw_ao(:,:,:)
    real(kind=dp), intent(out) :: du_ao(:,:,:), dv_ao(:,:,:)

    real(kind=dp), contiguous, pointer :: c(:,:), eps(:)
    real(kind=dp), allocatable :: um(:,:), vm(:,:), zm(:,:), tm(:,:), pm(:,:), wm(:,:)
    real(kind=dp), allocatable :: dum(:,:), dvm(:,:), dzm(:,:), dtm(:,:,:), dpm(:,:,:), dwm(:,:,:)
    real(kind=dp), allocatable :: hp(:,:), hm(:,:), ht(:,:), dhp(:,:,:), dhm(:,:,:), dht(:,:,:)
    real(kind=dp), allocatable :: gxp(:,:), dgxp(:,:,:)
    real(kind=dp), allocatable :: dpz_ao(:,:,:)
    real(kind=dp), allocatable :: mu(:,:), mv(:,:), dm_u(:,:,:), dm_v(:,:,:)
    real(kind=dp), allocatable :: work(:,:)
    integer :: k, nbf, ncoord, nocc, nvir, nexc

    nbf=infos%basis%nbf; nocc=infos%mol_prop%nocc; nvir=nbf-nocc
    nexc=nocc*nvir; ncoord=size(umat,3)
    if (size(z0)/=nexc .or. any(shape(dz)/=[nexc,ncoord]) .or. &
        any(shape(dp_relaxed)/=[nbf,nbf,ncoord]) .or. &
        any(shape(dw_ao)/=[nbf,nbf,ncoord])) &
      call show_message('TDHF relaxed-density derivative arrays have incompatible shapes.',WITH_ABORT)
    call tagarray_get_data(infos%dat,OQP_VEC_MO_A,c)
    call tagarray_get_data(infos%dat,OQP_E_MO_A,eps)
    allocate(um(nocc,nvir),vm(nocc,nvir),zm(nocc,nvir),tm(nbf,nbf),pm(nbf,nbf),wm(nbf,nbf))
    allocate(dum(nocc,nvir),dvm(nocc,nvir),dzm(nocc,nvir),dtm(nbf,nbf,ncoord), &
             dpm(nbf,nbf,ncoord),dwm(nbf,nbf,ncoord),source=0.0_dp)
    um=reshape(u0,[nocc,nvir]); vm=reshape(v0,[nocc,nvir]); zm=reshape(z0,[nocc,nvir])
    tm=0.0_dp
    call trans_square(tm,um,um,0.5_dp); call trans_square(tm,vm,vm,0.5_dp)
    do k=1,ncoord
      dum=reshape(du(:,k),[nocc,nvir]); dvm=reshape(dv(:,k),[nocc,nvir])
      call trans_square(dtm(:,:,k),dum,um,1.0_dp)
      call trans_square(dtm(:,:,k),dvm,vm,1.0_dp)
    end do
    ! The production gradient forms the energy-weighted density with H+[P],
    ! where P is the relaxed difference density T+Z (not H+[T] alone).
    ! iatogen/symmetrize/0.5 stores one half of Z in both OV blocks.
    pm=tm
    pm(1:nocc,nocc+1:)=pm(1:nocc,nocc+1:)+0.5_dp*zm
    pm(nocc+1:,1:nocc)=pm(nocc+1:,1:nocc)+0.5_dp*transpose(zm)
    dpm=dtm
    do k=1,ncoord
      dzm=reshape(dz(:,k),[nocc,nvir])
      dpm(1:nocc,nocc+1:,k)=dpm(1:nocc,nocc+1:,k)+0.5_dp*dzm
      dpm(nocc+1:,1:nocc,k)=dpm(nocc+1:,1:nocc,k)+0.5_dp*transpose(dzm)
    end do
    allocate(mu(nbf,nbf),mv(nbf,nbf),dm_u(nbf,nbf,ncoord),dm_v(nbf,nbf,ncoord),source=0.0_dp)
    mu(1:nocc,nocc+1:)=um; mv(1:nocc,nocc+1:)=vm
    do k=1,ncoord
      dm_u(1:nocc,nocc+1:,k)=reshape(du(:,k),[nocc,nvir])
      dm_v(1:nocc,nocc+1:,k)=reshape(dv(:,k),[nocc,nvir])
    end do
    allocate(hp(nbf,nbf),hm(nbf,nbf),ht(nbf,nbf),dhp(nbf,nbf,ncoord), &
             dhm(nbf,nbf,ncoord),dht(nbf,nbf,ncoord))
    call differentiated_channel(infos,c,umat,mu,dm_u,+1,hp,dhp)
    call differentiated_channel(infos,c,umat,mv,dm_v,-1,hm,dhm)
    call differentiated_channel(infos,c,umat,pm,dpm,+1,ht,dht)
    allocate(gxp(nbf,nbf),dgxp(nbf,nbf,ncoord),source=0.0_dp)
    if (infos%control%hamilton==20) &
      call build_gxc_and_derivative(infos,c,umat,um,du,gxp,dgxp)
    call form_w_and_derivative(um,vm,zm,du,dv,dz,eps,eps_deriv,omega,domega, &
                               hp,hm,ht,dhp,dhm,dht,gxp,dgxp,wm,dwm)

    allocate(work(nbf,nbf),dpz_ao(nbf,nbf,ncoord))
    do k=1,ncoord
      dum=reshape(du(:,k),[nocc,nvir]); dvm=reshape(dv(:,k),[nocc,nvir]); &
      dzm=reshape(dz(:,k),[nocc,nvir])
      call density_derivative(c,umat(:,:,k),tm,dtm(:,:,k),dp_relaxed(:,:,k),work)
      call nonsymmetric_vo_derivative(c,umat(:,:,k),zm,dzm,dpz_ao(:,:,k),work)
      dp_relaxed(:,:,k)=dp_relaxed(:,:,k)+dpz_ao(:,:,k)
      call symmetric_ov_derivative(c,umat(:,:,k),um,dum,du_ao(:,:,k),work)
      call nonsymmetric_ov_derivative(c,umat(:,:,k),vm,dvm,dv_ao(:,:,k),work)
      call density_derivative(c,umat(:,:,k),wm,dwm(:,:,k),dw_ao(:,:,k),work)
      dw_ao(:,:,k)=0.5_dp*dw_ao(:,:,k)
    end do
    deallocate(um,vm,zm,tm,pm,wm,dum,dvm,dzm,dtm,dpm,dwm,hp,hm,ht,dhp,dhm,dht,gxp,dgxp, &
               mu,mv,dm_u,dm_v,work,dpz_ao)
  end subroutine

  subroutine form_w_and_derivative(u,v,z,duv,dvv,dzv,eps,deps,omega,domega, &
                                    hp,hm,ht,dhp,dhm,dht,gxp,dgxp,w,dw)
    real(dp),intent(in)::u(:,:),v(:,:),z(:,:),duv(:,:),dvv(:,:),dzv(:,:),eps(:),deps(:,:)
    real(dp),intent(in)::omega,domega(:),hp(:,:),hm(:,:),ht(:,:),dhp(:,:,:),dhm(:,:,:),dht(:,:,:)
    real(dp),intent(in)::gxp(:,:),dgxp(:,:,:)
    real(dp),intent(out)::w(:,:),dw(:,:,:)
    real(dp),allocatable::du(:,:),dv(:,:),dz(:,:)
    integer::a,b,i,j,k,no,nv
    no=size(u,1); nv=size(u,2); w=0.0_dp; dw=0.0_dp
    allocate(du(no,nv),dv(no,nv),dz(no,nv))
    do j=1,no; do i=1,no
      w(i,j)=omega*(dot_product(u(i,:),v(j,:))+dot_product(v(i,:),u(j,:)))+ht(i,j)
      do a=1,nv; w(i,j)=w(i,j)-eps(no+a)*(u(i,a)*u(j,a)+v(i,a)*v(j,a)); end do
    end do; end do
    do b=1,nv; do a=1,nv
      w(no+a,no+b)=omega*(dot_product(u(:,a),v(:,b))+dot_product(v(:,a),u(:,b)))
      do i=1,no; w(no+a,no+b)=w(no+a,no+b)+eps(i)*(u(i,a)*u(i,b)+v(i,a)*v(i,b)); end do
    end do; end do
    w(1:no,no+1:)=matmul(transpose(hp(1:no,1:no)),u)+matmul(transpose(hm(1:no,1:no)),v)
    do i=1,no; w(i,no+1:)=w(i,no+1:)+eps(i)*z(i,:); end do
    w(no+1:,1:no)=transpose(w(1:no,no+1:))
    ! The production TDDFT gradient adds 2*Gxc[X+Y,X+Y] to W_oo.
    ! Its complete derivative contains the explicit grid derivative, density
    ! response, and both inner MO rotations (GAMESS dGXM_K).
    w(1:no,1:no)=w(1:no,1:no)+2.0_dp*gxp(1:no,1:no)
    do k=1,size(dw,3)
      du=reshape(duv(:,k),[no,nv]); dv=reshape(dvv(:,k),[no,nv]); dz=reshape(dzv(:,k),[no,nv])
      do j=1,no; do i=1,no
        dw(i,j,k)=domega(k)*(dot_product(u(i,:),v(j,:))+dot_product(v(i,:),u(j,:)))+dht(i,j,k) &
          +omega*(dot_product(du(i,:),v(j,:))+dot_product(u(i,:),dv(j,:)) &
                 +dot_product(dv(i,:),u(j,:))+dot_product(v(i,:),du(j,:)))
        do a=1,nv
          dw(i,j,k)=dw(i,j,k)-deps(no+a,k)*(u(i,a)*u(j,a)+v(i,a)*v(j,a)) &
            -eps(no+a)*(du(i,a)*u(j,a)+u(i,a)*du(j,a)+dv(i,a)*v(j,a)+v(i,a)*dv(j,a))
        end do
      end do; end do
      do b=1,nv; do a=1,nv
        dw(no+a,no+b,k)=domega(k)*(dot_product(u(:,a),v(:,b))+dot_product(v(:,a),u(:,b))) &
          +omega*(dot_product(du(:,a),v(:,b))+dot_product(u(:,a),dv(:,b)) &
                 +dot_product(dv(:,a),u(:,b))+dot_product(v(:,a),du(:,b)))
        do i=1,no
          dw(no+a,no+b,k)=dw(no+a,no+b,k)+deps(i,k)*(u(i,a)*u(i,b)+v(i,a)*v(i,b)) &
            +eps(i)*(du(i,a)*u(i,b)+u(i,a)*du(i,b)+dv(i,a)*v(i,b)+v(i,a)*dv(i,b))
        end do
      end do; end do
      dw(1:no,no+1:,k)=matmul(transpose(dhp(1:no,1:no,k)),u)+matmul(transpose(hp(1:no,1:no)),du) &
        +matmul(transpose(dhm(1:no,1:no,k)),v)+matmul(transpose(hm(1:no,1:no)),dv)
      do i=1,no; dw(i,no+1:,k)=dw(i,no+1:,k)+deps(i,k)*z(i,:)+eps(i)*dz(i,:); end do
      dw(no+1:,1:no,k)=transpose(dw(1:no,no+1:,k))
      dw(1:no,1:no,k)=dw(1:no,1:no,k)+2.0_dp*dgxp(1:no,1:no,k)
    end do
    deallocate(du,dv,dz)
  end subroutine

  pure subroutine trans_square(t,a,b,f)
    real(dp),intent(inout)::t(:,:); real(dp),intent(in)::a(:,:),b(:,:),f; integer::no
    no=size(a,1)
    t(:no,:no)=t(:no,:no)-0.5_dp*f*(matmul(a,transpose(b))+matmul(b,transpose(a)))
    t(no+1:,no+1:)=t(no+1:,no+1:)+0.5_dp*f*(matmul(transpose(a),b)+matmul(transpose(b),a))
  end subroutine

  subroutine density_derivative(c,u,m,dm,dao,w)
    real(dp),intent(in)::c(:,:),u(:,:),m(:,:),dm(:,:); real(dp),intent(out)::dao(:,:),w(:,:)
    w=matmul(c,u)
    dao=matmul(matmul(w,m),transpose(c))+matmul(matmul(c,dm),transpose(c)) &
       +matmul(matmul(c,m),transpose(w))
  end subroutine

  subroutine add_symmetric_ov_density(c,u,z,dz,dao,w,scale)
    real(dp),intent(in)::c(:,:),u(:,:),z(:,:),dz(:,:),scale; real(dp),intent(inout)::dao(:,:); real(dp),intent(out)::w(:,:)
    real(dp)::m(size(c,2),size(c,2)),dm(size(c,2),size(c,2)),x(size(c,1),size(c,1)); integer::no
    no=size(z,1); m=0; dm=0; m(:no,no+1:)=0.5_dp*z; m(no+1:,:no)=0.5_dp*transpose(z)
    dm(:no,no+1:)=0.5_dp*dz; dm(no+1:,:no)=0.5_dp*transpose(dz)
    call density_derivative(c,u,m,dm,x,w); dao=dao+scale*x
  end subroutine

  subroutine symmetric_ov_derivative(c,u,z,dz,dao,w)
    real(dp),intent(in)::c(:,:),u(:,:),z(:,:),dz(:,:); real(dp),intent(out)::dao(:,:),w(:,:)
    dao=0; call add_symmetric_ov_density(c,u,z,dz,dao,w,1.0_dp)
  end subroutine

  subroutine nonsymmetric_ov_derivative(c,u,z,dz,dao,w)
    real(dp),intent(in)::c(:,:),u(:,:),z(:,:),dz(:,:); real(dp),intent(out)::dao(:,:),w(:,:)
    real(dp)::m(size(c,2),size(c,2)),dm(size(c,2),size(c,2)); integer::no
    no=size(z,1); m=0; dm=0; m(:no,no+1:)=z; dm(:no,no+1:)=dz
    call density_derivative(c,u,m,dm,dao,w)
  end subroutine

  subroutine nonsymmetric_vo_derivative(c,u,z,dz,dao,w)
    real(dp),intent(in)::c(:,:),u(:,:),z(:,:),dz(:,:); real(dp),intent(out)::dao(:,:),w(:,:)
    real(dp)::m(size(c,2),size(c,2)),dm(size(c,2),size(c,2)); integer::no
    no=size(z,1); m=0; dm=0
    m(no+1:,:no)=transpose(z); dm(no+1:,:no)=transpose(dz)
    call density_derivative(c,u,m,dm,dao,w)
  end subroutine

end module tdhf_hessian_density_mod
