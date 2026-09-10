module mod_dft_gridint_tdlda_point
  use precision, only: fp
  use mod_dft_gridint_tdxc_deriv, only: lda_tdxc_fixed_value, &
    lda_tdxc_fixed_der1, lda_tdxc_fixed_der2, lda_tdxc_response_probe
  implicit none
  private
  integer, parameter :: hmap(3,3)=reshape([1,4,6,4,2,5,6,5,3],[3,3])
  public :: lda_density_nuclear_point
  public :: lda_add_owner_motion
  public :: lda_weighted_fixed_hessian
  public :: lda_weighted_response_row
contains

  pure subroutine lda_density_nuclear_point(density,ao_atom,aov,aog1,aog2, &
      d1,d2)
    real(fp),intent(in)::density(:,:),aov(:),aog1(:,:),aog2(:,:)
    integer,intent(in)::ao_atom(:)
    real(fp),intent(out)::d1(:,:),d2(:,:,:,:)
    integer::nao,nat,mu,nu,aa,bb,a,b
    real(fp)::p,qm,qn,rm,rn,qrm,qrn
    nao=size(aov);nat=size(d1,2)
    if(any(shape(density)/=[nao,nao]).or.size(ao_atom)/=nao.or. &
       any(shape(aog1)/=[nao,3]).or.any(shape(aog2)/=[nao,6]).or. &
       any(shape(d1)/=[3,nat]).or.any(shape(d2)/=[3,3,nat,nat])) &
      error stop 'lda_density_nuclear_point: shape mismatch'
    d1=0.0_fp;d2=0.0_fp
    do nu=1,nao;do mu=1,nao
      p=density(mu,nu);if(p==0.0_fp)cycle
      do aa=1,nat;do a=1,3
        qm=0.0_fp;qn=0.0_fp
        if(ao_atom(mu)==aa)qm=-aog1(mu,a)
        if(ao_atom(nu)==aa)qn=-aog1(nu,a)
        d1(a,aa)=d1(a,aa)+p*(qm*aov(nu)+aov(mu)*qn)
        do bb=1,nat;do b=1,3
          rm=0.0_fp;rn=0.0_fp;qrm=0.0_fp;qrn=0.0_fp
          if(ao_atom(mu)==bb)rm=-aog1(mu,b)
          if(ao_atom(nu)==bb)rn=-aog1(nu,b)
          if(ao_atom(mu)==aa.and.aa==bb)qrm=aog2(mu,hmap(a,b))
          if(ao_atom(nu)==aa.and.aa==bb)qrn=aog2(nu,hmap(a,b))
          d2(a,b,aa,bb)=d2(a,b,aa,bb)+p*(qrm*aov(nu)+qm*rn+rm*qn+aov(mu)*qrn)
        end do;end do
      end do;end do
    end do;end do
  end subroutine lda_density_nuclear_point

  pure subroutine lda_add_owner_motion(owner,f1,f2,d1,d2)
    integer,intent(in)::owner
    real(fp),intent(in)::f1(:,:),f2(:,:,:,:)
    real(fp),intent(out)::d1(:,:),d2(:,:,:,:)
    integer::nat,aa,bb,a,b
    nat=size(f1,2);d1=f1;d2=f2
    do aa=1,nat;do a=1,3
      if(aa==owner)d1(a,aa)=d1(a,aa)-sum(f1(a,:))
      do bb=1,nat;do b=1,3
        if(aa==owner)d2(a,b,aa,bb)=d2(a,b,aa,bb)-sum(f2(b,a,bb,:),dim=1)
        if(bb==owner)d2(a,b,aa,bb)=d2(a,b,aa,bb)-sum(f2(a,b,aa,:),dim=1)
        if(aa==owner.and.bb==owner)d2(a,b,aa,bb)=d2(a,b,aa,bb)+sum(f2(a,b,:,:))
      end do;end do
    end do;end do
  end subroutine lda_add_owner_motion

  subroutine lda_weighted_fixed_hessian(vxc,fxc,kxc,lxc,p,x,dr,dp,dx, &
      d2r,d2p,d2x,qscale,w,dw,d2w,h)
    real(fp),intent(in)::vxc,fxc,kxc,lxc,p,x,qscale,w
    real(fp),intent(in)::dr(:,:),dp(:,:),dx(:,:),d2r(:,:,:,:),d2p(:,:,:,:),d2x(:,:,:,:)
    real(fp),intent(in)::dw(:,:),d2w(:,:,:,:)
    real(fp),intent(inout)::h(:,:,:,:)
    integer::nat,aa,bb,a,b
    real(fp)::q,qa,qb,qab
    nat=size(dr,2);q=lda_tdxc_fixed_value(vxc,fxc,p,x)
    do aa=1,nat;do a=1,3;do bb=1,nat;do b=1,3
      qa=lda_tdxc_fixed_der1(vxc,fxc,kxc,p,x,dr(a,aa),dp(a,aa),dx(a,aa))
      qb=lda_tdxc_fixed_der1(vxc,fxc,kxc,p,x,dr(b,bb),dp(b,bb),dx(b,bb))
      qab=lda_tdxc_fixed_der2(vxc,fxc,kxc,lxc,p,x,dr(a,aa),dr(b,bb), &
        d2r(a,b,aa,bb),dp(a,aa),dp(b,bb),d2p(a,b,aa,bb), &
        dx(a,aa),dx(b,bb),d2x(a,b,aa,bb))
      h(a,b,aa,bb)=h(a,b,aa,bb)+qscale*(w*qab+dw(a,aa)*qb+dw(b,bb)*qa+d2w(a,aa,b,bb)*q)
    end do;end do;end do;end do
  end subroutine lda_weighted_fixed_hessian

  subroutine lda_weighted_response_row(fxc,kxc,lxc,p,x,rd,rx,dr,dp,dx, &
      drd,drx,qscale,w,dw,row)
    real(fp),intent(in)::fxc,kxc,lxc,p,x,rd,rx,qscale,w
    real(fp),intent(in)::dr(:,:),dp(:,:),dx(:,:),drd(:,:),drx(:,:),dw(:,:)
    real(fp),intent(inout)::row(:,:)
    integer::nat,aa,a
    real(fp)::fr,dfr
    nat=size(dr,2);fr=lda_tdxc_response_probe(fxc,kxc,p,x,rd,rx)
    do aa=1,nat;do a=1,3
      dfr=kxc*dr(a,aa)*rd*p+fxc*drd(a,aa)*p+fxc*rd*dp(a,aa) &
        +4.0_fp*(lxc*dr(a,aa)*rd*x*x+kxc*drd(a,aa)*x*x &
        +2.0_fp*kxc*rd*x*dx(a,aa)) &
        +8.0_fp*(kxc*dr(a,aa)*x*rx+fxc*dx(a,aa)*rx+fxc*x*drx(a,aa))
      row(a,aa)=row(a,aa)+qscale*(dw(a,aa)*fr+w*dfr)
    end do;end do
  end subroutine lda_weighted_response_row
end module mod_dft_gridint_tdlda_point
