module mod_dft_gridint_tdlda_driver
  use precision,only:fp
  use mod_dft_gridint_tdlda_point,only:lda_density_nuclear_point, &
    lda_add_owner_motion,lda_weighted_fixed_hessian,lda_weighted_response_row
  use mod_dft_partition_hessian,only:partition_weight_nuclear_derivatives
  use mod_dft_gridint,only:xc_engine_t,xc_consumer_t,xc_options_t,run_xc, &
    xc_der1,xc_der2_contr,xc_der3_contr,OQP_FUNTYP_LDA
  use functionals,only:functional_t
  implicit none
  private
  type,extends(xc_consumer_t)::xc_consumer_tdlda_t
    real(fp),allocatable::h(:,:,:,:,:),rows(:,:,:,:),p(:,:),x(:,:),dd(:,:,:),du(:,:,:)
    real(fp),allocatable::xyz(:,:),shift(:,:)
    integer,allocatable::ao_atom(:)
    logical,allocatable::dummy(:)
    integer::part_type=0
    logical::has_shift=.false.,response=.false.
    class(functional_t),pointer::functional=>null()
  contains
    procedure::parallel_start=>lda_start
    procedure::parallel_stop=>lda_stop
    procedure::update=>lda_update
    procedure::postUpdate=>lda_post
    procedure::clean=>lda_clean
  end type
  public::tddft_lda_fixed_hessian,tddft_lda_response_rows
contains
  subroutine tddft_lda_fixed_hessian(basis,grid,r,p,x,h,infos)
    use basis_tools,only:basis_set
    use mod_dft_molgrid,only:dft_grid_t
    use types,only:information
    type(basis_set),intent(in)::basis;type(dft_grid_t),target,intent(in)::grid
    real(fp),intent(in)::r(:,:),p(:,:),x(:,:);real(fp),intent(out)::h(:,:)
    type(information),target,intent(in)::infos
    type(xc_consumer_tdlda_t)::dat
    call setup_and_run(dat,basis,grid,r,p,x,infos,.false.)
    call export_hessian(dat,h);call dat%clean()
  end subroutine
  subroutine tddft_lda_response_rows(basis,grid,r,p,x,dd,du,rows,infos)
    use basis_tools,only:basis_set
    use mod_dft_molgrid,only:dft_grid_t
    use types,only:information
    type(basis_set),intent(in)::basis;type(dft_grid_t),target,intent(in)::grid
    real(fp),intent(in)::r(:,:),p(:,:),x(:,:),dd(:,:,:),du(:,:,:)
    real(fp),intent(out)::rows(:,:);type(information),target,intent(in)::infos
    type(xc_consumer_tdlda_t)::dat
    call setup_and_run(dat,basis,grid,r,p,x,infos,.true.,dd,du)
    call export_rows(dat,rows);call dat%clean()
  end subroutine
  subroutine setup_and_run(dat,basis,grid,r,p,x,infos,response,dd,du)
    use basis_tools,only:basis_set
    use mod_dft_molgrid,only:dft_grid_t
    use types,only:information
    type(xc_consumer_tdlda_t),intent(inout)::dat;type(basis_set),intent(in)::basis
    type(dft_grid_t),target,intent(in)::grid;type(information),target,intent(in)::infos
    real(fp),intent(in)::r(:,:),p(:,:),x(:,:);logical,intent(in)::response
    real(fp),intent(in),optional::dd(:,:,:),du(:,:,:)
    type(xc_options_t)::opts;real(fp),allocatable,target::rn(:,:)
    integer::n,nat,nc,i,k,s,first,last
    n=basis%nbf;nat=infos%mol_prop%natom;nc=3*nat;dat%response=response
    allocate(rn(n,n),dat%p(n,n),dat%x(n,n),dat%ao_atom(n))
    if(response)allocate(dat%dd(n,n,nc),dat%du(n,n,nc))
    do i=1,n
      rn(:,i)=r(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      dat%p(:,i)=2.0_fp*p(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      dat%x(:,i)=x(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      if(response)then;do k=1,nc
        dat%dd(:,i,k)=dd(:,i,k)*basis%bfnrm(:)*basis%bfnrm(i)
        dat%du(:,i,k)=du(:,i,k)*basis%bfnrm(:)*basis%bfnrm(i)
      end do;end if
    end do
    dat%ao_atom=0;do s=1,basis%nshell;first=basis%ao_offset(s);last=first+basis%naos(s)-1
      dat%ao_atom(first:last)=basis%origin(s);end do
    dat%xyz=infos%atoms%xyz;dat%dummy=grid%dummyAtom;dat%shift=grid%surfaceShift
    dat%part_type=grid%partFunType;dat%has_shift=grid%hasSurfaceShift;dat%functional=>infos%functional
    opts%isGGA=.false.;opts%needTau=.false.;opts%functional=>infos%functional
    opts%hasBeta=.false.;opts%isWFVecs=.false.;opts%numAOs=n;opts%maxPts=grid%maxSlicePts
    opts%limPts=grid%maxNRadTimesNAng;opts%numAtoms=nat;opts%maxAngMom=basis%mxam
    opts%nDer=2;opts%nXCDer=3;opts%numOccAlpha=infos%mol_prop%nelec_A
    opts%numOccBeta=infos%mol_prop%nelec_B;opts%wfAlpha=>rn
    opts%dft_threshold=infos%dft%grid_density_cutoff;opts%ao_threshold=infos%dft%grid_ao_threshold
    opts%ao_sparsity_ratio=0.0_fp;opts%molGrid=>grid
    call dat%pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi);call run_xc(opts,dat,basis)
  end subroutine
  subroutine export_hessian(dat,out)
    type(xc_consumer_tdlda_t),intent(in)::dat;real(fp),intent(out)::out(:,:)
    integer::a,b,i,j;do b=1,size(dat%h,4);do j=1,3;do a=1,size(dat%h,3);do i=1,3
      out(3*(a-1)+i,3*(b-1)+j)=dat%h(i,j,a,b,1)
    end do;end do;end do;end do
  end subroutine
  subroutine export_rows(dat,out)
    type(xc_consumer_tdlda_t),intent(in)::dat;real(fp),intent(out)::out(:,:)
    integer::a,i,k;do k=1,size(out,2);do a=1,size(dat%rows,2);do i=1,3
      out(3*(a-1)+i,k)=dat%rows(i,a,k,1)
    end do;end do;end do
  end subroutine
  subroutine lda_start(self,xce,nthreads)
    class(xc_consumer_tdlda_t),target,intent(inout)::self;class(xc_engine_t),intent(in)::xce;integer,intent(in)::nthreads
    if(self%response)then;allocate(self%rows(3,xce%numAtoms,3*xce%numAtoms,nthreads),source=0.0_fp)
    else;allocate(self%h(3,3,xce%numAtoms,xce%numAtoms,nthreads),source=0.0_fp);end if
  end subroutine
  subroutine lda_stop(self)
    class(xc_consumer_tdlda_t),intent(inout)::self
    if(self%response)then
      if(size(self%rows,4)>1)self%rows(:,:,:,1)=sum(self%rows,dim=4)
      call self%pe%allreduce(self%rows(:,:,:,1),size(self%rows(:,:,:,1)))
    else
      if(size(self%h,5)>1)self%h(:,:,:,:,1)=sum(self%h,dim=5)
      call self%pe%allreduce(self%h(:,:,:,:,1),size(self%h(:,:,:,:,1)))
    end if
  end subroutine
  subroutine lda_post(self,xce,mythread)
    class(xc_consumer_tdlda_t),intent(inout)::self;class(xc_engine_t),intent(in)::xce;integer::mythread
  end subroutine
  subroutine lda_clean(self)
    class(xc_consumer_tdlda_t),intent(inout)::self
    if(allocated(self%h))deallocate(self%h);if(allocated(self%rows))deallocate(self%rows)
    if(allocated(self%p))deallocate(self%p);if(allocated(self%x))deallocate(self%x)
    if(allocated(self%dd))deallocate(self%dd);if(allocated(self%du))deallocate(self%du)
    if(allocated(self%xyz))deallocate(self%xyz);if(allocated(self%shift))deallocate(self%shift)
    if(allocated(self%ao_atom))deallocate(self%ao_atom);if(allocated(self%dummy))deallocate(self%dummy)
    nullify(self%functional)
  end subroutine
  subroutine lda_update(self,xce,mythread)
    class(xc_consumer_tdlda_t),intent(inout)::self;class(xc_engine_t),intent(in)::xce;integer::mythread
    real(fp),allocatable::p(:,:),x(:,:),dd(:,:,:),du(:,:,:)
    integer,allocatable::atoms(:);real(fp),allocatable::v(:),f(:),k(:),l(:)
    integer::n,nc,i,status;n=xce%numAOs_p;nc=3*xce%numAtoms
    allocate(p(n,n),x(n,n),atoms(n),v(xce%numPts),f(xce%numPts),k(xce%numPts),l(xce%numPts))
    if(self%response)allocate(dd(n,n,nc),du(n,n,nc))
    if(xce%skip_p)then;p=self%p;x=self%x;atoms=self%ao_atom;if(self%response)then;dd=self%dd;du=self%du;end if
    else;p=self%p(xce%indices_p(1:n),xce%indices_p(1:n));x=self%x(xce%indices_p(1:n),xce%indices_p(1:n))
      atoms=self%ao_atom(xce%indices_p(1:n));if(self%response)then
        dd=self%dd(xce%indices_p(1:n),xce%indices_p(1:n),:);du=self%du(xce%indices_p(1:n),xce%indices_p(1:n),:);end if
    end if
    call lda_kernels(self,xce,v,f,k,l)
    do i=1,xce%numPts
      if(self%response)then
        call response_point(v(i),f(i),k(i),l(i),xce%wfAlpha_p,p,x,dd,du,atoms,xce%aoV(:,i), &
          xce%aoG1(:,i,:),xce%aoG2(:,i,:),xce%xyzw(i,:),xce%currAtom,self,mythread,status)
      else
        call fixed_point(v(i),f(i),k(i),l(i),xce%wfAlpha_p,p,x,atoms,xce%aoV(:,i), &
          xce%aoG1(:,i,:),xce%aoG2(:,i,:),xce%xyzw(i,:),xce%currAtom,self,mythread,status)
      end if
      if(status/=0)error stop 'lda_update: partition derivative failure'
    end do
  end subroutine
  subroutine fixed_point(v,f,k,l,r,p,x,atoms,aov,g1,g2,xyzw,owner,self,thread,status)
    real(fp),intent(in)::v,f,k,l,r(:,:),p(:,:),x(:,:),aov(:),g1(:,:),g2(:,:),xyzw(4)
    integer,intent(in)::atoms(:),owner,thread;class(xc_consumer_tdlda_t),intent(inout)::self;integer,intent(out)::status
    integer::nat;real(fp)::rv,pv,xv,scale;real(fp),allocatable::fr(:,:),fpd(:,:),fx(:,:),sr(:,:,:,:),sp(:,:,:,:),sx(:,:,:,:)
    real(fp),allocatable::dr(:,:),dp(:,:),dx(:,:),d2r(:,:,:,:),d2p(:,:,:,:),d2x(:,:,:,:),w(:),dw(:,:,:),d2w(:,:,:,:,:)
    nat=size(self%xyz,2);allocate(fr(3,nat),fpd(3,nat),fx(3,nat),sr(3,3,nat,nat),sp(3,3,nat,nat),sx(3,3,nat,nat), &
      dr(3,nat),dp(3,nat),dx(3,nat),d2r(3,3,nat,nat),d2p(3,3,nat,nat),d2x(3,3,nat,nat),w(nat),dw(3,nat,nat),d2w(3,nat,3,nat,nat))
    call value(r,aov,rv);if(rv<1e-12_fp)then;status=0;return;end if;call value(p,aov,pv);call value(x,aov,xv)
    call lda_density_nuclear_point(r,atoms,aov,g1,g2,fr,sr);call lda_add_owner_motion(owner,fr,sr,dr,d2r)
    call lda_density_nuclear_point(p,atoms,aov,g1,g2,fpd,sp);call lda_add_owner_motion(owner,fpd,sp,dp,d2p)
    call lda_density_nuclear_point(x,atoms,aov,g1,g2,fx,sx);call lda_add_owner_motion(owner,fx,sx,dx,d2x)
    call partition_weight_nuclear_derivatives(self%xyz,xyzw(1:3),owner,self%dummy,self%part_type,self%has_shift,self%shift,w,dw,d2w,status)
    if(status/=0.or.w(owner)<=sqrt(tiny(1.0_fp)))return;scale=xyzw(4)/w(owner)
    call lda_weighted_fixed_hessian(v,f,k,l,pv,xv,dr,dp,dx,d2r,d2p,d2x,scale,w(owner),dw(:,:,owner),d2w(:,:,:,:,owner),self%h(:,:,:,:,thread))
  end subroutine
  subroutine response_point(v,f,k,l,r,p,x,rd,ru,atoms,aov,g1,g2,xyzw,owner,self,thread,status)
    real(fp),intent(in)::v,f,k,l,r(:,:),p(:,:),x(:,:),rd(:,:,:),ru(:,:,:),aov(:),g1(:,:),g2(:,:),xyzw(4)
    integer,intent(in)::atoms(:),owner,thread;class(xc_consumer_tdlda_t),intent(inout)::self;integer,intent(out)::status
    integer::nat,nc,j;real(fp)::rv,pv,xv,rdv,ruv,scale;real(fp),allocatable::fr(:,:),fpd(:,:),fx(:,:),fd(:,:),fu(:,:),s(:,:,:,:)
    real(fp),allocatable::dr(:,:),dp(:,:),dx(:,:),drd(:,:),dru(:,:),junk(:,:,:,:),w(:),dw(:,:,:),d2w(:,:,:,:,:),row(:,:)
    nat=size(self%xyz,2);nc=3*nat;allocate(fr(3,nat),fpd(3,nat),fx(3,nat),fd(3,nat),fu(3,nat),s(3,3,nat,nat),dr(3,nat),dp(3,nat),dx(3,nat),drd(3,nat),dru(3,nat),junk(3,3,nat,nat),w(nat),dw(3,nat,nat),d2w(3,nat,3,nat,nat),row(3,nat))
    call value(r,aov,rv);if(rv<1e-12_fp)then;status=0;return;end if;call value(p,aov,pv);call value(x,aov,xv)
    call lda_density_nuclear_point(r,atoms,aov,g1,g2,fr,s);call lda_add_owner_motion(owner,fr,s,dr,junk)
    call lda_density_nuclear_point(p,atoms,aov,g1,g2,fpd,s);call lda_add_owner_motion(owner,fpd,s,dp,junk)
    call lda_density_nuclear_point(x,atoms,aov,g1,g2,fx,s);call lda_add_owner_motion(owner,fx,s,dx,junk)
    call partition_weight_nuclear_derivatives(self%xyz,xyzw(1:3),owner,self%dummy,self%part_type,self%has_shift,self%shift,w,dw,d2w,status)
    if(status/=0.or.w(owner)<=sqrt(tiny(1.0_fp)))return;scale=xyzw(4)/w(owner)
    do j=1,nc;call value(rd(:,:,j),aov,rdv);call value(ru(:,:,j),aov,ruv)
      call lda_density_nuclear_point(rd(:,:,j),atoms,aov,g1,g2,fd,s);call lda_add_owner_motion(owner,fd,s,drd,junk)
      call lda_density_nuclear_point(ru(:,:,j),atoms,aov,g1,g2,fu,s);call lda_add_owner_motion(owner,fu,s,dru,junk)
      row=0;call lda_weighted_response_row(f,k,l,pv,xv,rdv,ruv,dr,dp,dx,drd,dru,scale,w(owner),dw(:,:,owner),row)
      self%rows(:,:,j,thread)=self%rows(:,:,j,thread)+row
    end do
  end subroutine
  subroutine value(d,ao,r);real(fp),intent(in)::d(:,:),ao(:);real(fp),intent(out)::r;integer::i,j
    r=0;do j=1,size(ao);do i=1,size(ao);r=r+d(i,j)*ao(i)*ao(j);end do;end do
  end subroutine
  subroutine lda_kernels(self,xce,v,f,k,l)
    class(xc_consumer_tdlda_t)::self;class(xc_engine_t)::xce;real(fp),intent(out)::v(:),f(:),k(:),l(:)
    type(xc_engine_t)::probe;real(fp),allocatable::ones(:);real(fp)::vp,fplus,kp,vm,fm,km,rho,wt;integer::i,n
    real(fp),parameter::eps=1e-4_fp;n=xce%numPts;allocate(probe%xclib);probe%funTyp=OQP_FUNTYP_LDA
    call probe%xclib%init(.false.,.false.,.false.,.false.,2*n,3);call probe%xclib%setPts(2*n)
    do i=1,n;probe%xclib%rho(:,i)=xce%xclib%rho(:,i)*(1+eps);probe%xclib%rho(:,n+i)=xce%xclib%rho(:,i)*(1-eps);end do
    allocate(ones(2*n),source=1.0_fp);call probe%xclib%compute(self%functional,ones)
    do i=1,n;wt=xce%xyzw(i,4);call lda_kernel(xce,i,v(i),f(i),k(i));v(i)=v(i)/wt;f(i)=f(i)/wt;k(i)=k(i)/wt
      call lda_kernel(probe,i,vp,fplus,kp);call lda_kernel(probe,n+i,vm,fm,km);rho=sum(xce%xclib%rho(:,i))
      if(rho<1e-12_fp)then;l(i)=0.0_fp;else;l(i)=(kp-km)/(2*eps*rho);end if
    end do
  end subroutine
  subroutine lda_kernel(xce,i,v,f,k)
    class(xc_engine_t)::xce;integer,intent(in)::i;real(fp),intent(out)::v,f,k
    real(fp),parameter::cr(2)=[.5_fp,.5_fp];real(fp)::zr(2),zs(3),zt(2),dr(2),ds(3),dt(2),fr(2),fs(3),ft(2),ffs(3),gr(2),gs(3),gt(2)
    zr=0;zs=0;zt=0;call xc_der1(xce,.true.,i,dr,ds,dt);v=dot_product(cr,dr)
    call xc_der2_contr(xce,.true.,i,cr,zs,zt,fr,fs,ft);f=dot_product(cr,fr)
    call xc_der3_contr(xce,i,cr,zs,zt,zs,ffs,gr,gs,gt);k=dot_product(cr,gr)
  end subroutine
end module mod_dft_gridint_tdlda_driver
