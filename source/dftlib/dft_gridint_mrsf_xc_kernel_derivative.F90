! Analytic nuclear derivative of the spin-polarized semilocal XC-kernel
! action needed by the differentiated MRSF orbital-adjoint equations.
module mod_dft_gridint_mrsf_xc_kernel_derivative

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t,xc_consumer_t,xc_options_t,run_xc, &
    xc_der1,xc_der2_contr,xc_der3_contr
  use mod_dft_partition_hessian, only: partition_weight_nuclear_derivatives
  use mod_dft_gridint_mrsf_xc_fock_deriv_point, only: &
    moving_ao_pair_derivative,moving_density_derivative, &
    spin_fock_point_derivative
  use mod_dft_gridint_mrsf_xc_hessian_point, only: &
    mrsf_xc_kernel_fock_coefficients

  implicit none
  private

  type, extends(xc_consumer_t) :: kernel_derivative_consumer_t
    real(fp), allocatable :: derivative_a(:,:,:,:),derivative_b(:,:,:,:)
    real(fp), allocatable :: probe_a(:,:),probe_b(:,:)
    real(fp), allocatable :: dground_a(:,:,:),dground_b(:,:,:)
    real(fp), allocatable :: dprobe_a(:,:,:),dprobe_b(:,:,:)
    real(fp), allocatable :: atom_xyz(:,:),surface_shift(:,:)
    integer, allocatable :: ao_atom(:)
    logical, allocatable :: dummy_atom(:)
    real(fp), allocatable :: worker_error(:)
    integer :: part_fun_type=0
    logical :: has_surface_shift=.false.,is_gga=.false.
  contains
    procedure :: parallel_start => kernel_start
    procedure :: parallel_stop => kernel_stop
    procedure :: update => kernel_update
    procedure :: postUpdate => kernel_post
    procedure :: clean => kernel_clean
  end type kernel_derivative_consumer_t

  public :: mrsf_xc_kernel_fock_total_derivative

contains

!> Differentiate K_xc[D](Q) for physical alpha/beta reference and probe
!> densities D and Q.  dD and dQ are their AO-coefficient responses.  Moving
!> AO, atom-centred grid, partition-weight, fxc, and kxc terms are included.
!> No displaced geometry or electronic-state expansion is used.
  subroutine mrsf_xc_kernel_fock_total_derivative(basis,mol_grid,da,db, &
      qa,qb,dda,ddb,dqa,dqb,derivative_a,derivative_b,infos,status,threshold)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set), intent(in) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: da(:,:),db(:,:),qa(:,:),qb(:,:)
    real(fp), intent(in) :: dda(:,:,:),ddb(:,:,:),dqa(:,:,:),dqb(:,:,:)
    real(fp), intent(out) :: derivative_a(:,:,:),derivative_b(:,:,:)
    type(information), target, intent(in) :: infos
    integer, intent(out) :: status
    real(fp), intent(in), optional :: threshold

    type(kernel_derivative_consumer_t) :: dat
    type(xc_options_t) :: opts
    real(fp), allocatable, target :: da_normalized(:,:),db_normalized(:,:)
    real(fp) :: grid_threshold
    integer :: first,i,k,last,natom,nbf,ncart,shell

    nbf=basis%nbf
    natom=infos%mol_prop%natom
    ncart=3*natom
    status=0
    derivative_a=0.0_fp
    derivative_b=0.0_fp
    if(nbf<=0 .or. natom<=0 .or. any(shape(da)/=[nbf,nbf]) .or. &
       any(shape(db)/=[nbf,nbf]) .or. any(shape(qa)/=[nbf,nbf]) .or. &
       any(shape(qb)/=[nbf,nbf]) .or. any(shape(dda)/=[nbf,nbf,ncart]) .or. &
       any(shape(ddb)/=[nbf,nbf,ncart]) .or. &
       any(shape(dqa)/=[nbf,nbf,ncart]) .or. &
       any(shape(dqb)/=[nbf,nbf,ncart]) .or. &
       any(shape(derivative_a)/=[nbf,nbf,ncart]) .or. &
       any(shape(derivative_b)/=[nbf,nbf,ncart])) then
      status=-1
      return
    end if
    if(infos%functional%needTau .or. infos%functional%needLapl .or. &
       infos%dft%cam_flag) then
      status=-2
      return
    end if
    grid_threshold=infos%dft%grid_density_cutoff
    if(present(threshold)) grid_threshold=threshold
    if(grid_threshold<0.0_fp) then
      status=-3
      return
    end if

    allocate(da_normalized(nbf,nbf),db_normalized(nbf,nbf), &
      dat%probe_a(nbf,nbf),dat%probe_b(nbf,nbf), &
      dat%dground_a(nbf,nbf,ncart),dat%dground_b(nbf,nbf,ncart), &
      dat%dprobe_a(nbf,nbf,ncart),dat%dprobe_b(nbf,nbf,ncart), &
      dat%ao_atom(nbf))
    do i=1,nbf
      da_normalized(:,i)=da(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      db_normalized(:,i)=db(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      dat%probe_a(:,i)=qa(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      dat%probe_b(:,i)=qb(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      do k=1,ncart
        dat%dground_a(:,i,k)=dda(:,i,k)*basis%bfnrm(:)*basis%bfnrm(i)
        dat%dground_b(:,i,k)=ddb(:,i,k)*basis%bfnrm(:)*basis%bfnrm(i)
        dat%dprobe_a(:,i,k)=dqa(:,i,k)*basis%bfnrm(:)*basis%bfnrm(i)
        dat%dprobe_b(:,i,k)=dqb(:,i,k)*basis%bfnrm(:)*basis%bfnrm(i)
      end do
    end do
    dat%ao_atom=0
    do shell=1,basis%nshell
      first=basis%ao_offset(shell)
      last=first+basis%naos(shell)-1
      dat%ao_atom(first:last)=basis%origin(shell)
    end do
    if(any(dat%ao_atom<1) .or. any(dat%ao_atom>natom)) then
      status=-4
      call dat%clean()
      deallocate(da_normalized,db_normalized)
      return
    end if
    dat%atom_xyz=infos%atoms%xyz
    dat%dummy_atom=mol_grid%dummyAtom
    dat%surface_shift=mol_grid%surfaceShift
    dat%part_fun_type=mol_grid%partFunType
    dat%has_surface_shift=mol_grid%hasSurfaceShift
    dat%is_gga=infos%functional%needGrd

    opts%isGGA=dat%is_gga
    opts%needTau=.false.
    opts%functional=>infos%functional
    opts%hasBeta=.true.
    opts%isWFVecs=.false.
    opts%numAOs=nbf
    opts%maxPts=mol_grid%maxSlicePts
    opts%limPts=mol_grid%maxNRadTimesNAng
    opts%numAtoms=natom
    opts%maxAngMom=basis%mxam
    opts%nDer=merge(1,2,dat%is_gga)
    opts%nXCDer=3
    opts%numOccAlpha=infos%mol_prop%nelec_A
    opts%numOccBeta=infos%mol_prop%nelec_B
    opts%wfAlpha=>da_normalized
    opts%wfBeta=>db_normalized
    opts%dft_threshold=grid_threshold
    opts%ao_threshold=0.0_fp
    opts%ao_sparsity_ratio=0.0_fp
    opts%molGrid=>mol_grid

    call dat%pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    call run_xc(opts,dat,basis)
    if(allocated(dat%worker_error)) then
      if(dat%worker_error(1)>0.0_fp) status=-5
    end if
    if(status==0) then
      do k=1,ncart
        do i=1,nbf
          derivative_a(:,i,k)=dat%derivative_a(:,i,k,1)* &
            basis%bfnrm(:)*basis%bfnrm(i)
          derivative_b(:,i,k)=dat%derivative_b(:,i,k,1)* &
            basis%bfnrm(:)*basis%bfnrm(i)
        end do
      end do
    end if
    call dat%clean()
    deallocate(da_normalized,db_normalized)
  end subroutine mrsf_xc_kernel_fock_total_derivative

!-----------------------------------------------------------------------------

  subroutine kernel_start(self,xce,nthreads)
    class(kernel_derivative_consumer_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: ncart
    ncart=3*xce%numAtoms
    allocate(self%derivative_a(xce%numAOs,xce%numAOs,ncart,nthreads), &
      self%derivative_b(xce%numAOs,xce%numAOs,ncart,nthreads), &
      self%worker_error(nthreads),source=0.0_fp)
  end subroutine kernel_start

  subroutine kernel_stop(self)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    if(size(self%derivative_a,4)>1) then
      self%derivative_a(:,:,:,1)=sum(self%derivative_a,dim=4)
      self%derivative_b(:,:,:,1)=sum(self%derivative_b,dim=4)
      self%worker_error(1)=sum(self%worker_error)
    end if
    call self%pe%allreduce(self%derivative_a(:,:,:,1), &
      size(self%derivative_a(:,:,:,1)))
    call self%pe%allreduce(self%derivative_b(:,:,:,1), &
      size(self%derivative_b(:,:,:,1)))
    call self%pe%allreduce(self%worker_error(1),1)
  end subroutine kernel_stop

  subroutine kernel_post(self,xce,mythread)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
  end subroutine kernel_post

  subroutine kernel_clean(self)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    if(allocated(self%derivative_a)) deallocate(self%derivative_a)
    if(allocated(self%derivative_b)) deallocate(self%derivative_b)
    if(allocated(self%probe_a)) deallocate(self%probe_a,self%probe_b)
    if(allocated(self%dground_a)) &
      deallocate(self%dground_a,self%dground_b,self%dprobe_a,self%dprobe_b)
    if(allocated(self%atom_xyz)) deallocate(self%atom_xyz)
    if(allocated(self%surface_shift)) deallocate(self%surface_shift)
    if(allocated(self%ao_atom)) deallocate(self%ao_atom)
    if(allocated(self%dummy_atom)) deallocate(self%dummy_atom)
    if(allocated(self%worker_error)) deallocate(self%worker_error)
  end subroutine kernel_clean

!-----------------------------------------------------------------------------

  subroutine kernel_update(self,xce,mythread)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(fp), allocatable :: ga(:,:),gb(:,:),qa(:,:),qb(:,:)
    real(fp), allocatable :: dga(:,:,:),dgb(:,:,:),dqa(:,:,:),dqb(:,:,:)
    integer, allocatable :: atoms(:)
    integer :: n,ipt,status

    n=xce%numAOs_p
    allocate(ga(n,n),gb(n,n),qa(n,n),qb(n,n), &
      dga(n,n,3*xce%numAtoms),dgb(n,n,3*xce%numAtoms), &
      dqa(n,n,3*xce%numAtoms),dqb(n,n,3*xce%numAtoms),atoms(n))
    ga=xce%wfAlpha_p
    gb=xce%wfBeta_p
    if(xce%skip_p) then
      qa=self%probe_a
      qb=self%probe_b
      dga=self%dground_a
      dgb=self%dground_b
      dqa=self%dprobe_a
      dqb=self%dprobe_b
      atoms=self%ao_atom
    else
      qa=self%probe_a(xce%indices_p(1:n),xce%indices_p(1:n))
      qb=self%probe_b(xce%indices_p(1:n),xce%indices_p(1:n))
      dga=self%dground_a(xce%indices_p(1:n),xce%indices_p(1:n),:)
      dgb=self%dground_b(xce%indices_p(1:n),xce%indices_p(1:n),:)
      dqa=self%dprobe_a(xce%indices_p(1:n),xce%indices_p(1:n),:)
      dqb=self%dprobe_b(xce%indices_p(1:n),xce%indices_p(1:n),:)
      atoms=self%ao_atom(xce%indices_p(1:n))
    end if
    do ipt=1,xce%numPts
      call accumulate_kernel_point(self,xce,mythread,ipt,ga,gb,qa,qb, &
        dga,dgb,dqa,dqb,atoms,status)
      if(status/=0) then
        self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
        exit
      end if
    end do
    deallocate(ga,gb,qa,qb,dga,dgb,dqa,dqb,atoms)
  end subroutine kernel_update

!-----------------------------------------------------------------------------

  subroutine accumulate_kernel_point(self,xce,mythread,ipt,ga,gb,qa,qb, &
      dga,dgb,dqa,dqb,atoms,status)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: mythread,ipt,atoms(:)
    real(fp), intent(in) :: ga(:,:),gb(:,:),qa(:,:),qb(:,:)
    real(fp), intent(in) :: dga(:,:,:),dgb(:,:,:),dqa(:,:,:),dqb(:,:,:)
    integer, intent(out) :: status

    real(fp), allocatable :: dg_rho(:,:),dp_rho(:,:),dg_grad(:,:,:)
    real(fp), allocatable :: dp_grad(:,:,:),weights(:),dweights(:,:,:)
    real(fp), allocatable :: d2weights(:,:,:,:,:),ground(:),probe(:)
    real(fp), allocatable :: dground(:,:),dprobe(:,:),first(:),second(:,:)
    real(fp), allocatable :: third(:,:,:),v_r(:),dv_r(:,:)
    real(fp), allocatable :: dweight_flat(:),dpair(:),dgrad_pair(:,:)
    real(fp), allocatable :: point_derivative(:,:),coefficient(:,:)
    real(fp), allocatable :: dcoefficient(:,:,:)
    real(fp) :: rho(2),prho(2),grad_rho(3,2),grad_probe(3,2)
    real(fp) :: pair,grad_pair(3),finite_weight,quadrature_scale
    integer :: atom,cart,coordinate,global_mu,global_nu,i,j,k,local_status
    integer :: mu,n,nat,ncart,nvar,nu,owner,spin

    n=size(ga,1)
    nat=xce%numAtoms
    ncart=3*nat
    nvar=merge(5,2,self%is_gga)
    owner=xce%currAtom
    status=0
    finite_weight=xce%xyzw(ipt,4)
    if(abs(finite_weight)<=tiny(1.0_fp)) return
    allocate(dg_rho(2,ncart),dp_rho(2,ncart),dg_grad(3,2,ncart), &
      dp_grad(3,2,ncart),weights(nat),dweights(3,nat,nat), &
      d2weights(3,nat,3,nat,nat),ground(nvar),probe(nvar), &
      dground(nvar,ncart),dprobe(nvar,ncart),first(nvar), &
      second(nvar,nvar),third(nvar,nvar,nvar),v_r(2),dv_r(2,ncart), &
      dweight_flat(ncart), &
      dpair(ncart),dgrad_pair(3,ncart),point_derivative(2,ncart), &
      coefficient(3,2),dcoefficient(3,2,ncart))

    call field_value_gradient(ga,xce%aoV(:,ipt),xce%aoG1(:,ipt,:), &
      rho(1),grad_rho(:,1))
    call field_value_gradient(gb,xce%aoV(:,ipt),xce%aoG1(:,ipt,:), &
      rho(2),grad_rho(:,2))
    call field_value_gradient(qa,xce%aoV(:,ipt),xce%aoG1(:,ipt,:), &
      prho(1),grad_probe(:,1))
    call field_value_gradient(qb,xce%aoV(:,ipt),xce%aoG1(:,ipt,:), &
      prho(2),grad_probe(:,2))
    call total_density_derivative(ga,dga,atoms,owner,xce,ipt, &
      dg_rho(1,:),dg_grad(:,1,:),local_status)
    if(local_status==0) call total_density_derivative(gb,dgb,atoms,owner, &
      xce,ipt,dg_rho(2,:),dg_grad(:,2,:),local_status)
    if(local_status==0) call total_density_derivative(qa,dqa,atoms,owner, &
      xce,ipt,dp_rho(1,:),dp_grad(:,1,:),local_status)
    if(local_status==0) call total_density_derivative(qb,dqb,atoms,owner, &
      xce,ipt,dp_rho(2,:),dp_grad(:,2,:),local_status)
    if(local_status/=0) then
      status=-1
      call cleanup()
      return
    end if
    call build_density_variables(self%is_gga,rho,grad_rho,prho,grad_probe, &
      dg_rho,dg_grad,dp_rho,dp_grad,ground,probe,dground,dprobe)
    call build_unweighted_kernels(xce,ipt,nvar,finite_weight,first,second, &
      third)
    call mrsf_xc_kernel_fock_coefficients(first,second,third,probe,dground, &
      dprobe,self%is_gga,grad_rho,grad_probe,dg_grad,dp_grad,v_r, &
      coefficient,dv_r,dcoefficient,local_status)
    if(local_status/=0) then
      status=-5
      call cleanup()
      return
    end if

    call partition_weight_nuclear_derivatives(self%atom_xyz, &
      xce%xyzw(ipt,1:3),owner,self%dummy_atom,self%part_fun_type, &
      self%has_surface_shift,self%surface_shift,weights,dweights,d2weights, &
      local_status)
    if(local_status/=0 .or. weights(owner)<=sqrt(tiny(1.0_fp))) then
      status=-2
      call cleanup()
      return
    end if
    quadrature_scale=finite_weight/weights(owner)
    do atom=1,nat
      do cart=1,3
        coordinate=3*(atom-1)+cart
        dweight_flat(coordinate)=dweights(cart,atom,owner)
      end do
    end do

    do nu=1,n
      global_nu=nu
      if(.not.xce%skip_p) global_nu=xce%indices_p(nu)
      do mu=1,n
        global_mu=mu
        if(.not.xce%skip_p) global_mu=xce%indices_p(mu)
        call moving_ao_pair_derivative(owner,atoms(mu),atoms(nu), &
          xce%aoV(mu,ipt),xce%aoV(nu,ipt),xce%aoG1(mu,ipt,:), &
          xce%aoG1(nu,ipt,:),xce%aoG2(mu,ipt,:),xce%aoG2(nu,ipt,:), &
          pair,grad_pair,dpair,dgrad_pair,local_status)
        if(local_status/=0) then
          status=-3
          call cleanup()
          return
        end if
        call spin_fock_point_derivative(quadrature_scale,weights(owner), &
          dweight_flat,v_r,coefficient,dv_r, &
          dcoefficient,pair,grad_pair,dpair,dgrad_pair,point_derivative, &
          local_status)
        if(local_status/=0) then
          status=-4
          call cleanup()
          return
        end if
        do coordinate=1,ncart
          self%derivative_a(global_mu,global_nu,coordinate,mythread)= &
            self%derivative_a(global_mu,global_nu,coordinate,mythread)+ &
            point_derivative(1,coordinate)
          self%derivative_b(global_mu,global_nu,coordinate,mythread)= &
            self%derivative_b(global_mu,global_nu,coordinate,mythread)+ &
            point_derivative(2,coordinate)
        end do
      end do
    end do
    call cleanup()

  contains

    subroutine cleanup()
      if(allocated(dg_rho)) deallocate(dg_rho,dp_rho,dg_grad,dp_grad)
      if(allocated(weights)) deallocate(weights,dweights,d2weights)
      if(allocated(ground)) deallocate(ground,probe,dground,dprobe)
      if(allocated(first)) deallocate(first,second,third,v_r,dv_r)
      if(allocated(dweight_flat)) deallocate(dweight_flat,dpair,dgrad_pair, &
        point_derivative,coefficient,dcoefficient)
    end subroutine cleanup

  end subroutine accumulate_kernel_point

!-----------------------------------------------------------------------------

  subroutine total_density_derivative(density,response,atoms,owner,xce,ipt, &
      drho,dgrad,status)
    real(fp), intent(in) :: density(:,:),response(:,:,:)
    integer, intent(in) :: atoms(:),owner,ipt
    class(xc_engine_t), intent(in) :: xce
    real(fp), intent(out) :: drho(:),dgrad(:,:)
    integer, intent(out) :: status
    real(fp) :: value,gradient(3)
    integer :: coordinate

    call moving_density_derivative(density,atoms,owner,xce%aoV(:,ipt), &
      xce%aoG1(:,ipt,:),xce%aoG2(:,ipt,:),.true.,drho,dgrad,status)
    if(status/=0) return
    do coordinate=1,size(drho)
      call field_value_gradient(response(:,:,coordinate),xce%aoV(:,ipt), &
        xce%aoG1(:,ipt,:),value,gradient)
      drho(coordinate)=drho(coordinate)+value
      dgrad(:,coordinate)=dgrad(:,coordinate)+gradient
    end do
  end subroutine total_density_derivative

  pure subroutine field_value_gradient(density,aov,aog1,value,gradient)
    real(fp), intent(in) :: density(:,:),aov(:),aog1(:,:)
    real(fp), intent(out) :: value,gradient(3)
    integer :: mu,nu
    value=0.0_fp
    gradient=0.0_fp
    do nu=1,size(aov)
      do mu=1,size(aov)
        value=value+density(mu,nu)*aov(mu)*aov(nu)
        gradient=gradient+density(mu,nu)*( &
          aog1(mu,:)*aov(nu)+aov(mu)*aog1(nu,:))
      end do
    end do
  end subroutine field_value_gradient

  pure subroutine build_density_variables(is_gga,rho,grho,prho,pgrho, &
      drho,dgrho,dprho,dpgrho,ground,probe,dground,dprobe)
    logical, intent(in) :: is_gga
    real(fp), intent(in) :: rho(2),grho(3,2),prho(2),pgrho(3,2)
    real(fp), intent(in) :: drho(:,:),dgrho(:,:,:),dprho(:,:),dpgrho(:,:,:)
    real(fp), intent(out) :: ground(:),probe(:),dground(:,:),dprobe(:,:)
    integer :: coordinate

    ground=0.0_fp
    probe=0.0_fp
    dground=0.0_fp
    dprobe=0.0_fp
    ground(1:2)=rho
    probe(1:2)=prho
    dground(1:2,:)=drho
    dprobe(1:2,:)=dprho
    if(.not.is_gga) return
    ground(3)=dot_product(grho(:,1),grho(:,1))
    ground(4)=dot_product(grho(:,2),grho(:,2))
    ground(5)=dot_product(grho(:,1),grho(:,2))
    probe(3)=2.0_fp*dot_product(grho(:,1),pgrho(:,1))
    probe(4)=2.0_fp*dot_product(grho(:,2),pgrho(:,2))
    probe(5)=dot_product(grho(:,1),pgrho(:,2))+ &
      dot_product(grho(:,2),pgrho(:,1))
    do coordinate=1,size(drho,2)
      dground(3,coordinate)=2.0_fp*dot_product(grho(:,1), &
        dgrho(:,1,coordinate))
      dground(4,coordinate)=2.0_fp*dot_product(grho(:,2), &
        dgrho(:,2,coordinate))
      dground(5,coordinate)=dot_product(dgrho(:,1,coordinate),grho(:,2))+ &
        dot_product(grho(:,1),dgrho(:,2,coordinate))
      dprobe(3,coordinate)=2.0_fp*( &
        dot_product(dgrho(:,1,coordinate),pgrho(:,1))+ &
        dot_product(grho(:,1),dpgrho(:,1,coordinate)))
      dprobe(4,coordinate)=2.0_fp*( &
        dot_product(dgrho(:,2,coordinate),pgrho(:,2))+ &
        dot_product(grho(:,2),dpgrho(:,2,coordinate)))
      dprobe(5,coordinate)= &
        dot_product(dgrho(:,1,coordinate),pgrho(:,2))+ &
        dot_product(grho(:,1),dpgrho(:,2,coordinate))+ &
        dot_product(dgrho(:,2,coordinate),pgrho(:,1))+ &
        dot_product(grho(:,2),dpgrho(:,1,coordinate))
    end do
  end subroutine build_density_variables

  subroutine build_unweighted_kernels(xce,ipt,nvar,finite_weight,first, &
      second,third)
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: ipt,nvar
    real(fp), intent(in) :: finite_weight
    real(fp), intent(out) :: first(:),second(:,:),third(:,:,:)
    real(fp) :: dr(2),ds(3),dt(2),fr(2),fs(3),ft(2),ffs(3)
    real(fp) :: gr(2),gs(3),gt(2)
    real(fp), allocatable :: square(:,:),mixed(:)
    integer :: j,k

    call xc_der1(xce,.true.,ipt,dr,ds,dt)
    call join_direction(nvar,dr,ds,first)
    first=first/finite_weight
    second=0.0_fp
    do j=1,nvar
      call split_basis(nvar,j,dr,ds)
      dt=0.0_fp
      call xc_der2_contr(xce,.true.,ipt,dr,ds,dt,fr,fs,ft)
      call join_direction(nvar,fr,fs,second(:,j))
    end do
    second=0.5_fp*(second+transpose(second))/finite_weight
    allocate(square(nvar,nvar),mixed(nvar))
    do j=1,nvar
      call split_basis(nvar,j,dr,ds)
      dt=0.0_fp
      call xc_der3_contr(xce,ipt,dr,ds,dt,[0.0_fp,0.0_fp,0.0_fp], &
        ffs,gr,gs,gt)
      call join_direction(nvar,gr,gs,square(:,j))
    end do
    third=0.0_fp
    do j=1,nvar
      third(:,j,j)=square(:,j)/finite_weight
      do k=j+1,nvar
        call split_pair(nvar,j,k,dr,ds)
        dt=0.0_fp
        call xc_der3_contr(xce,ipt,dr,ds,dt,[0.0_fp,0.0_fp,0.0_fp], &
          ffs,gr,gs,gt)
        call join_direction(nvar,gr,gs,mixed)
        mixed=0.5_fp*(mixed-square(:,j)-square(:,k))/finite_weight
        third(:,j,k)=mixed
        third(:,k,j)=mixed
      end do
    end do
    deallocate(square,mixed)
  end subroutine build_unweighted_kernels

  pure subroutine split_basis(nvar,index,dr,ds)
    integer, intent(in) :: nvar,index
    real(fp), intent(out) :: dr(2),ds(3)
    dr=0.0_fp
    ds=0.0_fp
    if(index<=2) then
      dr(index)=1.0_fp
    else if(nvar==5) then
      ds(index-2)=1.0_fp
    end if
  end subroutine split_basis

  pure subroutine split_pair(nvar,index_a,index_b,dr,ds)
    integer, intent(in) :: nvar,index_a,index_b
    real(fp), intent(out) :: dr(2),ds(3)
    real(fp) :: dra(2),drb(2),dsa(3),dsb(3)
    call split_basis(nvar,index_a,dra,dsa)
    call split_basis(nvar,index_b,drb,dsb)
    dr=dra+drb
    ds=dsa+dsb
  end subroutine split_pair

  pure subroutine join_direction(nvar,dr,ds,result)
    integer, intent(in) :: nvar
    real(fp), intent(in) :: dr(2),ds(3)
    real(fp), intent(out) :: result(:)
    result=0.0_fp
    result(1:2)=dr
    if(nvar==5) result(3:5)=ds
  end subroutine join_direction

end module mod_dft_gridint_mrsf_xc_kernel_derivative
