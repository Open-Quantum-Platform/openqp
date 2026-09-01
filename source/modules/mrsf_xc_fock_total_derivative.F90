! Analytic total nuclear derivative of the spin-resolved semilocal XC Fock
! matrices needed by MRSF-TDDFT and ROKS response equations.
module mrsf_xc_fock_total_derivative_mod

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t,xc_consumer_t,xc_options_t,run_xc, &
    xc_der1,xc_der2_contr
  use mod_dft_gridint_fxc, only: utddft_fxc
  use mod_dft_partition_hessian, only: partition_weight_nuclear_derivatives
  use mod_dft_gridint_mrsf_xc_fock_deriv_point, only: &
    lda_spin_fock_point_derivative,gga_spin_fock_point_derivative, &
    moving_ao_pair_derivative,moving_density_derivative

  implicit none
  private

  type, extends(xc_consumer_t) :: xc_fock_derivative_consumer_t
    real(fp), allocatable :: derivative_a(:,:,:,:)
    real(fp), allocatable :: derivative_b(:,:,:,:)
    real(fp), allocatable :: atom_xyz(:,:),surface_shift(:,:)
    integer, allocatable :: ao_atom(:)
    logical, allocatable :: dummy_atom(:)
    real(fp), allocatable :: worker_error(:)
    integer :: part_fun_type=0
    logical :: has_surface_shift=.false.
    logical :: is_gga=.false.
  contains
    procedure :: parallel_start => derivative_parallel_start
    procedure :: parallel_stop => derivative_parallel_stop
    procedure :: update => derivative_update
    procedure :: postUpdate => derivative_post_update
    procedure :: clean => derivative_clean
  end type xc_fock_derivative_consumer_t

  public :: mrsf_xc_fock_total_derivative
  public :: mrsf_xc_fock_total_derivative_from_dmo

contains

!> Spin-resolved semilocal XC Fock total nuclear derivative.
!>
!> Methodological starting point: Hiroya Nakata's analytical TDHF/TDDFT
!> Hessian formulation.  The first derivative is separated exactly as
!>
!>   dVxc_s = (dVxc_s)_moving-AO/grid + sum_t fxc_st[dP_t].
!>
!> The first term is evaluated by analytic atom-centred grid integrals below;
!> the second reuses the production unrestricted `utddft_fxc` kernel action.
!> Nuclear finite differences and displaced geometries are not used.
!>
!> pa/pb and dpa/dpb are spin density matrices in the AO coefficient
!> representation.  The output coordinate order is x1,y1,z1,x2,y2,z2,... .
!> This density-only interface makes no determinant or Slater construction.
!> Global hybrids are accepted because their semilocal XC part has the same
!> derivative; the exact-exchange derivative is assembled by its own Fock
!> path.  Meta-GGA and range-separated/CAM cases fail closed.
  subroutine mrsf_xc_fock_total_derivative(basis,mol_grid,pa,pb,dpa,dpb, &
      derivative_a,derivative_b,infos,status,threshold)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: pa(:,:),pb(:,:),dpa(:,:,:),dpb(:,:,:)
    real(fp), intent(out) :: derivative_a(:,:,:),derivative_b(:,:,:)
    type(information), target, intent(in) :: infos
    integer, intent(out) :: status
    real(fp), intent(in), optional :: threshold

    type(xc_fock_derivative_consumer_t) :: dat
    type(xc_options_t) :: opts
    real(fp), allocatable, target :: pa_normalized(:,:),pb_normalized(:,:)
    real(fp), allocatable :: dpa_work(:,:,:),dpb_work(:,:,:)
    real(fp), allocatable :: kernel_a(:,:,:),kernel_b(:,:,:)
    real(fp) :: grid_threshold
    integer :: first,i,k,last,natom,nbf,ncart,shell

    nbf=basis%nbf
    natom=infos%mol_prop%natom
    ncart=3*natom
    status=0
    derivative_a=0.0_fp
    derivative_b=0.0_fp
    if(nbf<=0 .or. natom<=0 .or. any(shape(pa)/=[nbf,nbf]) .or. &
       any(shape(pb)/=[nbf,nbf]) .or. any(shape(dpa)/=[nbf,nbf,ncart]) .or. &
       any(shape(dpb)/=[nbf,nbf,ncart]) .or. &
       any(shape(derivative_a)/=[nbf,nbf,ncart]) .or. &
       any(shape(derivative_b)/=[nbf,nbf,ncart])) then
      status=-1
      return
    end if
    if(infos%functional%needTau .or. infos%functional%needLapl) then
      status=-2
      return
    end if
    if(infos%dft%cam_flag) then
      status=-3
      return
    end if

    grid_threshold=infos%dft%grid_density_cutoff
    if(present(threshold)) grid_threshold=threshold
    if(grid_threshold<0.0_fp) then
      status=-4
      return
    end if

    allocate(pa_normalized(nbf,nbf),pb_normalized(nbf,nbf), &
      dpa_work(nbf,nbf,ncart),dpb_work(nbf,nbf,ncart), &
      kernel_a(nbf,nbf,ncart),kernel_b(nbf,nbf,ncart), &
      dat%ao_atom(nbf))
    do i=1,nbf
      pa_normalized(:,i)=pa(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      pb_normalized(:,i)=pb(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
    end do
    dat%ao_atom=0
    do shell=1,basis%nshell
      first=basis%ao_offset(shell)
      last=first+basis%naos(shell)-1
      dat%ao_atom(first:last)=basis%origin(shell)
    end do
    if(any(dat%ao_atom<1) .or. any(dat%ao_atom>natom)) then
      status=-5
      call dat%clean()
      deallocate(pa_normalized,pb_normalized,dpa_work,dpb_work,kernel_a,kernel_b)
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
    ! run_xc adds one derivative for GGA.  Both branches therefore collocate
    ! AO values through G2, as required by d grad(phi_mu phi_nu)/dR.
    opts%nDer=merge(1,2,dat%is_gga)
    opts%nXCDer=2
    opts%numOccAlpha=infos%mol_prop%nelec_A
    opts%numOccBeta=infos%mol_prop%nelec_B
    opts%wfAlpha=>pa_normalized
    opts%wfBeta=>pb_normalized
    opts%dft_threshold=grid_threshold
    ! A geometry-dependent AO-pruning boundary is not admissible in an
    ! analytic nuclear derivative; retain only the density-screening rule.
    opts%ao_threshold=0.0_fp
    opts%ao_sparsity_ratio=0.0_fp
    opts%molGrid=>mol_grid

    call dat%pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    call run_xc(opts,dat,basis)
    if(allocated(dat%worker_error)) then
      if(dat%worker_error(1)>0.0_fp) status=-6
    end if
    if(status==0) then
      do k=1,ncart
        do i=1,nbf
          derivative_a(:,i,k)=dat%derivative_a(:,i,k,1) &
            *basis%bfnrm(:)*basis%bfnrm(i)
          derivative_b(:,i,k)=dat%derivative_b(:,i,k,1) &
            *basis%bfnrm(:)*basis%bfnrm(i)
        end do
      end do
    end if
    call dat%clean()

    if(status==0) then
      ! utddft_fxc scales its density arguments temporarily, so use private
      ! copies even though it restores them before returning.
      dpa_work=dpa
      dpb_work=dpb
      kernel_a=0.0_fp
      kernel_b=0.0_fp
      call utddft_fxc(basis=basis,molGrid=mol_grid,isVecs=.false., &
        wfa=pa,wfb=pb,fxa=kernel_a,fxb=kernel_b,dxa=dpa_work,dxb=dpb_work, &
        nMtx=ncart,threshold=grid_threshold,infos=infos)
      derivative_a=derivative_a+kernel_a
      derivative_b=derivative_b+kernel_b
    end if
    deallocate(pa_normalized,pb_normalized,dpa_work,dpb_work,kernel_a,kernel_b)
  end subroutine mrsf_xc_fock_total_derivative

!-------------------------------------------------------------------------------

!> Orbital-response convenience interface.  Occupation vectors permit ROKS
!> and other spin-resolved ensembles without introducing configuration-state
!> or determinant objects.
  subroutine mrsf_xc_fock_total_derivative_from_dmo(basis,mol_grid,pa,pb, &
      mo_a,mo_b,dmo_a,dmo_b,occupation_a,occupation_b,derivative_a, &
      derivative_b,infos,status,threshold)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: pa(:,:),pb(:,:),mo_a(:,:),mo_b(:,:)
    real(fp), intent(in) :: dmo_a(:,:,:),dmo_b(:,:,:)
    real(fp), intent(in) :: occupation_a(:),occupation_b(:)
    real(fp), intent(out) :: derivative_a(:,:,:),derivative_b(:,:,:)
    type(information), target, intent(in) :: infos
    integer, intent(out) :: status
    real(fp), intent(in), optional :: threshold

    real(fp), allocatable :: dpa(:,:,:),dpb(:,:,:)
    integer :: i,j,k,ncart,nbf,orbital

    nbf=size(pa,1)
    ncart=size(dmo_a,3)
    status=0
    derivative_a=0.0_fp
    derivative_b=0.0_fp
    if(nbf<=0 .or. size(pa,2)/=nbf .or. any(shape(pb)/=[nbf,nbf]) .or. &
       size(mo_a,1)/=nbf .or. size(mo_b,1)/=nbf .or. &
       size(dmo_a,1)/=nbf .or. size(dmo_b,1)/=nbf .or. &
       size(dmo_a,2)/=size(mo_a,2) .or. size(dmo_b,2)/=size(mo_b,2) .or. &
       size(dmo_b,3)/=ncart .or. size(occupation_a)/=size(mo_a,2) .or. &
       size(occupation_b)/=size(mo_b,2)) then
      status=-1
      return
    end if
    allocate(dpa(nbf,nbf,ncart),dpb(nbf,nbf,ncart),source=0.0_fp)
    do k=1,ncart
      do orbital=1,size(mo_a,2)
        if(abs(occupation_a(orbital))<=tiny(1.0_fp)) cycle
        do j=1,nbf
          do i=1,nbf
            dpa(i,j,k)=dpa(i,j,k)+occupation_a(orbital)*( &
              dmo_a(i,orbital,k)*mo_a(j,orbital) &
              +mo_a(i,orbital)*dmo_a(j,orbital,k))
          end do
        end do
      end do
      do orbital=1,size(mo_b,2)
        if(abs(occupation_b(orbital))<=tiny(1.0_fp)) cycle
        do j=1,nbf
          do i=1,nbf
            dpb(i,j,k)=dpb(i,j,k)+occupation_b(orbital)*( &
              dmo_b(i,orbital,k)*mo_b(j,orbital) &
              +mo_b(i,orbital)*dmo_b(j,orbital,k))
          end do
        end do
      end do
    end do
    call mrsf_xc_fock_total_derivative(basis,mol_grid,pa,pb,dpa,dpb, &
      derivative_a,derivative_b,infos,status,threshold)
    deallocate(dpa,dpb)
  end subroutine mrsf_xc_fock_total_derivative_from_dmo

!-------------------------------------------------------------------------------

  subroutine derivative_parallel_start(self,xce,nthreads)
    class(xc_fock_derivative_consumer_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: ncart
    ncart=3*xce%numAtoms
    allocate(self%derivative_a(xce%numAOs,xce%numAOs,ncart,nthreads), &
      self%derivative_b(xce%numAOs,xce%numAOs,ncart,nthreads), &
      self%worker_error(nthreads),source=0.0_fp)
  end subroutine derivative_parallel_start

  subroutine derivative_parallel_stop(self)
    class(xc_fock_derivative_consumer_t), intent(inout) :: self
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
  end subroutine derivative_parallel_stop

  subroutine derivative_post_update(self,xce,mythread)
    class(xc_fock_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
    ! Contributions are accumulated directly by update.
  end subroutine derivative_post_update

  subroutine derivative_clean(self)
    class(xc_fock_derivative_consumer_t), intent(inout) :: self
    if(allocated(self%derivative_a)) deallocate(self%derivative_a)
    if(allocated(self%derivative_b)) deallocate(self%derivative_b)
    if(allocated(self%atom_xyz)) deallocate(self%atom_xyz)
    if(allocated(self%surface_shift)) deallocate(self%surface_shift)
    if(allocated(self%ao_atom)) deallocate(self%ao_atom)
    if(allocated(self%dummy_atom)) deallocate(self%dummy_atom)
    if(allocated(self%worker_error)) deallocate(self%worker_error)
  end subroutine derivative_clean

!-------------------------------------------------------------------------------

  subroutine derivative_update(self,xce,mythread)
    class(xc_fock_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(fp), allocatable :: pa(:,:),pb(:,:),drho(:,:),dgrad_rho(:,:,:)
    real(fp), allocatable :: dweight_flat(:),dvr(:,:),dvs(:,:)
    real(fp), allocatable :: dpair(:),dgrad_pair(:,:),point_derivative(:,:)
    real(fp), allocatable :: weights(:),dweights(:,:,:),d2weights(:,:,:,:,:)
    integer, allocatable :: atoms(:)
    real(fp) :: dr(2),ds(3),dt(2),fr(2),fs(3),ft(2)
    real(fp) :: vr(2),vs(3),vt(2),grad_rho(3,2)
    real(fp) :: pair,grad_pair(3),quadrature_scale,finite_weight
    integer :: atom,cart,coordinate,global_mu,global_nu,ipt,mu,n,nat,ncart
    integer :: nu,owner,pair_status,partition_status,spin_status

    n=xce%numAOs_p
    nat=xce%numAtoms
    ncart=3*nat
    owner=xce%currAtom
    if(owner<1 .or. owner>nat) then
      self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
      return
    end if
    allocate(pa(n,n),pb(n,n),atoms(n),drho(2,ncart), &
      dgrad_rho(3,2,ncart),dweight_flat(ncart),dvr(2,ncart), &
      dvs(3,ncart),dpair(ncart),dgrad_pair(3,ncart), &
      point_derivative(2,ncart),weights(nat),dweights(3,nat,nat), &
      d2weights(3,nat,3,nat,nat))
    if(xce%skip_p) then
      pa=xce%wfAlpha_p
      pb=xce%wfBeta_p
      atoms=self%ao_atom
    else
      pa=xce%wfAlpha_p
      pb=xce%wfBeta_p
      atoms=self%ao_atom(xce%indices_p(1:n))
    end if

    do ipt=1,xce%numPts
      finite_weight=xce%xyzw(ipt,4)
      if(abs(finite_weight)<=tiny(1.0_fp)) cycle
      call partition_weight_nuclear_derivatives(self%atom_xyz, &
        xce%xyzw(ipt,1:3),owner,self%dummy_atom,self%part_fun_type, &
        self%has_surface_shift,self%surface_shift,weights,dweights, &
        d2weights,partition_status)
      if(partition_status/=0 .or. &
         weights(owner)<=sqrt(tiny(1.0_fp))) then
        self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
        return
      end if
      quadrature_scale=finite_weight/weights(owner)
      do atom=1,nat
        do cart=1,3
          coordinate=3*(atom-1)+cart
          dweight_flat(coordinate)=dweights(cart,atom,owner)
        end do
      end do

      call moving_density_derivative(pa,atoms,owner,xce%aoV(:,ipt), &
        xce%aoG1(:,ipt,:),xce%aoG2(:,ipt,:),self%is_gga,drho(1,:), &
        dgrad_rho(:,1,:),spin_status)
      if(spin_status/=0) then
        self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
        return
      end if
      call moving_density_derivative(pb,atoms,owner,xce%aoV(:,ipt), &
        xce%aoG1(:,ipt,:),xce%aoG2(:,ipt,:),self%is_gga,drho(2,:), &
        dgrad_rho(:,2,:),spin_status)
      if(spin_status/=0) then
        self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
        return
      end if

      call xc_der1(xce,.true.,ipt,vr,vs,vt)
      vr=vr/finite_weight
      vs=vs/finite_weight
      grad_rho(:,1)=xce%xclib%drho(1:3,ipt)
      grad_rho(:,2)=xce%xclib%drho(4:6,ipt)
      do coordinate=1,ncart
        dr=drho(:,coordinate)
        ds=0.0_fp
        if(self%is_gga) then
          ds(1)=2.0_fp*dot_product(grad_rho(:,1), &
            dgrad_rho(:,1,coordinate))
          ds(2)=2.0_fp*dot_product(grad_rho(:,2), &
            dgrad_rho(:,2,coordinate))
          ds(3)=dot_product(dgrad_rho(:,1,coordinate),grad_rho(:,2)) &
            +dot_product(grad_rho(:,1),dgrad_rho(:,2,coordinate))
        end if
        dt=0.0_fp
        call xc_der2_contr(xce,.true.,ipt,dr,ds,dt,fr,fs,ft)
        dvr(:,coordinate)=fr/finite_weight
        dvs(:,coordinate)=fs/finite_weight
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
            pair,grad_pair,dpair,dgrad_pair,pair_status)
          if(pair_status/=0) then
            self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
            return
          end if
          if(self%is_gga) then
            call gga_spin_fock_point_derivative(quadrature_scale, &
              weights(owner),dweight_flat,vr,vs,dvr,dvs,grad_rho, &
              dgrad_rho,pair,grad_pair,dpair,dgrad_pair,point_derivative, &
              spin_status)
          else
            call lda_spin_fock_point_derivative(quadrature_scale, &
              weights(owner),dweight_flat,vr,dvr,pair,dpair, &
              point_derivative,spin_status)
          end if
          if(spin_status/=0) then
            self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
            return
          end if
          self%derivative_a(global_mu,global_nu,:,mythread)= &
            self%derivative_a(global_mu,global_nu,:,mythread) &
            +point_derivative(1,:)
          self%derivative_b(global_mu,global_nu,:,mythread)= &
            self%derivative_b(global_mu,global_nu,:,mythread) &
            +point_derivative(2,:)
        end do
      end do
    end do
  end subroutine derivative_update

end module mrsf_xc_fock_total_derivative_mod
