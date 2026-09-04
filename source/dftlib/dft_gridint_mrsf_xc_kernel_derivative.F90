! Analytic nuclear derivative of the spin-polarized semilocal XC-kernel
! action needed by the differentiated MRSF orbital-adjoint equations.
module mod_dft_gridint_mrsf_xc_kernel_derivative

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t,xc_consumer_t,xc_options_t,run_xc, &
    xc_der1,xc_der2_contr,xc_der3_contr
  use mod_dft_partition_hessian, only: &
    partition_weight_nuclear_first_derivatives
  use mod_dft_gridint_mrsf_xc_fock_deriv_point, only: &
    moving_ao_pair_derivative, &
    spin_fock_point_derivative
  use mod_dft_gga_nuclear_point, only: &
    gga_density_nuclear_point_first_batch
  use mod_dft_gridint_tdgga_response, only: gga_add_owner_motion_first
  use mod_dft_gridint_mrsf_xc_hessian_point, only: &
    mrsf_xc_kernel_fock_coefficients

  implicit none
  private

  type :: kernel_point_workspace_t
    real(fp), allocatable :: dg_rho(:,:),dp_rho(:,:)
    real(fp), allocatable :: dg_grad(:,:,:),dp_grad(:,:,:)
    real(fp), allocatable :: fixed_d(:,:,:,:),fixed_g(:,:,:,:,:)
    real(fp), allocatable :: total_d(:,:,:,:),total_g(:,:,:,:,:)
    real(fp), allocatable :: response_value(:,:,:)
    real(fp), allocatable :: response_gradient(:,:,:,:)
    real(fp), allocatable :: weights(:),dweights(:,:,:)
    real(fp), allocatable :: ground(:),probe(:),dground(:,:),dprobe(:,:)
    real(fp), allocatable :: first(:),second(:,:),third(:,:,:)
    real(fp), allocatable :: v_r(:),dv_r(:,:)
    real(fp), allocatable :: dweight_flat(:),dpair(:),dgrad_pair(:,:)
    real(fp), allocatable :: point_derivative(:,:),coefficient(:,:)
    real(fp), allocatable :: dcoefficient(:,:,:)
  contains
    procedure :: init => kernel_workspace_init
    procedure :: clean => kernel_workspace_clean
  end type kernel_point_workspace_t

  type, extends(xc_consumer_t) :: kernel_derivative_consumer_t
    real(fp), allocatable :: derivative_a(:,:,:,:,:),derivative_b(:,:,:,:,:)
    type(kernel_point_workspace_t), allocatable :: workspace(:)
    real(fp), allocatable :: probe_a(:,:,:),probe_b(:,:,:)
    real(fp), allocatable :: dground_a(:,:,:),dground_b(:,:,:)
    real(fp), allocatable :: dprobe_a(:,:,:,:),dprobe_b(:,:,:,:)
    real(fp), allocatable :: atom_xyz(:,:),surface_shift(:,:)
    integer, allocatable :: ao_atom(:)
    logical, allocatable :: dummy_atom(:)
    real(fp), allocatable :: worker_error(:)
    integer :: part_fun_type=0
    integer :: nprobe=0
    logical :: has_surface_shift=.false.,is_gga=.false.
  contains
    procedure :: parallel_start => kernel_start
    procedure :: parallel_stop => kernel_stop
    procedure :: update => kernel_update
    procedure :: postUpdate => kernel_post
    procedure :: clean => kernel_clean
  end type kernel_derivative_consumer_t

  public :: mrsf_xc_kernel_fock_total_derivative
  public :: mrsf_xc_kernel_fock_total_derivative_batch

contains

!-----------------------------------------------------------------------------

  subroutine kernel_workspace_init(self,nat,nprobe,is_gga)
    class(kernel_point_workspace_t), intent(inout) :: self
    integer, intent(in) :: nat,nprobe
    logical, intent(in) :: is_gga
    integer :: ncart,nvar

    call self%clean()
    ncart=3*nat
    nvar=merge(5,2,is_gga)
    allocate(self%dg_rho(2,ncart),self%dp_rho(2,ncart), &
      self%dg_grad(3,2,ncart),self%dp_grad(3,2,ncart), &
      self%fixed_d(3,nat,2,nprobe+1), &
      self%fixed_g(3,3,nat,2,nprobe+1), &
      self%total_d(3,nat,2,nprobe+1), &
      self%total_g(3,3,nat,2,nprobe+1), &
      self%response_value(2,nprobe+1,ncart), &
      self%response_gradient(3,2,nprobe+1,ncart), &
      self%weights(nat),self%dweights(3,nat,nat),self%ground(nvar), &
      self%probe(nvar),self%dground(nvar,ncart),self%dprobe(nvar,ncart), &
      self%first(nvar),self%second(nvar,nvar),self%third(nvar,nvar,nvar), &
      self%v_r(2),self%dv_r(2,ncart),self%dweight_flat(ncart), &
      self%dpair(ncart),self%dgrad_pair(3,ncart), &
      self%point_derivative(2,ncart),self%coefficient(3,2), &
      self%dcoefficient(3,2,ncart))
  end subroutine kernel_workspace_init

  subroutine kernel_workspace_clean(self)
    class(kernel_point_workspace_t), intent(inout) :: self
    if(allocated(self%dg_rho)) &
      deallocate(self%dg_rho,self%dp_rho,self%dg_grad,self%dp_grad)
    if(allocated(self%fixed_d)) &
      deallocate(self%fixed_d,self%fixed_g,self%total_d,self%total_g, &
        self%response_value,self%response_gradient)
    if(allocated(self%weights)) &
      deallocate(self%weights,self%dweights)
    if(allocated(self%ground)) &
      deallocate(self%ground,self%probe,self%dground,self%dprobe)
    if(allocated(self%first)) &
      deallocate(self%first,self%second,self%third,self%v_r,self%dv_r)
    if(allocated(self%dweight_flat)) &
      deallocate(self%dweight_flat,self%dpair,self%dgrad_pair, &
        self%point_derivative,self%coefficient,self%dcoefficient)
  end subroutine kernel_workspace_clean

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

    real(fp), allocatable :: qa_batch(:,:,:),qb_batch(:,:,:)
    real(fp), allocatable :: dqa_batch(:,:,:,:),dqb_batch(:,:,:,:)
    real(fp), allocatable :: derivative_a_batch(:,:,:,:), &
      derivative_b_batch(:,:,:,:)
    integer :: nbf,ncart

    nbf=basis%nbf
    ncart=3*infos%mol_prop%natom
    allocate(qa_batch(nbf,nbf,1),qb_batch(nbf,nbf,1), &
      dqa_batch(nbf,nbf,ncart,1),dqb_batch(nbf,nbf,ncart,1), &
      derivative_a_batch(nbf,nbf,ncart,1), &
      derivative_b_batch(nbf,nbf,ncart,1),source=0.0_fp)
    qa_batch(:,:,1)=qa; qb_batch(:,:,1)=qb
    dqa_batch(:,:,:,1)=dqa; dqb_batch(:,:,:,1)=dqb
    call mrsf_xc_kernel_fock_total_derivative_batch(basis,mol_grid,da,db, &
      qa_batch,qb_batch,dda,ddb,dqa_batch,dqb_batch,derivative_a_batch, &
      derivative_b_batch,infos,status,threshold)
    derivative_a=derivative_a_batch(:,:,:,1)
    derivative_b=derivative_b_batch(:,:,:,1)
    deallocate(qa_batch,qb_batch,dqa_batch,dqb_batch, &
      derivative_a_batch,derivative_b_batch)
  end subroutine mrsf_xc_kernel_fock_total_derivative

!> Exact multi-probe derivative of K_xc[D](Q_p).  The physical reference D
!> and its nuclear response are common to every probe.  All probes share one
!> molecular-grid traversal and the same AO values; no grid, density, or
!> functional approximation is introduced.
  subroutine mrsf_xc_kernel_fock_total_derivative_batch(basis,mol_grid,da,db, &
      qa,qb,dda,ddb,dqa,dqb,derivative_a,derivative_b,infos,status,threshold)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set), intent(in) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: da(:,:),db(:,:),qa(:,:,:),qb(:,:,:)
    real(fp), intent(in) :: dda(:,:,:),ddb(:,:,:),dqa(:,:,:,:),dqb(:,:,:,:)
    real(fp), intent(out) :: derivative_a(:,:,:,:),derivative_b(:,:,:,:)
    type(information), target, intent(in) :: infos
    integer, intent(out) :: status
    real(fp), intent(in), optional :: threshold

    type(kernel_derivative_consumer_t) :: dat
    type(xc_options_t) :: opts
    real(fp), allocatable, target :: da_normalized(:,:),db_normalized(:,:)
    real(fp) :: grid_threshold
    integer :: first,i,k,last,natom,nbf,ncart,nprobe,probe,shell

    nbf=basis%nbf
    natom=infos%mol_prop%natom
    ncart=3*natom
    nprobe=size(qa,3)
    status=0
    derivative_a=0.0_fp
    derivative_b=0.0_fp
    if(nbf<=0 .or. natom<=0 .or. nprobe<=0 .or. &
       any(shape(da)/=[nbf,nbf]) .or. any(shape(db)/=[nbf,nbf]) .or. &
       any(shape(qa)/=[nbf,nbf,nprobe]) .or. &
       any(shape(qb)/=[nbf,nbf,nprobe]) .or. &
       any(shape(dda)/=[nbf,nbf,ncart]) .or. &
       any(shape(ddb)/=[nbf,nbf,ncart]) .or. &
       any(shape(dqa)/=[nbf,nbf,ncart,nprobe]) .or. &
       any(shape(dqb)/=[nbf,nbf,ncart,nprobe]) .or. &
       any(shape(derivative_a)/=[nbf,nbf,ncart,nprobe]) .or. &
       any(shape(derivative_b)/=[nbf,nbf,ncart,nprobe])) then
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
      dat%probe_a(nbf,nbf,nprobe),dat%probe_b(nbf,nbf,nprobe), &
      dat%dground_a(nbf,nbf,ncart),dat%dground_b(nbf,nbf,ncart), &
      dat%dprobe_a(nbf,nbf,ncart,nprobe), &
      dat%dprobe_b(nbf,nbf,ncart,nprobe), &
      dat%ao_atom(nbf))
    do i=1,nbf
      da_normalized(:,i)=da(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      db_normalized(:,i)=db(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      do k=1,ncart
        dat%dground_a(:,i,k)=dda(:,i,k)*basis%bfnrm(:)*basis%bfnrm(i)
        dat%dground_b(:,i,k)=ddb(:,i,k)*basis%bfnrm(:)*basis%bfnrm(i)
      end do
      do probe=1,nprobe
        dat%probe_a(:,i,probe)=qa(:,i,probe)*basis%bfnrm(:)*basis%bfnrm(i)
        dat%probe_b(:,i,probe)=qb(:,i,probe)*basis%bfnrm(:)*basis%bfnrm(i)
        do k=1,ncart
          dat%dprobe_a(:,i,k,probe)=dqa(:,i,k,probe)* &
            basis%bfnrm(:)*basis%bfnrm(i)
          dat%dprobe_b(:,i,k,probe)=dqb(:,i,k,probe)* &
            basis%bfnrm(:)*basis%bfnrm(i)
        end do
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
    dat%nprobe=nprobe
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
      do probe=1,nprobe
        do k=1,ncart
          do i=1,nbf
            derivative_a(:,i,k,probe)=dat%derivative_a(:,i,k,probe,1)* &
              basis%bfnrm(:)*basis%bfnrm(i)
            derivative_b(:,i,k,probe)=dat%derivative_b(:,i,k,probe,1)* &
              basis%bfnrm(:)*basis%bfnrm(i)
          end do
        end do
      end do
    end if
    call dat%clean()
    deallocate(da_normalized,db_normalized)
  end subroutine mrsf_xc_kernel_fock_total_derivative_batch

!-----------------------------------------------------------------------------

  subroutine kernel_start(self,xce,nthreads)
    class(kernel_derivative_consumer_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: ncart,thread
    ncart=3*xce%numAtoms
    allocate(self%derivative_a(xce%numAOs,xce%numAOs,ncart,self%nprobe,nthreads), &
      self%derivative_b(xce%numAOs,xce%numAOs,ncart,self%nprobe,nthreads), &
      self%worker_error(nthreads),source=0.0_fp)
    allocate(self%workspace(nthreads))
    do thread=1,nthreads
      call self%workspace(thread)%init(xce%numAtoms,self%nprobe,self%is_gga)
    end do
  end subroutine kernel_start

  subroutine kernel_stop(self)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    if(size(self%derivative_a,5)>1) then
      self%derivative_a(:,:,:,:,1)=sum(self%derivative_a,dim=5)
      self%derivative_b(:,:,:,:,1)=sum(self%derivative_b,dim=5)
      self%worker_error(1)=sum(self%worker_error)
    end if
    call self%pe%allreduce(self%derivative_a(:,:,:,:,1), &
      size(self%derivative_a(:,:,:,:,1)))
    call self%pe%allreduce(self%derivative_b(:,:,:,:,1), &
      size(self%derivative_b(:,:,:,:,1)))
    call self%pe%allreduce(self%worker_error(1),1)
  end subroutine kernel_stop

  subroutine kernel_post(self,xce,mythread)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
  end subroutine kernel_post

  subroutine kernel_clean(self)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    integer :: thread
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
    if(allocated(self%workspace)) then
      do thread=1,size(self%workspace)
        call self%workspace(thread)%clean()
      end do
      deallocate(self%workspace)
    end if
  end subroutine kernel_clean

!-----------------------------------------------------------------------------

  subroutine kernel_update(self,xce,mythread)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(fp), allocatable :: ga(:,:),gb(:,:),qa(:,:,:),qb(:,:,:)
    real(fp), allocatable :: base_density(:,:,:,:)
    real(fp), allocatable :: dga(:,:,:),dgb(:,:,:),dqa(:,:,:,:),dqb(:,:,:,:)
    real(fp), allocatable :: v_r_all(:,:),coefficient_all(:,:,:), &
      dv_r_all(:,:,:),dcoefficient_all(:,:,:,:)
    integer, allocatable :: atoms(:)
    integer :: n,ipt,probe,status

    n=xce%numAOs_p
    allocate(base_density(self%nprobe+1,2,n,n), &
      v_r_all(2,self%nprobe),coefficient_all(3,2,self%nprobe), &
      dv_r_all(2,3*xce%numAtoms,self%nprobe), &
      dcoefficient_all(3,2,3*xce%numAtoms,self%nprobe),atoms(n))
    if(xce%skip_p) then
      atoms=self%ao_atom
      base_density(1,1,:,:)=xce%wfAlpha_p
      base_density(1,2,:,:)=xce%wfBeta_p
      do probe=1,self%nprobe
        base_density(probe+1,1,:,:)=self%probe_a(:,:,probe)
        base_density(probe+1,2,:,:)=self%probe_b(:,:,probe)
      end do
      do ipt=1,xce%numPts
        call accumulate_kernel_point_batch(self,xce,mythread,ipt, &
          base_density,self%probe_a,self%probe_b,self%dground_a, &
          self%dground_b,self%dprobe_a,self%dprobe_b,atoms, &
          self%workspace(mythread),v_r_all,coefficient_all,dv_r_all, &
          dcoefficient_all,status)
        if(status/=0) then
          self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
          exit
        end if
      end do
    else
      allocate(ga(n,n),gb(n,n),qa(n,n,self%nprobe),qb(n,n,self%nprobe), &
        dga(n,n,3*xce%numAtoms),dgb(n,n,3*xce%numAtoms), &
        dqa(n,n,3*xce%numAtoms,self%nprobe), &
        dqb(n,n,3*xce%numAtoms,self%nprobe))
      ga=xce%wfAlpha_p
      gb=xce%wfBeta_p
      qa=self%probe_a(xce%indices_p(1:n),xce%indices_p(1:n),:)
      qb=self%probe_b(xce%indices_p(1:n),xce%indices_p(1:n),:)
      dga=self%dground_a(xce%indices_p(1:n),xce%indices_p(1:n),:)
      dgb=self%dground_b(xce%indices_p(1:n),xce%indices_p(1:n),:)
      dqa=self%dprobe_a(xce%indices_p(1:n),xce%indices_p(1:n),:,:)
      dqb=self%dprobe_b(xce%indices_p(1:n),xce%indices_p(1:n),:,:)
      atoms=self%ao_atom(xce%indices_p(1:n))
      base_density(1,1,:,:)=ga
      base_density(1,2,:,:)=gb
      do probe=1,self%nprobe
        base_density(probe+1,1,:,:)=qa(:,:,probe)
        base_density(probe+1,2,:,:)=qb(:,:,probe)
      end do
      do ipt=1,xce%numPts
        call accumulate_kernel_point_batch(self,xce,mythread,ipt, &
          base_density,qa,qb,dga,dgb,dqa,dqb,atoms, &
          self%workspace(mythread),v_r_all,coefficient_all,dv_r_all, &
          dcoefficient_all,status)
        if(status/=0) then
          self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
          exit
        end if
      end do
      deallocate(ga,gb,qa,qb,dga,dgb,dqa,dqb)
    end if
    deallocate(base_density,v_r_all,coefficient_all,dv_r_all, &
      dcoefficient_all,atoms)
  end subroutine kernel_update

!-----------------------------------------------------------------------------

  subroutine accumulate_kernel_point_batch(self,xce,mythread,ipt,base_density, &
      qa,qb,dga,dgb,dqa,dqb,atoms,workspace,v_r_all,coefficient_all, &
      dv_r_all,dcoefficient_all,status)
    class(kernel_derivative_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: mythread,ipt,atoms(:)
    real(fp), intent(inout) :: base_density(:,:,:,:)
    real(fp), intent(in) :: qa(:,:,:),qb(:,:,:),dga(:,:,:),dgb(:,:,:), &
      dqa(:,:,:,:),dqb(:,:,:,:)
    type(kernel_point_workspace_t), intent(inout) :: workspace
    real(fp), intent(out) :: v_r_all(:,:),coefficient_all(:,:,:), &
      dv_r_all(:,:,:),dcoefficient_all(:,:,:,:)
    integer, intent(out) :: status

    real(fp) :: rho(2),prho(2),grad_rho(3,2),grad_probe(3,2)
    real(fp) :: value(2,size(qa,3)+1),gradient(3,2,size(qa,3)+1)
    real(fp) :: pair,grad_pair(3),finite_weight,quadrature_scale
    integer :: atom,cart,coordinate,field,global_mu,global_nu,i,j,k,local_status
    integer :: mu,n,nat,ncart,nprobe,nvar,nu,owner,probe_index,spin

    n=size(base_density,3)
    nat=xce%numAtoms
    ncart=3*nat
    nprobe=size(qa,3)
    nvar=merge(5,2,self%is_gga)
    owner=xce%currAtom
    status=0
    finite_weight=xce%xyzw(ipt,4)
    if(abs(finite_weight)<=tiny(1.0_fp)) return
    associate(dg_rho=>workspace%dg_rho,dp_rho=>workspace%dp_rho, &
      dg_grad=>workspace%dg_grad,dp_grad=>workspace%dp_grad, &
      fixed_d=>workspace%fixed_d,fixed_g=>workspace%fixed_g, &
      total_d=>workspace%total_d,total_g=>workspace%total_g, &
      response_value=>workspace%response_value, &
      response_gradient=>workspace%response_gradient, &
      weights=>workspace%weights,dweights=>workspace%dweights, &
      ground=>workspace%ground, &
      probe=>workspace%probe,dground=>workspace%dground, &
      dprobe=>workspace%dprobe,first=>workspace%first, &
      second=>workspace%second,third=>workspace%third,v_r=>workspace%v_r, &
      dv_r=>workspace%dv_r,dweight_flat=>workspace%dweight_flat, &
      dpair=>workspace%dpair,dgrad_pair=>workspace%dgrad_pair, &
      point_derivative=>workspace%point_derivative, &
      coefficient=>workspace%coefficient,dcoefficient=>workspace%dcoefficient)
    call partition_weight_nuclear_first_derivatives(self%atom_xyz, &
      xce%xyzw(ipt,1:3),owner,self%dummy_atom,self%part_fun_type, &
      self%has_surface_shift,self%surface_shift,weights,dweights, &
      local_status)
    if(local_status/=0 .or. weights(owner)<=sqrt(tiny(1.0_fp))) then
      status=-2
      return
    end if
    quadrature_scale=finite_weight/weights(owner)
    do atom=1,nat
      do cart=1,3
        coordinate=3*(atom-1)+cart
        dweight_flat(coordinate)=dweights(cart,atom,owner)
      end do
    end do

    call build_unweighted_kernels(xce,ipt,nvar,finite_weight,first,second, &
      third)
    call field_value_gradient_batch(base_density,xce%aoV(:,ipt), &
      xce%aoG1(:,ipt,:),value,gradient)
    call gga_density_nuclear_point_first_batch(base_density,atoms, &
      xce%aoV(:,ipt),xce%aoG1(:,ipt,:),xce%aoG2(:,ipt,:),fixed_d,fixed_g)
    do field=1,nprobe+1
      do spin=1,2
        call gga_add_owner_motion_first(owner,fixed_d(:,:,spin,field), &
          fixed_g(:,:,:,spin,field),total_d(:,:,spin,field), &
          total_g(:,:,:,spin,field))
      end do
    end do
    call response_field_value_gradient_batch(dga,dgb,dqa,dqb, &
      xce%aoV(:,ipt),xce%aoG1(:,ipt,:),response_value,response_gradient)
    do probe_index=1,nprobe
      rho=value(:,1)
      prho=value(:,probe_index+1)
      grad_rho=gradient(:,:,1)
      grad_probe=gradient(:,:,probe_index+1)
      do field=1,2
        do spin=1,2
          if(field==1) then
            dg_rho(spin,:)=reshape(total_d(:,:,spin,1),[ncart])+ &
              response_value(spin,1,:)
            do coordinate=1,ncart
              atom=(coordinate-1)/3+1
              cart=coordinate-3*(atom-1)
              dg_grad(:,spin,coordinate)=total_g(:,cart,atom,spin,1)+ &
                response_gradient(:,spin,1,coordinate)
            end do
          else
            dp_rho(spin,:)=reshape(total_d(:,:,spin,probe_index+1),[ncart])+ &
              response_value(spin,probe_index+1,:)
            do coordinate=1,ncart
              atom=(coordinate-1)/3+1
              cart=coordinate-3*(atom-1)
              dp_grad(:,spin,coordinate)= &
                total_g(:,cart,atom,spin,probe_index+1)+ &
                response_gradient(:,spin,probe_index+1,coordinate)
            end do
          end if
        end do
      end do
      call build_density_variables(self%is_gga,rho,grad_rho,prho,grad_probe, &
        dg_rho,dg_grad,dp_rho,dp_grad,ground,probe,dground,dprobe)
      call mrsf_xc_kernel_fock_coefficients(first,second,third,probe,dground, &
        dprobe,self%is_gga,grad_rho,grad_probe,dg_grad,dp_grad,v_r, &
        coefficient,dv_r,dcoefficient,local_status)
      if(local_status/=0) then
        status=-5
        return
      end if
      v_r_all(:,probe_index)=v_r
      coefficient_all(:,:,probe_index)=coefficient
      dv_r_all(:,:,probe_index)=dv_r
      dcoefficient_all(:,:,:,probe_index)=dcoefficient
    end do

    do nu=1,n
      global_nu=nu
      if(.not.xce%skip_p) global_nu=xce%indices_p(nu)
      do mu=1,nu
        global_mu=mu
        if(.not.xce%skip_p) global_mu=xce%indices_p(mu)
        call moving_ao_pair_derivative(owner,atoms(mu),atoms(nu), &
          xce%aoV(mu,ipt),xce%aoV(nu,ipt),xce%aoG1(mu,ipt,:), &
          xce%aoG1(nu,ipt,:),xce%aoG2(mu,ipt,:),xce%aoG2(nu,ipt,:), &
          pair,grad_pair,dpair,dgrad_pair,local_status)
        if(local_status/=0) then
          status=-3
          return
        end if
        do probe_index=1,nprobe
          call spin_fock_point_derivative(quadrature_scale,weights(owner), &
            dweight_flat,v_r_all(:,probe_index), &
            coefficient_all(:,:,probe_index),dv_r_all(:,:,probe_index), &
            dcoefficient_all(:,:,:,probe_index),pair,grad_pair,dpair, &
            dgrad_pair,point_derivative,local_status)
          if(local_status/=0) then
            status=-4
            return
          end if
          do coordinate=1,ncart
            self%derivative_a(global_mu,global_nu,coordinate,probe_index,mythread)= &
              self%derivative_a(global_mu,global_nu,coordinate,probe_index,mythread)+ &
              point_derivative(1,coordinate)
            self%derivative_b(global_mu,global_nu,coordinate,probe_index,mythread)= &
              self%derivative_b(global_mu,global_nu,coordinate,probe_index,mythread)+ &
              point_derivative(2,coordinate)
            if(global_mu/=global_nu) then
              self%derivative_a(global_nu,global_mu,coordinate,probe_index,mythread)= &
                self%derivative_a(global_nu,global_mu,coordinate,probe_index,mythread)+ &
                point_derivative(1,coordinate)
              self%derivative_b(global_nu,global_mu,coordinate,probe_index,mythread)= &
                self%derivative_b(global_nu,global_mu,coordinate,probe_index,mythread)+ &
                point_derivative(2,coordinate)
            end if
          end do
        end do
      end do
    end do
    end associate

  end subroutine accumulate_kernel_point_batch

!-----------------------------------------------------------------------------

  pure subroutine field_value_gradient_batch(density,aov,aog1,value,gradient)
    real(fp), intent(in) :: density(:,:,:,:),aov(:),aog1(:,:)
    real(fp), intent(out) :: value(:,:),gradient(:,:,:)
    real(fp) :: pair,grad_pair(3),density_pair(size(density,1),2)
    integer :: field,mu,nu,spin,nao,nfield

    nao=size(aov)
    nfield=size(density,1)
    if(any(shape(density)/=[nfield,2,nao,nao]) .or. &
       any(shape(aog1)/=[nao,3]) .or. any(shape(value)/=[2,nfield]) .or. &
       any(shape(gradient)/=[3,2,nfield])) &
      error stop 'field_value_gradient_batch: shape mismatch'
    value=0.0_fp
    gradient=0.0_fp
    do nu=1,nao
      do mu=1,nu
        pair=aov(mu)*aov(nu)
        grad_pair=aog1(mu,:)*aov(nu)+aov(mu)*aog1(nu,:)
        density_pair=density(:,:,mu,nu)
        if(mu/=nu) density_pair=density_pair+density(:,:,nu,mu)
        do field=1,nfield
          do spin=1,2
            value(spin,field)=value(spin,field)+ &
              density_pair(field,spin)*pair
            gradient(:,spin,field)=gradient(:,spin,field)+ &
              density_pair(field,spin)*grad_pair
          end do
        end do
      end do
    end do
  end subroutine field_value_gradient_batch

  pure subroutine response_field_value_gradient_batch(dga,dgb,dqa,dqb, &
      aov,aog1,value,gradient)
    real(fp), intent(in) :: dga(:,:,:),dgb(:,:,:),dqa(:,:,:,:),dqb(:,:,:,:)
    real(fp), intent(in) :: aov(:),aog1(:,:)
    real(fp), intent(out) :: value(:,:,:),gradient(:,:,:,:)
    real(fp) :: pair,grad_pair(3)
    real(fp) :: pair_density(2,size(dqa,4)+1,size(dga,3))
    integer :: coordinate,field,mu,nu,nao,ncart,nfield

    nao=size(aov)
    ncart=size(dga,3)
    nfield=size(dqa,4)+1
    if(any(shape(dga)/=[nao,nao,ncart]) .or. &
       any(shape(dgb)/=[nao,nao,ncart]) .or. &
       any(shape(dqa)/=[nao,nao,ncart,nfield-1]) .or. &
       any(shape(dqb)/=[nao,nao,ncart,nfield-1]) .or. &
       any(shape(aog1)/=[nao,3]) .or. &
       any(shape(value)/=[2,nfield,ncart]) .or. &
       any(shape(gradient)/=[3,2,nfield,ncart])) &
      error stop 'response_field_value_gradient_batch: shape mismatch'
    value=0.0_fp
    gradient=0.0_fp
    do nu=1,nao
      do mu=1,nu
        pair=aov(mu)*aov(nu)
        grad_pair=aog1(mu,:)*aov(nu)+aov(mu)*aog1(nu,:)
        pair_density(1,1,:)=dga(mu,nu,:)
        pair_density(2,1,:)=dgb(mu,nu,:)
        do field=2,nfield
          pair_density(1,field,:)=dqa(mu,nu,:,field-1)
          pair_density(2,field,:)=dqb(mu,nu,:,field-1)
        end do
        if(mu/=nu) then
          pair_density(1,1,:)=pair_density(1,1,:)+dga(nu,mu,:)
          pair_density(2,1,:)=pair_density(2,1,:)+dgb(nu,mu,:)
          do field=2,nfield
            pair_density(1,field,:)=pair_density(1,field,:)+ &
              dqa(nu,mu,:,field-1)
            pair_density(2,field,:)=pair_density(2,field,:)+ &
              dqb(nu,mu,:,field-1)
          end do
        end if
        do coordinate=1,ncart
          do field=1,nfield
            value(:,field,coordinate)=value(:,field,coordinate)+ &
              pair_density(:,field,coordinate)*pair
            gradient(:,1,field,coordinate)= &
              gradient(:,1,field,coordinate)+ &
              pair_density(1,field,coordinate)*grad_pair
            gradient(:,2,field,coordinate)= &
              gradient(:,2,field,coordinate)+ &
              pair_density(2,field,coordinate)*grad_pair
          end do
        end do
      end do
    end do
  end subroutine response_field_value_gradient_batch

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
    real(fp) :: square(5,5),mixed(5)
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
    square=0.0_fp
    mixed=0.0_fp
    do j=1,nvar
      call split_basis(nvar,j,dr,ds)
      dt=0.0_fp
      call xc_der3_contr(xce,ipt,dr,ds,dt,[0.0_fp,0.0_fp,0.0_fp], &
        ffs,gr,gs,gt)
      call join_direction(nvar,gr,gs,square(1:nvar,j))
    end do
    third=0.0_fp
    do j=1,nvar
      third(:,j,j)=square(1:nvar,j)/finite_weight
      do k=j+1,nvar
        call split_pair(nvar,j,k,dr,ds)
        dt=0.0_fp
        call xc_der3_contr(xce,ipt,dr,ds,dt,[0.0_fp,0.0_fp,0.0_fp], &
          ffs,gr,gs,gt)
        call join_direction(nvar,gr,gs,mixed(1:nvar))
        mixed(1:nvar)=0.5_fp*(mixed(1:nvar)-square(1:nvar,j)- &
          square(1:nvar,k))/finite_weight
        third(:,j,k)=mixed(1:nvar)
        third(:,k,j)=mixed(1:nvar)
      end do
    end do
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
