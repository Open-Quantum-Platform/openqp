! Analytic spin-polarized LDA/GGA quadrature contribution to the founding
! two-SOMO MRSF-TDDFT Hessian.
module mod_dft_gridint_mrsf_xc_hessian

  use precision, only: fp
  use mod_dft_gridint, only: xc_engine_t,xc_consumer_t,xc_options_t,run_xc, &
    xc_der1,xc_der2_contr,xc_der3_contr
  use mod_dft_gga_nuclear_point, only: gga_density_nuclear_point_batch, &
    gga_density_nuclear_point_first_batch
  use mod_dft_gridint_tdgga_hessian, only: gga_add_owner_motion
  use mod_dft_gridint_tdgga_response, only: &
    gga_add_owner_motion_first
  use mod_dft_partition_hessian, only: partition_weight_nuclear_derivatives, &
    partition_derivative_workspace_t
  use mod_dft_gridint_mrsf_xc_hessian_point, only: &
    mrsf_xc_weighted_fixed_hessian,mrsf_xc_weighted_response_rows

  implicit none
  private

  type :: mrsf_xc_point_workspace_t
    type(partition_derivative_workspace_t) :: partition
    real(fp), allocatable :: fixed_d1(:,:,:),fixed_g1(:,:,:,:)
    real(fp), allocatable :: fixed_d2(:,:,:,:,:),fixed_g2(:,:,:,:,:,:)
    real(fp), allocatable :: d1(:,:,:),g1(:,:,:,:)
    real(fp), allocatable :: d2(:,:,:,:,:),g2(:,:,:,:,:,:)
    real(fp), allocatable :: weights(:),dweights(:,:,:),d2weights(:,:,:,:,:)
    real(fp), allocatable :: ug(:),up(:),dug(:,:,:),dup(:,:,:)
    real(fp), allocatable :: d2ug(:,:,:,:,:),d2up(:,:,:,:,:)
    real(fp), allocatable :: first(:),second(:,:),third(:,:,:)
    real(fp), allocatable :: response_value(:,:,:)
    real(fp), allocatable :: response_gradient(:,:,:,:)
    real(fp), allocatable :: response_fixed_d(:,:,:,:,:)
    real(fp), allocatable :: response_fixed_g(:,:,:,:,:,:)
    real(fp), allocatable :: response_d1(:,:,:,:,:)
    real(fp), allocatable :: response_g1(:,:,:,:,:,:)
    real(fp), allocatable :: rg(:,:),rp(:,:),drg(:,:,:,:),drp(:,:,:,:), &
      row(:,:,:)
  contains
    procedure :: init => mrsf_xc_workspace_init
    procedure :: clean => mrsf_xc_workspace_clean
  end type mrsf_xc_point_workspace_t

  type, extends(xc_consumer_t) :: mrsf_xc_hessian_consumer_t
    real(fp), allocatable :: hessian(:,:,:,:,:),rows(:,:,:,:)
    type(mrsf_xc_point_workspace_t), allocatable :: workspace(:)
    real(fp), allocatable :: base_density(:,:,:)
    real(fp), allocatable :: response_density(:,:,:,:)
    real(fp), allocatable :: atom_xyz(:,:),surface_shift(:,:)
    integer, allocatable :: ao_atom(:)
    logical, allocatable :: dummy_atom(:)
    real(fp), allocatable :: worker_error(:)
    integer :: part_fun_type=0
    logical :: has_surface_shift=.false.,is_gga=.false.
  contains
    procedure :: parallel_start => mrsf_xc_start
    procedure :: parallel_stop => mrsf_xc_stop
    procedure :: update => mrsf_xc_update
    procedure :: postUpdate => mrsf_xc_post
    procedure :: clean => mrsf_xc_clean
  end type mrsf_xc_hessian_consumer_t

  public :: mrsf_tddft_xc_hessian_rows

contains

!> Evaluate the fixed-density XC Hessian and all density-response rows.
!>
!> density_d and density_p are the alpha/beta reference and relaxed MRSF
!> densities.  response_d and response_p are their total first nuclear
!> derivatives in the AO coefficient representation.  Only these physical
!> spin densities enter; no state expansion is formed.
  subroutine mrsf_tddft_xc_hessian_rows(basis,mol_grid,density_d,density_p, &
      response_d,response_p,hessian,rows,infos,status)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information

    type(basis_set), intent(in) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: density_d(:,:,:),density_p(:,:,:)
    real(fp), intent(in) :: response_d(:,:,:,:),response_p(:,:,:,:)
    real(fp), intent(out) :: hessian(:,:),rows(:,:)
    type(information), target, intent(in) :: infos
    integer, intent(out) :: status

    type(mrsf_xc_hessian_consumer_t) :: dat
    type(xc_options_t) :: opts
    real(fp), allocatable, target :: da(:,:),db(:,:)
    integer :: nbf,natom,ncart,i,k,shell,first,last,atom_a,atom_b,cart_a,cart_b

    nbf=basis%nbf
    natom=infos%mol_prop%natom
    ncart=3*natom
    status=0
    hessian=0.0_fp
    rows=0.0_fp
    if(nbf<=0 .or. natom<=0 .or. &
       any(shape(density_d)/=[nbf,nbf,2]) .or. &
       any(shape(density_p)/=[nbf,nbf,2]) .or. &
       any(shape(response_d)/=[nbf,nbf,2,ncart]) .or. &
       any(shape(response_p)/=[nbf,nbf,2,ncart]) .or. &
       any(shape(hessian)/=[ncart,ncart]) .or. &
       any(shape(rows)/=[ncart,ncart])) then
      status=-1
      return
    end if
    if(infos%functional%needTau .or. infos%functional%needLapl .or. &
       infos%dft%cam_flag) then
      status=-2
      return
    end if
    allocate(da(nbf,nbf),db(nbf,nbf),dat%base_density(4,nbf,nbf), &
      dat%response_density(2*ncart,2,nbf,nbf),dat%ao_atom(nbf))
    do i=1,nbf
      da(:,i)=density_d(:,i,1)*basis%bfnrm(:)*basis%bfnrm(i)
      db(:,i)=density_d(:,i,2)*basis%bfnrm(:)*basis%bfnrm(i)
      dat%base_density(1,:,i)=da(:,i)
      dat%base_density(2,:,i)=db(:,i)
      dat%base_density(3,:,i)=density_p(:,i,1)* &
        basis%bfnrm(:)*basis%bfnrm(i)
      dat%base_density(4,:,i)=density_p(:,i,2)* &
        basis%bfnrm(:)*basis%bfnrm(i)
      do k=1,ncart
        dat%response_density(k,1,:,i)=response_d(:,i,1,k)* &
          basis%bfnrm(:)*basis%bfnrm(i)
        dat%response_density(k,2,:,i)=response_d(:,i,2,k)* &
          basis%bfnrm(:)*basis%bfnrm(i)
        dat%response_density(ncart+k,1,:,i)=response_p(:,i,1,k)* &
          basis%bfnrm(:)*basis%bfnrm(i)
        dat%response_density(ncart+k,2,:,i)=response_p(:,i,2,k)* &
          basis%bfnrm(:)*basis%bfnrm(i)
      end do
    end do
    dat%ao_atom=0
    do shell=1,basis%nshell
      first=basis%ao_offset(shell)
      last=first+basis%naos(shell)-1
      dat%ao_atom(first:last)=basis%origin(shell)
    end do
    if(any(dat%ao_atom<1) .or. any(dat%ao_atom>natom)) then
      status=-3
      call dat%clean()
      deallocate(da,db)
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
    ! GGA adds one AO derivative internally.  Both paths therefore collocate
    ! through G3, which is required by the fixed-density second derivative.
    opts%nDer=merge(2,3,dat%is_gga)
    opts%nXCDer=3
    opts%numOccAlpha=infos%mol_prop%nelec_A
    opts%numOccBeta=infos%mol_prop%nelec_B
    opts%wfAlpha=>da
    opts%wfBeta=>db
    opts%dft_threshold=0.0_fp
    opts%ao_threshold=0.0_fp
    opts%ao_sparsity_ratio=0.0_fp
    opts%molGrid=>mol_grid

    call dat%pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    call run_xc(opts,dat,basis)
    if(allocated(dat%worker_error)) then
      if(dat%worker_error(1)>0.0_fp) status=-4
    end if
    if(status==0) then
      do atom_b=1,natom
        do cart_b=1,3
          do atom_a=1,natom
            do cart_a=1,3
              hessian(3*(atom_a-1)+cart_a,3*(atom_b-1)+cart_b)= &
                dat%hessian(cart_a,cart_b,atom_a,atom_b,1)
            end do
          end do
        end do
      end do
      do k=1,ncart
        do atom_a=1,natom
          do cart_a=1,3
            rows(3*(atom_a-1)+cart_a,k)=dat%rows(cart_a,atom_a,k,1)
          end do
        end do
      end do
    end if
    call dat%clean()
    deallocate(da,db)
  end subroutine mrsf_tddft_xc_hessian_rows

!-----------------------------------------------------------------------------

  subroutine mrsf_xc_workspace_init(self,nat,is_gga)
    class(mrsf_xc_point_workspace_t), intent(inout) :: self
    integer, intent(in) :: nat
    logical, intent(in) :: is_gga
    integer :: ncart,nvar

    call self%clean()
    ncart=3*nat
    nvar=merge(5,2,is_gga)
    allocate(self%fixed_d1(3,nat,4),self%fixed_g1(3,3,nat,4), &
      self%fixed_d2(3,3,nat,nat,4), &
      self%fixed_g2(3,3,3,nat,nat,4), &
      self%d1(3,nat,4),self%g1(3,3,nat,4), &
      self%d2(3,3,nat,nat,4),self%g2(3,3,3,nat,nat,4), &
      self%weights(nat),self%dweights(3,nat,nat), &
      self%d2weights(3,nat,3,nat,nat),self%ug(nvar),self%up(nvar), &
      self%dug(nvar,3,nat),self%dup(nvar,3,nat), &
      self%d2ug(nvar,3,3,nat,nat),self%d2up(nvar,3,3,nat,nat), &
      self%first(nvar),self%second(nvar,nvar),self%third(nvar,nvar,nvar), &
      self%response_value(2,2,ncart), &
      self%response_gradient(3,2,2,ncart), &
      self%response_fixed_d(3,nat,2,2,ncart), &
      self%response_fixed_g(3,3,nat,2,2,ncart), &
      self%response_d1(3,nat,2,2,ncart), &
      self%response_g1(3,3,nat,2,2,ncart), &
      self%rg(nvar,ncart),self%rp(nvar,ncart), &
      self%drg(nvar,3,nat,ncart),self%drp(nvar,3,nat,ncart), &
      self%row(3,nat,ncart))
  end subroutine mrsf_xc_workspace_init

  subroutine mrsf_xc_workspace_clean(self)
    class(mrsf_xc_point_workspace_t), intent(inout) :: self
    call self%partition%clean()
    if(allocated(self%fixed_d1)) &
      deallocate(self%fixed_d1,self%fixed_g1,self%fixed_d2,self%fixed_g2)
    if(allocated(self%d1)) deallocate(self%d1,self%g1,self%d2,self%g2)
    if(allocated(self%weights)) &
      deallocate(self%weights,self%dweights,self%d2weights)
    if(allocated(self%ug)) &
      deallocate(self%ug,self%up,self%dug,self%dup,self%d2ug,self%d2up)
    if(allocated(self%first)) &
      deallocate(self%first,self%second,self%third)
    if(allocated(self%response_value)) &
      deallocate(self%response_value,self%response_gradient, &
        self%response_fixed_d,self%response_fixed_g,self%response_d1, &
        self%response_g1)
    if(allocated(self%rg)) &
      deallocate(self%rg,self%rp,self%drg,self%drp,self%row)
  end subroutine mrsf_xc_workspace_clean

!-----------------------------------------------------------------------------

  subroutine mrsf_xc_start(self,xce,nthreads)
    class(mrsf_xc_hessian_consumer_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    integer :: thread
    allocate(self%hessian(3,3,xce%numAtoms,xce%numAtoms,nthreads), &
      self%rows(3,xce%numAtoms,3*xce%numAtoms,nthreads), &
      self%worker_error(nthreads),source=0.0_fp)
    allocate(self%workspace(nthreads))
    do thread=1,nthreads
      call self%workspace(thread)%init(xce%numAtoms,self%is_gga)
    end do
  end subroutine mrsf_xc_start

  subroutine mrsf_xc_stop(self)
    class(mrsf_xc_hessian_consumer_t), intent(inout) :: self
    if(size(self%hessian,5)>1) then
      self%hessian(:,:,:,:,1)=sum(self%hessian,dim=5)
      self%rows(:,:,:,1)=sum(self%rows,dim=4)
      self%worker_error(1)=sum(self%worker_error)
    end if
    call self%pe%allreduce(self%hessian(:,:,:,:,1), &
      size(self%hessian(:,:,:,:,1)))
    call self%pe%allreduce(self%rows(:,:,:,1),size(self%rows(:,:,:,1)))
    call self%pe%allreduce(self%worker_error(1),1)
  end subroutine mrsf_xc_stop

  subroutine mrsf_xc_post(self,xce,mythread)
    class(mrsf_xc_hessian_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
    ! Contributions are accumulated directly in update.
  end subroutine mrsf_xc_post

  subroutine mrsf_xc_clean(self)
    class(mrsf_xc_hessian_consumer_t), intent(inout) :: self
    integer :: thread
    if(allocated(self%hessian)) deallocate(self%hessian)
    if(allocated(self%rows)) deallocate(self%rows)
    if(allocated(self%base_density)) deallocate(self%base_density)
    if(allocated(self%response_density)) deallocate(self%response_density)
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
  end subroutine mrsf_xc_clean

!-----------------------------------------------------------------------------

  subroutine mrsf_xc_update(self,xce,mythread)
    class(mrsf_xc_hessian_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    real(fp), allocatable :: base_density(:,:,:),response_combined(:,:,:,:), &
      response_fixed_combined(:,:,:,:),response_fixed_gradient_combined(:,:,:,:,:)
    integer, allocatable :: atoms(:)
    integer :: n,ncart,ipt,status

    n=xce%numAOs_p
    ncart=3*xce%numAtoms
    allocate(response_fixed_combined(3,xce%numAtoms,2,2*ncart), &
      response_fixed_gradient_combined(3,3,xce%numAtoms,2,2*ncart),atoms(n))
    if(xce%skip_p) then
      atoms=self%ao_atom
      do ipt=1,xce%numPts
        call accumulate_point(self,xce,mythread,ipt,self%base_density, &
          self%response_density,response_fixed_combined, &
          response_fixed_gradient_combined,atoms, &
          self%workspace(mythread),status)
        if(status/=0) then
          self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
          exit
        end if
      end do
    else
      allocate(base_density(4,n,n),response_combined(2*ncart,2,n,n))
      base_density=self%base_density(:,xce%indices_p(1:n), &
        xce%indices_p(1:n))
      response_combined=self%response_density(:,:,xce%indices_p(1:n), &
        xce%indices_p(1:n))
      atoms=self%ao_atom(xce%indices_p(1:n))
      do ipt=1,xce%numPts
        call accumulate_point(self,xce,mythread,ipt,base_density, &
          response_combined,response_fixed_combined, &
          response_fixed_gradient_combined,atoms, &
          self%workspace(mythread),status)
        if(status/=0) then
          self%worker_error(mythread)=self%worker_error(mythread)+1.0_fp
          exit
        end if
      end do
      deallocate(base_density,response_combined)
    end if
    deallocate(response_fixed_combined,response_fixed_gradient_combined,atoms)
  end subroutine mrsf_xc_update

!-----------------------------------------------------------------------------

  subroutine accumulate_point(self,xce,mythread,ipt,base_density, &
      response_combined, &
      response_fixed_combined,response_fixed_gradient_combined,atoms, &
      workspace,status)
    class(mrsf_xc_hessian_consumer_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: mythread,ipt,atoms(:)
    real(fp), intent(in) :: base_density(:,:,:),response_combined(:,:,:,:)
    real(fp), intent(inout) :: response_fixed_combined(:,:,:,:), &
      response_fixed_gradient_combined(:,:,:,:,:)
    type(mrsf_xc_point_workspace_t), intent(inout) :: workspace
    integer, intent(out) :: status

    integer :: nat,ncart,nvar,field,spin,k,matrix,owner,local_status
    real(fp) :: finite_weight,quadrature_scale,exc
    real(fp) :: value(2,2),gradient(3,2,2)

    nat=xce%numAtoms
    ncart=3*nat
    nvar=merge(5,2,self%is_gga)
    owner=xce%currAtom
    status=0
    finite_weight=xce%xyzw(ipt,4)
    if(abs(finite_weight)<=tiny(1.0_fp)) return
    associate(fixed_d1=>workspace%fixed_d1,fixed_g1=>workspace%fixed_g1, &
      fixed_d2=>workspace%fixed_d2,fixed_g2=>workspace%fixed_g2, &
      d1=>workspace%d1,g1=>workspace%g1,d2=>workspace%d2,g2=>workspace%g2, &
      weights=>workspace%weights,dweights=>workspace%dweights, &
      d2weights=>workspace%d2weights,ug=>workspace%ug,up=>workspace%up, &
      dug=>workspace%dug,dup=>workspace%dup,d2ug=>workspace%d2ug, &
      d2up=>workspace%d2up,first=>workspace%first,second=>workspace%second, &
      third=>workspace%third,response_value=>workspace%response_value, &
      response_gradient=>workspace%response_gradient, &
      response_fixed_d=>workspace%response_fixed_d, &
      response_fixed_g=>workspace%response_fixed_g, &
      response_d1=>workspace%response_d1,response_g1=>workspace%response_g1, &
      rg=>workspace%rg,rp=>workspace%rp,drg=>workspace%drg, &
      drp=>workspace%drp,row=>workspace%row)
    fixed_d1=0.0_fp;fixed_g1=0.0_fp;fixed_d2=0.0_fp;fixed_g2=0.0_fp
    d1=0.0_fp;g1=0.0_fp;d2=0.0_fp;g2=0.0_fp
    call gga_density_nuclear_point_batch(base_density,atoms,xce%aoV(:,ipt), &
      xce%aoG1(:,ipt,:),xce%aoG2(:,ipt,:),xce%aoG3(:,ipt,:), &
      fixed_d1,fixed_g1,fixed_d2,fixed_g2)
    call field_value_gradient_batch(base_density,xce%aoV(:,ipt), &
      xce%aoG1(:,ipt,:),value,gradient)
    do field=1,2
      do spin=1,2
        matrix=spin+2*(field-1)
        call gga_add_owner_motion(owner,fixed_d1(:,:,matrix), &
          fixed_g1(:,:,:,matrix),fixed_d2(:,:,:,:,matrix), &
          fixed_g2(:,:,:,:,:,matrix),d1(:,:,matrix), &
          g1(:,:,:,matrix),d2(:,:,:,:,matrix),g2(:,:,:,:,:,matrix))
      end do
    end do
    call build_variables(self%is_gga,value(:,1),gradient(:,:,1), &
      value(:,2),gradient(:,:,2),d1(:,:,1:2),g1(:,:,:,1:2), &
      d1(:,:,3:4),g1(:,:,:,3:4),d2(:,:,:,:,1:2), &
      g2(:,:,:,:,:,1:2),d2(:,:,:,:,3:4),g2(:,:,:,:,:,3:4), &
      ug,up,dug,dup,d2ug,d2up)
    call partition_weight_nuclear_derivatives(self%atom_xyz, &
      xce%xyzw(ipt,1:3),owner,self%dummy_atom,self%part_fun_type, &
      self%has_surface_shift,self%surface_shift,weights,dweights,d2weights, &
      local_status,workspace%partition)
    if(local_status/=0 .or. weights(owner)<=sqrt(tiny(1.0_fp))) then
      status=-1
      return
    end if
    quadrature_scale=finite_weight/weights(owner)
    call build_unweighted_kernels(xce,ipt,nvar,finite_weight,first,second, &
      third)
    exc=xce%xclib%exc(ipt)
    call mrsf_xc_weighted_fixed_hessian(exc,first,second,third,up,dug,dup, &
      d2ug,d2up,quadrature_scale,weights(owner),dweights(:,:,owner), &
      d2weights(:,:,:,:,owner),self%hessian(:,:,:,:,mythread),local_status)
    if(local_status/=0) then
      status=-2
      return
    end if
    call response_field_value_gradient_batch(response_combined, &
      xce%aoV(:,ipt),xce%aoG1(:,ipt,:),response_value,response_gradient)
    call gga_density_nuclear_point_first_batch(response_combined,atoms, &
      xce%aoV(:,ipt),xce%aoG1(:,ipt,:),xce%aoG2(:,ipt,:), &
      response_fixed_combined,response_fixed_gradient_combined)
    do k=1,ncart
      do field=1,2
        do spin=1,2
          matrix=k+(field-1)*ncart
          call gga_add_owner_motion_first(owner, &
            response_fixed_combined(:,:,spin,matrix), &
            response_fixed_gradient_combined(:,:,:,spin,matrix), &
            response_d1(:,:,spin,field,k), &
            response_g1(:,:,:,spin,field,k))
        end do
      end do
      call build_response_variables(self%is_gga,value(:,1),gradient(:,:,1), &
        value(:,2),gradient(:,:,2),d1(:,:,1:2),g1(:,:,:,1:2), &
        d1(:,:,3:4),g1(:,:,:,3:4),response_value(:,1,k), &
        response_gradient(:,:,1,k),response_value(:,2,k), &
        response_gradient(:,:,2,k),response_d1(:,:,:,1,k), &
        response_g1(:,:,:,:,1,k),response_d1(:,:,:,2,k), &
        response_g1(:,:,:,:,2,k),rg(:,k),rp(:,k),drg(:,:,:,k),drp(:,:,:,k))
    end do
    row=0.0_fp
    call mrsf_xc_weighted_response_rows(first,second,third,up,dug,dup,rg, &
      rp,drg,drp,quadrature_scale,weights(owner),dweights(:,:,owner),row, &
      local_status)
    if(local_status/=0) then
      status=-3
      return
    end if
    self%rows(:,:,:,mythread)=self%rows(:,:,:,mythread)+row
    end associate

  end subroutine accumulate_point

!-----------------------------------------------------------------------------

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

!-----------------------------------------------------------------------------

  pure subroutine field_value_gradient_batch(density,aov,aog1,value,gradient)
    real(fp), intent(in) :: density(:,:,:),aov(:),aog1(:,:)
    real(fp), intent(out) :: value(:,:),gradient(:,:,:)
    real(fp) :: pair,grad_pair(3),pair_density(4)
    integer :: field,matrix,mu,nu,spin,nao

    nao=size(aov)
    if(any(shape(density)/=[4,nao,nao]) .or. &
       any(shape(aog1)/=[nao,3]) .or. any(shape(value)/=[2,2]) .or. &
       any(shape(gradient)/=[3,2,2])) &
      error stop 'field_value_gradient_batch: shape mismatch'
    value=0.0_fp
    gradient=0.0_fp
    do nu=1,nao
      do mu=1,nu
        pair=aov(mu)*aov(nu)
        grad_pair=aog1(mu,:)*aov(nu)+aov(mu)*aog1(nu,:)
        pair_density=density(:,mu,nu)
        if(mu/=nu) pair_density=pair_density+density(:,nu,mu)
        do field=1,2
          do spin=1,2
            matrix=spin+2*(field-1)
            value(spin,field)=value(spin,field)+pair_density(matrix)*pair
            gradient(:,spin,field)=gradient(:,spin,field)+ &
              pair_density(matrix)*grad_pair
          end do
        end do
      end do
    end do
  end subroutine field_value_gradient_batch

!-----------------------------------------------------------------------------

  pure subroutine response_field_value_gradient_batch(density,aov,aog1, &
      value,gradient)
    real(fp), intent(in) :: density(:,:,:,:)
    real(fp), intent(in) :: aov(:),aog1(:,:)
    real(fp), intent(out) :: value(:,:,:),gradient(:,:,:,:)
    real(fp) :: qvalue,qgradient(3), &
      pair_ground(size(density,1)/2,2),pair_probe(size(density,1)/2,2)
    integer :: mu,nu,spin,space,nao,ncart

    nao=size(aov)
    ncart=size(density,1)/2
    if(any(shape(density)/=[2*ncart,2,nao,nao]) .or. &
       any(shape(aog1)/=[nao,3]) .or. &
       any(shape(value)/=[2,2,ncart]) .or. &
       any(shape(gradient)/=[3,2,2,ncart])) &
      error stop 'response_field_value_gradient_batch: shape mismatch'
    value=0.0_fp
    gradient=0.0_fp
    do nu=1,nao
      do mu=1,nu
        qvalue=aov(mu)*aov(nu)
        qgradient=aog1(mu,:)*aov(nu)+aov(mu)*aog1(nu,:)
        pair_ground=density(1:ncart,:,mu,nu)
        pair_probe=density(ncart+1:2*ncart,:,mu,nu)
        if(mu/=nu) then
          pair_ground=pair_ground+density(1:ncart,:,nu,mu)
          pair_probe=pair_probe+density(ncart+1:2*ncart,:,nu,mu)
        end if
        do spin=1,2
          value(spin,1,:)=value(spin,1,:)+ &
            pair_ground(:,spin)*qvalue
          value(spin,2,:)=value(spin,2,:)+ &
            pair_probe(:,spin)*qvalue
          do space=1,3
            gradient(space,spin,1,:)=gradient(space,spin,1,:)+ &
              pair_ground(:,spin)*qgradient(space)
            gradient(space,spin,2,:)=gradient(space,spin,2,:)+ &
              pair_probe(:,spin)*qgradient(space)
          end do
        end do
      end do
    end do
  end subroutine response_field_value_gradient_batch

  pure subroutine build_variables(is_gga,rho,grad_rho,prho,grad_prho, &
      drho,dgrad_rho,dprho,dgrad_prho,d2rho,d2grad_rho,d2prho, &
      d2grad_prho,ground,probe,dground,dprobe,d2ground,d2probe)
    logical, intent(in) :: is_gga
    real(fp), intent(in) :: rho(2),grad_rho(3,2),prho(2),grad_prho(3,2)
    real(fp), intent(in) :: drho(:,:,:),dgrad_rho(:,:,:,:)
    real(fp), intent(in) :: dprho(:,:,:),dgrad_prho(:,:,:,:)
    real(fp), intent(in) :: d2rho(:,:,:,:,:),d2grad_rho(:,:,:,:,:,:)
    real(fp), intent(in) :: d2prho(:,:,:,:,:),d2grad_prho(:,:,:,:,:,:)
    real(fp), intent(out) :: ground(:),probe(:),dground(:,:,:),dprobe(:,:,:)
    real(fp), intent(out) :: d2ground(:,:,:,:,:),d2probe(:,:,:,:,:)
    integer :: nat,atom_a,atom_b,cart_a,cart_b

    nat=size(drho,2)
    ground=0.0_fp;probe=0.0_fp;dground=0.0_fp;dprobe=0.0_fp
    d2ground=0.0_fp;d2probe=0.0_fp
    ground(1:2)=rho
    probe(1:2)=prho
    do atom_a=1,nat
      do cart_a=1,3
        dground(1:2,cart_a,atom_a)=drho(cart_a,atom_a,1:2)
        dprobe(1:2,cart_a,atom_a)=dprho(cart_a,atom_a,1:2)
        do atom_b=1,nat
          do cart_b=1,3
            d2ground(1:2,cart_a,cart_b,atom_a,atom_b)= &
              d2rho(cart_a,cart_b,atom_a,atom_b,1:2)
            d2probe(1:2,cart_a,cart_b,atom_a,atom_b)= &
              d2prho(cart_a,cart_b,atom_a,atom_b,1:2)
          end do
        end do
      end do
    end do
    if(.not.is_gga) return
    ground(3)=dot_product(grad_rho(:,1),grad_rho(:,1))
    ground(4)=dot_product(grad_rho(:,2),grad_rho(:,2))
    ground(5)=dot_product(grad_rho(:,1),grad_rho(:,2))
    probe(3)=2.0_fp*dot_product(grad_rho(:,1),grad_prho(:,1))
    probe(4)=2.0_fp*dot_product(grad_rho(:,2),grad_prho(:,2))
    probe(5)=dot_product(grad_rho(:,1),grad_prho(:,2))+ &
      dot_product(grad_rho(:,2),grad_prho(:,1))
    do atom_a=1,nat
      do cart_a=1,3
        call first_sigma(grad_rho,dgrad_rho(:,:,atom_a,:),grad_prho, &
          dgrad_prho(:,:,atom_a,:),cart_a,dground(3:5,cart_a,atom_a), &
          dprobe(3:5,cart_a,atom_a))
        do atom_b=1,nat
          do cart_b=1,3
            call second_sigma(grad_rho,dgrad_rho(:,:,atom_a,:), &
              dgrad_rho(:,:,atom_b,:), &
              d2grad_rho(:,:,:,atom_a,atom_b,:), &
              grad_prho,dgrad_prho(:,:,atom_a,:), &
              dgrad_prho(:,:,atom_b,:), &
              d2grad_prho(:,:,:,atom_a,atom_b,:),cart_a,cart_b, &
              d2ground(3:5,cart_a,cart_b,atom_a,atom_b), &
              d2probe(3:5,cart_a,cart_b,atom_a,atom_b))
          end do
        end do
      end do
    end do
  end subroutine build_variables

  pure subroutine first_sigma(g,dg,pg,dpg,cart,ds,dps)
    real(fp), intent(in) :: g(3,2),dg(3,3,2),pg(3,2),dpg(3,3,2)
    integer, intent(in) :: cart
    real(fp), intent(out) :: ds(3),dps(3)
    ds(1)=2.0_fp*dot_product(g(:,1),dg(:,cart,1))
    ds(2)=2.0_fp*dot_product(g(:,2),dg(:,cart,2))
    ds(3)=dot_product(dg(:,cart,1),g(:,2))+ &
      dot_product(g(:,1),dg(:,cart,2))
    dps(1)=2.0_fp*(dot_product(dg(:,cart,1),pg(:,1))+ &
      dot_product(g(:,1),dpg(:,cart,1)))
    dps(2)=2.0_fp*(dot_product(dg(:,cart,2),pg(:,2))+ &
      dot_product(g(:,2),dpg(:,cart,2)))
    dps(3)=dot_product(dg(:,cart,1),pg(:,2))+ &
      dot_product(g(:,1),dpg(:,cart,2))+ &
      dot_product(dg(:,cart,2),pg(:,1))+ &
      dot_product(g(:,2),dpg(:,cart,1))
  end subroutine first_sigma

  pure subroutine second_sigma(g,dga,dgb,d2g,pg,dpga,dpgb,d2pg, &
      cart_a,cart_b,d2s,d2ps)
    real(fp), intent(in) :: g(3,2),dga(3,3,2),dgb(3,3,2),d2g(3,3,3,2)
    real(fp), intent(in) :: pg(3,2),dpga(3,3,2),dpgb(3,3,2), &
      d2pg(3,3,3,2)
    integer, intent(in) :: cart_a,cart_b
    real(fp), intent(out) :: d2s(3),d2ps(3)
    d2s(1)=2.0_fp*(dot_product(dga(:,cart_a,1),dgb(:,cart_b,1))+ &
      dot_product(g(:,1),d2g(:,cart_a,cart_b,1)))
    d2s(2)=2.0_fp*(dot_product(dga(:,cart_a,2),dgb(:,cart_b,2))+ &
      dot_product(g(:,2),d2g(:,cart_a,cart_b,2)))
    d2s(3)=dot_product(d2g(:,cart_a,cart_b,1),g(:,2))+ &
      dot_product(dga(:,cart_a,1),dgb(:,cart_b,2))+ &
      dot_product(dgb(:,cart_b,1),dga(:,cart_a,2))+ &
      dot_product(g(:,1),d2g(:,cart_a,cart_b,2))
    d2ps(1)=2.0_fp*( &
      dot_product(d2g(:,cart_a,cart_b,1),pg(:,1))+ &
      dot_product(dga(:,cart_a,1),dpgb(:,cart_b,1))+ &
      dot_product(dgb(:,cart_b,1),dpga(:,cart_a,1))+ &
      dot_product(g(:,1),d2pg(:,cart_a,cart_b,1)))
    d2ps(2)=2.0_fp*( &
      dot_product(d2g(:,cart_a,cart_b,2),pg(:,2))+ &
      dot_product(dga(:,cart_a,2),dpgb(:,cart_b,2))+ &
      dot_product(dgb(:,cart_b,2),dpga(:,cart_a,2))+ &
      dot_product(g(:,2),d2pg(:,cart_a,cart_b,2)))
    d2ps(3)= &
      dot_product(d2g(:,cart_a,cart_b,1),pg(:,2))+ &
      dot_product(dga(:,cart_a,1),dpgb(:,cart_b,2))+ &
      dot_product(dgb(:,cart_b,1),dpga(:,cart_a,2))+ &
      dot_product(g(:,1),d2pg(:,cart_a,cart_b,2))+ &
      dot_product(d2g(:,cart_a,cart_b,2),pg(:,1))+ &
      dot_product(dga(:,cart_a,2),dpgb(:,cart_b,1))+ &
      dot_product(dgb(:,cart_b,2),dpga(:,cart_a,1))+ &
      dot_product(g(:,2),d2pg(:,cart_a,cart_b,1))
  end subroutine second_sigma

!-----------------------------------------------------------------------------

  pure subroutine build_response_variables(is_gga,rho,g,prho,pg,drho,dg, &
      dprho,dpg,rrho,rg,rprho,rpg,drrho,drg,drprho,drpg, &
      response_ground,response_probe,dresponse_ground,dresponse_probe)
    logical, intent(in) :: is_gga
    real(fp), intent(in) :: rho(2),g(3,2),prho(2),pg(3,2)
    real(fp), intent(in) :: drho(:,:,:),dg(:,:,:,:),dprho(:,:,:), &
      dpg(:,:,:,:)
    real(fp), intent(in) :: rrho(2),rg(3,2),rprho(2),rpg(3,2)
    real(fp), intent(in) :: drrho(:,:,:),drg(:,:,:,:),drprho(:,:,:), &
      drpg(:,:,:,:)
    real(fp), intent(out) :: response_ground(:),response_probe(:)
    real(fp), intent(out) :: dresponse_ground(:,:,:), &
      dresponse_probe(:,:,:)
    integer :: nat,atom,cart

    nat=size(drho,2)
    response_ground=0.0_fp;response_probe=0.0_fp
    dresponse_ground=0.0_fp;dresponse_probe=0.0_fp
    response_ground(1:2)=rrho
    response_probe(1:2)=rprho
    do atom=1,nat
      do cart=1,3
        dresponse_ground(1:2,cart,atom)=drrho(cart,atom,1:2)
        dresponse_probe(1:2,cart,atom)=drprho(cart,atom,1:2)
      end do
    end do
    if(.not.is_gga) return
    response_ground(3)=2.0_fp*dot_product(g(:,1),rg(:,1))
    response_ground(4)=2.0_fp*dot_product(g(:,2),rg(:,2))
    response_ground(5)=dot_product(rg(:,1),g(:,2))+dot_product(g(:,1),rg(:,2))
    response_probe(3)=2.0_fp*(dot_product(rg(:,1),pg(:,1))+ &
      dot_product(g(:,1),rpg(:,1)))
    response_probe(4)=2.0_fp*(dot_product(rg(:,2),pg(:,2))+ &
      dot_product(g(:,2),rpg(:,2)))
    response_probe(5)=dot_product(rg(:,1),pg(:,2))+ &
      dot_product(g(:,1),rpg(:,2))+dot_product(rg(:,2),pg(:,1))+ &
      dot_product(g(:,2),rpg(:,1))
    do atom=1,nat
      do cart=1,3
        dresponse_ground(3,cart,atom)=2.0_fp*( &
          dot_product(dg(:,cart,atom,1),rg(:,1))+ &
          dot_product(g(:,1),drg(:,cart,atom,1)))
        dresponse_ground(4,cart,atom)=2.0_fp*( &
          dot_product(dg(:,cart,atom,2),rg(:,2))+ &
          dot_product(g(:,2),drg(:,cart,atom,2)))
        dresponse_ground(5,cart,atom)= &
          dot_product(drg(:,cart,atom,1),g(:,2))+ &
          dot_product(rg(:,1),dg(:,cart,atom,2))+ &
          dot_product(dg(:,cart,atom,1),rg(:,2))+ &
          dot_product(g(:,1),drg(:,cart,atom,2))
        dresponse_probe(3,cart,atom)=2.0_fp*( &
          dot_product(drg(:,cart,atom,1),pg(:,1))+ &
          dot_product(rg(:,1),dpg(:,cart,atom,1))+ &
          dot_product(dg(:,cart,atom,1),rpg(:,1))+ &
          dot_product(g(:,1),drpg(:,cart,atom,1)))
        dresponse_probe(4,cart,atom)=2.0_fp*( &
          dot_product(drg(:,cart,atom,2),pg(:,2))+ &
          dot_product(rg(:,2),dpg(:,cart,atom,2))+ &
          dot_product(dg(:,cart,atom,2),rpg(:,2))+ &
          dot_product(g(:,2),drpg(:,cart,atom,2)))
        dresponse_probe(5,cart,atom)= &
          dot_product(drg(:,cart,atom,1),pg(:,2))+ &
          dot_product(rg(:,1),dpg(:,cart,atom,2))+ &
          dot_product(dg(:,cart,atom,1),rpg(:,2))+ &
          dot_product(g(:,1),drpg(:,cart,atom,2))+ &
          dot_product(drg(:,cart,atom,2),pg(:,1))+ &
          dot_product(rg(:,2),dpg(:,cart,atom,1))+ &
          dot_product(dg(:,cart,atom,2),rpg(:,1))+ &
          dot_product(g(:,2),drpg(:,cart,atom,1))
      end do
    end do
  end subroutine build_response_variables

!-----------------------------------------------------------------------------

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
    dr=0.0_fp;ds=0.0_fp
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

end module mod_dft_gridint_mrsf_xc_hessian
