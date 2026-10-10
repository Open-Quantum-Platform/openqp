module mod_dft_gridint_tdgga_response_driver

  use precision, only: fp
  use mod_dft_gridint_tdxc_deriv, only: gga_tdxc_kernel_t, &
    gga_tdxc_variables_t, gga_fourth_from_third
  use mod_dft_gridint_tdgga_response, only: &
    gga_build_response_direction, &
    gga_build_response_direction_derivative, &
    gga_weighted_response_row, gga_add_owner_motion_first, &
    gga_density_nuclear_point_first
  use mod_dft_partition_hessian, only: partition_weight_nuclear_derivatives
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t, xc_options_t, &
    run_xc, xc_der1, xc_der2_contr, xc_der3_contr, OQP_FUNTYP_GGA
  use functionals, only: functional_t

  implicit none
  private

  type, extends(xc_consumer_t) :: xc_consumer_tdgga_response_t
    real(fp), allocatable :: rows(:,:,:,:)
    real(fp), allocatable :: density_p(:,:), density_v(:,:)
    real(fp), allocatable :: density_d(:,:,:), density_u(:,:,:)
    real(fp), allocatable :: atom_xyz(:,:), surface_shift(:,:)
    integer, allocatable :: ao_atom(:)
    logical, allocatable :: dummy_atom(:)
    integer :: part_fun_type = 0
    logical :: has_surface_shift = .false.
    class(functional_t), pointer :: functional => null()
  contains
    procedure :: parallel_start => response_parallel_start
    procedure :: parallel_stop => response_parallel_stop
    procedure :: update => response_update
    procedure :: postUpdate => response_post_update
    procedure :: clean => response_clean
  end type xc_consumer_tdgga_response_t

  public :: tddft_gga_response_rows

contains

!> Direct restricted-GGA analogue of GAMESS TDHXR1G/TDHXRPG.
  subroutine tddft_gga_response_rows(basis, mol_grid, density_r, density_p, &
      density_v, density_d, density_u, rows, infos)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information
    type(basis_set), intent(in) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: density_r(:,:), density_p(:,:), density_v(:,:)
    real(fp), intent(in) :: density_d(:,:,:), density_u(:,:,:)
    real(fp), intent(out) :: rows(:,:)
    type(information), target, intent(in) :: infos

    type(xc_consumer_tdgga_response_t) :: dat
    type(xc_options_t) :: opts
    real(fp), allocatable, target :: dnorm(:,:)
    integer :: nbf, nat, ncart, i, k, atom, cart, shell, first, last

    nbf = basis%nbf
    nat = infos%mol_prop%natom
    ncart = 3*nat
    if (any(shape(density_r) /= [nbf,nbf]) .or. &
        any(shape(density_p) /= [nbf,nbf]) .or. &
        any(shape(density_v) /= [nbf,nbf]) .or. &
        any(shape(density_d) /= [nbf,nbf,ncart]) .or. &
        any(shape(density_u) /= [nbf,nbf,ncart]) .or. &
        any(shape(rows) /= [ncart,ncart])) &
      error stop 'tddft_gga_response_rows: shape mismatch'

    allocate(dnorm(nbf,nbf),dat%density_p(nbf,nbf), &
      dat%density_v(nbf,nbf),dat%density_d(nbf,nbf,ncart), &
      dat%density_u(nbf,nbf,ncart),dat%ao_atom(nbf))
    do i = 1, nbf
      dnorm(:,i) = density_r(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      ! OQP_TD_P is one half of GAMESS PEFF.  D, PV, dD_K, and dPV_K
      ! already use the spin-summed/DAF-943 conventions required here.
      dat%density_p(:,i) = 2.0_fp*density_p(:,i) &
        *basis%bfnrm(:)*basis%bfnrm(i)
      dat%density_v(:,i) = density_v(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      do k = 1, ncart
        dat%density_d(:,i,k) = density_d(:,i,k) &
          *basis%bfnrm(:)*basis%bfnrm(i)
        dat%density_u(:,i,k) = density_u(:,i,k) &
          *basis%bfnrm(:)*basis%bfnrm(i)
      end do
    end do
    dat%ao_atom = 0
    do shell = 1, basis%nshell
      first = basis%ao_offset(shell)
      last = first+basis%naos(shell)-1
      dat%ao_atom(first:last) = basis%origin(shell)
    end do
    if (any(dat%ao_atom == 0)) &
      error stop 'tddft_gga_response_rows: AO owner missing'
    dat%atom_xyz = infos%atoms%xyz
    dat%dummy_atom = mol_grid%dummyAtom
    dat%surface_shift = mol_grid%surfaceShift
    dat%part_fun_type = mol_grid%partFunType
    dat%has_surface_shift = mol_grid%hasSurfaceShift
    dat%functional => infos%functional

    opts%isGGA = .true.
    opts%needTau = .false.
    opts%functional => infos%functional
    opts%hasBeta = .false.
    opts%isWFVecs = .false.
    opts%numAOs = nbf
    opts%maxPts = mol_grid%maxSlicePts
    opts%limPts = mol_grid%maxNRadTimesNAng
    opts%numAtoms = nat
    opts%maxAngMom = basis%mxam
    opts%nDer = 1
    opts%nXCDer = 3
    opts%numOccAlpha = infos%mol_prop%nelec_A
    opts%numOccBeta = infos%mol_prop%nelec_B
    opts%wfAlpha => dnorm
    opts%dft_threshold = infos%dft%grid_density_cutoff
    opts%ao_threshold = infos%dft%grid_ao_threshold
    opts%ao_sparsity_ratio = 0.0_fp
    opts%molGrid => mol_grid

    call dat%pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    call run_xc(opts,dat,basis)
    rows = 0.0_fp
    do k = 1, ncart
      do atom = 1, nat
        do cart = 1, 3
          rows(3*(atom-1)+cart,k) = dat%rows(cart,atom,k,1)
        end do
      end do
    end do
    call dat%clean()
  end subroutine tddft_gga_response_rows

  subroutine response_parallel_start(self,xce,nthreads)
    class(xc_consumer_tdgga_response_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    if (allocated(self%rows)) deallocate(self%rows)
    allocate(self%rows(3,xce%numAtoms,3*xce%numAtoms,nthreads),source=0.0_fp)
  end subroutine response_parallel_start

  subroutine response_parallel_stop(self)
    class(xc_consumer_tdgga_response_t), intent(inout) :: self
    if (size(self%rows,4) > 1) self%rows(:,:,:,1) = sum(self%rows,dim=4)
    call self%pe%allreduce(self%rows(:,:,:,1),size(self%rows(:,:,:,1)))
  end subroutine response_parallel_stop

  subroutine response_post_update(self,xce,mythread)
    class(xc_consumer_tdgga_response_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
  end subroutine response_post_update

  subroutine response_clean(self)
    class(xc_consumer_tdgga_response_t), intent(inout) :: self
    if (allocated(self%rows)) deallocate(self%rows)
    if (allocated(self%density_p)) deallocate(self%density_p)
    if (allocated(self%density_v)) deallocate(self%density_v)
    if (allocated(self%density_d)) deallocate(self%density_d)
    if (allocated(self%density_u)) deallocate(self%density_u)
    if (allocated(self%atom_xyz)) deallocate(self%atom_xyz)
    if (allocated(self%surface_shift)) deallocate(self%surface_shift)
    if (allocated(self%ao_atom)) deallocate(self%ao_atom)
    if (allocated(self%dummy_atom)) deallocate(self%dummy_atom)
    nullify(self%functional)
  end subroutine response_clean

  subroutine response_update(self,xce,mythread)
    class(xc_consumer_tdgga_response_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    type(gga_tdxc_kernel_t), allocatable :: kernels(:)
    real(fp), allocatable :: p(:,:), v(:,:), dd(:,:,:), du(:,:,:)
    integer, allocatable :: atoms(:)
    integer :: i, n, ncart, status

    n = xce%numAOs_p
    ncart = 3*xce%numAtoms
    allocate(p(n,n),v(n,n),dd(n,n,ncart),du(n,n,ncart),atoms(n), &
      kernels(xce%numPts))
    if (xce%skip_p) then
      p = self%density_p
      v = self%density_v
      dd = self%density_d
      du = self%density_u
      atoms = self%ao_atom
    else
      p = self%density_p(xce%indices_p(1:n),xce%indices_p(1:n))
      v = self%density_v(xce%indices_p(1:n),xce%indices_p(1:n))
      dd = self%density_d(xce%indices_p(1:n),xce%indices_p(1:n),:)
      du = self%density_u(xce%indices_p(1:n),xce%indices_p(1:n),:)
      atoms = self%ao_atom(xce%indices_p(1:n))
    end if
    call build_unweighted_gga_kernels(self,xce,kernels)
    do i = 1, xce%numPts
      call response_quadrature_point(kernels(i),xce%wfAlpha_p,p,v,dd,du, &
        atoms,xce%aoV(:,i),xce%aoG1(:,i,:),xce%aoG2(:,i,:), &
        self%atom_xyz,xce%xyzw(i,1:3),xce%currAtom,self%dummy_atom, &
        self%part_fun_type,self%has_surface_shift,self%surface_shift, &
        xce%xyzw(i,4),self%rows(:,:,:,mythread),status)
      if (status /= 0) &
        error stop 'response_update: partition derivative failure'
    end do
  end subroutine response_update

  subroutine response_quadrature_point(kernel,density_r,density_p,density_v, &
      density_d,density_u,ao_atom,aov,aog1,aog2,atom_xyz,point,owner, &
      dummy_atom,part_fun_type,has_surface_shift,surface_shift,finite_weight, &
      rows,status)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    real(fp), intent(in) :: density_r(:,:), density_p(:,:), density_v(:,:)
    real(fp), intent(in) :: density_d(:,:,:), density_u(:,:,:)
    integer, intent(in) :: ao_atom(:), owner, part_fun_type
    real(fp), intent(in) :: aov(:), aog1(:,:), aog2(:,:)
    real(fp), intent(in) :: atom_xyz(:,:), point(3), surface_shift(:,:)
    logical, intent(in) :: dummy_atom(:), has_surface_shift
    real(fp), intent(in) :: finite_weight
    real(fp), intent(inout) :: rows(:,:,:)
    integer, intent(out) :: status

    integer :: nat, ncart, k, atom, cart
    real(fp), parameter :: rho_cutoff=1.0e-8_fp
    real(fp) :: r,p,v,rd,ru,gr(3),gp(3),gv(3),gd(3),gu(3)
    real(fp) :: duk(7), dduk1(7), qscale
    real(fp), allocatable :: fr(:,:),fp_fixed(:,:),fv(:,:),fd(:,:),fu(:,:)
    real(fp), allocatable :: fgr(:,:,:),fgp(:,:,:),fgv(:,:,:)
    real(fp), allocatable :: fgd(:,:,:),fgu(:,:,:)
    real(fp), allocatable :: dr(:,:),dp(:,:),dv(:,:),dd(:,:),duf(:,:)
    real(fp), allocatable :: dgr(:,:,:),dgp(:,:,:),dgv(:,:,:)
    real(fp), allocatable :: dgd(:,:,:),dgu(:,:,:)
    real(fp), allocatable :: weights(:),dweights(:,:,:),d2weights(:,:,:,:,:)
    real(fp), allocatable :: base_du(:,:,:),response_ddu(:,:,:),row(:,:)
    type(gga_tdxc_variables_t) :: u

    nat = size(atom_xyz,2)
    ncart = 3*nat
    status = 0
    allocate(fr(3,nat),fp_fixed(3,nat),fv(3,nat),fd(3,nat),fu(3,nat), &
      fgr(3,3,nat),fgp(3,3,nat),fgv(3,3,nat),fgd(3,3,nat), &
      fgu(3,3,nat),dr(3,nat),dp(3,nat),dv(3,nat),dd(3,nat), &
      duf(3,nat),dgr(3,3,nat),dgp(3,3,nat),dgv(3,3,nat), &
      dgd(3,3,nat),dgu(3,3,nat),weights(nat),dweights(3,nat,nat), &
      d2weights(3,nat,3,nat,nat),base_du(7,3,nat), &
      response_ddu(7,3,nat),row(3,nat))

    call field_value_gradient(density_r,aov,aog1,r,gr)
    if (r < rho_cutoff) return
    call field_value_gradient(density_p,aov,aog1,p,gp)
    call field_value_gradient(density_v,aov,aog1,v,gv)
    call first_field_derivatives(density_r,fr,fgr,dr,dgr)
    call first_field_derivatives(density_p,fp_fixed,fgp,dp,dgp)
    call first_field_derivatives(density_v,fv,fgv,dv,dgv)

    u%r=r; u%s=dot_product(gr,gr); u%p=p; u%pg=dot_product(gr,gp)
    u%v=v; u%gv=dot_product(gr,gv); u%w=dot_product(gv,gv)
    do atom = 1, nat
      do cart = 1, 3
        base_du(1,cart,atom)=dr(cart,atom)
        base_du(2,cart,atom)=2.0_fp*dot_product(gr,dgr(:,cart,atom))
        base_du(3,cart,atom)=dp(cart,atom)
        base_du(4,cart,atom)=dot_product(dgr(:,cart,atom),gp) &
          +dot_product(gr,dgp(:,cart,atom))
        base_du(5,cart,atom)=dv(cart,atom)
        base_du(6,cart,atom)=dot_product(dgr(:,cart,atom),gv) &
          +dot_product(gr,dgv(:,cart,atom))
        base_du(7,cart,atom)=2.0_fp*dot_product(gv,dgv(:,cart,atom))
      end do
    end do

    call partition_weight_nuclear_derivatives(atom_xyz,point,owner,dummy_atom, &
      part_fun_type,has_surface_shift,surface_shift,weights,dweights, &
      d2weights,status)
    if (status /= 0) return
    if (weights(owner) <= sqrt(tiny(1.0_fp))) return
    qscale = finite_weight/weights(owner)

    do k = 1, ncart
      call field_value_gradient(density_d(:,:,k),aov,aog1,rd,gd)
      call field_value_gradient(density_u(:,:,k),aov,aog1,ru,gu)
      call first_field_derivatives(density_d(:,:,k),fd,fgd,dd,dgd)
      call first_field_derivatives(density_u(:,:,k),fu,fgu,duf,dgu)
      call gga_build_response_direction(gr,gp,gv,rd,gd,ru,gu,duk)
      do atom = 1, nat
        do cart = 1, 3
          call gga_build_response_direction_derivative(gr,gp,gv,gd,gu, &
            dgr(:,cart,atom),dgp(:,cart,atom),dgv(:,cart,atom), &
            dd(cart,atom),dgd(:,cart,atom),duf(cart,atom), &
            dgu(:,cart,atom),dduk1)
          response_ddu(:,cart,atom)=dduk1
        end do
      end do
      row=0.0_fp
      call gga_weighted_response_row(kernel,u,base_du,duk,response_ddu, &
        qscale,weights(owner),dweights(:,:,owner),row)
      rows(:,:,k)=rows(:,:,k)+row
    end do

  contains

    subroutine first_field_derivatives(density,fixed_r,fixed_g,r1,g1)
      real(fp), intent(in) :: density(:,:)
      real(fp), intent(out) :: fixed_r(:,:),fixed_g(:,:,:),r1(:,:),g1(:,:,:)
      call gga_density_nuclear_point_first(density,ao_atom,aov,aog1,aog2, &
        fixed_r,fixed_g)
      call gga_add_owner_motion_first(owner,fixed_r,fixed_g,r1,g1)
    end subroutine first_field_derivatives

    subroutine field_value_gradient(density,values,gradients,value,gradient)
      real(fp), intent(in) :: density(:,:),values(:),gradients(:,:)
      real(fp), intent(out) :: value,gradient(3)
      integer :: mu,nu
      value=0.0_fp; gradient=0.0_fp
      do nu=1,size(values)
        do mu=1,size(values)
          value=value+density(mu,nu)*values(mu)*values(nu)
          gradient=gradient+density(mu,nu)*(gradients(mu,:)*values(nu) &
            +values(mu)*gradients(nu,:))
        end do
      end do
    end subroutine field_value_gradient
  end subroutine response_quadrature_point

!> Recover unweighted total-density kernels and obtain the five fourth
!> derivatives by the TDHXFD4 local rho/sigma finite differences.
  subroutine build_unweighted_gga_kernels(self,xce,kernels)
    class(xc_consumer_tdgga_response_t) :: self
    class(xc_engine_t) :: xce
    type(gga_tdxc_kernel_t), intent(out) :: kernels(:)
    type(xc_engine_t) :: probe
    type(gga_tdxc_kernel_t) :: krp,krm,ksp,ksm
    real(fp), allocatable :: ones(:)
    real(fp), parameter :: eps_r=1.0e-3_fp,eps_s=1.0e-3_fp
    real(fp), parameter :: rho_cutoff=1.0e-8_fp
    real(fp) :: rho,sigma,wt
    integer :: i,n

    n=xce%numPts
    allocate(probe%xclib)
    probe%funTyp=OQP_FUNTYP_GGA
    call probe%xclib%init(.true.,.false.,.false.,.false.,4*n,3)
    call probe%xclib%setPts(4*n)
    probe%xclib%drho=0.0_fp; probe%xclib%tau=0.0_fp; probe%xclib%lapl=0.0_fp
    do i=1,n
      probe%xclib%rho(:,i)=xce%xclib%rho(:,i)*(1.0_fp+eps_r)
      probe%xclib%rho(:,n+i)=xce%xclib%rho(:,i)*(1.0_fp-eps_r)
      probe%xclib%rho(:,2*n+i)=xce%xclib%rho(:,i)
      probe%xclib%rho(:,3*n+i)=xce%xclib%rho(:,i)
      probe%xclib%sig(:,i)=xce%xclib%sig(:,i)
      probe%xclib%sig(:,n+i)=xce%xclib%sig(:,i)
      probe%xclib%sig(:,2*n+i)=xce%xclib%sig(:,i)*(1.0_fp+eps_s)
      probe%xclib%sig(:,3*n+i)=xce%xclib%sig(:,i)*(1.0_fp-eps_s)
    end do
    allocate(ones(4*n),source=1.0_fp)
    call probe%xclib%compute(self%functional,ones)
    do i=1,n
      wt=xce%xyzw(i,4)
      call restricted_gga_kernel(xce,i,kernels(i)); call scale_kernel(kernels(i),1.0_fp/wt)
      call restricted_gga_kernel(probe,i,krp)
      call restricted_gga_kernel(probe,n+i,krm)
      call restricted_gga_kernel(probe,2*n+i,ksp)
      call restricted_gga_kernel(probe,3*n+i,ksm)
      rho=sum(xce%xclib%rho(:,i)); sigma=4.0_fp*xce%xclib%sig(xce%xclib%ids%ga,i)
      call gga_fourth_from_third(rho,sigma,krp,krm,ksp,ksm, &
        eps_r,eps_s,rho_cutoff,kernels(i))
    end do
  end subroutine build_unweighted_gga_kernels

  subroutine scale_kernel(k,s)
    type(gga_tdxc_kernel_t), intent(inout) :: k
    real(fp), intent(in) :: s
    k%fr=k%fr*s; k%fs=k%fs*s; k%frr=k%frr*s; k%frs=k%frs*s
    k%fss=k%fss*s; k%frrr=k%frrr*s; k%frrs=k%frrs*s
    k%frss=k%frss*s; k%fsss=k%fsss*s
  end subroutine scale_kernel

  subroutine restricted_gga_kernel(xce,ipt,kernel)
    class(xc_engine_t) :: xce
    integer, intent(in) :: ipt
    type(gga_tdxc_kernel_t), intent(out) :: kernel
    real(fp), parameter :: cr(2)=[0.5_fp,0.5_fp]
    real(fp), parameter :: cs(3)=[0.25_fp,0.25_fp,0.25_fp]
    real(fp) :: zr(2),zs(3),zt(2),dr(2),ds(3),dt(2),fr(2),fs(3),ft(2)
    real(fp) :: ffs(3),gr(2),gs(3),gt(2)
    zr=0;zs=0;zt=0;kernel=gga_tdxc_kernel_t()
    call xc_der1(xce,.true.,ipt,dr,ds,dt)
    kernel%fr=dot_product(cr,dr);kernel%fs=dot_product(cs,ds)
    call xc_der2_contr(xce,.true.,ipt,cr,zs,zt,fr,fs,ft)
    kernel%frr=dot_product(cr,fr);kernel%frs=dot_product(cs,fs)
    call xc_der2_contr(xce,.true.,ipt,zr,cs,zt,fr,fs,ft)
    kernel%fss=dot_product(cs,fs)
    call xc_der3_contr(xce,ipt,cr,zs,zt,zs,ffs,gr,gs,gt)
    kernel%frrr=dot_product(cr,gr);kernel%frrs=dot_product(cs,gs)
    call xc_der3_contr(xce,ipt,zr,cs,zt,zs,ffs,gr,gs,gt)
    kernel%frss=dot_product(cr,gr);kernel%fsss=dot_product(cs,gs)
  end subroutine restricted_gga_kernel

end module mod_dft_gridint_tdgga_response_driver
