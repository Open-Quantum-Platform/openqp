module mod_dft_gridint_tdgga_hessian_driver
  use precision, only: fp
  use mod_dft_gridint_tdxc_deriv, only: gga_tdxc_kernel_t, gga_fourth_from_third
  use mod_dft_gridint_tdgga_hessian, only: gga_fixed_density_quadrature_point
  use mod_dft_gridint, only: xc_engine_t, xc_consumer_t, xc_options_t, &
    run_xc, xc_der1, xc_der2_contr, xc_der3_contr, OQP_FUNTYP_GGA
  use functionals, only: functional_t
  implicit none
  private

  type, extends(xc_consumer_t) :: xc_consumer_tdgga_hess_t
    real(fp), allocatable :: hessian(:,:,:,:,:)
    real(fp), allocatable :: density_p(:,:), density_v(:,:)
    real(fp), allocatable :: atom_xyz(:,:), surface_shift(:,:)
    integer, allocatable :: ao_atom(:)
    logical, allocatable :: dummy_atom(:)
    integer :: part_fun_type = 0
    logical :: has_surface_shift = .false.
    class(functional_t), pointer :: functional => null()
  contains
    procedure :: parallel_start => tdgga_hess_parallel_start
    procedure :: parallel_stop => tdgga_hess_parallel_stop
    procedure :: update => tdgga_hess_update
    procedure :: postUpdate => tdgga_hess_post_update
    procedure :: clean => tdgga_hess_clean
  end type xc_consumer_tdgga_hess_t

  public :: tddft_gga_fixed_hessian

contains

!> Evaluate the direct fixed-density restricted-GGA contribution on the real
!> molecular grid.  Density matrices use the same spin-summed restricted
!> convention as dmatd_density_blk.
  subroutine tddft_gga_fixed_hessian(basis, mol_grid, density_r, density_p, &
      density_v, hessian, infos)
    use basis_tools, only: basis_set
    use mod_dft_molgrid, only: dft_grid_t
    use types, only: information
    type(basis_set), intent(in) :: basis
    type(dft_grid_t), target, intent(in) :: mol_grid
    real(fp), intent(in) :: density_r(:,:), density_p(:,:), density_v(:,:)
    real(fp), intent(out) :: hessian(:,:)
    type(information), target, intent(in) :: infos

    type(xc_consumer_tdgga_hess_t) :: dat
    type(xc_options_t) :: opts
    real(fp), allocatable, target :: dnorm(:,:)
    integer :: nbf, nat, i, j, atom_a, atom_b, shell, first, last

    nbf = basis%nbf
    nat = infos%mol_prop%natom
    if (any(shape(density_r) /= [nbf,nbf]) .or. &
        any(shape(density_p) /= [nbf,nbf]) .or. &
        any(shape(density_v) /= [nbf,nbf]) .or. &
        any(shape(hessian) /= [3*nat,3*nat])) &
      error stop 'tddft_gga_fixed_hessian: shape mismatch'

    allocate(dnorm(nbf,nbf),dat%density_p(nbf,nbf),dat%density_v(nbf,nbf), &
      dat%ao_atom(nbf))
    do i = 1, nbf
      dnorm(:,i) = density_r(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      ! Restricted TDDFT Lagrangian: 2*Peff and 4*X*Fxc*X.  The point
      ! function already contains the factor four for X, so apply only the
      ! remaining Peff spin degeneracy while normalizing the AO matrix.
      dat%density_p(:,i) = 2.0_fp*density_p(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
      dat%density_v(:,i) = density_v(:,i)*basis%bfnrm(:)*basis%bfnrm(i)
    end do
    dat%ao_atom = 0
    do shell = 1, basis%nshell
      first = basis%ao_offset(shell)
      last = first+basis%naos(shell)-1
      dat%ao_atom(first:last) = basis%origin(shell)
    end do
    if (any(dat%ao_atom == 0)) error stop 'tddft_gga_fixed_hessian: AO owner missing'
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
    opts%nDer = 2                 ! GGA adds one: AO values through G3
    opts%nXCDer = 3              ! fourth kernel is local FD of Libxc third
    opts%numOccAlpha = infos%mol_prop%nelec_A
    opts%numOccBeta = infos%mol_prop%nelec_B
    opts%wfAlpha => dnorm
    opts%dft_threshold = infos%dft%grid_density_cutoff
    opts%ao_threshold = infos%dft%grid_ao_threshold
    opts%ao_sparsity_ratio = 0.0_fp
    opts%molGrid => mol_grid

    call dat%pe%init(infos%mpiinfo%comm,infos%mpiinfo%usempi)
    call run_xc(opts,dat,basis)
    do atom_b = 1, nat
      do j = 1, 3
        do atom_a = 1, nat
          do i = 1, 3
            hessian(3*(atom_a-1)+i,3*(atom_b-1)+j) = &
              dat%hessian(i,j,atom_a,atom_b,1)
          end do
        end do
      end do
    end do
    call dat%clean()
  end subroutine tddft_gga_fixed_hessian

  subroutine tdgga_hess_parallel_start(self,xce,nthreads)
    class(xc_consumer_tdgga_hess_t), target, intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer, intent(in) :: nthreads
    if (allocated(self%hessian)) deallocate(self%hessian)
    allocate(self%hessian(3,3,xce%numAtoms,xce%numAtoms,nthreads),source=0.0_fp)
  end subroutine tdgga_hess_parallel_start

  subroutine tdgga_hess_parallel_stop(self)
    class(xc_consumer_tdgga_hess_t), intent(inout) :: self
    if (size(self%hessian,5) > 1) &
      self%hessian(:,:,:,:,1) = sum(self%hessian,dim=5)
    call self%pe%allreduce(self%hessian(:,:,:,:,1), &
      size(self%hessian(:,:,:,:,1)))
  end subroutine tdgga_hess_parallel_stop

  subroutine tdgga_hess_clean(self)
    class(xc_consumer_tdgga_hess_t), intent(inout) :: self
    if (allocated(self%hessian)) deallocate(self%hessian)
    if (allocated(self%density_p)) deallocate(self%density_p)
    if (allocated(self%density_v)) deallocate(self%density_v)
    if (allocated(self%atom_xyz)) deallocate(self%atom_xyz)
    if (allocated(self%surface_shift)) deallocate(self%surface_shift)
    if (allocated(self%ao_atom)) deallocate(self%ao_atom)
    if (allocated(self%dummy_atom)) deallocate(self%dummy_atom)
    nullify(self%functional)
  end subroutine tdgga_hess_clean

  subroutine tdgga_hess_post_update(self,xce,mythread)
    class(xc_consumer_tdgga_hess_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread
    ! All work is accumulated directly by update.
  end subroutine tdgga_hess_post_update

  subroutine tdgga_hess_update(self,xce,mythread)
    class(xc_consumer_tdgga_hess_t), intent(inout) :: self
    class(xc_engine_t), intent(in) :: xce
    integer :: mythread

    type(gga_tdxc_kernel_t), allocatable :: kernels(:)
    real(fp), allocatable :: p(:,:), v(:,:)
    integer, allocatable :: atoms(:)
    integer :: i, n, status

    n = xce%numAOs_p
    allocate(p(n,n),v(n,n),atoms(n),kernels(xce%numPts))
    if (xce%skip_p) then
      p = self%density_p
      v = self%density_v
      atoms = self%ao_atom
    else
      p = self%density_p(xce%indices_p(1:n),xce%indices_p(1:n))
      v = self%density_v(xce%indices_p(1:n),xce%indices_p(1:n))
      atoms = self%ao_atom(xce%indices_p(1:n))
    end if
    call build_unweighted_gga_kernels(self,xce,kernels)
    do i = 1, xce%numPts
      call gga_fixed_density_quadrature_point(kernels(i),xce%wfAlpha_p,p,v, &
        atoms,xce%aoV(:,i),xce%aoG1(:,i,:),xce%aoG2(:,i,:),xce%aoG3(:,i,:), &
        self%atom_xyz,xce%xyzw(i,1:3),xce%currAtom,self%dummy_atom, &
        self%part_fun_type,self%has_surface_shift,self%surface_shift, &
        xce%xyzw(i,4),self%hessian(:,:,:,:,mythread),status)
      if (status /= 0) error stop 'tdgga_hess_update: partition derivative failure'
    end do
  end subroutine tdgga_hess_update

!> Current Libxc arrays are quadrature weighted.  Recover the unweighted
!> first-through-third kernel and evaluate the fourth derivatives by the same
!> local symmetric perturbation of third derivatives used by GAMESS TDHXFD4.
  subroutine build_unweighted_gga_kernels(self,xce,kernels)
    class(xc_consumer_tdgga_hess_t) :: self
    class(xc_engine_t) :: xce
    type(gga_tdxc_kernel_t), intent(out) :: kernels(:)
    type(xc_engine_t) :: probe
    type(gga_tdxc_kernel_t) :: krp, krm, ksp, ksm
    real(fp), allocatable :: ones(:)
    real(fp), parameter :: eps_r=1.0e-4_fp, eps_s=1.0e-4_fp
    real(fp) :: rho, sigma, wt
    integer :: i, n, irp, irm, isp, ism

    n = xce%numPts
    if (size(kernels) /= n) error stop 'build_unweighted_gga_kernels: shape mismatch'
    allocate(probe%xclib)
    probe%funTyp = OQP_FUNTYP_GGA
    call probe%xclib%init(.true.,.false.,.false.,.false.,4*n,3)
    call probe%xclib%setPts(4*n)
    probe%xclib%drho = 0.0_fp
    probe%xclib%tau = 0.0_fp
    probe%xclib%lapl = 0.0_fp
    do i = 1, n
      irp=i; irm=n+i; isp=2*n+i; ism=3*n+i
      probe%xclib%rho(:,irp)=xce%xclib%rho(:,i)*(1.0_fp+eps_r)
      probe%xclib%rho(:,irm)=xce%xclib%rho(:,i)*(1.0_fp-eps_r)
      probe%xclib%rho(:,isp)=xce%xclib%rho(:,i)
      probe%xclib%rho(:,ism)=xce%xclib%rho(:,i)
      probe%xclib%sig(:,irp)=xce%xclib%sig(:,i)
      probe%xclib%sig(:,irm)=xce%xclib%sig(:,i)
      probe%xclib%sig(:,isp)=xce%xclib%sig(:,i)*(1.0_fp+eps_s)
      probe%xclib%sig(:,ism)=xce%xclib%sig(:,i)*(1.0_fp-eps_s)
    end do
    allocate(ones(4*n),source=1.0_fp)
    call probe%xclib%compute(self%functional,ones)
    do i = 1, n
      wt=xce%xyzw(i,4)
      call restricted_gga_kernel(xce,i,kernels(i))
      call scale_kernel(kernels(i),1.0_fp/wt)
      call restricted_gga_kernel(probe,i,krp)
      call restricted_gga_kernel(probe,n+i,krm)
      call restricted_gga_kernel(probe,2*n+i,ksp)
      call restricted_gga_kernel(probe,3*n+i,ksm)
      rho=sum(xce%xclib%rho(:,i))
      sigma=4.0_fp*xce%xclib%sig(xce%xclib%ids%ga,i)
      call gga_fourth_from_third(rho,sigma,krp,krm,ksp,ksm, &
        eps_r,eps_s,1.0e-14_fp,kernels(i))
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
    real(fp), parameter :: cr(2)=[0.5_fp,0.5_fp], cs(3)=[0.25_fp,0.25_fp,0.25_fp]
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

end module mod_dft_gridint_tdgga_hessian_driver
