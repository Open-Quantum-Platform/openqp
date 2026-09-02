! Pointwise analytic nuclear derivatives of spin-resolved semilocal XC Fock
! matrix elements.  This module is deliberately independent of Libxc and the
! molecular-grid driver so its chain rule can be verified with exact local
! polynomial oracles.
module mod_dft_gridint_mrsf_xc_fock_deriv_point

  use precision, only: fp

  implicit none
  private

  integer, parameter :: hmap(3,3) = reshape([1,4,6, 4,2,5, 6,5,3],[3,3])

  public :: lda_spin_fock_point_derivative
  public :: gga_spin_fock_point_derivative
  public :: spin_fock_point_derivative
  public :: moving_ao_pair_derivative
  public :: moving_density_derivative

contains

!> Differentiate a spin-resolved AO operator with independently supplied
!> scalar-pair and gradient-pair coefficients.  This is the general form
!> needed for the derivative of the spin-polarized XC-kernel action.
  pure subroutine spin_fock_point_derivative(quadrature_scale,weight, &
      dweight,v_r,coefficient,dv_r,dcoefficient,pair,grad_pair,dpair, &
      dgrad_pair,derivative,status)
    real(fp), intent(in) :: quadrature_scale,weight,dweight(:)
    real(fp), intent(in) :: v_r(2),coefficient(3,2),dv_r(:,:), &
      dcoefficient(:,:,:)
    real(fp), intent(in) :: pair,grad_pair(3),dpair(:),dgrad_pair(:,:)
    real(fp), intent(out) :: derivative(:,:)
    integer, intent(out) :: status

    real(fp) :: integrand,dintegrand
    integer :: coordinate,ncoord,spin

    ncoord=size(dweight)
    status=0
    derivative=0.0_fp
    if(ncoord<=0 .or. any(shape(dv_r)/=[2,ncoord]) .or. &
       any(shape(dcoefficient)/=[3,2,ncoord]) .or. &
       size(dpair)/=ncoord .or. any(shape(dgrad_pair)/=[3,ncoord]) .or. &
       any(shape(derivative)/=[2,ncoord])) then
      status=-1
      return
    end if
    do coordinate=1,ncoord
      do spin=1,2
        integrand=v_r(spin)*pair+dot_product(coefficient(:,spin),grad_pair)
        dintegrand=dv_r(spin,coordinate)*pair+ &
          v_r(spin)*dpair(coordinate)+ &
          dot_product(dcoefficient(:,spin,coordinate),grad_pair)+ &
          dot_product(coefficient(:,spin),dgrad_pair(:,coordinate))
        derivative(spin,coordinate)=quadrature_scale*( &
          dweight(coordinate)*integrand+weight*dintegrand)
      end do
    end do
  end subroutine spin_fock_point_derivative

!-------------------------------------------------------------------------------

!> Differentiate one spin-resolved LDA AO-matrix integrand.
!>
!> quadrature_scale excludes the normalized owner-cell weight.  Thus the
!> returned derivative contains both the partition derivative and the
!> derivative of the local potential and AO pair.
  pure subroutine lda_spin_fock_point_derivative(quadrature_scale,weight, &
      dweight,v_r,dv_r,pair,dpair,derivative,status)
    real(fp), intent(in) :: quadrature_scale,weight,dweight(:)
    real(fp), intent(in) :: v_r(2),dv_r(:,:),pair,dpair(:)
    real(fp), intent(out) :: derivative(:,:)
    integer, intent(out) :: status

    integer :: coordinate,spin,ncoord

    ncoord=size(dweight)
    status=0
    derivative=0.0_fp
    if(ncoord<=0 .or. any(shape(dv_r)/=[2,ncoord]) .or. &
       size(dpair)/=ncoord .or. any(shape(derivative)/=[2,ncoord])) then
      status=-1
      return
    end if
    do coordinate=1,ncoord
      do spin=1,2
        derivative(spin,coordinate)=quadrature_scale*( &
          dweight(coordinate)*v_r(spin)*pair+weight*( &
          dv_r(spin,coordinate)*pair+v_r(spin)*dpair(coordinate)))
      end do
    end do
  end subroutine lda_spin_fock_point_derivative

!-------------------------------------------------------------------------------

!> Differentiate one spin-resolved GGA AO-matrix integrand.
!>
!> The sigma order is (aa,bb,ab), matching xc_der1/xc_der2_contr.  For spin a,
!> c_a=2 e_sigma_aa grad(rho_a)+e_sigma_ab grad(rho_b), and analogously for b.
  pure subroutine gga_spin_fock_point_derivative(quadrature_scale,weight, &
      dweight,v_r,v_sigma,dv_r,dv_sigma,grad_rho,dgrad_rho,pair,grad_pair, &
      dpair,dgrad_pair,derivative,status)
    real(fp), intent(in) :: quadrature_scale,weight,dweight(:)
    real(fp), intent(in) :: v_r(2),v_sigma(3),dv_r(:,:),dv_sigma(:,:)
    real(fp), intent(in) :: grad_rho(3,2),dgrad_rho(:,:,:)
    real(fp), intent(in) :: pair,grad_pair(3),dpair(:),dgrad_pair(:,:)
    real(fp), intent(out) :: derivative(:,:)
    integer, intent(out) :: status

    real(fp) :: coefficient(3,2),dcoefficient(3,2)
    real(fp) :: integrand,dintegrand
    integer :: coordinate,ncoord,spin

    ncoord=size(dweight)
    status=0
    derivative=0.0_fp
    if(ncoord<=0 .or. any(shape(dv_r)/=[2,ncoord]) .or. &
       any(shape(dv_sigma)/=[3,ncoord]) .or. &
       any(shape(dgrad_rho)/=[3,2,ncoord]) .or. &
       size(dpair)/=ncoord .or. any(shape(dgrad_pair)/=[3,ncoord]) .or. &
       any(shape(derivative)/=[2,ncoord])) then
      status=-1
      return
    end if

    coefficient(:,1)=2.0_fp*v_sigma(1)*grad_rho(:,1) &
      +v_sigma(3)*grad_rho(:,2)
    coefficient(:,2)=2.0_fp*v_sigma(2)*grad_rho(:,2) &
      +v_sigma(3)*grad_rho(:,1)
    do coordinate=1,ncoord
      dcoefficient(:,1)=2.0_fp*dv_sigma(1,coordinate)*grad_rho(:,1) &
        +2.0_fp*v_sigma(1)*dgrad_rho(:,1,coordinate) &
        +dv_sigma(3,coordinate)*grad_rho(:,2) &
        +v_sigma(3)*dgrad_rho(:,2,coordinate)
      dcoefficient(:,2)=2.0_fp*dv_sigma(2,coordinate)*grad_rho(:,2) &
        +2.0_fp*v_sigma(2)*dgrad_rho(:,2,coordinate) &
        +dv_sigma(3,coordinate)*grad_rho(:,1) &
        +v_sigma(3)*dgrad_rho(:,1,coordinate)
      do spin=1,2
        integrand=v_r(spin)*pair+dot_product(coefficient(:,spin),grad_pair)
        dintegrand=dv_r(spin,coordinate)*pair &
          +v_r(spin)*dpair(coordinate) &
          +dot_product(dcoefficient(:,spin),grad_pair) &
          +dot_product(coefficient(:,spin),dgrad_pair(:,coordinate))
        derivative(spin,coordinate)=quadrature_scale*( &
          dweight(coordinate)*integrand+weight*dintegrand)
      end do
    end do
  end subroutine gga_spin_fock_point_derivative

!-------------------------------------------------------------------------------

!> Nuclear derivatives of phi_mu phi_nu and its electronic gradient at a
!> quadrature point translated rigidly with its owner atom.
!>
!> D_A phi_mu=(delta_A,owner-delta_A,center(mu))*grad(phi_mu).  This includes
!> both moving AO centers and moving atom-centred grid points and therefore
!> obeys exact rigid-translation cancellation point by point.
  pure subroutine moving_ao_pair_derivative(owner,atom_mu,atom_nu,aov_mu, &
      aov_nu,aog1_mu,aog1_nu,aog2_mu,aog2_nu,pair,grad_pair,dpair, &
      dgrad_pair,status)
    integer, intent(in) :: owner,atom_mu,atom_nu
    real(fp), intent(in) :: aov_mu,aov_nu,aog1_mu(:),aog1_nu(:)
    real(fp), intent(in) :: aog2_mu(:),aog2_nu(:)
    real(fp), intent(out) :: pair,grad_pair(3),dpair(:),dgrad_pair(:,:)
    integer, intent(out) :: status

    real(fp) :: dmu,dnu,dgmu(3),dgnu(3),factor_mu,factor_nu
    integer :: atom,cart,index,nat,ncoord,space,center_index,ncenter
    integer :: centers(3)

    ncoord=size(dpair)
    nat=ncoord/3
    status=0
    pair=0.0_fp
    grad_pair=0.0_fp
    dpair=0.0_fp
    dgrad_pair=0.0_fp
    if(ncoord<=0 .or. 3*nat/=ncoord .or. owner<1 .or. owner>nat .or. &
       atom_mu<1 .or. atom_mu>nat .or. atom_nu<1 .or. atom_nu>nat .or. &
       size(aog1_mu)/=3 .or. size(aog1_nu)/=3 .or. &
       size(aog2_mu)/=6 .or. size(aog2_nu)/=6 .or. &
       any(shape(dgrad_pair)/=[3,ncoord])) then
      status=-1
      return
    end if
    pair=aov_mu*aov_nu
    grad_pair=aog1_mu*aov_nu+aov_mu*aog1_nu
    centers(1)=owner
    ncenter=1
    if(atom_mu/=owner) then
      ncenter=ncenter+1
      centers(ncenter)=atom_mu
    end if
    if(atom_nu/=owner .and. atom_nu/=atom_mu) then
      ncenter=ncenter+1
      centers(ncenter)=atom_nu
    end if
    do center_index=1,ncenter
      atom=centers(center_index)
      factor_mu=real(merge(1,0,atom==owner)-merge(1,0,atom==atom_mu),fp)
      factor_nu=real(merge(1,0,atom==owner)-merge(1,0,atom==atom_nu),fp)
      do cart=1,3
        index=3*(atom-1)+cart
        dmu=factor_mu*aog1_mu(cart)
        dnu=factor_nu*aog1_nu(cart)
        do space=1,3
          dgmu(space)=factor_mu*aog2_mu(hmap(cart,space))
          dgnu(space)=factor_nu*aog2_nu(hmap(cart,space))
        end do
        dpair(index)=dmu*aov_nu+aov_mu*dnu
        dgrad_pair(:,index)=dgmu*aov_nu+aog1_mu*dnu &
          +dmu*aog1_nu+aov_mu*dgnu
      end do
    end do
  end subroutine moving_ao_pair_derivative

!-------------------------------------------------------------------------------

!> Moving-grid nuclear derivatives of a density and its electronic gradient.
!> No orbital, determinant, or Slater representation is introduced: density
!> matrices are the sole state-dependent input.
  pure subroutine moving_density_derivative(density,ao_atom,owner,aov,aog1, &
      aog2,need_gradient,drho,dgrad_rho,status)
    real(fp), intent(in) :: density(:,:),aov(:),aog1(:,:),aog2(:,:)
    integer, intent(in) :: ao_atom(:),owner
    logical, intent(in) :: need_gradient
    real(fp), intent(out) :: drho(:),dgrad_rho(:,:)
    integer, intent(out) :: status

    real(fp) :: pair,grad_pair(3),dpair(size(drho)),pair_density
    real(fp) :: dgrad_pair(3,size(drho))
    integer :: mu,nu,nao,ncoord,pair_status

    nao=size(aov)
    ncoord=size(drho)
    status=0
    drho=0.0_fp
    dgrad_rho=0.0_fp
    if(nao<=0 .or. any(shape(density)/=[nao,nao]) .or. &
       size(ao_atom)/=nao .or. any(shape(aog1)/=[nao,3]) .or. &
       any(shape(aog2)/=[nao,6]) .or. any(shape(dgrad_rho)/=[3,ncoord])) then
      status=-1
      return
    end if
    do nu=1,nao
      do mu=1,nu
        pair_density=density(mu,nu)
        if(mu/=nu) pair_density=pair_density+density(nu,mu)
        if(abs(pair_density)<=tiny(1.0_fp)) cycle
        call moving_ao_pair_derivative(owner,ao_atom(mu),ao_atom(nu), &
          aov(mu),aov(nu),aog1(mu,:),aog1(nu,:),aog2(mu,:),aog2(nu,:), &
          pair,grad_pair,dpair,dgrad_pair,pair_status)
        if(pair_status/=0) then
          status=-2
          return
        end if
        drho=drho+pair_density*dpair
        if(need_gradient) &
          dgrad_rho=dgrad_rho+pair_density*dgrad_pair
      end do
    end do
  end subroutine moving_density_derivative

end module mod_dft_gridint_mrsf_xc_fock_deriv_point
