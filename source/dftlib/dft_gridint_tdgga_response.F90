! Pointwise restricted-GGA TDDFT response-row chain rule.
module mod_dft_gridint_tdgga_response

  use precision, only: fp
  use mod_dft_gridint_tdxc_deriv, only: gga_tdxc_kernel_t, &
    gga_tdxc_variables_t, gga_tdxc_fixed_derivatives

  implicit none
  private

  integer, parameter :: hmap(3,3) = reshape([1,4,6, 4,2,5, 6,5,3], [3,3])

  public :: gga_build_response_direction
  public :: gga_build_response_direction_derivative
  public :: gga_weighted_response_row
  public :: gga_add_owner_motion_first
  public :: gga_density_nuclear_point_first

contains

!> Build the TDHXRPG response direction in
!> u=(rho,sigma,rho_P,Gamma_P,rho_v,Gamma_v,|grad rho_v|2).
!> Peff is frozen in this direction.
  pure subroutine gga_build_response_direction(grad_r, grad_p, grad_v, &
      rho_d, grad_d, rho_u, grad_u, duk)
    real(fp), intent(in) :: grad_r(3), grad_p(3), grad_v(3)
    real(fp), intent(in) :: rho_d, grad_d(3), rho_u, grad_u(3)
    real(fp), intent(out) :: duk(7)

    duk(1) = rho_d
    duk(2) = 2.0_fp*dot_product(grad_r,grad_d)
    duk(3) = 0.0_fp
    duk(4) = dot_product(grad_d,grad_p)
    duk(5) = rho_u
    duk(6) = dot_product(grad_d,grad_v)+dot_product(grad_r,grad_u)
    duk(7) = 2.0_fp*dot_product(grad_v,grad_u)
  end subroutine gga_build_response_direction

!> Differentiate the response direction with respect to one skeleton
!> coordinate while all five AO density matrices remain frozen.
  pure subroutine gga_build_response_direction_derivative(grad_r, grad_p, &
      grad_v, grad_d, grad_u, dgrad_r, dgrad_p, dgrad_v, d_rho_d, &
      dgrad_d, d_rho_u, dgrad_u, dduk)
    real(fp), intent(in) :: grad_r(3), grad_p(3), grad_v(3)
    real(fp), intent(in) :: grad_d(3), grad_u(3)
    real(fp), intent(in) :: dgrad_r(3), dgrad_p(3), dgrad_v(3)
    real(fp), intent(in) :: d_rho_d, dgrad_d(3), d_rho_u, dgrad_u(3)
    real(fp), intent(out) :: dduk(7)

    dduk(1) = d_rho_d
    dduk(2) = 2.0_fp*(dot_product(dgrad_r,grad_d) &
      +dot_product(grad_r,dgrad_d))
    dduk(3) = 0.0_fp
    dduk(4) = dot_product(dgrad_d,grad_p) &
      +dot_product(grad_d,dgrad_p)
    dduk(5) = d_rho_u
    dduk(6) = dot_product(dgrad_d,grad_v) &
      +dot_product(grad_d,dgrad_v) &
      +dot_product(dgrad_r,grad_u) &
      +dot_product(grad_r,dgrad_u)
    dduk(7) = 2.0_fp*(dot_product(dgrad_v,grad_u) &
      +dot_product(grad_v,dgrad_u))
  end subroutine gga_build_response_direction_derivative

!> Add one response direction at one grid point to a complete skeleton row.
!> quadrature_scale is the radial/angular factor, excluding the normalized
!> owner-cell partition weight.
  subroutine gga_weighted_response_row(kernel, u, du, duk, dduk, &
      quadrature_scale, weight, dweight, row)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    type(gga_tdxc_variables_t), intent(in) :: u
    real(fp), intent(in) :: du(:,:,:), duk(7), dduk(:,:,:)
    real(fp), intent(in) :: quadrature_scale, weight
    real(fp), intent(in) :: dweight(:,:)
    real(fp), intent(inout) :: row(:,:)

    real(fp) :: fa(7), fab(7,7), gb(7), fr, dfr
    integer :: nat, atom, cart

    nat = size(dweight,2)
    if (any(shape(du) /= [7,3,nat]) .or. &
        any(shape(dduk) /= [7,3,nat]) .or. &
        any(shape(row) /= [3,nat])) &
      error stop 'gga_weighted_response_row: shape mismatch'
    call gga_tdxc_fixed_derivatives(kernel,u,fa,fab)
    fr = dot_product(fa,duk)
    gb = matmul(transpose(fab),duk)
    do atom = 1, nat
      do cart = 1, 3
        dfr = dot_product(fa,dduk(:,cart,atom)) &
          +dot_product(gb,du(:,cart,atom))
        row(cart,atom) = row(cart,atom)+quadrature_scale*( &
          dweight(cart,atom)*fr+weight*dfr)
      end do
    end do
  end subroutine gga_weighted_response_row

!> Convert fixed-laboratory-grid first nuclear derivatives to derivatives at
!> a point translated rigidly with its owner atom.
  pure subroutine gga_add_owner_motion_first(owner, fixed_d1, &
      fixed_grad_d1, d1, grad_d1)
    integer, intent(in) :: owner
    real(fp), intent(in) :: fixed_d1(:,:), fixed_grad_d1(:,:,:)
    real(fp), intent(out) :: d1(:,:), grad_d1(:,:,:)

    integer :: nat, cart

    nat = size(fixed_d1,2)
    if (owner < 1 .or. owner > nat .or. &
        any(shape(fixed_d1) /= [3,nat]) .or. &
        any(shape(fixed_grad_d1) /= [3,3,nat]) .or. &
        any(shape(d1) /= [3,nat]) .or. &
        any(shape(grad_d1) /= [3,3,nat])) &
      error stop 'gga_add_owner_motion_first: shape mismatch'
    d1 = fixed_d1
    grad_d1 = fixed_grad_d1
    do cart = 1, 3
      d1(cart,owner) = d1(cart,owner)-sum(fixed_d1(cart,:))
      grad_d1(:,cart,owner) = grad_d1(:,cart,owner) &
        -sum(fixed_grad_d1(:,cart,:),dim=2)
    end do
  end subroutine gga_add_owner_motion_first

!> Fixed-grid first nuclear derivatives of rho and grad(rho).  This is the
!> first-order subset of gga_density_nuclear_point and needs AO values only
!> through G2, as does GAMESS TDHXRDG.
  pure subroutine gga_density_nuclear_point_first(density,ao_atom,aov,aog1, &
      aog2,drho,dgrho)
    real(fp), intent(in) :: density(:,:),aov(:),aog1(:,:),aog2(:,:)
    integer, intent(in) :: ao_atom(:)
    real(fp), intent(out) :: drho(:,:),dgrho(:,:,:)
    integer :: nao,nat,mu,nu,atom,cart,space
    real(fp) :: p,qmu,qnu,qgmu,qgnu

    nao=size(aov);nat=size(drho,2)
    if (any(shape(density) /= [nao,nao]) .or. size(ao_atom) /= nao .or. &
        any(shape(aog1) /= [nao,3]) .or. any(shape(aog2) /= [nao,6]) .or. &
        any(shape(drho) /= [3,nat]) .or. &
        any(shape(dgrho) /= [3,3,nat])) &
      error stop 'gga_density_nuclear_point_first: shape mismatch'
    drho=0.0_fp;dgrho=0.0_fp
    do nu=1,nao
      do mu=1,nao
        p=density(mu,nu)
        if (p == 0.0_fp) cycle
        do atom=1,nat
          do cart=1,3
            qmu=0.0_fp;qnu=0.0_fp
            if (ao_atom(mu) == atom) qmu=-aog1(mu,cart)
            if (ao_atom(nu) == atom) qnu=-aog1(nu,cart)
            drho(cart,atom)=drho(cart,atom) &
              +p*(qmu*aov(nu)+aov(mu)*qnu)
            do space=1,3
              qgmu=0.0_fp;qgnu=0.0_fp
              if (ao_atom(mu) == atom) qgmu=-aog2(mu,hmap(cart,space))
              if (ao_atom(nu) == atom) qgnu=-aog2(nu,hmap(cart,space))
              dgrho(space,cart,atom)=dgrho(space,cart,atom)+p*( &
                qgmu*aov(nu)+aog1(mu,space)*qnu &
                +qmu*aog1(nu,space)+aov(mu)*qgnu)
            end do
          end do
        end do
      end do
    end do
  end subroutine gga_density_nuclear_point_first

end module mod_dft_gridint_tdgga_response
