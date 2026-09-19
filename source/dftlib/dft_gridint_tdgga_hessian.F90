! Pointwise restricted-GGA fixed-density TDDFT Hessian assembly.
!
! This module consumes the AO-centre derivatives from
! mod_dft_gga_nuclear_point and the normalized partition-weight derivatives
! from mod_dft_partition_hessian.  It implements the TDHXFAB/TDHXD2G chain
! rule without finite differences of a molecular XC gradient.
module mod_dft_gridint_tdgga_hessian

  use precision, only: fp
  use mod_dft_gridint_tdxc_deriv, only: gga_tdxc_kernel_t, &
    gga_tdxc_variables_t, gga_tdxc_fixed_value, &
    gga_tdxc_fixed_derivatives
  use mod_dft_gga_nuclear_point, only: gga_density_nuclear_point
  use mod_dft_partition_hessian, only: partition_weight_nuclear_derivatives

  implicit none
  private

  public :: gga_add_owner_motion
  public :: gga_build_variable_derivatives
  public :: gga_weighted_point_hessian
  public :: gga_fixed_density_quadrature_point

contains

!> Complete contribution of one owner-centred quadrature point.
!>
!> density_r, density_p, and density_v generate respectively the ground,
!> effective-difference, and singlet transition density fields.  AO values
!> and electronic-coordinate derivatives must be evaluated at point.  The
!> routine adds the resulting contribution to hessian and returns a nonzero
!> status only when the partition-weight construction fails.
  subroutine gga_fixed_density_quadrature_point(kernel, density_r, density_p, &
      density_v, ao_atom, aov, aog1, aog2, aog3, atom_xyz, point, owner, &
      dummy_atom, part_fun_type, has_surface_shift, surface_shift, &
      finite_weight, hessian, status)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    real(fp), intent(in) :: density_r(:,:), density_p(:,:), density_v(:,:)
    integer, intent(in) :: ao_atom(:), owner, part_fun_type
    real(fp), intent(in) :: aov(:), aog1(:,:), aog2(:,:), aog3(:,:)
    real(fp), intent(in) :: atom_xyz(:,:), point(3), surface_shift(:,:)
    logical, intent(in) :: dummy_atom(:), has_surface_shift
    real(fp), intent(in) :: finite_weight
    real(fp), intent(inout) :: hessian(:,:,:,:)
    integer, intent(out) :: status

    integer :: nat
    real(fp), parameter :: rho_cutoff=1.0e-8_fp
    real(fp) :: r, p, v, grad_r(3), grad_p(3), grad_v(3)
    real(fp), allocatable :: fr1(:,:), fp1(:,:), fv1(:,:)
    real(fp), allocatable :: fgr1(:,:,:), fgp1(:,:,:), fgv1(:,:,:)
    real(fp), allocatable :: fr2(:,:,:,:), fp2(:,:,:,:), fv2(:,:,:,:)
    real(fp), allocatable :: fgr2(:,:,:,:,:), fgp2(:,:,:,:,:)
    real(fp), allocatable :: fgv2(:,:,:,:,:)
    real(fp), allocatable :: dr(:,:), dp(:,:), dv(:,:)
    real(fp), allocatable :: dgr(:,:,:), dgp(:,:,:), dgv(:,:,:)
    real(fp), allocatable :: d2r(:,:,:,:), d2p(:,:,:,:), d2v(:,:,:,:)
    real(fp), allocatable :: d2gr(:,:,:,:,:), d2gp(:,:,:,:,:)
    real(fp), allocatable :: d2gv(:,:,:,:,:)
    real(fp), allocatable :: weights(:), dweights(:,:,:), d2weights(:,:,:,:,:)
    real(fp), allocatable :: du(:,:,:), d2u(:,:,:,:,:)
    type(gga_tdxc_variables_t) :: u

    nat = size(atom_xyz,2)
    status = 0
    if (size(atom_xyz,1) /= 3 .or. size(dummy_atom) /= nat .or. &
        any(shape(surface_shift) /= [nat,nat]) .or. &
        any(shape(hessian) /= [3,3,nat,nat])) then
      status = 10
      return
    end if
    allocate(fr1(3,nat),fp1(3,nat),fv1(3,nat), &
      fgr1(3,3,nat),fgp1(3,3,nat),fgv1(3,3,nat), &
      fr2(3,3,nat,nat),fp2(3,3,nat,nat),fv2(3,3,nat,nat), &
      fgr2(3,3,3,nat,nat),fgp2(3,3,3,nat,nat), &
      fgv2(3,3,3,nat,nat),dr(3,nat),dp(3,nat),dv(3,nat), &
      dgr(3,3,nat),dgp(3,3,nat),dgv(3,3,nat), &
      d2r(3,3,nat,nat),d2p(3,3,nat,nat),d2v(3,3,nat,nat), &
      d2gr(3,3,3,nat,nat),d2gp(3,3,3,nat,nat), &
      d2gv(3,3,3,nat,nat),weights(nat),dweights(3,nat,nat), &
      d2weights(3,nat,3,nat,nat),du(7,3,nat),d2u(7,3,3,nat,nat))

    call field_value_gradient(density_r,aov,aog1,r,grad_r)
    ! Match GAMESS TDHXD2G/TDHXDPG: the complete point contribution is
    ! omitted below RCUTG, not only the finite-difference fourth kernel.
    if (r < rho_cutoff) return
    call field_value_gradient(density_p,aov,aog1,p,grad_p)
    call field_value_gradient(density_v,aov,aog1,v,grad_v)
    call gga_density_nuclear_point(density_r,ao_atom,aov,aog1,aog2,aog3, &
      fr1,fgr1,fr2,fgr2)
    call gga_density_nuclear_point(density_p,ao_atom,aov,aog1,aog2,aog3, &
      fp1,fgp1,fp2,fgp2)
    call gga_density_nuclear_point(density_v,ao_atom,aov,aog1,aog2,aog3, &
      fv1,fgv1,fv2,fgv2)
    call gga_add_owner_motion(owner,fr1,fgr1,fr2,fgr2,dr,dgr,d2r,d2gr)
    call gga_add_owner_motion(owner,fp1,fgp1,fp2,fgp2,dp,dgp,d2p,d2gp)
    call gga_add_owner_motion(owner,fv1,fgv1,fv2,fgv2,dv,dgv,d2v,d2gv)
    call gga_build_variable_derivatives(r,p,v,grad_r,grad_p,grad_v, &
      dr,dp,dv,dgr,dgp,dgv,d2r,d2p,d2v,d2gr,d2gp,d2gv,u,du,d2u)
    call partition_weight_nuclear_derivatives(atom_xyz,point,owner,dummy_atom, &
      part_fun_type,has_surface_shift,surface_shift,weights,dweights, &
      d2weights,status)
    if (status /= 0) return
    if (weights(owner) <= sqrt(tiny(1.0_fp))) return
    call gga_weighted_point_hessian(kernel,u,du,d2u, &
      finite_weight/weights(owner),weights(owner),dweights(:,:,owner), &
      d2weights(:,:,:,:,owner),hessian)

  contains

    subroutine field_value_gradient(density,values,gradients,value,gradient)
      real(fp), intent(in) :: density(:,:), values(:), gradients(:,:)
      real(fp), intent(out) :: value, gradient(3)
      integer :: mu, nu
      if (size(density,1) /= size(values) .or. &
          size(density,2) /= size(values) .or. &
          any(shape(gradients) /= [size(values),3])) &
        error stop 'gga_fixed_density_quadrature_point: AO shape mismatch'
      value = 0.0_fp
      gradient = 0.0_fp
      do nu = 1, size(values)
        do mu = 1, size(values)
          value = value+density(mu,nu)*values(mu)*values(nu)
          gradient = gradient+density(mu,nu)*( &
            gradients(mu,:)*values(nu)+values(mu)*gradients(nu,:))
        end do
      end do
    end subroutine field_value_gradient

  end subroutine gga_fixed_density_quadrature_point

!> Add the rigid translation of an atom-centred integration point.
!>
!> The input derivatives are at fixed laboratory grid position.  Spatial
!> derivatives not present explicitly in the AO helper follow from
!> translational invariance: grad(f)=-sum_A d_A(f),
!> grad[d_B(f)]=-sum_A d_A d_B(f), and
!> grad grad(f)=sum_AB d_A d_B(f).  Consequently AO G3 is sufficient for the
!> second derivative of grad(f), including owner-grid motion.
  subroutine gga_add_owner_motion(owner, fixed_d1, fixed_grad_d1, &
      fixed_d2, fixed_grad_d2, d1, grad_d1, d2, grad_d2)
    integer, intent(in) :: owner
    real(fp), intent(in) :: fixed_d1(:,:), fixed_grad_d1(:,:,:)
    real(fp), intent(in) :: fixed_d2(:,:,:,:), fixed_grad_d2(:,:,:,:,:)
    real(fp), intent(out) :: d1(:,:), grad_d1(:,:,:)
    real(fp), intent(out) :: d2(:,:,:,:), grad_d2(:,:,:,:,:)

    integer :: nat, a, b, atom_a, atom_b, c, atom_c
    real(fp) :: hess(3,3), hess_from_d2(3,3), third(3,3,3)
    real(fp) :: spatial_fixed_a, spatial_fixed_b

    nat = size(fixed_d1,2)
    call check_field_shapes(nat)
    if (owner < 1 .or. owner > nat) &
      error stop 'gga_add_owner_motion: owner out of range'

    hess = 0.0_fp
    hess_from_d2 = 0.0_fp
    third = 0.0_fp
    do atom_a = 1, nat
      hess = hess-fixed_grad_d1(:,:,atom_a)
      do atom_b = 1, nat
        hess_from_d2 = hess_from_d2+fixed_d2(:,:,atom_a,atom_b)
        third = third+fixed_grad_d2(:,:,:,atom_a,atom_b)
      end do
    end do

    d1 = fixed_d1
    grad_d1 = fixed_grad_d1
    d2 = fixed_d2
    grad_d2 = fixed_grad_d2
    do atom_a = 1, nat
      do a = 1, 3
        if (atom_a == owner) then
          d1(a,atom_a) = d1(a,atom_a)-sum(fixed_d1(a,:))
          grad_d1(:,a,atom_a) = grad_d1(:,a,atom_a)+hess(:,a)
        end if
        do atom_b = 1, nat
          do b = 1, 3
            spatial_fixed_a = 0.0_fp
            spatial_fixed_b = 0.0_fp
            do atom_c = 1, nat
              spatial_fixed_a = spatial_fixed_a &
                -fixed_d2(b,a,atom_b,atom_c)
              spatial_fixed_b = spatial_fixed_b &
                -fixed_d2(a,b,atom_a,atom_c)
            end do
            if (atom_a == owner) &
              d2(a,b,atom_a,atom_b) = d2(a,b,atom_a,atom_b)+spatial_fixed_a
            if (atom_b == owner) &
              d2(a,b,atom_a,atom_b) = d2(a,b,atom_a,atom_b)+spatial_fixed_b
            if (atom_a == owner .and. atom_b == owner) &
              d2(a,b,atom_a,atom_b) = &
                d2(a,b,atom_a,atom_b)+hess_from_d2(a,b)

            do c = 1, 3
              if (atom_a == owner) grad_d2(c,a,b,atom_a,atom_b) = &
                grad_d2(c,a,b,atom_a,atom_b) &
                -sum(fixed_grad_d2(c,b,a,atom_b,:),dim=1)
              if (atom_b == owner) grad_d2(c,a,b,atom_a,atom_b) = &
                grad_d2(c,a,b,atom_a,atom_b) &
                -sum(fixed_grad_d2(c,a,b,atom_a,:),dim=1)
              if (atom_a == owner .and. atom_b == owner) &
                grad_d2(c,a,b,atom_a,atom_b) = &
                grad_d2(c,a,b,atom_a,atom_b)+third(c,a,b)
            end do
          end do
        end do
      end do
    end do

  contains

    subroutine check_field_shapes(n)
      integer, intent(in) :: n
      if (any(shape(fixed_d1) /= [3,n]) .or. &
          any(shape(fixed_grad_d1) /= [3,3,n]) .or. &
          any(shape(fixed_d2) /= [3,3,n,n]) .or. &
          any(shape(fixed_grad_d2) /= [3,3,3,n,n]) .or. &
          any(shape(d1) /= [3,n]) .or. &
          any(shape(grad_d1) /= [3,3,n]) .or. &
          any(shape(d2) /= [3,3,n,n]) .or. &
          any(shape(grad_d2) /= [3,3,3,n,n])) &
        error stop 'gga_add_owner_motion: shape mismatch'
    end subroutine check_field_shapes

  end subroutine gga_add_owner_motion

!> Form derivatives of u=(rho,sigma,rho_P,Gamma_P,rho_v,Gamma_v,W_v).
!>
!> Gamma_P=grad(rho).grad(rho_P), Gamma_v=grad(rho).grad(rho_v),
!> sigma=|grad(rho)|2, and W_v=|grad(rho_v)|2.
  subroutine gga_build_variable_derivatives(r, p, v, grad_r, grad_p, grad_v, &
      dr, dp, dv, dgrad_r, dgrad_p, dgrad_v, &
      d2r, d2p, d2v, d2grad_r, d2grad_p, d2grad_v, u, du, d2u)
    real(fp), intent(in) :: r, p, v, grad_r(3), grad_p(3), grad_v(3)
    real(fp), intent(in) :: dr(:,:), dp(:,:), dv(:,:)
    real(fp), intent(in) :: dgrad_r(:,:,:), dgrad_p(:,:,:), dgrad_v(:,:,:)
    real(fp), intent(in) :: d2r(:,:,:,:), d2p(:,:,:,:), d2v(:,:,:,:)
    real(fp), intent(in) :: d2grad_r(:,:,:,:,:), d2grad_p(:,:,:,:,:)
    real(fp), intent(in) :: d2grad_v(:,:,:,:,:)
    type(gga_tdxc_variables_t), intent(out) :: u
    real(fp), intent(out) :: du(:,:,:), d2u(:,:,:,:,:)

    integer :: nat, a, b, atom_a, atom_b

    nat = size(dr,2)
    call check_shapes(nat)
    u%r = r
    u%s = dot_product(grad_r,grad_r)
    u%p = p
    u%pg = dot_product(grad_r,grad_p)
    u%v = v
    u%gv = dot_product(grad_r,grad_v)
    u%w = dot_product(grad_v,grad_v)
    du = 0.0_fp
    d2u = 0.0_fp

    do atom_a = 1, nat
      do a = 1, 3
        du(1,a,atom_a) = dr(a,atom_a)
        du(2,a,atom_a) = 2.0_fp*dot_product(grad_r,dgrad_r(:,a,atom_a))
        du(3,a,atom_a) = dp(a,atom_a)
        du(4,a,atom_a) = dot_product(dgrad_r(:,a,atom_a),grad_p) &
          +dot_product(grad_r,dgrad_p(:,a,atom_a))
        du(5,a,atom_a) = dv(a,atom_a)
        du(6,a,atom_a) = dot_product(dgrad_r(:,a,atom_a),grad_v) &
          +dot_product(grad_r,dgrad_v(:,a,atom_a))
        du(7,a,atom_a) = 2.0_fp*dot_product(grad_v,dgrad_v(:,a,atom_a))
        do atom_b = 1, nat
          do b = 1, 3
            d2u(1,a,b,atom_a,atom_b) = d2r(a,b,atom_a,atom_b)
            d2u(2,a,b,atom_a,atom_b) = 2.0_fp*( &
              dot_product(dgrad_r(:,a,atom_a),dgrad_r(:,b,atom_b)) &
              +dot_product(grad_r,d2grad_r(:,a,b,atom_a,atom_b)))
            d2u(3,a,b,atom_a,atom_b) = d2p(a,b,atom_a,atom_b)
            d2u(4,a,b,atom_a,atom_b) = &
              dot_product(d2grad_r(:,a,b,atom_a,atom_b),grad_p) &
              +dot_product(dgrad_r(:,a,atom_a),dgrad_p(:,b,atom_b)) &
              +dot_product(dgrad_r(:,b,atom_b),dgrad_p(:,a,atom_a)) &
              +dot_product(grad_r,d2grad_p(:,a,b,atom_a,atom_b))
            d2u(5,a,b,atom_a,atom_b) = d2v(a,b,atom_a,atom_b)
            d2u(6,a,b,atom_a,atom_b) = &
              dot_product(d2grad_r(:,a,b,atom_a,atom_b),grad_v) &
              +dot_product(dgrad_r(:,a,atom_a),dgrad_v(:,b,atom_b)) &
              +dot_product(dgrad_r(:,b,atom_b),dgrad_v(:,a,atom_a)) &
              +dot_product(grad_r,d2grad_v(:,a,b,atom_a,atom_b))
            d2u(7,a,b,atom_a,atom_b) = 2.0_fp*( &
              dot_product(dgrad_v(:,a,atom_a),dgrad_v(:,b,atom_b)) &
              +dot_product(grad_v,d2grad_v(:,a,b,atom_a,atom_b)))
          end do
        end do
      end do
    end do

  contains

    subroutine check_shapes(n)
      integer, intent(in) :: n
      if (any(shape(dp) /= [3,n]) .or. any(shape(dv) /= [3,n]) .or. &
          any(shape(dgrad_r) /= [3,3,n]) .or. &
          any(shape(dgrad_p) /= [3,3,n]) .or. &
          any(shape(dgrad_v) /= [3,3,n]) .or. &
          any(shape(d2r) /= [3,3,n,n]) .or. &
          any(shape(d2p) /= [3,3,n,n]) .or. &
          any(shape(d2v) /= [3,3,n,n]) .or. &
          any(shape(d2grad_r) /= [3,3,3,n,n]) .or. &
          any(shape(d2grad_p) /= [3,3,3,n,n]) .or. &
          any(shape(d2grad_v) /= [3,3,3,n,n]) .or. &
          any(shape(du) /= [7,3,n]) .or. &
          any(shape(d2u) /= [7,3,3,n,n])) &
        error stop 'gga_build_variable_derivatives: shape mismatch'
    end subroutine check_shapes

  end subroutine gga_build_variable_derivatives

!> Add one quadrature point to the fixed-density GGA Hessian.
!>
!> quadrature_scale contains the radial/angular factor.  weight and its
!> derivatives are the selected normalized atom-cell partition weight.
  subroutine gga_weighted_point_hessian(kernel, u, du, d2u, &
      quadrature_scale, weight, dweight, d2weight, hessian, point_value)
    type(gga_tdxc_kernel_t), intent(in) :: kernel
    type(gga_tdxc_variables_t), intent(in) :: u
    real(fp), intent(in) :: du(:,:,:), d2u(:,:,:,:,:)
    real(fp), intent(in) :: quadrature_scale, weight
    real(fp), intent(in) :: dweight(:,:), d2weight(:,:,:,:)
    real(fp), intent(inout) :: hessian(:,:,:,:)
    real(fp), intent(out), optional :: point_value

    real(fp) :: fa(7), fab(7,7), q, qa, qb, qab
    integer :: nat, a, b, atom_a, atom_b, i, j

    nat = size(dweight,2)
    if (any(shape(du) /= [7,3,nat]) .or. &
        any(shape(d2u) /= [7,3,3,nat,nat]) .or. &
        any(shape(d2weight) /= [3,nat,3,nat]) .or. &
        any(shape(hessian) /= [3,3,nat,nat])) &
      error stop 'gga_weighted_point_hessian: shape mismatch'

    q = gga_tdxc_fixed_value(kernel,u)
    call gga_tdxc_fixed_derivatives(kernel,u,fa,fab)
    if (present(point_value)) point_value = quadrature_scale*weight*q
    do atom_a = 1, nat
      do a = 1, 3
        qa = dot_product(fa,du(:,a,atom_a))
        do atom_b = 1, nat
          do b = 1, 3
            qb = dot_product(fa,du(:,b,atom_b))
            qab = dot_product(fa,d2u(:,a,b,atom_a,atom_b))
            do i = 1, 7
              do j = 1, 7
                qab = qab+fab(i,j)*du(i,a,atom_a)*du(j,b,atom_b)
              end do
            end do
            hessian(a,b,atom_a,atom_b) = hessian(a,b,atom_a,atom_b) &
              +quadrature_scale*(weight*qab+dweight(a,atom_a)*qb &
              +dweight(b,atom_b)*qa+d2weight(a,atom_a,b,atom_b)*q)
          end do
        end do
      end do
    end do
  end subroutine gga_weighted_point_hessian

end module mod_dft_gridint_tdgga_hessian
