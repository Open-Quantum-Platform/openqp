!> @brief Native one-dimensional-projection (ODP) umbrella bias.
!>
!> The raw collective variables q are mapped into an explicitly scaled,
!> dimensionless space X_i = scale_i*q_i before the reactant/product
!> projection is formed.  The C ABI returns the conservative atomic force,
!> not the gradient, so callers can add it directly to an MD force.
module odp_umbrella_mod

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private

  integer(c_int), parameter, public :: ODP_CV_DISTANCE = 1_c_int
  integer(c_int), parameter, public :: ODP_CV_ASYMMETRIC_DISTANCE = 2_c_int
  integer(c_int), parameter, public :: ODP_CV_ANGLE = 3_c_int

  real(c_double), parameter :: distance_epsilon = 1.0e-12_c_double
  real(c_double), parameter :: angle_sine_epsilon = 1.0e-10_c_double
  real(c_double), parameter :: path_epsilon_squared = 1.0e-24_c_double

  public :: odp_umbrella_evaluate

contains

!> @brief Value and Cartesian Jacobian of one continuous internal CV.
  subroutine odp_cv_value_jacobian(natom, cv_type, atoms, coordinates, &
      value, jacobian, info)
    integer, intent(in) :: natom
    integer(c_int), intent(in) :: cv_type
    integer(c_int64_t), intent(in) :: atoms(4)
    real(c_double), intent(in) :: coordinates(3, natom)
    real(c_double), intent(out) :: value
    real(c_double), intent(out) :: jacobian(3, natom)
    integer(c_int), intent(out) :: info
    integer :: a, b, c, d
    real(c_double) :: vab(3), vcd(3), unit_ab(3), unit_cd(3)
    real(c_double) :: rab, rcd, cosine, sine, grad_a(3), grad_c(3)

    info = 0_c_int
    value = 0.0_c_double
    jacobian = 0.0_c_double

    select case (cv_type)
    case (ODP_CV_DISTANCE)
      a = int(atoms(1)) + 1
      b = int(atoms(2)) + 1
      if (a < 1 .or. a > natom .or. b < 1 .or. b > natom .or. a == b) then
        info = 2_c_int
        return
      end if
      vab = coordinates(:, a) - coordinates(:, b)
      rab = sqrt(dot_product(vab, vab))
      if (.not. ieee_is_finite(rab) .or. rab <= distance_epsilon) then
        info = 4_c_int
        return
      end if
      unit_ab = vab/rab
      value = rab
      jacobian(:, a) = unit_ab
      jacobian(:, b) = -unit_ab

    case (ODP_CV_ASYMMETRIC_DISTANCE)
      a = int(atoms(1)) + 1
      b = int(atoms(2)) + 1
      c = int(atoms(3)) + 1
      d = int(atoms(4)) + 1
      if (min(a, b, c, d) < 1 .or. max(a, b, c, d) > natom .or. &
          a == b .or. c == d) then
        info = 2_c_int
        return
      end if
      vab = coordinates(:, a) - coordinates(:, b)
      vcd = coordinates(:, c) - coordinates(:, d)
      rab = sqrt(dot_product(vab, vab))
      rcd = sqrt(dot_product(vcd, vcd))
      if (.not. ieee_is_finite(rab) .or. .not. ieee_is_finite(rcd) .or. &
          rab <= distance_epsilon .or. rcd <= distance_epsilon) then
        info = 4_c_int
        return
      end if
      unit_ab = vab/rab
      unit_cd = vcd/rcd
      value = rab - rcd
      jacobian(:, a) = jacobian(:, a) + unit_ab
      jacobian(:, b) = jacobian(:, b) - unit_ab
      jacobian(:, c) = jacobian(:, c) - unit_cd
      jacobian(:, d) = jacobian(:, d) + unit_cd

    case (ODP_CV_ANGLE)
      a = int(atoms(1)) + 1
      b = int(atoms(2)) + 1
      c = int(atoms(3)) + 1
      if (min(a, b, c) < 1 .or. max(a, b, c) > natom .or. &
          a == b .or. b == c .or. a == c) then
        info = 2_c_int
        return
      end if
      vab = coordinates(:, a) - coordinates(:, b)
      vcd = coordinates(:, c) - coordinates(:, b)
      rab = sqrt(dot_product(vab, vab))
      rcd = sqrt(dot_product(vcd, vcd))
      if (.not. ieee_is_finite(rab) .or. .not. ieee_is_finite(rcd) .or. &
          rab <= distance_epsilon .or. rcd <= distance_epsilon) then
        info = 5_c_int
        return
      end if
      cosine = dot_product(vab, vcd)/(rab*rcd)
      cosine = max(-1.0_c_double, min(1.0_c_double, cosine))
      sine = sqrt(max(0.0_c_double, 1.0_c_double - cosine*cosine))
      if (.not. ieee_is_finite(sine) .or. sine <= angle_sine_epsilon) then
        info = 6_c_int
        return
      end if
      value = acos(cosine)
      grad_a = -(vcd/(rab*rcd) - cosine*vab/(rab*rab))/sine
      grad_c = -(vab/(rab*rcd) - cosine*vcd/(rcd*rcd))/sine
      jacobian(:, a) = grad_a
      jacobian(:, c) = grad_c
      jacobian(:, b) = -grad_a - grad_c

    case default
      info = 2_c_int
    end select
  end subroutine odp_cv_value_jacobian

!> @brief Evaluate a scaled ODP umbrella energy and atomic force.
!>
!> cv_atoms uses zero-based C/Python atom indices.  Each row contains two
!> indices for distance, four for asymmetric distance d(a,b)-d(c,d), or three
!> for angle a-b-c.  Angles and their references are in radians.
!>
!> Status: 0 success; 1 invalid dimensions/force constants; 2 invalid CV;
!> 3 invalid scale/reference; 4 singular distance; 5 zero angle bond;
!> 6 collinear angle; 7 degenerate R/P path; 8 non-finite input/output.
  function odp_umbrella_evaluate(natom_c, ncv_c, cv_types, cv_atoms, &
      cv_scales, reference_r, reference_p, center, k_parallel, &
      k_perpendicular, coordinates, energy, force, xi, cv_raw, cv_scaled, &
      cv_perpendicular, perpendicular_norm, energy_parallel, &
      energy_perpendicular) result(info) &
      bind(C, name="oqp_odp_umbrella_evaluate")
    integer(c_int64_t), value, intent(in) :: natom_c, ncv_c
    integer(c_int), intent(in) :: cv_types(*)
    integer(c_int64_t), intent(in) :: cv_atoms(4, *)
    real(c_double), intent(in) :: cv_scales(*), reference_r(*), reference_p(*)
    real(c_double), value, intent(in) :: center, k_parallel, k_perpendicular
    real(c_double), intent(in) :: coordinates(3, *)
    real(c_double), intent(out) :: energy, force(3, *), xi
    real(c_double), intent(out) :: cv_raw(*), cv_scaled(*)
    real(c_double), intent(out) :: cv_perpendicular(*)
    real(c_double), intent(out) :: perpendicular_norm
    real(c_double), intent(out) :: energy_parallel, energy_perpendicular
    integer(c_int) :: info
    integer :: natom, ncv, atom_index, cv_index
    real(c_double), allocatable :: jacobian(:, :), direction(:), &
                                   displacement(:), dudx(:)
    real(c_double) :: path_length_squared, projection_length
    real(c_double) :: cv_value, coefficient
    integer(c_int) :: cv_info

    info = 1_c_int
    energy = 0.0_c_double
    xi = 0.0_c_double
    perpendicular_norm = 0.0_c_double
    energy_parallel = 0.0_c_double
    energy_perpendicular = 0.0_c_double
    if (natom_c <= 0_c_int64_t .or. ncv_c <= 0_c_int64_t) return
    if (.not. ieee_is_finite(center) .or. &
        .not. ieee_is_finite(k_parallel) .or. &
        .not. ieee_is_finite(k_perpendicular) .or. &
        k_parallel < 0.0_c_double .or. k_perpendicular < 0.0_c_double) return

    natom = int(natom_c)
    ncv = int(ncv_c)
    force(:, 1:natom) = 0.0_c_double
    cv_raw(1:ncv) = 0.0_c_double
    cv_scaled(1:ncv) = 0.0_c_double
    cv_perpendicular(1:ncv) = 0.0_c_double
    do atom_index = 1, natom
      if (any(.not. ieee_is_finite(coordinates(:, atom_index)))) then
        info = 8_c_int
        return
      end if
    end do

    allocate(jacobian(3, natom), direction(ncv), displacement(ncv), dudx(ncv))
    direction = 0.0_c_double
    displacement = 0.0_c_double
    dudx = 0.0_c_double

    do cv_index = 1, ncv
      if (.not. ieee_is_finite(cv_scales(cv_index)) .or. &
          cv_scales(cv_index) <= 0.0_c_double .or. &
          .not. ieee_is_finite(reference_r(cv_index)) .or. &
          .not. ieee_is_finite(reference_p(cv_index))) then
        info = 3_c_int
        return
      end if
      call odp_cv_value_jacobian(natom, cv_types(cv_index), &
          cv_atoms(:, cv_index), coordinates(:, 1:natom), cv_value, &
          jacobian, cv_info)
      if (cv_info /= 0_c_int) then
        info = cv_info
        return
      end if
      cv_raw(cv_index) = cv_value
      cv_scaled(cv_index) = cv_scales(cv_index)*cv_value
      direction(cv_index) = cv_scales(cv_index)*( &
          reference_p(cv_index) - reference_r(cv_index))
      displacement(cv_index) = cv_scaled(cv_index) - &
          cv_scales(cv_index)*reference_r(cv_index)
    end do

    path_length_squared = dot_product(direction, direction)
    if (.not. ieee_is_finite(path_length_squared) .or. &
        path_length_squared <= path_epsilon_squared) then
      info = 7_c_int
      return
    end if
    projection_length = dot_product(displacement, direction)
    xi = projection_length/path_length_squared
    cv_perpendicular(1:ncv) = displacement - xi*direction
    perpendicular_norm = sqrt(dot_product( &
        cv_perpendicular(1:ncv), cv_perpendicular(1:ncv)))
    energy_parallel = 0.5_c_double*k_parallel*(xi - center)**2
    energy_perpendicular = 0.5_c_double*k_perpendicular* &
        perpendicular_norm**2
    energy = energy_parallel + energy_perpendicular

    dudx = k_parallel*(xi - center)*direction/path_length_squared + &
           k_perpendicular*cv_perpendicular(1:ncv)
    do cv_index = 1, ncv
      call odp_cv_value_jacobian(natom, cv_types(cv_index), &
          cv_atoms(:, cv_index), coordinates(:, 1:natom), cv_value, &
          jacobian, cv_info)
      if (cv_info /= 0_c_int) then
        info = cv_info
        return
      end if
      coefficient = cv_scales(cv_index)*dudx(cv_index)
      force(:, 1:natom) = force(:, 1:natom) - coefficient*jacobian
    end do

    if (.not. ieee_is_finite(energy) .or. .not. ieee_is_finite(xi) .or. &
        .not. ieee_is_finite(perpendicular_norm) .or. &
        any(.not. ieee_is_finite(force(:, 1:natom)))) then
      info = 8_c_int
      return
    end if
    info = 0_c_int
  end function odp_umbrella_evaluate

end module odp_umbrella_mod
