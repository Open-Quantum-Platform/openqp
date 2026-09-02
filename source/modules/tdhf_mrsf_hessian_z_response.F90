module tdhf_mrsf_hessian_z_response_mod

  use precision, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use tdhf_sf_lib, only: sfrorhs, sfrolhs
  use tdhf_hessian_response_mod, only: solve_mrsf_z_response_matrix_free, &
    solve_mrsf_z_response_batch_matrix_free

  implicit none

  private

  integer, parameter, public :: MRSF_Z_RESPONSE_SUCCESS = 0
  integer, parameter, public :: MRSF_Z_RESPONSE_INVALID_SHAPE = -1
  integer, parameter, public :: MRSF_Z_RESPONSE_NOT_TWO_SOMO = -2
  integer, parameter, public :: MRSF_Z_RESPONSE_UNSUPPORTED_MULTIPLICITY = -3
  integer, parameter, public :: MRSF_Z_RESPONSE_UMRSF_UNSUPPORTED = -4
  integer, parameter, public :: MRSF_Z_RESPONSE_CAM_UNSUPPORTED = -5
  integer, parameter, public :: MRSF_Z_RESPONSE_NONFINITE = -6

  public :: build_mrsf_z_rhs_derivative
  public :: build_mrsf_z_operator_derivative_action
  public :: solve_mrsf_z_response_from_mo_derivatives

  ! Hiroya Nakata's TDHF/TDDFT Hessian formulation is the methodological
  ! starting point for this differentiated MRSF Z-vector equation.
  !
  ! The supplied quantities must already represent the complete spin-adapted
  ! CO(+), OV(+), CV(+), OO, CV(-), OV(-), CO(-) response topology.  This
  ! routine is restricted to ordinary two-SOMO singlet and triplet targets.

  abstract interface
    subroutine mrsf_z_orbital_hessian_action(vector, result, status)
      import dp
      real(kind=dp), intent(in) :: vector(:)
      real(kind=dp), intent(out) :: result(:)
      integer, intent(out) :: status
    end subroutine mrsf_z_orbital_hessian_action

    subroutine mrsf_z_orbital_hessian_batch_action(vectors, results, status)
      import dp
      real(kind=dp), intent(in) :: vectors(:,:)
      real(kind=dp), intent(out) :: results(:,:)
      integer, intent(out) :: status
    end subroutine mrsf_z_orbital_hessian_batch_action
  end interface

contains

!###############################################################################

  subroutine validate_mrsf_z_method(target_multiplicity, is_umrsf, cam_flag, &
                                    nocca, noccb, status)
    integer, intent(in) :: target_multiplicity, nocca, noccb
    logical, intent(in) :: is_umrsf, cam_flag
    integer, intent(out) :: status

    status = MRSF_Z_RESPONSE_SUCCESS
    if (is_umrsf) then
      status = MRSF_Z_RESPONSE_UMRSF_UNSUPPORTED
      return
    end if
    if (cam_flag) then
      status = MRSF_Z_RESPONSE_CAM_UNSUPPORTED
      return
    end if
    if (target_multiplicity /= 1 .and. target_multiplicity /= 3) then
      status = MRSF_Z_RESPONSE_UNSUPPORTED_MULTIPLICITY
      return
    end if
    if (nocca <= 0 .or. noccb < 0 .or. noccb >= nocca) then
      status = MRSF_Z_RESPONSE_INVALID_SHAPE
      return
    end if
    if (nocca - noccb /= 2) then
      status = MRSF_Z_RESPONSE_NOT_TWO_SOMO
      return
    end if
  end subroutine validate_mrsf_z_method

!###############################################################################

  subroutine evaluate_mrsf_z_rhs(hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, &
                                 tij, tab, fa, fb, nocca, noccb, rhs)
    real(kind=dp), intent(in) :: hxa(:,:), hxb(:,:)
    real(kind=dp), intent(in) :: rhs_ab1_mo_a(:,:), rhs_ab1_mo_b(:,:)
    real(kind=dp), intent(in) :: tij(:,:), tab(:,:), fa(:,:), fb(:,:)
    integer, intent(in) :: nocca, noccb
    real(kind=dp), intent(out) :: rhs(:)

    real(kind=dp), allocatable :: hxa_work(:,:), hxb_work(:,:)

    allocate(hxa_work, source=hxa)
    allocate(hxb_work, source=hxb)
    call sfrorhs(rhs, hxa_work, hxb_work, rhs_ab1_mo_a, rhs_ab1_mo_b, &
                 tij, tab, fa, fb, nocca, noccb)
    deallocate(hxa_work, hxb_work)
  end subroutine evaluate_mrsf_z_rhs

!###############################################################################

  subroutine build_mrsf_z_rhs_derivative(target_multiplicity, is_umrsf, &
      cam_flag, nocca, noccb, fa, fb, hxa, hxb, rhs_ab1_mo_a, &
      rhs_ab1_mo_b, tij, tab, dfa, dfb, dhxa, dhxb, drhs_ab1_mo_a, &
      drhs_ab1_mo_b, dtij, dtab, drhs, status)
    ! Differentiate the spin-adapted MRSF Z-vector right-hand side by exact
    ! algebraic polarization of the variables entering sfrorhs.
    !
    ! The direct response, Fock, and difference-density directions are
    ! polarized separately.  Thus products such as (dF)(dT), which are second
    ! order in a nuclear displacement, never enter the first derivative.

    integer, intent(in) :: target_multiplicity, nocca, noccb
    logical, intent(in) :: is_umrsf, cam_flag
    real(kind=dp), intent(in) :: fa(:,:), fb(:,:), hxa(:,:), hxb(:,:)
    real(kind=dp), intent(in) :: rhs_ab1_mo_a(:,:), rhs_ab1_mo_b(:,:)
    real(kind=dp), intent(in) :: tij(:,:), tab(:,:)
    real(kind=dp), intent(in) :: dfa(:,:,:), dfb(:,:,:)
    real(kind=dp), intent(in) :: dhxa(:,:,:), dhxb(:,:,:)
    real(kind=dp), intent(in) :: drhs_ab1_mo_a(:,:,:)
    real(kind=dp), intent(in) :: drhs_ab1_mo_b(:,:,:)
    real(kind=dp), intent(in) :: dtij(:,:,:), dtab(:,:,:)
    real(kind=dp), intent(out) :: drhs(:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: plus(:), minus(:)
    integer :: coordinate, lzdim, nbf, ncoord, nvira, nvirb

    drhs = 0.0_dp
    call validate_mrsf_z_method(target_multiplicity, is_umrsf, cam_flag, &
                                nocca, noccb, status)
    if (status /= MRSF_Z_RESPONSE_SUCCESS) return

    nbf = size(fa, 1)
    ncoord = size(drhs, 2)
    nvira = nbf - nocca
    nvirb = nbf - noccb
    lzdim = noccb*(nocca - noccb + nvira) + (nocca - noccb)*nvira
    if (nbf <= nocca .or. ncoord <= 0 .or. &
        any(shape(fa) /= [nbf, nbf]) .or. &
        any(shape(fb) /= [nbf, nbf]) .or. &
        any(shape(hxa) /= [nbf, nocca]) .or. &
        any(shape(hxb) /= [nbf, nbf]) .or. &
        any(shape(rhs_ab1_mo_a) /= [nocca, nvira]) .or. &
        any(shape(rhs_ab1_mo_b) /= [noccb, nvirb]) .or. &
        any(shape(tij) /= [nocca, nocca]) .or. &
        any(shape(tab) /= [nvirb, nvirb]) .or. &
        any(shape(dfa) /= [nbf, nbf, ncoord]) .or. &
        any(shape(dfb) /= [nbf, nbf, ncoord]) .or. &
        any(shape(dhxa) /= [nbf, nocca, ncoord]) .or. &
        any(shape(dhxb) /= [nbf, nbf, ncoord]) .or. &
        any(shape(drhs_ab1_mo_a) /= [nocca, nvira, ncoord]) .or. &
        any(shape(drhs_ab1_mo_b) /= [noccb, nvirb, ncoord]) .or. &
        any(shape(dtij) /= [nocca, nocca, ncoord]) .or. &
        any(shape(dtab) /= [nvirb, nvirb, ncoord]) .or. &
        any(shape(drhs) /= [lzdim, ncoord])) then
      status = MRSF_Z_RESPONSE_INVALID_SHAPE
      return
    end if
    if (any(.not. ieee_is_finite(fa)) .or. &
        any(.not. ieee_is_finite(fb)) .or. &
        any(.not. ieee_is_finite(hxa)) .or. &
        any(.not. ieee_is_finite(hxb)) .or. &
        any(.not. ieee_is_finite(rhs_ab1_mo_a)) .or. &
        any(.not. ieee_is_finite(rhs_ab1_mo_b)) .or. &
        any(.not. ieee_is_finite(tij)) .or. &
        any(.not. ieee_is_finite(tab)) .or. &
        any(.not. ieee_is_finite(dfa)) .or. &
        any(.not. ieee_is_finite(dfb)) .or. &
        any(.not. ieee_is_finite(dhxa)) .or. &
        any(.not. ieee_is_finite(dhxb)) .or. &
        any(.not. ieee_is_finite(drhs_ab1_mo_a)) .or. &
        any(.not. ieee_is_finite(drhs_ab1_mo_b)) .or. &
        any(.not. ieee_is_finite(dtij)) .or. &
        any(.not. ieee_is_finite(dtab))) then
      status = MRSF_Z_RESPONSE_NONFINITE
      return
    end if

    allocate(plus(lzdim), minus(lzdim))
    do coordinate = 1, ncoord
      ! Direct response polarization: Hx and the RHS response image.
      call evaluate_mrsf_z_rhs(hxa + dhxa(:,:,coordinate), &
        hxb + dhxb(:,:,coordinate), &
        rhs_ab1_mo_a + drhs_ab1_mo_a(:,:,coordinate), &
        rhs_ab1_mo_b + drhs_ab1_mo_b(:,:,coordinate), &
        tij, tab, fa, fb, nocca, noccb, plus)
      call evaluate_mrsf_z_rhs(hxa - dhxa(:,:,coordinate), &
        hxb - dhxb(:,:,coordinate), &
        rhs_ab1_mo_a - drhs_ab1_mo_a(:,:,coordinate), &
        rhs_ab1_mo_b - drhs_ab1_mo_b(:,:,coordinate), &
        tij, tab, fa, fb, nocca, noccb, minus)
      drhs(:,coordinate) = 0.5_dp*(plus - minus)

      ! Fock polarization at fixed baseline difference densities.
      call evaluate_mrsf_z_rhs(hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, &
        tij, tab, fa + dfa(:,:,coordinate), fb + dfb(:,:,coordinate), &
        nocca, noccb, plus)
      call evaluate_mrsf_z_rhs(hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, &
        tij, tab, fa - dfa(:,:,coordinate), fb - dfb(:,:,coordinate), &
        nocca, noccb, minus)
      drhs(:,coordinate) = drhs(:,coordinate) + 0.5_dp*(plus - minus)

      ! Difference-density polarization at fixed baseline Fock matrices.
      call evaluate_mrsf_z_rhs(hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, &
        tij + dtij(:,:,coordinate), tab + dtab(:,:,coordinate), &
        fa, fb, nocca, noccb, plus)
      call evaluate_mrsf_z_rhs(hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, &
        tij - dtij(:,:,coordinate), tab - dtab(:,:,coordinate), &
        fa, fb, nocca, noccb, minus)
      drhs(:,coordinate) = drhs(:,coordinate) + 0.5_dp*(plus - minus)
    end do
    deallocate(plus, minus)

    if (any(.not. ieee_is_finite(drhs))) then
      drhs = 0.0_dp
      status = MRSF_Z_RESPONSE_NONFINITE
    end if
  end subroutine build_mrsf_z_rhs_derivative

!###############################################################################

  subroutine evaluate_mrsf_z_operator_action(z, mo_energy, fa, fb, &
      z_ab1_mo_a, z_ab1_mo_b, nocca, noccb, action)
    real(kind=dp), intent(in) :: z(:), mo_energy(:), fa(:,:), fb(:,:)
    real(kind=dp), intent(in) :: z_ab1_mo_a(:,:), z_ab1_mo_b(:,:)
    integer, intent(in) :: nocca, noccb
    real(kind=dp), intent(out) :: action(:)

    call sfrolhs(action, z, mo_energy, fa, fb, z_ab1_mo_a, z_ab1_mo_b, &
                 nocca, noccb)
  end subroutine evaluate_mrsf_z_operator_action

!###############################################################################

  subroutine build_mrsf_z_operator_derivative_action( &
      target_multiplicity, is_umrsf, cam_flag, nocca, noccb, z, &
      mo_energy, fa, fb, z_ab1_mo_a, z_ab1_mo_b, dmo_energy, dfa, &
      dfb, dz_ab1_mo_a, dz_ab1_mo_b, operator_derivative_z, status)
    ! Form (dH_orb)Z with the baseline Z held fixed.  The derivative response
    ! images dz_ab1_mo_a/b must likewise be evaluated at fixed Z.  The unknown
    ! dZ is obtained only after this quantity enters the response equation.

    integer, intent(in) :: target_multiplicity, nocca, noccb
    logical, intent(in) :: is_umrsf, cam_flag
    real(kind=dp), intent(in) :: z(:), mo_energy(:), fa(:,:), fb(:,:)
    real(kind=dp), intent(in) :: z_ab1_mo_a(:,:), z_ab1_mo_b(:,:)
    real(kind=dp), intent(in) :: dmo_energy(:,:), dfa(:,:,:), dfb(:,:,:)
    real(kind=dp), intent(in) :: dz_ab1_mo_a(:,:,:), dz_ab1_mo_b(:,:,:)
    real(kind=dp), intent(out) :: operator_derivative_z(:,:)
    integer, intent(out) :: status

    real(kind=dp), allocatable :: plus(:), minus(:)
    integer :: coordinate, lzdim, nbf, ncoord, nvira, nvirb

    operator_derivative_z = 0.0_dp
    call validate_mrsf_z_method(target_multiplicity, is_umrsf, cam_flag, &
                                nocca, noccb, status)
    if (status /= MRSF_Z_RESPONSE_SUCCESS) return

    nbf = size(fa, 1)
    ncoord = size(operator_derivative_z, 2)
    nvira = nbf - nocca
    nvirb = nbf - noccb
    lzdim = noccb*(nocca - noccb + nvira) + (nocca - noccb)*nvira
    if (nbf <= nocca .or. ncoord <= 0 .or. &
        size(z) /= lzdim .or. size(mo_energy) /= nbf .or. &
        any(shape(fa) /= [nbf, nbf]) .or. &
        any(shape(fb) /= [nbf, nbf]) .or. &
        any(shape(z_ab1_mo_a) /= [nocca, nvira]) .or. &
        any(shape(z_ab1_mo_b) /= [noccb, nvirb]) .or. &
        any(shape(dmo_energy) /= [nbf, ncoord]) .or. &
        any(shape(dfa) /= [nbf, nbf, ncoord]) .or. &
        any(shape(dfb) /= [nbf, nbf, ncoord]) .or. &
        any(shape(dz_ab1_mo_a) /= [nocca, nvira, ncoord]) .or. &
        any(shape(dz_ab1_mo_b) /= [noccb, nvirb, ncoord]) .or. &
        any(shape(operator_derivative_z) /= [lzdim, ncoord])) then
      status = MRSF_Z_RESPONSE_INVALID_SHAPE
      return
    end if
    if (any(.not. ieee_is_finite(z)) .or. &
        any(.not. ieee_is_finite(mo_energy)) .or. &
        any(.not. ieee_is_finite(fa)) .or. &
        any(.not. ieee_is_finite(fb)) .or. &
        any(.not. ieee_is_finite(z_ab1_mo_a)) .or. &
        any(.not. ieee_is_finite(z_ab1_mo_b)) .or. &
        any(.not. ieee_is_finite(dmo_energy)) .or. &
        any(.not. ieee_is_finite(dfa)) .or. &
        any(.not. ieee_is_finite(dfb)) .or. &
        any(.not. ieee_is_finite(dz_ab1_mo_a)) .or. &
        any(.not. ieee_is_finite(dz_ab1_mo_b))) then
      status = MRSF_Z_RESPONSE_NONFINITE
      return
    end if

    allocate(plus(lzdim), minus(lzdim))
    do coordinate = 1, ncoord
      ! Orbital-energy polarization at fixed Z and fixed remaining operator.
      call evaluate_mrsf_z_operator_action(z, &
        mo_energy + dmo_energy(:,coordinate), fa, fb, &
        z_ab1_mo_a, z_ab1_mo_b, nocca, noccb, plus)
      call evaluate_mrsf_z_operator_action(z, &
        mo_energy - dmo_energy(:,coordinate), fa, fb, &
        z_ab1_mo_a, z_ab1_mo_b, nocca, noccb, minus)
      operator_derivative_z(:,coordinate) = 0.5_dp*(plus - minus)

      ! Fock polarization at fixed Z.
      call evaluate_mrsf_z_operator_action(z, mo_energy, &
        fa + dfa(:,:,coordinate), fb + dfb(:,:,coordinate), &
        z_ab1_mo_a, z_ab1_mo_b, nocca, noccb, plus)
      call evaluate_mrsf_z_operator_action(z, mo_energy, &
        fa - dfa(:,:,coordinate), fb - dfb(:,:,coordinate), &
        z_ab1_mo_a, z_ab1_mo_b, nocca, noccb, minus)
      operator_derivative_z(:,coordinate) = &
        operator_derivative_z(:,coordinate) + 0.5_dp*(plus - minus)

      ! Response-image polarization at fixed Z.
      call evaluate_mrsf_z_operator_action(z, mo_energy, fa, fb, &
        z_ab1_mo_a + dz_ab1_mo_a(:,:,coordinate), &
        z_ab1_mo_b + dz_ab1_mo_b(:,:,coordinate), &
        nocca, noccb, plus)
      call evaluate_mrsf_z_operator_action(z, mo_energy, fa, fb, &
        z_ab1_mo_a - dz_ab1_mo_a(:,:,coordinate), &
        z_ab1_mo_b - dz_ab1_mo_b(:,:,coordinate), &
        nocca, noccb, minus)
      operator_derivative_z(:,coordinate) = &
        operator_derivative_z(:,coordinate) + 0.5_dp*(plus - minus)
    end do
    deallocate(plus, minus)

    if (any(.not. ieee_is_finite(operator_derivative_z))) then
      operator_derivative_z = 0.0_dp
      status = MRSF_Z_RESPONSE_NONFINITE
    end if
  end subroutine build_mrsf_z_operator_derivative_action

!###############################################################################

  subroutine solve_mrsf_z_response_from_mo_derivatives( &
      apply_orbital_hessian, target_multiplicity, is_umrsf, cam_flag, &
      nocca, noccb, mo_energy, fa, fb, hxa, hxb, rhs_ab1_mo_a, &
      rhs_ab1_mo_b, z_ab1_mo_a, z_ab1_mo_b, tij, tab, z, dmo_energy, &
      dfa, dfb, dhxa, dhxb, drhs_ab1_mo_a, drhs_ab1_mo_b, &
      dz_ab1_mo_a, dz_ab1_mo_b, dtij, dtab, drhs, &
      operator_derivative_z, dz, residual_max, status, tol, maxit, restart, &
      apply_orbital_hessian_batch,apply_orbital_preconditioner_batch)
    procedure(mrsf_z_orbital_hessian_action) :: apply_orbital_hessian
    procedure(mrsf_z_orbital_hessian_batch_action), optional :: &
      apply_orbital_hessian_batch
    procedure(mrsf_z_orbital_hessian_batch_action), optional :: &
      apply_orbital_preconditioner_batch
    integer, intent(in) :: target_multiplicity, nocca, noccb
    logical, intent(in) :: is_umrsf, cam_flag
    real(kind=dp), intent(in) :: mo_energy(:), fa(:,:), fb(:,:)
    real(kind=dp), intent(in) :: hxa(:,:), hxb(:,:)
    real(kind=dp), intent(in) :: rhs_ab1_mo_a(:,:), rhs_ab1_mo_b(:,:)
    real(kind=dp), intent(in) :: z_ab1_mo_a(:,:), z_ab1_mo_b(:,:)
    real(kind=dp), intent(in) :: tij(:,:), tab(:,:), z(:)
    real(kind=dp), intent(in) :: dmo_energy(:,:), dfa(:,:,:), dfb(:,:,:)
    real(kind=dp), intent(in) :: dhxa(:,:,:), dhxb(:,:,:)
    real(kind=dp), intent(in) :: drhs_ab1_mo_a(:,:,:)
    real(kind=dp), intent(in) :: drhs_ab1_mo_b(:,:,:)
    real(kind=dp), intent(in) :: dz_ab1_mo_a(:,:,:)
    real(kind=dp), intent(in) :: dz_ab1_mo_b(:,:,:)
    real(kind=dp), intent(in) :: dtij(:,:,:), dtab(:,:,:)
    real(kind=dp), intent(out) :: drhs(:,:), operator_derivative_z(:,:)
    real(kind=dp), intent(out) :: dz(:,:), residual_max
    integer, intent(out) :: status
    real(kind=dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxit, restart

    real(kind=dp) :: solve_tolerance
    integer :: solve_iterations, solve_restart

    dz = 0.0_dp
    residual_max = 0.0_dp
    call build_mrsf_z_rhs_derivative(target_multiplicity, is_umrsf, &
      cam_flag, nocca, noccb, fa, fb, hxa, hxb, rhs_ab1_mo_a, &
      rhs_ab1_mo_b, tij, tab, dfa, dfb, dhxa, dhxb, &
      drhs_ab1_mo_a, drhs_ab1_mo_b, dtij, dtab, drhs, status)
    if (status /= MRSF_Z_RESPONSE_SUCCESS) return

    call build_mrsf_z_operator_derivative_action(target_multiplicity, &
      is_umrsf, cam_flag, nocca, noccb, z, mo_energy, fa, fb, &
      z_ab1_mo_a, z_ab1_mo_b, dmo_energy, dfa, dfb, &
      dz_ab1_mo_a, dz_ab1_mo_b, operator_derivative_z, status)
    if (status /= MRSF_Z_RESPONSE_SUCCESS) return

    solve_tolerance = 1.0e-11_dp
    if (present(tol)) solve_tolerance = tol
    solve_iterations = max(200, 4*size(z))
    if (present(maxit)) solve_iterations = maxit
    solve_restart = min(max(12, size(z)/4), max(1, size(z)))
    if (present(restart)) solve_restart = restart
    if(present(apply_orbital_hessian_batch)) then
      if(present(apply_orbital_preconditioner_batch)) then
        call solve_mrsf_z_response_batch_matrix_free( &
          apply_orbital_hessian_batch,drhs,operator_derivative_z,dz, &
          residual_max,status,tol=solve_tolerance,maxit=solve_iterations, &
          apply_preconditioner=apply_orbital_preconditioner_batch)
      else
        call solve_mrsf_z_response_batch_matrix_free( &
          apply_orbital_hessian_batch,drhs,operator_derivative_z,dz, &
          residual_max,status,tol=solve_tolerance,maxit=solve_iterations)
      end if
    else
      call solve_mrsf_z_response_matrix_free(apply_orbital_hessian, drhs, &
        operator_derivative_z, dz, residual_max, status, &
        tol=solve_tolerance, maxit=solve_iterations, restart=solve_restart)
    end if
  end subroutine solve_mrsf_z_response_from_mo_derivatives

end module tdhf_mrsf_hessian_z_response_mod
