program test_mrsf_hessian_z_response

  use precision, only: dp
  use tdhf_sf_lib, only: sfrorhs, sfrolhs
  use tdhf_mrsf_hessian_z_response_mod, only: &
    MRSF_Z_RESPONSE_SUCCESS, MRSF_Z_RESPONSE_NOT_TWO_SOMO, &
    MRSF_Z_RESPONSE_UNSUPPORTED_MULTIPLICITY, &
    MRSF_Z_RESPONSE_UMRSF_UNSUPPORTED, MRSF_Z_RESPONSE_CAM_UNSUPPORTED, &
    build_mrsf_z_rhs_derivative, &
    build_mrsf_z_operator_derivative_action, &
    solve_mrsf_z_response_from_mo_derivatives

  implicit none

  integer, parameter :: nbf = 5, nocca = 3, noccb = 1
  integer, parameter :: nvira = nbf - nocca, nvirb = nbf - noccb
  integer, parameter :: nsomo = nocca - noccb
  integer, parameter :: lzdim = noccb*(nsomo + nvira) + nsomo*nvira
  integer, parameter :: ncoord = 2
  real(kind=dp), parameter :: polarization_scale = 0.37_dp

  real(kind=dp) :: mo_energy(nbf), fa(nbf,nbf), fb(nbf,nbf)
  real(kind=dp) :: hxa(nbf,nocca), hxb(nbf,nbf)
  real(kind=dp) :: rhs_ab1_mo_a(nocca,nvira)
  real(kind=dp) :: rhs_ab1_mo_b(noccb,nvirb)
  real(kind=dp) :: z_ab1_mo_a(nocca,nvira)
  real(kind=dp) :: z_ab1_mo_b(noccb,nvirb)
  real(kind=dp) :: tij(nocca,nocca), tab(nvirb,nvirb), z(lzdim)
  real(kind=dp) :: dmo_energy(nbf,ncoord)
  real(kind=dp) :: dfa(nbf,nbf,ncoord), dfb(nbf,nbf,ncoord)
  real(kind=dp) :: dhxa(nbf,nocca,ncoord), dhxb(nbf,nbf,ncoord)
  real(kind=dp) :: drhs_ab1_mo_a(nocca,nvira,ncoord)
  real(kind=dp) :: drhs_ab1_mo_b(noccb,nvirb,ncoord)
  real(kind=dp) :: dz_ab1_mo_a(nocca,nvira,ncoord)
  real(kind=dp) :: dz_ab1_mo_b(noccb,nvirb,ncoord)
  real(kind=dp) :: dtij(nocca,nocca,ncoord)
  real(kind=dp) :: dtab(nvirb,nvirb,ncoord)
  real(kind=dp) :: drhs(lzdim,ncoord), dhz(lzdim,ncoord)
  real(kind=dp) :: triplet_drhs(lzdim,ncoord), triplet_dhz(lzdim,ncoord)
  real(kind=dp) :: expected_drhs(lzdim,ncoord)
  real(kind=dp) :: expected_dhz(lzdim,ncoord)
  real(kind=dp) :: plus(lzdim), minus(lzdim)
  real(kind=dp) :: solved_drhs(lzdim,ncoord)
  real(kind=dp) :: solved_dhz(lzdim,ncoord)
  real(kind=dp) :: dz(lzdim,ncoord), expected_dz(lzdim,ncoord)
  real(kind=dp) :: orbital_diagonal(lzdim), residual
  integer :: coordinate, i, j, status

  call initialize_inputs()

  call build_mrsf_z_rhs_derivative(1, .false., .false., nocca, noccb, &
    fa, fb, hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, tij, tab, dfa, dfb, &
    dhxa, dhxb, drhs_ab1_mo_a, drhs_ab1_mo_b, dtij, dtab, drhs, status)
  if (status /= MRSF_Z_RESPONSE_SUCCESS) &
    error stop 'singlet RHS derivative rejected valid quantities'

  call build_mrsf_z_operator_derivative_action(1, .false., .false., &
    nocca, noccb, z, mo_energy, fa, fb, z_ab1_mo_a, z_ab1_mo_b, &
    dmo_energy, dfa, dfb, dz_ab1_mo_a, dz_ab1_mo_b, dhz, status)
  if (status /= MRSF_Z_RESPONSE_SUCCESS) &
    error stop 'singlet operator derivative rejected valid quantities'

  do coordinate = 1, ncoord
    call evaluate_rhs_direction(polarization_scale, coordinate, plus)
    call evaluate_rhs_direction(-polarization_scale, coordinate, minus)
    expected_drhs(:,coordinate) = &
      (plus - minus)/(2.0_dp*polarization_scale)

    call evaluate_operator_direction(polarization_scale, coordinate, plus)
    call evaluate_operator_direction(-polarization_scale, coordinate, minus)
    expected_dhz(:,coordinate) = &
      (plus - minus)/(2.0_dp*polarization_scale)
  end do
  if (maxval(abs(drhs - expected_drhs)) > 2.0e-12_dp) &
    error stop 'RHS derivative differs from the directional algebraic oracle'
  if (maxval(abs(dhz - expected_dhz)) > 2.0e-12_dp) &
    error stop 'operator derivative differs from the directional algebraic oracle'

  call build_mrsf_z_rhs_derivative(3, .false., .false., nocca, noccb, &
    fa, fb, hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, tij, tab, dfa, dfb, &
    dhxa, dhxb, drhs_ab1_mo_a, drhs_ab1_mo_b, dtij, dtab, &
    triplet_drhs, status)
  if (status /= MRSF_Z_RESPONSE_SUCCESS) &
    error stop 'triplet RHS derivative rejected valid quantities'
  call build_mrsf_z_operator_derivative_action(3, .false., .false., &
    nocca, noccb, z, mo_energy, fa, fb, z_ab1_mo_a, z_ab1_mo_b, &
    dmo_energy, dfa, dfb, dz_ab1_mo_a, dz_ab1_mo_b, triplet_dhz, status)
  if (status /= MRSF_Z_RESPONSE_SUCCESS) &
    error stop 'triplet operator derivative rejected valid quantities'
  if (maxval(abs(triplet_drhs - drhs)) > 1.0e-14_dp) &
    error stop 'singlet and triplet orbital RHS conventions differ'
  if (maxval(abs(triplet_dhz - dhz)) > 1.0e-14_dp) &
    error stop 'singlet and triplet orbital-operator conventions differ'

  do i = 1, lzdim
    orbital_diagonal(i) = (-1.0_dp)**i*(1.5_dp + 0.2_dp*i)
  end do
  do coordinate = 1, ncoord
    expected_dz(:,coordinate) = &
      (expected_drhs(:,coordinate) - expected_dhz(:,coordinate)) &
      / orbital_diagonal
  end do
  call solve_mrsf_z_response_from_mo_derivatives( &
    apply_orbital_hessian, 1, .false., .false., nocca, noccb, &
    mo_energy, fa, fb, hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, &
    z_ab1_mo_a, z_ab1_mo_b, tij, tab, z, dmo_energy, dfa, dfb, &
    dhxa, dhxb, drhs_ab1_mo_a, drhs_ab1_mo_b, dz_ab1_mo_a, &
    dz_ab1_mo_b, dtij, dtab, solved_drhs, solved_dhz, dz, residual, &
    status, tol=1.0e-13_dp, maxit=100, restart=lzdim)
  if (status /= MRSF_Z_RESPONSE_SUCCESS) &
    error stop 'matrix-free derivative Z-vector equation did not converge'
  if (maxval(abs(solved_drhs - expected_drhs)) > 2.0e-12_dp) &
    error stop 'combined solver changed the RHS derivative'
  if (maxval(abs(solved_dhz - expected_dhz)) > 2.0e-12_dp) &
    error stop 'combined solver changed the operator derivative'
  if (maxval(abs(dz - expected_dz)) > 2.0e-11_dp) &
    error stop 'matrix-free derivative Z vector is incorrect'
  if (residual > 2.0e-11_dp) &
    error stop 'matrix-free derivative Z residual is too large'

  call build_mrsf_z_rhs_derivative(5, .false., .false., nocca, noccb, &
    fa, fb, hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, tij, tab, dfa, dfb, &
    dhxa, dhxb, drhs_ab1_mo_a, drhs_ab1_mo_b, dtij, dtab, drhs, status)
  if (status /= MRSF_Z_RESPONSE_UNSUPPORTED_MULTIPLICITY) &
    error stop 'unsupported multiplicity did not fail closed'

  call build_mrsf_z_rhs_derivative(1, .true., .false., nocca, noccb, &
    fa, fb, hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, tij, tab, dfa, dfb, &
    dhxa, dhxb, drhs_ab1_mo_a, drhs_ab1_mo_b, dtij, dtab, drhs, status)
  if (status /= MRSF_Z_RESPONSE_UMRSF_UNSUPPORTED) &
    error stop 'unrestricted MRSF did not fail closed'

  call build_mrsf_z_rhs_derivative(1, .false., .true., nocca, noccb, &
    fa, fb, hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, tij, tab, dfa, dfb, &
    dhxa, dhxb, drhs_ab1_mo_a, drhs_ab1_mo_b, dtij, dtab, drhs, status)
  if (status /= MRSF_Z_RESPONSE_CAM_UNSUPPORTED) &
    error stop 'range-separated response did not fail closed'

  call build_mrsf_z_rhs_derivative(1, .false., .false., 2, 1, &
    fa, fb, hxa, hxb, rhs_ab1_mo_a, rhs_ab1_mo_b, tij, tab, dfa, dfb, &
    dhxa, dhxb, drhs_ab1_mo_a, drhs_ab1_mo_b, dtij, dtab, drhs, status)
  if (status /= MRSF_Z_RESPONSE_NOT_TWO_SOMO) &
    error stop 'non-two-SOMO reference did not fail closed'

contains

  subroutine initialize_inputs()
    mo_energy = [-1.10_dp, -0.63_dp, -0.41_dp, 0.27_dp, 0.82_dp]
    z = [(0.06_dp*i - 0.23_dp, i=1,lzdim)]
    do j = 1, nbf
      do i = 1, nbf
        fa(i,j) = 0.07_dp*i - 0.03_dp*j + 0.004_dp*i*j
        fb(i,j) = -0.02_dp*i + 0.05_dp*j - 0.003_dp*i*j
        hxb(i,j) = 0.011_dp*i + 0.017_dp*j
        do coordinate = 1, ncoord
          dfa(i,j,coordinate) = &
            0.005_dp*i - 0.002_dp*j + 0.001_dp*coordinate
          dfb(i,j,coordinate) = &
            -0.003_dp*i + 0.004_dp*j - 0.0007_dp*coordinate
          dhxb(i,j,coordinate) = &
            0.002_dp*i - 0.001_dp*j + 0.0003_dp*coordinate
        end do
      end do
    end do
    do j = 1, nocca
      do i = 1, nbf
        hxa(i,j) = -0.013_dp*i + 0.009_dp*j
        do coordinate = 1, ncoord
          dhxa(i,j,coordinate) = &
            0.0015_dp*i + 0.0008_dp*j - 0.0002_dp*coordinate
        end do
      end do
    end do
    do j = 1, nocca
      do i = 1, nocca
        tij(i,j) = 0.021_dp*i - 0.016_dp*j
        do coordinate = 1, ncoord
          dtij(i,j,coordinate) = &
            0.001_dp*i + 0.0005_dp*j + 0.0002_dp*coordinate
        end do
      end do
    end do
    do j = 1, nvirb
      do i = 1, nvirb
        tab(i,j) = -0.018_dp*i + 0.012_dp*j
        do coordinate = 1, ncoord
          dtab(i,j,coordinate) = &
            -0.0009_dp*i + 0.0006_dp*j - 0.0001_dp*coordinate
        end do
      end do
    end do
    do j = 1, nvira
      do i = 1, nocca
        rhs_ab1_mo_a(i,j) = 0.031_dp*i - 0.014_dp*j
        z_ab1_mo_a(i,j) = -0.022_dp*i + 0.019_dp*j
        do coordinate = 1, ncoord
          drhs_ab1_mo_a(i,j,coordinate) = &
            0.0013_dp*i - 0.0004_dp*j + 0.0002_dp*coordinate
          dz_ab1_mo_a(i,j,coordinate) = &
            -0.0008_dp*i + 0.0007_dp*j - 0.0001_dp*coordinate
        end do
      end do
    end do
    do j = 1, nvirb
      rhs_ab1_mo_b(1,j) = -0.017_dp + 0.008_dp*j
      z_ab1_mo_b(1,j) = 0.026_dp - 0.006_dp*j
      do coordinate = 1, ncoord
        drhs_ab1_mo_b(1,j,coordinate) = &
          -0.0006_dp + 0.0003_dp*j + 0.0001_dp*coordinate
        dz_ab1_mo_b(1,j,coordinate) = &
          0.0005_dp - 0.0002_dp*j + 0.00015_dp*coordinate
      end do
    end do
    do coordinate = 1, ncoord
      do i = 1, nbf
        dmo_energy(i,coordinate) = &
          0.003_dp*i - 0.001_dp*coordinate
      end do
    end do
  end subroutine initialize_inputs

  subroutine evaluate_rhs_direction(scale, direction, value)
    real(kind=dp), intent(in) :: scale
    integer, intent(in) :: direction
    real(kind=dp), intent(out) :: value(lzdim)
    real(kind=dp) :: hxa_work(nbf,nocca), hxb_work(nbf,nbf)

    hxa_work = hxa + scale*dhxa(:,:,direction)
    hxb_work = hxb + scale*dhxb(:,:,direction)
    call sfrorhs(value, hxa_work, hxb_work, &
      rhs_ab1_mo_a + scale*drhs_ab1_mo_a(:,:,direction), &
      rhs_ab1_mo_b + scale*drhs_ab1_mo_b(:,:,direction), &
      tij + scale*dtij(:,:,direction), &
      tab + scale*dtab(:,:,direction), &
      fa + scale*dfa(:,:,direction), fb + scale*dfb(:,:,direction), &
      nocca, noccb)
  end subroutine evaluate_rhs_direction

  subroutine evaluate_operator_direction(scale, direction, value)
    real(kind=dp), intent(in) :: scale
    integer, intent(in) :: direction
    real(kind=dp), intent(out) :: value(lzdim)

    call sfrolhs(value, z, mo_energy + scale*dmo_energy(:,direction), &
      fa + scale*dfa(:,:,direction), fb + scale*dfb(:,:,direction), &
      z_ab1_mo_a + scale*dz_ab1_mo_a(:,:,direction), &
      z_ab1_mo_b + scale*dz_ab1_mo_b(:,:,direction), nocca, noccb)
  end subroutine evaluate_operator_direction

  subroutine apply_orbital_hessian(vector, result, operator_status)
    real(kind=dp), intent(in) :: vector(:)
    real(kind=dp), intent(out) :: result(:)
    integer, intent(out) :: operator_status

    if (size(vector) /= lzdim .or. size(result) /= lzdim) then
      result = 0.0_dp
      operator_status = -1
      return
    end if
    result = orbital_diagonal*vector
    operator_status = 0
  end subroutine apply_orbital_hessian

end program test_mrsf_hessian_z_response
