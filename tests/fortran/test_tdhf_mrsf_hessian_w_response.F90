program test_tdhf_mrsf_hessian_w_response

  use precision, only: dp
  use tdhf_mrsf_hessian_w_response_mod, only: &
      MRSF_W_SUCCESS, MRSF_W_BAD_REFERENCE, MRSF_W_BAD_TARGET, &
      MRSF_W_UMRSF_UNSUPPORTED, MRSF_W_CAM_UNSUPPORTED, &
      build_mrsf_w_ao, build_mrsf_w_ao_derivative

  implicit none

  integer, parameter :: nbf = 6, noca = 4, nocb = 2, npert = 2
  integer, parameter :: nz = nocb*(2 + nbf - noca) + 2*(nbf - noca)

  type :: state_input_t
    real(kind=dp) :: mo(nbf,nbf), mo_energy(nbf)
    real(kind=dp) :: fa(nbf,nbf), fb(nbf,nbf), z(nz)
    real(kind=dp) :: hxa(nbf,noca), hxb(nbf,nbf)
    real(kind=dp) :: ppija(noca,noca), ppijb(nocb,nocb)
  end type state_input_t

  type(state_input_t) :: baseline
  real(kind=dp) :: dmo(nbf,nbf,npert), dmo_energy(nbf,npert)
  real(kind=dp) :: dfa(nbf,nbf,npert), dfb(nbf,nbf,npert), dz(nz,npert)
  real(kind=dp) :: dhxa(nbf,noca,npert), dhxb(nbf,nbf,npert)
  real(kind=dp) :: dppija(noca,noca,npert), dppijb(nocb,nocb,npert)
  real(kind=dp) :: w_ao(nbf,nbf), dw_ao(nbf,nbf,npert)
  real(kind=dp) :: triplet_dw(nbf,nbf,npert), reference(nbf,nbf)
  real(kind=dp) :: frozen_w_ao(nbf,nbf)
  integer :: perturbation, status

  call initialize_inputs()
  frozen_w_ao = reshape([ &
      -3.43814924999999915e-03_dp, -7.20363649999999911e-03_dp, &
      -1.94349674999999972e-02_dp, -1.85025297499999998e-02_dp, &
       3.87753300000000028e-03_dp,  7.41365824999999992e-03_dp, &
      -7.20363649999999911e-03_dp, -7.96893249999999949e-03_dp, &
      -2.16014909999999970e-02_dp, -2.23503120000000076e-02_dp, &
       6.41192000000000281e-04_dp,  3.13914599999999804e-03_dp, &
      -1.94349674999999972e-02_dp, -2.16014909999999970e-02_dp, &
      -1.99129457500000009e-02_dp, -2.40575692500000042e-02_dp, &
      -1.19229275000000325e-03_dp,  6.47946249999994603e-04_dp, &
      -1.85025297499999998e-02_dp, -2.23503120000000076e-02_dp, &
      -2.40575692500000042e-02_dp, -2.86374015000000065e-02_dp, &
      -6.33588375000000096e-03_dp, -5.43401600000000505e-03_dp, &
       3.87753300000000028e-03_dp,  6.41192000000000281e-04_dp, &
      -1.19229275000000325e-03_dp, -6.33588375000000096e-03_dp, &
      -2.27292475000000129e-03_dp,  1.23077174999999730e-03_dp, &
       7.41365824999999992e-03_dp,  3.13914599999999804e-03_dp, &
       6.47946249999994603e-04_dp, -5.43401600000000505e-03_dp, &
       1.23077174999999730e-03_dp,  8.29344999999927327e-05_dp ], &
      [nbf,nbf])

  call evaluate_state(baseline, 1, w_ao, status)
  call require(status == MRSF_W_SUCCESS, 'singlet baseline was rejected')
  call require(maxval(abs(w_ao - frozen_w_ao)) < 5.0e-15_dp, &
               'baseline differs from the frozen mrsfrowcal oracle')

  call evaluate_derivative(.true., 3, 1, .false., .false., dw_ao, status)
  call require(status == MRSF_W_SUCCESS, 'singlet derivative was rejected')
  do perturbation = 1, npert
    reference = 0.0_dp
    call accumulate_separated_polarization(perturbation, reference)
    call require(maxval(abs(dw_ao(:,:,perturbation) - reference)) < 2.0e-13_dp, &
                 'analytic derivative differs from separated polarization')
  end do

  call evaluate_derivative(.true., 3, 3, .false., .false., triplet_dw, status)
  call require(status == MRSF_W_SUCCESS, 'triplet derivative was rejected')
  call require(maxval(abs(triplet_dw - dw_ao)) < 2.0e-14_dp, &
               'common singlet/triplet row mapping changed')

  call evaluate_derivative(.false., 3, 1, .false., .false., triplet_dw, status)
  call require(status == MRSF_W_BAD_REFERENCE, 'non-ROHF reference was accepted')
  call require(maxval(abs(triplet_dw)) < tiny(1.0_dp), &
               'rejected reference returned data')

  call evaluate_derivative(.true., 3, 5, .false., .false., triplet_dw, status)
  call require(status == MRSF_W_BAD_TARGET, 'unsupported target was accepted')
  call require(maxval(abs(triplet_dw)) < tiny(1.0_dp), &
               'rejected target returned data')

  call evaluate_derivative(.true., 3, 1, .true., .false., triplet_dw, status)
  call require(status == MRSF_W_UMRSF_UNSUPPORTED, 'UMRSF case was accepted')
  call require(maxval(abs(triplet_dw)) < tiny(1.0_dp), &
               'rejected UMRSF case returned data')

  call evaluate_derivative(.true., 3, 1, .false., .true., triplet_dw, status)
  call require(status == MRSF_W_CAM_UNSUPPORTED, 'CAM case was accepted')
  call require(maxval(abs(triplet_dw)) < tiny(1.0_dp), &
               'rejected CAM case returned data')

contains

  subroutine initialize_inputs()
    integer :: i, j, k, p

    do j = 1, nbf
      do i = 1, nbf
        baseline%mo(i,j) = 0.03_dp*real(i,dp) - 0.02_dp*real(j,dp)
        if (i == j) baseline%mo(i,j) = baseline%mo(i,j) + 1.0_dp
        baseline%fa(i,j) = 0.04_dp*real(i,dp) + 0.015_dp*real(j,dp)
        baseline%fb(i,j) = -0.025_dp*real(i,dp) + 0.035_dp*real(j,dp)
        baseline%hxb(i,j) = 0.012_dp*real(i,dp) - 0.009_dp*real(j,dp)
        do p = 1, npert
          dmo(i,j,p) = real(p,dp)*(0.001_dp*real(i,dp) - 0.0007_dp*real(j,dp))
          dfa(i,j,p) = real(p,dp)*(0.0013_dp*real(i,dp) + 0.0004_dp*real(j,dp))
          dfb(i,j,p) = real(p,dp)*(-0.0008_dp*real(i,dp) + 0.0006_dp*real(j,dp))
          dhxb(i,j,p) = real(p,dp)*(0.0003_dp*real(i,dp) - 0.0002_dp*real(j,dp))
        end do
      end do
    end do

    do i = 1, nbf
      baseline%mo_energy(i) = -0.7_dp + 0.16_dp*real(i,dp)
      do p = 1, npert
        dmo_energy(i,p) = real(p,dp)*(0.0009_dp - 0.0001_dp*real(i,dp))
      end do
    end do
    do k = 1, nz
      baseline%z(k) = 0.011_dp*real(k,dp) - 0.05_dp
      do p = 1, npert
        dz(k,p) = real(p,dp)*(0.0005_dp*real(k,dp) - 0.001_dp)
      end do
    end do

    do j = 1, noca
      do i = 1, nbf
        baseline%hxa(i,j) = -0.008_dp*real(i,dp) + 0.013_dp*real(j,dp)
        do p = 1, npert
          dhxa(i,j,p) = real(p,dp)*(-0.0002_dp*real(i,dp) + 0.0007_dp*real(j,dp))
        end do
      end do
      do i = 1, noca
        baseline%ppija(i,j) = 0.006_dp*real(i,dp) + 0.004_dp*real(j,dp)
        do p = 1, npert
          dppija(i,j,p) = real(p,dp)*(0.0002_dp*real(i+j,dp) - 0.0003_dp)
        end do
      end do
    end do
    do j = 1, nocb
      do i = 1, nocb
        baseline%ppijb(i,j) = -0.005_dp*real(i,dp) + 0.003_dp*real(j,dp)
        do p = 1, npert
          dppijb(i,j,p) = real(p,dp)*(-0.0001_dp*real(i,dp) + 0.0002_dp*real(j,dp))
        end do
      end do
    end do
  end subroutine initialize_inputs


  subroutine evaluate_state(state, target_multiplicity, output, state_status)
    type(state_input_t), intent(in) :: state
    integer, intent(in) :: target_multiplicity
    real(kind=dp), intent(out) :: output(nbf,nbf)
    integer, intent(out) :: state_status

    call build_mrsf_w_ao( &
        state%mo, state%mo_energy, state%fa, state%fb, state%z, &
        state%hxa, state%hxb, state%ppija, state%ppijb, noca, nocb, &
        .true., 3, target_multiplicity, .false., .false., output, state_status)
  end subroutine evaluate_state


  subroutine evaluate_derivative(rohf, reference_mult, target_mult, umrsf, &
                                 cam_flag, output, derivative_status)
    logical, intent(in) :: rohf, umrsf, cam_flag
    integer, intent(in) :: reference_mult, target_mult
    real(kind=dp), intent(out) :: output(nbf,nbf,npert)
    integer, intent(out) :: derivative_status

    call build_mrsf_w_ao_derivative( &
        baseline%mo, dmo, baseline%mo_energy, dmo_energy, &
        baseline%fa, dfa, baseline%fb, dfb, baseline%z, dz, &
        baseline%hxa, dhxa, baseline%hxb, dhxb, &
        baseline%ppija, dppija, baseline%ppijb, dppijb, noca, nocb, &
        rohf, reference_mult, target_mult, umrsf, cam_flag, &
        output, derivative_status)
  end subroutine evaluate_derivative


  subroutine accumulate_separated_polarization(perturbation, total)
    integer, intent(in) :: perturbation
    real(kind=dp), intent(inout) :: total(nbf,nbf)

    type(state_input_t) :: positive, negative
    real(kind=dp) :: w_positive(nbf,nbf), w_negative(nbf,nbf)
    integer :: group, positive_status, negative_status

    do group = 1, 9
      positive = baseline
      negative = baseline
      select case (group)
      case (1)
        positive%mo = positive%mo + dmo(:,:,perturbation)
        negative%mo = negative%mo - dmo(:,:,perturbation)
      case (2)
        positive%mo_energy = positive%mo_energy + dmo_energy(:,perturbation)
        negative%mo_energy = negative%mo_energy - dmo_energy(:,perturbation)
      case (3)
        positive%fa = positive%fa + dfa(:,:,perturbation)
        negative%fa = negative%fa - dfa(:,:,perturbation)
      case (4)
        positive%fb = positive%fb + dfb(:,:,perturbation)
        negative%fb = negative%fb - dfb(:,:,perturbation)
      case (5)
        positive%z = positive%z + dz(:,perturbation)
        negative%z = negative%z - dz(:,perturbation)
      case (6)
        positive%hxa = positive%hxa + dhxa(:,:,perturbation)
        negative%hxa = negative%hxa - dhxa(:,:,perturbation)
      case (7)
        positive%hxb = positive%hxb + dhxb(:,:,perturbation)
        negative%hxb = negative%hxb - dhxb(:,:,perturbation)
      case (8)
        positive%ppija = positive%ppija + dppija(:,:,perturbation)
        negative%ppija = negative%ppija - dppija(:,:,perturbation)
      case (9)
        positive%ppijb = positive%ppijb + dppijb(:,:,perturbation)
        negative%ppijb = negative%ppijb - dppijb(:,:,perturbation)
      end select

      call evaluate_state(positive, 1, w_positive, positive_status)
      call evaluate_state(negative, 1, w_negative, negative_status)
      call require(positive_status == MRSF_W_SUCCESS, 'positive polarization failed')
      call require(negative_status == MRSF_W_SUCCESS, 'negative polarization failed')
      total = total + 0.5_dp*(w_positive - w_negative)
    end do
  end subroutine accumulate_separated_polarization


  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_tdhf_mrsf_hessian_w_response
