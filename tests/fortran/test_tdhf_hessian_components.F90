program test_tdhf_hessian_components

  use precision, only: dp
  use tdhf_hessian_components_mod, only: assemble_tdhf_cartesian_hessian, &
    tdhf_hessian_is_applicable

  implicit none

  real(kind=dp) :: hfixed(2,2), hxc(2,2), rows(2,2)
  real(kind=dp) :: htotal(2,2), expected(2,2), asymmetry
  integer :: status

  hfixed = reshape([1.11_dp, 2.22_dp, 2.22_dp, 4.44_dp], [2,2])
  hxc = 0.0_dp
  rows = reshape([0.0_dp, 2.0_dp, 4.0_dp, 0.0_dp], [2,2])
  expected = reshape([1.11_dp, 5.22_dp, 5.22_dp, 4.44_dp], [2,2])

  call assemble_tdhf_cartesian_hessian(hfixed, hxc, rows, htotal, &
                                        asymmetry, status)
  if (status /= 0) error stop 'valid component shapes were rejected'
  if (maxval(abs(htotal - expected)) > 1.0e-14_dp) &
    error stop 'TDHF Hessian component sum is incorrect'
  if (abs(asymmetry - 2.0_dp) > 1.0e-14_dp) &
    error stop 'directional-row asymmetry is incorrect'

  if (.not. tdhf_hessian_is_applicable(1, 1, .false., .false., 1)) &
    error stop 'verified closed-shell full-response case was rejected'
  if (tdhf_hessian_is_applicable(2, 1, .false., .false., 1)) &
    error stop 'open-shell reference was accepted'
  if (tdhf_hessian_is_applicable(1, 1, .true., .false., 1)) &
    error stop 'Tamm-Dancoff response was accepted'
  if (tdhf_hessian_is_applicable(1, 3, .false., .false., 1)) &
    error stop 'triplet target was accepted'
  if (tdhf_hessian_is_applicable(1, 1, .false., .true., 1)) &
    error stop 'kinetic-energy-density functional was accepted'
  if (tdhf_hessian_is_applicable(1, 1, .false., .false., 2)) &
    error stop 'unverified multi-rank evaluation was accepted'

end program test_tdhf_hessian_components
