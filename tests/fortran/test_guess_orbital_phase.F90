program test_guess_orbital_phase
  use precision, only: dp
  use orbital_phase_test_m, only: canonicalize_orbital_phases
  implicit none

  real(dp) :: orbitals(4,4), once(4,4)

  orbitals(:,1) = [0.1_dp, -0.7_dp, 0.3_dp, 0.0_dp]
  orbitals(:,2) = [-0.5_dp, 0.5_dp + epsilon(1.0_dp), 0.1_dp, 0.0_dp]
  orbitals(:,3) = 0.0_dp
  orbitals(:,4) = [0.5_dp, -0.5_dp - epsilon(1.0_dp), 0.0_dp, 0.0_dp]

  call canonicalize_orbital_phases(orbitals)

  if (maxval(abs(orbitals(:,1) - &
      [-0.1_dp, 0.7_dp, -0.3_dp, 0.0_dp])) > 1.0e-15_dp) error stop 1
  ! Both near-ties select AO 1, rather than following the round-off-level peak.
  if (orbitals(1,2) <= 0.0_dp) error stop 2
  if (orbitals(1,4) <= 0.0_dp) error stop 3
  ! A zero orbital has no sign to choose.
  if (any(orbitals(:,3) /= 0.0_dp)) error stop 4

  once = orbitals
  call canonicalize_orbital_phases(orbitals)
  if (any(orbitals /= once)) error stop 5

  write(*,'(a)') 'orbital phase canonicalization selftest passed'
end program test_guess_orbital_phase
