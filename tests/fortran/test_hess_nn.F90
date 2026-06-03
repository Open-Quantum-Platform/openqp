program test_hess_nn
  ! Link test for the native nuclear-repulsion Hessian (grd1::hess_nn).
  ! Validates the analytic (3N,3N) nuclear-repulsion Hessian against a central
  ! finite difference of the nuclear-repulsion gradient (grd1::grad_nn), using
  ! the production routines linked from liboqp.
  use precision, only: dp
  use atomic_structure_m, only: atomic_structure
  use grd1, only: grad_nn, hess_nn
  implicit none

  integer, parameter :: nat = 3
  type(atomic_structure) :: atoms
  integer :: ecp(nat)
  real(dp) :: hess_an(3*nat,3*nat), hess_fd(3*nat,3*nat)
  real(dp) :: gp(3,nat), gm(3,nat)
  real(dp) :: h, maxerr
  integer :: a, k, ka
  logical :: ok

  ok = atoms%init(nat)
  atoms%xyz(:,1) = [ 0.0_dp,  0.0_dp,  0.20_dp]
  atoms%xyz(:,2) = [ 0.0_dp,  1.45_dp, -0.95_dp]
  atoms%xyz(:,3) = [ 0.0_dp, -1.45_dp, -0.95_dp]
  atoms%zn = [8.0_dp, 1.0_dp, 1.0_dp]
  ecp = 0

  hess_an = 0.0_dp
  call hess_nn(atoms, ecp, hess_an)

  h = 1.0e-5_dp
  do k = 1, nat
    do a = 1, 3
      ka = 3*(k-1) + a
      atoms%xyz(a,k) = atoms%xyz(a,k) + h
      atoms%grad = 0.0_dp; call grad_nn(atoms, ecp); gp = atoms%grad
      atoms%xyz(a,k) = atoms%xyz(a,k) - 2*h
      atoms%grad = 0.0_dp; call grad_nn(atoms, ecp); gm = atoms%grad
      atoms%xyz(a,k) = atoms%xyz(a,k) + h
      hess_fd(:, ka) = reshape((gp - gm) / (2*h), [3*nat])
    end do
  end do

  maxerr = maxval(abs(hess_an - hess_fd))
  print '(a,es12.4)', 'max |analytic - FD| = ', maxerr
  ! Hessian symmetry of the analytic result
  print '(a,es12.4)', 'max asymmetry       = ', maxval(abs(hess_an - transpose(hess_an)))

  if (maxerr < 1.0e-6_dp) then
    print '(a)', 'PASS: hess_nn matches finite-difference nuclear-repulsion gradient'
  else
    print '(a)', 'FAIL: hess_nn disagrees with finite difference'
    call exit(1)
  end if

end program test_hess_nn
