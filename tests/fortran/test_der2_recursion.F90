program test_der2_recursion
  ! Validate the production 1D second-derivative recursion der2_kinovl_xyz by
  ! checking it equals the bra-center first-derivative recursion der_kinovl_xyz
  ! composed with itself, on a random input array. This pins the closed-form
  !   d2[j,i] = 4 ai^2 [j,i+2] - 2 ai (2i+1) [j,i] + i(i-1) [j,i-2]
  ! against d(d[.]).
  use, intrinsic :: iso_fortran_env, only: real64
  use mod_1e_primitives, only: der_kinovl_xyz, der2_kinovl_xyz
  implicit none

  integer, parameter :: ni = 4, nj = 3
  real(real64) :: xyz(0:nj,0:ni+2,3)
  real(real64) :: d1(0:nj,0:ni+2,3)
  real(real64) :: dd(0:nj,0:ni,3)
  real(real64) :: d2(0:nj,0:ni,3)
  real(real64) :: ai, maxerr
  integer :: seed_size
  integer, allocatable :: seed(:)

  call random_seed(size=seed_size)
  allocate(seed(seed_size)); seed = 20260529
  call random_seed(put=seed)
  call random_number(xyz)
  ai = 1.37_real64

  ! compose first derivative twice: d1 up to bra index ni+1, then dd up to ni
  call der_kinovl_xyz(d1, xyz, ni+1, nj, ai)
  call der_kinovl_xyz(dd, d1, ni, nj, ai)

  ! closed-form second derivative
  call der2_kinovl_xyz(d2, xyz, ni, nj, ai)

  maxerr = maxval(abs(dd(0:nj,0:ni,:) - d2(0:nj,0:ni,:)))
  print '(a,es12.4)', 'max |der(der) - der2| = ', maxerr
  if (maxerr < 1.0e-10_real64) then
    print '(a)', 'PASS: der2_kinovl_xyz matches composed first derivatives'
  else
    print '(a)', 'FAIL: der2_kinovl_xyz disagrees with composed first derivatives'
    call exit(1)
  end if

end program test_der2_recursion
