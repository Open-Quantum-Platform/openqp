program test_der2_recursion
  ! Validate the production 1D second-derivative recursion der2_kinovl_xyz by
  ! checking it equals the bra-center first-derivative recursion der_kinovl_xyz
  ! composed with itself, on a random input array. This pins the closed-form
  !   d2[j,i] = 4 ai^2 [j,i+2] - 2 ai (2i+1) [j,i] + i(i-1) [j,i-2]
  ! against d(d[.]).
  use, intrinsic :: iso_fortran_env, only: real64
  use mod_1e_primitives, only: der_kinovl_xyz, der2_kinovl_xyz, &
                               der_coul_xyz, der2_coul_xyz
  implicit none

  integer, parameter :: ni = 4, nj = 3
  real(real64) :: xyz(0:nj,0:ni+2,3)
  real(real64) :: d1(0:nj,0:ni+2,3)
  real(real64) :: dd(0:nj,0:ni,3)
  real(real64) :: d2(0:nj,0:ni,3)
  real(real64) :: ai, maxerr
  integer :: seed_size
  integer, allocatable :: seed(:)

  ! Coulomb (Rys) arrays carry an extra root dimension
  integer, parameter :: nr = 2
  real(real64) :: cxyz(0:nj,0:ni+2,3,nr)
  real(real64) :: cd1(0:nj,0:ni+2,3,nr)
  real(real64) :: cdd(0:nj,0:ni,3,nr)
  real(real64) :: cd2(0:nj,0:ni,3,nr)
  real(real64) :: cerr

  call random_seed(size=seed_size)
  allocate(seed(seed_size)); seed = 20260529
  call random_seed(put=seed)
  call random_number(xyz)
  call random_number(cxyz)
  ai = 1.37_real64

  ! --- overlap/kinetic recursion: der2 == der(der) ---
  call der_kinovl_xyz(d1, xyz, ni+1, nj, ai)
  call der_kinovl_xyz(dd, d1, ni, nj, ai)
  call der2_kinovl_xyz(d2, xyz, ni, nj, ai)
  maxerr = maxval(abs(dd(0:nj,0:ni,:) - d2(0:nj,0:ni,:)))
  print '(a,es12.4)', 'kinovl  max |der(der) - der2| = ', maxerr

  ! --- Coulomb (Rys) recursion: der2 == der(der), per root ---
  call der_coul_xyz(cd1, cxyz, ni+1, nj, ai, nr)
  call der_coul_xyz(cdd, cd1, ni, nj, ai, nr)
  call der2_coul_xyz(cd2, cxyz, ni, nj, ai, nr)
  cerr = maxval(abs(cdd(0:nj,0:ni,:,:) - cd2(0:nj,0:ni,:,:)))
  print '(a,es12.4)', 'coulomb max |der(der) - der2| = ', cerr

  if (max(maxerr, cerr) < 1.0e-10_real64) then
    print '(a)', 'PASS: der2 routines match composed first derivatives'
  else
    print '(a)', 'FAIL: der2 routines disagree with composed first derivatives'
    call exit(1)
  end if

end program test_der2_recursion
