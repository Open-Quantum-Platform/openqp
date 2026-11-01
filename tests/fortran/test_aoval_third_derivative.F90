program test_aoval_third_derivative

  use, intrinsic :: iso_c_binding, only: c_int64_t
  use precision, only: fp
  use atomic_structure_m, only: atomic_structure
  use basis_tools, only: basis_set

  implicit none

  type(atomic_structure), target :: atoms
  type(basis_set) :: basis
  real(fp), parameter :: h = 1.0e-5_fp
  real(fp) :: point(3), shifted(3), maxerr
  real(fp) :: value(3), g1(3,3), g2(3,6), g3(3,10)
  real(fp) :: vp(3), vm(3), gp(3,3), gm(3,3), hp(3,6), hm(3,6)
  real(fp) :: fd3(3,10)
  integer :: axis, ierr, nnz

  ierr = atoms%init(1_c_int64_t)
  if (ierr /= 0) error stop 'failed to allocate test atom'
  atoms%xyz = 0.0_fp
  atoms%zn = 1.0_fp
  atoms%mass = 1.0_fp
  atoms%grad = 0.0_fp

  ! One unnormalized Cartesian p shell is sufficient to exercise pure and
  ! mixed third derivatives without depending on basis-library I/O.
  call basis%reserve(1, 1, 3)
  basis%atoms => atoms
  basis%nshell = 1
  basis%nprim = 1
  basis%nbf = 3
  basis%mxcontr = 1
  basis%mxam = 1
  basis%g_offset(1) = 1
  basis%origin(1) = 1
  basis%am(1) = 1
  basis%harmonic(1) = 0
  basis%ncontr(1) = 1
  basis%ao_offset(1) = 1
  basis%naos(1) = 3
  basis%ex(1) = 0.73_fp
  basis%cc(1) = 1.17_fp
  basis%bfnrm = 1.0_fp
  call basis%set_screening(50.0_fp)

  point = [0.31_fp, -0.27_fp, 0.19_fp]
  call eval3(point, value, g1, g2, g3)

  do axis = 1, 3
    shifted = point
    shifted(axis) = shifted(axis) + h
    call eval2(shifted, vp, gp, hp)
    shifted(axis) = point(axis) - h
    call eval2(shifted, vm, gm, hm)
    select case (axis)
    case (1)
      fd3(:,1) = (hp(:,1)-hm(:,1))/(2*h)   ! XXX
      fd3(:,6) = (hp(:,2)-hm(:,2))/(2*h)   ! YYX
      fd3(:,8) = (hp(:,3)-hm(:,3))/(2*h)   ! ZZX
    case (2)
      fd3(:,2) = (hp(:,2)-hm(:,2))/(2*h)   ! YYY
      fd3(:,4) = (hp(:,1)-hm(:,1))/(2*h)   ! XXY
      fd3(:,9) = (hp(:,3)-hm(:,3))/(2*h)   ! ZZY
    case (3)
      fd3(:,3) = (hp(:,3)-hm(:,3))/(2*h)   ! ZZZ
      fd3(:,5) = (hp(:,1)-hm(:,1))/(2*h)   ! XXZ
      fd3(:,7) = (hp(:,2)-hm(:,2))/(2*h)   ! YYZ
      fd3(:,10) = (hp(:,4)-hm(:,4))/(2*h)  ! XYZ = d(XY)/dZ
    end select
  end do

  maxerr = maxval(abs(g3-fd3))
  print '(a,es12.4)', 'AO G3 max |analytic - FD(G2)| = ', maxerr
  if (maxerr > 2.0e-8_fp) error stop 'AO third derivatives failed finite-difference check'

contains

  subroutine eval2(xyz, v, g, hess)
    real(fp), intent(in) :: xyz(3)
    real(fp), intent(out) :: v(3), g(3,3), hess(3,6)
    call basis%aoval(xyz, nnz, v, g(:,1), g(:,2), g(:,3), &
      hess(:,1), hess(:,2), hess(:,3), hess(:,4), hess(:,5), hess(:,6))
  end subroutine eval2

  subroutine eval3(xyz, v, g, hess, third)
    real(fp), intent(in) :: xyz(3)
    real(fp), intent(out) :: v(3), g(3,3), hess(3,6), third(3,10)
    call basis%aoval(xyz, nnz, v, g(:,1), g(:,2), g(:,3), &
      hess(:,1), hess(:,2), hess(:,3), hess(:,4), hess(:,5), hess(:,6), &
      third(:,1), third(:,2), third(:,3), third(:,4), third(:,5), &
      third(:,6), third(:,7), third(:,8), third(:,9), third(:,10))
  end subroutine eval3

end program test_aoval_third_derivative
