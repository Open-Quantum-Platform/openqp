program test_mrsf_hessian_space

  use precision, only: dp
  use tdhf_mrsf_conventions_mod, only: mrsf_raw_spc_multiplier
  use tdhf_mrsf_hessian_space_mod, only: build_mrsf_packed_transform, &
    mrsf_physical_dimensions

  implicit none

  real(kind=dp) :: singlet(9,8),triplet(9,6)
  real(kind=dp) :: gram_s(8,8),gram_t(6,6)
  integer :: expanded,physical,status

  if (mrsf_raw_spc_multiplier(1) /= 1 .or. &
      mrsf_raw_spc_multiplier(3) /= -1 .or. &
      mrsf_raw_spc_multiplier(5) /= 0) &
    error stop 'raw-to-physical SPC sign conversion is incorrect'

  call mrsf_physical_dimensions(4,3,1,1,expanded,physical,status)
  if (status /= 0 .or. expanded /= 9 .or. physical /= 8) &
    error stop 'singlet MRSF dimensions are incorrect'
  call build_mrsf_packed_transform(4,3,1,1,singlet,status)
  if (status /= 0) error stop 'singlet MRSF packed transformation failed'
  gram_s = matmul(transpose(singlet),singlet)
  if (maxval(abs(gram_s-identity8())) > 1.0e-14_dp) &
    error stop 'singlet MRSF transformation is not orthonormal'
  if (abs(singlet(2,2)-1.0_dp) > 1.0e-14_dp .or. &
      maxval(abs(singlet(6,:))) > 1.0e-14_dp) &
    error stop 'singlet packed L amplitude was folded twice'

  call mrsf_physical_dimensions(4,3,1,3,expanded,physical,status)
  if (status /= 0 .or. expanded /= 9 .or. physical /= 6) &
    error stop 'triplet MRSF dimensions are incorrect'
  call build_mrsf_packed_transform(4,3,1,3,triplet,status)
  if (status /= 0) error stop 'triplet MRSF packed transformation failed'
  gram_t = matmul(transpose(triplet),triplet)
  if (maxval(abs(gram_t-identity6())) > 1.0e-14_dp) &
    error stop 'triplet MRSF transformation is not orthonormal'
  if (abs(triplet(2,2)-1.0_dp) > 1.0e-14_dp .or. &
      maxval(abs(triplet(6,:))) > 1.0e-14_dp) &
    error stop 'triplet packed L amplitude was folded twice'
  if (maxval(abs(triplet([3,5],:))) > 1.0e-14_dp) &
    error stop 'triplet transformation retains G or D'

contains

  pure function identity8() result(matrix)
    real(kind=dp) :: matrix(8,8)
    integer :: i
    matrix=0.0_dp
    do i=1,8
      matrix(i,i)=1.0_dp
    end do
  end function identity8

  pure function identity6() result(matrix)
    real(kind=dp) :: matrix(6,6)
    integer :: i
    matrix=0.0_dp
    do i=1,6
      matrix(i,i)=1.0_dp
    end do
  end function identity6

end program test_mrsf_hessian_space
