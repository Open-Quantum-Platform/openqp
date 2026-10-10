! module precision
!>    @author  Vladimir Mironov
!>
!>    @brief   Contains constants for floating
!>             point number precision
!
!     REVISION HISTORY:
!>    @date _May, 2016_ Initial release
!
module precision

  use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64, real128
  implicit none

  public

  integer, parameter :: I1B = int8
  integer, parameter :: I2B = int16
  integer, parameter :: I4B = int32
  integer, parameter :: I8B = int64

  integer, parameter :: SP = real32
  integer, parameter :: DP = real64
!   Quad precision:
  integer, parameter :: QP = real128
!   Default floating pcision:
  integer, parameter :: FP = DP

  integer, parameter :: SPC = real32
  integer, parameter :: DPC = real64
  integer, parameter :: LGC = kind(.true.)

end module precision
