! @brief This module contains Boys function code
module boys

  use precision, only: dp, qp
  use boys_lut, only: fgrid, xgrid, rfinc, rmr, tlgm, tmax, rxinc, nord, igrid, irgrd
  implicit none

  private
  public boysf

  public fgrid
  public xgrid
  public rfinc
  public rmr
  public rxinc
  public tlgm
  public tmax
  public nord
  public igrid
  public irgrd

  real(kind=8), parameter :: &
    pi = 4*atan(1.0_dp), &
    sqrtpi = sqrt(pi), &
    halfsqrtpi = 0.5d0*sqrtpi

contains

!    subroutine boysf(n, tt, ftout)
  pure subroutine boysf(n, tt, ft)
!!$omp declare simd(boysf) uniform(n) linear(ref(tt))

    implicit none

    real(kind=8), intent(in) :: tt
    integer, intent(in) :: n

    !real(kind=8), intent(out) :: ftout(0:)

    real(kind=8) :: ftf
    !real(kind=8) :: ft(0:16)
    real(kind=8), intent(out) :: ft(0:*)

    real(kind=8) :: tv, tx, fx, et, t2
    integer :: m, ip, ix, ifxgrd, iftgrd

    if (tt > tmax) then

      ftf = halfsqrtpi/sqrt(tt)
      t2 = tt*2
      do m = 0, n
        ft(m) = tlgm(m)*ftf
        ftf = ftf/t2
      end do

    else

      ifxgrd = igrid(n)
      iftgrd = irgrd(n)

      tv = tt*rfinc(ifxgrd)
      tx = tt*rxinc
      ip = nint(tv)
      ix = nint(tx)

      fx = 0
      et = 0
      do m = nord, 0, -1
        fx = (fx*tv+fgrid(m, ip, ifxgrd))
        et = (et*tx+xgrid(m, ix))
      end do

      ft(iftgrd) = fx

      t2 = tt+tt
      do m = iftgrd, 1, -1
        ft(m-1) = (t2*ft(m)+et)*rmr(m)
      end do

    end if

    !ftout(0:n) = ft(0:n)

  end subroutine boysf

end module boys
