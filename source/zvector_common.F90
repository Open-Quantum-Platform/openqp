!> Shared helpers for the TDHF/SF/MRSF z-vector (CPHF/CPKS) solvers.
!>
!> These routines were previously duplicated, one copy per response module
!> (`tdhf_z_vector`, `tdhf_sf_z_vector`, `tdhf_mrsf_z_vector`). They are
!> collected here so the three z-vector drivers share a single, tested
!> implementation.
module zvector_common

  use precision, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public :: sanitize_zvector_preconditioner

contains

!> @brief Build a finite diagonal (Jacobi) preconditioner `xminv = 1/xm`.
!>
!> Each denominator that is non-finite or smaller in magnitude than `floor`
!> is replaced by `+/-floor` (sign preserved), so the preconditioner can never
!> introduce a NaN/Inf or an overflow. The number of regularized entries is
!> reported on `log_unit`, tagged with `tag` (e.g. "RHF", "SF", "MRSF").
!>
!> `floor` is supplied by the caller because the response modules use slightly
!> different thresholds (1e-12 for RHF/SF, 1e-14 for MRSF).
  subroutine sanitize_zvector_preconditioner(xm, xminv, log_unit, floor, tag)
    real(kind=dp),    intent(in)  :: xm(:)
    real(kind=dp),    intent(out) :: xminv(:)
    integer,          intent(in)  :: log_unit
    real(kind=dp),    intent(in)  :: floor
    character(len=*), intent(in), optional :: tag

    integer :: i, regularized
    real(kind=dp) :: denom
    character(len=16) :: prefix

    prefix = ''
    if (present(tag)) prefix = trim(tag)//' '

    regularized = 0
    do i = 1, size(xm)
      denom = xm(i)
      if (.not. ieee_is_finite(denom) .or. abs(denom) < floor) then
        if (ieee_is_finite(denom) .and. denom < 0.0_dp) then
          denom = -floor
        else
          denom =  floor
        end if
        regularized = regularized + 1
      end if
      xminv(i) = 1.0_dp / denom
    end do

    if (regularized > 0) then
      write(log_unit,'(1x,A,"z-vector preconditioner regularized ",I0," denominator(s)")') &
            trim(prefix), regularized
      call flush(log_unit)
    end if

  end subroutine sanitize_zvector_preconditioner

end module zvector_common
