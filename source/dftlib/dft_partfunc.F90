module mod_dft_partfunc

  use precision, only: fp
  implicit none

  !> @brief Interface for func: double -> double
  abstract interface
    pure real(KIND=fp) function func_d_d(x)
      import
      real(KIND=fp), intent(IN) :: x
    end function func_d_d
  end interface

!> @brief Type for partition function calculation
  type partition_function
    !< Values beyond (-limit,+limit) interval considered as 0.0 or 1.0
    real(KIND=fp) :: limit = 1.0_fp
    !< Compute partition function value
    procedure(func_d_d), nopass, pointer :: eval
    !< Compute partition function derivative
    procedure(func_d_d), nopass, pointer :: deriv
  contains
    procedure :: set => set_partition_function
  end type

  integer, parameter :: PTYPE_SSF = 0
  integer, parameter :: PTYPE_ERF = 1
  integer, parameter :: PTYPE_BECKE4 = 2
  integer, parameter :: PTYPE_SMSTP2 = 3
  integer, parameter :: PTYPE_SMSTP3 = 4
  integer, parameter :: PTYPE_SMSTP4 = 5
  integer, parameter :: PTYPE_SMSTP5 = 6

  private
  public partition_function
  public PTYPE_BECKE4
  public PTYPE_SSF
  public PTYPE_ERF
  public PTYPE_SMSTP2
  public PTYPE_SMSTP3
  public PTYPE_SMSTP4
  public PTYPE_SMSTP5

!*******************************************************************************
contains

!> @brief Becke's partition function (4th order)
  pure function partf_eval_becke4(x) result(f)
    real(KIND=fp) :: f
    real(KIND=fp), intent(IN) :: x
!    REAL(KIND=fp), PARAMETER :: LIMIT = 1.0, SCALEF = 1.0
    integer :: i
    f = x
    do i = 1, 4
      f = 0.5_fp*f*(3.0_fp-f*f)
    end do
    f = 0.5_fp-0.5_fp*f
  end function

!> @brief Becke's partition function (4th order) derivative
  pure function partf_diff_becke4(x) result(df)
    real(KIND=fp) :: df
    real(KIND=fp), intent(IN) :: x
!    REAL(KIND=fp), PARAMETER :: LIMIT = 1.0, SCALEF = 1.0
    real(KIND=fp), parameter :: FACTOR = 81.0_fp/32.0_fp
    real(KIND=fp) :: f
    integer :: i
    f = x
    df = 1.0_fp
    do i = 1, 4
      df = df*(1.0_fp-f*f)
      f = 0.5_fp*f*(3.0_fp-f*f)
    end do
    df = -FACTOR*df
  end function

!-------------------------------------------------------------------------------

!> @brief SSF partition function
  pure function partf_eval_ssf(x) result(f)
    real(KIND=fp) :: f
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp), parameter :: LIMIT = 0.64_fp, SCALEF = 1.0_fp/LIMIT
    real(KIND=fp) :: f1, f2, f4

    if (abs(x) > LIMIT) then
      f = 0.5_fp-sign(0.5_fp, x)
      return
    end if
    f1 = x*SCALEF
    f2 = f1*f1
    f4 = f2*f2
    f = 0.0625_fp*f1*((35.0_fp-35.0_fp*f2)+f4*(21.0_fp-5.0_fp*f2))
    f = 0.5_fp-0.5_fp*f
  end function

!> @brief SSF partition function derivative
  pure function partf_diff_ssf(x) result(df)
    real(KIND=fp) :: df
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp), parameter :: LIMIT = 0.64_fp, SCALEF = 1.0_fp/LIMIT
    real(KIND=fp) :: f1, f2, f4

    if (abs(x) > LIMIT) then
      df = 0.0_fp
      return
    end if
    f1 = x*SCALEF
    f2 = f1*f1
    f4 = f2*f2
    df = 0.0625_fp*((35.0_fp-105.0_fp*f2)+ &
                    f4*(105.0_fp-35.0_fp*f2))
    df = -0.5_fp*SCALEF*df
  end function

!-------------------------------------------------------------------------------

!> @brief Modified SSF with erf(a*x/(1-x^2)) scaling (like in NWChem)
  pure function partf_eval_erf(x) result(f)
    real(KIND=fp) :: f
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp), parameter :: LIMIT = 0.725_fp, SCALEF = 1.0_fp/0.3_fp

    if (abs(x) > LIMIT) then
      f = 0.5_fp-sign(0.5_fp, x)
      return
    end if
    f = erf(x/(1.0_fp-x**2)*SCALEF)
    f = 0.5_fp-0.5_fp*f
  end function

!> @brief Modified SSF with erf(a*x/(1-x^2)) scaling (like in NWChem) derivative
  pure function partf_diff_erf(x) result(df)
    real(KIND=fp) :: df
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp), parameter :: LIMIT = 0.725_fp, SCALEF = 1.0_fp/0.3_fp
    real(KIND=fp), parameter :: PI = 3.141592653589793d0
    real(KIND=fp), parameter :: FACTOR = SCALEF/sqrt(PI)
    real(KIND=fp) :: ex, frac

    if (abs(x) > LIMIT) then
      df = 0.0_fp
      return
    end if
    frac = 1.0_fp/(1.0_fp-x*x)
    ex = exp(-(SCALEF*SCALEF*x*x*frac*frac))
    df = -FACTOR*ex*(1.0_fp+x*x)*frac
  end function

!-------------------------------------------------------------------------------

!> @brief Smoothstep function (2nd order)
  pure function partf_eval_smoothstep2(x) result(f)
    real(KIND=fp) :: f
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp) :: f1, f2
    real(KIND=fp), parameter :: LIMIT = 0.55_fp, SCALEF = 0.5_fp/LIMIT
    if (abs(x) > LIMIT) then
      f = 0.5_fp-sign(0.5_fp, x)
      return
    else
      f1 = 0.5_fp-x*SCALEF
      f2 = f1*f1
      f = f2*f1*(10.0_fp-15.0_fp*f1+6.0_fp*f2)
    end if
  end function

!> @brief Smoothstep function (2nd order) derivative
  pure function partf_diff_smoothstep2(x) result(df)
    real(KIND=fp) :: df
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp) :: f1, f2
    real(KIND=fp), parameter :: LIMIT = 0.55_fp, SCALEF = 0.5_fp/LIMIT
    if (abs(x) > LIMIT) then
      df = 0.0_fp
      return
    else
      f1 = 0.5_fp-x*SCALEF
      f2 = f1*f1
      df = f2*(30.0_fp-60.0_fp*f1+30.0_fp*f2)
    end if
  end function

!-------------------------------------------------------------------------------

!> @brief Smoothstep function (3th order)
  pure function partf_eval_smoothstep3(x) result(f)
    real(KIND=fp) :: f
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp) :: f1, f2, f4
    real(KIND=fp), parameter :: LIMIT = 0.62_fp, SCALEF = 0.5_fp/LIMIT
    if (abs(x) > LIMIT) then
      f = 0.5_fp-sign(0.5_fp, x)
      return
    else
      f1 = 0.5_fp-x*SCALEF
      f2 = f1*f1
      f4 = f2*f2
      f = f4*((35.0_fp-84.0_fp*f1)+ &
              f2*(70.0_fp-20.0_fp*f1))
    end if
  end function

!> @brief Smoothstep function (3th order) derivative
  pure function partf_diff_smoothstep3(x) result(df)
    real(KIND=fp) :: df
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp) :: f1, f2, f4
    real(KIND=fp), parameter :: LIMIT = 0.62_fp, SCALEF = 0.5_fp/LIMIT
    if (abs(x) > LIMIT) then
      df = 0.0_fp
      return
    else
      f1 = 0.5_fp-x*SCALEF
      f2 = f1*f1
      f4 = f2*f2
      df = f2*f1*((140.0_fp-420.0_fp*f1)+ &
                  f2*(420.0_fp-140.0_fp*f1))
    end if
  end function

!-------------------------------------------------------------------------------

!> @brief Smoothstep function (4th order)
  pure function partf_eval_smoothstep4(x) result(f)
    real(KIND=fp) :: f
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp) :: f1, f2, f4
    real(KIND=fp), parameter :: LIMIT = 0.69_fp, SCALEF = 0.5_fp/LIMIT
    if (abs(x) > LIMIT) then
      f = 0.5_fp-sign(0.5_fp, x)
      return
    else
      f1 = 0.5_fp-x*SCALEF
      f2 = f1*f1
      f4 = f2*f2
      f = f4*f1*((126.0_fp-420.0_fp*f1)+ &
                 f2*(540.0_fp-315.0_fp*f1)+ &
                 f4*70.0_fp)
    end if
  end function

!> @brief Smoothstep function (4th order) derivative
  pure function partf_diff_smoothstep4(x) result(df)
    real(KIND=fp) :: df
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp) :: f1, f2, f4
    real(KIND=fp), parameter :: LIMIT = 0.69_fp, SCALEF = 0.5_fp/LIMIT
    if (abs(x) > LIMIT) then
      df = 0.0_fp
      return
    else
      f1 = 0.5_fp-x*SCALEF
      f2 = f1*f1
      f4 = f2*f2
      df = f4*((630.0_fp-2520.0_fp*f1)+ &
               f2*(3780.0_fp-2520.0_fp*f1)+ &
               f4*630.0_fp)
    end if
  end function

!-------------------------------------------------------------------------------

!> @brief Smoothstep function (5th order)
  pure function partf_eval_smoothstep5(x) result(f)
    real(KIND=fp) :: f
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp) :: f1, f2, f4
    real(KIND=fp), parameter :: LIMIT = 0.73_fp, SCALEF = 0.5_fp/LIMIT
    if (abs(x) > LIMIT) then
      f = 0.5_fp-sign(0.5_fp, x)
      return
    else
      f1 = 0.5_fp-x*SCALEF
      f2 = f1*f1
      f4 = f2*f2
      f = f4*f2*((462.0_fp-1980.0_fp*f1)+ &
                 f2*(3465.0_fp-3080.0_fp*f1)+ &
                 f4*(1386.0_fp-252.0_fp*f1))
    end if
  end function

!> @brief Smoothstep function (5th order) derivative
  pure function partf_diff_smoothstep5(x) result(df)
    real(KIND=fp) :: df
    real(KIND=fp), intent(IN) :: x
    real(KIND=fp) :: f1, f2, f4
    real(KIND=fp), parameter :: LIMIT = 0.73_fp, SCALEF = 0.5_fp/LIMIT
    if (abs(x) > LIMIT) then
      df = 0.0_fp
      return
    else
      f1 = 0.5_fp-x*SCALEF
      f2 = f1*f1
      f4 = f2*f2
      df = f4*f1*((2772.0_fp-13860.0_fp*f1)+ &
                  f2*(27720.0_fp-27720.0_fp*f1)+ &
                  f4*(13860.0_fp-2772.0_fp*f1))
    end if
  end function

!*******************************************************************************

!> @brief Set up partition function parameters
  subroutine set_partition_function(partfunc, ptype)
    class(partition_function), intent(INOUT) :: partfunc
    integer, intent(IN) :: ptype

    select case (ptype)
    case (PTYPE_BECKE4)
      partfunc%limit = 1.0_fp
      partfunc%eval => partf_eval_becke4
      partfunc%deriv => partf_diff_becke4

    case (PTYPE_SSF)
      partfunc%limit = 0.64_fp
      partfunc%eval => partf_eval_ssf
      partfunc%deriv => partf_diff_ssf

    case (PTYPE_ERF)
      partfunc%limit = 0.725_fp
      partfunc%eval => partf_eval_erf
      partfunc%deriv => partf_diff_erf

    case (PTYPE_SMSTP2)
      partfunc%limit = 0.55_fp
      partfunc%eval => partf_eval_smoothstep2
      partfunc%deriv => partf_diff_smoothstep2

    case (PTYPE_SMSTP3)
      partfunc%limit = 0.62_fp
      partfunc%eval => partf_eval_smoothstep3
      partfunc%deriv => partf_diff_smoothstep3

    case (PTYPE_SMSTP4)
      partfunc%limit = 0.69_fp
      partfunc%eval => partf_eval_smoothstep4
      partfunc%deriv => partf_diff_smoothstep4

    case (PTYPE_SMSTP5)
      partfunc%limit = 0.74_fp
      partfunc%eval => partf_eval_smoothstep5
      partfunc%deriv => partf_diff_smoothstep5

    case DEFAULT
      partfunc%limit = 0.62_fp
      partfunc%eval => partf_eval_smoothstep3
      partfunc%deriv => partf_diff_smoothstep3
    end select

  end subroutine

!-------------------------------------------------------------------------------

! SUBROUTINE toupper(s, u)
!     CHARACTER(LEN=*), INTENT(IN)  :: s
!     CHARACTER(LEN=*), INTENT(OUT) :: u
!     INTEGER, PARAMETER :: ISHIFT = (iachar("A")-iachar("a"))
!     INTEGER :: i, code
!
!     DO i = 1, min(len(s), len(u))
!         SELECT CASE (s(i:i))
!         CASE ( "a" : "z" )
!             code = iachar(s(i:i))
!             u(i:i) = achar(code+ISHIFT)
!         CASE DEFAULT
!             u(i:i) = s(i:i)
!         END SELECT
!     END DO
!
! END SUBROUTINE

end module mod_dft_partfunc
