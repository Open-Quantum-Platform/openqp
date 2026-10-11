!*MODULE messages
!> @brief   This module provides routines where the output is
!> @details Mostly, this file is needed for simplifying of usage
!>            the LibXC interface in different software
!>          This file can be expanded for additional OpenQP messages.
!>            For example, aborting with printing custom message
!> @author  Igor S. Gerasimov
!> @date    July, 2021 - Initial release -
!> @todo    Remove it
!> @params  WITH_ABORT    - logical key for stopmode; stop will be
!> @params  WITHOUT_ABORT - logical key for stopmode; stop will not be
module messages
  use precision, only: dp
#ifdef OQP
  use io_constants, only: write_unit => iw
#else
  use comm_IOFILE, only: write_unit => IW
  use comm_PAR, only: master_worker => MASWRK
#endif
  implicit none
  private
  logical, parameter :: WITH_ABORT = .true.
  logical, parameter :: WITHOUT_ABORT = .false.
  interface show_message
    module procedure show_message_text, &
      show_message_with_integer, &
      show_message_with_double, &
      show_message_with_double_and_text, &
      show_message_with_integer_and_text, &
      show_message_with_keys
  end interface show_message
  public show_message, with_abort, without_abort
contains
  !> @brief   Unit that shows a message to the user
  !> @details OpenQP modules open the log on write_unit only while they
  !>          write to it.  A message issued while it is closed went to a
  !>          file fort.6 that gfortran creates in the working directory,
  !>          where nobody looks (the analytic TD Hessian aborted that way).
  !>          Use standard error until the log is open again.
  integer function message_unit() result(unit)
    use, intrinsic :: iso_fortran_env, only: error_unit
    logical :: connected
    inquire (unit=write_unit, opened=connected)
    unit = write_unit
    if (.not. connected) unit = error_unit
  end function message_unit
  !> @brief   Print simple message
  !> @author  Igor S. Gerasimov
  !> @date    July, 2021 - Initial release -
  !> @params  message  (in)           - message for displaying
  !> @params  stopmode (in, optional) - is aborting required?
  subroutine show_message_text(message, stopmode)
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: stopmode
    logical :: stopmode_
    integer :: u
    if (.not. present(stopmode)) then
      stopmode_ = WITHOUT_ABORT
    else
      stopmode_ = stopmode
    end if
    u = message_unit()
#ifndef OQP
    if (master_worker) &
#endif
      write (u, "(A)") message
    call abort(stopmode_)
  end subroutine show_message_text
  !> @brief   Print simple message
  !> @details write( ,format) message, value
  !> @author  Igor S. Gerasimov
  !> @date    July, 2021 - Initial release -
  !> @params  format   (in)           - format for displaying
  !> @params  message  (in)           - message for displaying
  !> @params  value    (in)           - value for displaying
  !> @params  stopmode (in, optional) - is aborting required?
  subroutine show_message_with_integer(format, message, value, stopmode)
    character(len=*), intent(in) :: format
    character(len=*), intent(in) :: message
    integer, intent(in) :: value
    logical, intent(in), optional :: stopmode
    logical :: stopmode_
    integer :: u
    if (.not. present(stopmode)) then
      stopmode_ = WITHOUT_ABORT
    else
      stopmode_ = stopmode
    end if
    u = message_unit()
#ifndef OQP
    if (master_worker) &
#endif
      write (u, format) message, value
    call abort(stopmode_)
  end subroutine show_message_with_integer
  !> @brief   Print simple message
  !> @details write( ,format) message, value
  !> @author  Igor S. Gerasimov
  !> @date    July, 2021 - Initial release -
  !> @params  format   (in)           - format for displaying
  !> @params  message  (in)           - message for displaying
  !> @params  value    (in)           - value for displaying
  !> @params  stopmode (in, optional) - is aborting required?
  subroutine show_message_with_double(format, message, value, stopmode)
    character(len=*), intent(in) :: format
    character(len=*), intent(in) :: message
    real(dp), intent(in) :: value
    logical, intent(in), optional :: stopmode
    logical :: stopmode_
    integer :: u
    if (.not. present(stopmode)) then
      stopmode_ = WITHOUT_ABORT
    else
      stopmode_ = stopmode
    end if
    u = message_unit()
#ifndef OQP
    if (master_worker) &
#endif
      write (u, format) message, value
    call abort(stopmode_)
  end subroutine show_message_with_double
  !> @brief   Print simple message
  !> @details write( ,format) message1, value, message2
  !> @author  Igor S. Gerasimov
  !> @date    July, 2021 - Initial release -
  !> @params  format   (in)           - format for displaying
  !> @params  message1 (in)           - message for displaying
  !> @params  value    (in)           - value for displaying
  !> @params  message2 (in)           - message for displaying
  !> @params  stopmode (in, optional) - is aborting required?
  subroutine show_message_with_double_and_text(format, message1, value, message2, stopmode)
    character(len=*), intent(in) :: format
    character(len=*), intent(in) :: message1, message2
    real(kind=dp), intent(in) :: value
    logical, intent(in), optional :: stopmode
    logical :: stopmode_
    integer :: u
    if (.not. present(stopmode)) then
      stopmode_ = WITHOUT_ABORT
    else
      stopmode_ = stopmode
    end if
    u = message_unit()
#ifndef OQP
    if (master_worker) &
#endif
      write (u, format) message1, value, message2
    call abort(stopmode_)
  end subroutine show_message_with_double_and_text
  !> @brief   Print simple message
  !> @details write( ,format) message1, value, message2
  !> @author  Igor S. Gerasimov
  !> @date    July, 2021 - Initial release -
  !> @params  format   (in)           - format for displaying
  !> @params  message1 (in)           - message for displaying
  !> @params  value    (in)           - value for displaying
  !> @params  message2 (in)           - message for displaying
  !> @params  stopmode (in, optional) - is aborting required?
  subroutine show_message_with_integer_and_text(format, message1, value, message2, stopmode)
    character(len=*), intent(in) :: format
    character(len=*), intent(in) :: message1, message2
    integer, intent(in) :: value
    logical, intent(in), optional :: stopmode
    logical :: stopmode_
    integer :: u
    if (.not. present(stopmode)) then
      stopmode_ = WITHOUT_ABORT
    else
      stopmode_ = stopmode
    end if
    u = message_unit()
#ifndef OQP
    if (master_worker) &
#endif
      write (u, format) message1, value, message2
    call abort(stopmode_)
  end subroutine show_message_with_integer_and_text
  !> @brief   Print simple message
  !> @details write( ,'(A)') message
  !>          write( ,format) keys
  !> @author  Igor S. Gerasimov
  !> @date    July, 2021 - Initial release -
  !> @params  message  (in)           - message for displaying
  !> @params  format   (in)           - format for displaying
  !> @params  keys     (in)           - keys for displaying
  !> @params  stopmode (in, optional) - is aborting required?
  subroutine show_message_with_keys(message, format, keys, stopmode)
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: format
    character(len=*), intent(in), optional :: keys(:)
    logical, intent(in), optional :: stopmode
    logical :: stopmode_
    integer :: u
    if (.not. present(stopmode)) then
      stopmode_ = WITHOUT_ABORT
    else
      stopmode_ = stopmode
    end if
    u = message_unit()
#ifndef OQP
    if (master_worker) then
#endif
      write (u, "(A)") message
      write (u, format) keys
#ifndef OQP
    end if
#endif
    call abort(stopmode_)
  end subroutine show_message_with_keys
  subroutine abort(stopmode)
    logical, intent(in) :: stopmode
    integer :: u
    u = message_unit()
    flush(u)
    if (stopmode) then
#ifdef OQP
#ifdef __GFORTRAN__
      call backtrace
#endif
      error stop "See above"
#else
      call abrt
#endif
    end if
  end subroutine abort
end module messages
