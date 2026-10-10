module logger
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  private
  type, public :: logger_t
    integer :: log_unit = error_unit
    character(:), allocatable :: log_file_name
  contains
    procedure :: log_open
    procedure :: log_close
    procedure :: timestamp
    final :: finalize
  end type
contains
  subroutine timestamp(this, message)
    class(logger_t) :: this
    character(*), intent(in) :: message
    character(8) :: date
    character(10) :: time
    character(5) :: zone
    call date_and_time(date, time, zone)
    write (this%log_unit, '("[ ",2(a,"-"),a,x,2(a,":"),a," ]  ",a)') &
      date(1:4), date(5:6), date(7:8), &
      time(1:2), time(3:4), time(5:6), message
!        write(this%log_unit,'(x,a,2x,"[ ",2(a,"-"),a,x,2(a,":"),a,"(",a,") ]")') &
!                message, date(1:4), date(5:6), date(7:8), &
!                time(1:2), time(3:4), time(5:), zone
  end subroutine

  function log_open(this, fname) result(res)
    class(logger_t) :: this
    character(*), intent(in) :: fname
    integer :: res
    integer :: iunit
    character(:), allocatable :: trim_fname
    character(7) :: is_ro

    res = this%log_close()
    trim_fname = trim(adjustl(fname))

    inquire (file=trim_fname, number=iunit, read=is_ro)
    if (iunit /= -1) then
      open (file=trim_fname, newunit=this%log_unit, iostat=res)
      call move_alloc(trim_fname, this%log_file_name)
    else
      if (is_ro /= 'YES') then
        this%log_unit = iunit
        call move_alloc(trim_fname, this%log_file_name)
      else
        write (error_unit, *) "File: '", trim_fname, "', is already opened for -reading-"
        write (error_unit, *) "Close this file prior to opening it again"
        write (error_unit, *) "Log unit unchanged"
        res = 1
      end if
    end if
  end function

  function log_close(this) result(res)
    class(logger_t) :: this
    integer :: res
    res = 0
    if (this%log_unit /= error_unit) then
      close (this%log_unit, iostat=res)
      if (res == 0) deallocate (this%log_file_name)
      this%log_unit = error_unit
    end if
  end function

  subroutine finalize(this)
    type(logger_t) :: this
    integer :: res
    res = this%log_close()
    if (res /= 0) then
      write (error_unit, *) "Cannot close file: '", this%log_file_name, "'"
      error stop "See above"
    end if
  end subroutine

end module
