module messages
  implicit none
  integer, parameter :: with_abort = 1
contains
  subroutine show_message(msg, mode)
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: mode
    write(*,'(A)') 'show_message: '//trim(msg)
    if (present(mode)) then
      if (mode == with_abort) stop 1
    end if
  end subroutine show_message
end module messages
