!> @brief Where self-test and opt-in diagnostic routines write their reports.
!> @details A fixed path such as /tmp/hess1_selftest.out is shared by every
!>   user and every concurrent test run on a host: a second user cannot
!>   replace a file the first one created, and two suites overwrite each
!>   other's results.  Reports go to $OQP_SELFTEST_DIR when it is set (the
!>   pytest drivers point it at a per-test directory), otherwise to the
!>   current working directory.
module selftest_output
  implicit none
  private
  public :: selftest_output_path

contains

  !> @brief Path of the report file @p name.
  function selftest_output_path(name) result(path)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: path
    character(len=:), allocatable :: dir
    integer :: length, status

    call get_environment_variable('OQP_SELFTEST_DIR', length=length, status=status)
    if (status /= 0 .or. length == 0) then
      path = name
      return
    end if
    allocate(character(len=length) :: dir)
    call get_environment_variable('OQP_SELFTEST_DIR', dir)
    path = dir//'/'//name
  end function selftest_output_path

end module selftest_output
