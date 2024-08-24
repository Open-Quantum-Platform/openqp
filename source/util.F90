module util
  use precision, only: dp

  implicit none

  character(len=*), parameter :: module_name = "util"

  private
  public :: measure_time
  public :: e_charge_repulsion

contains
  subroutine measure_time(print_total, log_unit)
    use iso_fortran_env, only: int64, real64
    implicit none
    integer, intent(in) :: print_total, log_unit
    integer(int64), save :: start_clock, previous_clock
    integer(int64) :: clock_rate, current_clock
    real(real64), save :: start_cpu_time, previous_cpu_time
    real(real64) :: current_cpu_time, elapsed_cpu_time, elapsed_wall_time
    logical, save :: first_call = .true.

    if (first_call) then
      call system_clock(start_clock, clock_rate)
      call cpu_time(start_cpu_time)
      previous_clock = start_clock
      previous_cpu_time = start_cpu_time
      first_call = .false.
    end if

    call system_clock(current_clock, clock_rate)
    call cpu_time(current_cpu_time)

    elapsed_wall_time = real(current_clock - previous_clock, &
                             real64) / real(clock_rate, real64)
    elapsed_cpu_time = current_cpu_time - previous_cpu_time

    write(log_unit, "(3X, A, F10.3, 3X, A, F10.3)") &
            "Step  CPU time (seconds): ", elapsed_cpu_time, &
            "Wall time (seconds): ", elapsed_wall_time

    previous_clock = current_clock
    previous_cpu_time = current_cpu_time

    if (print_total /= 0) then
      elapsed_wall_time = real(current_clock-start_clock, &
                               real64) / real(clock_rate, real64)
      elapsed_cpu_time = current_cpu_time - start_cpu_time

      write(log_unit, "(3X, A, F10.3, 3X, A, F10.3)") &
              "Total CPU time (seconds): ", elapsed_cpu_time, &
              "Wall time (seconds): ", elapsed_wall_time
    end if
  end subroutine measure_time

  pure function e_charge_repulsion(xyz, q) result(enuc)
    use precision, only: dp
    implicit none
    real(kind=dp) :: enuc
    real(kind=dp), intent(in) :: xyz(:,:), q(:)
    real(kind=dp), parameter :: disttol = 1.0d-10
    integer :: i, j, nat
    real(kind=dp) :: rr

    nat = min(ubound(q,1), ubound(xyz,2))
    enuc = 0.0d0
    do i = 2, nat
      do j = 1, i-1
        rr = sum((xyz(:,i)-xyz(:,j))**2)
        if (rr >= disttol) enuc = enuc + q(i)*q(j)/sqrt(rr)
      end do
    end do
  end function e_charge_repulsion

end module util
