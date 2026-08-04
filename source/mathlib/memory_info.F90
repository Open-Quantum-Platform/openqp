!> @brief What this process may actually allocate, and a guard built on it.
!>
!> @details Modules that allocate one dominant array up front should refuse
!>   before the allocator does, and say what was needed against what was
!>   there.  A ceiling compiled into the source cannot do that -- the same
!>   binary runs on a laptop and on a 500 GB node, and a fixed cap is wrong
!>   in both directions: it lets a laptop start a job that dies in the
!>   allocator, and it refuses a node that had the memory all along.
!>
!>   `oqp_available_memory_gb` reports the tightest of physical RAM, Linux
!>   MemAvailable, and the cgroup limit -- the last being what actually
!>   binds under SLURM or a container.  `OQP_MEMORY_LIMIT_GB` overrides it.
!>
!>   It returns a negative value when nothing could be determined.  Callers
!>   must treat that as "unknown" and let the job proceed: refusing every
!>   calculation because the probe failed would be worse than the failure
!>   mode it protects against.
module memory_info

  use precision, only: dp
  use, intrinsic :: iso_c_binding, only: c_int64_t

  implicit none

  private
  public :: oqp_available_memory_gb
  public :: oqp_memory_check
  public :: oqp_mem_str

  !> Fraction of the probed memory a single estimate may claim.  The
  !> estimates cover a module's own dominant arrays, not the SCF data
  !> already resident, the integral engine's buffers, or the allocator's
  !> fragmentation, so leaving headroom is what makes the guard useful
  !> rather than merely optimistic.
  real(dp), parameter, public :: OQP_MEMORY_SAFETY_FRACTION = 0.9_dp

  interface
    function oqp_available_memory_bytes() bind(c, name="oqp_available_memory_bytes") &
        result(n)
      import :: c_int64_t
      integer(c_int64_t) :: n
    end function
  end interface

contains

!> @brief Format a size in GB for humans, picking the unit from the magnitude.
!>
!> A fixed "F10.2 GB" renders every small figure as `0.00 GB`, which turns
!> "needed X, had Y" into "needed 0.00, had 0.00".  The guard's whole value is
!> in those two numbers, so they have to survive the printing.
  function oqp_mem_str(gb) result(s)
    real(dp), intent(in) :: gb
    character(len=24) :: s
    real(dp) :: mb

    if (gb >= 1024.0_dp) then
      write(s,'(F0.2,A)') gb/1024.0_dp, ' TB'
    else if (gb >= 1.0_dp) then
      write(s,'(F0.2,A)') gb, ' GB'
    else
      mb = gb*1024.0_dp
      if (mb >= 1.0_dp) then
        write(s,'(F0.1,A)') mb, ' MB'
      else
        write(s,'(F0.1,A)') mb*1024.0_dp, ' kB'
      end if
    end if
    s = adjustl(s)
  end function oqp_mem_str

!> @brief Memory this process may use, in GB; negative if it cannot be probed.
  function oqp_available_memory_gb() result(gb)
    real(dp) :: gb
    integer(c_int64_t) :: bytes

    bytes = oqp_available_memory_bytes()
    if (bytes <= 0_c_int64_t) then
      gb = -1.0_dp
    else
      gb = real(bytes, dp)/1.073741824e9_dp
    end if
  end function oqp_available_memory_gb

!> @brief Refuse a calculation whose estimated peak does not fit this machine.
!>
!> @param[in] need_gb  estimated peak for the caller's dominant arrays
!> @param[in] what     short description used in the message, e.g. "CCSD(T)"
!> @param[in] advice   what the user can change, e.g. "freeze more core orbitals"
!> @param[in] iw       unit to report on
!>
!> Reports the estimate and the probe on @p iw either way, so a run that
!> proceeds still records what it expected to need.  Aborts only when the
!> probe succeeded and the estimate does not fit.
  subroutine oqp_memory_check(need_gb, what, advice, iw)

    use messages, only: show_message, with_abort

    real(dp), intent(in) :: need_gb
    character(len=*), intent(in) :: what, advice
    integer, intent(in) :: iw

    real(dp) :: avail_gb, budget_gb
    character(len=512) :: msg

    avail_gb = oqp_available_memory_gb()

    if (avail_gb < 0.0_dp) then
      write(iw,'(2X,A)') trim(what)//': estimated peak memory ~'// &
          trim(oqp_mem_str(need_gb))
      write(iw,'(2X,A)') trim(what)//': cannot determine available memory on ' // &
          'this machine; proceeding without a check.'
      write(iw,'(2X,A)') trim(what)//': set OQP_MEMORY_LIMIT_GB to impose a budget.'
      return
    end if

    budget_gb = OQP_MEMORY_SAFETY_FRACTION*avail_gb

    write(iw,'(2X,A)') trim(what)//': estimated peak memory ~'// &
        trim(oqp_mem_str(need_gb))//' of '//trim(oqp_mem_str(avail_gb))//' available'

    if (need_gb > budget_gb) then
      msg = trim(what)//' needs about '//trim(oqp_mem_str(need_gb))// &
            ' but only '//trim(oqp_mem_str(avail_gb))// &
            ' is available on this machine (usable budget '// &
            trim(oqp_mem_str(budget_gb))//').'
      write(iw,'(/2X,A)') trim(msg)
      write(iw,'(2X,A)') 'Reduce the cost -- '//trim(advice)//' -- or run ' // &
          'on a machine with more memory.'
      write(iw,'(2X,A)') 'If this machine really has more memory than reported ' // &
          '(some batch systems do not expose it), set OQP_MEMORY_LIMIT_GB.'
      ! No close(iw) before the abort: show_message writes through this same
      ! unit, and closing it first would discard the one-line summary in the
      ! exact low-memory scenario the message exists for.
      call show_message(trim(msg)//' Reduce the cost ('//trim(advice)// &
          '), use a larger machine, or set OQP_MEMORY_LIMIT_GB if the ' // &
          'probe under-reports.', with_abort)
    end if

  end subroutine oqp_memory_check

end module memory_info
