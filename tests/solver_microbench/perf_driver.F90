program driver
  use precision, only: dp
  use iso_c_binding, only: c_loc
  use bench_common
  use pcg_mod
  implicit none

  integer, parameter :: n = 200000        ! large vectors -> per-iter overhead visible
  integer, parameter :: n_repeat = 300    ! many solves -> stable timing
  integer, parameter :: mxit = 1000
  real(dp), parameter :: tol = 1.0e-10_dp

  integer, target :: dummy
  real(dp), allocatable :: b(:), x(:), xtrue(:)
  type(pcg_t) :: s
  integer :: rep, k, iters, last_iters
  integer(8) :: c0, c1, crate
  real(dp) :: t, resid, maxabs

  allocate(b(n), x(n), xtrue(n))

  ! ---- PERFORMANCE ---------------------------------------------------------
  ! Fixed RHS; identical for old and new builds.
  do k = 1, n
    b(k) = 1.0_dp + 0.5_dp*sin(real(k,dp))
  end do

  call system_clock(count_rate=crate)
  call reset_counters()
  call system_clock(c0)
  last_iters = 0
  do rep = 1, n_repeat
    call s%init(b=b, update=matvec, precond=precond, dat=dummy, tol=tol)
    iters = 0
    if (s%errcode == PCG_OK) then
      do k = 1, mxit
        call s%step()
        if (s%errcode == PCG_CONVERGED) exit
        if (s%errcode /= PCG_OK) exit
        iters = iters + 1
      end do
    end if
    x = s%x
    resid = s%error
    call s%clean()
    last_iters = iters
  end do
  call system_clock(c1)
  t = real(c1 - c0, dp) / real(crate, dp)

  write(*,'(A,I0,A,I0,A,I0,A,ES12.5,A,F9.4)') &
    'PERF n=', n, ' repeats=', n_repeat, ' iters/solve=', last_iters, &
    ' final_resid=', resid, ' total_time_s=', t

  ! ---- STABILITY: exact initial solution converges at init ----------------
  inject_mode = 0
  ! Solve once to get a reference solution, then feed it back as x0.
  call solve_once(b, x, iters, resid)
  xtrue = x
  call s%init(b=b, update=matvec, precond=precond, dat=dummy, x0=xtrue, tol=tol)
  write(*,'(A,I0,A,L1)') 'STAB exact_x0: errcode=', int(s%errcode), &
       '  converged_at_init=', (s%errcode == PCG_CONVERGED)
  call s%clean()

  ! ---- STABILITY: NaN injected into operator at iteration 3 ---------------
  inject_mode = 1; inject_at = 3
  call reset_counters()
  call s%init(b=b, update=matvec, precond=precond, dat=dummy, tol=tol)
  iters = 0
  do k = 1, mxit
    call s%step()
    if (s%errcode /= PCG_OK) exit
    iters = iters + 1
  end do
  maxabs = solmax(s%x)
  write(*,'(A,I0,A,I0,A,L1,A,ES10.3)') &
    'STAB nan_in_Ax: errcode=', int(s%errcode), ' breakdown_iter=', iters, &
    '  is_breakdown=', (s%errcode == PCG_BREAKDOWN), '  max|x|=', maxabs
  call s%clean()

  ! ---- STABILITY: +Inf injected into preconditioner at iteration 3 --------
  inject_mode = 2; inject_at = 3
  call reset_counters()
  call s%init(b=b, update=matvec, precond=precond, dat=dummy, tol=tol)
  iters = 0
  do k = 1, mxit
    call s%step()
    if (s%errcode /= PCG_OK) exit
    iters = iters + 1
  end do
  maxabs = solmax(s%x)
  write(*,'(A,I0,A,I0,A,L1,A,ES10.3)') &
    'STAB inf_in_Minv: errcode=', int(s%errcode), ' breakdown_iter=', iters, &
    '  is_breakdown=', (s%errcode == PCG_BREAKDOWN), '  max|x|=', maxabs
  call s%clean()
  inject_mode = 0

contains

  subroutine solve_once(b, x, iters, resid)
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: iters
    real(dp), intent(out) :: resid
    type(pcg_t) :: ss
    integer :: kk
    call ss%init(b=b, update=matvec, precond=precond, dat=dummy, tol=tol)
    iters = 0
    do kk = 1, mxit
      call ss%step()
      if (ss%errcode == PCG_CONVERGED) exit
      if (ss%errcode /= PCG_OK) exit
      iters = iters + 1
    end do
    x = ss%x
    resid = ss%error
    call ss%clean()
  end subroutine solve_once

  real(dp) function solmax(v)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(dp), intent(in) :: v(:)
    integer :: i
    solmax = 0.0_dp
    do i = 1, size(v)
      if (ieee_is_finite(v(i))) then
        solmax = max(solmax, abs(v(i)))
      else
        solmax = huge(1.0_dp)   ! flag: a non-finite entry leaked into the solution
        return
      end if
    end do
  end function solmax

end program driver
