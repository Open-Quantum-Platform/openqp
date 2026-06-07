! Synthetic SPD operator + Jacobi preconditioner for the PCG micro-benchmark,
! plus controllable NaN/Inf injection to exercise the breakdown guards.
module bench_common
  use precision, only: dp
  use iso_c_binding, only: c_ptr
  implicit none

  ! Operator:  A x (i) = diag*x(i) - offd*( x(i-1) + x(i+1) )   (Dirichlet ends)
  ! Symmetric, positive-definite, diagonally dominant for offd < diag/2.
  real(dp) :: a_diag = 2.0_dp
  real(dp) :: a_offd = 0.49_dp        ! cond ~ (diag+2offd)/(diag-2offd)
  real(dp) :: op_shift = 0.0_dp       ! A <- A - op_shift*I  (>0 -> indefinite)
  logical  :: precond_identity = .false.  ! .true. keeps M = I (SPD) for indefinite A
  logical  :: indef_diag = .false.    ! .true. -> A = diag(+1,-1,+1,-1,...) (indefinite)

  ! Injection control (for stability tests):
  !   inject_mode = 0  none
  !               = 1  put NaN into A*x output at/after inject_at-th matvec call
  !               = 2  put +Inf into the preconditioner output at/after inject_at
  integer :: inject_mode = 0
  integer :: inject_at   = 1
  integer :: matvec_calls = 0
  integer :: precond_calls = 0

contains

  subroutine reset_counters()
    matvec_calls = 0
    precond_calls = 0
  end subroutine reset_counters

  subroutine matvec(y, x, dat)
    real(dp) :: x(:)
    real(dp) :: y(:)
    type(c_ptr) :: dat
    integer :: i, n
    n = size(x)
    matvec_calls = matvec_calls + 1
    if (indef_diag) then
      do i = 1, n
        y(i) = merge(1.0_dp, -1.0_dp, mod(i,2) == 1) * x(i)
      end do
      if (inject_mode == 1 .and. matvec_calls >= inject_at) y(1) = ieee_nan()
      return
    end if
    do i = 1, n
      y(i) = (a_diag - op_shift)*x(i)
      if (i > 1) y(i) = y(i) - a_offd*x(i-1)
      if (i < n) y(i) = y(i) - a_offd*x(i+1)
    end do
    if (inject_mode == 1 .and. matvec_calls >= inject_at) then
      y(1) = ieee_nan()
    end if
  end subroutine matvec

  subroutine precond(y, x, dat)
    real(dp) :: x(:)
    real(dp) :: y(:)
    type(c_ptr) :: dat
    precond_calls = precond_calls + 1
    if (precond_identity) then
      y = x                             ! M = I (SPD) for the indefinite case
    else
      y = x / a_diag                    ! Jacobi (SPD): M^-1 = diag(1/a_diag)
    end if
    if (inject_mode == 2 .and. precond_calls >= inject_at) then
      y(1) = ieee_posinf()
    end if
  end subroutine precond

  real(dp) function ieee_nan()
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    ieee_nan = ieee_value(1.0_dp, ieee_quiet_nan)
  end function ieee_nan

  real(dp) function ieee_posinf()
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    ieee_posinf = ieee_value(1.0_dp, ieee_positive_inf)
  end function ieee_posinf

end module bench_common
