! Preconditioned MINRES for symmetric (possibly indefinite) A with an SPD
! preconditioner M.  Paige & Saunders (1975) Lanczos + Givens formulation,
! same module/API shape as pcg_mod so it can slot into the z-vector selection.
!
! Unlike CG/PCG, MINRES does NOT require A to be positive-definite: it minimizes
! the residual over the Krylov space and stays stable on indefinite systems
! (e.g. near singlet/triplet instabilities) at CG-like 3-term-recurrence cost.
module minres_mod

  use precision, only: dp
  use iso_c_binding, only: c_ptr, c_loc, c_null_ptr
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public MINRES_CONVERGED, MINRES_OK, MINRES_NOT_INITIALIZED
  public MINRES_BAD_ARGUMENT, MINRES_BREAKDOWN
  public minres_matvec, minres_t, minres_optimize

  integer, parameter :: MINRES_CONVERGED       = -1
  integer, parameter :: MINRES_OK              = 0
  integer, parameter :: MINRES_NOT_INITIALIZED = 1
  integer, parameter :: MINRES_BAD_ARGUMENT    = 2
  integer, parameter :: MINRES_BREAKDOWN       = 3
  real(kind=dp), parameter :: MINRES_DENOMINATOR_FLOOR = 1.0d-24

  interface
    subroutine minres_matvec(y, x, dat)
      import
      real(kind=dp) :: x(:)
      real(kind=dp) :: y(:)
      type(c_ptr)   :: dat
    end subroutine
  end interface

  !> @brief MINRES solver state for A x = b (A symmetric, M SPD).
  type :: minres_t
    logical :: initialized = .false.
    integer(kind=8) :: errcode = 0
    ! Krylov / Lanczos work vectors
    real(kind=dp), allocatable :: b(:), x(:)
    real(kind=dp), allocatable :: r1(:), r2(:), y(:), v(:), av(:)
    real(kind=dp), allocatable :: w(:), w1(:), w2(:)
    ! Lanczos / Givens scalars
    real(kind=dp) :: beta1 = 0.0_dp, beta = 0.0_dp, oldb = 0.0_dp
    real(kind=dp) :: dbar = 0.0_dp, epsln = 0.0_dp, phibar = 0.0_dp
    real(kind=dp) :: cs = -1.0_dp, sn = 0.0_dp
    real(kind=dp) :: error = huge(1.0_dp), tol = 0.0_dp
    integer :: iter = 0
    procedure(minres_matvec), nopass, pointer :: precond => null()
    procedure(minres_matvec), nopass, pointer :: update  => null()
    type(c_ptr) :: dat = c_null_ptr
  contains
    procedure :: init  => minres_init
    procedure :: clean => minres_clean
    procedure :: step  => minres_step
  end type

contains

!#################################################################

  logical function safe_denominator(value)
    real(kind=dp), intent(in) :: value
    safe_denominator = ieee_is_finite(value) .and. abs(value) >= MINRES_DENOMINATOR_FLOOR
  end function safe_denominator

!#################################################################

  subroutine minres_init(this, b, update, precond, dat, x0, tol)
    class(minres_t), intent(inout) :: this
    real(kind=dp), intent(in) :: b(:)
    procedure(minres_matvec) :: update, precond
    real(kind=dp), optional, intent(in) :: x0(:)
    real(kind=dp), optional, intent(in) :: tol
    type(*), target :: dat
    integer :: n

    if (size(b) <= 0) then
      this%errcode = MINRES_BAD_ARGUMENT
      return
    end if
    if (present(x0)) then
      if (size(x0) /= size(b)) then
        this%errcode = MINRES_BAD_ARGUMENT
        return
      end if
    end if
    if (.not. all(ieee_is_finite(b))) then
      this%errcode = MINRES_BREAKDOWN
      return
    end if

    n = size(b)
    allocate(this%b(n), this%x(n), this%r1(n), this%r2(n), this%y(n), &
             this%v(n), this%av(n), this%w(n), this%w1(n), this%w2(n), &
             source=0.0_dp)

    this%precond => precond
    this%update  => update
    this%b = b
    if (present(x0)) this%x = x0
    if (present(tol)) this%tol = tol
    this%dat = c_loc(dat)

    ! r1 = b - A x
    call this%update(this%av, this%x, this%dat)
    if (.not. all(ieee_is_finite(this%av))) then
      this%errcode = MINRES_BREAKDOWN
      return
    end if
    this%r1 = this%b - this%av
    this%r2 = this%r1

    ! y = M^-1 r1 ; beta1 = sqrt(r1 . y)  (>0 since M SPD)
    call this%precond(this%y, this%r1, this%dat)
    if (.not. all(ieee_is_finite(this%y))) then
      this%errcode = MINRES_BREAKDOWN
      return
    end if
    this%beta1 = dot_product(this%r1, this%y)
    if (.not. ieee_is_finite(this%beta1) .or. this%beta1 < 0.0_dp) then
      this%errcode = MINRES_BREAKDOWN
      return
    end if
    this%beta1 = sqrt(this%beta1)

    this%beta   = this%beta1
    this%oldb   = 0.0_dp
    this%dbar   = 0.0_dp
    this%epsln  = 0.0_dp
    this%phibar = this%beta1
    this%cs     = -1.0_dp
    this%sn     = 0.0_dp
    this%iter   = 0
    this%w  = 0.0_dp
    this%w1 = 0.0_dp
    this%w2 = 0.0_dp

    this%error = this%beta1
    this%initialized = .true.
    if (this%beta1 <= this%tol) this%errcode = MINRES_CONVERGED
  end subroutine

!#################################################################

  subroutine minres_clean(this)
    class(minres_t), intent(inout) :: this
    if (allocated(this%b))  deallocate(this%b)
    if (allocated(this%x))  deallocate(this%x)
    if (allocated(this%r1)) deallocate(this%r1)
    if (allocated(this%r2)) deallocate(this%r2)
    if (allocated(this%y))  deallocate(this%y)
    if (allocated(this%v))  deallocate(this%v)
    if (allocated(this%av)) deallocate(this%av)
    if (allocated(this%w))  deallocate(this%w)
    if (allocated(this%w1)) deallocate(this%w1)
    if (allocated(this%w2)) deallocate(this%w2)
    nullify(this%precond); nullify(this%update)
    this%dat = c_null_ptr
    this%initialized = .false.
    this%errcode = 0
    this%error = huge(1.0_dp)
    this%tol = 0.0_dp
    this%iter = 0
  end subroutine

!#################################################################

  subroutine minres_step(this)
    class(minres_t), intent(inout) :: this
    real(kind=dp) :: s, alpha, gamma, delta, gbar, oldeps, phi, denom

    if (.not. this%initialized) then
      this%errcode = MINRES_NOT_INITIALIZED
      return
    end if

    associate(x => this%x, r1 => this%r1, r2 => this%r2, y => this%y, &
              v => this%v, av => this%av, w => this%w, w1 => this%w1, w2 => this%w2)

      this%iter = this%iter + 1

      ! --- Lanczos step: generate next vector ------------------------------
      if (.not. safe_denominator(this%beta)) then
        this%errcode = MINRES_BREAKDOWN
        return
      end if
      s = 1.0_dp / this%beta
      v = s * y

      call this%update(av, v, this%dat)        ! av = A v
      y = av
      if (this%iter >= 2) y = y - (this%beta / this%oldb) * r1
      alpha = dot_product(v, y)
      y = y - (alpha / this%beta) * r2
      r1 = r2
      r2 = y
      call this%precond(y, r2, this%dat)       ! y = M^-1 r2

      this%oldb = this%beta
      this%beta = dot_product(r2, y)
      if (.not. ieee_is_finite(this%beta) .or. this%beta < 0.0_dp) then
        this%errcode = MINRES_BREAKDOWN
        return
      end if
      this%beta = sqrt(this%beta)

      ! --- apply previous Givens rotation ----------------------------------
      oldeps     = this%epsln
      delta      = this%cs * this%dbar + this%sn * alpha
      gbar       = this%sn * this%dbar - this%cs * alpha
      this%epsln =                       this%sn * this%beta
      this%dbar  =                      -this%cs * this%beta

      ! --- compute and apply new Givens rotation ---------------------------
      gamma = sqrt(gbar*gbar + this%beta*this%beta)
      if (.not. safe_denominator(gamma)) then
        this%errcode = MINRES_BREAKDOWN
        return
      end if
      this%cs    = gbar / gamma
      this%sn    = this%beta / gamma
      phi        = this%cs * this%phibar
      this%phibar = this%sn * this%phibar

      ! --- update solution -------------------------------------------------
      denom = 1.0_dp / gamma
      w1 = w2
      w2 = w
      w  = (v - oldeps*w1 - delta*w2) * denom
      x  = x + phi * w

      this%error = abs(this%phibar)
      if (.not. ieee_is_finite(this%error)) then
        this%errcode = MINRES_BREAKDOWN
        return
      end if
      if (this%error <= this%tol) then
        if (.not. all(ieee_is_finite(x))) then
          this%errcode = MINRES_BREAKDOWN
          return
        end if
        this%errcode = MINRES_CONVERGED
      end if
    end associate
  end subroutine

!#################################################################

  subroutine minres_optimize(b, update, precond, dat, mxit, x0, tol, err, iters)
    real(kind=dp), intent(inout) :: b(:)
    procedure(minres_matvec) :: update, precond
    type(*), intent(in) :: dat
    integer, intent(in) :: mxit
    real(kind=dp), optional, intent(in) :: x0(:)
    real(kind=dp), intent(in) :: tol
    real(kind=dp), optional, intent(out) :: err
    integer, optional, intent(out) :: iters
    type(minres_t) :: m
    integer :: it

    if (present(iters)) iters = 0
    call m%init(b=b, update=update, precond=precond, dat=dat, x0=x0, tol=tol)
    select case (m%errcode)
    case (MINRES_OK)
      do it = 1, mxit
        if (present(iters)) iters = it
        call m%step()
        if (m%errcode /= MINRES_OK) exit
      end do
    case (MINRES_CONVERGED)
      ! initial guess already solves the system
    end select

    if (m%errcode == MINRES_CONVERGED .or. m%errcode == MINRES_OK) b = m%x
    if (present(err)) err = m%error
    call m%clean()
  end subroutine

end module minres_mod
