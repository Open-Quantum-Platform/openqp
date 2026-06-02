module pcg_mod

  use precision, only: dp
  use iso_c_binding, only: c_ptr, c_loc, c_null_ptr, c_f_pointer
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

!#################################################################

  private
  public PCG_CONVERGED
  public PCG_OK
  public PCG_NOT_INITIALIZED
  public PCG_BAD_ARGUMENT
  public PCG_BREAKDOWN
  public pcg_matvec
  public pcg_t
  public pcg_optimize

!#################################################################

  integer, parameter :: PCG_CONVERGED       = -1
  integer, parameter :: PCG_OK              = 0
  integer, parameter :: PCG_NOT_INITIALIZED = 1
  integer, parameter :: PCG_BAD_ARGUMENT    = 2
  integer, parameter :: PCG_BREAKDOWN       = 3
  real(kind=dp), parameter :: PCG_DENOMINATOR_FLOOR = 1.0d-24

  integer, parameter :: msglen = 32

  character(len=msglen), parameter :: &
    errmsg(-1:*) = [ &
        character(len=msglen) :: &
        "PCG_CONVERGED" &
      , "PCG_OK" &
      , "PCG_NOT_INITIALIZED" &
      , "PCG_BAD_ARGUMENT" &
      , "PCG_BREAKDOWN" &
    ]

  interface
    subroutine pcg_matvec(y, x, dat)
      import
      real(kind=dp) :: x(:)
      real(kind=dp) :: y(:)
      type(c_ptr) :: dat
    end subroutine
  end interface

!> @brief PCG solver for equation Ax=b
  type :: pcg_t
    logical :: initialized = .false.
    integer(kind=8) :: errcode = 0
    real(kind=dp), allocatable :: b(:)
    real(kind=dp), allocatable :: x(:)
    real(kind=dp), allocatable :: Ap(:)
    real(kind=dp), allocatable :: p(:)
    real(kind=dp), allocatable :: r(:)
    real(kind=dp), allocatable :: y(:)
    real(kind=dp) :: error = huge(1.0_dp)
    real(kind=dp) :: tol = 0.0_dp
    procedure(pcg_matvec), nopass, pointer :: precond => null()
    procedure(pcg_matvec), nopass, pointer :: update => null()
    type(c_ptr) :: dat = c_null_ptr

  contains

    procedure :: init  => pcg_init
    procedure :: clean => pcg_clean
    procedure :: step  => pcg_step

  end type

!#################################################################

contains

!#################################################################

  subroutine pcg_init(this, b, update, precond, dat, x0, tol)
    implicit none
    class(pcg_t), intent(inout) :: this
    real(kind=dp), intent(in) :: b(:)
    procedure(pcg_matvec) :: update
    procedure(pcg_matvec) :: precond
    real(kind=dp), optional, intent(in) :: x0(:)
    real(kind=dp), optional, intent(in) :: tol
    type(*), target :: dat

    integer :: veclen

    if (size(b) <= 0) then
      this%errcode = PCG_BAD_ARGUMENT
      return
    end if

    if (present(x0)) then
      if (size(x0) /= size(b)) then
        this%errcode = PCG_BAD_ARGUMENT
        return
      end if
    end if

    if (.not. all(ieee_is_finite(b))) then
      this%errcode = PCG_BREAKDOWN
      return
    end if
    if (present(x0)) then
      if (.not. all(ieee_is_finite(x0))) then
        this%errcode = PCG_BREAKDOWN
        return
      end if
    end if

    veclen = ubound(b,1)

    allocate(this%x(veclen), &
             this%Ap(veclen), &
             this%b(veclen), &
             this%p(veclen), &
             this%r(veclen), &
             this%y(veclen), &
             source=0.0_dp)

    this%precond => precond
    this%update => update


    this%b = b
    if (present(x0)) this%x = x0
    if (present(tol)) this%tol = tol

    this%dat = c_loc(dat)

    call this%update(this%Ap, this%x, this%dat)
    if (any(.not. ieee_is_finite(this%Ap))) then
      this%errcode = PCG_BREAKDOWN
      return
    end if
    this%r(:) = this%b - this%Ap
    if (any(.not. ieee_is_finite(this%r))) then
      this%errcode = PCG_BREAKDOWN
      return
    end if
    call this%precond(this%y, this%r, this%dat)
    if (any(.not. ieee_is_finite(this%y))) then
      this%errcode = PCG_BREAKDOWN
      return
    end if
    this%p(:) = this%y

    this%error = norm2(this%r)
    if (.not. ieee_is_finite(this%error)) then
      this%errcode = PCG_BREAKDOWN
      return
    end if

    this%initialized = .true.

  end subroutine

!#################################################################

  subroutine pcg_clean(this)
    implicit none
    class(pcg_t), intent(inout) :: this

    if ( allocated(this%x )) deallocate(this%x)
    if ( allocated(this%Ap)) deallocate(this%Ap)
    if ( allocated(this%b))  deallocate(this%b)
    if ( allocated(this%p )) deallocate(this%p)
    if ( allocated(this%r )) deallocate(this%r)
    if ( allocated(this%y )) deallocate(this%y)

    nullify(this%precond)
    nullify(this%update)
    this%dat = c_null_ptr

    this%error = huge(1.0_dp)
    this%tol = 0.0_dp

    this%errcode = 0
    this%initialized = .false.

  end subroutine

!#################################################################

  subroutine pcg_step(this)
    implicit none
    class(pcg_t), intent(inout) :: this

    real(kind=dp) :: rz, rz_new, pap, alpha, beta

    if (.not.this%initialized) then
      this%errcode = PCG_NOT_INITIALIZED
      return
    end if

    if (any(.not. ieee_is_finite(this%x)) .or. any(.not. ieee_is_finite(this%b)) .or. &
        any(.not. ieee_is_finite(this%p)) .or. any(.not. ieee_is_finite(this%r)) .or. &
        any(.not. ieee_is_finite(this%y))) then
      this%errcode = PCG_BREAKDOWN
      return
    end if

    associate(x  => this%x, Ap => this%Ap, &
              p  => this%p, r  => this%r, y  => this%y, &
              error => this%error)

      call this%update(Ap, p, this%dat)
      if (any(.not. ieee_is_finite(Ap))) then
        this%errcode = PCG_BREAKDOWN
        return
      end if

      if (any(.not. ieee_is_finite(p))) then
        this%errcode = PCG_BREAKDOWN
        return
      end if

      rz = dot_product(r, y)
      pap = dot_product(p, Ap)
      if (.not. ieee_is_finite(pap) .or. .not. ieee_is_finite(rz)) then
        this%errcode = PCG_BREAKDOWN
        return
      end if
      if (.not. pcg_safe_positive_denominator(pap) .or. &
          .not. pcg_safe_positive_denominator(rz)) then
        this%errcode = PCG_BREAKDOWN
        return
      end if

      alpha = rz / pap
      if (.not. ieee_is_finite(alpha)) then
        this%errcode = PCG_BREAKDOWN
        return
      end if

      x(:) = x(:) + alpha*p(:)
      if (.not. ieee_is_finite(x(1))) then
        this%errcode = PCG_BREAKDOWN
        return
      end if
      if (.not. all(ieee_is_finite(x))) then
        this%errcode = PCG_BREAKDOWN
        return
      end if
      r(:) = r(:) - alpha*Ap(:)
      if (any(.not. ieee_is_finite(r))) then
        this%errcode = PCG_BREAKDOWN
        return
      end if

      error = norm2(r)
      if (.not. ieee_is_finite(error)) then
        this%errcode = PCG_BREAKDOWN
        return
      end if

      if (error<this%tol) then
        this%errcode = PCG_CONVERGED
        return
      end if

      if (.not. pcg_safe_positive_denominator(error)) then
        this%errcode = PCG_BREAKDOWN
        return
      end if

      call this%precond(y, r, this%dat)
      if (any(.not. ieee_is_finite(y))) then
        this%errcode = PCG_BREAKDOWN
        return
      end if
      rz_new = dot_product(r, y)
      if (.not. ieee_is_finite(rz_new)) then
        this%errcode = PCG_BREAKDOWN
        return
      end if
      if (.not. pcg_safe_positive_denominator(rz_new) .or. &
          .not. ieee_is_finite(rz_new) .or. any(.not. ieee_is_finite(y))) then
        this%errcode = PCG_BREAKDOWN
        return
      end if
      beta = rz_new / rz
      if (.not. ieee_is_finite(beta)) then
        this%errcode = PCG_BREAKDOWN
        return
      end if
      p(:) = y(:) + beta*p(:)
      if (any(.not. ieee_is_finite(p))) then
        this%errcode = PCG_BREAKDOWN
        return
      end if

    end associate

  end subroutine

!#################################################################

  logical function pcg_safe_positive_denominator(value)
    implicit none
    real(kind=dp), intent(in) :: value

    pcg_safe_positive_denominator = ieee_is_finite(value) .and. &
      abs(value) >= PCG_DENOMINATOR_FLOOR

  end function pcg_safe_positive_denominator

!#################################################################

  subroutine pcg_optimize(b, update, precond, dat, mxit, x0, tol, err, cgiters)

    use messages, only: show_message, with_abort
    implicit none

    real(kind=dp), intent(inout) :: b(:)
    procedure(pcg_matvec) :: update
    procedure(pcg_matvec) :: precond
    real(kind=dp), optional, intent(in) :: x0(:)
    integer, intent(in) :: mxit
    real(kind=dp), intent(in) :: tol
    type(*), intent(in) :: dat
    real(kind=dp), optional, intent(out) :: err
    real(kind=dp), optional, intent(out) :: cgiters

    type(pcg_t) :: pcg
    integer :: iter

    if (present(cgiters)) cgiters = 0

    call pcg%init(b=b, update=update, precond=precond, dat=dat, x0=x0, tol=tol)
    if (pcg%errcode /= PCG_OK) goto 9999

    do iter = 1, mxit
      if (present(cgiters)) cgiters = iter
      call pcg%step()
      select case (pcg%errcode)
        case (PCG_OK)
          continue
        case (PCG_CONVERGED)
          exit
        case default
          goto 9999
      end select
    end do

    if (pcg%errcode == PCG_CONVERGED .or. pcg%errcode == PCG_OK) then
      b = pcg%x
      if (present(err)) err = pcg%error
    end if

    return
    9999 continue

    call show_message('PCG: an error has occured, ' // &
                      trim(errmsg(pcg%errcode)), WITH_ABORT)

  end subroutine


!#################################################################

end module
