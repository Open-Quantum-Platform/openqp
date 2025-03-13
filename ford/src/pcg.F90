module pcg_mod

  use precision, only: dp
  use iso_c_binding, only: c_ptr, c_loc, c_null_ptr, c_f_pointer

  implicit none

!#################################################################

  private
  public PCG_CONVERGED
  public PCG_OK
  public PCG_NOT_INITIALIZED
  public PCG_BAD_ARGUMENT
  public pcg_matvec
  public pcg_t
  public pcg_optimize

!#################################################################

  integer, parameter :: PCG_CONVERGED       = -1
  integer, parameter :: PCG_OK              = 0
  integer, parameter :: PCG_NOT_INITIALIZED = 1
  integer, parameter :: PCG_BAD_ARGUMENT    = 2

  integer, parameter :: msglen = 32

  character(len=msglen), parameter :: &
    errmsg(-1:*) = [ &
        character(len=msglen) :: &
        "PCG_CONVERGED" &
      , "PCG_OK" &
      , "PCG_NOT_INITIALIZED" &
      , "PCG_BAD_ARGUMENT" &
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

    if (present(x0)) then
      if (size(x0) /= size(b)) then
        this%errcode = PCG_BAD_ARGUMENT
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

    this%r(:) = this%b - this%Ap
    call this%precond(this%y, this%r, this%dat)
    this%p(:) = this%y

    this%error = norm2(this%r)

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

    real(kind=dp) :: rz, alpha, beta

    if (.not.this%initialized) then
      this%error = PCG_NOT_INITIALIZED
      return
    end if

    associate(x  => this%x, Ap => this%Ap, &
              p  => this%p, r  => this%r, y  => this%y, &
              error => this%error)

      call this%update(Ap, p, this%dat)

      rz = dot_product(r, y)
      alpha = rz / dot_product(p, Ap)

      x(:) = x(:) + alpha*p(:)
      r(:) = r(:) - alpha*Ap(:)

      error = norm2(r)

      if (error<this%tol) then
        this%errcode = PCG_CONVERGED
        return
      end if

      call this%precond(y, r, this%dat)
      beta = dot_product(r, y) / rz
      p(:) = y(:) + beta*p(:)

    end associate

  end subroutine

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
      err = pcg%error
    end if

    return
    9999 continue

    call show_message('PCG: an error has occured, ' // &
                      trim(errmsg(pcg%errcode)), WITH_ABORT)

  end subroutine


!#################################################################

end module
