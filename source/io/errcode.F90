module errcode
  implicit none

  private
  public assignment(=), operator(==)
  public errcode_t

  integer, parameter, private :: &
    OQP_TOPIC_MASK = int(z'FF00'), &
    OQP_VALUE_MASK = int(z'00FF')

  type :: errcode_t
    integer :: topic
    integer :: code
  contains
    procedure :: getcode => errcode_t_getcode
    procedure :: explain => errcode_t_explain

    procedure, private :: errcode_t_set
  end type

  interface errcode_t
    module procedure errcode_t_init
  end interface

  interface assignment(=)
    module procedure errcode_t_set, errcode_t_set_i
  end interface

  interface operator(==)
    module procedure errcode_t_compare_ee, &
      errcode_t_compare_ie, errcode_t_compare_ei
  end interface
contains

!###############################################################
! error codes
!###############################################################

  function errcode_t_init(val) result(res)
    integer, intent(in) :: val
    type(errcode_t) :: res
    res = val
  end function

  subroutine errcode_t_explain(this)
    class(errcode_t), intent(in) :: this
    write (*, '(2(A,I0.4))') &
      "ErrMsg topic=", this%topic, &
      ", code=", this%code
  end subroutine

  subroutine errcode_t_set(this, val)
    class(errcode_t), intent(inout) :: this
    integer, intent(in) :: val
    this%topic = ishft(val, -8)
    this%code = iand(val, OQP_VALUE_MASK)
  end subroutine

  function errcode_t_getcode(this) result(res)
    class(errcode_t), intent(in) :: this
    integer :: res
    res = ishft(this%topic, 8)+this%code
  end function

  subroutine errcode_t_set_i(val, this)
    integer, intent(out) :: val
    class(errcode_t), intent(in) :: this
    val = this%getcode()
  end subroutine

  function errcode_t_compare_ee(this, another) result(res)
    class(errcode_t), intent(in) :: this
    class(errcode_t), intent(in) :: another
    logical :: res
    res = this%getcode() == another%getcode()
  end function

  function errcode_t_compare_ei(this, another) result(res)
    class(errcode_t), intent(in) :: this
    integer, intent(in) :: another
    logical :: res
    res = this%getcode() == another
  end function

  function errcode_t_compare_ie(this, another) result(res)
    integer, intent(in) :: this
    class(errcode_t), intent(in) :: another
    logical :: res
    res = this == another%getcode()
  end function

end module
