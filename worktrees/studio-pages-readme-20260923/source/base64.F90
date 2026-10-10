module base64
  use, intrinsic :: iso_c_binding, only: c_char, c_int64_t, c_long_long, c_ptr, c_loc, c_null_char
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public b64_encode, b64_decode

  character(*), parameter :: BASE64_TABLE = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

  interface
    integer(c_long_long) function base64_encode(src, dst, nbytes) bind(c)
      use, intrinsic :: iso_c_binding, only: c_ptr, c_long_long
      type(c_ptr), value :: src
      type(c_ptr), value :: dst
      integer(c_long_long), value :: nbytes
    end function
    integer(c_long_long) function base64_decode(src, dst) bind(c)
      use, intrinsic :: iso_c_binding, only: c_ptr, c_long_long
      type(c_ptr), value :: src
      type(c_ptr), value :: dst
    end function
    integer(c_size_t) function strlen(str) bind(c)
      use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
      type(c_ptr), value :: str
    end function
  end interface

  interface b64_encode
    module procedure b64_encode_int32, b64_encode_int64, b64_encode_real32, b64_encode_real64, b64_encode_char
  end interface

  interface b64_decode
    module procedure b64_decode_int32, b64_decode_int64, b64_decode_real32, b64_decode_real64, b64_decode_char
  end interface
contains

  function c_to_f_string(str) result(res)
    character(c_char), target :: str(*)
    character(:), allocatable :: res
    integer :: slen, i
    slen = 0
    do
      if (str(slen+1) == c_null_char) exit
      slen = slen+1
    end do

    if (slen == 0) then
      res = ''
    else
      allocate (character(len=slen) :: res)
      do i = 1, slen
        res(i:i) = str(i)
      end do
    end if
  end function

  function f_to_c_string(str) result(res)
    character(*) :: str
    character(c_char), allocatable :: res(:)
    integer :: slen, i
    slen = len(str)+1
    allocate (res(slen))
    do i = 1, slen-1
      res(i) = str(i:i)
    end do
    res(slen) = c_null_char
  end function

  subroutine string_fix_c_length(string)
    character(*), target :: string
    integer :: slen
#ifdef __GFORTRAN__
    logical, parameter :: need_fix = .true.
#else
    logical, parameter :: need_fix = .false.
#endif
    if (need_fix) then
      slen = len(string)
      if (slen < strlen(c_loc(string))) string(slen+1:slen+1) = c_null_char
    end if
  end subroutine

  function b64_encode_int32(src) result(res)
    integer(int32), target :: src(:)
    character(:), allocatable :: res
    integer :: slen, rlen
    integer(c_long_long) :: retlen
    character(c_char), target, allocatable :: cres(:)
    slen = ubound(src, 1)
    rlen = (slen*storage_size(src)/8+2)/3*4
    allocate (cres(rlen+1))
    retlen = base64_encode(c_loc(src), c_loc(cres), &
                           int(slen*storage_size(src)/8, c_long_long))
    cres(rlen+1) = c_null_char
    res = c_to_f_string(cres)
  end function

  function b64_encode_int64(src) result(res)
    integer(int64), target :: src(:)
    character(:), target, allocatable :: res
    integer :: slen, rlen
    integer(c_long_long) :: retlen
    character(c_char), target, allocatable :: cres(:)
    slen = ubound(src, 1)
    rlen = (slen*storage_size(src)/8+2)/3*4
    allocate (cres(rlen+1))
    retlen = base64_encode(c_loc(src), c_loc(cres), &
                           int(slen*storage_size(src)/8, c_long_long))
    cres(rlen+1) = c_null_char
    res = c_to_f_string(cres)
  end function

  function b64_encode_real32(src) result(res)
    real(real32), target :: src(:)
    character(:), target, allocatable :: res
    integer :: slen, rlen
    integer(c_long_long) :: retlen
    character(c_char), target, allocatable :: cres(:)
    slen = ubound(src, 1)
    rlen = (slen*storage_size(src)/8+2)/3*4
    allocate (cres(rlen+1))
    retlen = base64_encode(c_loc(src), c_loc(cres), &
                           int(slen*storage_size(src)/8, c_long_long))
    cres(rlen+1) = c_null_char
    res = c_to_f_string(cres)
  end function

  function b64_encode_real64(src) result(res)
    real(real64), target :: src(:)
    character(:), target, allocatable :: res
    integer :: slen, rlen
    integer(c_long_long) :: retlen
    character(c_char), target, allocatable :: cres(:)
    slen = ubound(src, 1)
    rlen = (slen*storage_size(src)/8+2)/3*4
    allocate (cres(rlen+1))
    retlen = base64_encode(c_loc(src), c_loc(cres), &
                           int(slen*storage_size(src)/8, c_long_long))
    cres(rlen+1) = c_null_char
    res = c_to_f_string(cres)
  end function

  function b64_encode_char(src) result(res)
    character(*), target :: src
    character(:), target, allocatable :: res
    character(c_char), target, allocatable :: csrc(:), cres(:)
    integer :: slen, rlen
    integer(c_long_long) :: retlen
    csrc = f_to_c_string(src)
    slen = len(src)
    rlen = (slen*storage_size('a')/8+2)/3*4
    allocate (cres(rlen+1))
    retlen = base64_encode(c_loc(csrc), c_loc(cres), &
                           int(slen, c_long_long))
    cres(rlen+1) = c_null_char
    res = c_to_f_string(cres)
  end function

  subroutine b64_decode_int32(src, res)
    character(*), target :: src
    integer(int32), target, allocatable, intent(inout) :: res(:)
    integer :: slen, rlen
    integer(kind=c_int64_t) :: nbytes
    character(kind=c_char), allocatable, target :: csrc(:)

    slen = len(src)
    csrc = f_to_c_string(src)
    rlen = ((slen+3)/4)*3/(storage_size(res)/8)
    if (allocated(res)) deallocate (res)
    allocate (res(rlen))
    nbytes = base64_decode(c_loc(csrc), c_loc(res))
    deallocate (csrc)
  end subroutine

  subroutine b64_decode_int64(src, res)
    character(*), target :: src
    integer(int64), target, allocatable, intent(inout) :: res(:)
    integer :: slen, rlen
    integer(kind=c_int64_t) :: nbytes
    character(kind=c_char), allocatable, target :: csrc(:)

    slen = len(src)
    csrc = f_to_c_string(src)
    rlen = ((slen+3)/4)*3/(storage_size(res)/8)
    if (allocated(res)) deallocate (res)
    allocate (res(rlen))
    nbytes = base64_decode(c_loc(csrc), c_loc(res))
    deallocate (csrc)
  end subroutine

  subroutine b64_decode_real32(src, res)
    character(*), target :: src
    real(real32), target, allocatable, intent(inout) :: res(:)
    integer :: slen, rlen
    integer(kind=c_int64_t) :: nbytes
    character(kind=c_char), allocatable, target :: csrc(:)

    slen = len(src)
    csrc = f_to_c_string(src)
    rlen = ((slen+3)/4)*3/(storage_size(res)/8)
    if (allocated(res)) deallocate (res)
    allocate (res(rlen))
    nbytes = base64_decode(c_loc(csrc), c_loc(res))
    deallocate (csrc)
  end subroutine

  subroutine b64_decode_real64(src, res)
    character(*), target :: src
    real(real64), target, allocatable, intent(inout) :: res(:)
    integer :: slen, rlen
    integer(kind=c_int64_t) :: nbytes
    character(kind=c_char), allocatable, target :: csrc(:)

    slen = len(src)
    csrc = f_to_c_string(src)
    rlen = ((slen+3)/4)*3/(storage_size(res)/8)
    if (allocated(res)) deallocate (res)
    allocate (res(rlen))
    nbytes = base64_decode(c_loc(csrc), c_loc(res))
    deallocate (csrc)
  end subroutine

  subroutine b64_decode_char(src, res)
    character(*), target :: src
    character(:), target, allocatable, intent(inout) :: res
    integer :: slen, rlen
    integer(kind=c_int64_t) :: nbytes
    integer :: npad
    character(kind=c_char), allocatable, target :: csrc(:), cres(:)

    slen = len(src)
    csrc = f_to_c_string(src)
    npad = scan(src, BASE64_TABLE, back=.true.)
    if (npad /= 0) npad = slen-npad
    rlen = (slen+3)/4*3/(storage_size('a')/8)-npad
    if (allocated(res)) deallocate (res)
    allocate (cres(rlen+1))
    cres(rlen+1) = c_null_char
    nbytes = base64_decode(c_loc(csrc), c_loc(cres))
    res = c_to_f_string(cres)
    deallocate (csrc, cres)
  end subroutine

end module base64
