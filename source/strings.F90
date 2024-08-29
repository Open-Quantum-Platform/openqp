module strings
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t
  implicit none
  private
  character(*), parameter, public :: &
    ALPHANUM = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_' &
    , AN_DASH = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-'
  integer, parameter :: &
    A_LOWER = iachar('a'), &
    Z_LOWER = iachar('z'), &
    A_UPPER = iachar('A'), &
    Z_UPPER = iachar('Z'), &
    UPPER_TO_LOWER = A_LOWER-A_UPPER
  public :: to_upper, to_lower, count_substring, remove_spaces, index_ith, fstring
  public :: f_c_char
  public :: c_f_char
  character, target :: char_empty = ''

  type, public :: tokenizer_t
    character(:), pointer :: str => null()
    character(:), pointer :: delims => null()
  contains
    procedure :: get => next_token
  end type

  type, public, bind(C) :: Cstring
    integer(c_int64_t) :: length
    type(c_ptr) :: string
  end type Cstring

contains

!> @brief  Split string on tokens, delimiting by specified set of symbols
!> @author Vladimir Mironov
!> @date   Sept, 2021
!> @param  this   [inout] `tokenizer_t` instance
!> @param  str    [in]    optional, input string
!> @param  delim  [in]    optional, string containing delimiter symbols, default=''
!> @param  result         pointer to the next token
!> @note   Delimiter instance saves its state after `str` and `delim` has been set up.
!>         At each successfull call user can redefine `delim` string
!>         The behavior is similar to C `strtok` function.
!>         Returns null pointer when line ends and resets `tokenizer_t` internal state.
!>         Example usage:
!>         ```fortran
!>         type(tokenizer_t) :: tokz
!>         character(:), pointer :: ptr
!>         ptr => tokz%get(str, ' ;,')
!>         do while (associated(ptr))
!>           print *, ptr
!>           ptr => tokz%get()
!>         end do
!>         ```
!> @note   Do not use for csv reading!
  function next_token(this, str, delim) result(res)
    class(tokenizer_t) :: this
    character(*), target, intent(in), optional :: str
    character(*), target, intent(in), optional :: delim
    character(:), pointer :: res
    integer :: pos_start, res_len, pos_end
    character(:), pointer :: temp_str

    res => null()

    if (present(str)) then
      nullify (this%str)
      this%str => str
    end if

    if (present(delim)) then
      if (len(delim) /= 0) then
        this%delims => delim
      end if
    end if

    if (len(this%str) == 0) then
      nullify (this%str)
      this%delims => char_empty
      return
    end if

    if (.not. associated(this%str)) return
    if (.not. associated(this%delims)) then
      this%delims => char_empty
    end if

    pos_start = verify(this%str, this%delims)

    if (pos_start == 0) then
      nullify (this%str)
      this%delims => char_empty
      return
    end if

    res_len = scan(this%str(pos_start:), this%delims)-1
    temp_str => this%str(pos_start:)
    if (res_len < 0) res_len = len(temp_str)
    nullify(temp_str)

    pos_end = pos_start+res_len

    res => this%str(pos_start:pos_end-1)

    this%str => this%str(pos_end+1:)

  end function next_token

  !> @brief  return string in upper case
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019 --Initial release--
  !> @date   May, 2021 Moved to strings
  !> @date   Sep, 2021 subroutine -> function
  !> @param  line - (in)
  pure function to_upper(input_string) result(output_string)
    character(*), intent(in) :: input_string
    character(:), allocatable :: output_string
    integer :: i, ic
    output_string = input_string
    do i = 1, len(output_string)
      ic = iachar(output_string(i:i))
      if (ic <= Z_LOWER .and. ic >= A_LOWER) then
        output_string(i:i) = achar(ic-UPPER_TO_LOWER)
      end if
    end do
  end function to_upper
  !> @brief  return string in lower case
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2021 --Initial release--
  !> @param  line - (inout)
  pure function to_lower(input_string) result(output_string)
    character(*), intent(in) :: input_string
    character(:), allocatable :: output_string
    integer :: i, ic
    output_string = input_string
    do i = 1, len(output_string)
      ic = iachar(output_string(i:i))
      if (ic <= Z_UPPER .and. ic >= A_UPPER) then
        output_string(i:i) = achar(ic+UPPER_TO_LOWER)
      end if
    end do
  end function to_lower
  !> @brief  This function return count of substring in string
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019 --Initial release--
  !> @date   May, 2021 Moved to strings
  !> @param  substring - (in)
  !> @param  string    - (in)
  pure integer function count_substring(substring, string) result(res)
    character(len=*), intent(in) :: substring, string
    ! internal variables
    character(len=:), allocatable :: tmp_string
    res = 0
    tmp_string = string
    do
      if (index(tmp_string, substring) == 0) exit
      res = res+1
      tmp_string = tmp_string(index(tmp_string, substring)+1:)
    end do
    return
  end function count_substring
  !> @brief  This routine remove not needed spaces from section line
  !> @detail This routine used revert reading of lines.
  !>         Example, that this routine do:
  !>       > DFTTYP = PBE0 BASNAM=APC4 , ACC5
  !>       < DFTTYP=PBE0 BASNAM=APC4,ACC5
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019 --Initial release--
  !> @date   May, 2021 Moved to strings
  !> @param  line - (inout) worked line
  pure subroutine remove_spaces(line)
    character(len=:), allocatable, intent(inout) :: line
    ! internal variables
    character(len=:), allocatable :: tmp_line, res_line
    integer :: i, ind, ind_end
    logical :: skip, first
    tmp_line = line
    line = ""
    do
      res_line = ""
      ind_end = index(tmp_line, "=")
      first = .true.
      skip = .false.
      ind = ind_end
      if (ind == 0) ind = len(tmp_line)
      do i = ind, 1, -1
        if (tmp_line(i:i) /= " ") then
          res_line = tmp_line(i:i)//res_line
          if (tmp_line(i:i) /= "=" .and. first .and. ind_end /= 0) then
            skip = .true.
            first = .false.
          end if
        else if (skip) then
          res_line = " "//res_line
          skip = .false.
        end if
      end do
      line = line//trim(adjustl(res_line))
      if (ind_end == 0) exit
      tmp_line = tmp_line(ind+1:)
    end do
  end subroutine remove_spaces
  !> @brief   This function return index of i'th substring in string
  !> @details if ind is negative, backward search will be
  !>          result is equal 0 if search was failed
  !> @author  Igor S. Gerasimov
  !> @date    May, 2021 --Initial release--
  !> @param   substring - (in)
  !> @param   string    - (in)
  !> @param   ith       - (in) index of needed substring
  pure integer function index_ith(substring, string, ith) result(res)
    character(len=*), intent(in) :: substring, string
    integer, intent(in) :: ith
    ! internal variables
    character(len=:), allocatable :: tmp_string
    integer :: i, diff
    res = 0
    tmp_string = string
    if (ith > 0) then
      do i = 1, ith
        diff = index(tmp_string, substring)
        res = res+diff
        if (diff == 0) then
          res = 0
          exit
        end if
        tmp_string = tmp_string(index(tmp_string, substring)+1:)
      end do
    else if (ith < 0) then
      do i = -1, ith, -1
        res = index(tmp_string, substring, back=.true.)
        if (res == 0) exit
        tmp_string = tmp_string(:index(tmp_string, substring, back=.true.)-1)
      end do
    else
      res = 0
    end if
  end function index_ith
  !> @brief   This function return fortran allocatable string
  !> @author  Igor S. Gerasimov
  !> @date    April, 2022 --Initial release--
  !> @param   string    - (in) C-like string
  function fstring(string) result(res)
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_char
    type(Cstring), intent(in) :: string
    character(len=:), allocatable :: res
    character(len=1,kind=c_char), pointer :: fpstring(:)
    ! internal variables
    integer :: i
    allocate(character(len=string%length) :: res)
    call c_f_pointer(string%string, fpstring, shape=[string%length])
    do i = 1, string%length
      res(i:i) = fpstring(i)
    end do
  end function fstring

!> @brief Convert Fortran string to a raw null-terminated C-string
!> @note No real boundary checks for c-string are performed, use at your own risk!
!> @note C-string should be already allocated
!> @param[in]   from    fortran string
!> @param[out]  to      null-terminated c-string, len(to) should be at least numchar+1
!> @param[in]   numchar maximum number of characters to pass from fortran to c-string
  subroutine f_c_char(from, to, numchar)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    character(len=*), intent(in) :: from
    character(kind=c_char, len=1) :: to(*)
    integer, intent(in) :: numchar
    integer :: i
    do i = 1, min(len(from), numchar)
      to(i) = from(i:i)
    end do
    to(i) = c_null_char
  end subroutine f_c_char

!> @brief Convert raw null-terminated C-string to Fortran string
!> @note No boundary checks for c-string are performed, use at your own risk!
  function c_f_char(cchar) result(fstr)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    character(kind=c_char), intent(in) :: cchar(*)
    character(:), allocatable :: fstr
    integer :: i, strlen

    strlen = 0
    do
        if (cchar(strlen+1) == c_null_char) exit
        strlen = strlen + 1
    end do

    allocate(character(strlen) :: fstr)

    do i = 1, strlen
        fstr(i:i) = cchar(i)
    end do
  end function c_f_char
end module strings
