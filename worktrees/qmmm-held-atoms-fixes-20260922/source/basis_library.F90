module basis_library

  use, intrinsic :: iso_fortran_env, only: real64
  use elements, only: MAX_ELEMENT_Z, get_element_id
  use strings, only: to_upper
  use constants, only: ANGULAR_LABEL, NUM_CART_BF
  use io_constants, only: IW

  implicit none

  private

  type, public :: atom_basis_t
    integer :: nshells = 0
    integer :: nbfs = 0
    integer :: nprims = 0
    integer, allocatable :: ang(:)
    integer, allocatable :: ncontract(:)
    real(real64), allocatable :: ex(:)
    real(real64), allocatable :: cc(:)
  contains
    procedure, private :: reserve => atom_basis_reserve
  end type

  type, public :: basis_library_t
    type(atom_basis_t) :: atoms(MAX_ELEMENT_Z)
  contains
    procedure :: from_file => read_basis_library
    procedure :: calc_req_storage
    procedure :: echo => basis_library_echo
  end type

contains

!-------------------------------------------------------------------------------

  subroutine calc_req_storage(this, z, nshell, ngauss, nbasis)
    class(basis_library_t), target, intent(in) :: this
    integer, intent(in) :: z(:)
    integer, intent(out) :: nshell, ngauss, nbasis

    nshell = sum(this%atoms(z)%nshells)
    nbasis = sum(this%atoms(z)%nbfs)
    ngauss = sum(this%atoms(z)%nprims)

  end subroutine

!-------------------------------------------------------------------------------

  subroutine read_basis_library(this, fname)
    class(basis_library_t), target, intent(inout) :: this
    character(*), intent(in) :: fname
    integer :: iunit
    integer :: stat
    character(len=1024) :: msg

    open (newunit=iunit, file=trim(adjustl(fname)), status='old', action='read', &
          iostat=stat, iomsg=msg)

    if (stat /= 0) then
      write (IW, *) 'Can''t open basis library file ''', trim(adjustl(fname)), ''''
      write (IW, *) 'The error is:'
      write (IW, *) trim(msg)
      flush (IW)
      call abort
    end if

    call read_basis_library_file(this, iunit)

    close (iunit)

  end subroutine

!-------------------------------------------------------------------------------

  subroutine read_basis_library_file(this, iunit)
    class(basis_library_t), target, intent(inout) :: this
    integer, intent(in) :: iunit
    type(atom_basis_t), pointer :: atom
    integer :: z
    integer :: ig
    integer :: ln
    integer :: state
    character(1) :: label
    integer :: ng
    character(len=1024) :: error
    integer :: iost

    integer :: id

    character(1024) :: line

    integer :: l
    integer, parameter :: &
      READ_ATOM = 0, &
      READ_SHELL = 1, &
      READ_GAUSS = 2

    state = READ_ATOM
    ig = 0
    ln = 0

    do
      read (unit=iunit, fmt='(A)', end=100) line
      ln = ln+1
!          Check if line is commented out
      if (skip_string(line)) cycle

      select case (state)
      case (READ_ATOM)
        z = get_element_id(trim(adjustl(line)))

        if (z > 0) then
          state = READ_SHELL
          atom => this%atoms(z)
        end if

      case (READ_SHELL)
        if (len_trim(line) == 0) then
          state = READ_ATOM
          ig = 0
          cycle
        end if

        read (line, *, err=200, iostat=iost, iomsg=error) label, ng

        atom%nshells = atom%nshells+1
        call atom%reserve(nshell=atom%nshells, ngauss=ig+ng)

        l = scan(ANGULAR_LABEL, label)
        atom%ang(atom%nshells) = l-1
        atom%nbfs = atom%nbfs+NUM_CART_BF(l-1)

        atom%ncontract(atom%nshells) = ng
        atom%nprims = atom%nprims+ng

        if (atom%ang(atom%nshells) < 0) then
          write (IW, *) "Unknown basis function type: '", label, "'"
          goto 200
        end if

        state = READ_GAUSS

      case (READ_GAUSS)
        ig = ig+1
        read (line, fmt=*, err=200, iostat=iost, iomsg=error) id, atom%ex(ig), atom%cc(ig)

        ng = ng-1
        if (ng == 0) state = READ_SHELL

      end select
    end do

100 continue
    return

200 continue
    write (IW, '("Error while reading basis set file at line ",I0,":")') ln
    write (IW, '(A)') trim(line)
    call abort

  end subroutine

!-------------------------------------------------------------------------------

  subroutine basis_library_echo(this)
    use elements, only: MAX_ELEMENT_Z, ELEMENTS_LONG_NAME
    class(basis_library_t), intent(in) :: this
    integer :: ia, sh, ig, g

    do ia = 1, MAX_ELEMENT_Z
      associate (at => this%atoms(ia))
        if (at%nshells == 0) cycle
        write (*, *) trim(ELEMENTS_LONG_NAME(ia)), at%nshells, at%nprims, at%nbfs
        ig = 0
        do sh = 1, at%nshells
          write (*, '(A1,i10)') ANGULAR_LABEL(at%ang(sh):at%ang(sh)), at%ncontract(sh)
          do g = 1, at%ncontract(sh)
            ig = ig+1
            write (*, '(i4,2ES25.15)') g, at%ex(ig), at%cc(ig)
          end do
        end do
      end associate
    end do
  end subroutine

!-------------------------------------------------------------------------------

  subroutine atom_basis_reserve(this, nshell, ngauss)
    class(atom_basis_t), intent(inout) :: this
    integer, optional, intent(in) :: nshell, ngauss

    if (present(nshell)) then
      call reserve_array_int(this%ang, nshell)
      call reserve_array_int(this%ncontract, nshell)
    end if

    if (present(ngauss)) then
      call reserve_array_real(this%ex, ngauss)
      call reserve_array_real(this%cc, ngauss)
    end if
  end subroutine

  subroutine reserve_array_int(a, n)
    integer, intent(in) :: n
    integer, allocatable, intent(inout) :: a(:)
    integer, allocatable :: tmp(:)
    integer :: lda
    lda = 0
    if (allocated(a)) lda = ubound(a, 1)
    if (n <= lda) return
    allocate (tmp(max(n, lda+16)))
    if (allocated(a)) tmp(:) = a(:)
    call move_alloc(from=tmp, to=a)
  end subroutine

  subroutine reserve_array_real(a, n)
    integer, intent(in) :: n
    real(real64), allocatable, intent(inout) :: a(:)
    real(real64), allocatable :: tmp(:)
    integer :: lda
    lda = 0
    if (allocated(a)) lda = ubound(a, 1)
    if (n <= lda) return
    allocate (tmp(max(n, lda+16)))
    if (allocated(a)) tmp(:lda) = a(:)
    call move_alloc(from=tmp, to=a)
  end subroutine

  logical function skip_string(str) result(res)
    character(*) :: str
    character(*), parameter :: TO_SKIP = '!#&$/\'
    integer :: pos
    res = .false.
    pos = verify(str, ' ')
    if (pos == 0) return
    res = scan(TO_SKIP, str(pos:pos)) /= 0
  end function

end module
