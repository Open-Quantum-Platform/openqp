module atomic_structure_m
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  implicit none
!  Structured Data types for basis set index
  type, public :: atomic_structure
    real(c_double), allocatable :: zn(:)      !< atomic number or nuclear charge
    real(c_double), allocatable :: mass(:)    !< atomic mass
    real(c_double), allocatable :: grad(:,:)  !< Gradient
    real(c_double), allocatable :: xyz(:,:)   !< Atomic coordinates
!      character(len=2) :: Symbol   !< Atomic symbol
!      character(len=8) :: SHTYPS   !< Shell type of basis set
  contains
    procedure, non_overridable :: init => atomic_structure_init
    procedure, non_overridable :: clean => atomic_structure_clean
    procedure, non_overridable :: center => atomic_structure_center
  end type atomic_structure
contains
  function atomic_structure_init(self, natoms) result(ok)
    class(atomic_structure) :: self
    integer :: natoms
    integer(c_int) :: ok
    ok = self%clean()
    if (ok /= 0) return
    allocate( self%zn(natoms) &
            , self%mass(natoms) &
            , self%grad(3,natoms) &
            , self%xyz(3,natoms) &
            , stat=ok)
!      character(len=2) :: Symbol   !< Atomic symbol
  end function atomic_structure_init

  function atomic_structure_clean(self) result(ok)
    class(atomic_structure) :: self
    integer(c_int) :: ok

    ok = 0

    if (allocated(self%zn   ))  deallocate(self%zn   , stat=ok)
    if (ok /= 0) return
    if (allocated(self%mass ))  deallocate(self%mass , stat=ok)
    if (ok /= 0) return
    if (allocated(self%grad ))  deallocate(self%grad , stat=ok)
    if (ok /= 0) return
    if (allocated(self%xyz))  deallocate(self%xyz, stat=ok)
  end function atomic_structure_clean

  function atomic_structure_center(self, weight) result(r)
    use strings, only: to_upper
    implicit none
    class(atomic_structure) :: self
    character(len=*), optional :: weight
    real(c_double) :: r(3)
    character(len=8) :: wtype
    character(len=8), parameter :: WTYPE_NONE = 'NONE'
    character(len=8), parameter :: WTYPE_MASS = 'MASS'

    wtype = WTYPE_NONE

    if (present(weight)) wtype = to_upper(weight)

    r = 0

    select case (wtype)

    case(WTYPE_NONE)
      r = sum(self%xyz,2)/ubound(self%xyz,2)

    case(WTYPE_MASS)
      r = matmul(self%xyz,self%mass)/sum(self%mass)

    case default
      error stop "Unknown weight type in atomic_structure_center: "//wtype
    end select

  end function atomic_structure_center

end module atomic_structure_m
