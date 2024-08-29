!> @brief Module to store data of DFT atomic quadratures
!> @author Vladimir Mironov
module mod_grid_storage

  use precision, only: fp
  implicit none

!*******************************************************************************

!> @brief Type to store 3D quadrature grid
!> @author Vladimir Mironov
  type :: grid_3d_t
    !< X coordinates
    real(KIND=fp), allocatable :: x(:)
    !< Y coordinates
    real(KIND=fp), allocatable :: y(:)
    !< Z coordinates
    real(KIND=fp), allocatable :: z(:)
    !< quadrature weights
    real(KIND=fp), allocatable :: w(:)
    !< number of points
    integer :: nPts = 0
    !< grid array index, used for reverse search
    integer :: idGrid = 0
    !< tesselation data: (number or points per zone, tesselation degree)
    integer(KIND=2), allocatable :: izones(:, :)
  contains
    procedure, pass :: set => setGrid
    procedure, pass :: get => getGrid
    procedure, pass :: check => checkGrid
    procedure, pass :: clear => clearGrid
  end type

!*******************************************************************************

!> @brief Pointer to 3D grid container (to be used in arrays)
!> @author Vladimir Mironov
  type :: grid_3d_pt
    type(grid_3d_t), pointer :: p
  end type

!*******************************************************************************

!> @brief Basic array list type to store 3d grids
!> @note DO not use any method except `p => list%get`
!>  in a performance-critical loop!
!> @author Vladimir Mironov
  type :: list_grid_3d_t
    private
    !< Current number of grids stored
    integer, public :: nGrids = 0
    !< Maximum number of grids
    integer :: maxGrids = 0
    !< Grid data
    type(grid_3d_pt), allocatable :: elem(:)
  contains
    procedure, non_overridable :: get_pts => get_grid_pts
    procedure, non_overridable :: set_pts => set_grid_pts
    procedure, non_overridable :: findID => findIDListGrid
    procedure, non_overridable :: getByID => getByIDListGrid
    procedure, non_overridable :: push => pushListGrid
    procedure, non_overridable :: pop => popListGrid

    procedure :: get => getListGrid
    procedure :: set => setListGrid
    procedure :: init => initListGrid
    procedure :: clear => clearListGrid
    procedure :: delete => deleteListGrid
    procedure, private :: extend => extendListGrid
  end type

  type :: atomic_grid_t
    type(list_grid_3d_t), pointer :: spherical_grids
    integer, allocatable       :: sph_npts(:)
    integer                    :: idAtm
    real(kind=fp)              :: rAtm
    real(KIND=fp), allocatable :: sph_radii(:)
    real(KIND=fp), allocatable :: rad_pts(:)
    real(KIND=fp), allocatable :: rad_wts(:)
  end type

  integer, parameter :: DEFAULT_GRID_CHUNK = 32

!*******************************************************************************
  private

  public grid_3d_t
  public grid_3d_pt
  public list_grid_3d_t
  public atomic_grid_t

contains

!*******************************************************************************

!-------------------------------------------------------------------------------

!> @brief Get angular grid data
!> @param[in]   npts  number of points
!> @param[out]  x        x coordinates
!> @param[out]  y        y coordinates
!> @param[out]  z        z coordinates
!> @param[out]  w        weights
!> @param[out]  izones   tesselation data
!> @param[out]  found    whether the grid was found in the dataset
!> @author Vladimir Mironov
  subroutine get_grid_pts(grids, npts, x, y, z, w, izones, found)
    class(list_grid_3d_t) :: grids
    integer, intent(OUT) :: izones(:, :)
    real(KIND=fp), intent(OUT) :: x(:), y(:), z(:), w(:)
    integer, intent(IN) :: npts
    logical, intent(OUT) :: found
    type(grid_3d_t), pointer :: tmp

    tmp => grids%get(npts)
    found = associated(tmp)

    if (found) then
      x = tmp%x
      y = tmp%y
      z = tmp%z
      w = tmp%w
      izones = tmp%izones
    end if

  end subroutine

!-------------------------------------------------------------------------------

!> @brief Set angular grid data
!> @param[in]  npts  number of points
!> @param[in]  x        x coordinates
!> @param[in]  y        y coordinates
!> @param[in]  z        z coordinates
!> @param[in]  w        weights
!> @param[in]  izones   tesselation data
!> @author Vladimir Mironov
  subroutine set_grid_pts(grids, npts, x, y, z, w, izones)
    class(list_grid_3d_t) :: grids
    real(KIND=fp), intent(IN) :: x(:), y(:), z(:), w(:)
    integer, intent(IN) :: izones(:, :)
    integer, intent(IN) :: npts
    type(grid_3d_t), pointer :: tmp
    !TYPE(grid_3d_t) :: grid
    !CALL grid%set(npts,x,y,z,w,izones)
    !tmp => grids%set(grid)
    tmp => grids%set(grid_3d_t(x=x, y=y, z=z, w=w, &
                                    npts=npts, izones=int(izones,2)))
  end subroutine

!*******************************************************************************

!> @brief Get a pointer to a grid stored in a list
!> @param[in]  npts  number of points
!> @author Vladimir Mironov
  function getListGrid(list, npts) result(res)
    class(list_grid_3d_t), target, intent(IN) :: list
    type(grid_3d_t), pointer :: res

    integer, intent(IN) :: npts
    integer :: i
    logical :: found

    found = .false.

    do i = 1, list%nGrids
      res => list%elem(i)%p
      found = res%check(npts)
      if (found) exit
    end do

    if (.not. found) res => null()

  end function getListGrid

!-------------------------------------------------------------------------------

!> @brief Get index of a grid stored in a list
!> @param[in]  npts  number of points
!> @author Vladimir Mironov
  function findIDListGrid(list, npts) result(id)
    class(list_grid_3d_t), target, intent(IN) :: list
    integer :: id

    integer, intent(IN) :: npts
    logical :: found

    found = .false.

    do id = 1, list%nGrids
      found = list%elem(id)%p%check(npts)
      if (found) exit
    end do

    if (.not. found) id = 0

  end function findIDListGrid

!-------------------------------------------------------------------------------

!> @brief Get a pointer to a grid stored in a list by its index
!> @param[in]  id    grid index
!> @author Vladimir Mironov
  function getByIDListGrid(list, id) result(res)
    class(list_grid_3d_t), target, intent(IN) :: list
    type(grid_3d_t), pointer :: res
    integer, intent(IN) :: id

    res => null()
    if (id >= 1 .and. id <= list%nGrids) res => list%elem(id)%p

  end function getByIDListGrid

!-------------------------------------------------------------------------------

!> @brief Set a grid in the list. Return a pointer to the list entry.
!> @param[in]  grid  grid data
!> @author Vladimir Mironov
  function setListGrid(list, grid) result(ptr)
    class(list_grid_3d_t), intent(INOUT) :: list
    type(grid_3d_t), intent(IN) :: grid
    type(grid_3d_t), pointer :: ptr

    ptr => list%get(grid%npts)

    if (associated(ptr)) then
      ptr = grid
    else
      ptr => list%push(grid)
    end if

  end function setListGrid

!-------------------------------------------------------------------------------

!> @brief Add a new grid to the list. Return a pointer to the list entry.
!> @param[in]  grid  grid data
!> @author Vladimir Mironov
  function pushListGrid(list, grid) result(ptr)
    class(list_grid_3d_t), intent(INOUT) :: list
    type(grid_3d_t), intent(IN) :: grid
    type(grid_3d_t), pointer :: ptr

    if (list%nGrids == list%maxGrids) call list%extend

    list%nGrids = list%nGrids+1
    allocate (list%elem(list%nGrids)%p)!, source=grid)
    list%elem(list%nGrids)%p = grid

    ptr => list%elem(list%nGrids)%p

    ptr%idGrid = list%nGrids

  end function pushListGrid

!-------------------------------------------------------------------------------

!> @brief Pop a grid from the top of the list.
!> @author Vladimir Mironov
  function popListGrid(list) result(grid)
    class(list_grid_3d_t), target, intent(INOUT) :: list
    type(grid_3d_t) :: grid

    if (list%nGrids /= 0) then
      grid = list%elem(list%nGrids)%p
      deallocate (list%elem(list%nGrids)%p)
      list%nGrids = list%nGrids-1
    end if

  end function popListGrid

!-------------------------------------------------------------------------------

!> @brief Initialize grid list, possibly specifying its size
!> @param[in]  iSize   desired list size
!> @author Vladimir Mironov
  subroutine initListGrid(list, iSize)
    class(list_grid_3d_t) :: list
    integer, optional, intent(IN) :: iSize
    integer :: isz

    isz = DEFAULT_GRID_CHUNK
    if (present(iSize)) isz = iSize

    if (allocated(list%elem)) then
      call list%clear
      if (list%maxGrids < isz) then
        deallocate (list%elem)
      end if
    end if

    if (.not. allocated(list%elem)) then
      allocate (list%elem(isz))
    end if

    list%nGrids = 0
    list%maxGrids = max(isz, list%maxGrids)

  end subroutine initListGrid

!-------------------------------------------------------------------------------

!> @brief Clear grid list
!> @author Vladimir Mironov
  subroutine clearListGrid(list)
    class(list_grid_3d_t) :: list
    integer :: i
    do i = 1, list%nGrids
      if (associated(list%elem(i)%p)) deallocate (list%elem(i)%p)
    end do
    list%nGrids = 0
  end subroutine clearListGrid

!-------------------------------------------------------------------------------

!> @brief Destroy grid list
!> @author Vladimir Mironov
  subroutine deleteListGrid(list)
    class(list_grid_3d_t) :: list
    call list%clear
    deallocate (list%elem)
    list%maxGrids = 0
  end subroutine deleteListGrid

!-------------------------------------------------------------------------------

!> @brief Extend grid list to a new size
!> @param[in] iSize  size to extend the list
!> @author Vladimir Mironov
  subroutine extendListGrid(list, iSize)
    class(list_grid_3d_t), intent(INOUT) :: list
    type(grid_3d_pt), allocatable :: elem_new(:)

    integer, optional, intent(IN) :: iSize
    integer :: isz

    isz = DEFAULT_GRID_CHUNK
    if (present(iSize)) isz = iSize

    allocate (elem_new(list%maxGrids+isz), source=list%elem)
    call move_alloc(from=elem_new, to=list%elem)
    list%maxGrids = list%maxGrids+isz

  end subroutine

!*******************************************************************************

!> @brief Set 3D grid data
!> @param[in]  npts  number of points
!> @param[in]  x        x coordinates
!> @param[in]  y        y coordinates
!> @param[in]  z        z coordinates
!> @param[in]  w        weights
!> @param[in]  izones   tesselation data
!> @author Vladimir Mironov
  subroutine setGrid(grid, npts, x, y, z, w, izones)
    class(grid_3d_t), intent(INOUT) :: grid
    real(KIND=fp), intent(IN) :: x(:), y(:), z(:), w(:)
    integer, intent(IN) :: izones(:, :)
    integer, intent(IN) :: npts

    grid%npts = npts
    grid%x = x
    grid%y = y
    grid%z = z
    grid%w = w
    grid%izones = int(izones, kind=2)

  end subroutine setGrid

!-------------------------------------------------------------------------------

!> @brief Get 3D grid data
!> @param[in]   npts  number of points
!> @param[out]  x        x coordinates
!> @param[out]  y        y coordinates
!> @param[out]  z        z coordinates
!> @param[out]  w        weights
!> @param[out]  izones   tesselation data
!> @author Vladimir Mironov
  subroutine getGrid(grid, npts, x, y, z, w, izones)
    class(grid_3d_t), intent(IN) :: grid
    real(KIND=fp), intent(OUT) :: x(:), y(:), z(:), w(:)
    integer, intent(OUT) :: izones(:, :)
    integer, intent(OUT) :: npts

    npts = grid%npts
    x = grid%x
    y = grid%y
    z = grid%z
    w = grid%w
    izones = grid%izones

  end subroutine getGrid

!-------------------------------------------------------------------------------

!> @brief Check if it is the grid you are looking for.
!> @param[in]   npts  number of points
!> @author Vladimir Mironov
  function checkGrid(grid, npts) result(ok)
    class(grid_3d_t) :: grid
    integer, intent(IN) :: npts
    logical :: ok
    ok = .false.
    if (grid%npts == npts) ok = .true.
  end function checkGrid

!-------------------------------------------------------------------------------

!> @brief Destroy grid
!> @author Vladimir Mironov
  subroutine clearGrid(grid)
    class(grid_3d_t), intent(INOUT) :: grid
    grid%npts = 0
    deallocate (grid%x, grid%y, grid%z, grid%w)
    deallocate (grid%izones)
  end subroutine clearGrid

!-------------------------------------------------------------------------------

end module mod_grid_storage
