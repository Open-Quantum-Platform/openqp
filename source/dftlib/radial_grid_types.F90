module dft_radial_grid_types
  use precision, only: dp
  implicit none

  private
  public :: get_radial_grid
  public :: dft_radial_grid_mhl
  public :: dft_radial_grid_mk3
  public :: dft_radial_grid_ta
  public :: dft_radial_grid_becke
  public :: dft_radial_grid_map_uniform
  public :: dft_radial_grid_map_cheb2

!> No grid
  integer, parameter :: dft_radial_grid_none = -1

!> Murray-Handy-Laming grid
  integer, parameter :: dft_radial_grid_mhl = 0

!> Mura-Knowles Log3 grid
  integer, parameter :: dft_radial_grid_mk3 = 1

!> Treutler-Ahlrichs grid
  integer, parameter :: dft_radial_grid_ta = 2

!> Becke grid
  integer, parameter :: dft_radial_grid_becke = 3

!> Use default mapping
  integer, parameter :: dft_radial_grid_map_default = 0

!> Map the uniform grid
  integer, parameter :: dft_radial_grid_map_uniform = 1

!> Map the Chebyshev 2nd kind grid
  integer, parameter :: dft_radial_grid_map_cheb2 = 2

!> @brief Base class for radial grids calculation
  type, abstract :: rad_grid_transform_t
    real(kind=dp) :: interval(2)
    integer :: map_grid = 0
  contains
    procedure, nopass :: empty => rad_grid_transform_t_empty
    procedure(set_radial_grid), deferred :: set
    procedure, nopass  :: query => base_query
    procedure(int_1d_transformation), deferred :: transform
  end type

  abstract interface
    subroutine set_radial_grid(this, map_grid, opt)
      import
      class(rad_grid_transform_t) :: this
      integer, optional :: map_grid
      real(kind=dp), optional :: opt
    end subroutine

    function int_1d_transformation(this, r) result(rdr)
      import
      class(rad_grid_transform_t) :: this
      real(kind=dp), intent(inout) :: r !< root
      real(kind=dp) :: rdr(2) !< transformed root and the derivative
    end function
  end interface

  type, extends(rad_grid_transform_t) :: mhl_grid
    real(kind=dp) :: mr
  contains
    procedure, pass :: set => mhl_set
    procedure, nopass :: query => mhl_query
    procedure, pass :: transform => mhl_transform
  end type

  type, extends(rad_grid_transform_t) :: mk3_grid
    real(kind=dp) :: alpha0
  contains
    procedure :: set => mk3_set
    procedure, nopass :: query => mk3_query
    procedure :: transform => mk3_transform
  end type

  type, extends(rad_grid_transform_t) :: ta_grid
    real(kind=dp) :: pwr
  contains
    procedure :: set => ta_set
    procedure, nopass :: query => ta_query
    procedure :: transform => ta_transform
  end type

  type, extends(rad_grid_transform_t) :: becke_grid
  contains
    procedure :: set => becke_set
    procedure, nopass :: query => becke_query
    procedure :: transform => becke_transform
  end type

contains

!> @brief Set the radial points and weights
!> @param ptrad [inout]  quadrature points
!> @param wtrad [inout]  quadrature weights
!> @param rad_grid_type [in]  type of radial grid
!> @param rad_grid_type [in]  which quadrature to map on the selected grid
  subroutine get_radial_grid(ptrad, wtrad, nrad, rad_grid_type, map_grid)
    real(kind=dp), intent(inout) :: ptrad(nrad), wtrad(nrad)
    integer, intent(in) :: nrad, rad_grid_type
    integer, intent(in), optional :: map_grid

    real(kind=dp) :: rdr(2)
    integer :: i

    class(rad_grid_transform_t), allocatable :: rad_grid_transform

    select case (rad_grid_type)
    case (dft_radial_grid_mhl)
      allocate (mhl_grid :: rad_grid_transform)
    case (dft_radial_grid_mk3)
      allocate (mk3_grid :: rad_grid_transform)
    case (dft_radial_grid_ta)
      allocate (ta_grid :: rad_grid_transform)
    case (dft_radial_grid_becke)
      allocate (becke_grid :: rad_grid_transform)
!   Unknown grid type
    case default
      write (*, *) 'unknown radial grid type =', rad_grid_type
      call abort
    end select

    call rad_grid_transform%set(map_grid)

    select case (rad_grid_transform%map_grid)
    case (dft_radial_grid_map_uniform)
      call getuniform(nrad, ptrad, wtrad, rad_grid_transform%interval)
    case (dft_radial_grid_map_cheb2)
      call getcheb2(nrad, ptrad, wtrad, rad_grid_transform%interval)
!   Unknown map_grid type
    case default
      write (*, *) 'unknown map grid type =', rad_grid_transform%map_grid
      call abort
    end select

!   Apply variable transformation and scale weights of the original
!   quadrature
    do i = 1, nrad
      rdr = rad_grid_transform%transform(ptrad(i))
      ptrad(i) = rdr(1)
      wtrad(i) = rdr(2)*wtrad(i)
    end do

!   Include spherical r**2 factor
    wtrad = wtrad*ptrad**2

  end subroutine

!----------------------------------------------------------------------

!> @brief  Murray-Handy-Laming grid variable transformation.
!> @detail  This grid is based on the following variable transformation:
!>   ri = R0 * (r / (1-r) )**m_r
!>   to map interval (0, +1) to (0, +inf).
!>   R0 is a scaling coefficient a.k.a. Bragg-Slater radius
  function mhl_transform(this, r) result(rdr)
    class(mhl_grid) :: this
    real(kind=dp), intent(inout) :: r !< root
    real(kind=dp) :: rdr(2)
    rdr(1) = (r/(1.0d0-r))**this%mr ! r
    rdr(2) = this%mr*r**(this%mr-1)/(1.0d0-r)**(this%mr+1) ! dr
  end function

!> @brief  Murray-Handy-Laming grid init
  subroutine mhl_set(this, map_grid, opt)
    class(mhl_grid) :: this
    integer, optional :: map_grid
    real(kind=dp), optional :: opt
!   Murray-Handy-Laming grid parameter `m_r`:
    integer, parameter :: m_r = 2

    this%mr = m_r
    if (present(opt)) this%mr = opt

!   Default uniform grid
    this%map_grid = dft_radial_grid_map_uniform
    if (present(map_grid)) then
        if (map_grid/=dft_radial_grid_map_default) this%map_grid = map_grid
    end if

    this%interval = [0.0d0, 1.0d0]
  end subroutine

  subroutine mhl_query(grid_type, map_default)
    integer, optional :: grid_type
    integer, optional :: map_default
    grid_type = dft_radial_grid_mhl
    map_default = dft_radial_grid_map_uniform
  end subroutine

!----------------------------------------------------------------------

!> @brief  Mura-Knowles Log-3 grid
!> @detail  This grid is based on the following variable transformation:
!>   r_i = -alpha*log(1-(x_i)**3),
!>   Alpha is a scaling factor
!>   In the original paper authors suggest alpha=7.0 for 1st and 2nd groups
!>   of periodic table and alpha=5.0 otherwise. Here different value
!>   is used:
!>   alpha=alpha0*Rbs,
!>   where Rbs is Bragg-Slater radius. Similar approach is used in NWChem.
  function mk3_transform(this, r) result(rdr)
    class(mk3_grid) :: this
    real(kind=dp), intent(inout) :: r !< root
    real(kind=dp) :: rdr(2)
    real(kind=dp) :: y, dy

    y = 1.0d0-r**3
    dy = 3.0d0*r**2
    rdr(1) = -this%alpha0*log(y) ! r
    rdr(2) = this%alpha0*dy/y ! dr
  end function

!> @brief  Mura-Knowles Log-3 grid init
  subroutine mk3_set(this, map_grid, opt)
    class(mk3_grid) :: this
    integer, optional :: map_grid
    real(kind=dp), optional :: opt

    this%alpha0 = 3.95d0
    if (present(opt)) this%alpha0 = opt

!   Default uniform grid
    this%map_grid = dft_radial_grid_map_uniform
    if (present(map_grid)) then
        if (map_grid/=dft_radial_grid_map_default) this%map_grid = map_grid
    end if

    this%interval = [0.0d0, 1.0d0]
  end subroutine

  subroutine mk3_query(grid_type, map_default)
    integer, optional :: grid_type
    integer, optional :: map_default
    grid_type = dft_radial_grid_mhl
    map_default = dft_radial_grid_map_uniform
  end subroutine


!----------------------------------------------------------------------

!> @brief  Treutler-Ahlrichs grid variable transformation.
!> @detail  This grid is based on the following variable transformation:
!>   ri = R0/log(2) * (1+(xi))**ta_pow * log(2/(1-xi))
!>   to map interval (-1, +1) to (0, +inf).
!>   In the original paper, Chebyshev 2nd kind grid is
!>   used. R0 are per-atom scaling coefficients.
  function ta_transform(this, r) result(rdr)
    class(ta_grid) :: this
    real(kind=dp), intent(inout) :: r !< root
    real(kind=dp) :: rdr(2)

    real(kind=dp), parameter :: log2m1 = 1.0d0/log(2.0d0)
    real(kind=dp) :: tpow, tlog

    tpow = (1.0d0+r)**this%pwr
    tlog = log(2.0d0/(1.0d0-r))
    rdr(1) = log2m1*tpow*tlog ! r
    rdr(2) = log2m1*(this%pwr*tpow*tlog/(1.0d0+r)+tpow/(1.0d0-r)) ! dr
  end function

!> @brief  Treutler-Ahlrichs grid init
  subroutine ta_set(this, map_grid, opt)
    class(ta_grid) :: this
    integer, optional :: map_grid
    real(kind=dp), optional :: opt
!   Treutler-Ahlrichs grid parameter:
    real(kind=dp), parameter :: ta_pow = 0.6d0

    this%pwr = ta_pow
    if (present(opt)) this%pwr = opt

!   Default Chebyshev grid
    this%map_grid = dft_radial_grid_map_cheb2
    if (present(map_grid)) then
        if (map_grid/=dft_radial_grid_map_default) this%map_grid = map_grid
    end if

    this%interval = [-1.0d0, 1.0d0]
  end subroutine

  subroutine ta_query(grid_type, map_default)
    integer, optional :: grid_type
    integer, optional :: map_default
    grid_type = dft_radial_grid_mhl
    map_default = dft_radial_grid_map_cheb2
  end subroutine

!----------------------------------------------------------------------

!> @brief  Becke grid variable transformation.
!> @detail  This grid is based on the following variable transformation:
!>   ri = R0*(1+xi)/(1-xi)
!>   to map interval (-1, +1) to (0, +inf).
!>   In the original paper, Chebyshev 2nd kind grid is
!>   used. R0 are per-atom scaling coefficients.
  function becke_transform(this, r) result(rdr)
    class(becke_grid) :: this
    real(kind=dp), intent(inout) :: r !< root
    real(kind=dp) :: rdr(2)

    call this%empty !dummy, avoid warnings

    rdr(1) = (1.0d0+r)/(1.0d0-r) ! r
    rdr(2) = 2.0d0/(1.0d0-r)**2 ! dr
  end function

!> @brief  Becke grid init
  subroutine becke_set(this, map_grid, opt)
    class(becke_grid) :: this
    integer, optional :: map_grid
    real(kind=dp), optional :: opt

    if (present(opt)) continue

!   Default Chebyshev grid
    this%map_grid = dft_radial_grid_map_cheb2
    if (present(map_grid)) then
        if (map_grid/=dft_radial_grid_map_default) this%map_grid = map_grid
    end if

    this%interval = [-1.0d0, 1.0d0]
  end subroutine

  subroutine becke_query(grid_type, map_default)
    integer, optional :: grid_type
    integer, optional :: map_default
    grid_type = dft_radial_grid_mhl
    map_default = dft_radial_grid_map_cheb2
  end subroutine

!----------------------------------------------------------------------
! Helper subroutines
!----------------------------------------------------------------------

!> @brief Get Chebyshev 2nd kind grid for numerical integration
!> @param n [in] order
!> @param r [out] roots
!> @param w [out] weights
!> @param d [in] optional initial interval
  subroutine getcheb2(n, r, w, d)
    integer :: n !< order
    real(kind=dp), intent(out) :: r(:) !< roots
    real(kind=dp), intent(out) :: w(:) !< weights
    real(kind=dp), intent(in), optional :: d(2) !< initial interval
    real(kind=dp), parameter :: d0(2) = [-1.0d0, 1.0d0]
    real(kind=dp), parameter :: pi = 3.141592653589793238d0

    integer :: i

    r = [(cos(pi*(n-i+1.0d0)/(n+1)), i=1, n)] ! to ensure ascending root order
    w = [(pi/(n+1)*sin(pi*(n-i+1.0d0)/(n+1)), i=1, n)]

    if (present(d)) call linear_map(r, w, d0, d)

  end subroutine

!----------------------------------------------------------------------

!> @brief Get uniform grid for numerical integration
!> @param n [in] order
!> @param r [out] roots
!> @param w [out] weights
!> @param d [in] optional initial interval
  subroutine getuniform(n, r, w, d)
    integer :: n !< order
    real(kind=dp), intent(out) :: r(:) !< roots
    real(kind=dp), intent(out) :: w(:) !< weights
    real(kind=dp), intent(in), optional :: d(2) !< initial interval
    real(kind=dp), parameter :: d0(2) = [0.0d0, 1.0d0]

    integer :: i

    r = [(i*1.0d0/(n+1), i=1, n)]
    w = [(1.0d0/(n+1), i=1, n)]

    if (present(d)) call linear_map(r, w, d0, d)

  end subroutine

!----------------------------------------------------------------------

!> @brief Linear map of roots and weights to a new interval
!> @param r [out] roots
!> @param w [out] weights
!> @param d0 [in] interval to map from
!> @param d [in] interval to map to
  subroutine linear_map(r, w, d0, d)
    real(kind=dp), intent(inout) :: r(:) !< roots
    real(kind=dp), intent(inout) :: w(:) !< weights
    real(kind=dp), intent(in) :: d0(2)   !< initial interval
    real(kind=dp), intent(in) :: d(2)    !< new interval

    r = (d(2)-d(1))/(d0(2)-d0(1))*(r-d0(1))+d(1)
    w = (d(2)-d(1))/(d0(2)-d0(1))*w
  end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------

  subroutine rad_grid_transform_t_empty
  end subroutine

  subroutine base_query(grid_type, map_default)
    integer, optional :: grid_type
    integer, optional :: map_default
    grid_type = dft_radial_grid_none
    map_default = dft_radial_grid_map_default
  end subroutine
end module
