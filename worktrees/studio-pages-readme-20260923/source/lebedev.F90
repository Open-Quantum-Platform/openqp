module lebedev
  use precision, only: dp
  implicit none

  private
  public :: lebedev_get_grid
  public :: spherical_grid_type_lebedev
  public :: spherical_grid_type_o
  public :: spherical_grid_type_oh
  public :: spherical_grid_type_i
  public :: lebedev_npts
  public :: oct_npts
  public :: oh_npts
  public :: i_npts
  public :: lebedev_orders
  public :: oct_orders
  public :: oh_orders
  public :: i_orders

  integer, parameter :: spherical_grid_type_lebedev = 0
  integer, parameter :: spherical_grid_type_o = 1
  integer, parameter :: spherical_grid_type_oh = 2
  integer, parameter :: spherical_grid_type_i = 3

  integer, parameter :: &
    lebedev_npts(*) = [ &
      6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, &
      590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, &
      4802, 5294, 5810]
  integer, parameter :: oct_npts(*) = [246, 264, 342, 432]
  integer, parameter :: oh_npts(*) = [350, 398]
  integer, parameter :: i_npts(*) = [132, 152, 180, 192, 212, 242]

  integer, parameter :: &
    lebedev_orders(*) = [ &
      3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, &
      65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131]
  integer, parameter :: oct_orders(*) = [26, 27, 31, 35]
  integer, parameter :: oh_orders(*) = [31, 33]
  integer, parameter :: i_orders(*) = [19, 20, 21, 23, 24, 25]

contains

  subroutine lebedev_get_grid(npts, xyz, w, grid_type)
    use messages, only: show_message, WITH_ABORT
    integer, intent(in) :: npts
    integer, optional, intent(in) :: grid_type
    real(kind=dp), intent(inout) :: xyz(npts, *), w(*)

    character(:), allocatable :: errmsg
    integer :: n
    integer :: grdtype
    grdtype = 0
    if (present(grid_type)) grdtype = grid_type

    select case (grdtype)
!   Lebedev-Laikov grids
    case (spherical_grid_type_lebedev)
      select case (npts)
      case (6); call ld0006(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (14); call ld0014(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (26); call ld0026(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (38); call ld0038(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (50); call ld0050(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (74); call ld0074(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (86); call ld0086(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (110); call ld0110(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (146); call ld0146(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (170); call ld0170(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (194); call ld0194(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (230); call ld0230(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (266); call ld0266(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (302); call ld0302(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (350); call ld0350(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (434); call ld0434(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (590); call ld0590(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (770); call ld0770(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (974); call ld0974(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (1202); call ld1202(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (1454); call ld1454(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (1730); call ld1730(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (2030); call ld2030(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (2354); call ld2354(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (2702); call ld2702(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (3074); call ld3074(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (3470); call ld3470(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (3890); call ld3890(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (4334); call ld4334(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (4802); call ld4802(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (5294); call ld5294(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (5810); call ld5810(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case default
        write (errmsg, "(' Lebedev grid with ',i4,' points not available.')") npts
        call show_message(errmsg, WITH_ABORT)
      end select

!   Class of grids developed by A.S.Popov.
!   Rotational octahedral symmetry grids (without inversion):
    case (spherical_grid_type_o)
      select case (npts)
      case (246); call od0246(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (264); call od0264(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (342); call od0342(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (432); call od0432(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case default
        write (errmsg, "(' O-symmetry grid with ',i4,' points not available.')") npts
        call show_message(errmsg, WITH_ABORT)
      end select

!   Full octahedral symmetry grids (same symmetry as Lebedev grids):
    case (spherical_grid_type_oh)
      select case (npts)
      case (350); call ohd0350(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (398); call ohd0398(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case default
        write (errmsg, "(' Oh-symmetry grid with ',i4,' points not available.')") npts
        call show_message(errmsg, WITH_ABORT)
      end select

!   Rotational icosahedral symmetry grids:
    case (spherical_grid_type_i)
      select case (npts)
      case (132); call id0132(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (152); call id0152(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (180); call id0180(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (192); call id0192(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (212); call id0212(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case (242); call id0242(xyz(:, 1), xyz(:, 2), xyz(:, 3), w, n)
      case default
        write (errmsg, "(' I-symmetry grid with ',i4,' points not available.')") npts
        call show_message(errmsg, WITH_ABORT)
      end select

    case default
      call show_message("Unknown spherical grid type", WITH_ABORT)

    end select

!   Final check
    if (n /= npts) then
      write (errmsg, '(2(a,i8))') 'Error, the computed number of points N=', n, &
        ' does not match requested NPTS=', npts
      call show_message(errmsg, WITH_ABORT)
    end if
  end subroutine

!> @brief Given a point on a sphere (specified by `a` and `b`), generate all
!>    the equivalent points under octahedral (Oh or O) symmetry,
!>    making grid points with weight `v`.
!> @detail
!>   The variable `num` is increased by the number of different points
!>   generated.
!>
!>   Depending on code, there are 6...48 different but equivalent
!>   points.
!>   Full octahedral symmetry points:
!>   code=1:   (0,0,1) etc                                (  6 points)
!>   code=2:   (0,a,a) etc, a=1/sqrt(2)                   ( 12 points)
!>   code=3:   (a,a,a) etc, a=1/sqrt(3)                   (  8 points)
!>   code=4:   (a,a,b) etc, b=sqrt(1-2 a^2)               ( 24 points)
!>   code=5:   (a,b,0) etc, b=sqrt(1-a^2), a input        ( 24 points)
!>   code=6:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input  ( 48 points)
!>   General octahedral points:
!>   code=7:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input  ( 24 points)

!> @note This routine is part of a set of routines that generate
!>  Lebedev grids [1-6] for integration on a sphere. The original
!>  C-code [1] was kindly provided by Dr. Dmitri N. Laikov and
!>  translated into fortran by Dr. Christoph van Wuellen.
!>  This routine was translated from C to Fortran77 by hand.
!>
!>  The original subroutine was heavily modified and extended
!>  by Dr. Vladimir Mironov.
!>
!>  Users of this code are asked to include reference [1] in their
!>  publications, and in the user- and programmers-manuals
!>  describing their codes.
!>
!>  This code was distributed through CCL (http://www.ccl.net/).
!>
!>  [1] V.I. Lebedev, and D.N. Laikov
!>      "A quadrature formula for the sphere of the 131st
!>       algebraic order of accuracy"
!>      Doklady Mathematics, vol. 59, no. 3, 1999, pp. 477-481.
!>
!>  [2] V.I. Lebedev
!>      "A quadrature formula for the sphere of 59th algebraic
!>       order of accuracy"
!>      russian acad. sci. dokl. math., vol. 50, 1995, pp. 283-286.
!>
!>  [3] V.I. Lebedev, and A.L. Skorokhodov
!>      "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!>      Russian Acad. Sci. Dokl. Math., vol. 45, 1992, pp. 587-592.
!>
!>  [4] V.I. Lebedev
!>      "Spherical quadrature formulas exact to orders 25-29"
!>      Siberian Mathematical Journal, vol. 18, 1977, pp. 99-107.
!>
!>  [5] V.I. Lebedev
!>      "Quadratures on a sphere"
!>      Computational Mathematics and Mathematical Physics, vol. 16,
!>      1976, pp. 10-24.
!>
!>  [6] V.I. Lebedev
!>      "Values of the nodes and weights of ninth to seventeenth
!>       order gauss-markov quadrature formulae invariant under the
!>       octahedron group with inversion"
!>      Computational Mathematics and Mathematical Physics, vol. 15,
!>      1975, pp. 44-51.
!
  subroutine gen_oh(code, num, x, y, z, w, a, b, v)
    use messages, only: show_message, WITH_ABORT
    real(kind=dp), intent(out) :: x(*), y(*), z(*), w(*)
    real(kind=dp), intent(inout) :: a, b, v
    integer, intent(in) :: code
    integer, intent(inout) :: num

    real(kind=dp) :: c

    select case (code)
    case (1)
!     (0,0,1) etc,   6 points
!     Octahedron vertices
      a = 1.0d0
      x(1:6) = 0.0d0
      y(1:6) = 0.0d0
      z(1:6) = 0.0d0

      x(1) =  a
      x(2) = -a

      y(3) =  a
      y(4) = -a

      z(5) =  a
      z(6) = -a

      w(1:6) = v

      num = num+6

    case (2)
!     (0,a,a) etc, a=1/sqrt(2), 12 points
      a = sqrt(0.5d0)
      call gen_td(num, x, y, z, w, 0.0d0, a, a, v)

    case (3)
!     (a,a,a) etc, a=1/sqrt(3), 8 points
!     Cube vertices
      a = sqrt(1.0d0/3.0d0)

      x(1) =  a; y(1) =  a; z(1) =  a
      x(2) = -a; y(2) =  a; z(2) =  a
      x(3) =  a; y(3) = -a; z(3) =  a
      x(4) = -a; y(4) = -a; z(4) =  a
      x(5) =  a; y(5) =  a; z(5) = -a
      x(6) = -a; y(6) =  a; z(6) = -a
      x(7) =  a; y(7) = -a; z(7) = -a
      x(8) = -a; y(8) = -a; z(8) = -a

      w(1:8) = v

      num = num+8

    case (4)
!     (a,a,b) etc, b=sqrt(1-2 a^2), 24 points
      b = sqrt(1.0d0-2.0d0*a*a)
      call gen_td(num, x, y, z, w, a, a, b, v)
      call gen_td(num, x(13), y(13), z(13), w(13), a, a, -b, v)

    case (5)
!     (a,b,0) etc, b=sqrt(1-a^2), a input, 24 points
      b = sqrt(1.0d0-a*a)
      call gen_td(num, x, y, z, w, a, 0.0d0, b, v)
      call gen_td(num, x(13), y(13), z(13), w(13), 0.0d0, a, b, v)

    case (6)
!     (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input, 48 points
      c = sqrt(1.0d0-a*a-b*b)
      call gen_td(num, x, y, z, w, a, b, c, v)
      call gen_td(num, x(13), y(13), z(13), w(13), a, b, -c, v)
      call gen_td(num, x(25), y(25), z(25), w(25), b, a, c, v)
      call gen_td(num, x(37), y(37), z(37), w(37), b, a, -c, v)

    case (7)
!     (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input, 24 points
      c = sqrt(1.0d0-a*a-b*b)

      call gen_td(num, x(1), y(1), z(1), w(1), a, b, c, v)
      call gen_td(num, x(13), y(13), z(13), w(13), a, -c, b, v)

    case default
      call show_message('GEN_OH: INVALID CODE', WITH_ABORT)

    end select

  end subroutine

!---------------------------------------------------------------------

!> @brief Given a point on a sphere (specified by `a` and `b`), generate all
!>    the equivalent points under Ih symmetry, making grid points with
!>    weight `v`.
!> @author Vladimir Mironov
  subroutine gen_ih(code, num, x, y, z, w, a, b, v)
    use messages, only: show_message, WITH_ABORT
    implicit none
    real(kind=dp), intent(inout) :: x(*), y(*), z(*), w(*)
    real(kind=dp), intent(inout) :: a, b, v
    real(kind=dp) :: c
    integer code
    integer num
    integer n
    double precision g, h, t
    parameter(g=(sqrt(5.0d0)+1.0d0)*0.25d0)
    parameter(h=(sqrt(5.0d0)-1.0d0)*0.25d0)
    parameter(t=0.5d0)

    n = 1

    select case (code)
    case (1)
      a = sqrt((5.0d0+sqrt(5.0d0))/10.0d0)
      b = sqrt((5.0d0-sqrt(5.0d0))/10.0d0)
      c = 0.0d0
      call gen_td(n, x(n), y(n), z(n), w(n), a, b, c, v)
      num = num+12

    case (2)
      a = sqrt((3.0d0-sqrt(5.0d0))/6.0d0)
      b = sqrt((3.0d0+sqrt(5.0d0))/6.0d0)
      c = 0.0d0
      call gen_td(n, x(n), y(n), z(n), w(n), a, b, c, v)

      c = 1.0d0/sqrt(3.0d0)
      call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
      num = num+20

    case (3)
      call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
      a = g
      b = h
      c = t
      call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
      num = num+30

    case (4)
      c = sqrt(1.0d0-a*a-b*b)
      call gen_td(n, x(n), y(n), z(n), w(n), a, b, c, v)
      call ih_tr(a, b, c)
      call gen_td(n, x(n), y(n), z(n), w(n), a, b, c, v)
      call ih_tr(a, b, c)
      call gen_td(n, x(n), y(n), z(n), w(n), a, b, c, v)
      call ih_tr(a, b, c)
      call gen_td(n, x(n), y(n), z(n), w(n), a, b, c, v)
      call ih_tr(a, b, c)
      call gen_td(n, x(n), y(n), z(n), w(n), a, b, c, v)
      num = num+60

    case default
      call show_message('GEN_IH: INVALID CODE', WITH_ABORT)

    end select

  end subroutine

!---------------------------------------------------------------------

!> @brief Given a point on a sphere (specified by `x`, `y`, and `z`),
!>    permute its coordinates, so the new point belongs to the next
!>    tetrahedron within icosahedron.
!> @note This is supplementary subroutine to the others
!> @author Vladimir Mironov
  subroutine ih_tr(x, y, z)
    real(kind=dp), intent(inout) :: x, y, z
    real(kind=dp) :: x_tmp, y_tmp, z_tmp
    real(kind=dp), parameter :: g = (sqrt(5.0d0)+1.0d0)*0.25d0
    real(kind=dp), parameter :: h = (sqrt(5.0d0)-1.0d0)*0.25d0
    real(kind=dp), parameter :: t = 0.5d0

! New point:
! |x|   |g  h -t|   |x|
! |y| = |h  t  g| x |y|
! |z|   |t -g  h|   |z|

    x_tmp = g*x+h*y-t*z
    y_tmp = h*x+t*y+g*z
    z_tmp = t*x-g*y+h*z

    x = x_tmp
    y = y_tmp
    z = z_tmp
  end subroutine

!---------------------------------------------------------------------

!> @brief Given a point on a sphere (specified by `a` and `b`), generate all
!>    the equivalent points under Td symmetry, making grid points with
!>    weight `v`.
!> @note This is supplementary subroutine to the others
!> @author Vladimir Mironov
  subroutine gen_td(num, x, y, z, w, a, b, c, v)
    real(kind=dp), intent(inout) :: x(*), y(*), z(*), w(*)
    real(kind=dp), intent(in) :: a, b, c, v
    integer, intent(inout) :: num
!     a*a + b*b + c*c = 1
    x(1)  =  a ;   y(1)  =  b ;   z(1)  =  c
    x(2)  = -a ;   y(2)  = -b ;   z(2)  =  c
    x(3)  = -a ;   y(3)  =  b ;   z(3)  = -c
    x(4)  =  a ;   y(4)  = -b ;   z(4)  = -c
    x(5)  =  c ;   y(5)  =  a ;   z(5)  =  b
    x(6)  =  c ;   y(6)  = -a ;   z(6)  = -b
    x(7)  = -c ;   y(7)  = -a ;   z(7)  =  b
    x(8)  = -c ;   y(8)  =  a ;   z(8)  = -b
    x(9)  =  b ;   y(9)  =  c ;   z(9)  =  a
    x(10) = -b ;   y(10) =  c ;   z(10) = -a
    x(11) =  b ;   y(11) = -c ;   z(11) = -a
    x(12) = -b ;   y(12) = -c ;   z(12) =  a

    w(1:12) = v

    num = num+12
  end subroutine

!---------------------------------------------------------------------

  subroutine ld0006(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 6-POINT ANGULAR GRID
!
!      This routine is part of a set of routines that generate
!      Lebedev grids [1-6] for integration on a sphere. The original
!      C-code [1] was kindly provided by Dr. Dmitri N. Laikov and
!      translated into Fortran by Dr. Christoph van Wüllen.
!      This routine was translated using a C to Fortran77 conversion
!      tool written by Dr. Christoph van Wüllen.
!
!      Users of this code are asked to include reference [1] in their
!      publications, and in the user and programmer manuals
!      describing their codes.
!
!      This code was distributed through CCL (http://www.ccl.net/).
!
!      References:
!
!      [1] V.I. Lebedev and D.N. Laikov,
!          "A Quadrature Formula for the Sphere of the 131st
!          Algebraic Order of Accuracy,"
!          Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
!      [2] V.I. Lebedev,
!          "A Quadrature Formula for the Sphere of 59th Algebraic
!          Order of Accuracy,"
!          Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.
!
!      [3] V.I. Lebedev and A.L. Skorokhodov,
!          "Quadrature Formulas of Orders 41, 47, and 53 for the Sphere,"
!          Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.
!
!      [4] V.I. Lebedev,
!          "Spherical Quadrature Formulas Exact to Orders 25-29,"
!          Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.
!
!      [5] V.I. Lebedev,
!          "Quadratures on a Sphere,"
!          Computational Mathematics and Mathematical Physics, Vol. 16,
!          1976, pp. 10-24.
!
!      [6] V.I. Lebedev,
!          "Values of the Nodes and Weights of Ninth to Seventeenth
!          Order Gauss-Markov Quadrature Formulae Invariant under the
!          Octahedron Group with Inversion,"
!          Computational Mathematics and Mathematical Physics, Vol. 15,
!          1975, pp. 44-51.
!
    n = 1
    v = 0.1666666666666667d+00
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0014(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 14-POINT ANGULAR GRID
!
!      This routine is part of a set of routines that generate
!      Lebedev grids [1-6] for integration on a sphere. The original
!      C-code [1] was kindly provided by Dr. Dmitri N. Laikov and
!      translated into Fortran by Dr. Christoph van Wüllen.
!      This routine was translated using a C to Fortran77 conversion
!      tool written by Dr. Christoph van Wüllen.
!
!      Users of this code are asked to include reference [1] in their
!      publications, and in the user and programmer manuals
!      describing their codes.
!
!      This code was distributed through CCL (http://www.ccl.net/).
!
!      References:
!
!      [1] V.I. Lebedev and D.N. Laikov,
!          "A Quadrature Formula for the Sphere of the 131st
!          Algebraic Order of Accuracy,"
!          Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
!      [2] V.I. Lebedev,
!          "A Quadrature Formula for the Sphere of 59th Algebraic
!          Order of Accuracy,"
!          Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.
!
!      [3] V.I. Lebedev and A.L. Skorokhodov,
!          "Quadrature Formulas of Orders 41, 47, and 53 for the Sphere,"
!          Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.
!
!      [4] V.I. Lebedev,
!          "Spherical Quadrature Formulas Exact to Orders 25-29,"
!          Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.
!
!      [5] V.I. Lebedev,
!          "Quadratures on a Sphere,"
!          Computational Mathematics and Mathematical Physics, Vol. 16,
!          1976, pp. 10-24.
!
!      [6] V.I. Lebedev,
!          "Values of the Nodes and Weights of Ninth to Seventeenth
!          Order Gauss-Markov Quadrature Formulae Invariant under the
!          Octahedron Group with Inversion,"
!          Computational Mathematics and Mathematical Physics, Vol. 15,
!          1975, pp. 44-51.
!
    n = 1
    v = 0.6666666666666667d-01
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.7500000000000000d-01
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0026(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV   26-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.4761904761904762d-01
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3809523809523810d-01
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3214285714285714d-01
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0038(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV   38-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.9523809523809524d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3214285714285714d-01
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4597008433809831d+00
    v = 0.2857142857142857d-01
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0050(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV   50-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.1269841269841270d-01
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2257495590828924d-01
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2109375000000000d-01
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3015113445777636d+00
    v = 0.2017333553791887d-01
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0074(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV   74-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.5130671797338464d-03
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1660406956574204d-01
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = -0.2958603896103896d-01
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4803844614152614d+00
    v = 0.2657620708215946d-01
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3207726489807764d+00
    v = 0.1652217099371571d-01
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0086(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV   86-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.1154401154401154d-01
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1194390908585628d-01
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3696028464541502d+00
    v = 0.1111055571060340d-01
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6943540066026664d+00
    v = 0.1187650129453714d-01
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3742430390903412d+00
    v = 0.1181230374690448d-01
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0110(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  110-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.3828270494937162d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.9793737512487512d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1851156353447362d+00
    v = 0.8211737283191111d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6904210483822922d+00
    v = 0.9942814891178103d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3956894730559419d+00
    v = 0.9595471336070963d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4783690288121502d+00
    v = 0.9694996361663028d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0146(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  146-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.5996313688621381d-03
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.7372999718620756d-02
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.7210515360144488d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6764410400114264d+00
    v = 0.7116355493117555d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4174961227965453d+00
    v = 0.6753829486314477d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1574676672039082d+00
    v = 0.7574394159054034d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1403553811713183d+00
    b = 0.4493328323269557d+00
    v = 0.6991087353303262d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0170(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  170-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.5544842902037365d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.6071332770670752d-02
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.6383674773515093d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2551252621114134d+00
    v = 0.5183387587747790d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6743601460362766d+00
    v = 0.6317929009813725d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4318910696719410d+00
    v = 0.6201670006589077d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2613931360335988d+00
    v = 0.5477143385137348d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4990453161796037d+00
    b = 0.1446630744325115d+00
    v = 0.5968383987681156d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0194(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  194-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.1782340447244611d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.5716905949977102d-02
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.5573383178848738d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6712973442695226d+00
    v = 0.5608704082587997d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2892465627575439d+00
    v = 0.5158237711805383d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4446933178717437d+00
    v = 0.5518771467273614d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1299335447650067d+00
    v = 0.4106777028169394d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3457702197611283d+00
    v = 0.5051846064614808d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1590417105383530d+00
    b = 0.8360360154824589d+00
    v = 0.5530248916233094d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0230(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  230-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = -0.5522639919727325d-01
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4450274607445226d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4492044687397611d+00
    v = 0.4496841067921404d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2520419490210201d+00
    v = 0.5049153450478750d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6981906658447242d+00
    v = 0.3976408018051883d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6587405243460960d+00
    v = 0.4401400650381014d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4038544050097660d-01
    v = 0.1724544350544401d-01
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5823842309715585d+00
    v = 0.4231083095357343d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3545877390518688d+00
    v = 0.5198069864064399d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2272181808998187d+00
    b = 0.4864661535886647d+00
    v = 0.4695720972568883d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0266(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  266-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = -0.1313769127326952d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = -0.2522728704859336d-02
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4186853881700583d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7039373391585475d+00
    v = 0.5315167977810885d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1012526248572414d+00
    v = 0.4047142377086219d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4647448726420539d+00
    v = 0.4112482394406990d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3277420654971629d+00
    v = 0.3595584899758782d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6620338663699974d+00
    v = 0.4256131351428158d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8506508083520399d+00
    v = 0.4229582700647240d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3233484542692899d+00
    b = 0.1153112011009701d+00
    v = 0.4080914225780505d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2314790158712601d+00
    b = 0.5244939240922365d+00
    v = 0.4071467593830964d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0302(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  302-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.8545911725128148d-03
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3599119285025571d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3515640345570105d+00
    v = 0.3449788424305883d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6566329410219612d+00
    v = 0.3604822601419882d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4729054132581005d+00
    v = 0.3576729661743367d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9618308522614784d-01
    v = 0.2352101413689164d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2219645236294178d+00
    v = 0.3108953122413675d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7011766416089545d+00
    v = 0.3650045807677255d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2644152887060663d+00
    v = 0.2982344963171804d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5718955891878961d+00
    v = 0.3600820932216460d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2510034751770465d+00
    b = 0.8000727494073952d+00
    v = 0.3571540554273387d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1233548532583327d+00
    b = 0.4127724083168531d+00
    v = 0.3392312205006170d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0350(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  350-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.3006796749453936d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3050627745650771d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7068965463912316d+00
    v = 0.1621104600288991d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4794682625712025d+00
    v = 0.3005701484901752d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1927533154878019d+00
    v = 0.2990992529653774d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6930357961327123d+00
    v = 0.2982170644107595d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3608302115520091d+00
    v = 0.2721564237310992d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6498486161496169d+00
    v = 0.3033513795811141d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1932945013230339d+00
    v = 0.3007949555218533d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3800494919899303d+00
    v = 0.2881964603055307d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2899558825499574d+00
    b = 0.7934537856582316d+00
    v = 0.2958357626535696d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9684121455103957d-01
    b = 0.8280801506686862d+00
    v = 0.3036020026407088d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1833434647041659d+00
    b = 0.9074658265305127d+00
    v = 0.2832187403926303d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0434(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  434-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.5265897968224436d-03
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2548219972002607d-02
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2512317418927307d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6909346307509111d+00
    v = 0.2530403801186355d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1774836054609158d+00
    v = 0.2014279020918528d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4914342637784746d+00
    v = 0.2501725168402936d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6456664707424256d+00
    v = 0.2513267174597564d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2861289010307638d+00
    v = 0.2302694782227416d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7568084367178018d-01
    v = 0.1462495621594614d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3927259763368002d+00
    v = 0.2445373437312980d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8818132877794288d+00
    v = 0.2417442375638981d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9776428111182649d+00
    v = 0.1910951282179532d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2054823696403044d+00
    b = 0.8689460322872412d+00
    v = 0.2416930044324775d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5905157048925271d+00
    b = 0.7999278543857286d+00
    v = 0.2512236854563495d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5550152361076807d+00
    b = 0.7717462626915901d+00
    v = 0.2496644054553086d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9371809858553722d+00
    b = 0.3344363145343455d+00
    v = 0.2236607760437849d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0590(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  590-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.3095121295306187d-03
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1852379698597489d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7040954938227469d+00
    v = 0.1871790639277744d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6807744066455243d+00
    v = 0.1858812585438317d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6372546939258752d+00
    v = 0.1852028828296213d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5044419707800358d+00
    v = 0.1846715956151242d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4215761784010967d+00
    v = 0.1818471778162769d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3317920736472123d+00
    v = 0.1749564657281154d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2384736701421887d+00
    v = 0.1617210647254411d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1459036449157763d+00
    v = 0.1384737234851692d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6095034115507196d-01
    v = 0.9764331165051050d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6116843442009876d+00
    v = 0.1857161196774078d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3964755348199858d+00
    v = 0.1705153996395864d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1724782009907724d+00
    v = 0.1300321685886048d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5610263808622060d+00
    b = 0.3518280927733519d+00
    v = 0.1842866472905286d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4742392842551980d+00
    b = 0.2634716655937950d+00
    v = 0.1802658934377451d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5984126497885380d+00
    b = 0.1816640840360209d+00
    v = 0.1849830560443660d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3791035407695563d+00
    b = 0.1720795225656878d+00
    v = 0.1713904507106709d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2778673190586244d+00
    b = 0.8213021581932511d-01
    v = 0.1555213603396808d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5033564271075117d+00
    b = 0.8999205842074875d-01
    v = 0.1802239128008525d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0770(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  770-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.2192942088181184d-03
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1436433617319080d-02
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1421940344335877d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5087204410502360d-01
    v = 0.6798123511050502d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1228198790178831d+00
    v = 0.9913184235294912d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2026890814408786d+00
    v = 0.1180207833238949d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2847745156464294d+00
    v = 0.1296599602080921d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3656719078978026d+00
    v = 0.1365871427428316d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4428264886713469d+00
    v = 0.1402988604775325d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5140619627249735d+00
    v = 0.1418645563595609d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6306401219166803d+00
    v = 0.1421376741851662d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6716883332022612d+00
    v = 0.1423996475490962d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6979792685336881d+00
    v = 0.1431554042178567d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1446865674195309d+00
    v = 0.9254401499865368d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3390263475411216d+00
    v = 0.1250239995053509d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5335804651263506d+00
    v = 0.1394365843329230d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6944024393349413d-01
    b = 0.2355187894242326d+00
    v = 0.1127089094671749d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2269004109529460d+00
    b = 0.4102182474045730d+00
    v = 0.1345753760910670d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8025574607775339d-01
    b = 0.6214302417481605d+00
    v = 0.1424957283316783d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1467999527896572d+00
    b = 0.3245284345717394d+00
    v = 0.1261523341237750d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1571507769824727d+00
    b = 0.5224482189696630d+00
    v = 0.1392547106052696d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2365702993157246d+00
    b = 0.6017546634089558d+00
    v = 0.1418761677877656d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7714815866765732d-01
    b = 0.4346575516141163d+00
    v = 0.1338366684479554d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3062936666210730d+00
    b = 0.4908826589037616d+00
    v = 0.1393700862676131d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3822477379524787d+00
    b = 0.5648768149099500d+00
    v = 0.1415914757466932d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld0974(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV  974-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.1438294190527431d-03
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1125772288287004d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4292963545341347d-01
    v = 0.4948029341949241d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1051426854086404d+00
    v = 0.7357990109125470d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1750024867623087d+00
    v = 0.8889132771304384d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2477653379650257d+00
    v = 0.9888347838921435d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3206567123955957d+00
    v = 0.1053299681709471d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3916520749849983d+00
    v = 0.1092778807014578d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4590825874187624d+00
    v = 0.1114389394063227d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5214563888415861d+00
    v = 0.1123724788051555d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6253170244654199d+00
    v = 0.1125239325243814d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6637926744523170d+00
    v = 0.1126153271815905d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6910410398498301d+00
    v = 0.1130286931123841d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7052907007457760d+00
    v = 0.1134986534363955d-02
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1236686762657990d+00
    v = 0.6823367927109931d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2940777114468387d+00
    v = 0.9454158160447096d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4697753849207649d+00
    v = 0.1074429975385679d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6334563241139567d+00
    v = 0.1129300086569132d-02
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5974048614181342d-01
    b = 0.2029128752777523d+00
    v = 0.8436884500901954d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1375760408473636d+00
    b = 0.4602621942484054d+00
    v = 0.1075255720448885d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3391016526336286d+00
    b = 0.5030673999662036d+00
    v = 0.1108577236864462d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1271675191439820d+00
    b = 0.2817606422442134d+00
    v = 0.9566475323783357d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2693120740413512d+00
    b = 0.4331561291720157d+00
    v = 0.1080663250717391d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1419786452601918d+00
    b = 0.6256167358580814d+00
    v = 0.1126797131196295d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6709284600738255d-01
    b = 0.3798395216859157d+00
    v = 0.1022568715358061d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7057738183256172d-01
    b = 0.5517505421423520d+00
    v = 0.1108960267713108d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2783888477882155d+00
    b = 0.6029619156159187d+00
    v = 0.1122790653435766d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1979578938917407d+00
    b = 0.3589606329589096d+00
    v = 0.1032401847117460d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2087307061103274d+00
    b = 0.5348666438135476d+00
    v = 0.1107249382283854d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4055122137872836d+00
    b = 0.5674997546074373d+00
    v = 0.1121780048519972d-02
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld1202(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 1202-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.1105189233267572d-03
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.9205232738090741d-03
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.9133159786443561d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3712636449657089d-01
    v = 0.3690421898017899d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9140060412262223d-01
    v = 0.5603990928680660d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1531077852469906d+00
    v = 0.6865297629282609d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2180928891660612d+00
    v = 0.7720338551145630d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2839874532200175d+00
    v = 0.8301545958894795d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3491177600963764d+00
    v = 0.8686692550179628d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4121431461444309d+00
    v = 0.8927076285846890d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4718993627149127d+00
    v = 0.9060820238568219d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5273145452842337d+00
    v = 0.9119777254940867d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6209475332444019d+00
    v = 0.9128720138604181d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6569722711857291d+00
    v = 0.9130714935691735d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6841788309070143d+00
    v = 0.9152873784554116d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7012604330123631d+00
    v = 0.9187436274321654d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1072382215478166d+00
    v = 0.5176977312965694d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2582068959496968d+00
    v = 0.7331143682101417d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4172752955306717d+00
    v = 0.8463232836379928d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5700366911792503d+00
    v = 0.9031122694253992d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9827986018263947d+00
    b = 0.1771774022615325d+00
    v = 0.6485778453163257d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9624249230326228d+00
    b = 0.2475716463426288d+00
    v = 0.7435030910982369d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9402007994128811d+00
    b = 0.3354616289066489d+00
    v = 0.7998527891839054d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9320822040143202d+00
    b = 0.3173615246611977d+00
    v = 0.8101731497468018d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9043674199393299d+00
    b = 0.4090268427085357d+00
    v = 0.8483389574594331d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8912407560074747d+00
    b = 0.3854291150669224d+00
    v = 0.8556299257311812d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8676435628462708d+00
    b = 0.4932221184851285d+00
    v = 0.8803208679738260d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8581979986041619d+00
    b = 0.4785320675922435d+00
    v = 0.8811048182425720d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8396753624049856d+00
    b = 0.4507422593157064d+00
    v = 0.8850282341265444d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8165288564022188d+00
    b = 0.5632123020762100d+00
    v = 0.9021342299040653d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8015469370783529d+00
    b = 0.5434303569693900d+00
    v = 0.9010091677105086d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7773563069070351d+00
    b = 0.5123518486419871d+00
    v = 0.9022692938426915d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7661621213900394d+00
    b = 0.6394279634749102d+00
    v = 0.9158016174693465d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7553584143533510d+00
    b = 0.6269805509024392d+00
    v = 0.9131578003189435d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7344305757559503d+00
    b = 0.6031161693096310d+00
    v = 0.9107813579482705d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7043837184021765d+00
    b = 0.5693702498468441d+00
    v = 0.9105760258970126d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld1454(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 1454-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.7777160743261247d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.7557646413004701d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3229290663413854d-01
    v = 0.2841633806090617d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8036733271462222d-01
    v = 0.4374419127053555d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1354289960531653d+00
    v = 0.5417174740872172d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1938963861114426d+00
    v = 0.6148000891358593d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2537343715011275d+00
    v = 0.6664394485800705d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3135251434752570d+00
    v = 0.7025039356923220d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3721558339375338d+00
    v = 0.7268511789249627d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4286809575195696d+00
    v = 0.7422637534208629d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4822510128282994d+00
    v = 0.7509545035841214d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5320679333566263d+00
    v = 0.7548535057718401d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6172998195394274d+00
    v = 0.7554088969774001d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6510679849127481d+00
    v = 0.7553147174442808d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6777315251687360d+00
    v = 0.7564767653292297d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6963109410648741d+00
    v = 0.7587991808518730d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7058935009831749d+00
    v = 0.7608261832033027d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9955546194091857d+00
    v = 0.4021680447874916d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9734115901794209d+00
    v = 0.5804871793945964d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9275693732388626d+00
    v = 0.6792151955945159d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8568022422795103d+00
    v = 0.7336741211286294d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7623495553719372d+00
    v = 0.7581866300989608d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5707522908892223d+00
    b = 0.4387028039889501d+00
    v = 0.7538257859800743d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5196463388403083d+00
    b = 0.3858908414762617d+00
    v = 0.7483517247053123d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4646337531215351d+00
    b = 0.3301937372343854d+00
    v = 0.7371763661112059d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4063901697557691d+00
    b = 0.2725423573563777d+00
    v = 0.7183448895756934d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3456329466643087d+00
    b = 0.2139510237495250d+00
    v = 0.6895815529822191d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2831395121050332d+00
    b = 0.1555922309786647d+00
    v = 0.6480105801792886d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2197682022925330d+00
    b = 0.9892878979686097d-01
    v = 0.5897558896594636d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1564696098650355d+00
    b = 0.4598642910675510d-01
    v = 0.5095708849247346d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6027356673721295d+00
    b = 0.3376625140173426d+00
    v = 0.7536906428909755d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5496032320255096d+00
    b = 0.2822301309727988d+00
    v = 0.7472505965575118d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4921707755234567d+00
    b = 0.2248632342592540d+00
    v = 0.7343017132279698d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4309422998598483d+00
    b = 0.1666224723456479d+00
    v = 0.7130871582177445d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3664108182313672d+00
    b = 0.1086964901822169d+00
    v = 0.6817022032112776d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2990189057758436d+00
    b = 0.5251989784120085d-01
    v = 0.6380941145604121d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6268724013144998d+00
    b = 0.2297523657550023d+00
    v = 0.7550381377920310d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5707324144834607d+00
    b = 0.1723080607093800d+00
    v = 0.7478646640144802d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5096360901960365d+00
    b = 0.1140238465390513d+00
    v = 0.7335918720601220d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4438729938312456d+00
    b = 0.5611522095882537d-01
    v = 0.7110120527658118d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6419978471082389d+00
    b = 0.1164174423140873d+00
    v = 0.7571363978689501d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5817218061802611d+00
    b = 0.5797589531445219d-01
    v = 0.7489908329079234d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld1730(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 1730-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.6309049437420976d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.6398287705571748d-03
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.6357185073530720d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2860923126194662d-01
    v = 0.2221207162188168d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7142556767711522d-01
    v = 0.3475784022286848d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1209199540995559d+00
    v = 0.4350742443589804d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1738673106594379d+00
    v = 0.4978569136522127d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2284645438467734d+00
    v = 0.5435036221998053d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2834807671701512d+00
    v = 0.5765913388219542d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3379680145467339d+00
    v = 0.6001200359226003d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3911355454819537d+00
    v = 0.6162178172717512d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4422860353001403d+00
    v = 0.6265218152438485d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4907781568726057d+00
    v = 0.6323987160974212d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5360006153211468d+00
    v = 0.6350767851540569d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6142105973596603d+00
    v = 0.6354362775297107d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6459300387977504d+00
    v = 0.6352302462706235d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6718056125089225d+00
    v = 0.6358117881417972d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6910888533186254d+00
    v = 0.6373101590310117d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7030467416823252d+00
    v = 0.6390428961368665d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8354951166354646d-01
    v = 0.3186913449946576d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2050143009099486d+00
    v = 0.4678028558591711d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3370208290706637d+00
    v = 0.5538829697598626d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4689051484233963d+00
    v = 0.6044475907190476d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5939400424557334d+00
    v = 0.6313575103509012d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1394983311832261d+00
    b = 0.4097581162050343d-01
    v = 0.4078626431855630d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1967999180485014d+00
    b = 0.8851987391293348d-01
    v = 0.4759933057812725d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2546183732548967d+00
    b = 0.1397680182969819d+00
    v = 0.5268151186413440d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3121281074713875d+00
    b = 0.1929452542226526d+00
    v = 0.5643048560507316d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3685981078502492d+00
    b = 0.2467898337061562d+00
    v = 0.5914501076613073d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4233760321547856d+00
    b = 0.3003104124785409d+00
    v = 0.6104561257874195d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4758671236059246d+00
    b = 0.3526684328175033d+00
    v = 0.6230252860707806d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5255178579796463d+00
    b = 0.4031134861145713d+00
    v = 0.6305618761760796d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5718025633734589d+00
    b = 0.4509426448342351d+00
    v = 0.6343092767597889d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2686927772723415d+00
    b = 0.4711322502423248d-01
    v = 0.5176268945737826d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3306006819904809d+00
    b = 0.9784487303942695d-01
    v = 0.5564840313313692d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3904906850594983d+00
    b = 0.1505395810025273d+00
    v = 0.5856426671038980d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4479957951904390d+00
    b = 0.2039728156296050d+00
    v = 0.6066386925777091d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5027076848919780d+00
    b = 0.2571529941121107d+00
    v = 0.6208824962234458d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5542087392260217d+00
    b = 0.3092191375815670d+00
    v = 0.6296314297822907d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6020850887375187d+00
    b = 0.3593807506130276d+00
    v = 0.6340423756791859d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4019851409179594d+00
    b = 0.5063389934378671d-01
    v = 0.5829627677107342d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4635614567449800d+00
    b = 0.1032422269160612d+00
    v = 0.6048693376081110d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5215860931591575d+00
    b = 0.1566322094006254d+00
    v = 0.6202362317732461d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5758202499099271d+00
    b = 0.2098082827491099d+00
    v = 0.6299005328403779d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6259893683876795d+00
    b = 0.2618824114553391d+00
    v = 0.6347722390609353d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5313795124811891d+00
    b = 0.5263245019338556d-01
    v = 0.6203778981238834d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5893317955931995d+00
    b = 0.1061059730982005d+00
    v = 0.6308414671239979d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6426246321215801d+00
    b = 0.1594171564034221d+00
    v = 0.6362706466959498d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6511904367376113d+00
    b = 0.5354789536565540d-01
    v = 0.6375414170333233d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld2030(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 2030-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.4656031899197431d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.5421549195295507d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2540835336814348d-01
    v = 0.1778522133346553d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6399322800504915d-01
    v = 0.2811325405682796d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1088269469804125d+00
    v = 0.3548896312631459d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1570670798818287d+00
    v = 0.4090310897173364d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2071163932282514d+00
    v = 0.4493286134169965d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2578914044450844d+00
    v = 0.4793728447962723d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3085687558169623d+00
    v = 0.5015415319164265d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3584719706267024d+00
    v = 0.5175127372677937d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4070135594428709d+00
    v = 0.5285522262081019d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4536618626222638d+00
    v = 0.5356832703713962d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4979195686463577d+00
    v = 0.5397914736175170d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5393075111126999d+00
    v = 0.5416899441599930d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6115617676843916d+00
    v = 0.5419308476889938d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6414308435160159d+00
    v = 0.5416936902030596d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6664099412721607d+00
    v = 0.5419544338703164d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6859161771214913d+00
    v = 0.5428983656630975d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6993625593503890d+00
    v = 0.5442286500098193d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7062393387719380d+00
    v = 0.5452250345057301d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7479028168349763d-01
    v = 0.2568002497728530d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1848951153969366d+00
    v = 0.3827211700292145d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3059529066581305d+00
    v = 0.4579491561917824d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4285556101021362d+00
    v = 0.5042003969083574d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5468758653496526d+00
    v = 0.5312708889976025d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6565821978343439d+00
    v = 0.5438401790747117d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1253901572367117d+00
    b = 0.3681917226439641d-01
    v = 0.3316041873197344d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1775721510383941d+00
    b = 0.7982487607213301d-01
    v = 0.3899113567153771d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2305693358216114d+00
    b = 0.1264640966592335d+00
    v = 0.4343343327201309d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2836502845992063d+00
    b = 0.1751585683418957d+00
    v = 0.4679415262318919d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3361794746232590d+00
    b = 0.2247995907632670d+00
    v = 0.4930847981631031d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3875979172264824d+00
    b = 0.2745299257422246d+00
    v = 0.5115031867540091d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4374019316999074d+00
    b = 0.3236373482441118d+00
    v = 0.5245217148457367d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4851275843340022d+00
    b = 0.3714967859436741d+00
    v = 0.5332041499895321d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5303391803806868d+00
    b = 0.4175353646321745d+00
    v = 0.5384583126021542d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5726197380596287d+00
    b = 0.4612084406355461d+00
    v = 0.5411067210798852d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2431520732564863d+00
    b = 0.4258040133043952d-01
    v = 0.4259797391468714d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3002096800895869d+00
    b = 0.8869424306722721d-01
    v = 0.4604931368460021d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3558554457457432d+00
    b = 0.1368811706510655d+00
    v = 0.4871814878255202d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4097782537048887d+00
    b = 0.1860739985015033d+00
    v = 0.5072242910074885d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4616337666067458d+00
    b = 0.2354235077395853d+00
    v = 0.5217069845235350d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5110707008417874d+00
    b = 0.2842074921347011d+00
    v = 0.5315785966280310d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5577415286163795d+00
    b = 0.3317784414984102d+00
    v = 0.5376833708758905d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6013060431366950d+00
    b = 0.3775299002040700d+00
    v = 0.5408032092069521d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3661596767261781d+00
    b = 0.4599367887164592d-01
    v = 0.4842744917904866d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4237633153506581d+00
    b = 0.9404893773654421d-01
    v = 0.5048926076188130d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4786328454658452d+00
    b = 0.1431377109091971d+00
    v = 0.5202607980478373d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5305702076789774d+00
    b = 0.1924186388843570d+00
    v = 0.5309932388325743d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5793436224231788d+00
    b = 0.2411590944775190d+00
    v = 0.5377419770895208d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6247069017094747d+00
    b = 0.2886871491583605d+00
    v = 0.5411696331677717d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4874315552535204d+00
    b = 0.4804978774953206d-01
    v = 0.5197996293282420d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5427337322059053d+00
    b = 0.9716857199366665d-01
    v = 0.5311120836622945d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5943493747246700d+00
    b = 0.1465205839795055d+00
    v = 0.5384309319956951d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6421314033564943d+00
    b = 0.1953579449803574d+00
    v = 0.5421859504051886d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6020628374713980d+00
    b = 0.4916375015738108d-01
    v = 0.5390948355046314d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6529222529856881d+00
    b = 0.9861621540127005d-01
    v = 0.5433312705027845d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld2354(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 2354-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.3922616270665292d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4703831750854424d-03
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4678202801282136d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2290024646530589d-01
    v = 0.1437832228979900d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5779086652271284d-01
    v = 0.2303572493577644d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9863103576375984d-01
    v = 0.2933110752447454d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1428155792982185d+00
    v = 0.3402905998359838d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1888978116601463d+00
    v = 0.3759138466870372d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2359091682970210d+00
    v = 0.4030638447899798d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2831228833706171d+00
    v = 0.4236591432242211d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3299495857966693d+00
    v = 0.4390522656946746d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3758840802660796d+00
    v = 0.4502523466626247d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4204751831009480d+00
    v = 0.4580577727783541d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4633068518751051d+00
    v = 0.4631391616615899d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5039849474507313d+00
    v = 0.4660928953698676d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5421265793440747d+00
    v = 0.4674751807936953d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6092660230557310d+00
    v = 0.4676414903932920d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6374654204984869d+00
    v = 0.4674086492347870d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6615136472609892d+00
    v = 0.4674928539483207d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6809487285958127d+00
    v = 0.4680748979686447d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6952980021665196d+00
    v = 0.4690449806389040d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7041245497695400d+00
    v = 0.4699877075860818d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6744033088306065d-01
    v = 0.2099942281069176d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1678684485334166d+00
    v = 0.3172269150712804d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2793559049539613d+00
    v = 0.3832051358546523d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3935264218057639d+00
    v = 0.4252193818146985d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5052629268232558d+00
    v = 0.4513807963755000d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6107905315437531d+00
    v = 0.4657797469114178d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1135081039843524d+00
    b = 0.3331954884662588d-01
    v = 0.2733362800522836d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1612866626099378d+00
    b = 0.7247167465436538d-01
    v = 0.3235485368463559d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2100786550168205d+00
    b = 0.1151539110849745d+00
    v = 0.3624908726013453d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2592282009459942d+00
    b = 0.1599491097143677d+00
    v = 0.3925540070712828d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3081740561320203d+00
    b = 0.2058699956028027d+00
    v = 0.4156129781116235d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3564289781578164d+00
    b = 0.2521624953502911d+00
    v = 0.4330644984623263d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4035587288240703d+00
    b = 0.2982090785797674d+00
    v = 0.4459677725921312d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4491671196373903d+00
    b = 0.3434762087235733d+00
    v = 0.4551593004456795d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4928854782917489d+00
    b = 0.3874831357203437d+00
    v = 0.4613341462749918d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5343646791958988d+00
    b = 0.4297814821746926d+00
    v = 0.4651019618269806d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5732683216530990d+00
    b = 0.4699402260943537d+00
    v = 0.4670249536100625d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2214131583218986d+00
    b = 0.3873602040643895d-01
    v = 0.3549555576441708d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2741796504750071d+00
    b = 0.8089496256902013d-01
    v = 0.3856108245249010d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3259797439149485d+00
    b = 0.1251732177620872d+00
    v = 0.4098622845756882d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3765441148826891d+00
    b = 0.1706260286403185d+00
    v = 0.4286328604268950d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4255773574530558d+00
    b = 0.2165115147300408d+00
    v = 0.4427802198993945d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4727795117058430d+00
    b = 0.2622089812225259d+00
    v = 0.4530473511488561d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5178546895819012d+00
    b = 0.3071721431296201d+00
    v = 0.4600805475703138d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5605141192097460d+00
    b = 0.3508998998801138d+00
    v = 0.4644599059958017d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6004763319352512d+00
    b = 0.3929160876166931d+00
    v = 0.4667274455712508d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3352842634946949d+00
    b = 0.4202563457288019d-01
    v = 0.4069360518020356d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3891971629814670d+00
    b = 0.8614309758870850d-01
    v = 0.4260442819919195d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4409875565542281d+00
    b = 0.1314500879380001d+00
    v = 0.4408678508029063d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4904893058592484d+00
    b = 0.1772189657383859d+00
    v = 0.4518748115548597d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5375056138769549d+00
    b = 0.2228277110050294d+00
    v = 0.4595564875375116d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5818255708669969d+00
    b = 0.2677179935014386d+00
    v = 0.4643988774315846d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6232334858144959d+00
    b = 0.3113675035544165d+00
    v = 0.4668827491646946d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4489485354492058d+00
    b = 0.4409162378368174d-01
    v = 0.4400541823741973d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5015136875933150d+00
    b = 0.8939009917748489d-01
    v = 0.4514512890193797d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5511300550512623d+00
    b = 0.1351806029383365d+00
    v = 0.4596198627347549d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5976720409858000d+00
    b = 0.1808370355053196d+00
    v = 0.4648659016801781d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6409956378989354d+00
    b = 0.2257852192301602d+00
    v = 0.4675502017157673d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5581222330827514d+00
    b = 0.4532173421637160d-01
    v = 0.4598494476455523d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6074705984161695d+00
    b = 0.9117488031840314d-01
    v = 0.4654916955152048d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6532272537379033d+00
    b = 0.1369294213140155d+00
    v = 0.4684709779505137d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6594761494500487d+00
    b = 0.4589901487275583d-01
    v = 0.4691445539106986d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld2702(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 2702-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.2998675149888161d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4077860529495355d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2065562538818703d-01
    v = 0.1185349192520667d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5250918173022379d-01
    v = 0.1913408643425751d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8993480082038376d-01
    v = 0.2452886577209897d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1306023924436019d+00
    v = 0.2862408183288702d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1732060388531418d+00
    v = 0.3178032258257357d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2168727084820249d+00
    v = 0.3422945667633690d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2609528309173586d+00
    v = 0.3612790520235922d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3049252927938952d+00
    v = 0.3758638229818521d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3483484138084404d+00
    v = 0.3868711798859953d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3908321549106406d+00
    v = 0.3949429933189938d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4320210071894814d+00
    v = 0.4006068107541156d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4715824795890053d+00
    v = 0.4043192149672723d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5091984794078453d+00
    v = 0.4064947495808078d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5445580145650803d+00
    v = 0.4075245619813152d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6072575796841768d+00
    v = 0.4076423540893566d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6339484505755803d+00
    v = 0.4074280862251555d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6570718257486958d+00
    v = 0.4074163756012244d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6762557330090709d+00
    v = 0.4077647795071246d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6911161696923790d+00
    v = 0.4084517552782530d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7012841911659961d+00
    v = 0.4092468459224052d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7064559272410020d+00
    v = 0.4097872687240906d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6123554989894765d-01
    v = 0.1738986811745028d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1533070348312393d+00
    v = 0.2659616045280191d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2563902605244206d+00
    v = 0.3240596008171533d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3629346991663361d+00
    v = 0.3621195964432943d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4683949968987538d+00
    v = 0.3868838330760539d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5694479240657952d+00
    v = 0.4018911532693111d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6634465430993955d+00
    v = 0.4089929432983252d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1033958573552305d+00
    b = 0.3034544009063584d-01
    v = 0.2279907527706409d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1473521412414395d+00
    b = 0.6618803044247135d-01
    v = 0.2715205490578897d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1924552158705967d+00
    b = 0.1054431128987715d+00
    v = 0.3057917896703976d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2381094362890328d+00
    b = 0.1468263551238858d+00
    v = 0.3326913052452555d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2838121707936760d+00
    b = 0.1894486108187886d+00
    v = 0.3537334711890037d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3291323133373415d+00
    b = 0.2326374238761579d+00
    v = 0.3700567500783129d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3736896978741460d+00
    b = 0.2758485808485768d+00
    v = 0.3825245372589122d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4171406040760013d+00
    b = 0.3186179331996921d+00
    v = 0.3918125171518296d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4591677985256915d+00
    b = 0.3605329796303794d+00
    v = 0.3984720419937579d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4994733831718418d+00
    b = 0.4012147253586509d+00
    v = 0.4029746003338211d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5377731830445096d+00
    b = 0.4403050025570692d+00
    v = 0.4057428632156627d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5737917830001331d+00
    b = 0.4774565904277483d+00
    v = 0.4071719274114857d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2027323586271389d+00
    b = 0.3544122504976147d-01
    v = 0.2990236950664119d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2516942375187273d+00
    b = 0.7418304388646328d-01
    v = 0.3262951734212878d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3000227995257181d+00
    b = 0.1150502745727186d+00
    v = 0.3482634608242413d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3474806691046342d+00
    b = 0.1571963371209364d+00
    v = 0.3656596681700892d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3938103180359209d+00
    b = 0.1999631877247100d+00
    v = 0.3791740467794218d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4387519590455703d+00
    b = 0.2428073457846535d+00
    v = 0.3894034450156905d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4820503960077787d+00
    b = 0.2852575132906155d+00
    v = 0.3968600245508371d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5234573778475101d+00
    b = 0.3268884208674639d+00
    v = 0.4019931351420050d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5627318647235282d+00
    b = 0.3673033321675939d+00
    v = 0.4052108801278599d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5996390607156954d+00
    b = 0.4061211551830290d+00
    v = 0.4068978613940934d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3084780753791947d+00
    b = 0.3860125523100059d-01
    v = 0.3454275351319704d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3589988275920223d+00
    b = 0.7928938987104867d-01
    v = 0.3629963537007920d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4078628415881973d+00
    b = 0.1212614643030087d+00
    v = 0.3770187233889873d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4549287258889735d+00
    b = 0.1638770827382693d+00
    v = 0.3878608613694378d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5000278512957279d+00
    b = 0.2065965798260176d+00
    v = 0.3959065270221274d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5429785044928199d+00
    b = 0.2489436378852235d+00
    v = 0.4015286975463570d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5835939850491711d+00
    b = 0.2904811368946891d+00
    v = 0.4050866785614717d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6216870353444856d+00
    b = 0.3307941957666609d+00
    v = 0.4069320185051913d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4151104662709091d+00
    b = 0.4064829146052554d-01
    v = 0.3760120964062763d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4649804275009218d+00
    b = 0.8258424547294755d-01
    v = 0.3870969564418064d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5124695757009662d+00
    b = 0.1251841962027289d+00
    v = 0.3955287790534055d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5574711100606224d+00
    b = 0.1679107505976331d+00
    v = 0.4015361911302668d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5998597333287227d+00
    b = 0.2102805057358715d+00
    v = 0.4053836986719548d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6395007148516600d+00
    b = 0.2518418087774107d+00
    v = 0.4073578673299117d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5188456224746252d+00
    b = 0.4194321676077518d-01
    v = 0.3954628379231406d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5664190707942778d+00
    b = 0.8457661551921499d-01
    v = 0.4017645508847530d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6110464353283153d+00
    b = 0.1273652932519396d+00
    v = 0.4059030348651293d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6526430302051563d+00
    b = 0.1698173239076354d+00
    v = 0.4080565809484880d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6167551880377548d+00
    b = 0.4266398851548864d-01
    v = 0.4063018753664651d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6607195418355383d+00
    b = 0.8551925814238349d-01
    v = 0.4087191292799671d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld3074(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 3074-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.2599095953754734d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3603134089687541d-03
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3586067974412447d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1886108518723392d-01
    v = 0.9831528474385880d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4800217244625303d-01
    v = 0.1605023107954450d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8244922058397242d-01
    v = 0.2072200131464099d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1200408362484023d+00
    v = 0.2431297618814187d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1595773530809965d+00
    v = 0.2711819064496707d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2002635973434064d+00
    v = 0.2932762038321116d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2415127590139982d+00
    v = 0.3107032514197368d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2828584158458477d+00
    v = 0.3243808058921213d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3239091015338138d+00
    v = 0.3349899091374030d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3643225097962194d+00
    v = 0.3430580688505218d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4037897083691802d+00
    v = 0.3490124109290343d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4420247515194127d+00
    v = 0.3532148948561955d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4787572538464938d+00
    v = 0.3559862669062833d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5137265251275234d+00
    v = 0.3576224317551411d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5466764056654611d+00
    v = 0.3584050533086076d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6054859420813535d+00
    v = 0.3584903581373224d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6308106701764562d+00
    v = 0.3582991879040586d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6530369230179584d+00
    v = 0.3582371187963125d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6718609524611158d+00
    v = 0.3584353631122350d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6869676499894013d+00
    v = 0.3589120166517785d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6980467077240748d+00
    v = 0.3595445704531601d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7048241721250522d+00
    v = 0.3600943557111074d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5591105222058232d-01
    v = 0.1456447096742039d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1407384078513916d+00
    v = 0.2252370188283782d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2364035438976309d+00
    v = 0.2766135443474897d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3360602737818170d+00
    v = 0.3110729491500851d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4356292630054665d+00
    v = 0.3342506712303391d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5321569415256174d+00
    v = 0.3491981834026860d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6232956305040554d+00
    v = 0.3576003604348932d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9469870086838469d-01
    b = 0.2778748387309470d-01
    v = 0.1921921305788564d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1353170300568141d+00
    b = 0.6076569878628364d-01
    v = 0.2301458216495632d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1771679481726077d+00
    b = 0.9703072762711040d-01
    v = 0.2604248549522893d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2197066664231751d+00
    b = 0.1354112458524762d+00
    v = 0.2845275425870697d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2624783557374927d+00
    b = 0.1750996479744100d+00
    v = 0.3036870897974840d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3050969521214442d+00
    b = 0.2154896907449802d+00
    v = 0.3188414832298066d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3472252637196021d+00
    b = 0.2560954625740152d+00
    v = 0.3307046414722089d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3885610219026360d+00
    b = 0.2965070050624096d+00
    v = 0.3398330969031360d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4288273776062765d+00
    b = 0.3363641488734497d+00
    v = 0.3466757899705373d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4677662471302948d+00
    b = 0.3753400029836788d+00
    v = 0.3516095923230054d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5051333589553359d+00
    b = 0.4131297522144286d+00
    v = 0.3549645184048486d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5406942145810492d+00
    b = 0.4494423776081795d+00
    v = 0.3570415969441392d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5742204122576457d+00
    b = 0.4839938958841502d+00
    v = 0.3581251798496118d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1865407027225188d+00
    b = 0.3259144851070796d-01
    v = 0.2543491329913348d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2321186453689432d+00
    b = 0.6835679505297343d-01
    v = 0.2786711051330776d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2773159142523882d+00
    b = 0.1062284864451989d+00
    v = 0.2985552361083679d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3219200192237254d+00
    b = 0.1454404409323047d+00
    v = 0.3145867929154039d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3657032593944029d+00
    b = 0.1854018282582510d+00
    v = 0.3273290662067609d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4084376778363622d+00
    b = 0.2256297412014750d+00
    v = 0.3372705511943501d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4499004945751427d+00
    b = 0.2657104425000896d+00
    v = 0.3448274437851510d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4898758141326335d+00
    b = 0.3052755487631557d+00
    v = 0.3503592783048583d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5281547442266309d+00
    b = 0.3439863920645423d+00
    v = 0.3541854792663162d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5645346989813992d+00
    b = 0.3815229456121914d+00
    v = 0.3565995517909428d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5988181252159848d+00
    b = 0.4175752420966734d+00
    v = 0.3578802078302898d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2850425424471603d+00
    b = 0.3562149509862536d-01
    v = 0.2958644592860982d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3324619433027876d+00
    b = 0.7330318886871096d-01
    v = 0.3119548129116835d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3785848333076282d+00
    b = 0.1123226296008472d+00
    v = 0.3250745225005984d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4232891028562115d+00
    b = 0.1521084193337708d+00
    v = 0.3355153415935208d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4664287050829722d+00
    b = 0.1921844459223610d+00
    v = 0.3435847568549328d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5078458493735726d+00
    b = 0.2321360989678303d+00
    v = 0.3495786831622488d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5473779816204180d+00
    b = 0.2715886486360520d+00
    v = 0.3537767805534621d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5848617133811376d+00
    b = 0.3101924707571355d+00
    v = 0.3564459815421428d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6201348281584888d+00
    b = 0.3476121052890973d+00
    v = 0.3578464061225468d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3852191185387871d+00
    b = 0.3763224880035108d-01
    v = 0.3239748762836212d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4325025061073423d+00
    b = 0.7659581935637135d-01
    v = 0.3345491784174287d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4778486229734490d+00
    b = 0.1163381306083900d+00
    v = 0.3429126177301782d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5211663693009000d+00
    b = 0.1563890598752899d+00
    v = 0.3492420343097421d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5623469504853703d+00
    b = 0.1963320810149200d+00
    v = 0.3537399050235257d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6012718188659246d+00
    b = 0.2357847407258738d+00
    v = 0.3566209152659172d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6378179206390117d+00
    b = 0.2743846121244060d+00
    v = 0.3581084321919782d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4836936460214534d+00
    b = 0.3895902610739024d-01
    v = 0.3426522117591512d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5293792562683797d+00
    b = 0.7871246819312640d-01
    v = 0.3491848770121379d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5726281253100033d+00
    b = 0.1187963808202981d+00
    v = 0.3539318235231476d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6133658776169068d+00
    b = 0.1587914708061787d+00
    v = 0.3570231438458694d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6515085491865307d+00
    b = 0.1983058575227646d+00
    v = 0.3586207335051714d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5778692716064976d+00
    b = 0.3977209689791542d-01
    v = 0.3541196205164025d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6207904288086192d+00
    b = 0.7990157592981152d-01
    v = 0.3574296911573953d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6608688171046802d+00
    b = 0.1199671308754309d+00
    v = 0.3591993279818963d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6656263089489130d+00
    b = 0.4015955957805969d-01
    v = 0.3595855034661997d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld3470(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 3470-POINT ANGULAR GRID
!
!      This routine is part of a set of routines that generate
!      Lebedev grids [1-6] for integration on a sphere. The original
!      C-code [1] was kindly provided by Dr. Dmitri N. Laikov and
!      translated into Fortran by Dr. Christoph van Wüllen.
!      This routine was translated using a C to Fortran77 conversion
!      tool written by Dr. Christoph van Wüllen.
!
!      Users of this code are asked to include reference [1] in their
!      publications, and in the user and programmer manuals
!      describing their codes.
!
!      This code was distributed through CCL (http://www.ccl.net/).
!
!      References:
!
!      [1] V.I. Lebedev and D.N. Laikov,
!          "A Quadrature Formula for the Sphere of the 131st
!          Algebraic Order of Accuracy,"
!          Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
!      [2] V.I. Lebedev,
!          "A Quadrature Formula for the Sphere of 59th Algebraic
!          Order of Accuracy,"
!          Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.
!
!      [3] V.I. Lebedev and A.L. Skorokhodov,
!          "Quadrature Formulas of Orders 41, 47, and 53 for the Sphere,"
!          Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.
!
!      [4] V.I. Lebedev,
!          "Spherical Quadrature Formulas Exact to Orders 25-29,"
!          Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.
!
!      [5] V.I. Lebedev,
!          "Quadratures on a Sphere,"
!          Computational Mathematics and Mathematical Physics, Vol. 16,
!          1976, pp. 10-24.
!
!      [6] V.I. Lebedev,
!          "Values of the Nodes and Weights of Ninth to Seventeenth
!          Order Gauss-Markov Quadrature Formulae Invariant under the
!          Octahedron Group with Inversion,"
!          Computational Mathematics and Mathematical Physics, Vol. 15,
!          1975, pp. 44-51.
!
    n = 1
    v = 0.2040382730826330d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3178149703889544d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1721420832906233d-01
    v = 0.8288115128076110d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4408875374981770d-01
    v = 0.1360883192522954d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7594680813878681d-01
    v = 0.1766854454542662d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1108335359204799d+00
    v = 0.2083153161230153d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1476517054388567d+00
    v = 0.2333279544657158d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1856731870860615d+00
    v = 0.2532809539930247d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2243634099428821d+00
    v = 0.2692472184211158d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2633006881662727d+00
    v = 0.2819949946811885d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3021340904916283d+00
    v = 0.2920953593973030d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3405594048030089d+00
    v = 0.2999889782948352d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3783044434007372d+00
    v = 0.3060292120496902d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4151194767407910d+00
    v = 0.3105109167522192d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4507705766443257d+00
    v = 0.3136902387550312d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4850346056573187d+00
    v = 0.3157984652454632d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5176950817792470d+00
    v = 0.3170516518425422d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5485384240820989d+00
    v = 0.3176568425633755d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6039117238943308d+00
    v = 0.3177198411207062d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6279956655573113d+00
    v = 0.3175519492394733d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6493636169568952d+00
    v = 0.3174654952634756d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6677644117704504d+00
    v = 0.3175676415467654d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6829368572115624d+00
    v = 0.3178923417835410d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6946195818184121d+00
    v = 0.3183788287531909d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7025711542057026d+00
    v = 0.3188755151918807d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7066004767140119d+00
    v = 0.3191916889313849d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5132537689946062d-01
    v = 0.1231779611744508d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1297994661331225d+00
    v = 0.1924661373839880d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2188852049401307d+00
    v = 0.2380881867403424d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3123174824903457d+00
    v = 0.2693100663037885d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4064037620738195d+00
    v = 0.2908673382834366d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4984958396944782d+00
    v = 0.3053914619381535d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5864975046021365d+00
    v = 0.3143916684147777d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6686711634580175d+00
    v = 0.3187042244055363d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8715738780835950d-01
    b = 0.2557175233367578d-01
    v = 0.1635219535869790d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1248383123134007d+00
    b = 0.5604823383376681d-01
    v = 0.1968109917696070d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1638062693383378d+00
    b = 0.8968568601900765d-01
    v = 0.2236754342249974d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2035586203373176d+00
    b = 0.1254086651976279d+00
    v = 0.2453186687017181d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2436798975293774d+00
    b = 0.1624780150162012d+00
    v = 0.2627551791580541d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2838207507773806d+00
    b = 0.2003422342683208d+00
    v = 0.2767654860152220d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3236787502217692d+00
    b = 0.2385628026255263d+00
    v = 0.2879467027765895d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3629849554840691d+00
    b = 0.2767731148783578d+00
    v = 0.2967639918918702d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4014948081992087d+00
    b = 0.3146542308245309d+00
    v = 0.3035900684660351d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4389818379260225d+00
    b = 0.3519196415895088d+00
    v = 0.3087338237298308d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4752331143674377d+00
    b = 0.3883050984023654d+00
    v = 0.3124608838860167d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5100457318374018d+00
    b = 0.4235613423908649d+00
    v = 0.3150084294226743d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5432238388954868d+00
    b = 0.4574484717196220d+00
    v = 0.3165958398598402d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5745758685072442d+00
    b = 0.4897311639255524d+00
    v = 0.3174320440957372d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1723981437592809d+00
    b = 0.3010630597881105d-01
    v = 0.2182188909812599d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2149553257844597d+00
    b = 0.6326031554204694d-01
    v = 0.2399727933921445d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2573256081247422d+00
    b = 0.9848566980258631d-01
    v = 0.2579796133514652d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2993163751238106d+00
    b = 0.1350835952384266d+00
    v = 0.2727114052623535d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3407238005148000d+00
    b = 0.1725184055442181d+00
    v = 0.2846327656281355d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3813454978483264d+00
    b = 0.2103559279730725d+00
    v = 0.2941491102051334d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4209848104423343d+00
    b = 0.2482278774554860d+00
    v = 0.3016049492136107d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4594519699996300d+00
    b = 0.2858099509982883d+00
    v = 0.3072949726175648d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4965640166185930d+00
    b = 0.3228075659915428d+00
    v = 0.3114768142886460d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5321441655571562d+00
    b = 0.3589459907204151d+00
    v = 0.3143823673666223d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5660208438582166d+00
    b = 0.3939630088864310d+00
    v = 0.3162269764661535d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5980264315964364d+00
    b = 0.4276029922949089d+00
    v = 0.3172164663759821d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2644215852350733d+00
    b = 0.3300939429072552d-01
    v = 0.2554575398967435d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3090113743443063d+00
    b = 0.6803887650078501d-01
    v = 0.2701704069135677d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3525871079197808d+00
    b = 0.1044326136206709d+00
    v = 0.2823693413468940d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3950418005354029d+00
    b = 0.1416751597517679d+00
    v = 0.2922898463214289d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4362475663430163d+00
    b = 0.1793408610504821d+00
    v = 0.3001829062162428d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4760661812145854d+00
    b = 0.2170630750175722d+00
    v = 0.3062890864542953d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5143551042512103d+00
    b = 0.2545145157815807d+00
    v = 0.3108328279264746d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5509709026935597d+00
    b = 0.2913940101706601d+00
    v = 0.3140243146201245d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5857711030329428d+00
    b = 0.3274169910910705d+00
    v = 0.3160638030977130d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6186149917404392d+00
    b = 0.3623081329317265d+00
    v = 0.3171462882206275d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3586894569557064d+00
    b = 0.3497354386450040d-01
    v = 0.2812388416031796d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4035266610019441d+00
    b = 0.7129736739757095d-01
    v = 0.2912137500288045d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4467775312332510d+00
    b = 0.1084758620193165d+00
    v = 0.2993241256502206d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4883638346608543d+00
    b = 0.1460915689241772d+00
    v = 0.3057101738983822d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5281908348434601d+00
    b = 0.1837790832369980d+00
    v = 0.3105319326251432d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5661542687149311d+00
    b = 0.2212075390874021d+00
    v = 0.3139565514428167d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6021450102031452d+00
    b = 0.2580682841160985d+00
    v = 0.3161543006806366d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6360520783610050d+00
    b = 0.2940656362094121d+00
    v = 0.3172985960613294d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4521611065087196d+00
    b = 0.3631055365867002d-01
    v = 0.2989400336901431d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4959365651560963d+00
    b = 0.7348318468484350d-01
    v = 0.3054555883947677d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5376815804038283d+00
    b = 0.1111087643812648d+00
    v = 0.3104764960807702d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5773314480243768d+00
    b = 0.1488226085145408d+00
    v = 0.3141015825977616d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6148113245575056d+00
    b = 0.1862892274135151d+00
    v = 0.3164520621159896d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6500407462842380d+00
    b = 0.2231909701714456d+00
    v = 0.3176652305912204d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5425151448707213d+00
    b = 0.3718201306118944d-01
    v = 0.3105097161023939d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5841860556907931d+00
    b = 0.7483616335067346d-01
    v = 0.3143014117890550d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6234632186851500d+00
    b = 0.1125990834266120d+00
    v = 0.3168172866287200d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6602934551848843d+00
    b = 0.1501303813157619d+00
    v = 0.3181401865570968d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6278573968375105d+00
    b = 0.3767559930245720d-01
    v = 0.3170663659156037d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6665611711264577d+00
    b = 0.7548443301360158d-01
    v = 0.3185447944625510d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld3890(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 3890-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.1807395252196920d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2848008782238827d-03
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2836065837530581d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1587876419858352d-01
    v = 0.7013149266673816d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4069193593751206d-01
    v = 0.1162798021956766d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7025888115257997d-01
    v = 0.1518728583972105d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1027495450028704d+00
    v = 0.1798796108216934d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1371457730893426d+00
    v = 0.2022593385972785d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1727758532671953d+00
    v = 0.2203093105575464d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2091492038929037d+00
    v = 0.2349294234299855d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2458813281751915d+00
    v = 0.2467682058747003d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2826545859450066d+00
    v = 0.2563092683572224d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3191957291799622d+00
    v = 0.2639253896763318d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3552621469299578d+00
    v = 0.2699137479265108d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3906329503406230d+00
    v = 0.2745196420166739d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4251028614093031d+00
    v = 0.2779529197397593d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4584777520111870d+00
    v = 0.2803996086684265d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4905711358710193d+00
    v = 0.2820302356715842d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5212011669847385d+00
    v = 0.2830056747491068d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5501878488737995d+00
    v = 0.2834808950776839d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6025037877479342d+00
    v = 0.2835282339078929d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6254572689549016d+00
    v = 0.2833819267065800d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6460107179528248d+00
    v = 0.2832858336906784d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6639541138154251d+00
    v = 0.2833268235451244d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6790688515667495d+00
    v = 0.2835432677029253d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6911338580371512d+00
    v = 0.2839091722743049d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6999385956126490d+00
    v = 0.2843308178875841d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7053037748656896d+00
    v = 0.2846703550533846d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4732224387180115d-01
    v = 0.1051193406971900d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1202100529326803d+00
    v = 0.1657871838796974d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2034304820664855d+00
    v = 0.2064648113714232d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2912285643573002d+00
    v = 0.2347942745819741d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3802361792726768d+00
    v = 0.2547775326597726d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4680598511056146d+00
    v = 0.2686876684847025d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5528151052155599d+00
    v = 0.2778665755515867d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6329386307803041d+00
    v = 0.2830996616782929d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8056516651369069d-01
    b = 0.2363454684003124d-01
    v = 0.1403063340168372d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1156476077139389d+00
    b = 0.5191291632545936d-01
    v = 0.1696504125939477d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1520473382760421d+00
    b = 0.8322715736994519d-01
    v = 0.1935787242745390d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1892986699745931d+00
    b = 0.1165855667993712d+00
    v = 0.2130614510521968d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2270194446777792d+00
    b = 0.1513077167409504d+00
    v = 0.2289381265931048d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2648908185093273d+00
    b = 0.1868882025807859d+00
    v = 0.2418630292816186d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3026389259574136d+00
    b = 0.2229277629776224d+00
    v = 0.2523400495631193d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3400220296151384d+00
    b = 0.2590951840746235d+00
    v = 0.2607623973449605d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3768217953335510d+00
    b = 0.2951047291750847d+00
    v = 0.2674441032689209d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4128372900921884d+00
    b = 0.3307019714169930d+00
    v = 0.2726432360343356d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4478807131815630d+00
    b = 0.3656544101087634d+00
    v = 0.2765787685924545d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4817742034089257d+00
    b = 0.3997448951939695d+00
    v = 0.2794428690642224d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5143472814653344d+00
    b = 0.4327667110812024d+00
    v = 0.2814099002062895d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5454346213905650d+00
    b = 0.4645196123532293d+00
    v = 0.2826429531578994d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5748739313170252d+00
    b = 0.4948063555703345d+00
    v = 0.2832983542550884d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1599598738286342d+00
    b = 0.2792357590048985d-01
    v = 0.1886695565284976d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1998097412500951d+00
    b = 0.5877141038139065d-01
    v = 0.2081867882748234d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2396228952566202d+00
    b = 0.9164573914691377d-01
    v = 0.2245148680600796d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2792228341097746d+00
    b = 0.1259049641962687d+00
    v = 0.2380370491511872d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3184251107546741d+00
    b = 0.1610594823400863d+00
    v = 0.2491398041852455d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3570481164426244d+00
    b = 0.1967151653460898d+00
    v = 0.2581632405881230d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3949164710492144d+00
    b = 0.2325404606175168d+00
    v = 0.2653965506227417d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4318617293970503d+00
    b = 0.2682461141151439d+00
    v = 0.2710857216747087d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4677221009931678d+00
    b = 0.3035720116011973d+00
    v = 0.2754434093903659d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5023417939270955d+00
    b = 0.3382781859197439d+00
    v = 0.2786579932519380d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5355701836636128d+00
    b = 0.3721383065625942d+00
    v = 0.2809011080679474d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5672608451328771d+00
    b = 0.4049346360466055d+00
    v = 0.2823336184560987d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5972704202540162d+00
    b = 0.4364538098633802d+00
    v = 0.2831101175806309d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2461687022333596d+00
    b = 0.3070423166833368d-01
    v = 0.2221679970354546d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2881774566286831d+00
    b = 0.6338034669281885d-01
    v = 0.2356185734270703d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3293963604116978d+00
    b = 0.9742862487067941d-01
    v = 0.2469228344805590d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3697303822241377d+00
    b = 0.1323799532282290d+00
    v = 0.2562726348642046d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4090663023135127d+00
    b = 0.1678497018129336d+00
    v = 0.2638756726753028d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4472819355411712d+00
    b = 0.2035095105326114d+00
    v = 0.2699311157390862d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4842513377231437d+00
    b = 0.2390692566672091d+00
    v = 0.2746233268403837d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5198477629962928d+00
    b = 0.2742649818076149d+00
    v = 0.2781225674454771d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5539453011883145d+00
    b = 0.3088503806580094d+00
    v = 0.2805881254045684d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5864196762401251d+00
    b = 0.3425904245906614d+00
    v = 0.2821719877004913d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6171484466668390d+00
    b = 0.3752562294789468d+00
    v = 0.2830222502333124d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3350337830565727d+00
    b = 0.3261589934634747d-01
    v = 0.2457995956744870d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3775773224758284d+00
    b = 0.6658438928081572d-01
    v = 0.2551474407503706d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4188155229848973d+00
    b = 0.1014565797157954d+00
    v = 0.2629065335195311d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4586805892009344d+00
    b = 0.1368573320843822d+00
    v = 0.2691900449925075d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4970895714224235d+00
    b = 0.1724614851951608d+00
    v = 0.2741275485754276d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5339505133960747d+00
    b = 0.2079779381416412d+00
    v = 0.2778530970122595d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5691665792531440d+00
    b = 0.2431385788322288d+00
    v = 0.2805010567646741d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6026387682680377d+00
    b = 0.2776901883049853d+00
    v = 0.2822055834031040d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6342676150163307d+00
    b = 0.3113881356386632d+00
    v = 0.2831016901243473d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4237951119537067d+00
    b = 0.3394877848664351d-01
    v = 0.2624474901131803d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4656918683234929d+00
    b = 0.6880219556291447d-01
    v = 0.2688034163039377d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5058857069185980d+00
    b = 0.1041946859721635d+00
    v = 0.2738932751287636d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5443204666713996d+00
    b = 0.1398039738736393d+00
    v = 0.2777944791242523d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5809298813759742d+00
    b = 0.1753373381196155d+00
    v = 0.2806011661660987d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6156416039447128d+00
    b = 0.2105215793514010d+00
    v = 0.2824181456597460d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6483801351066604d+00
    b = 0.2450953312157051d+00
    v = 0.2833585216577828d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5103616577251688d+00
    b = 0.3485560643800719d-01
    v = 0.2738165236962878d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5506738792580681d+00
    b = 0.7026308631512033d-01
    v = 0.2778365208203180d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5889573040995292d+00
    b = 0.1059035061296403d+00
    v = 0.2807852940418966d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6251641589516930d+00
    b = 0.1414823925236026d+00
    v = 0.2827245949674705d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6592414921570178d+00
    b = 0.1767207908214530d+00
    v = 0.2837342344829828d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5930314017533384d+00
    b = 0.3542189339561672d-01
    v = 0.2809233907610981d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6309812253390175d+00
    b = 0.7109574040369549d-01
    v = 0.2829930809742694d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6666296011353230d+00
    b = 0.1067259792282730d+00
    v = 0.2841097874111479d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6703715271049922d+00
    b = 0.3569455268820809d-01
    v = 0.2843455206008783d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld4334(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 4334-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.1449063022537883d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2546377329828424d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1462896151831013d-01
    v = 0.6018432961087496d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3769840812493139d-01
    v = 0.1002286583263673d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6524701904096891d-01
    v = 0.1315222931028093d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9560543416134648d-01
    v = 0.1564213746876724d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1278335898929198d+00
    v = 0.1765118841507736d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1613096104466031d+00
    v = 0.1928737099311080d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1955806225745371d+00
    v = 0.2062658534263270d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2302935218498028d+00
    v = 0.2172395445953787d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2651584344113027d+00
    v = 0.2262076188876047d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2999276825183209d+00
    v = 0.2334885699462397d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3343828669718798d+00
    v = 0.2393355273179203d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3683265013750518d+00
    v = 0.2439559200468863d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4015763206518108d+00
    v = 0.2475251866060002d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4339612026399770d+00
    v = 0.2501965558158773d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4653180651114582d+00
    v = 0.2521081407925925d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4954893331080803d+00
    v = 0.2533881002388081d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5243207068924930d+00
    v = 0.2541582900848261d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5516590479041704d+00
    v = 0.2545365737525860d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6012371927804176d+00
    v = 0.2545726993066799d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6231574466449819d+00
    v = 0.2544456197465555d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6429416514181271d+00
    v = 0.2543481596881064d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6604124272943595d+00
    v = 0.2543506451429194d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6753851470408250d+00
    v = 0.2544905675493763d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6876717970626160d+00
    v = 0.2547611407344429d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6970895061319234d+00
    v = 0.2551060375448869d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7034746912553310d+00
    v = 0.2554291933816039d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7067017217542295d+00
    v = 0.2556255710686343d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4382223501131123d-01
    v = 0.9041339695118195d-04
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1117474077400006d+00
    v = 0.1438426330079022d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1897153252911440d+00
    v = 0.1802523089820518d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2724023009910331d+00
    v = 0.2060052290565496d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3567163308709902d+00
    v = 0.2245002248967466d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4404784483028087d+00
    v = 0.2377059847731150d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5219833154161411d+00
    v = 0.2468118955882525d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5998179868977553d+00
    v = 0.2525410872966528d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6727803154548222d+00
    v = 0.2553101409933397d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7476563943166086d-01
    b = 0.2193168509461185d-01
    v = 0.1212879733668632d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1075341482001416d+00
    b = 0.4826419281533887d-01
    v = 0.1472872881270931d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1416344885203259d+00
    b = 0.7751191883575742d-01
    v = 0.1686846601010828d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1766325315388586d+00
    b = 0.1087558139247680d+00
    v = 0.1862698414660208d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2121744174481514d+00
    b = 0.1413661374253096d+00
    v = 0.2007430956991861d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2479669443408145d+00
    b = 0.1748768214258880d+00
    v = 0.2126568125394796d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2837600452294113d+00
    b = 0.2089216406612073d+00
    v = 0.2224394603372113d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3193344933193984d+00
    b = 0.2431987685545972d+00
    v = 0.2304264522673135d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3544935442438745d+00
    b = 0.2774497054377770d+00
    v = 0.2368854288424087d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3890571932288154d+00
    b = 0.3114460356156915d+00
    v = 0.2420352089461772d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4228581214259090d+00
    b = 0.3449806851913012d+00
    v = 0.2460597113081295d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4557387211304052d+00
    b = 0.3778618641248256d+00
    v = 0.2491181912257687d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4875487950541643d+00
    b = 0.4099086391698978d+00
    v = 0.2513528194205857d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5181436529962997d+00
    b = 0.4409474925853973d+00
    v = 0.2528943096693220d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5473824095600661d+00
    b = 0.4708094517711291d+00
    v = 0.2538660368488136d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5751263398976174d+00
    b = 0.4993275140354637d+00
    v = 0.2543868648299022d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1489515746840028d+00
    b = 0.2599381993267017d-01
    v = 0.1642595537825183d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1863656444351767d+00
    b = 0.5479286532462190d-01
    v = 0.1818246659849308d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2238602880356348d+00
    b = 0.8556763251425254d-01
    v = 0.1966565649492420d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2612723375728160d+00
    b = 0.1177257802267011d+00
    v = 0.2090677905657991d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2984332990206190d+00
    b = 0.1508168456192700d+00
    v = 0.2193820409510504d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3351786584663333d+00
    b = 0.1844801892177727d+00
    v = 0.2278870827661928d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3713505522209120d+00
    b = 0.2184145236087598d+00
    v = 0.2348283192282090d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4067981098954663d+00
    b = 0.2523590641486229d+00
    v = 0.2404139755581477d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4413769993687534d+00
    b = 0.2860812976901373d+00
    v = 0.2448227407760734d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4749487182516394d+00
    b = 0.3193686757808996d+00
    v = 0.2482110455592573d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5073798105075426d+00
    b = 0.3520226949547602d+00
    v = 0.2507192397774103d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5385410448878654d+00
    b = 0.3838544395667890d+00
    v = 0.2524765968534880d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5683065353670530d+00
    b = 0.4146810037640963d+00
    v = 0.2536052388539425d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5965527620663510d+00
    b = 0.4443224094681121d+00
    v = 0.2542230588033068d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2299227700856157d+00
    b = 0.2865757664057584d-01
    v = 0.1944817013047896d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2695752998553267d+00
    b = 0.5923421684485993d-01
    v = 0.2067862362746635d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3086178716611389d+00
    b = 0.9117817776057715d-01
    v = 0.2172440734649114d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3469649871659077d+00
    b = 0.1240593814082605d+00
    v = 0.2260125991723423d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3845153566319655d+00
    b = 0.1575272058259175d+00
    v = 0.2332655008689523d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4211600033403215d+00
    b = 0.1912845163525413d+00
    v = 0.2391699681532458d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4567867834329882d+00
    b = 0.2250710177858171d+00
    v = 0.2438801528273928d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4912829319232061d+00
    b = 0.2586521303440910d+00
    v = 0.2475370504260665d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5245364793303812d+00
    b = 0.2918112242865407d+00
    v = 0.2502707235640574d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5564369788915756d+00
    b = 0.3243439239067890d+00
    v = 0.2522031701054241d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5868757697775287d+00
    b = 0.3560536787835351d+00
    v = 0.2534511269978784d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6157458853519617d+00
    b = 0.3867480821242581d+00
    v = 0.2541284914955151d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3138461110672113d+00
    b = 0.3051374637507278d-01
    v = 0.2161509250688394d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3542495872050569d+00
    b = 0.6237111233730755d-01
    v = 0.2248778513437852d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3935751553120181d+00
    b = 0.9516223952401907d-01
    v = 0.2322388803404617d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4317634668111147d+00
    b = 0.1285467341508517d+00
    v = 0.2383265471001355d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4687413842250821d+00
    b = 0.1622318931656033d+00
    v = 0.2432476675019525d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5044274237060283d+00
    b = 0.1959581153836453d+00
    v = 0.2471122223750674d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5387354077925727d+00
    b = 0.2294888081183837d+00
    v = 0.2500291752486870d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5715768898356105d+00
    b = 0.2626031152713945d+00
    v = 0.2521055942764682d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6028627200136111d+00
    b = 0.2950904075286713d+00
    v = 0.2534472785575503d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6325039812653463d+00
    b = 0.3267458451113286d+00
    v = 0.2541599713080121d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3981986708423407d+00
    b = 0.3183291458749821d-01
    v = 0.2317380975862936d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4382791182133300d+00
    b = 0.6459548193880908d-01
    v = 0.2378550733719775d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4769233057218166d+00
    b = 0.9795757037087952d-01
    v = 0.2428884456739118d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5140823911194238d+00
    b = 0.1316307235126655d+00
    v = 0.2469002655757292d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5496977833862983d+00
    b = 0.1653556486358704d+00
    v = 0.2499657574265851d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5837047306512727d+00
    b = 0.1988931724126510d+00
    v = 0.2521676168486082d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6160349566926879d+00
    b = 0.2320174581438950d+00
    v = 0.2535935662645334d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6466185353209440d+00
    b = 0.2645106562168662d+00
    v = 0.2543356743363214d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4810835158795404d+00
    b = 0.3275917807743992d-01
    v = 0.2427353285201535d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5199925041324341d+00
    b = 0.6612546183967181d-01
    v = 0.2468258039744386d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5571717692207494d+00
    b = 0.9981498331474143d-01
    v = 0.2500060956440310d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5925789250836378d+00
    b = 0.1335687001410374d+00
    v = 0.2523238365420979d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6261658523859670d+00
    b = 0.1671444402896463d+00
    v = 0.2538399260252846d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6578811126669331d+00
    b = 0.2003106382156076d+00
    v = 0.2546255927268069d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5609624612998100d+00
    b = 0.3337500940231335d-01
    v = 0.2500583360048449d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5979959659984670d+00
    b = 0.6708750335901803d-01
    v = 0.2524777638260203d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6330523711054002d+00
    b = 0.1008792126424850d+00
    v = 0.2540951193860656d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6660960998103972d+00
    b = 0.1345050343171794d+00
    v = 0.2549524085027472d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6365384364585819d+00
    b = 0.3372799460737052d-01
    v = 0.2542569507009158d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6710994302899275d+00
    b = 0.6755249309678028d-01
    v = 0.2552114127580376d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld4802(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 4802-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.9687521879420705d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2307897895367918d-03
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2297310852498558d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2335728608887064d-01
    v = 0.7386265944001919d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4352987836550653d-01
    v = 0.8257977698542210d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6439200521088801d-01
    v = 0.9706044762057630d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9003943631993181d-01
    v = 0.1302393847117003d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1196706615548473d+00
    v = 0.1541957004600968d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1511715412838134d+00
    v = 0.1704459770092199d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1835982828503801d+00
    v = 0.1827374890942906d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2165081259155405d+00
    v = 0.1926360817436107d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2496208720417563d+00
    v = 0.2008010239494833d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2827200673567900d+00
    v = 0.2075635983209175d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3156190823994346d+00
    v = 0.2131306638690909d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3481476793749115d+00
    v = 0.2176562329937335d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3801466086947226d+00
    v = 0.2212682262991018d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4114652119634011d+00
    v = 0.2240799515668565d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4419598786519751d+00
    v = 0.2261959816187525d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4714925949329543d+00
    v = 0.2277156368808855d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4999293972879466d+00
    v = 0.2287351772128336d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5271387221431248d+00
    v = 0.2293490814084085d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5529896780837761d+00
    v = 0.2296505312376273d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6000856099481712d+00
    v = 0.2296793832318756d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6210562192785175d+00
    v = 0.2295785443842974d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6401165879934240d+00
    v = 0.2295017931529102d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6571144029244334d+00
    v = 0.2295059638184868d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6718910821718863d+00
    v = 0.2296232343237362d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6842845591099010d+00
    v = 0.2298530178740771d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6941353476269816d+00
    v = 0.2301579790280501d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7012965242212991d+00
    v = 0.2304690404996513d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7056471428242644d+00
    v = 0.2307027995907102d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4595557643585895d-01
    v = 0.9312274696671092d-04
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1049316742435023d+00
    v = 0.1199919385876926d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1773548879549274d+00
    v = 0.1598039138877690d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2559071411236127d+00
    v = 0.1822253763574900d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3358156837985898d+00
    v = 0.1988579593655040d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4155835743763893d+00
    v = 0.2112620102533307d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4937894296167472d+00
    v = 0.2201594887699007d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5691569694793316d+00
    v = 0.2261622590895036d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6405840854894251d+00
    v = 0.2296458453435705d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7345133894143348d-01
    b = 0.2177844081486067d-01
    v = 0.1006006990267000d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1009859834044931d+00
    b = 0.4590362185775188d-01
    v = 0.1227676689635876d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1324289619748758d+00
    b = 0.7255063095690877d-01
    v = 0.1467864280270117d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1654272109607127d+00
    b = 0.1017825451960684d+00
    v = 0.1644178912101232d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1990767186776461d+00
    b = 0.1325652320980364d+00
    v = 0.1777664890718961d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2330125945523278d+00
    b = 0.1642765374496765d+00
    v = 0.1884825664516690d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2670080611108287d+00
    b = 0.1965360374337889d+00
    v = 0.1973269246453848d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3008753376294316d+00
    b = 0.2290726770542238d+00
    v = 0.2046767775855328d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3344475596167860d+00
    b = 0.2616645495370823d+00
    v = 0.2107600125918040d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3675709724070786d+00
    b = 0.2941150728843141d+00
    v = 0.2157416362266829d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4001000887587812d+00
    b = 0.3262440400919066d+00
    v = 0.2197557816920721d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4318956350436028d+00
    b = 0.3578835350611916d+00
    v = 0.2229192611835437d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4628239056795531d+00
    b = 0.3888751854043678d+00
    v = 0.2253385110212775d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4927563229773636d+00
    b = 0.4190678003222840d+00
    v = 0.2271137107548774d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5215687136707969d+00
    b = 0.4483151836883852d+00
    v = 0.2283414092917525d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5491402346984905d+00
    b = 0.4764740676087880d+00
    v = 0.2291161673130077d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5753520160126075d+00
    b = 0.5034021310998277d+00
    v = 0.2295313908576598d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1388326356417754d+00
    b = 0.2435436510372806d-01
    v = 0.1438204721359031d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1743686900537244d+00
    b = 0.5118897057342652d-01
    v = 0.1607738025495257d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2099737037950268d+00
    b = 0.8014695048539634d-01
    v = 0.1741483853528379d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2454492590908548d+00
    b = 0.1105117874155699d+00
    v = 0.1851918467519151d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2807219257864278d+00
    b = 0.1417950531570966d+00
    v = 0.1944628638070613d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3156842271975842d+00
    b = 0.1736604945719597d+00
    v = 0.2022495446275152d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3502090945177752d+00
    b = 0.2058466324693981d+00
    v = 0.2087462382438514d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3841684849519686d+00
    b = 0.2381284261195919d+00
    v = 0.2141074754818308d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4174372367906016d+00
    b = 0.2703031270422569d+00
    v = 0.2184640913748162d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4498926465011892d+00
    b = 0.3021845683091309d+00
    v = 0.2219309165220329d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4814146229807701d+00
    b = 0.3335993355165720d+00
    v = 0.2246123118340624d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5118863625734701d+00
    b = 0.3643833735518232d+00
    v = 0.2266062766915125d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5411947455119144d+00
    b = 0.3943789541958179d+00
    v = 0.2280072952230796d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5692301500357246d+00
    b = 0.4234320144403542d+00
    v = 0.2289082025202583d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5958857204139576d+00
    b = 0.4513897947419260d+00
    v = 0.2294012695120025d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2156270284785766d+00
    b = 0.2681225755444491d-01
    v = 0.1722434488736947d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2532385054909710d+00
    b = 0.5557495747805614d-01
    v = 0.1830237421455091d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2902564617771537d+00
    b = 0.8569368062950249d-01
    v = 0.1923855349997633d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3266979823143256d+00
    b = 0.1167367450324135d+00
    v = 0.2004067861936271d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3625039627493614d+00
    b = 0.1483861994003304d+00
    v = 0.2071817297354263d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3975838937548699d+00
    b = 0.1803821503011405d+00
    v = 0.2128250834102103d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4318396099009774d+00
    b = 0.2124962965666424d+00
    v = 0.2174513719440102d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4651706555732742d+00
    b = 0.2445221837805913d+00
    v = 0.2211661839150214d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4974752649620969d+00
    b = 0.2762701224322987d+00
    v = 0.2240665257813102d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5286517579627517d+00
    b = 0.3075627775211328d+00
    v = 0.2262439516632620d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5586001195731895d+00
    b = 0.3382311089826877d+00
    v = 0.2277874557231869d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5872229902021319d+00
    b = 0.3681108834741399d+00
    v = 0.2287854314454994d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6144258616235123d+00
    b = 0.3970397446872839d+00
    v = 0.2293268499615575d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2951676508064861d+00
    b = 0.2867499538750441d-01
    v = 0.1912628201529828d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3335085485472725d+00
    b = 0.5867879341903510d-01
    v = 0.1992499672238701d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3709561760636381d+00
    b = 0.8961099205022284d-01
    v = 0.2061275533454027d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4074722861667498d+00
    b = 0.1211627927626297d+00
    v = 0.2119318215968572d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4429923648839117d+00
    b = 0.1530748903554898d+00
    v = 0.2167416581882652d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4774428052721736d+00
    b = 0.1851176436721877d+00
    v = 0.2206430730516600d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5107446539535904d+00
    b = 0.2170829107658179d+00
    v = 0.2237186938699523d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5428151370542935d+00
    b = 0.2487786689026271d+00
    v = 0.2260480075032884d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5735699292556964d+00
    b = 0.2800239952795016d+00
    v = 0.2277098884558542d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6029253794562866d+00
    b = 0.3106445702878119d+00
    v = 0.2287845715109671d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6307998987073145d+00
    b = 0.3404689500841194d+00
    v = 0.2293547268236294d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3752652273692719d+00
    b = 0.2997145098184479d-01
    v = 0.2056073839852528d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4135383879344028d+00
    b = 0.6086725898678011d-01
    v = 0.2114235865831876d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4506113885153907d+00
    b = 0.9238849548435643d-01
    v = 0.2163175629770551d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4864401554606072d+00
    b = 0.1242786603851851d+00
    v = 0.2203392158111650d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5209708076611709d+00
    b = 0.1563086731483386d+00
    v = 0.2235473176847839d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5541422135830122d+00
    b = 0.1882696509388506d+00
    v = 0.2260024141501235d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5858880915113817d+00
    b = 0.2199672979126059d+00
    v = 0.2277675929329182d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6161399390603444d+00
    b = 0.2512165482924867d+00
    v = 0.2289102112284834d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6448296482255090d+00
    b = 0.2818368701871888d+00
    v = 0.2295027954625118d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4544796274917948d+00
    b = 0.3088970405060312d-01
    v = 0.2161281589879992d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4919389072146628d+00
    b = 0.6240947677636835d-01
    v = 0.2201980477395102d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5279313026985183d+00
    b = 0.9430706144280313d-01
    v = 0.2234952066593166d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5624169925571135d+00
    b = 0.1263547818770374d+00
    v = 0.2260540098520838d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5953484627093287d+00
    b = 0.1583430788822594d+00
    v = 0.2279157981899988d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6266730715339185d+00
    b = 0.1900748462555988d+00
    v = 0.2291296918565571d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6563363204278871d+00
    b = 0.2213599519592567d+00
    v = 0.2297533752536649d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5314574716585696d+00
    b = 0.3152508811515374d-01
    v = 0.2234927356465995d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5674614932298185d+00
    b = 0.6343865291465561d-01
    v = 0.2261288012985219d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6017706004970264d+00
    b = 0.9551503504223951d-01
    v = 0.2280818160923688d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6343471270264178d+00
    b = 0.1275440099801196d+00
    v = 0.2293773295180159d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6651494599127802d+00
    b = 0.1593252037671960d+00
    v = 0.2300528767338634d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6050184986005704d+00
    b = 0.3192538338496105d-01
    v = 0.2281893855065666d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6390163550880400d+00
    b = 0.6402824353962306d-01
    v = 0.2295720444840727d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6711199107088448d+00
    b = 0.9609805077002909d-01
    v = 0.2303227649026753d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6741354429572275d+00
    b = 0.3211853196273233d-01
    v = 0.2304831913227114d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld5294(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 5294-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.9080510764308163d-04
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2084824361987793d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2303261686261450d-01
    v = 0.5011105657239616d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3757208620162394d-01
    v = 0.5942520409683854d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5821912033821852d-01
    v = 0.9564394826109721d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8403127529194872d-01
    v = 0.1185530657126338d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1122927798060578d+00
    v = 0.1364510114230331d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1420125319192987d+00
    v = 0.1505828825605415d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1726396437341978d+00
    v = 0.1619298749867023d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2038170058115696d+00
    v = 0.1712450504267789d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2352849892876508d+00
    v = 0.1789891098164999d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2668363354312461d+00
    v = 0.1854474955629795d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2982941279900452d+00
    v = 0.1908148636673661d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3295002922087076d+00
    v = 0.1952377405281833d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3603094918363593d+00
    v = 0.1988349254282232d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3905857895173920d+00
    v = 0.2017079807160050d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4202005758160837d+00
    v = 0.2039473082709094d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4490310061597227d+00
    v = 0.2056360279288953d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4769586160311491d+00
    v = 0.2068525823066865d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5038679887049750d+00
    v = 0.2076724877534488d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5296454286519961d+00
    v = 0.2081694278237885d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5541776207164850d+00
    v = 0.2084157631219326d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5990467321921213d+00
    v = 0.2084381531128593d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6191467096294587d+00
    v = 0.2083476277129307d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6375251212901849d+00
    v = 0.2082686194459732d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6540514381131168d+00
    v = 0.2082475686112415d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6685899064391510d+00
    v = 0.2083139860289915d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6810013009681648d+00
    v = 0.2084745561831237d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6911469578730340d+00
    v = 0.2087091313375890d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6988956915141736d+00
    v = 0.2089718413297697d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7041335794868720d+00
    v = 0.2092003303479793d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7067754398018567d+00
    v = 0.2093336148263241d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3840368707853623d-01
    v = 0.7591708117365267d-04
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9835485954117399d-01
    v = 0.1083383968169186d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1665774947612998d+00
    v = 0.1403019395292510d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2405702335362910d+00
    v = 0.1615970179286436d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3165270770189046d+00
    v = 0.1771144187504911d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3927386145645443d+00
    v = 0.1887760022988168d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4678825918374656d+00
    v = 0.1973474670768214d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5408022024266935d+00
    v = 0.2033787661234659d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6104967445752438d+00
    v = 0.2072343626517331d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6760910702685738d+00
    v = 0.2091177834226918d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6655644120217392d-01
    b = 0.1936508874588424d-01
    v = 0.9316684484675566d-04
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9446246161270182d-01
    b = 0.4252442002115869d-01
    v = 0.1116193688682976d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1242651925452509d+00
    b = 0.6806529315354374d-01
    v = 0.1298623551559414d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1553438064846751d+00
    b = 0.9560957491205369d-01
    v = 0.1450236832456426d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1871137110542670d+00
    b = 0.1245931657452888d+00
    v = 0.1572719958149914d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2192612628836257d+00
    b = 0.1545385828778978d+00
    v = 0.1673234785867195d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2515682807206955d+00
    b = 0.1851004249723368d+00
    v = 0.1756860118725188d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2838535866287290d+00
    b = 0.2160182608272384d+00
    v = 0.1826776290439367d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3159578817528521d+00
    b = 0.2470799012277111d+00
    v = 0.1885116347992865d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3477370882791392d+00
    b = 0.2781014208986402d+00
    v = 0.1933457860170574d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3790576960890540d+00
    b = 0.3089172523515731d+00
    v = 0.1973060671902064d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4097938317810200d+00
    b = 0.3393750055472244d+00
    v = 0.2004987099616311d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4398256572859637d+00
    b = 0.3693322470987730d+00
    v = 0.2030170909281499d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4690384114718480d+00
    b = 0.3986541005609877d+00
    v = 0.2049461460119080d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4973216048301053d+00
    b = 0.4272112491408562d+00
    v = 0.2063653565200186d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5245681526132446d+00
    b = 0.4548781735309936d+00
    v = 0.2073507927381027d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5506733911803888d+00
    b = 0.4815315355023251d+00
    v = 0.2079764593256122d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5755339829522475d+00
    b = 0.5070486445801855d+00
    v = 0.2083150534968778d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1305472386056362d+00
    b = 0.2284970375722366d-01
    v = 0.1262715121590664d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1637327908216477d+00
    b = 0.4812254338288384d-01
    v = 0.1414386128545972d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1972734634149637d+00
    b = 0.7531734457511935d-01
    v = 0.1538740401313898d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2308694653110130d+00
    b = 0.1039043639882017d+00
    v = 0.1642434942331432d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2643899218338160d+00
    b = 0.1334526587117626d+00
    v = 0.1729790609237496d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2977171599622171d+00
    b = 0.1636414868936382d+00
    v = 0.1803505190260828d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3307293903032310d+00
    b = 0.1942195406166568d+00
    v = 0.1865475350079657d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3633069198219073d+00
    b = 0.2249752879943753d+00
    v = 0.1917182669679069d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3953346955922727d+00
    b = 0.2557218821820032d+00
    v = 0.1959851709034382d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4267018394184914d+00
    b = 0.2862897925213193d+00
    v = 0.1994529548117882d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4573009622571704d+00
    b = 0.3165224536636518d+00
    v = 0.2022138911146548d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4870279559856109d+00
    b = 0.3462730221636496d+00
    v = 0.2043518024208592d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5157819581450322d+00
    b = 0.3754016870282835d+00
    v = 0.2059450313018110d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5434651666465393d+00
    b = 0.4037733784993613d+00
    v = 0.2070685715318472d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5699823887764627d+00
    b = 0.4312557784139123d+00
    v = 0.2077955310694373d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5952403350947741d+00
    b = 0.4577175367122110d+00
    v = 0.2081980387824712d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2025152599210369d+00
    b = 0.2520253617719557d-01
    v = 0.1521318610377956d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2381066653274425d+00
    b = 0.5223254506119000d-01
    v = 0.1622772720185755d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2732823383651612d+00
    b = 0.8060669688588620d-01
    v = 0.1710498139420709d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3080137692611118d+00
    b = 0.1099335754081255d+00
    v = 0.1785911149448736d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3422405614587601d+00
    b = 0.1399120955959857d+00
    v = 0.1850125313687736d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3758808773890420d+00
    b = 0.1702977801651705d+00
    v = 0.1904229703933298d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4088458383438932d+00
    b = 0.2008799256601680d+00
    v = 0.1949259956121987d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4410450550841152d+00
    b = 0.2314703052180836d+00
    v = 0.1986161545363960d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4723879420561312d+00
    b = 0.2618972111375892d+00
    v = 0.2015790585641370d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5027843561874343d+00
    b = 0.2920013195600270d+00
    v = 0.2038934198707418d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5321453674452458d+00
    b = 0.3216322555190551d+00
    v = 0.2056334060538251d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5603839113834030d+00
    b = 0.3506456615934198d+00
    v = 0.2068705959462289d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5874150706875146d+00
    b = 0.3789007181306267d+00
    v = 0.2076753906106002d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6131559381660038d+00
    b = 0.4062580170572782d+00
    v = 0.2081179391734803d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2778497016394506d+00
    b = 0.2696271276876226d-01
    v = 0.1700345216228943d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3143733562261912d+00
    b = 0.5523469316960465d-01
    v = 0.1774906779990410d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3501485810261827d+00
    b = 0.8445193201626464d-01
    v = 0.1839659377002642d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3851430322303653d+00
    b = 0.1143263119336083d+00
    v = 0.1894987462975169d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4193013979470415d+00
    b = 0.1446177898344475d+00
    v = 0.1941548809452595d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4525585960458567d+00
    b = 0.1751165438438091d+00
    v = 0.1980078427252384d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4848447779622947d+00
    b = 0.2056338306745660d+00
    v = 0.2011296284744488d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5160871208276894d+00
    b = 0.2359965487229226d+00
    v = 0.2035888456966776d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5462112185696926d+00
    b = 0.2660430223139146d+00
    v = 0.2054516325352142d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5751425068101757d+00
    b = 0.2956193664498032d+00
    v = 0.2067831033092635d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6028073872853596d+00
    b = 0.3245763905312779d+00
    v = 0.2076485320284876d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6291338275278409d+00
    b = 0.3527670026206972d+00
    v = 0.2081141439525255d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3541797528439391d+00
    b = 0.2823853479435550d-01
    v = 0.1834383015469222d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3908234972074657d+00
    b = 0.5741296374713106d-01
    v = 0.1889540591777677d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4264408450107590d+00
    b = 0.8724646633650199d-01
    v = 0.1936677023597375d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4609949666553286d+00
    b = 0.1175034422915616d+00
    v = 0.1976176495066504d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4944389496536006d+00
    b = 0.1479755652628428d+00
    v = 0.2008536004560983d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5267194884346086d+00
    b = 0.1784740659484352d+00
    v = 0.2034280351712291d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5577787810220990d+00
    b = 0.2088245700431244d+00
    v = 0.2053944466027758d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5875563763536670d+00
    b = 0.2388628136570763d+00
    v = 0.2068077642882360d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6159910016391269d+00
    b = 0.2684308928769185d+00
    v = 0.2077250949661599d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6430219602956268d+00
    b = 0.2973740761960252d+00
    v = 0.2082062440705320d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4300647036213646d+00
    b = 0.2916399920493977d-01
    v = 0.1934374486546626d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4661486308935531d+00
    b = 0.5898803024755659d-01
    v = 0.1974107010484300d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5009658555287261d+00
    b = 0.8924162698525409d-01
    v = 0.2007129290388658d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5344824270447704d+00
    b = 0.1197185199637321d+00
    v = 0.2033736947471293d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5666575997416371d+00
    b = 0.1502300756161382d+00
    v = 0.2054287125902493d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5974457471404752d+00
    b = 0.1806004191913564d+00
    v = 0.2069184936818894d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6267984444116886d+00
    b = 0.2106621764786252d+00
    v = 0.2078883689808782d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6546664713575417d+00
    b = 0.2402526932671914d+00
    v = 0.2083886366116359d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5042711004437253d+00
    b = 0.2982529203607657d-01
    v = 0.2006593275470817d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5392127456774380d+00
    b = 0.6008728062339922d-01
    v = 0.2033728426135397d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5726819437668618d+00
    b = 0.9058227674571398d-01
    v = 0.2055008781377608d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6046469254207278d+00
    b = 0.1211219235803400d+00
    v = 0.2070651783518502d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6350716157434952d+00
    b = 0.1515286404791580d+00
    v = 0.2080953335094320d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6639177679185454d+00
    b = 0.1816314681255552d+00
    v = 0.2086284998988521d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5757276040972253d+00
    b = 0.3026991752575440d-01
    v = 0.2055549387644668d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6090265823139755d+00
    b = 0.6078402297870770d-01
    v = 0.2071871850267654d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6406735344387661d+00
    b = 0.9135459984176636d-01
    v = 0.2082856600431965d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6706397927793709d+00
    b = 0.1218024155966590d+00
    v = 0.2088705858819358d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6435019674426665d+00
    b = 0.3052608357660639d-01
    v = 0.2083995867536322d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6747218676375681d+00
    b = 0.6112185773983089d-01
    v = 0.2090509712889637d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

  subroutine ld5810(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
!
!      LEBEDEV 5810-POINT ANGULAR GRID
!
!
!      THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
!      LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
!      C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
!      TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
!      THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
!      TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
!
!      USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
!      PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
!      DESCRIBING THEIR CODES.
!
!      THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
!
!      [1] V.I. LEBEDEV, AND D.N. LAIKOV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
!           ALGEBRAIC ORDER OF ACCURACY"
!          DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
!
!      [2] V.I. LEBEDEV
!          "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
!           ORDER OF ACCURACY"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
!
!      [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
!          "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
!          RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
!
!      [4] V.I. LEBEDEV
!          "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
!          SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
!
!      [5] V.I. LEBEDEV
!          "QUADRATURES ON A SPHERE"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
!          1976, PP. 10-24.
!
!      [6] V.I. LEBEDEV
!          "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
!           ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
!           OCTAHEDRON GROUP WITH INVERSION"
!          COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
!          1975, PP. 44-51.
!
    n = 1
    v = 0.9735347946175486d-05
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1907581241803167d-03
    call gen_oh(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1901059546737578d-03
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1182361662400277d-01
    v = 0.3926424538919212d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3062145009138958d-01
    v = 0.6667905467294382d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5329794036834243d-01
    v = 0.8868891315019135d-04
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7848165532862220d-01
    v = 0.1066306000958872d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1054038157636201d+00
    v = 0.1214506743336128d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1335577797766211d+00
    v = 0.1338054681640871d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1625769955502252d+00
    v = 0.1441677023628504d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1921787193412792d+00
    v = 0.1528880200826557d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2221340534690548d+00
    v = 0.1602330623773609d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2522504912791132d+00
    v = 0.1664102653445244d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2823610860679697d+00
    v = 0.1715845854011323d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3123173966267560d+00
    v = 0.1758901000133069d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3419847036953789d+00
    v = 0.1794382485256736d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3712386456999758d+00
    v = 0.1823238106757407d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3999627649876828d+00
    v = 0.1846293252959976d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4280466458648093d+00
    v = 0.1864284079323098d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4553844360185711d+00
    v = 0.1877882694626914d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4818736094437834d+00
    v = 0.1887716321852025d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5074138709260629d+00
    v = 0.1894381638175673d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5319061304570707d+00
    v = 0.1898454899533629d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5552514978677286d+00
    v = 0.1900497929577815d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5981009025246183d+00
    v = 0.1900671501924092d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6173990192228116d+00
    v = 0.1899837555533510d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6351365239411131d+00
    v = 0.1899014113156229d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6512010228227200d+00
    v = 0.1898581257705106d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6654758363948120d+00
    v = 0.1898804756095753d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6778410414853370d+00
    v = 0.1899793610426402d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6881760887484110d+00
    v = 0.1901464554844117d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6963645267094598d+00
    v = 0.1903533246259542d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7023010617153579d+00
    v = 0.1905556158463228d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.7059004636628753d+00
    v = 0.1907037155663528d-03
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3552470312472575d-01
    v = 0.5992997844249967d-04
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.9151176620841283d-01
    v = 0.9749059382456978d-04
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1566197930068980d+00
    v = 0.1241680804599158d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2265467599271907d+00
    v = 0.1437626154299360d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2988242318581361d+00
    v = 0.1584200054793902d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3717482419703886d+00
    v = 0.1694436550982744d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4440094491758889d+00
    v = 0.1776617014018108d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5145337096756642d+00
    v = 0.1836132434440077d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5824053672860230d+00
    v = 0.1876494727075983d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6468283961043370d+00
    v = 0.1899906535336482d-03
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6095964259104373d-01
    b = 0.1787828275342931d-01
    v = 0.8143252820767350d-04
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.8811962270959388d-01
    b = 0.3953888740792096d-01
    v = 0.9998859890887728d-04
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1165936722428831d+00
    b = 0.6378121797722990d-01
    v = 0.1156199403068359d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1460232857031785d+00
    b = 0.8985890813745037d-01
    v = 0.1287632092635513d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1761197110181755d+00
    b = 0.1172606510576162d+00
    v = 0.1398378643365139d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2066471190463718d+00
    b = 0.1456102876970995d+00
    v = 0.1491876468417391d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2374076026328152d+00
    b = 0.1746153823011775d+00
    v = 0.1570855679175456d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2682305474337051d+00
    b = 0.2040383070295584d+00
    v = 0.1637483948103775d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2989653312142369d+00
    b = 0.2336788634003698d+00
    v = 0.1693500566632843d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3294762752772209d+00
    b = 0.2633632752654219d+00
    v = 0.1740322769393633d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3596390887276086d+00
    b = 0.2929369098051601d+00
    v = 0.1779126637278296d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3893383046398812d+00
    b = 0.3222592785275512d+00
    v = 0.1810908108835412d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4184653789358347d+00
    b = 0.3512004791195743d+00
    v = 0.1836529132600190d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4469172319076166d+00
    b = 0.3796385677684537d+00
    v = 0.1856752841777379d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4745950813276976d+00
    b = 0.4074575378263879d+00
    v = 0.1872270566606832d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5014034601410262d+00
    b = 0.4345456906027828d+00
    v = 0.1883722645591307d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5272493404551239d+00
    b = 0.4607942515205134d+00
    v = 0.1891714324525297d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5520413051846366d+00
    b = 0.4860961284181720d+00
    v = 0.1896827480450146d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5756887237503077d+00
    b = 0.5103447395342790d+00
    v = 0.1899628417059528d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1225039430588352d+00
    b = 0.2136455922655793d-01
    v = 0.1123301829001669d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1539113217321372d+00
    b = 0.4520926166137188d-01
    v = 0.1253698826711277d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1856213098637712d+00
    b = 0.7086468177864818d-01
    v = 0.1366266117678531d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2174998728035131d+00
    b = 0.9785239488772918d-01
    v = 0.1462736856106918d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2494128336938330d+00
    b = 0.1258106396267210d+00
    v = 0.1545076466685412d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2812321562143480d+00
    b = 0.1544529125047001d+00
    v = 0.1615096280814007d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3128372276456111d+00
    b = 0.1835433512202753d+00
    v = 0.1674366639741759d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3441145160177973d+00
    b = 0.2128813258619585d+00
    v = 0.1724225002437900d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3749567714853510d+00
    b = 0.2422913734880829d+00
    v = 0.1765810822987288d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4052621732015610d+00
    b = 0.2716163748391453d+00
    v = 0.1800104126010751d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4349335453522385d+00
    b = 0.3007127671240280d+00
    v = 0.1827960437331284d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4638776641524965d+00
    b = 0.3294470677216479d+00
    v = 0.1850140300716308d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4920046410462687d+00
    b = 0.3576932543699155d+00
    v = 0.1867333507394938d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5192273554861704d+00
    b = 0.3853307059757764d+00
    v = 0.1880178688638289d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5454609081136522d+00
    b = 0.4122425044452694d+00
    v = 0.1889278925654758d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5706220661424140d+00
    b = 0.4383139587781027d+00
    v = 0.1895213832507346d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5946286755181518d+00
    b = 0.4634312536300553d+00
    v = 0.1898548277397420d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.1905370790924295d+00
    b = 0.2371311537781979d-01
    v = 0.1349105935937341d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2242518717748009d+00
    b = 0.4917878059254806d-01
    v = 0.1444060068369326d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2577190808025936d+00
    b = 0.7595498960495142d-01
    v = 0.1526797390930008d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2908724534927187d+00
    b = 0.1036991083191100d+00
    v = 0.1598208771406474d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3236354020056219d+00
    b = 0.1321348584450234d+00
    v = 0.1659354368615331d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3559267359304543d+00
    b = 0.1610316571314789d+00
    v = 0.1711279910946440d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3876637123676956d+00
    b = 0.1901912080395707d+00
    v = 0.1754952725601440d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4187636705218842d+00
    b = 0.2194384950137950d+00
    v = 0.1791247850802529d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4491449019883107d+00
    b = 0.2486155334763858d+00
    v = 0.1820954300877716d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4787270932425445d+00
    b = 0.2775768931812335d+00
    v = 0.1844788524548449d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5074315153055574d+00
    b = 0.3061863786591120d+00
    v = 0.1863409481706220d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5351810507738336d+00
    b = 0.3343144718152556d+00
    v = 0.1877433008795068d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5619001025975381d+00
    b = 0.3618362729028427d+00
    v = 0.1887444543705232d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5875144035268046d+00
    b = 0.3886297583620408d+00
    v = 0.1894009829375006d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6119507308734495d+00
    b = 0.4145742277792031d+00
    v = 0.1897683345035198d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2619733870119463d+00
    b = 0.2540047186389353d-01
    v = 0.1517327037467653d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.2968149743237949d+00
    b = 0.5208107018543989d-01
    v = 0.1587740557483543d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3310451504860488d+00
    b = 0.7971828470885599d-01
    v = 0.1649093382274097d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3646215567376676d+00
    b = 0.1080465999177927d+00
    v = 0.1701915216193265d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3974916785279360d+00
    b = 0.1368413849366629d+00
    v = 0.1746847753144065d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4295967403772029d+00
    b = 0.1659073184763559d+00
    v = 0.1784555512007570d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4608742854473447d+00
    b = 0.1950703730454614d+00
    v = 0.1815687562112174d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4912598858949903d+00
    b = 0.2241721144376724d+00
    v = 0.1840864370663302d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5206882758945558d+00
    b = 0.2530655255406489d+00
    v = 0.1860676785390006d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5490940914019819d+00
    b = 0.2816118409731066d+00
    v = 0.1875690583743703d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5764123302025542d+00
    b = 0.3096780504593238d+00
    v = 0.1886453236347225d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6025786004213506d+00
    b = 0.3371348366394987d+00
    v = 0.1893501123329645d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6275291964794956d+00
    b = 0.3638547827694396d+00
    v = 0.1897366184519868d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3348189479861771d+00
    b = 0.2664841935537443d-01
    v = 0.1643908815152736d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.3699515545855295d+00
    b = 0.5424000066843495d-01
    v = 0.1696300350907768d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4042003071474669d+00
    b = 0.8251992715430854d-01
    v = 0.1741553103844483d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4375320100182624d+00
    b = 0.1112695182483710d+00
    v = 0.1780015282386092d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4699054490335947d+00
    b = 0.1402964116467816d+00
    v = 0.1812116787077125d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5012739879431952d+00
    b = 0.1694275117584291d+00
    v = 0.1838323158085421d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5315874883754966d+00
    b = 0.1985038235312689d+00
    v = 0.1859113119837737d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5607937109622117d+00
    b = 0.2273765660020893d+00
    v = 0.1874969220221698d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5888393223495521d+00
    b = 0.2559041492849764d+00
    v = 0.1886375612681076d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6156705979160163d+00
    b = 0.2839497251976899d+00
    v = 0.1893819575809276d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6412338809078123d+00
    b = 0.3113791060500690d+00
    v = 0.1897794748256767d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4076051259257167d+00
    b = 0.2757792290858463d-01
    v = 0.1738963926584846d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4423788125791520d+00
    b = 0.5584136834984293d-01
    v = 0.1777442359873466d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4760480917328258d+00
    b = 0.8457772087727143d-01
    v = 0.1810010815068719d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5085838725946297d+00
    b = 0.1135975846359248d+00
    v = 0.1836920318248129d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5399513637391218d+00
    b = 0.1427286904765053d+00
    v = 0.1858489473214328d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5701118433636380d+00
    b = 0.1718112740057635d+00
    v = 0.1875079342496592d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5990240530606021d+00
    b = 0.2006944855985351d+00
    v = 0.1887080239102310d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6266452685139695d+00
    b = 0.2292335090598907d+00
    v = 0.1894905752176822d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6529320971415942d+00
    b = 0.2572871512353714d+00
    v = 0.1898991061200695d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.4791583834610126d+00
    b = 0.2826094197735932d-01
    v = 0.1809065016458791d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5130373952796940d+00
    b = 0.5699871359683649d-01
    v = 0.1836297121596799d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5456252429628476d+00
    b = 0.8602712528554394d-01
    v = 0.1858426916241869d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5768956329682385d+00
    b = 0.1151748137221281d+00
    v = 0.1875654101134641d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6068186944699046d+00
    b = 0.1442811654136362d+00
    v = 0.1888240751833503d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6353622248024907d+00
    b = 0.1731930321657680d+00
    v = 0.1896497383866979d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6624927035731797d+00
    b = 0.2017619958756061d+00
    v = 0.1900775530219121d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5484933508028488d+00
    b = 0.2874219755907391d-01
    v = 0.1858525041478814d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.5810207682142106d+00
    b = 0.5778312123713695d-01
    v = 0.1876248690077947d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6120955197181352d+00
    b = 0.8695262371439526d-01
    v = 0.1889404439064607d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6416944284294319d+00
    b = 0.1160893767057166d+00
    v = 0.1898168539265290d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6697926391731260d+00
    b = 0.1450378826743251d+00
    v = 0.1902779940661772d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6147594390585488d+00
    b = 0.2904957622341456d-01
    v = 0.1890125641731815d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6455390026356783d+00
    b = 0.5823809152617197d-01
    v = 0.1899434637795751d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6747258588365477d+00
    b = 0.8740384899884715d-01
    v = 0.1904520856831751d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    a = 0.6772135750395347d+00
    b = 0.2919946135808105d-01
    v = 0.1905534498734563d-03
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

!=====================================================================
!      Lebedev-Laikov grids end here
!=====================================================================

!
!      Octahedral symmetry (w/o inversion) 246-point angular grid
!      Order: 26
!
!      [1] A.S. POPOV
!          "THE SEARCH FOR THE SPHERE OF THE BEST CUBATURE FORMULAE
!           INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS"
!          SIB. ZH. VYCHISL. MAT., VOL. 5, NO. 4, 2002, PP. 367–372
!
  subroutine od0246(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.3897235138937088d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3493936461315880d-02
    a = 0.9566807992266356d+00
    b = 0.2833386636495198d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3591650554371342d-02
    a = 0.9757117308401445d+00
    b = 0.1167317879067011d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3967953538531415d-02
    a = 0.9001427110863452d+00
    b = 0.1004016871311435d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4044472616286787d-02
    a = 0.9016081884868434d+00
    b = 0.3046563962638665d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4158726826302798d-02
    a = 0.8657185254531113d+00
    b = 0.4789871418590519d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4223724523384502d-02
    a = 0.7958577819235757d+00
    b = 0.2817783845311169d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4287698313089224d-02
    a = 0.7828854398832512d+00
    b = 0.3632315411734175d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4295521563270202d-02
    a = 0.7797755454487808d+00
    b = 0.4876721073724419d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4306264958684786d-02
    a = 0.6482368665122326d+00
    b = 0.4557085192757556d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4322408526695458d-02
    a = 0.7237785649544816d+00
    b = 0.6543863740765325d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

!
!      Octahedral symmetry (w/o inversion) 264-point angular grid
!      Order: 27
!
!      [1] A.S. POPOV
!          "THE SEARCH FOR THE SPHERE OF THE BEST CUBATURE FORMULAE
!           INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS"
!          SIB. ZH. VYCHISL. MAT., VOL. 5, NO. 4, 2002, PP. 367–372
!
  subroutine od0264(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 3.404743138426950d-03
    a = 7.501162447877730d-01
    b = 6.569874882943270d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 3.527346280748710d-03
    a = 8.810738528375910d-01
    b = 2.804539315326130d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 3.560705552891180d-03
    a = 6.985522168362600d-01
    b = 6.638349164875560d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 3.661015502447230d-03
    a = 9.884416293559190d-01
    b = 4.187763670001430d-02
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 3.734364518407020d-03
    a = 8.703452942381050d-01
    b = 4.831216465940420d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 3.786218818475490d-03
    a = 7.577370457280040d-01
    b = 3.737914001474930d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 3.812295213144400d-03
    a = 9.340175312890680d-01
    b = 7.308888363903200d-02
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 3.922282615281700d-03
    a = 6.702365796250490d-01
    b = 5.696274199623510d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 3.934290436352150d-03
    a = 8.218581364908940d-01
    b = 4.788479450589300d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 4.143508784124530d-03
    a = 9.465732078368030d-01
    b = 2.750687576636780d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 4.179895806367260d-03
    a = 8.172239405554010d-01
    b = 1.335491891058270d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

!
!      Octahedral symmetry (w/o inversion) 342-point angular grid
!      Order: 31
!
!      [1] A.S. POPOV
!          "CUBATURE FORMULAE FOR THE SPHERE
!           INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS" (IN RUSSIAN)
!          COMPUT. MATH. MATH. PHYS., VOL. 38, NO. 1, 1998, PP. 30-37
!
  subroutine od0342(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.3023648748223408d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2640009601341421d-02
    a = 0.8596638928344739d+00
    b = 0.4062242558266217d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2743272994096360d-02
    a = 0.9048298313140403d+00
    b = 0.2370253708753523d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2892274286735084d-02
    a = 0.7919900954039140d+00
    b = 0.3635572640746994d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2898956168072910d-02
    a = 0.9246930024023610d+00
    b = 0.5409950346576234d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2922211185898755d-02
    a = 0.7983845764909109d+00
    b = 0.5645514094465199d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2926667482335111d-02
    a = 0.9629723259634715d+00
    b = 0.2038117136355577d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2931773085197849d-02
    a = 0.9120669247729113d+00
    b = 0.3877575423799307d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2936869297379136d-02
    a = 0.7377581016329626d+00
    b = 0.5411494384156677d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2975927626183587d-02
    a = 0.9810590688132057d+00
    b = 0.1311454620366076d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2986111941772548d-02
    a = 0.6522834517550465d+00
    b = 0.4833339354284559d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2987229421263409d-02
    a = 0.6952847727766393d+00
    b = 0.2919213249111464d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2995506191909761d-02
    a = 0.8278662627227689d+00
    b = 0.1738700791494905d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3035869254195839d-02
    a = 0.8316507582159955d+00
    b = 0.5548939979903705d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3038075943229046d-02
    a = 0.7126773308879356d+00
    b = 0.9613108329184178d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

!      Oh symmetry 350-point angular grid
!      Order: 31
!
!      [1] A.S. POPOV
!          "THE SEARCH FOR THE SPHERE OF THE BEST CUBATURE FORMULAE
!           INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS
!           WITH INVERSION FOR A SPHERE"
!          SIB. ZH. VYCHISL. MAT., VOL. 8, NO. 2, 2005, PP. 143–148
!
  subroutine ohd0350(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.3022655957073096d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3053589782677049d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1606436947071299d-02
    a = 0.7240037770684867d+00
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2894462580088764d-02
    a = 0.9249389183479164d+00
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2992303278280904d-02
    a = 0.9810676749667857d+00
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2705550334445959d-02
    a = 0.3622473127982810d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2966242835900409d-02
    a = 0.1917690066919518d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2999330613486812d-02
    a = 0.4801685990565060d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3041905028826839d-02
    a = 0.6502682026896775d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3044469577443018d-02
    a = 0.6935461975932684d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2846671108207593d-02
    a = 0.9077761244847474d+00
    b = 0.3765168541264181d+00
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2954001654360567d-02
    a = 0.7932087174843227d+00
    b = 0.5354121187940431d+00
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3020546347912858d-02
    a = 0.8290201309029455d+00
    b = 0.5505983319551926d+00
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)

    n = n-1
  end subroutine

!
!      Oh symmetry 398-point angular grid
!      Order: 33
!
!      [1] A.S. POPOV
!          "THE SEARCH FOR THE SPHERE OF THE BEST CUBATURE FORMULAE
!           INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS
!           WITH INVERSION FOR A SPHERE"
!          SIB. ZH. VYCHISL. MAT., VOL. 8, NO. 2, 2005, PP. 143–148
!
  subroutine ohd0398(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.1822084579093247d-02
    call gen_oh(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2827377032685751d-02
    call gen_oh(3, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.1575484799105965d-02
    a = 0.7062503174494366d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2580087158579662d-02
    a = 0.1133790435116091d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2639846092802613d-02
    a = 0.6901038654954956d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2783844769874123d-02
    a = 0.6457774767552608d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2846001949687294d-02
    a = 0.4796087844363489d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2873971674037779d-02
    a = 0.3632128775162115d+00
    call gen_oh(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2553168380704661d-02
    a = 0.9648952557917490d+00
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2654875161144892d-02
    a = 0.9014928604201154d+00
    call gen_oh(5, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2220113527157410d-02
    a = 0.8797615447057731d+00
    b = 0.4380336735755111d+00
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2326306559741858d-02
    a = 0.9388803094389767d+00
    b = 0.2946263848640117d+00
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2611549635094245d-02
    a = 0.8091315291511537d+00
    b = 0.5793558553216747d+00
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2722733540537044d-02
    a = 0.7838817431841564d+00
    b = 0.5439488406950955d+00
    call gen_oh(6, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

!
!      Octahedral symmetry (w/o inversion) 432-point angular grid
!      Order: 25
!
!      [1] A.S. POPOV
!          "CUBATURE FORMULAE FOR THE SPHERE
!           INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS" (IN RUSSIAN)
!          COMPUT. MATH. MATH. PHYS., VOL. 38, NO. 1, 1998, PP. 30-37
!
  subroutine od0432(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.2073915331320519d-02
    a = 0.7489806258088988d+00
    b = 0.6108701068802251d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2107255030397068d-02
    a = 0.9161730444214968d+00
    b = 0.3197412043673950d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2115342322975654d-02
    a = 0.9641217858844791d+00
    b = 0.1775069067845282d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2139217879967886d-02
    a = 0.9932625280244080d+00
    b = 0.4909487312548763d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2213455223136293d-02
    a = 0.9350867590143857d+00
    b = 0.1021948630931257d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2220017128107630d-02
    a = 0.8328350332317131d+00
    b = 0.5345156354592039d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2255553046127313d-02
    a = 0.9113314027451442d+00
    b = 0.4001239125939150d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2317129463567652d-02
    a = 0.8726339023393756d+00
    b = 0.5449493821656796d-01
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2319753441076136d-02
    a = 0.9678242462766119d+00
    b = 0.2485055252796990d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2346020122932169d-02
    a = 0.7078265977814910d+00
    b = 0.1321992811820306d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2365479437733548d-02
    a = 0.7070053872543563d+00
    b = 0.3277599108679579d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2373641673707936d-02
    a = 0.7815225804060471d+00
    b = 0.6238577434353162d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2434041150670229d-02
    a = 0.7200100968665680d+00
    b = 0.5549556437629704d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2458995407627726d-02
    a = 0.8068208852018608d+00
    b = 0.1984078314629552d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2465111275848897d-02
    a = 0.8385089184456621d+00
    b = 0.4494208832677718d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2466377894866924d-02
    a = 0.8801102003441889d+00
    b = 0.2566114675058704d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2496086422248489d-02
    a = 0.7879155638776240d+00
    b = 0.3902320295286575d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2499274414354599d-02
    a = 0.6504900372238729d+00
    b = 0.4941000883702588d+00
    call gen_oh(7, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

! Icosahedral symmetry (w/o inversion) 132-point angular grid
! Order: 19
!
! [1] Popov, A. S. (1994). Cubature formulae for a sphere invariant under
!     cyclic rotation groups.
!     Russian Journal of Numerical Analysis and Mathematical Modelling, 9(6), 535–546.
  subroutine id0132(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.6359381359381359d-02
    call gen_ih(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.7540532344007263d-02
    a = 0.8222245632603329d+00
    b = 0.4844155031786994d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.7854258050783132d-02
    a = 0.6921548451959329d+00
    b = 0.4181234628202419d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

! Icosahedral symmetry (w/o inversion) 152-point angular grid
! Order: 20
!
! [1] Popov, A. S. (1994). Cubature formulae for a sphere invariant under
!     cyclic rotation groups.
!     Russian Journal of Numerical Analysis and Mathematical Modelling, 9(6), 535–546.
  subroutine id0152(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.6611961667098409d-02
    call gen_ih(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4648223955885233d-02
    call gen_ih(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.6409543073731093d-02
    a = 0.7539155496229768d+00
    b = 0.4457443002522118d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.7385323274220814d-02
    a = 0.7228854451847260d+00
    b = 0.6440377245132042d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

! Icosahedral symmetry (w/o inversion) 180-point angular grid
! Order: 21
!
! [1] Popov, A. S. (1994). Cubature formulae for a sphere invariant under
!     cyclic rotation groups.
!     Russian Journal of Numerical Analysis and Mathematical Modelling, 9(6), 535–546.
  subroutine id0180(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.5018022982483986d-02
    a = 0.8243592074145533d+00
    b = 0.5307071420708551d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.5592501087262595d-02
    a = 0.7425058274010514d+00
    b = 0.2794547353533492d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.6056142596920085d-02
    a = 0.6952402305312859d+00
    b = 0.5565510052111179d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

! Icosahedral symmetry (w/o inversion) 192-point angular grid
! Order: 23
!
! [1] Popov, A. S. (1994). Cubature formulae for a sphere invariant under
!     cyclic rotation groups.
!     Russian Journal of Numerical Analysis and Mathematical Modelling, 9(6), 535–546.
  subroutine id0192(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.4164880157580243d-02
    call gen_ih(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.5077679306001083d-02
    a = 0.8032626540631647d+00
    b = 0.5443960115691368d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.5349546802216597d-02
    a = 0.6978513452732168d+00
    b = 0.5282061903582911d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.5406464526932938d-02
    a = 0.7323722630950026d+00
    b = 0.2934089823272757d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

! Icosahedral symmetry (w/o inversion) 212-point angular grid
! Order: 24
!
! [1] Popov, A. S. (1994). Cubature formulae for a sphere invariant under
!     cyclic rotation groups.
!     Russian Journal of Numerical Analysis and Mathematical Modelling, 9(6), 535–546.
  subroutine id0212(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.3741936943619837d-02
    call gen_ih(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4935771576102468d-02
    call gen_ih(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4551001693496304d-02
    a = 0.8452548782475415d+00
    b = 0.4824136177694623d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4855145722266675d-02
    a = 0.6779115110470655d+00
    b = 0.6479557053061163d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4866874670145564d-02
    a = 0.7753898899839309d+00
    b = 0.4269810983100433d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

! Icosahedral symmetry (w/o inversion) 242-point angular grid
! Order: 25
!
! [1] Popov, A. S. (1994). Cubature formulae for a sphere invariant under
!     cyclic rotation groups.
!     Russian Journal of Numerical Analysis and Mathematical Modelling, 9(6), 535–546.
  subroutine id0242(x, y, z, w, n)
    real(kind=dp) :: x(*)
    real(kind=dp) :: y(*)
    real(kind=dp) :: z(*)
    real(kind=dp) :: w(*)
    integer :: n
    real(kind=dp) :: a, b, v
    n = 1
    v = 0.3773233796543805d-02
    call gen_ih(1, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4742703016240548d-02
    call gen_ih(2, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.2350527849974790d-02
    call gen_ih(3, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.3976391306805966d-02
    a = 0.7588635353153991d+00
    b = 0.4631578933426286d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4548524050513189d-02
    a = 0.8418750816582262d+00
    b = 0.4876800613655622d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    v = 0.4630939619637840d-02
    a = 0.6718027690441337d+00
    b = 0.6591034254439872d+00
    call gen_ih(4, n, x(n), y(n), z(n), w(n), a, b, v)
    n = n-1
  end subroutine

end module
