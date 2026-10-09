!> Numerical kernel for the fractional-derivative (Caputo) ingredient xi^alpha.
!>
!> K_n(x) = 1/Gamma(c) * int_0^1 (1-t)^(c-1) t^n exp(-x t^2) dt,   c = p - alpha > 0,
!> x = zeta*|r - R_i|^2 >= 0, n = 0..nmax (Cartesian power sum of the primitive, +2 for the
!> derivative branch).  Two regimes: Gauss-Jacobi on [0,1] with the (1-t)^(c-1) weight for
!> x <= T^2, and Gauss-Legendre on [0, T/sqrt(x)] (Gaussian-peak region, weight regular) for
!> x > T^2.  Validated to <1e-11 relative against the exact 2F2 form
!>   Gamma(n+1)/Gamma(n+1+c) 2F2((n+1)/2,(n+2)/2; (n+1+c)/2,(n+2+c)/2; -x)
!> for 0 < c <= 2 (see sessions/20261009_fracderiv_functional/kernel_proto.py).
module xi_alpha_kernel
  use precision, only: fp
  implicit none
  private
  public :: xi_kernel_t

  integer, parameter :: NJ_DEFAULT = 48, NL_DEFAULT = 32
  real(fp), parameter :: T_SPLIT = 7.0_fp

  type :: xi_kernel_t
    real(fp) :: alpha = 1.0_fp
    integer  :: p = 1          !< integer order of the inner derivative (0 or 1)
    real(fp) :: c = 0.0_fp     !< p - alpha
    real(fp) :: x0 = T_SPLIT**2
    integer  :: nj = 0, nl = 0
    real(fp), allocatable :: tj(:), wj(:)   !< Gauss-Jacobi nodes on [0,1], weights incl. 1/Gamma(c)
    real(fp), allocatable :: sl(:), wl(:)   !< Gauss-Legendre nodes/weights on [0,1]
    logical :: ready = .false.
  contains
    procedure :: init => xi_kernel_init
    procedure :: eval => xi_kernel_eval
  end type xi_kernel_t

contains

  !> Prepare quadrature rules for a given fractional order alpha (<1).
  subroutine xi_kernel_init(this, alpha, p, nj, nl)
    class(xi_kernel_t), intent(inout) :: this
    real(fp), intent(in) :: alpha
    integer, intent(in), optional :: p, nj, nl
    integer :: i
    real(fp) :: lg

    this%alpha = alpha
    if (present(p)) then
      this%p = p
    else
      this%p = max(0, ceiling(alpha))
    end if
    this%c = real(this%p, fp) - alpha
    if (this%c < 0.0_fp) error stop 'xi_alpha_kernel: p - alpha must be non-negative'
    this%ready = .true.
    if (this%c == 0.0_fp) return   ! integer order: K_n(x) = exp(-x), no quadrature needed
    this%nj = NJ_DEFAULT; if (present(nj)) this%nj = nj
    this%nl = NL_DEFAULT; if (present(nl)) this%nl = nl
    if (allocated(this%tj)) deallocate(this%tj, this%wj, this%sl, this%wl)
    allocate(this%tj(this%nj), this%wj(this%nj), this%sl(this%nl), this%wl(this%nl))

    call gauss_jacobi(this%nj, this%c - 1.0_fp, 0.0_fp, this%tj, this%wj)
    call gauss_legendre(this%nl, this%sl, this%wl)
    ! map [-1,1] -> [0,1]:  (1-t)^(c-1) dt = 2^(-c) (1-s)^(c-1) ds
    lg = log_gamma(this%c)
    do i = 1, this%nj
      this%tj(i) = 0.5_fp*(this%tj(i) + 1.0_fp)
      this%wj(i) = this%wj(i) * exp(-this%c*log(2.0_fp) - lg)
    end do
    do i = 1, this%nl
      this%sl(i) = 0.5_fp*(this%sl(i) + 1.0_fp)
      this%wl(i) = 0.5_fp*this%wl(i)
    end do
    this%ready = .true.
  end subroutine xi_kernel_init

  !> K_n(x) for n = 0..nmax.
  subroutine xi_kernel_eval(this, nmax, x, k)
    class(xi_kernel_t), intent(in) :: this
    integer, intent(in) :: nmax
    real(fp), intent(in) :: x
    real(fp), intent(out) :: k(0:nmax)
    integer :: i, n
    real(fp) :: tc, t, f, tn, lg

    k = 0.0_fp
    if (this%c == 0.0_fp) then
      k = exp(-x)
      return
    end if
    if (x <= this%x0) then
      do i = 1, this%nj
        t = this%tj(i)
        f = this%wj(i) * exp(-x*t*t)
        tn = f
        do n = 0, nmax
          k(n) = k(n) + tn
          tn = tn * t
        end do
      end do
    else
      tc = T_SPLIT / sqrt(x)
      lg = log_gamma(this%c)
      do i = 1, this%nl
        t = this%sl(i) * tc
        f = this%wl(i) * tc * exp((this%c - 1.0_fp)*log(1.0_fp - t) - x*t*t - lg)
        tn = f
        do n = 0, nmax
          k(n) = k(n) + tn
          tn = tn * t
        end do
      end do
    end if
  end subroutine xi_kernel_eval

  !> Gauss-Jacobi nodes/weights on [-1,1] for weight (1-s)^a (1+s)^b by Golub-Welsch:
  !> eigen-decomposition of the symmetric tridiagonal Jacobi matrix (implicit QL, no LAPACK).
  subroutine gauss_jacobi(n, a, b, x, w)
    integer, intent(in) :: n
    real(fp), intent(in) :: a, b
    real(fp), intent(out) :: x(n), w(n)
    real(fp) :: d(n), e(n), z(n, n), ab, k, mu0
    integer :: i

    ab = a + b
    do i = 1, n
      k = real(i - 1, fp)
      if (i == 1) then
        d(i) = (b - a)/(ab + 2.0_fp)
      else
        d(i) = (b*b - a*a)/((2.0_fp*k + ab)*(2.0_fp*k + ab + 2.0_fp))
      end if
    end do
    e(1) = 0.0_fp
    do i = 2, n
      k = real(i - 1, fp)
      if (i == 2) then
        e(i) = sqrt(4.0_fp*(1.0_fp + a)*(1.0_fp + b)/((ab + 2.0_fp)**2*(ab + 3.0_fp)))
      else
        e(i) = sqrt(4.0_fp*k*(k + a)*(k + b)*(k + ab)/ &
                    ((2.0_fp*k + ab)**2*(2.0_fp*k + ab + 1.0_fp)*(2.0_fp*k + ab - 1.0_fp)))
      end if
    end do
    z = 0.0_fp
    do i = 1, n
      z(i, i) = 1.0_fp
    end do
    call tridiag_ql(n, d, e, z)
    mu0 = exp((ab + 1.0_fp)*log(2.0_fp) + log_gamma(a + 1.0_fp) + log_gamma(b + 1.0_fp) - log_gamma(ab + 2.0_fp))
    do i = 1, n
      x(i) = d(i)
      w(i) = mu0*z(1, i)**2
    end do
  end subroutine gauss_jacobi

  !> Implicit QL with eigenvectors for a symmetric tridiagonal matrix
  !> (diagonal d, sub-diagonal e(2:n)); on exit d holds ascending eigenvalues
  !> and the columns of z the eigenvectors.
  subroutine tridiag_ql(n, d, e, z)
    integer, intent(in) :: n
    real(fp), intent(inout) :: d(n), e(n), z(n, n)
    integer :: i, iter, k, l, m
    real(fp) :: b, c, dd, f, g, p, r, s, t

    do i = 2, n
      e(i - 1) = e(i)
    end do
    e(n) = 0.0_fp
    do l = 1, n
      iter = 0
      do
        do m = l, n - 1
          dd = abs(d(m)) + abs(d(m + 1))
          if (abs(e(m)) <= epsilon(1.0_fp)*dd) exit
        end do
        if (m == l) exit
        iter = iter + 1
        if (iter > 60) error stop 'xi_alpha_kernel: tridiag_ql did not converge'
        g = (d(l + 1) - d(l))/(2.0_fp*e(l))
        r = hypot(g, 1.0_fp)
        g = d(m) - d(l) + e(l)/(g + sign(r, g))
        s = 1.0_fp; c = 1.0_fp; p = 0.0_fp
        do i = m - 1, l, -1
          f = s*e(i); b = c*e(i)
          r = hypot(f, g)
          e(i + 1) = r
          if (r == 0.0_fp) then
            d(i + 1) = d(i + 1) - p
            e(m) = 0.0_fp
            exit
          end if
          s = f/r; c = g/r
          g = d(i + 1) - p
          r = (d(i) - g)*s + 2.0_fp*c*b
          p = s*r
          d(i + 1) = g + p
          g = c*r - b
          do k = 1, n
            f = z(k, i + 1)
            z(k, i + 1) = s*z(k, i) + c*f
            z(k, i) = c*z(k, i) - s*f
          end do
        end do
        if (r == 0.0_fp .and. i >= l) cycle
        d(l) = d(l) - p
        e(l) = g
        e(m) = 0.0_fp
      end do
    end do
    ! sort ascending
    do i = 1, n - 1
      k = i; p = d(i)
      do m = i + 1, n
        if (d(m) < p) then
          k = m; p = d(m)
        end if
      end do
      if (k /= i) then
        d(k) = d(i); d(i) = p
        do m = 1, n
          t = z(m, i); z(m, i) = z(m, k); z(m, k) = t
        end do
      end if
    end do
  end subroutine tridiag_ql

  !> Gauss-Legendre nodes/weights on [-1,1].
  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    real(fp), intent(out) :: x(n), w(n)
    real(fp), parameter :: pi = 3.14159265358979323846_fp, eps = 1.0e-15_fp
    integer :: i, j, m, its
    real(fp) :: z, z1, p1, p2, p3, pp
    m = (n + 1)/2
    do i = 1, m
      z = cos(pi*(i - 0.25_fp)/(n + 0.5_fp))
      do its = 1, 100
        p1 = 1.0_fp; p2 = 0.0_fp
        do j = 1, n
          p3 = p2; p2 = p1
          p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j
        end do
        pp = n*(z*p1 - p2)/(z*z - 1.0_fp)
        z1 = z; z = z1 - p1/pp
        if (abs(z - z1) <= eps) exit
      end do
      x(i) = -z; x(n + 1 - i) = z
      w(i) = 2.0_fp/((1.0_fp - z*z)*pp*pp); w(n + 1 - i) = w(i)
    end do
  end subroutine gauss_legendre

end module xi_alpha_kernel
