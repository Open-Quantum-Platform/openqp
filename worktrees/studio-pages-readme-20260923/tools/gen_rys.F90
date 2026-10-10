module rys_prepare

  implicit none
  public

  integer, parameter :: dp = 8

  real(kind=dp), parameter :: pi      = 4.0_dp * atan(1.0_dp)
  real(kind=dp), parameter :: sqrt_pi = sqrt(pi)
  real(kind=dp), parameter :: pio4    = pi / 4.0_dp
  real(kind=dp), parameter :: eps_rp = 1.0e-14_dp

  integer, parameter :: eri_mxrys = 16, maux = 55
  integer, parameter :: nauxs(eri_mxrys) = [ &
    20, 25, 30, 30, 35, &
    40, 40, 40, 45, 50, &
    50, 55, 55, 0, 0, 0], &

    maprys(eri_mxrys) = &
    [1, 2, 3, 3, 4, 5, 5, 5, &
     6, 7, 7, 8, 8, 0, 0, 0]

  real(kind=dp), parameter :: xasymp(eri_mxrys) = [ &
    29.0_dp, 37.0_dp, 43.0_dp, 49.0_dp, &
    55.0_dp, 60.0_dp, 65.0_dp, 71.0_dp, &
    76.0_dp, 81.0_dp, 86.0_dp, 91.0_dp, &
    96.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

  real(kind=dp):: rtsasy(eri_mxrys, eri_mxrys), wtsasy(eri_mxrys, eri_mxrys)
  real(kind=dp) :: rtsaux(maux,8), wtsaux(maux,8)

contains

!-----------------------------------------------------------------

  subroutine egr_rysgw(n, alpha, beta, eps, roots, weight, ierr, wrk)

    integer, intent(IN) :: &
      n

    real(KIND=dp), dimension(n), intent(IN) :: &
      alpha, beta

    real(KIND=dp), dimension(n), intent(OUT) :: &
      roots, weight

    real(KIND=dp), dimension(n), intent(INOUT) :: &
      wrk

    integer, intent(OUT) :: &
      ierr

    real(KIND=dp), intent(IN) :: &
      eps

    integer :: &
      i, j, k, l, m, mml, ii

    real(KIND=dp) :: &
      b, c, f, g, p, r, s

    if (n < 1) then
      ierr = -1
      return
    end if

    ierr = 0
    roots(1:n) = alpha(1:n)
    do k = 1, n
      if (beta(k) < 0.0e+00_dp) then
        ierr = -2
        return
      end if
    end do
    if (n == 1) then
      weight(1) = beta(1)
      return
    end if

    wrk(1:n-1) = sqrt(beta(2:n))
    !CALL vdsqrt(n-1,beta(2:n),wrk(1:n-1))
    wrk(n) = 0.0e+00_dp
    weight(1) = 1.0e+00_dp
    weight(2:n) = 0.0e+00_dp

    large: do l = 1, n
      j = 0

! look for a small subdiagonal element.

      iter: do
        do m = l, n
          if (m == n) exit
          if (abs(wrk(m)) <= &
              eps*(abs(roots(m))+abs(roots(m+1)))) exit
        end do
        p = roots(l)
        if (m == l) exit iter
        if (j == 30) then
          ierr = l
          return
        end if
        j = j+1

! form shift.

        g = (roots(l+1)-p)/(2.0e+00_dp*wrk(l))
        r = sqrt(g*g+1.0e+00_dp)
        g = roots(m)-p+wrk(l)/(g+sign(r, g))
        s = 1.0e+00_dp
        c = 1.0e+00_dp
        p = 0.0e+00_dp
        mml = m-l

! for i=m-1 step -1 until l do ...

        do ii = 1, mml
          i = m-ii
          f = s*wrk(i)
          b = c*wrk(i)
          if (abs(f) < abs(g)) then
            s = f/g
            r = sqrt(s*s+1.0e+00_dp)
            wrk(i+1) = g*r
            c = 1.0e+00_dp/r
            s = s*c
          else
            c = g/f
            r = sqrt(c*c+1.0e+00_dp)
            wrk(i+1) = f*r
            s = 1.0e+00_dp/r
            c = c*s
          end if
          g = roots(i+1)-p
          r = (roots(i)-g)*s+2.0e+00_dp*c*b
          p = s*r
          roots(i+1) = g+p
          g = c*r-b

! form first component of vector.

          f = weight(i+1)
          weight(i+1) = s*weight(i)+c*f
          weight(i) = c*weight(i)-s*f
        end do
        roots(l) = roots(l)-p
        wrk(l) = g
        wrk(m) = 0.0e+00_dp
      end do iter
    end do large

    weight(1:n) = beta(1)*weight(1:n)*weight(1:n)

  end subroutine egr_rysgw

!-----------------------------------------------------------------

  subroutine egr_setrys

    real(KIND=dp), dimension(:) :: &
      rts(maux), wts(maux), wrk(maux, maux), &
      alpha(0:maux-1), beta(0:maux-1)

    integer :: &
      i, j, m, n, ierr, nauxsv, naux, igrid
!
!     ----- initialize the rys quadrature procedure -----
!
!     numerical parameters were determined by experimentation:
!     a) the -xx- value at which we can go to the asymptotic formula
!        was found by matching roughly 13 digits in roots/weights.
!     b) density of the auxiliary grid for the discretization was
!        determined by eyeball, at the -xx- value where we can begin
!        to use the asymptotic formula.  for a given -nroots-, the
!        required grid density gradually increases with -xx-.
!
!
!     set up asymptotic root computation, by generating
!     the roots of the hermite polynomials of order 2n.
!
!
    do i = 1, eri_mxrys
!            note that maux must be at least twice eri_mxrys, from next line
      n = 2*i

      alpha(0:n-1) = 0.0d0
      beta(0) = sqrt_pi
      do j = 1, n-1
        beta(j) = j/2.0d0
      end do
      call egr_rysgw(n, alpha, beta, eps_rp, rts, wts, ierr, wrk)
      call egr_sort(n, rts, wts)

      rtsasy(1:i, i) = rts(i+1:i+i)*rts(i+1:i+i)
      wtsasy(1:i, i) = wts(i+1:i+i)
    end do
!
!        generate auxiliary grids, at 8 distinct point densities
!        the auxiliary quadrature is taken to be "shifted legendre"
!
    nauxsv = 0
    igrid = 0
    do m = 1, eri_mxrys
      naux = nauxs(m)
      if (naux == nauxsv) cycle
      if (naux == 0) exit
      igrid = igrid+1
      nauxsv = naux

      alpha(0:naux-1) = 0.5d0

      beta(0) = 1.0d0
      do i = 1, naux-1
        beta(i) = 0.25d0/(4.0d0-1.0d0/(i*i))
      end do

      call egr_rysgw(naux, alpha(0:naux-1), beta(0:naux-1), &
                     eps_rp, rts(1:naux), wts(1:naux), ierr, wrk)

      rtsaux(1:naux, igrid) = rts(1:naux)**2
      wtsaux(1:naux, igrid) = wts(1:naux)

    end do

  end subroutine egr_setrys

  subroutine egr_sort(n, roots, weights)
    implicit none
    integer :: n
    real(kind=dp), intent(inout) :: roots(*), weights(*)
    integer :: ii, i, j, k
    real(kind=dp) :: p

    do ii = 2, n
      i = ii-1
      k = i
      p = roots(i)
      do j = ii, n
        if (roots(j) < p) then
          k = j
          p = roots(j)
        end if
      end do
      if (k /= i) then
        roots(k) = roots(i)
        roots(i) = p
        p = weights(i)
        weights(i) = weights(k)
        weights(k) = p
      end if
    end do
  end subroutine egr_sort

end module rys_prepare

program gen_rys
  use rys_prepare
  character(*), parameter :: int_name = 'integer, public, parameter'
  character(*), parameter :: real_name = 'real(kind=8), public, parameter'
  character(*), parameter :: int_fmt = '(2XA/ *( 4x, 20(I0,:,", "), "&", /) )'
  character(*), parameter :: real_fmt = '(2XA/ *( 4x, 3(D23.16,:,", "), "&", /) )'

  call egr_setrys

  write(*,'(G0)') 'module rys_lut'
  write(*,'(2XA)') 'implicit none'
  write(*,'(2XA)') 'private'

  write(*,'(2XA,I0)') int_name//' :: eri_mxrys  = ', eri_mxrys
  write(*,'(2XA,I0)') int_name//' :: maux = ', maux

  write(*,advance='no',fmt=int_fmt) int_name//' :: nauxs(eri_mxrys) = [ &', nauxs
  write(*,*) ']'

  write(*,advance='no',fmt=int_fmt) int_name//' :: maprys(eri_mxrys) = [ &', maprys
  write(*,*) ']'

  write(*,advance='no',fmt=real_fmt) real_name//' :: xasymp(eri_mxrys) = [ &', xasymp
  write(*,*) ']'

  write(*,advance='no',fmt=real_fmt) real_name//' :: rtsasy(eri_mxrys,eri_mxrys) = reshape([ &', rtsasy
  write(*,*) '], shape=[eri_mxrys, eri_mxrys])'

  write(*,advance='no',fmt=real_fmt) real_name//' :: wtsasy(eri_mxrys,eri_mxrys) = reshape([ &', wtsasy
  write(*,*) '], shape=[eri_mxrys, eri_mxrys])'

  write(*,advance='no',fmt=real_fmt) real_name//' :: rtsaux(maux,8) = reshape([ &', rtsaux
  write(*,*) '], shape=[maux, 8])'

  write(*,advance='no',fmt=real_fmt) real_name//' :: wtsaux(maux,8) = reshape([ &', wtsaux
  write(*,*) '], shape=[maux, 8])'

  write(*,'(G0)') 'end module rys_lut'
end program gen_rys
