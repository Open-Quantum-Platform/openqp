! Reference DGEMM and DGEMV for the pure-Fortran oracle tests.
!
! The oracles compile production routines that call BLAS with default
! (LP64) integers.  CI's Linux runner ships only an ILP64 OpenBLAS and no
! libblas at all, so linking a system BLAS is both unavailable and the wrong
! integer ABI.  These are straight transcriptions of the Netlib reference
! loops; speed is irrelevant at oracle sizes.

subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  implicit none
  character, intent(in) :: transa, transb
  integer, intent(in) :: m, n, k, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  double precision, intent(in) :: a(lda, *), b(ldb, *)
  double precision, intent(inout) :: c(ldc, *)
  logical :: nota, notb
  integer :: i, j, l
  double precision :: temp

  nota = transa == 'n' .or. transa == 'N'
  notb = transb == 'n' .or. transb == 'N'
  if (m == 0 .or. n == 0) return
  do j = 1, n
    if (beta == 0.0d0) then
      c(1:m, j) = 0.0d0
    else if (beta /= 1.0d0) then
      c(1:m, j) = beta*c(1:m, j)
    end if
    if (alpha == 0.0d0) cycle
    do l = 1, k
      if (notb) then
        temp = alpha*b(l, j)
      else
        temp = alpha*b(j, l)
      end if
      if (nota) then
        do i = 1, m
          c(i, j) = c(i, j) + temp*a(i, l)
        end do
      else
        do i = 1, m
          c(i, j) = c(i, j) + temp*a(l, i)
        end do
      end if
    end do
  end do
end subroutine dgemm

subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
  implicit none
  character, intent(in) :: trans
  integer, intent(in) :: m, n, lda, incx, incy
  double precision, intent(in) :: alpha, beta
  double precision, intent(in) :: a(lda, *), x(*)
  double precision, intent(inout) :: y(*)
  logical :: notr
  integer :: i, j, lenx, leny, ix, iy, jx, jy, kx, ky
  double precision :: temp

  notr = trans == 'n' .or. trans == 'N'
  if (m == 0 .or. n == 0) return
  if (notr) then
    lenx = n
    leny = m
  else
    lenx = m
    leny = n
  end if
  kx = 1
  if (incx < 0) kx = 1 - (lenx - 1)*incx
  ky = 1
  if (incy < 0) ky = 1 - (leny - 1)*incy
  iy = ky
  do i = 1, leny
    if (beta == 0.0d0) then
      y(iy) = 0.0d0
    else if (beta /= 1.0d0) then
      y(iy) = beta*y(iy)
    end if
    iy = iy + incy
  end do
  if (alpha == 0.0d0) return
  if (notr) then
    jx = kx
    do j = 1, n
      temp = alpha*x(jx)
      iy = ky
      do i = 1, m
        y(iy) = y(iy) + temp*a(i, j)
        iy = iy + incy
      end do
      jx = jx + incx
    end do
  else
    jy = ky
    do j = 1, n
      temp = 0.0d0
      ix = kx
      do i = 1, m
        temp = temp + a(i, j)*x(ix)
        ix = ix + incx
      end do
      y(jy) = y(jy) + alpha*temp
      jy = jy + incy
    end do
  end if
end subroutine dgemv
