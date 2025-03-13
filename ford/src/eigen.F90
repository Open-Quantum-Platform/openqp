module eigen
  use precision, only: dp
  use mathlib_types, only: blas_int
  use messages, only: show_message, WITH_ABORT
  use oqp_linalg

  implicit none

  private
  public :: diag_symm_packed
  public :: diag_symm_full
  public :: schmd
  real(dp), parameter :: zero = 0.0_dp, two = 2.0_dp
contains

!>  @brief Find eigenvalues and eigenvectors of symmetric matrix
!>         in packed format
!>  @param[in]     mode   algorithm of diagonalization (not used now)
!>  @param[in]     n      matrix dimension
!>  @param[in]     ldvect leading dimension of eigenvector matrix
!>  @param[in]     nvect  required number of eigenvectors
!>  @param[in,out] h      matrix to be diagonalized
!>  @param[out]    eig    eigenvalues
!>  @param[out]    vector eigenvectors
!>  @param[out]    ierr   status
  subroutine diag_symm_packed(mode, ldvect, nvect, n, h, eig, vector, ierr)
    use messages, only: show_message, WITH_ABORT, WITHOUT_ABORT
!
    integer, intent(in) :: mode
    integer, intent(in) :: ldvect, nvect, n
    integer, optional, intent(out) :: ierr
    real(dp), intent(inout) :: h(*)
    real(kind=dp), intent(out) :: eig(*), vector(*)

    integer(blas_int), dimension(:), allocatable :: iwork, ifail
    integer(blas_int) :: ldvect_, n_, nvect_, info, ione
    integer :: iok
    real(dp), dimension(:), allocatable :: work
    real(dp) :: abstol, dlamch
    logical :: fatal

    character(16) :: driver

    ldvect_ = int(ldvect, kind=blas_int)
    n_      = int(n, kind=blas_int)
    nvect_  = int(nvect, kind=blas_int)
    ione    = 1

    fatal = WITH_ABORT
    if (present(ierr)) fatal = WITHOUT_ABORT

    allocate (work(n*8), iwork(n*5), ifail(n), stat=iok)
    if (iok /= 0) then
      if (present(ierr)) ierr = iok
      call show_message('Cannot allocate memory', fatal)
      return
    end if

    if (nvect == n .and. ldvect >= n) then
      driver = 'DSPEV'
      call dspev('V', 'U', n_, h, eig, vector, ldvect_, work, info)
    else
      abstol = two*DLAMCH('S')
      driver = 'DSPEVX'
      call dspevx('V', 'A', 'U', &
            ldvect_, h, zero, zero, ione, ione, abstol, n_, &
            eig, vector, nvect_, work, iwork, ifail, info)
    end if

    if (present(ierr)) ierr = info

    if (info /= 0) then
      call show_message('(A,I0)', &
              trim(driver)//' FAILED! INFO: ', int(info), fatal)
    end if

  end subroutine diag_symm_packed

!>  @brief Find eigenvalues and eigenvectors of symmetric matrix
!>         in full format
!>  @param[in]     mode   algorithm of diagonalization (not used now)
!>  @param[in]     n      matrix dimension
!>  @param[in,out] a      matrix to be diagonalized, overwritten by
!>                        the eigenvectors on the exit
!>  @param[in]     lda    leading dimension of the matrix
!>  @param[out]    eig    eigenvalues
!>  @param[out]    ierr   status
  subroutine diag_symm_full(mode, n, a, lda, eival, ierr)
    use messages, only: show_message, WITH_ABORT, WITHOUT_ABORT
!
    integer, intent(in) :: mode
    integer, intent(in) :: n, lda
    real(dp), intent(inout) :: a(*)
    real(kind=dp), intent(out) :: eival(*)
    integer, optional, intent(out) :: ierr

    integer(blas_int) :: lda_, n_, info, ione, lwork
    integer :: iok
    real(dp), dimension(:), allocatable :: work
    real(dp) :: rwork(1)
    logical :: fatal
    character(16) :: driver

    lda_    = int(lda, kind=blas_int)
    n_      = int(n, kind=blas_int)
    ione    = 1

    fatal = WITH_ABORT
    if (present(ierr)) fatal = WITHOUT_ABORT

    driver = 'DSYEV'
    call dsyev('V', 'U', n_, a, lda_, eival, rwork, -1_blas_int, info)
    lwork = int(nint(rwork(1)), blas_int)
    allocate (work(lwork), stat=iok)
    if (iok /= 0) then
      if (present(ierr)) ierr = iok
      call show_message('Cannot allocate memory', fatal)
      return
    end if
    call dsyev('V', 'U', n_, a, lda_, eival, work, lwork, info)


    if (present(ierr)) ierr = info

    if (info /= 0) then
      call show_message('(A,I0)', &
              trim(driver)//' FAILED! INFO: ', int(info), fatal)
    end if

  end subroutine diag_symm_full

  subroutine schmd(v, m, n, ldv, x)
    use, intrinsic :: iso_fortran_env, only: real64
    use messages, only: show_message, WITH_ABORT
    implicit none

    integer, intent(IN) :: ldv, m, n
    real(real64), intent(INOUT) :: v(ldv, n), x(n)

    real(real64), allocatable :: work(:)
    integer :: lwork
    integer :: info
    real(real64) :: wrksize(1)

    if (M > N) then
      call show_message("SCHMD: M > N", WITH_ABORT)
    end if
    if (N > LDV) then
      call show_message("SCHMD: N > LDV", WITH_ABORT)
    end if
!   Householder QR-based version using LAPACK:
    call dgeqrf(n, m, v, ldv, x, wrksize, -1, info)
    lwork = max(int(wrksize(1)), n)
    allocate (work(lwork))
    call dgeqrf(n, m, v, ldv, x, work, lwork, info)
    call dorgqr(n, n, m, v, ldv, x, work, lwork, info)
  end subroutine schmd
end module eigen
