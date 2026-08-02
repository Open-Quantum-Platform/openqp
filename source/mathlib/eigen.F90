module eigen
  use precision, only: dp
  use mathlib_types, only: blas_int
  use messages, only: show_message, WITH_ABORT
  use oqp_linalg
  use, intrinsic :: iso_c_binding, only: c_int64_t
  use blas_thread, only: blas_thread_count, blas_thread_set

  implicit none

  private
  public :: diag_symm_packed
  public :: diag_symm_full
  public :: schmd
  real(dp), parameter :: zero = 0.0_dp, two = 2.0_dp

!>  Rows of the matrix required before a dense eigensolve asks for one more
!>  BLAS thread.  See `eigen_blas_threads` for what this is buying.
  integer, parameter :: EIGEN_ROWS_PER_THREAD_DEFAULT = 512
contains

!>  @brief BLAS threads to use for a dense eigensolve of order `n`
!>
!>  LAPACK's symmetric eigensolvers do not thread the way GEMM does.  Inside
!>  OpenBLAS they are largely the netlib reference implementation, so they
!>  parallelise only through the BLAS calls they themselves make -- and the
!>  tridiagonal reduction that dominates DSYEVD is memory-bound level-2 work
!>  in many small calls.  Each of those calls forks and joins a full thread
!>  team, and below a certain size that coordination costs more than the
!>  arithmetic it is coordinating.
!>
!>  Measured on a 2 x 38-core Xeon 8368 against libopenblaso64 (OpenBLAS
!>  0.3.15, SkylakeX kernel), DSYEVD wall time in ms, best of several repeats:
!>
!>       n     1t      2t      4t      8t     16t     32t     64t    128t
!>     128   0.64    1.55    4.49    6.97    18.0    23.2    54.9   190.7
!>     256   3.11    5.10    8.00    12.4    37.7    49.4   105.8   509.6
!>     512   19.8    22.6    31.9    53.3   254.5   334.6   977.1  3630.8
!>    1024  131.4   122.7   152.1   319.0    1555    2092    5547   30527
!>    2048  935.9   741.8   583.3    1203    5497    7353   18552  102936
!>    4096   9279    5887    4302    4935   16308   21139   48123  259416
!>
!>  One thread is strictly best below n = 1024, nothing above four threads
!>  ever wins at any size tested, and the worst case is catastrophic: at
!>  n = 128, 128 threads is 299x slower than one.  For contrast, on the SAME
!>  library DGEMM scales properly -- 90.6 GFlop/s on one thread to 1775 on
!>  128 at n = 2048 -- so this is a property of the LAPACK layer, not of the
!>  BLAS, and GEMM-heavy paths are deliberately left alone.
!>
!>  The rule below is a single linear one rather than the table above, because
!>  the table is this machine and this OpenBLAS.  The mechanism generalises;
!>  the constants do not.  One thread per `EIGEN_ROWS_PER_THREAD` rows
!>  reproduces the measured optimum at every size except n = 4096, where it
!>  asks for 8 instead of 4 and gives up 15%.  It is also monotone and has no
!>  cliff, so on a BLAS that does thread these routines well it merely leaves
!>  some speed unclaimed for very large matrices -- where the absolute cost is
!>  already small next to the O(nbf^4) integral work around it.
!>
!>  Two properties matter for safety:
!>
!>    * it only ever LOWERS the count (`min` against what the caller had), so
!>      a caller that has already pinned BLAS -- `int2_twoei` does -- keeps
!>      its pin; and
!>    * it does nothing inside an OpenMP parallel region, where mutating a
!>      global thread count from several threads at once would be a race.
!>
!>  Override with OQP_EIGEN_ROWS_PER_THREAD; <= 0 disables the policy.
  function eigen_blas_threads(n, navail) result(nuse)
!$  use omp_lib, only: omp_in_parallel
    integer, intent(in) :: n
    integer(c_int64_t), intent(in) :: navail
    integer(c_int64_t) :: nuse
    character(32) :: env
    integer :: rows, stat_

    nuse = navail
    if (navail <= 1) return
!$  if (omp_in_parallel()) return

    rows = EIGEN_ROWS_PER_THREAD_DEFAULT
    call get_environment_variable('OQP_EIGEN_ROWS_PER_THREAD', env, status=stat_)
    if (stat_ == 0) then
      read(env, *, iostat=stat_) rows
      if (stat_ /= 0) rows = EIGEN_ROWS_PER_THREAD_DEFAULT
    end if
    if (rows <= 0) return          ! policy disabled

    nuse = max(1_c_int64_t, min(navail, int(n, c_int64_t) / int(rows, c_int64_t)))
  end function eigen_blas_threads

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
    integer(c_int64_t) :: nBlasSaved

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

!   Same reasoning as diag_symm_full: the packed drivers are, if anything,
!   less able to use a wide thread team than DSYEVD is.
    nBlasSaved = blas_thread_count()
    if (nBlasSaved > 0) call blas_thread_set(eigen_blas_threads(n, nBlasSaved))

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

    if (nBlasSaved > 0) call blas_thread_set(nBlasSaved)

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

    integer(blas_int) :: lda_, n_, info, lwork, liwork
    integer :: iok
    real(dp), dimension(:), allocatable :: work
    integer(blas_int), dimension(:), allocatable :: iwork
    real(dp) :: rwork(1)
    integer(blas_int) :: irwork(1)
    logical :: fatal
    character(16) :: driver
    integer(c_int64_t) :: nBlasSaved

    lda_    = int(lda, kind=blas_int)
    n_      = int(n, kind=blas_int)

    fatal = WITH_ABORT
    if (present(ierr)) fatal = WITHOUT_ABORT

!   Size the BLAS thread count to the eigenproblem; see eigen_blas_threads.
    nBlasSaved = blas_thread_count()
    if (nBlasSaved > 0) call blas_thread_set(eigen_blas_threads(n, nBlasSaved))

!   Divide-and-conquer driver: much faster than the QR-based DSYEV
!   for large matrices (it does most of its work in blocked level-3
!   BLAS), at the cost of a larger workspace. Same accuracy class.
    driver = 'DSYEVD'
    call dsyevd('V', 'U', n_, a, lda_, eival, rwork, -1_blas_int, &
                irwork, -1_blas_int, info)
    lwork = int(nint(rwork(1)), blas_int)
    liwork = irwork(1)
    allocate (work(lwork), iwork(liwork), stat=iok)
    if (iok /= 0) then
      if (nBlasSaved > 0) call blas_thread_set(nBlasSaved)
      if (present(ierr)) ierr = iok
      call show_message('Cannot allocate memory', fatal)
      return
    end if
    call dsyevd('V', 'U', n_, a, lda_, eival, work, lwork, iwork, liwork, info)

    if (nBlasSaved > 0) call blas_thread_set(nBlasSaved)

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
!   query both routines so that dorgqr can also run blocked
    call dgeqrf(n, m, v, ldv, x, wrksize, -1, info)
    lwork = max(int(wrksize(1)), n)
    call dorgqr(n, n, m, v, ldv, x, wrksize, -1, info)
    lwork = max(lwork, int(wrksize(1)))
    allocate (work(lwork))
    call dgeqrf(n, m, v, ldv, x, work, lwork, info)
    call dorgqr(n, n, m, v, ldv, x, work, lwork, info)
  end subroutine schmd
end module eigen
