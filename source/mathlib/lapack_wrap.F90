module lapack_wrap

  use precision, only: dp
  use mathlib_types, only: BLAS_INT, HUGE_BLAS_INT
  use messages, only: show_message, WITH_ABORT

  implicit none

  public
  private dp, blas_int, huge_blas_int, show_message, with_abort

  logical, parameter, private :: ARG_CHECK = .false.
  character(len=*), parameter, private :: BITNESS(2) = ["32", "64"]
  character(len=*), parameter, private :: ERRMSG = &
          "Integer is too big for " // BITNESS(BLAS_INT/4) // "bit BLAS/LAPACK"

contains

!----------------------------------------------------------------------

  subroutine oqp_dgeqrf_i64(m, n, a, lda, tau, work, lwork, info)
    integer :: info, lda, lwork, m, n
    real(dp) ::  a(lda,*), tau(*), work(*)

    integer(blas_int) :: info_, lda_, lwork_, m_, n_

    logical :: ok

    if (ARG_CHECK) then
    ok = .true.

    ok = ok .and. (lda   <= HUGE_BLAS_INT)
    ok = ok .and. (lwork <= HUGE_BLAS_INT)
    ok = ok .and. (m     <= HUGE_BLAS_INT)
    ok = ok .and. (n     <= HUGE_BLAS_INT)

    if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    lda_   = int(lda,   blas_int)
    lwork_ = int(lwork, blas_int)
    m_     = int(m,     blas_int)
    n_     = int(n,     blas_int)

    call dgeqrf(m_, n_, a, lda_, tau, work, lwork_, info_)

    info  = info_

  end subroutine

!----------------------------------------------------------------------

  subroutine oqp_dgesv_i64(n, nrhs, a, lda, ipiv, b, ldb, info)
    integer :: info, lda, ldb, n, nrhs
    integer :: ipiv(*)
    real(dp) ::  a(lda,*), b(ldb, *)

    integer(blas_int) :: info_, lda_, ldb_, n_, nrhs_
    integer(blas_int), allocatable :: ipiv_(:)

    logical :: ok

    if (ARG_CHECK) then
    ok = .true.

    ok = ok .and. (lda   <= HUGE_BLAS_INT)
    ok = ok .and. (ldb   <= HUGE_BLAS_INT)
    ok = ok .and. (n     <= HUGE_BLAS_INT)
    ok = ok .and. (nrhs  <= HUGE_BLAS_INT)

    if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    lda_   = int(lda,   blas_int)
    ldb_   = int(ldb,   blas_int)
    n_     = int(n,     blas_int)
    nrhs_  = int(nrhs,  blas_int)
    allocate(ipiv_(n))

    call dgesv(n_, nrhs_, a, lda_, ipiv, b, ldb_, info_)


    info  = info_
    ipiv(:n) = ipiv_

  end subroutine

!----------------------------------------------------------------------

  subroutine oqp_dsysv_i64(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    character(*) :: uplo
    integer :: info, lda, ldb, lwork, n, nrhs
    integer :: ipiv(*)
    real(dp) ::  a(*), b(*), work(*)

    integer(blas_int) :: info_, lda_, ldb_, lwork_, n_, nrhs_
    integer(blas_int), allocatable :: ipiv_(:)

    logical :: ok

    if (ARG_CHECK) then
    ok = .true.

    ok = ok .and. (lda   <= HUGE_BLAS_INT)
    ok = ok .and. (ldb   <= HUGE_BLAS_INT)
    ok = ok .and. (lwork <= HUGE_BLAS_INT)
    ok = ok .and. (n     <= HUGE_BLAS_INT)
    ok = ok .and. (nrhs  <= HUGE_BLAS_INT)

    if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    lda_   = int(lda,   blas_int)
    ldb_   = int(ldb,   blas_int)
    lwork_ = int(lwork, blas_int)
    n_     = int(n,     blas_int)
    nrhs_  = int(nrhs,  blas_int)
    if (lwork_ /= -1) allocate(ipiv_(n_))

    call dsysv(uplo, n_, nrhs_, a, lda_, ipiv_, b, ldb_, work, lwork_, info_)

    info  = info_
    if (lwork_ /= -1) ipiv(:n_) = ipiv_(:n_)

  end subroutine
!----------------------------------------------------------------------

  subroutine oqp_dgglse_i64(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    integer :: info, lda, ldb, lwork, m, n, p
    real(dp) ::  a(lda,*), b(ldb, *), c(*), d(*), work(*), x(*)

    integer(blas_int) :: info_, lda_, ldb_, lwork_, m_, n_, p_

    logical :: ok

    if (ARG_CHECK) then
    ok = .true.

    ok = ok .and. (info  <= HUGE_BLAS_INT)
    ok = ok .and. (lda   <= HUGE_BLAS_INT)
    ok = ok .and. (ldb   <= HUGE_BLAS_INT)
    ok = ok .and. (lwork <= HUGE_BLAS_INT)
    ok = ok .and. (m     <= HUGE_BLAS_INT)
    ok = ok .and. (n     <= HUGE_BLAS_INT)
    ok = ok .and. (p     <= HUGE_BLAS_INT)

    if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    lda_   = int(lda,   blas_int)
    ldb_   = int(ldb,   blas_int)
    lwork_ = int(lwork, blas_int)
    m_     = int(m,     blas_int)
    n_     = int(n,     blas_int)
    p_     = int(p,     blas_int)

    call dgglse(m_, n_, p_, a, lda_, b, ldb_, c, d, x, work, lwork_, info_)

    info  = info_

  end subroutine

!----------------------------------------------------------------------

  subroutine oqp_dorgqr_i64(m, n, k, a, lda, tau, work, lwork, info)
    integer :: info, k, lda, lwork, m, n
    real(dp) ::  a(lda,*), tau(*), work(*)

    integer(blas_int) :: info_, k_, lda_, lwork_, m_, n_

    logical :: ok

    if (ARG_CHECK) then
    ok = .true.

    ok = ok .and. (lda   <= HUGE_BLAS_INT)
    ok = ok .and. (lwork <= HUGE_BLAS_INT)
    ok = ok .and. (m     <= HUGE_BLAS_INT)
    ok = ok .and. (n     <= HUGE_BLAS_INT)
    ok = ok .and. (k     <= HUGE_BLAS_INT)

    if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    lda_   = int(lda,   blas_int)
    lwork_ = int(lwork, blas_int)
    m_     = int(m,     blas_int)
    n_     = int(n,     blas_int)
    k_     = int(k,     blas_int)

    call dorgqr(m_, n_, k_, a, lda_, tau, work, lwork_, info_)

    info  = info_

  end subroutine

!----------------------------------------------------------------------

  subroutine oqp_dormqr_i64(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    character(*) :: side, trans
    integer :: info, k, lda, ldc, lwork, m, n
    real(dp) ::  a(lda,*), c(ldc,*), tau(*), work(*)

    integer(blas_int) :: info_, k_, lda_, ldc_, lwork_, m_, n_

    logical :: ok

    if (ARG_CHECK) then
    ok = .true.

    ok = ok .and. (lda   <= HUGE_BLAS_INT)
    ok = ok .and. (ldc   <= HUGE_BLAS_INT)
    ok = ok .and. (lwork <= HUGE_BLAS_INT)
    ok = ok .and. (m     <= HUGE_BLAS_INT)
    ok = ok .and. (n     <= HUGE_BLAS_INT)
    ok = ok .and. (k     <= HUGE_BLAS_INT)

    if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    lda_   = int(lda,   blas_int)
    ldc_   = int(ldc,   blas_int)
    lwork_ = int(lwork, blas_int)
    m_     = int(m,     blas_int)
    n_     = int(n,     blas_int)
    k_     = int(k,     blas_int)

    call dormqr(side, trans, m_, n_, k_, a, lda_, tau, c, ldc_, work, lwork_, info_)

    info  = info_

  end subroutine

!----------------------------------------------------------------------

  subroutine oqp_dtpttr_i64(uplo, n, ap, a, lda, info)
    character(len=*) :: uplo
    integer :: info, n, lda
    real(dp) ::  a(lda,*), ap(*)

    integer(blas_int) :: info_, n_, lda_

    logical :: ok

    if (ARG_CHECK) then
    ok = .true.

    ok = ok .and. (n   <= HUGE_BLAS_INT)
    ok = ok .and. (lda <= HUGE_BLAS_INT)

    if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n,   blas_int)
    lda_ = int(lda, blas_int)

    call dtpttr(uplo, n_, ap, a, lda_, info_)

    info = info_

  end subroutine

!----------------------------------------------------------------------

  subroutine oqp_dtrttp_i64(uplo, n, a, lda, ap, info)
    character(len=*) :: uplo
    integer :: info, n, lda
    real(dp) ::  a(lda,*), ap(*)

    integer(blas_int) :: info_, n_, lda_

    logical :: ok

    if (ARG_CHECK) then
    ok = .true.

    ok = ok .and. (n   <= HUGE_BLAS_INT)
    ok = ok .and. (lda <= HUGE_BLAS_INT)

    if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n,   blas_int)
    lda_ = int(lda, blas_int)

    call dtrttp(uplo, n_, a, lda_, ap, info_)

    info = info_

  end subroutine

!----------------------------------------------------------------------

!> Compute the inverse of a real matrix using LAPACK routine dgetri.
!! @param n       Size of the matrix.
!! @param a       On entry, the matrix to be inverted. On exit, the inverted matrix.
!! @param lda     Leading dimension of the array `a`.
!! @param ipiv    Integer array of size `n` containing pivot indices computed by dgetrf.
!! @param work    Workspace array of size `lwork`.
!! @param lwork   Size of the workspace array.
!! @param info    On exit, info = 0 for successful exit. If info = -i, the i-th argument had an illegal value.
subroutine oqp_dgetri_i64(n, a, lda, ipiv, work, lwork, info)
    integer :: n, lda, ipiv(:), lwork, info
    real(dp) :: a(lda,*), work(lwork)

    integer(blas_int) :: n_, lda_, lwork_, info_
    logical :: ok

    ! Check arguments
    if (ARG_CHECK) then
        ok = .true.

        ! Check if the size of n, lda, and lwork is within acceptable bounds
        ok = ok .and. (n <= HUGE_BLAS_INT)
        ok = ok .and. (lda <= HUGE_BLAS_INT)
        ok = ok .and. (lwork <= HUGE_BLAS_INT)

        ! If any of the arguments are invalid, display error message and abort
        if (.not. ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    ! Convert integer arguments to blas_int type
    n_ = int(n, blas_int)
    lda_ = int(lda, blas_int)
    lwork_ = int(lwork, blas_int)

    ! Call LAPACK routine dgetri to compute the inverse of the matrix
    call dgetri(n_, a, lda_, ipiv, work, lwork_, info_)

    ! Assign the output info_ to the info argument
    info = info_
end subroutine

!----------------------------------------------------------------------

end module
