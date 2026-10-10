module blas_wrap

  use precision, only: dp
  use mathlib_types, only: BLAS_INT, HUGE_BLAS_INT
  use messages, only: show_message, WITH_ABORT

  implicit none

  logical, parameter :: ARG_CHECK = .false.

  character(len=*), parameter :: BITNESS(2) = ["32", "64"]
  character(len=*), parameter :: ERRMSG = &
          "Integer is too big for " // BITNESS(BLAS_INT/4) // "bit BLAS/LAPACK"

  private

  public oqp_caxpy_i64
  public oqp_ccopy_i64
  public oqp_cdotc_i64
  public oqp_cdotu_i64
  public oqp_cgbmv_i64
  public oqp_cgemm_i64
  public oqp_cgemv_i64
  public oqp_cgerc_i64
  public oqp_cgeru_i64
  public oqp_chbmv_i64
  public oqp_chemm_i64
  public oqp_chemv_i64
  public oqp_cher_i64
  public oqp_cher2_i64
  public oqp_cher2k_i64
  public oqp_cherk_i64
  public oqp_chpmv_i64
  public oqp_chpr_i64
  public oqp_chpr2_i64
  public oqp_cscal_i64
  public oqp_csrot_i64
  public oqp_csscal_i64
  public oqp_cswap_i64
  public oqp_csymm_i64
  public oqp_csyr2k_i64
  public oqp_csyrk_i64
  public oqp_ctbmv_i64
  public oqp_ctbsv_i64
  public oqp_ctpmv_i64
  public oqp_ctpsv_i64
  public oqp_ctrmm_i64
  public oqp_ctrmv_i64
  public oqp_ctrsm_i64
  public oqp_ctrsv_i64
  public oqp_dasum_i64
  public oqp_daxpy_i64
  public oqp_dcopy_i64
  public oqp_ddot_i64
  public oqp_dgbmv_i64
  public oqp_dgemm_i64
  public oqp_dgemv_i64
  public oqp_dger_i64
  public oqp_drot_i64
  public oqp_drotm_i64
  public oqp_dsbmv_i64
  public oqp_dscal_i64
  public oqp_dsdot_i64
  public oqp_dspmv_i64
  public oqp_dspr_i64
  public oqp_dspr2_i64
  public oqp_dswap_i64
  public oqp_dsymm_i64
  public oqp_dsymv_i64
  public oqp_dsyr_i64
  public oqp_dsyr2_i64
  public oqp_dsyr2k_i64
  public oqp_dsyrk_i64
  public oqp_dtbmv_i64
  public oqp_dtbsv_i64
  public oqp_dtpmv_i64
  public oqp_dtpsv_i64
  public oqp_dtrmm_i64
  public oqp_dtrmv_i64
  public oqp_dtrsm_i64
  public oqp_dtrsv_i64
  public oqp_dzasum_i64
  public oqp_icamax_i64
  public oqp_idamax_i64
  public oqp_isamax_i64
  public oqp_izamax_i64
  public oqp_sasum_i64
  public oqp_saxpy_i64
  public oqp_scasum_i64
  public oqp_scopy_i64
  public oqp_sdot_i64
  public oqp_sdsdot_i64
  public oqp_sgbmv_i64
  public oqp_sgemm_i64
  public oqp_sgemv_i64
  public oqp_sger_i64
  public oqp_srot_i64
  public oqp_srotm_i64
  public oqp_ssbmv_i64
  public oqp_sscal_i64
  public oqp_sspmv_i64
  public oqp_sspr_i64
  public oqp_sspr2_i64
  public oqp_sswap_i64
  public oqp_ssymm_i64
  public oqp_ssymv_i64
  public oqp_ssyr_i64
  public oqp_ssyr2_i64
  public oqp_ssyr2k_i64
  public oqp_ssyrk_i64
  public oqp_stbmv_i64
  public oqp_stbsv_i64
  public oqp_stpmv_i64
  public oqp_stpsv_i64
  public oqp_strmm_i64
  public oqp_strmv_i64
  public oqp_strsm_i64
  public oqp_strsv_i64
  public oqp_xerbla_i64
!  public oqp_xerbla_array_i64
  public oqp_zaxpy_i64
  public oqp_zcopy_i64
  public oqp_zdotc_i64
  public oqp_zdotu_i64
  public oqp_zdrot_i64
  public oqp_zdscal_i64
  public oqp_zgbmv_i64
  public oqp_zgemm_i64
  public oqp_zgemv_i64
  public oqp_zgerc_i64
  public oqp_zgeru_i64
  public oqp_zhbmv_i64
  public oqp_zhemm_i64
  public oqp_zhemv_i64
  public oqp_zher_i64
  public oqp_zher2_i64
  public oqp_zher2k_i64
  public oqp_zherk_i64
  public oqp_zhpmv_i64
  public oqp_zhpr_i64
  public oqp_zhpr2_i64
  public oqp_zscal_i64
  public oqp_zswap_i64
  public oqp_zsymm_i64
  public oqp_zsyr2k_i64
  public oqp_zsyrk_i64
  public oqp_ztbmv_i64
  public oqp_ztbsv_i64
  public oqp_ztpmv_i64
  public oqp_ztpsv_i64
  public oqp_ztrmm_i64
  public oqp_ztrmv_i64
  public oqp_ztrsm_i64
  public oqp_ztrsv_i64
  public oqp_dnrm2_i64
  public oqp_dznrm2_i64
  public oqp_scnrm2_i64
  public oqp_snrm2_i64

contains

  subroutine oqp_caxpy_i64(n, ca, cx, incx, cy, incy)
    complex :: ca
    integer :: incx
    integer :: incy
    integer :: n
    complex :: cx(*)
    complex :: cy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call caxpy(n_, ca, cx, incx_, cy, incy_)

  end subroutine oqp_caxpy_i64

  subroutine oqp_ccopy_i64(n, cx, incx, cy, incy)
    integer :: incx
    integer :: incy
    integer :: n
    complex :: cx(*)
    complex :: cy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call ccopy(n_, cx, incx_, cy, incy_)

  end subroutine oqp_ccopy_i64

  function oqp_cdotc_i64(n, cx, incx, cy, incy)
    complex, external :: cdotc
    complex :: oqp_cdotc_i64
    integer :: incx
    integer :: incy
    integer :: n
    complex :: cx(*)
    complex :: cy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    oqp_cdotc_i64 = cdotc(n_, cx, incx_, cy, incy_)

  end function oqp_cdotc_i64

  function oqp_cdotu_i64(n, cx, incx, cy, incy)
    complex, external :: cdotu
    complex :: oqp_cdotu_i64
    integer :: incx
    integer :: incy
    integer :: n
    complex :: cx(*)
    complex :: cy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    oqp_cdotu_i64 = cdotu(n_, cx, incx_, cy, incy_)

  end function oqp_cdotu_i64

  subroutine oqp_cgbmv_i64(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    complex :: alpha
    complex :: beta
    integer :: incx
    integer :: incy
    integer :: kl
    integer :: ku
    integer :: lda
    integer :: m
    integer :: n
    character :: trans
    complex :: a(lda,*)
    complex :: x(*)
    complex :: y(*)

    integer(blas_int) :: m_, n_, kl_, ku_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(kl   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ku   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    kl_   = int(kl   , blas_int)
    ku_   = int(ku   , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call cgbmv(trans, m_, n_, kl_, ku_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_cgbmv_i64

  subroutine oqp_cgemm_i64(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    complex :: alpha
    complex :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: transa
    character :: transb
    complex :: a(lda,*)
    complex :: b(ldb,*)
    complex :: c(ldc,*)

    integer(blas_int) :: m_, n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call cgemm(transa, transb, m_, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_cgemm_i64

  subroutine oqp_cgemv_i64(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    complex :: alpha
    complex :: beta
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    character :: trans
    complex :: a(lda,*)
    complex :: x(*)
    complex :: y(*)

    integer(blas_int) :: m_, n_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call cgemv(trans, m_, n_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_cgemv_i64

  subroutine oqp_cgerc_i64(m, n, alpha, x, incx, y, incy, a, lda)
    complex :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    complex :: a(lda,*)
    complex :: x(*)
    complex :: y(*)

    integer(blas_int) :: m_, n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call cgerc(m_, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_cgerc_i64

  subroutine oqp_cgeru_i64(m, n, alpha, x, incx, y, incy, a, lda)
    complex :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    complex :: a(lda,*)
    complex :: x(*)
    complex :: y(*)

    integer(blas_int) :: m_, n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call cgeru(m_, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_cgeru_i64

  subroutine oqp_chbmv_i64(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    complex :: alpha
    complex :: beta
    integer :: incx
    integer :: incy
    integer :: k
    integer :: lda
    integer :: n
    character :: uplo
    complex :: a(lda,*)
    complex :: x(*)
    complex :: y(*)

    integer(blas_int) :: n_, k_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call chbmv(uplo, n_, k_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_chbmv_i64

  subroutine oqp_chemm_i64(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    complex :: alpha
    complex :: beta
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: side
    character :: uplo
    complex :: a(lda,*)
    complex :: b(ldb,*)
    complex :: c(ldc,*)

    integer(blas_int) :: m_, n_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call chemm(side, uplo, m_, n_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_chemm_i64

  subroutine oqp_chemv_i64(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    complex :: alpha
    complex :: beta
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: n
    character :: uplo
    complex :: a(lda,*)
    complex :: x(*)
    complex :: y(*)

    integer(blas_int) :: n_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call chemv(uplo, n_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_chemv_i64

  subroutine oqp_cher_i64(uplo, n, alpha, x, incx, a, lda)
    real :: alpha
    integer :: incx
    integer :: lda
    integer :: n
    character :: uplo
    complex :: a(lda,*)
    complex :: x(*)

    integer(blas_int) :: n_, incx_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    lda_  = int(lda  , blas_int)

    call cher(uplo, n_, alpha, x, incx_, a, lda_)

  end subroutine oqp_cher_i64

  subroutine oqp_cher2_i64(uplo, n, alpha, x, incx, y, incy, a, lda)
    complex :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: n
    character :: uplo
    complex :: a(lda,*)
    complex :: x(*)
    complex :: y(*)

    integer(blas_int) :: n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call cher2(uplo, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_cher2_i64

  subroutine oqp_cher2k_i64(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    complex :: alpha
    real :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    complex :: a(lda,*)
    complex :: b(ldb,*)
    complex :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call cher2k(uplo, trans, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_cher2k_i64

  subroutine oqp_cherk_i64(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    real :: alpha
    real :: beta
    integer :: k
    integer :: lda
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    complex :: a(lda,*)
    complex :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldc_ = int(ldc , blas_int)

    call cherk(uplo, trans, n_, k_, alpha, a, lda_, beta, c, ldc_)

  end subroutine oqp_cherk_i64

  subroutine oqp_chpmv_i64(uplo, n, alpha, ap, x, incx, beta, y, incy)
    complex :: alpha
    complex :: beta
    integer :: incx
    integer :: incy
    integer :: n
    character :: uplo
    complex :: ap(*)
    complex :: x(*)
    complex :: y(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call chpmv(uplo, n_, alpha, ap, x, incx_, beta, y, incy_)

  end subroutine oqp_chpmv_i64

  subroutine oqp_chpr_i64(uplo, n, alpha, x, incx, ap)
    real :: alpha
    integer :: incx
    integer :: n
    character :: uplo
    complex :: ap(*)
    complex :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call chpr(uplo, n_, alpha, x, incx_, ap)

  end subroutine oqp_chpr_i64

  subroutine oqp_chpr2_i64(uplo, n, alpha, x, incx, y, incy, ap)
    complex :: alpha
    integer :: incx
    integer :: incy
    integer :: n
    character :: uplo
    complex :: ap(*)
    complex :: x(*)
    complex :: y(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call chpr2(uplo, n_, alpha, x, incx_, y, incy_, ap)

  end subroutine oqp_chpr2_i64

  subroutine oqp_cscal_i64(n, ca, cx, incx)
    complex :: ca
    integer :: incx
    integer :: n
    complex :: cx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call cscal(n_, ca, cx, incx_)

  end subroutine oqp_cscal_i64

  subroutine oqp_csrot_i64(n, cx, incx, cy, incy, c, s)
    integer :: incx
    integer :: incy
    integer :: n
    real :: c
    real :: s
    complex :: cx(*)
    complex :: cy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call csrot(n_, cx, incx_, cy, incy_, c, s)

  end subroutine oqp_csrot_i64

  subroutine oqp_csscal_i64(n, sa, cx, incx)
    real :: sa
    integer :: incx
    integer :: n
    complex :: cx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call csscal(n_, sa, cx, incx_)

  end subroutine oqp_csscal_i64

  subroutine oqp_cswap_i64(n, cx, incx, cy, incy)
    integer :: incx
    integer :: incy
    integer :: n
    complex :: cx(*)
    complex :: cy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call cswap(n_, cx, incx_, cy, incy_)

  end subroutine oqp_cswap_i64

  subroutine oqp_csymm_i64(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    complex :: alpha
    complex :: beta
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: side
    character :: uplo
    complex :: a(lda,*)
    complex :: b(ldb,*)
    complex :: c(ldc,*)

    integer(blas_int) :: m_, n_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call csymm(side, uplo, m_, n_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_csymm_i64

  subroutine oqp_csyr2k_i64(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    complex :: alpha
    complex :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    complex :: a(lda,*)
    complex :: b(ldb,*)
    complex :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call csyr2k(uplo, trans, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_csyr2k_i64

  subroutine oqp_csyrk_i64(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    complex :: alpha
    complex :: beta
    integer :: k
    integer :: lda
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    complex :: a(lda,*)
    complex :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldc_ = int(ldc , blas_int)

    call csyrk(uplo, trans, n_, k_, alpha, a, lda_, beta, c, ldc_)

  end subroutine oqp_csyrk_i64

  subroutine oqp_ctbmv_i64(uplo, trans, diag, n, k, a, lda, x, incx)
    integer :: incx
    integer :: k
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex :: a(lda,*)
    complex :: x(*)

    integer(blas_int) :: n_, k_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call ctbmv(uplo, trans, diag, n_, k_, a, lda_, x, incx_)

  end subroutine oqp_ctbmv_i64

  subroutine oqp_ctbsv_i64(uplo, trans, diag, n, k, a, lda, x, incx)
    integer :: incx
    integer :: k
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex :: a(lda,*)
    complex :: x(*)

    integer(blas_int) :: n_, k_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call ctbsv(uplo, trans, diag, n_, k_, a, lda_, x, incx_)

  end subroutine oqp_ctbsv_i64

  subroutine oqp_ctpmv_i64(uplo, trans, diag, n, ap, x, incx)
    integer :: incx
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex :: ap(*)
    complex :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call ctpmv(uplo, trans, diag, n_, ap, x, incx_)

  end subroutine oqp_ctpmv_i64

  subroutine oqp_ctpsv_i64(uplo, trans, diag, n, ap, x, incx)
    integer :: incx
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex :: ap(*)
    complex :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call ctpsv(uplo, trans, diag, n_, ap, x, incx_)

  end subroutine oqp_ctpsv_i64

  subroutine oqp_ctrmm_i64(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    complex :: alpha
    integer :: lda
    integer :: ldb
    integer :: m
    integer :: n
    character :: diag
    character :: side
    character :: transa
    character :: uplo
    complex :: a(lda,*)
    complex :: b(ldb,*)

    integer(blas_int) :: m_, n_, lda_, ldb_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)

    call ctrmm(side, uplo, transa, diag, m_, n_, alpha, a, lda_, b, ldb_)

  end subroutine oqp_ctrmm_i64

  subroutine oqp_ctrmv_i64(uplo, trans, diag, n, a, lda, x, incx)
    integer :: incx
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex :: a(lda,*)
    complex :: x(*)

    integer(blas_int) :: n_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call ctrmv(uplo, trans, diag, n_, a, lda_, x, incx_)

  end subroutine oqp_ctrmv_i64

  subroutine oqp_ctrsm_i64(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    complex :: alpha
    integer :: lda
    integer :: ldb
    integer :: m
    integer :: n
    character :: diag
    character :: side
    character :: transa
    character :: uplo
    complex :: a(lda,*)
    complex :: b(ldb,*)

    integer(blas_int) :: m_, n_, lda_, ldb_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)

    call ctrsm(side, uplo, transa, diag, m_, n_, alpha, a, lda_, b, ldb_)

  end subroutine oqp_ctrsm_i64

  subroutine oqp_ctrsv_i64(uplo, trans, diag, n, a, lda, x, incx)
    integer :: incx
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex :: a(lda,*)
    complex :: x(*)

    integer(blas_int) :: n_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call ctrsv(uplo, trans, diag, n_, a, lda_, x, incx_)

  end subroutine oqp_ctrsv_i64

  function oqp_dasum_i64(n, dx, incx)
    double precision, external :: dasum
    double precision :: oqp_dasum_i64
    integer :: incx
    integer :: n
    double precision :: dx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_dasum_i64 = dasum(n_, dx, incx_)

  end function oqp_dasum_i64

  subroutine oqp_daxpy_i64(n, da, dx, incx, dy, incy)
    double precision :: da
    integer :: incx
    integer :: incy
    integer :: n
    double precision :: dx(*)
    double precision :: dy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call daxpy(n_, da, dx, incx_, dy, incy_)

  end subroutine oqp_daxpy_i64

  subroutine oqp_dcopy_i64(n, dx, incx, dy, incy)
    integer :: incx
    integer :: incy
    integer :: n
    double precision :: dx(*)
    double precision :: dy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call dcopy(n_, dx, incx_, dy, incy_)

  end subroutine oqp_dcopy_i64

  function oqp_ddot_i64(n, dx, incx, dy, incy)
    double precision, external :: ddot
    double precision :: oqp_ddot_i64
    integer :: incx
    integer :: incy
    integer :: n
    double precision :: dx(*)
    double precision :: dy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    oqp_ddot_i64 = ddot(n_, dx, incx_, dy, incy_)

  end function oqp_ddot_i64

  subroutine oqp_dgbmv_i64(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    double precision :: alpha
    double precision :: beta
    integer :: incx
    integer :: incy
    integer :: kl
    integer :: ku
    integer :: lda
    integer :: m
    integer :: n
    character :: trans
    double precision :: a(lda,*)
    double precision :: x(*)
    double precision :: y(*)

    integer(blas_int) :: m_, n_, kl_, ku_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(kl   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ku   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    kl_   = int(kl   , blas_int)
    ku_   = int(ku   , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call dgbmv(trans, m_, n_, kl_, ku_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_dgbmv_i64

  subroutine oqp_dgemm_i64(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    double precision :: alpha
    double precision :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: transa
    character :: transb
    double precision :: a(lda,*)
    double precision :: b(ldb,*)
    double precision :: c(ldc,*)

    integer(blas_int) :: m_, n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call dgemm(transa, transb, m_, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_dgemm_i64

  subroutine oqp_dgemv_i64(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    double precision :: alpha
    double precision :: beta
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    character :: trans
    double precision :: a(lda,*)
    double precision :: x(*)
    double precision :: y(*)

    integer(blas_int) :: m_, n_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call dgemv(trans, m_, n_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_dgemv_i64

  subroutine oqp_dger_i64(m, n, alpha, x, incx, y, incy, a, lda)
    double precision :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    double precision :: a(lda,*)
    double precision :: x(*)
    double precision :: y(*)

    integer(blas_int) :: m_, n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call dger(m_, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_dger_i64

  subroutine oqp_drot_i64(n, dx, incx, dy, incy, c, s)
    double precision :: c
    double precision :: s
    integer :: incx
    integer :: incy
    integer :: n
    double precision :: dx(*)
    double precision :: dy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call drot(n_, dx, incx_, dy, incy_, c, s)

  end subroutine oqp_drot_i64

  subroutine oqp_drotm_i64(n, dx, incx, dy, incy, dparam)
    integer :: incx
    integer :: incy
    integer :: n
    double precision :: dparam(5)
    double precision :: dx(*)
    double precision :: dy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call drotm(n_, dx, incx_, dy, incy_, dparam)

  end subroutine oqp_drotm_i64

  subroutine oqp_dsbmv_i64(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    double precision :: alpha
    double precision :: beta
    integer :: incx
    integer :: incy
    integer :: k
    integer :: lda
    integer :: n
    character :: uplo
    double precision :: a(lda,*)
    double precision :: x(*)
    double precision :: y(*)

    integer(blas_int) :: n_, k_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call dsbmv(uplo, n_, k_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_dsbmv_i64

  subroutine oqp_dscal_i64(n, da, dx, incx)
    double precision :: da
    integer :: incx
    integer :: n
    double precision :: dx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call dscal(n_, da, dx, incx_)

  end subroutine oqp_dscal_i64

  function oqp_dsdot_i64(n, sx, incx, sy, incy)
    double precision, external :: dsdot
    double precision :: oqp_dsdot_i64
    integer :: incx
    integer :: incy
    integer :: n
    real :: sx(*)
    real :: sy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    oqp_dsdot_i64 = dsdot(n_, sx, incx_, sy, incy_)

  end function oqp_dsdot_i64

  subroutine oqp_dspmv_i64(uplo, n, alpha, ap, x, incx, beta, y, incy)
    double precision :: alpha
    double precision :: beta
    integer :: incx
    integer :: incy
    integer :: n
    character :: uplo
    double precision :: ap(*)
    double precision :: x(*)
    double precision :: y(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call dspmv(uplo, n_, alpha, ap, x, incx_, beta, y, incy_)

  end subroutine oqp_dspmv_i64

  subroutine oqp_dspr_i64(uplo, n, alpha, x, incx, ap)
    double precision :: alpha
    integer :: incx
    integer :: n
    character :: uplo
    double precision :: ap(*)
    double precision :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call dspr(uplo, n_, alpha, x, incx_, ap)

  end subroutine oqp_dspr_i64

  subroutine oqp_dspr2_i64(uplo, n, alpha, x, incx, y, incy, ap)
    double precision :: alpha
    integer :: incx
    integer :: incy
    integer :: n
    character :: uplo
    double precision :: ap(*)
    double precision :: x(*)
    double precision :: y(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call dspr2(uplo, n_, alpha, x, incx_, y, incy_, ap)

  end subroutine oqp_dspr2_i64

  subroutine oqp_dswap_i64(n, dx, incx, dy, incy)
    integer :: incx
    integer :: incy
    integer :: n
    double precision :: dx(*)
    double precision :: dy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call dswap(n_, dx, incx_, dy, incy_)

  end subroutine oqp_dswap_i64

  subroutine oqp_dsymm_i64(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    double precision :: alpha
    double precision :: beta
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: side
    character :: uplo
    double precision :: a(lda,*)
    double precision :: b(ldb,*)
    double precision :: c(ldc,*)

    integer(blas_int) :: m_, n_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call dsymm(side, uplo, m_, n_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_dsymm_i64

  subroutine oqp_dsymv_i64(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    double precision :: alpha
    double precision :: beta
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: n
    character :: uplo
    double precision :: a(lda,*)
    double precision :: x(*)
    double precision :: y(*)

    integer(blas_int) :: n_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call dsymv(uplo, n_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_dsymv_i64

  subroutine oqp_dsyr_i64(uplo, n, alpha, x, incx, a, lda)
    double precision :: alpha
    integer :: incx
    integer :: lda
    integer :: n
    character :: uplo
    double precision :: a(lda,*)
    double precision :: x(*)

    integer(blas_int) :: n_, incx_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    lda_  = int(lda  , blas_int)

    call dsyr(uplo, n_, alpha, x, incx_, a, lda_)

  end subroutine oqp_dsyr_i64

  subroutine oqp_dsyr2_i64(uplo, n, alpha, x, incx, y, incy, a, lda)
    double precision :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: n
    character :: uplo
    double precision :: a(lda,*)
    double precision :: x(*)
    double precision :: y(*)

    integer(blas_int) :: n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call dsyr2(uplo, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_dsyr2_i64

  subroutine oqp_dsyr2k_i64(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    double precision :: alpha
    double precision :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    double precision :: a(lda,*)
    double precision :: b(ldb,*)
    double precision :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call dsyr2k(uplo, trans, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_dsyr2k_i64

  subroutine oqp_dsyrk_i64(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    double precision :: alpha
    double precision :: beta
    integer :: k
    integer :: lda
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    double precision :: a(lda,*)
    double precision :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldc_ = int(ldc , blas_int)

    call dsyrk(uplo, trans, n_, k_, alpha, a, lda_, beta, c, ldc_)

  end subroutine oqp_dsyrk_i64

  subroutine oqp_dtbmv_i64(uplo, trans, diag, n, k, a, lda, x, incx)
    integer :: incx
    integer :: k
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    double precision :: a(lda,*)
    double precision :: x(*)

    integer(blas_int) :: n_, k_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call dtbmv(uplo, trans, diag, n_, k_, a, lda_, x, incx_)

  end subroutine oqp_dtbmv_i64

  subroutine oqp_dtbsv_i64(uplo, trans, diag, n, k, a, lda, x, incx)
    integer :: incx
    integer :: k
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    double precision :: a(lda,*)
    double precision :: x(*)

    integer(blas_int) :: n_, k_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call dtbsv(uplo, trans, diag, n_, k_, a, lda_, x, incx_)

  end subroutine oqp_dtbsv_i64

  subroutine oqp_dtpmv_i64(uplo, trans, diag, n, ap, x, incx)
    integer :: incx
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    double precision :: ap(*)
    double precision :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call dtpmv(uplo, trans, diag, n_, ap, x, incx_)

  end subroutine oqp_dtpmv_i64

  subroutine oqp_dtpsv_i64(uplo, trans, diag, n, ap, x, incx)
    integer :: incx
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    double precision :: ap(*)
    double precision :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call dtpsv(uplo, trans, diag, n_, ap, x, incx_)

  end subroutine oqp_dtpsv_i64

  subroutine oqp_dtrmm_i64(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    double precision :: alpha
    integer :: lda
    integer :: ldb
    integer :: m
    integer :: n
    character :: diag
    character :: side
    character :: transa
    character :: uplo
    double precision :: a(lda,*)
    double precision :: b(ldb,*)

    integer(blas_int) :: m_, n_, lda_, ldb_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)

    call dtrmm(side, uplo, transa, diag, m_, n_, alpha, a, lda_, b, ldb_)

  end subroutine oqp_dtrmm_i64

  subroutine oqp_dtrmv_i64(uplo, trans, diag, n, a, lda, x, incx)
    integer :: incx
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    double precision :: a(lda,*)
    double precision :: x(*)

    integer(blas_int) :: n_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call dtrmv(uplo, trans, diag, n_, a, lda_, x, incx_)

  end subroutine oqp_dtrmv_i64

  subroutine oqp_dtrsm_i64(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    double precision :: alpha
    integer :: lda
    integer :: ldb
    integer :: m
    integer :: n
    character :: diag
    character :: side
    character :: transa
    character :: uplo
    double precision :: a(lda,*)
    double precision :: b(ldb,*)

    integer(blas_int) :: m_, n_, lda_, ldb_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)

    call dtrsm(side, uplo, transa, diag, m_, n_, alpha, a, lda_, b, ldb_)

  end subroutine oqp_dtrsm_i64

  subroutine oqp_dtrsv_i64(uplo, trans, diag, n, a, lda, x, incx)
    integer :: incx
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    double precision :: a(lda,*)
    double precision :: x(*)

    integer(blas_int) :: n_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call dtrsv(uplo, trans, diag, n_, a, lda_, x, incx_)

  end subroutine oqp_dtrsv_i64

  function oqp_dzasum_i64(n, zx, incx)
    double precision, external :: dzasum
    double precision :: oqp_dzasum_i64
    integer :: incx
    integer :: n
    complex(kind=8) :: zx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_dzasum_i64 = dzasum(n_, zx, incx_)

  end function oqp_dzasum_i64

  function oqp_icamax_i64(n, cx, incx)
    integer(blas_int), external :: icamax
    integer :: oqp_icamax_i64
    integer :: incx
    integer :: n
    complex :: cx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_icamax_i64 = icamax(n_, cx, incx_)

  end function oqp_icamax_i64

  function oqp_idamax_i64(n, dx, incx)
    integer(blas_int), external :: idamax
    integer :: oqp_idamax_i64
    integer :: incx
    integer :: n
    double precision :: dx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_idamax_i64 = idamax(n_, dx, incx_)

  end function oqp_idamax_i64

  function oqp_isamax_i64(n, sx, incx)
    integer(blas_int), external :: isamax
    integer :: oqp_isamax_i64
    integer :: incx
    integer :: n
    real :: sx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_isamax_i64 = isamax(n_, sx, incx_)

  end function oqp_isamax_i64

  function oqp_izamax_i64(n, zx, incx)
    integer(blas_int), external :: izamax
    integer :: oqp_izamax_i64
    integer :: incx
    integer :: n
    complex(kind=8) :: zx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_izamax_i64 = izamax(n_, zx, incx_)

  end function oqp_izamax_i64

  function oqp_sasum_i64(n, sx, incx)
    real, external :: sasum
    real :: oqp_sasum_i64
    integer :: incx
    integer :: n
    real :: sx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_sasum_i64 = sasum(n_, sx, incx_)

  end function oqp_sasum_i64

  subroutine oqp_saxpy_i64(n, sa, sx, incx, sy, incy)
    real :: sa
    integer :: incx
    integer :: incy
    integer :: n
    real :: sx(*)
    real :: sy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call saxpy(n_, sa, sx, incx_, sy, incy_)

  end subroutine oqp_saxpy_i64

  function oqp_scasum_i64(n, cx, incx)
    real, external :: scasum
    real :: oqp_scasum_i64
    integer :: incx
    integer :: n
    complex :: cx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_scasum_i64 = scasum(n_, cx, incx_)

  end function oqp_scasum_i64

  subroutine oqp_scopy_i64(n, sx, incx, sy, incy)
    integer :: incx
    integer :: incy
    integer :: n
    real :: sx(*)
    real :: sy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call scopy(n_, sx, incx_, sy, incy_)

  end subroutine oqp_scopy_i64

  function oqp_sdot_i64(n, sx, incx, sy, incy)
    real, external :: sdot
    real :: oqp_sdot_i64
    integer :: incx
    integer :: incy
    integer :: n
    real :: sx(*)
    real :: sy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    oqp_sdot_i64 = sdot(n_, sx, incx_, sy, incy_)

  end function oqp_sdot_i64

  function oqp_sdsdot_i64(n, sb, sx, incx, sy, incy)
    real, external :: sdsdot
    real :: oqp_sdsdot_i64
    real :: sb
    integer :: incx
    integer :: incy
    integer :: n
    real :: sx(*)
    real :: sy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    oqp_sdsdot_i64 = sdsdot(n_, sb, sx, incx_, sy, incy_)

  end function oqp_sdsdot_i64

  subroutine oqp_sgbmv_i64(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    real :: alpha
    real :: beta
    integer :: incx
    integer :: incy
    integer :: kl
    integer :: ku
    integer :: lda
    integer :: m
    integer :: n
    character :: trans
    real :: a(lda,*)
    real :: x(*)
    real :: y(*)

    integer(blas_int) :: m_, n_, kl_, ku_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(kl   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ku   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    kl_   = int(kl   , blas_int)
    ku_   = int(ku   , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call sgbmv(trans, m_, n_, kl_, ku_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_sgbmv_i64

  subroutine oqp_sgemm_i64(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    real :: alpha
    real :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: transa
    character :: transb
    real :: a(lda,*)
    real :: b(ldb,*)
    real :: c(ldc,*)

    integer(blas_int) :: m_, n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call sgemm(transa, transb, m_, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_sgemm_i64

  subroutine oqp_sgemv_i64(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    real :: alpha
    real :: beta
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    character :: trans
    real :: a(lda,*)
    real :: x(*)
    real :: y(*)

    integer(blas_int) :: m_, n_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call sgemv(trans, m_, n_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_sgemv_i64

  subroutine oqp_sger_i64(m, n, alpha, x, incx, y, incy, a, lda)
    real :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    real :: a(lda,*)
    real :: x(*)
    real :: y(*)

    integer(blas_int) :: m_, n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call sger(m_, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_sger_i64

  subroutine oqp_srot_i64(n, sx, incx, sy, incy, c, s)
    real :: c
    real :: s
    integer :: incx
    integer :: incy
    integer :: n
    real :: sx(*)
    real :: sy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call srot(n_, sx, incx_, sy, incy_, c, s)

  end subroutine oqp_srot_i64

  subroutine oqp_srotm_i64(n, sx, incx, sy, incy, sparam)
    integer :: incx
    integer :: incy
    integer :: n
    real :: sparam(5)
    real :: sx(*)
    real :: sy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call srotm(n_, sx, incx_, sy, incy_, sparam)

  end subroutine oqp_srotm_i64

  subroutine oqp_ssbmv_i64(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    real :: alpha
    real :: beta
    integer :: incx
    integer :: incy
    integer :: k
    integer :: lda
    integer :: n
    character :: uplo
    real :: a(lda,*)
    real :: x(*)
    real :: y(*)

    integer(blas_int) :: n_, k_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call ssbmv(uplo, n_, k_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_ssbmv_i64

  subroutine oqp_sscal_i64(n, sa, sx, incx)
    real :: sa
    integer :: incx
    integer :: n
    real :: sx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call sscal(n_, sa, sx, incx_)

  end subroutine oqp_sscal_i64

  subroutine oqp_sspmv_i64(uplo, n, alpha, ap, x, incx, beta, y, incy)
    real :: alpha
    real :: beta
    integer :: incx
    integer :: incy
    integer :: n
    character :: uplo
    real :: ap(*)
    real :: x(*)
    real :: y(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call sspmv(uplo, n_, alpha, ap, x, incx_, beta, y, incy_)

  end subroutine oqp_sspmv_i64

  subroutine oqp_sspr_i64(uplo, n, alpha, x, incx, ap)
    real :: alpha
    integer :: incx
    integer :: n
    character :: uplo
    real :: ap(*)
    real :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call sspr(uplo, n_, alpha, x, incx_, ap)

  end subroutine oqp_sspr_i64

  subroutine oqp_sspr2_i64(uplo, n, alpha, x, incx, y, incy, ap)
    real :: alpha
    integer :: incx
    integer :: incy
    integer :: n
    character :: uplo
    real :: ap(*)
    real :: x(*)
    real :: y(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call sspr2(uplo, n_, alpha, x, incx_, y, incy_, ap)

  end subroutine oqp_sspr2_i64

  subroutine oqp_sswap_i64(n, sx, incx, sy, incy)
    integer :: incx
    integer :: incy
    integer :: n
    real :: sx(*)
    real :: sy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call sswap(n_, sx, incx_, sy, incy_)

  end subroutine oqp_sswap_i64

  subroutine oqp_ssymm_i64(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    real :: alpha
    real :: beta
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: side
    character :: uplo
    real :: a(lda,*)
    real :: b(ldb,*)
    real :: c(ldc,*)

    integer(blas_int) :: m_, n_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call ssymm(side, uplo, m_, n_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_ssymm_i64

  subroutine oqp_ssymv_i64(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    real :: alpha
    real :: beta
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: n
    character :: uplo
    real :: a(lda,*)
    real :: x(*)
    real :: y(*)

    integer(blas_int) :: n_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call ssymv(uplo, n_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_ssymv_i64

  subroutine oqp_ssyr_i64(uplo, n, alpha, x, incx, a, lda)
    real :: alpha
    integer :: incx
    integer :: lda
    integer :: n
    character :: uplo
    real :: a(lda,*)
    real :: x(*)

    integer(blas_int) :: n_, incx_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    lda_  = int(lda  , blas_int)

    call ssyr(uplo, n_, alpha, x, incx_, a, lda_)

  end subroutine oqp_ssyr_i64

  subroutine oqp_ssyr2_i64(uplo, n, alpha, x, incx, y, incy, a, lda)
    real :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: n
    character :: uplo
    real :: a(lda,*)
    real :: x(*)
    real :: y(*)

    integer(blas_int) :: n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call ssyr2(uplo, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_ssyr2_i64

  subroutine oqp_ssyr2k_i64(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    real :: alpha
    real :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    real :: a(lda,*)
    real :: b(ldb,*)
    real :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call ssyr2k(uplo, trans, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_ssyr2k_i64

  subroutine oqp_ssyrk_i64(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    real :: alpha
    real :: beta
    integer :: k
    integer :: lda
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    real :: a(lda,*)
    real :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldc_ = int(ldc , blas_int)

    call ssyrk(uplo, trans, n_, k_, alpha, a, lda_, beta, c, ldc_)

  end subroutine oqp_ssyrk_i64

  subroutine oqp_stbmv_i64(uplo, trans, diag, n, k, a, lda, x, incx)
    integer :: incx
    integer :: k
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    real :: a(lda,*)
    real :: x(*)

    integer(blas_int) :: n_, k_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call stbmv(uplo, trans, diag, n_, k_, a, lda_, x, incx_)

  end subroutine oqp_stbmv_i64

  subroutine oqp_stbsv_i64(uplo, trans, diag, n, k, a, lda, x, incx)
    integer :: incx
    integer :: k
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    real :: a(lda,*)
    real :: x(*)

    integer(blas_int) :: n_, k_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call stbsv(uplo, trans, diag, n_, k_, a, lda_, x, incx_)

  end subroutine oqp_stbsv_i64

  subroutine oqp_stpmv_i64(uplo, trans, diag, n, ap, x, incx)
    integer :: incx
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    real :: ap(*)
    real :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call stpmv(uplo, trans, diag, n_, ap, x, incx_)

  end subroutine oqp_stpmv_i64

  subroutine oqp_stpsv_i64(uplo, trans, diag, n, ap, x, incx)
    integer :: incx
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    real :: ap(*)
    real :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call stpsv(uplo, trans, diag, n_, ap, x, incx_)

  end subroutine oqp_stpsv_i64

  subroutine oqp_strmm_i64(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    real :: alpha
    integer :: lda
    integer :: ldb
    integer :: m
    integer :: n
    character :: diag
    character :: side
    character :: transa
    character :: uplo
    real :: a(lda,*)
    real :: b(ldb,*)

    integer(blas_int) :: m_, n_, lda_, ldb_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)

    call strmm(side, uplo, transa, diag, m_, n_, alpha, a, lda_, b, ldb_)

  end subroutine oqp_strmm_i64

  subroutine oqp_strmv_i64(uplo, trans, diag, n, a, lda, x, incx)
    integer :: incx
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    real :: a(lda,*)
    real :: x(*)

    integer(blas_int) :: n_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call strmv(uplo, trans, diag, n_, a, lda_, x, incx_)

  end subroutine oqp_strmv_i64

  subroutine oqp_strsm_i64(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    real :: alpha
    integer :: lda
    integer :: ldb
    integer :: m
    integer :: n
    character :: diag
    character :: side
    character :: transa
    character :: uplo
    real :: a(lda,*)
    real :: b(ldb,*)

    integer(blas_int) :: m_, n_, lda_, ldb_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)

    call strsm(side, uplo, transa, diag, m_, n_, alpha, a, lda_, b, ldb_)

  end subroutine oqp_strsm_i64

  subroutine oqp_strsv_i64(uplo, trans, diag, n, a, lda, x, incx)
    integer :: incx
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    real :: a(lda,*)
    real :: x(*)

    integer(blas_int) :: n_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call strsv(uplo, trans, diag, n_, a, lda_, x, incx_)

  end subroutine oqp_strsv_i64

  subroutine oqp_xerbla_i64(srname, info)
    character :: srname
    integer :: info

    integer(blas_int) :: info_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(info ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    info_ = int(info , blas_int)

    call xerbla(srname, info_)

  end subroutine oqp_xerbla_i64

!  subroutine oqp_xerbla_array_i64(srname_array, srname_len, info)
!    integer :: srname_len
!    integer :: info
!    character :: srname_array(srname_len)
!
!    integer(blas_int) :: srname_len_, info_
!    logical :: ok
!
!    if (ARG_CHECK) then
!      ok = .true.
!      ok = ok .and. abs(srname_len ) <= HUGE_BLAS_INT-1
!      ok = ok .and. abs(info       ) <= HUGE_BLAS_INT-1
!      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
!    end if
!
!    srname_len_ = int(srname_len , blas_int)
!    info_       = int(info       , blas_int)
!
!    call xerbla_array(srname_array, srname_len_, info_)
!
!  end subroutine oqp_xerbla_array_i64

  subroutine oqp_zaxpy_i64(n, za, zx, incx, zy, incy)
    complex(kind=8) :: za
    integer :: incx
    integer :: incy
    integer :: n
    complex(kind=8) :: zx(*)
    complex(kind=8) :: zy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zaxpy(n_, za, zx, incx_, zy, incy_)

  end subroutine oqp_zaxpy_i64

  subroutine oqp_zcopy_i64(n, zx, incx, zy, incy)
    integer :: incx
    integer :: incy
    integer :: n
    complex(kind=8) :: zx(*)
    complex(kind=8) :: zy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zcopy(n_, zx, incx_, zy, incy_)

  end subroutine oqp_zcopy_i64

  function oqp_zdotc_i64(n, zx, incx, zy, incy)
    complex(kind=8), external :: zdotc
    complex(kind=8) :: oqp_zdotc_i64
    integer :: incx
    integer :: incy
    integer :: n
    complex(kind=8) :: zx(*)
    complex(kind=8) :: zy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    oqp_zdotc_i64 = zdotc(n_, zx, incx_, zy, incy_)

  end function oqp_zdotc_i64

  function oqp_zdotu_i64(n, zx, incx, zy, incy)
    complex(kind=8), external :: zdotu
    complex(kind=8) :: oqp_zdotu_i64
    integer :: incx
    integer :: incy
    integer :: n
    complex(kind=8) :: zx(*)
    complex(kind=8) :: zy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    oqp_zdotu_i64 = zdotu(n_, zx, incx_, zy, incy_)

  end function oqp_zdotu_i64

  subroutine oqp_zdrot_i64(n, zx, incx, zy, incy, c, s)
    integer :: incx
    integer :: incy
    integer :: n
    double precision :: c
    double precision :: s
    complex(kind=8) :: zx(*)
    complex(kind=8) :: zy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zdrot(n_, zx, incx_, zy, incy_, c, s)

  end subroutine oqp_zdrot_i64

  subroutine oqp_zdscal_i64(n, da, zx, incx)
    double precision :: da
    integer :: incx
    integer :: n
    complex(kind=8) :: zx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call zdscal(n_, da, zx, incx_)

  end subroutine oqp_zdscal_i64

  subroutine oqp_zgbmv_i64(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: incx
    integer :: incy
    integer :: kl
    integer :: ku
    integer :: lda
    integer :: m
    integer :: n
    character :: trans
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)
    complex(kind=8) :: y(*)

    integer(blas_int) :: m_, n_, kl_, ku_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(kl   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ku   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    kl_   = int(kl   , blas_int)
    ku_   = int(ku   , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zgbmv(trans, m_, n_, kl_, ku_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_zgbmv_i64

  subroutine oqp_zgemm_i64(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: transa
    character :: transb
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: b(ldb,*)
    complex(kind=8) :: c(ldc,*)

    integer(blas_int) :: m_, n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call zgemm(transa, transb, m_, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_zgemm_i64

  subroutine oqp_zgemv_i64(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    character :: trans
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)
    complex(kind=8) :: y(*)

    integer(blas_int) :: m_, n_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zgemv(trans, m_, n_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_zgemv_i64

  subroutine oqp_zgerc_i64(m, n, alpha, x, incx, y, incy, a, lda)
    complex(kind=8) :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)
    complex(kind=8) :: y(*)

    integer(blas_int) :: m_, n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call zgerc(m_, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_zgerc_i64

  subroutine oqp_zgeru_i64(m, n, alpha, x, incx, y, incy, a, lda)
    complex(kind=8) :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: m
    integer :: n
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)
    complex(kind=8) :: y(*)

    integer(blas_int) :: m_, n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_    = int(m    , blas_int)
    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call zgeru(m_, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_zgeru_i64

  subroutine oqp_zhbmv_i64(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: incx
    integer :: incy
    integer :: k
    integer :: lda
    integer :: n
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)
    complex(kind=8) :: y(*)

    integer(blas_int) :: n_, k_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zhbmv(uplo, n_, k_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_zhbmv_i64

  subroutine oqp_zhemm_i64(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: side
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: b(ldb,*)
    complex(kind=8) :: c(ldc,*)

    integer(blas_int) :: m_, n_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call zhemm(side, uplo, m_, n_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_zhemm_i64

  subroutine oqp_zhemv_i64(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: n
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)
    complex(kind=8) :: y(*)

    integer(blas_int) :: n_, lda_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zhemv(uplo, n_, alpha, a, lda_, x, incx_, beta, y, incy_)

  end subroutine oqp_zhemv_i64

  subroutine oqp_zher_i64(uplo, n, alpha, x, incx, a, lda)
    double precision :: alpha
    integer :: incx
    integer :: lda
    integer :: n
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)

    integer(blas_int) :: n_, incx_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    lda_  = int(lda  , blas_int)

    call zher(uplo, n_, alpha, x, incx_, a, lda_)

  end subroutine oqp_zher_i64

  subroutine oqp_zher2_i64(uplo, n, alpha, x, incx, y, incy, a, lda)
    complex(kind=8) :: alpha
    integer :: incx
    integer :: incy
    integer :: lda
    integer :: n
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)
    complex(kind=8) :: y(*)

    integer(blas_int) :: n_, incx_, incy_, lda_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)
    lda_  = int(lda  , blas_int)

    call zher2(uplo, n_, alpha, x, incx_, y, incy_, a, lda_)

  end subroutine oqp_zher2_i64

  subroutine oqp_zher2k_i64(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    complex(kind=8) :: alpha
    double precision :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: b(ldb,*)
    complex(kind=8) :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call zher2k(uplo, trans, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_zher2k_i64

  subroutine oqp_zherk_i64(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    double precision :: alpha
    double precision :: beta
    integer :: k
    integer :: lda
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldc_ = int(ldc , blas_int)

    call zherk(uplo, trans, n_, k_, alpha, a, lda_, beta, c, ldc_)

  end subroutine oqp_zherk_i64

  subroutine oqp_zhpmv_i64(uplo, n, alpha, ap, x, incx, beta, y, incy)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: incx
    integer :: incy
    integer :: n
    character :: uplo
    complex(kind=8) :: ap(*)
    complex(kind=8) :: x(*)
    complex(kind=8) :: y(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zhpmv(uplo, n_, alpha, ap, x, incx_, beta, y, incy_)

  end subroutine oqp_zhpmv_i64

  subroutine oqp_zhpr_i64(uplo, n, alpha, x, incx, ap)
    double precision :: alpha
    integer :: incx
    integer :: n
    character :: uplo
    complex(kind=8) :: ap(*)
    complex(kind=8) :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call zhpr(uplo, n_, alpha, x, incx_, ap)

  end subroutine oqp_zhpr_i64

  subroutine oqp_zhpr2_i64(uplo, n, alpha, x, incx, y, incy, ap)
    complex(kind=8) :: alpha
    integer :: incx
    integer :: incy
    integer :: n
    character :: uplo
    complex(kind=8) :: ap(*)
    complex(kind=8) :: x(*)
    complex(kind=8) :: y(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zhpr2(uplo, n_, alpha, x, incx_, y, incy_, ap)

  end subroutine oqp_zhpr2_i64

  subroutine oqp_zscal_i64(n, za, zx, incx)
    complex(kind=8) :: za
    integer :: incx
    integer :: n
    complex(kind=8) :: zx(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call zscal(n_, za, zx, incx_)

  end subroutine oqp_zscal_i64

  subroutine oqp_zswap_i64(n, zx, incx, zy, incy)
    integer :: incx
    integer :: incy
    integer :: n
    complex(kind=8) :: zx(*)
    complex(kind=8) :: zy(*)

    integer(blas_int) :: n_, incx_, incy_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incy ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)
    incy_ = int(incy , blas_int)

    call zswap(n_, zx, incx_, zy, incy_)

  end subroutine oqp_zswap_i64

  subroutine oqp_zsymm_i64(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: m
    integer :: n
    character :: side
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: b(ldb,*)
    complex(kind=8) :: c(ldc,*)

    integer(blas_int) :: m_, n_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call zsymm(side, uplo, m_, n_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_zsymm_i64

  subroutine oqp_zsyr2k_i64(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: b(ldb,*)
    complex(kind=8) :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldb_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)
    ldc_ = int(ldc , blas_int)

    call zsyr2k(uplo, trans, n_, k_, alpha, a, lda_, b, ldb_, beta, c, ldc_)

  end subroutine oqp_zsyr2k_i64

  subroutine oqp_zsyrk_i64(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    complex(kind=8) :: alpha
    complex(kind=8) :: beta
    integer :: k
    integer :: lda
    integer :: ldc
    integer :: n
    character :: trans
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: c(ldc,*)

    integer(blas_int) :: n_, k_, lda_, ldc_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldc ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_   = int(n   , blas_int)
    k_   = int(k   , blas_int)
    lda_ = int(lda , blas_int)
    ldc_ = int(ldc , blas_int)

    call zsyrk(uplo, trans, n_, k_, alpha, a, lda_, beta, c, ldc_)

  end subroutine oqp_zsyrk_i64

  subroutine oqp_ztbmv_i64(uplo, trans, diag, n, k, a, lda, x, incx)
    integer :: incx
    integer :: k
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)

    integer(blas_int) :: n_, k_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call ztbmv(uplo, trans, diag, n_, k_, a, lda_, x, incx_)

  end subroutine oqp_ztbmv_i64

  subroutine oqp_ztbsv_i64(uplo, trans, diag, n, k, a, lda, x, incx)
    integer :: incx
    integer :: k
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)

    integer(blas_int) :: n_, k_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(k    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    k_    = int(k    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call ztbsv(uplo, trans, diag, n_, k_, a, lda_, x, incx_)

  end subroutine oqp_ztbsv_i64

  subroutine oqp_ztpmv_i64(uplo, trans, diag, n, ap, x, incx)
    integer :: incx
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex(kind=8) :: ap(*)
    complex(kind=8) :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call ztpmv(uplo, trans, diag, n_, ap, x, incx_)

  end subroutine oqp_ztpmv_i64

  subroutine oqp_ztpsv_i64(uplo, trans, diag, n, ap, x, incx)
    integer :: incx
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex(kind=8) :: ap(*)
    complex(kind=8) :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    call ztpsv(uplo, trans, diag, n_, ap, x, incx_)

  end subroutine oqp_ztpsv_i64

  subroutine oqp_ztrmm_i64(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    complex(kind=8) :: alpha
    integer :: lda
    integer :: ldb
    integer :: m
    integer :: n
    character :: diag
    character :: side
    character :: transa
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: b(ldb,*)

    integer(blas_int) :: m_, n_, lda_, ldb_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)

    call ztrmm(side, uplo, transa, diag, m_, n_, alpha, a, lda_, b, ldb_)

  end subroutine oqp_ztrmm_i64

  subroutine oqp_ztrmv_i64(uplo, trans, diag, n, a, lda, x, incx)
    integer :: incx
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)

    integer(blas_int) :: n_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call ztrmv(uplo, trans, diag, n_, a, lda_, x, incx_)

  end subroutine oqp_ztrmv_i64

  subroutine oqp_ztrsm_i64(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    complex(kind=8) :: alpha
    integer :: lda
    integer :: ldb
    integer :: m
    integer :: n
    character :: diag
    character :: side
    character :: transa
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: b(ldb,*)

    integer(blas_int) :: m_, n_, lda_, ldb_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(m   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(n   ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(ldb ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    m_   = int(m   , blas_int)
    n_   = int(n   , blas_int)
    lda_ = int(lda , blas_int)
    ldb_ = int(ldb , blas_int)

    call ztrsm(side, uplo, transa, diag, m_, n_, alpha, a, lda_, b, ldb_)

  end subroutine oqp_ztrsm_i64

  subroutine oqp_ztrsv_i64(uplo, trans, diag, n, a, lda, x, incx)
    integer :: incx
    integer :: lda
    integer :: n
    character :: diag
    character :: trans
    character :: uplo
    complex(kind=8) :: a(lda,*)
    complex(kind=8) :: x(*)

    integer(blas_int) :: n_, lda_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(lda  ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    lda_  = int(lda  , blas_int)
    incx_ = int(incx , blas_int)

    call ztrsv(uplo, trans, diag, n_, a, lda_, x, incx_)

  end subroutine oqp_ztrsv_i64

  function oqp_dnrm2_i64(n, x, incx)
    real(kind=8), external :: dnrm2
    real(kind=8) :: oqp_dnrm2_i64
    integer :: incx
    integer :: n
    real(kind=8) :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_dnrm2_i64 = dnrm2(n_, x, incx_)

  end function oqp_dnrm2_i64

  function oqp_dznrm2_i64(n, x, incx)
    real(kind=8), external :: dznrm2
    real(kind=8) :: oqp_dznrm2_i64
    integer :: incx
    integer :: n
    complex(kind=8) :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_dznrm2_i64 = dznrm2(n_, x, incx_)

  end function oqp_dznrm2_i64

  function oqp_scnrm2_i64(n, x, incx)
    real(kind=4), external :: scnrm2
    real(kind=4) :: oqp_scnrm2_i64
    integer :: incx
    integer :: n
    complex(kind=4) :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_scnrm2_i64 = scnrm2(n_, x, incx_)

  end function oqp_scnrm2_i64

  function oqp_snrm2_i64(n, x, incx)
    real(kind=4), external :: snrm2
    real(kind=4) :: oqp_snrm2_i64
    integer :: incx
    integer :: n
    real(kind=4) :: x(*)

    integer(blas_int) :: n_, incx_
    logical :: ok

    if (ARG_CHECK) then
      ok = .true.
      ok = ok .and. abs(n    ) <= HUGE_BLAS_INT-1
      ok = ok .and. abs(incx ) <= HUGE_BLAS_INT-1
      if (.not.ok) call show_message(ERRMSG, WITH_ABORT)
    end if

    n_    = int(n    , blas_int)
    incx_ = int(incx , blas_int)

    oqp_snrm2_i64 = snrm2(n_, x, incx_)

  end function oqp_snrm2_i64
end module
