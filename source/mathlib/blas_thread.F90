!> @brief Runtime control of BLAS-library internal threading
!> @details Thin Fortran interface to `blas_thread_ctl.c`, which resolves
!>   the thread-count setter/getter of the linked BLAS library (OpenBLAS,
!>   MKL, BLIS) at run time via dlsym().  Use it to switch BLAS to
!>   single-threaded mode around OpenMP-parallel regions that issue many
!>   small BLAS calls, where a BLAS-internal thread pool (e.g. pthread
!>   builds of OpenBLAS) oversubscribes the machine and serializes on its
!>   pool lock.
!>
!>   If the BLAS library is not recognized, `blas_thread_count` returns -1
!>   and `blas_thread_set` is a no-op, so the calls are always safe:
!>
!>     nSaved = blas_thread_count()
!>     call blas_thread_set(1_c_int64_t)
!>     !$omp parallel
!>     ...
!>     !$omp end parallel
!>     call blas_thread_set(nSaved)  ! no-op if nSaved == -1
module blas_thread
  use, intrinsic :: iso_c_binding, only: c_int64_t
  implicit none

  private
  public :: blas_thread_count
  public :: blas_thread_set

  interface
    !> @brief Current BLAS thread count, or -1 if the BLAS library is unknown
    function blas_thread_count() bind(c, name="oqp_blas_thread_count") result(n)
      import :: c_int64_t
      integer(c_int64_t) :: n
    end function

    !> @brief Set BLAS thread count; no-op for n < 1 or unknown BLAS library
    subroutine blas_thread_set(n) bind(c, name="oqp_blas_thread_set")
      import :: c_int64_t
      integer(c_int64_t), value :: n
    end subroutine
  end interface

end module blas_thread
