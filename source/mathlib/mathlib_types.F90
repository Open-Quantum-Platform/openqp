#ifndef OQP_BLAS_INT
#define OQP_BLAS_INT 4
#endif
module mathlib_types
  implicit none
  integer, parameter :: BLAS_INT = OQP_BLAS_INT
  integer, parameter :: HUGE_BLAS_INT = HUGE(1_BLAS_INT)
  logical, parameter :: INT_REQ_CONV = storage_size(1_blas_int) == storage_size(1)
end module mathlib_types
