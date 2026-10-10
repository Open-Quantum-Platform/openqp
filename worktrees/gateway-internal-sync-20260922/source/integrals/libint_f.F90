MODULE libint_f
   USE ISO_C_BINDING, ONLY: C_DOUBLE, C_PTR, C_NULL_PTR, C_INT, C_FUNPTR, C_F_POINTER, C_F_PROCPOINTER, C_SIZE_T

#ifdef OQP_LIBINT_ENABLE
#include "libint2/config.h"
#include "libint2/util/generated/libint2_params.h"
#include "fortran_incldefs.h"
#endif

   IMPLICIT NONE
   private
   public :: libint_t
   public :: libint2_static_init
   public :: libint2_static_cleanup
   public :: libint2_init_eri
   public :: libint2_cleanup_eri
   public :: libint2_build
#ifdef INCLUDE_ERI
   public :: libint2_build_eri
#if INCLUDE_ERI >= 1
   public :: libint2_build_eri1
#endif
#if INCLUDE_ERI >= 2
   public :: libint2_build_eri2
#endif
#endif
   public :: libint2_active

   logical, parameter :: libint2_active = &
#ifdef OQP_LIBINT_ENABLE
    .true.
#else
    .false.
#endif

#ifndef OQP_LIBINT_ENABLE
  type, bind(C) :: libint_t
    type(c_ptr) :: targets(1)
  end type
#endif

#ifdef OQP_LIBINT_ENABLE
#ifdef LIBINT2_MAX_AM
   INTEGER, PARAMETER :: libint2_max_am = LIBINT2_MAX_AM
#endif
#ifdef LIBINT2_MAX_AM_default
   INTEGER, PARAMETER :: libint2_max_am_default = LIBINT2_MAX_AM_default
#else
#  error "LIBINT2_MAX_AM_default is expected to be defined, libint2_params.h is misgenerated"
#endif
#ifdef LIBINT2_MAX_AM_default1
   INTEGER, PARAMETER :: libint2_max_am_default1 = LIBINT2_MAX_AM_default1
#else
   INTEGER, PARAMETER :: libint2_max_am_default1 = LIBINT2_MAX_AM_default
#endif
#ifdef LIBINT2_MAX_AM_default2
   INTEGER, PARAMETER :: libint2_max_am_default2 = LIBINT2_MAX_AM_default2
#else
   INTEGER, PARAMETER :: libint2_max_am_default2 = LIBINT2_MAX_AM_default
#endif
#ifdef LIBINT2_MAX_AM_eri
   INTEGER, PARAMETER :: libint2_max_am_eri = LIBINT2_MAX_AM_eri
#endif
#ifdef LIBINT2_MAX_AM_eri1
   INTEGER, PARAMETER :: libint2_max_am_eri1 = LIBINT2_MAX_AM_eri1
#endif
#ifdef LIBINT2_MAX_AM_eri2
   INTEGER, PARAMETER :: libint2_max_am_eri2 = LIBINT2_MAX_AM_eri2
#endif
#ifdef LIBINT2_MAX_AM_3eri
   INTEGER, PARAMETER :: libint2_max_am_3eri = LIBINT2_MAX_AM_3eri
#endif
#ifdef LIBINT2_MAX_AM_3eri1
   INTEGER, PARAMETER :: libint2_max_am_3eri1 = LIBINT2_MAX_AM_3eri1
#endif
#ifdef LIBINT2_MAX_AM_3eri2
   INTEGER, PARAMETER :: libint2_max_am_3eri2 = LIBINT2_MAX_AM_3eri2
#endif
#ifdef LIBINT2_MAX_AM_2eri
   INTEGER, PARAMETER :: libint2_max_am_2eri = LIBINT2_MAX_AM_2eri
#endif
#ifdef LIBINT2_MAX_AM_2eri1
   INTEGER, PARAMETER :: libint2_max_am_2eri1 = LIBINT2_MAX_AM_2eri1
#endif
#ifdef LIBINT2_MAX_AM_2eri2
   INTEGER, PARAMETER :: libint2_max_am_2eri2 = LIBINT2_MAX_AM_2eri2
#endif

   INTEGER, PARAMETER :: libint2_max_veclen = LIBINT2_MAX_VECLEN

#include "libint2_types_f.h"

#ifdef INCLUDE_ERI
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri, 0:libint2_max_am_eri, 0:libint2_max_am_eri, 0:libint2_max_am_eri), &
      BIND(C) :: libint2_build_eri
#if INCLUDE_ERI >= 1
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri1, 0:libint2_max_am_eri1, 0:libint2_max_am_eri1, 0:libint2_max_am_eri1), &
      BIND(C) :: libint2_build_eri1
#endif
#if INCLUDE_ERI >= 2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_eri2, 0:libint2_max_am_eri2, 0:libint2_max_am_eri2, 0:libint2_max_am_eri2), &
      BIND(C) :: libint2_build_eri2
#endif
#endif

#ifdef INCLUDE_ERI2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_2eri, 0:libint2_max_am_2eri), &
      BIND(C) :: libint2_build_2eri
#if INCLUDE_ERI2 >= 1
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_2eri1, 0:libint2_max_am_2eri1), &
      BIND(C) :: libint2_build_2eri1
#endif
#if INCLUDE_ERI2 >= 2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_2eri2, 0:libint2_max_am_2eri2), &
      BIND(C) :: libint2_build_2eri2
#endif
#endif

#ifdef INCLUDE_ERI3
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_default, 0:libint2_max_am_default, 0:libint2_max_am_3eri), &
      BIND(C) :: libint2_build_3eri
#if INCLUDE_ERI3 >= 1
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_default1, 0:libint2_max_am_default1, 0:libint2_max_am_3eri1), &
      BIND(C) :: libint2_build_3eri1
#endif
#if INCLUDE_ERI3 >= 2
   TYPE(C_FUNPTR), DIMENSION(0:libint2_max_am_default2, 0:libint2_max_am_default2, 0:libint2_max_am_3eri2), &
      BIND(C) :: libint2_build_3eri2
#endif
#endif

   INTERFACE
      SUBROUTINE libint2_static_init() BIND(C)
      END SUBROUTINE

      SUBROUTINE libint2_static_cleanup() BIND(C)
      END SUBROUTINE

#ifdef INCLUDE_ERI
      SUBROUTINE libint2_init_eri(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_eri(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_eri(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_eri
      END FUNCTION

#if INCLUDE_ERI >= 1
      SUBROUTINE libint2_init_eri1(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_eri1(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_eri1(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_eri1
      END FUNCTION
#endif

#if INCLUDE_ERI >= 2
      SUBROUTINE libint2_init_eri2(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_eri2(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_eri2(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_eri2
      END FUNCTION
#endif
#endif

#ifdef INCLUDE_ERI2
      SUBROUTINE libint2_init_2eri(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_2eri(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_2eri(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_2eri
      END FUNCTION

#if INCLUDE_ERI2 >= 1
      SUBROUTINE libint2_init_2eri1(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_2eri1(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_2eri1(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_2eri1
      END FUNCTION
#endif

#if INCLUDE_ERI2 >= 2
      SUBROUTINE libint2_init_2eri2(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_2eri2(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_2eri2(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_2eri2
      END FUNCTION
#endif
#endif

#ifdef INCLUDE_ERI3
      SUBROUTINE libint2_init_3eri(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_3eri(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_3eri(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_3eri
      END FUNCTION

#if INCLUDE_ERI3 >= 1
      SUBROUTINE libint2_init_3eri1(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_3eri1(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_3eri1(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_3eri1
      END FUNCTION
#endif

#if INCLUDE_ERI3 >= 2
      SUBROUTINE libint2_init_3eri2(libint, max_am, buf) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
         INTEGER(KIND=C_INT), VALUE :: max_am
         TYPE(C_PTR), VALUE :: buf
      END SUBROUTINE

      SUBROUTINE libint2_cleanup_3eri2(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE

      FUNCTION libint2_need_memory_3eri2(max_am) BIND(C)
         IMPORT
         INTEGER(KIND=C_INT), VALUE :: max_am
         INTEGER(KIND=C_SIZE_T) :: libint2_need_memory_3eri2
      END FUNCTION
#endif
#endif
   END INTERFACE

   ABSTRACT INTERFACE
      SUBROUTINE libint2_build(libint) BIND(C)
         IMPORT
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE
   END INTERFACE

#else
CONTAINS

  subroutine libint2_init_eri(libint, max_am, buf) bind(c)
     type(libint_t), dimension(*) :: libint
     integer(kind=c_int), value :: max_am
     type(c_ptr), value :: buf
  end subroutine

  subroutine libint2_cleanup_eri(libint) bind(c)
     type(libint_t), dimension(*) :: libint
  end subroutine

  subroutine libint2_build(libint) bind(c)
     type(libint_t), dimension(*) :: libint
  end subroutine

  subroutine libint2_static_init() bind(c)
  end subroutine

  subroutine libint2_static_cleanup() bind(c)
  end subroutine
#endif

END MODULE
