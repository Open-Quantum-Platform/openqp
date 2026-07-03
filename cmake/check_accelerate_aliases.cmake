# check_accelerate_aliases.cmake
# ---------------------------------------------------------------------------
# POST-BUILD guard for the Accelerate-ILP64 backend (macOS; OpenQP is ILP64-only).
#
# Fails the build if liboqp still references any *classic* (LP64, 4-byte
# integer) Accelerate BLAS/LAPACK symbol that was NOT interposed onto its
# modern 8-byte `<name>$NEWLAPACK$ILP64` variant. An un-aliased classic symbol
# would execute a 4-byte-integer routine while OpenQP passes 8-byte integers
# => silent numerical corruption. This is what keeps the curated symbol list
# in oqp_functions.cmake honest: adding a new BLAS/LAPACK call without adding
# it to OQP_ACCELERATE_ILP64_SYMS breaks the build here, loudly, with the
# missing name -- instead of shipping corrupted numerics.
#
# Invoked from the Apple-ILP64 link step in findLinearAlgebra():
#     add_custom_command(TARGET oqp POST_BUILD
#       COMMAND ${CMAKE_COMMAND} -DLIB=$<TARGET_FILE:oqp>
#               -P ${CMAKE_SOURCE_DIR}/cmake/check_accelerate_aliases.cmake)
#
# nm -m annotations this relies on (verified on macOS 15/26, Xcode ld):
#   correctly interposed classic symbol:
#       (indirect) external _dgesvd_ (for _dgesvd$NEWLAPACK$ILP64)
#   its ILP64 target (expected, OK):
#       (undefined) external _dgesvd$NEWLAPACK$ILP64 (from Accelerate)
#   a LEAK (classic symbol still bound to legacy LP64 Accelerate):
#       (undefined) external _dgesvd_ (from Accelerate)
#
# The pattern matches ANY lowercase Fortran-style symbol resolved from
# Accelerate (not just [sdcz]-prefixed): integer-returning routines (idamax,
# isamax, ...) and auxiliaries (lsame, xerbla, dlamch) are BLAS/LAPACK too,
# and a leak in any of them is just as corrupting.
# ---------------------------------------------------------------------------

if(NOT DEFINED LIB OR NOT EXISTS "${LIB}")
  message(FATAL_ERROR "check_accelerate_aliases: LIB not set or missing: '${LIB}'")
endif()

find_program(NM_EXE nm)
if(NOT NM_EXE)
  message(FATAL_ERROR "check_accelerate_aliases: 'nm' not found in PATH")
endif()

execute_process(COMMAND "${NM_EXE}" -m "${LIB}"
                OUTPUT_VARIABLE _nm
                RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "check_accelerate_aliases: 'nm -m ${LIB}' failed (rc=${_rc})")
endif()

# The trailing "_ " (underscore then space) excludes the $NEWLAPACK$ILP64
# targets themselves (they end in a digit) and C/CBLAS names (no trailing
# underscore before the annotation).
string(REGEX MATCHALL
       "\\(undefined\\) external _[a-z][a-z0-9]+_ \\(from Accelerate\\)"
       _leak_lines "${_nm}")

set(_bad "")
foreach(_line IN LISTS _leak_lines)
  string(REGEX MATCH "_[a-z][a-z0-9]+_" _sym "${_line}")
  list(APPEND _bad "${_sym}")
endforeach()

if(_bad)
  list(REMOVE_DUPLICATES _bad)
  list(SORT _bad)
  string(REPLACE ";" " " _badstr "${_bad}")
  message(FATAL_ERROR
    "Accelerate ILP64 alias check FAILED for:\n  ${LIB}\n"
    "These classic LP64 BLAS/LAPACK symbols are NOT routed to \$NEWLAPACK\$ILP64:\n"
    "    ${_badstr}\n"
    "Add each (base name, without the leading '_' / trailing '_') to "
    "OQP_ACCELERATE_ILP64_SYMS in cmake/oqp_functions.cmake and rebuild.\n"
    "Rationale: an un-aliased LP64 symbol executes a 4-byte-integer BLAS/LAPACK "
    "routine while OpenQP passes 8-byte integers => silent numerical corruption.")
endif()

message(STATUS "Accelerate ILP64 alias check OK: no LP64 BLAS/LAPACK leaks in ${LIB}")
