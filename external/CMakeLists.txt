include(ExternalProject)

set(_OQP_EXTERNALS_CACHE_REVISION "1")
set(_OQP_LIBINT2_VERSION "2.7.1.1-am4")
set(_OQP_NLOPT_VERSION "2.9.1")
set(_OQP_LIBXC_VERSION "7.0.0")
set(_OQP_TAGARRAY_VERSION "0.0.6")
set(_OQP_LIBECPINT_VERSION "1.0.7")
set(_OQP_LAPACK_VERSION "3.10.0")
set(_OQP_OPENTRUSTREGION_VERSION "2.0.0")
set(_OQP_FPRETTIFY_VERSION "0.3.7")
set(_OQP_DDX_VERSION "0.8.0")
set(_OQP_EXTERNALS_VERSION_KEY
    "ext${_OQP_EXTERNALS_CACHE_REVISION}-libint${_OQP_LIBINT2_VERSION}-nlopt${_OQP_NLOPT_VERSION}-libxc${_OQP_LIBXC_VERSION}-tag${_OQP_TAGARRAY_VERSION}-ecp${_OQP_LIBECPINT_VERSION}-lapack${_OQP_LAPACK_VERSION}-otr${_OQP_OPENTRUSTREGION_VERSION}-fmt${_OQP_FPRETTIFY_VERSION}-ddx${_OQP_DDX_VERSION}")

set(OQP_EXTERNAL_GENERATOR_ARGS
        -DCMAKE_GENERATOR:INTERNAL=${CMAKE_GENERATOR}
)
if(CMAKE_MAKE_PROGRAM)
    list(APPEND OQP_EXTERNAL_GENERATOR_ARGS -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM})
endif()

if(OQP_REUSE_EXTERNALS)
  if(OQP_EXTERNALS_ROOT)
    set(_OQP_EXTERNALS_ROOT "${OQP_EXTERNALS_ROOT}")
  elseif(APPLE AND DEFINED ENV{HOME} AND NOT "$ENV{HOME}" STREQUAL "")
    set(_OQP_EXTERNALS_ROOT "$ENV{HOME}/Library/Caches/openqp/externals")
  elseif(DEFINED ENV{XDG_CACHE_HOME} AND NOT "$ENV{XDG_CACHE_HOME}" STREQUAL "")
    set(_OQP_EXTERNALS_ROOT "$ENV{XDG_CACHE_HOME}/openqp/externals")
  elseif(DEFINED ENV{HOME} AND NOT "$ENV{HOME}" STREQUAL "")
    set(_OQP_EXTERNALS_ROOT "$ENV{HOME}/.cache/openqp/externals")
  else()
    set(_OQP_EXTERNALS_ROOT "${CMAKE_BINARY_DIR}/external-cache")
  endif()
  get_filename_component(_OQP_EXTERNALS_ROOT "${_OQP_EXTERNALS_ROOT}" ABSOLUTE)
  set(_oqp_external_linkage "static")
  if(BUILD_SHARED_LIBS)
    set(_oqp_external_linkage "shared")
  endif()
  set(_oqp_external_generator "${CMAKE_GENERATOR}")
  if(CMAKE_GENERATOR_PLATFORM)
    string(APPEND _oqp_external_generator "-${CMAKE_GENERATOR_PLATFORM}")
  endif()
  if(CMAKE_GENERATOR_TOOLSET)
    string(APPEND _oqp_external_generator "-${CMAKE_GENERATOR_TOOLSET}")
  endif()
  if(CMAKE_GENERATOR_INSTANCE)
    string(APPEND _oqp_external_generator "-${CMAKE_GENERATOR_INSTANCE}")
  endif()
  if(CMAKE_MAKE_PROGRAM)
    string(MD5 _oqp_external_make_program_hash "${CMAKE_MAKE_PROGRAM}")
    string(SUBSTRING "${_oqp_external_make_program_hash}" 0 12 _oqp_external_make_program_hash)
    string(APPEND _oqp_external_generator "-make${_oqp_external_make_program_hash}")
  endif()
  if(CMAKE_CONFIGURATION_TYPES)
    string(REPLACE ";" "+" _oqp_external_configurations "${CMAKE_CONFIGURATION_TYPES}")
    string(APPEND _oqp_external_generator "-${_oqp_external_configurations}")
  endif()
  # User compiler flags are forwarded into the external builds but are not
  # otherwise visible in the key, so hash them in: a build with e.g.
  # -DCMAKE_Fortran_FLAGS="-fcheck=all" must not reuse externals compiled
  # without those flags.
  set(_oqp_external_flag_blob
      "C=${CMAKE_C_FLAGS}|CXX=${CMAKE_CXX_FLAGS}|Fortran=${CMAKE_Fortran_FLAGS}")
  set(_oqp_external_flag_configs "${CMAKE_BUILD_TYPE}")
  if(CMAKE_CONFIGURATION_TYPES)
    set(_oqp_external_flag_configs "${CMAKE_CONFIGURATION_TYPES}")
  endif()
  foreach(_oqp_external_config IN LISTS _oqp_external_flag_configs)
    string(TOUPPER "${_oqp_external_config}" _oqp_external_config_upper)
    string(APPEND _oqp_external_flag_blob
        "|C_${_oqp_external_config_upper}=${CMAKE_C_FLAGS_${_oqp_external_config_upper}}|CXX_${_oqp_external_config_upper}=${CMAKE_CXX_FLAGS_${_oqp_external_config_upper}}|Fortran_${_oqp_external_config_upper}=${CMAKE_Fortran_FLAGS_${_oqp_external_config_upper}}")
  endforeach()
  string(MD5 _oqp_external_flags_hash "${_oqp_external_flag_blob}")
  string(SUBSTRING "${_oqp_external_flags_hash}" 0 12 _oqp_external_flags_hash)
  # _LINALG_LIB_TYPE collapses every non-MKL system BLAS to "other", but some
  # cached externals bake the concrete provider into their artifacts (the
  # shared libddx.so links BLAS/LAPACK at its own build time, and ddX receives
  # BLA_VENDOR), so hash the resolved provider and library paths into the key.
  # The bundled-NetLib route is normalized to a constant: it must produce the
  # same key whether netlib was chosen explicitly or reached by fallback, and
  # any BLAS variables left over from a failed system search must not split it.
  if(_LINALG_LIB_TYPE STREQUAL "NetLib")
    set(_oqp_external_blas_blob "netlib-bundled")
  else()
    set(_oqp_external_blas_blob
        "${LINALG_LIB}|${BLA_VENDOR}|${BLAS_LIBRARIES}|${LAPACK_LIBRARIES}")
  endif()
  string(MD5 _oqp_external_blas_hash "${_oqp_external_blas_blob}")
  string(SUBSTRING "${_oqp_external_blas_hash}" 0 12 _oqp_external_blas_hash)
  set(_oqp_external_key
      "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${_oqp_external_generator}-${CMAKE_C_COMPILER_ID}${CMAKE_C_COMPILER_VERSION}-${CMAKE_CXX_COMPILER_ID}${CMAKE_CXX_COMPILER_VERSION}-${CMAKE_Fortran_COMPILER_ID}${CMAKE_Fortran_COMPILER_VERSION}-${CMAKE_BUILD_TYPE}-${_oqp_external_linkage}-${_LINALG_LIB_TYPE}-blas${_oqp_external_blas_hash}-blasint${BLA_SIZEOF_INTEGER}-flags${_oqp_external_flags_hash}-${_OQP_EXTERNALS_VERSION_KEY}")
  string(REGEX REPLACE "[^A-Za-z0-9_.+-]" "_" _oqp_external_key "${_oqp_external_key}")
  set(OQP_EXTERNALS_SOURCE_ROOT "${_OQP_EXTERNALS_ROOT}/${_oqp_external_key}/src" CACHE INTERNAL "OQP_EXTERNALS_SOURCE_ROOT" FORCE)
  set(OQP_EXTERNALS_BUILD_ROOT "${_OQP_EXTERNALS_ROOT}/${_oqp_external_key}/build" CACHE INTERNAL "OQP_EXTERNALS_BUILD_ROOT" FORCE)
  set(OQP_EXTERNALS_INSTALL_ROOT "${_OQP_EXTERNALS_ROOT}/${_oqp_external_key}/install" CACHE INTERNAL "OQP_EXTERNALS_INSTALL_ROOT" FORCE)
  message(STATUS "Reusable bundled externals cache: ${_OQP_EXTERNALS_ROOT}")
  message(STATUS "Reusable bundled externals key: ${_oqp_external_key}")
endif()

function(oqp_set_external_paths var_prefix name default_prefix default_source default_build default_stamp default_install)
  if(OQP_REUSE_EXTERNALS)
    set(_prefix "${OQP_EXTERNALS_BUILD_ROOT}/${name}")
    set(_source "${OQP_EXTERNALS_SOURCE_ROOT}/${name}")
    set(_build "${_prefix}/build-${CMAKE_Fortran_COMPILER_ID}")
    set(_stamp "${_prefix}/stamp-${CMAKE_Fortran_COMPILER_ID}")
    set(_install "${OQP_EXTERNALS_INSTALL_ROOT}/${name}")
  else()
    set(_prefix "${default_prefix}")
    set(_source "${default_source}")
    set(_build "${default_build}")
    set(_stamp "${default_stamp}")
    set(_install "${default_install}")
  endif()
  set(${var_prefix}_PREFIX "${_prefix}" CACHE INTERNAL "${var_prefix}_PREFIX" FORCE)
  set(${var_prefix}_SOURCE_DIR "${_source}" CACHE INTERNAL "${var_prefix}_SOURCE_DIR" FORCE)
  set(${var_prefix}_BUILD_DIR "${_build}" CACHE INTERNAL "${var_prefix}_BUILD_DIR" FORCE)
  set(${var_prefix}_STAMP_DIR "${_stamp}" CACHE INTERNAL "${var_prefix}_STAMP_DIR" FORCE)
  set(${var_prefix}_INSTALL_DIR "${_install}" CACHE INTERNAL "${var_prefix}_INSTALL_DIR" FORCE)
endfunction()

# Reuse already-built bundled externals across fresh build directories.
#
# The externals are compiled into a persistent per-toolchain cache (see
# oqp_set_external_paths) so their artifacts survive `rm -rf build`. ExternalProject
# alone is not enough to *reuse* them, though: deleting the build directory drops
# build/.ninja_log, and on the next build ninja finds no recorded command line for
# any ExternalProject stamp and conservatively marks every step dirty
# ("command line not found in log"), recompiling all externals despite the cache.
#
# This helper short-circuits that. When OQP_REUSE_EXTERNALS is on and every cached
# artifact listed already exists, it creates a no-op stand-in target with the same
# name as the ExternalProject. oqp consumes the externals through -L/library-name
# and full-path links (not through these targets), which are only used for
# add_dependencies() ordering, so the stand-in keeps add_dependencies() valid while
# the real ExternalProject_Add is skipped (guarded by `if(NOT TARGET <name>)`).
# A cold build, or any incomplete cache, leaves the target undefined and falls
# through to a normal ExternalProject build.
function(oqp_reuse_or_build target)
  if(NOT OQP_REUSE_EXTERNALS)
    return()
  endif()
  foreach(_artifact IN LISTS ARGN)
    if(NOT EXISTS "${_artifact}")
      return()
    endif()
  endforeach()
  add_custom_target(${target})
  message(STATUS "Reusing cached external '${target}' (skipping rebuild)")
endfunction()

set(OQP_EXTERNAL_LIST_SEPARATOR "|")
function(oqp_external_cmake_list_arg out_var cache_var)
  set(_value "${ARGN}")
  string(REPLACE ";" "${OQP_EXTERNAL_LIST_SEPARATOR}" _value "${_value}")
  set(${out_var} "-D${cache_var}:STRING=${_value}" PARENT_SCOPE)
endfunction()

if(${USE_LIBINT})
oqp_set_external_paths(LIBINT2 libint
    "${CMAKE_BINARY_DIR}/external/libint"
    "${CMAKE_SOURCE_DIR}/external/libint"
    "${CMAKE_BINARY_DIR}/external/libint/build-${CMAKE_Fortran_COMPILER_ID}"
    "${CMAKE_BINARY_DIR}/external/libint/stamp-${CMAKE_Fortran_COMPILER_ID}"
    "${CMAKE_BINARY_DIR}/external/install/libint")
set(LIBINT2_LIBRARY "${LIBINT2_BUILD_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}int2${CMAKE_STATIC_LIBRARY_SUFFIX}" CACHE INTERNAL "LIBINT2_LIBRARY" FORCE)

oqp_reuse_or_build(libint2 "${LIBINT2_LIBRARY}" "${LIBINT2_BUILD_DIR}/fortran" "${LIBINT2_BUILD_DIR}/include")
if(NOT TARGET libint2)
ExternalProject_Add(libint2
    PREFIX ${LIBINT2_PREFIX}
    URL https://qchemlab.knu.ac.kr/vmironov/oqp-libint/-/archive/v${_OQP_LIBINT2_VERSION}/oqp-libint-v${_OQP_LIBINT2_VERSION}.tar.gz
    TLS_VERIFY TRUE
    SOURCE_DIR ${LIBINT2_SOURCE_DIR}
    BINARY_DIR ${LIBINT2_BUILD_DIR}
    STAMP_DIR  ${LIBINT2_STAMP_DIR}
    BUILD_BYPRODUCTS ${LIBINT2_LIBRARY}
    UPDATE_COMMAND ""
    INSTALL_COMMAND ""
    CMAKE_ARGS -DENABLE_FORTRAN=ON
               -DREQUIRE_CXX_API=OFF
               -DENABLE_MPFR=OFF
               -DLIBINT2_SHGAUSS_ORDERING=standard
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
               -DLIBINT2_BUILD_SHARED_AND_STATIC_LIBS=${BUILD_SHARED_LIBS}
               ${OQP_EXTERNAL_GENERATOR_ARGS}
)
endif()
endif()

oqp_set_external_paths(NLOPT nlopt
        "${CMAKE_BINARY_DIR}/external/nlopt"
        "${CMAKE_BINARY_DIR}/external/nlopt/src/nlopt"
        "${CMAKE_BINARY_DIR}/external/nlopt/build-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/nlopt/stamp-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/install/nlopt")
set(NLOPT_LIBRARY "${NLOPT_BUILD_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}nlopt${CMAKE_STATIC_LIBRARY_SUFFIX}" CACHE INTERNAL "NLOPT_LIBRARY" FORCE)
oqp_reuse_or_build(nlopt "${NLOPT_LIBRARY}" "${NLOPT_BUILD_DIR}/nlopt.f")
if(NOT TARGET nlopt)
ExternalProject_Add(nlopt
        PREFIX ${NLOPT_PREFIX}
        URL https://github.com/stevengj/nlopt/archive/v${_OQP_NLOPT_VERSION}.tar.gz
        URL_HASH SHA256=1e6c33f8cbdc4138d525f3326c231f14ed50d99345561e85285638c49b64ee93
        SOURCE_DIR ${NLOPT_SOURCE_DIR}
        BINARY_DIR ${NLOPT_BUILD_DIR}
        STAMP_DIR  ${NLOPT_STAMP_DIR}
        BUILD_BYPRODUCTS ${NLOPT_LIBRARY}
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
                   -DNLOPT_CXX=OFF
                   -DNLOPT_FORTRAN=ON
                   -DNLOPT_GUILE=OFF
                   -DNLOPT_MATLAB=OFF
                   -DNLOPT_OCTAVE=OFF
                   -DNLOPT_PYTHON=OFF
                   -DNLOPT_SWIG=OFF
                   -DNLOPT_TESTS=OFF
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                   ${OQP_EXTERNAL_GENERATOR_ARGS}
        )
endif()

set(LIBXC_CMAKE_ARGS
        -DDISABLE_VXC=OFF
        -DDISABLE_FXC=OFF
        -DDISABLE_KXC=OFF
        -DENABLE_FORTRAN=ON
        -DBUILD_TESTING=OFF
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_POLICY_VERSION_MINIMUM=3.25
        ${OQP_EXTERNAL_GENERATOR_ARGS}
)
set(OQP_LIBXC_C_FLAGS_RELEASE "" CACHE STRING "Optional release C flags for external libxc")
if(OQP_LIBXC_C_FLAGS_RELEASE)
    list(APPEND LIBXC_CMAKE_ARGS -DCMAKE_C_FLAGS_RELEASE=${OQP_LIBXC_C_FLAGS_RELEASE})
endif()

oqp_set_external_paths(LIBXC libxc
        "${CMAKE_SOURCE_DIR}/external/libxc"
        "${CMAKE_SOURCE_DIR}/external/libxc/src/libxc"
        "${CMAKE_BINARY_DIR}/external/libxc/src/libxc-build-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/libxc/src/libxc-stamp-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/install/libxc")
set(LIBXC_LIBRARY_DIR "${LIBXC_BUILD_DIR}" CACHE INTERNAL "LIBXC_LIBRARY_DIR" FORCE)
set(LIBXC_INCLUDE_DIR "${LIBXC_BUILD_DIR}" CACHE INTERNAL "LIBXC_INCLUDE_DIR" FORCE)
oqp_reuse_or_build(libxc
        "${LIBXC_BUILD_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}xc${CMAKE_STATIC_LIBRARY_SUFFIX}"
        "${LIBXC_BUILD_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}xcf03${CMAKE_STATIC_LIBRARY_SUFFIX}"
        "${LIBXC_BUILD_DIR}/xc_f03_lib_m.mod")
if(NOT TARGET libxc)
ExternalProject_Add(libxc
        PREFIX ${LIBXC_PREFIX}
	URL https://gitlab.com/libxc/libxc/-/archive/${_OQP_LIBXC_VERSION}/libxc-${_OQP_LIBXC_VERSION}.tar.gz
        URL_HASH SHA256=8d4e343041c9cd869833822f57744872076ae709a613c118d70605539fb13a77
        SOURCE_DIR ${LIBXC_SOURCE_DIR}
        BINARY_DIR ${LIBXC_BUILD_DIR}
        STAMP_DIR  ${LIBXC_STAMP_DIR}
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        CMAKE_ARGS ${LIBXC_CMAKE_ARGS}
        )
endif()

oqp_set_external_paths(TAGARRAY libtagarray
        "${CMAKE_SOURCE_DIR}/external/libtagarray"
        "${CMAKE_SOURCE_DIR}/external/libtagarray/src/libtagarray"
        "${CMAKE_BINARY_DIR}/external/libtagarray/src/libtagarray-build-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/libtagarray/src/libtagarray-stamp-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/install/libtagarray")
set(TAGARRAY_LIBRARY_DIR "${TAGARRAY_BUILD_DIR}" CACHE INTERNAL "TAGARRAY_LIBRARY_DIR" FORCE)
set(TAGARRAY_INCLUDE_DIRS
    "${TAGARRAY_BUILD_DIR}/include/tagarray;${TAGARRAY_SOURCE_DIR}/include"
    CACHE INTERNAL "TAGARRAY_INCLUDE_DIRS" FORCE)
oqp_reuse_or_build(libtagarray
        "${TAGARRAY_BUILD_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}tagarray${CMAKE_STATIC_LIBRARY_SUFFIX}"
        "${TAGARRAY_BUILD_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}tagarray_f${CMAKE_STATIC_LIBRARY_SUFFIX}"
        "${TAGARRAY_BUILD_DIR}/include/tagarray")
if(NOT TARGET libtagarray)
ExternalProject_Add(libtagarray
        PREFIX ${TAGARRAY_PREFIX}
        URL https://github.com/foxtran/libtagarray/archive/refs/tags/v${_OQP_TAGARRAY_VERSION}.tar.gz
        URL_HASH SHA256=9c844092f123276bb461e37e3451068c097b5a1fec9d65514ad3018baaeb2b32
#        GIT_REPOSITORY https://github.com/foxtran/libtagarray.git
#        GIT_TAG master
        SOURCE_DIR ${TAGARRAY_SOURCE_DIR}
        BINARY_DIR ${TAGARRAY_BUILD_DIR}
        STAMP_DIR  ${TAGARRAY_STAMP_DIR}
        UPDATE_COMMAND ""
        PATCH_COMMAND /bin/sh -c "grep -q reserve_data_default <SOURCE_DIR>/source/API/Fortran/container.f90 || patch -p1 -i ${CMAKE_SOURCE_DIR}/cmake/patches/libtagarray-v0.0.6-default-integer-shapes.patch"
        INSTALL_COMMAND ""
        CMAKE_ARGS -DENABLE_FORTRAN=ON
                   -DENABLE_TESTING=OFF
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                   ${OQP_EXTERNAL_GENERATOR_ARGS}
        )
endif()

oqp_set_external_paths(LIBECPINT libecpint
        "${CMAKE_SOURCE_DIR}/external/libecpint"
        "${CMAKE_SOURCE_DIR}/external/libecpint/src/libecpint"
        "${CMAKE_BINARY_DIR}/external/libecpint/src/libecpint-build-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/libecpint/src/libecpint-stamp-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/libecpint/install")
set(LIBECPINT_LIBRARY_DIR "${LIBECPINT_BUILD_DIR}/src" CACHE INTERNAL "LIBECPINT_LIBRARY_DIR" FORCE)
set(FADDEEVA_LIBRARY_DIR "${LIBECPINT_BUILD_DIR}/external/Faddeeva" CACHE INTERNAL "FADDEEVA_LIBRARY_DIR" FORCE)
set(LIBECPINT_INCLUDE_DIRS
    "${LIBECPINT_INSTALL_DIR}/include/libecpint;${LIBECPINT_INSTALL_DIR}/include"
    CACHE INTERNAL "LIBECPINT_INCLUDE_DIRS" FORCE)
oqp_reuse_or_build(libecpint
        "${LIBECPINT_LIBRARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}ecpint${CMAKE_STATIC_LIBRARY_SUFFIX}"
        "${FADDEEVA_LIBRARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}Faddeeva${CMAKE_STATIC_LIBRARY_SUFFIX}"
        "${LIBECPINT_INSTALL_DIR}/include/libecpint")
if(NOT TARGET libecpint)
ExternalProject_Add(libecpint
        PREFIX ${LIBECPINT_PREFIX}
        URL https://github.com/robashaw/libecpint/archive/refs/tags/v${_OQP_LIBECPINT_VERSION}.tar.gz
        URL_HASH SHA256=e9c60fddb2614f113ab59ec620799d961db73979845e6e637c4a6fb72aee51cc
        SOURCE_DIR ${LIBECPINT_SOURCE_DIR}
        BINARY_DIR ${LIBECPINT_BUILD_DIR}
        STAMP_DIR  ${LIBECPINT_STAMP_DIR}
	INSTALL_DIR ${LIBECPINT_INSTALL_DIR}
        UPDATE_COMMAND ""
        CMAKE_ARGS -DLIBECPINT_BUILD_DOCS=OFF
                   -DLIBECPINT_BUILD_TESTS=OFF
		   -DLIBECPINT_USE_PUGIXML=OFF
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                   ${OQP_EXTERNAL_GENERATOR_ARGS}
		   -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        INSTALL_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target install
        )
endif()
if(_LINALG_LIB_TYPE STREQUAL NetLib)
  if(LINALG_LIB_INT64)
    set(BLA_LIB_SUFFIX "64")
  else()
    set(BLA_LIB_SUFFIX "")
  endif()
  if(BUILD_SHARED_LIBS)
    set(la_fortran_flags "${CMAKE_Fortran_FLAGS} -fpic")
  else()
    set(la_fortran_flags ${CMAKE_Fortran_FLAGS})
  endif()
  oqp_set_external_paths(LAPACK lapack
        "${CMAKE_SOURCE_DIR}/external/lapack"
        "${CMAKE_SOURCE_DIR}/external/lapack/src/LAPACK"
        "${CMAKE_BINARY_DIR}/external/lapack/src/lapack-build"
        "${CMAKE_BINARY_DIR}/external/lapack/src/lapack-stamp"
        "${CMAKE_BINARY_DIR}/external/install/lapack")
  set(LIBLAPACK ${LAPACK_BUILD_DIR}/lib/liblapack${BLA_LIB_SUFFIX}.a CACHE INTERNAL "LIBLAPACK" FORCE)

  set(LIBBLAS   ${LAPACK_BUILD_DIR}/lib/libblas${BLA_LIB_SUFFIX}.a CACHE INTERNAL "LIBBLAS" FORCE)
  message(STATUS "BLAS_SUFFIX = ${BLA_LIB_SUFFIX}")

  oqp_reuse_or_build(LAPACK "${LIBLAPACK}" "${LIBBLAS}")
  if(NOT TARGET LAPACK)
  ExternalProject_Add(LAPACK
        PREFIX ${LAPACK_PREFIX}
        URL https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v${_OQP_LAPACK_VERSION}.tar.gz
        URL_HASH SHA256=328c1bea493a32cac5257d84157dc686cc3ab0b004e2bea22044e0a59f6f8a19
        SOURCE_DIR ${LAPACK_SOURCE_DIR}
        BINARY_DIR ${LAPACK_BUILD_DIR}
        STAMP_DIR  ${LAPACK_STAMP_DIR}
        BUILD_BYPRODUCTS ${LIBLAPACK}
                         ${LIBBLAS}
        UPDATE_COMMAND ""
        PATCH_COMMAND /usr/bin/perl -0pi -e "s/cmake_minimum_required\\(VERSION 3\\.2\\)/cmake_minimum_required(VERSION 3.5)/g" <SOURCE_DIR>/CMakeLists.txt <SOURCE_DIR>/INSTALL/CMakeLists.txt
        INSTALL_COMMAND ""
        CMAKE_ARGS -DBUILD_INDEX64=${LINALG_LIB_INT64}
                   -DBUILD_SINGLE=ON
                   -DBUILD_DOUBLE=ON
                   -DBUILD_COMPLEX=ON
                   -DBUILD_COMPLEX16=ON
                   -DUSE_OPTIMIZED_BLAS=OFF
                   -DUSE_OPTIMIZED_LAPACK=OFF
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                   -DCMAKE_Fortran_FLAGS=${la_fortran_flags}
                   -DCMAKE_POLICY_VERSION_MINIMUM=3.5
                   -DCMAKE_TRY_COMPILE_PLATFORM_VARIABLES=CMAKE_POLICY_VERSION_MINIMUM
                   ${OQP_EXTERNAL_GENERATOR_ARGS}
        )
  endif()
endif()
if(ENABLE_OPENTRAH)
oqp_set_external_paths(OPENTRUSTREGION opentrustregion
    "${CMAKE_SOURCE_DIR}/external/opentrustregion"
    "${CMAKE_SOURCE_DIR}/external/opentrustregion/src/libopentrustregion"
    "${CMAKE_BINARY_DIR}/external/opentrustregion/src/libopentrustregion-build-${CMAKE_Fortran_COMPILER_ID}"
    "${CMAKE_BINARY_DIR}/external/opentrustregion/src/libopentrustregion-stamp-${CMAKE_Fortran_COMPILER_ID}"
    "${CMAKE_BINARY_DIR}/external/install/opentrustregion")
set(otr_integer_flags "")
if(LINALG_LIB_INT64)
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
        set(otr_integer_flags "-fdefault-integer-8")
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang")
        set(otr_integer_flags "-fdefault-integer-8")
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
        set(otr_integer_flags "-i8")
    endif()
endif()
set(otr_fortran_flags "${CMAKE_Fortran_FLAGS} ${OTR_DEFS} ${otr_integer_flags}")
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
    set(otr_fortran_flags "${otr_fortran_flags} -fallow-argument-mismatch")
endif()
set(OTR_CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${OPENTRUSTREGION_INSTALL_DIR}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBLA_SIZEOF_INTEGER=${BLA_SIZEOF_INTEGER}
        -DINTEGER_SIZE=${BLA_SIZEOF_INTEGER}
        -DBUILD_SHARED_LIBS=OFF
        -DBUILD_TESTS=OFF
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_Fortran_FLAGS=${otr_fortran_flags}
        ${OQP_EXTERNAL_GENERATOR_ARGS}
)
if(CMAKE_PREFIX_PATH)
    oqp_external_cmake_list_arg(_otr_prefix_path_arg CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
    list(APPEND OTR_CMAKE_ARGS "${_otr_prefix_path_arg}")
endif()
if(_LINALG_LIB_TYPE STREQUAL NetLib)
    oqp_external_cmake_list_arg(_otr_blas_arg BLAS_LIBRARIES ${LIBBLAS})
    oqp_external_cmake_list_arg(_otr_lapack_arg LAPACK_LIBRARIES ${LIBLAPACK})
    list(APPEND OTR_CMAKE_ARGS
        "${_otr_blas_arg}"
        "${_otr_lapack_arg}"
    )
elseif(DEFINED BLAS_LIBRARIES AND DEFINED LAPACK_LIBRARIES)
    oqp_external_cmake_list_arg(_otr_blas_arg BLAS_LIBRARIES ${BLAS_LIBRARIES})
    oqp_external_cmake_list_arg(_otr_lapack_arg LAPACK_LIBRARIES ${LAPACK_LIBRARIES})
    list(APPEND OTR_CMAKE_ARGS
        "${_otr_blas_arg}"
        "${_otr_lapack_arg}"
    )
elseif(DEFINED BLA_VENDOR)
    list(APPEND OTR_CMAKE_ARGS -DBLA_VENDOR=${BLA_VENDOR})
    if(DEFINED BLA_PREFER_PKGCONFIG)
        list(APPEND OTR_CMAKE_ARGS -DBLA_PREFER_PKGCONFIG=${BLA_PREFER_PKGCONFIG})
    endif()
    if(DEFINED BLA_PKGCONFIG_BLAS)
        list(APPEND OTR_CMAKE_ARGS -DBLA_PKGCONFIG_BLAS=${BLA_PKGCONFIG_BLAS})
    endif()
    if(DEFINED BLA_PKGCONFIG_LAPACK)
        list(APPEND OTR_CMAKE_ARGS -DBLA_PKGCONFIG_LAPACK=${BLA_PKGCONFIG_LAPACK})
    endif()
endif()

# Portability: the bundled OpenTrustRegion hardcodes -march=native/-xHost in its
# own CMakeLists (CMAKE_Fortran_FLAGS_RELEASE). Two problems with leaving it:
#   * -march=native on an AVX-512 build host emits AVX-512 into level_shifted_davidson,
#     which SIGILLs (exit 132) on AVX2-only compute nodes (the original bug report);
#   * a newer gfortran driver resolves -march=native to a microarch *level* name
#     (e.g. x86-64-v3) that an OLDER backend rejects -- "bad value ('x86-64-v3') for
#     '-march=' switch" -- breaking the Docker build.
# Strip the host-tuning flags entirely so OTR builds with a generic -O2: portable
# across toolchains (old/new gfortran, x86, aarch64) and free of the SIGILL. OTR is a
# small solver, so losing -march tuning is negligible; main OpenQP uses no -march.
set(_OTR_ARCH_SED "s/-march=native//g; s/-xHost//g")
oqp_reuse_or_build(libopentrustregion
    "${OPENTRUSTREGION_BUILD_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}opentrustregion${OTR_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")
if(NOT TARGET libopentrustregion)
ExternalProject_Add(libopentrustregion
    PREFIX ${OPENTRUSTREGION_PREFIX}
    URL https://github.com/eriksen-lab/opentrustregion/archive/refs/tags/v${_OQP_OPENTRUSTREGION_VERSION}.tar.gz
    URL_HASH SHA256=682201cc9cc850a892d132e5e5fc11f6d87f8cd843fdbad88dafbf7277407f1c
    SOURCE_DIR ${OPENTRUSTREGION_SOURCE_DIR}
    UPDATE_COMMAND ""
    PATCH_COMMAND /usr/bin/perl -0pi -e "${_OTR_ARCH_SED}" <SOURCE_DIR>/CMakeLists.txt
    LIST_SEPARATOR ${OQP_EXTERNAL_LIST_SEPARATOR}
    CMAKE_ARGS ${OTR_CMAKE_ARGS}
    BINARY_DIR ${OPENTRUSTREGION_BUILD_DIR}
    STAMP_DIR  ${OPENTRUSTREGION_STAMP_DIR}
)
endif()
if(_LINALG_LIB_TYPE STREQUAL NetLib)
  add_dependencies(libopentrustregion LAPACK)
endif()
endif() # ENABLE_OPENTRAH

set(OTR_LIB_NAME "${CMAKE_STATIC_LIBRARY_PREFIX}opentrustregion${OTR_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")
if(ENABLE_Formatter)
  find_package(Python3 REQUIRED)
  if(OQP_REUSE_EXTERNALS)
    set(FPRETTIFY_PREFIX "${OQP_EXTERNALS_BUILD_ROOT}/fprettify" CACHE INTERNAL "FPRETTIFY_PREFIX" FORCE)
    set(FPRETTIFY_SOURCE_DIR "${OQP_EXTERNALS_SOURCE_ROOT}/fprettify" CACHE INTERNAL "FPRETTIFY_SOURCE_DIR" FORCE)
  else()
    set(FPRETTIFY_PREFIX "${CMAKE_SOURCE_DIR}/external/fprettify" CACHE INTERNAL "FPRETTIFY_PREFIX" FORCE)
    set(FPRETTIFY_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/fprettify/src/fprettify" CACHE INTERNAL "FPRETTIFY_SOURCE_DIR" FORCE)
  endif()
  oqp_reuse_or_build(fprettify "${FPRETTIFY_SOURCE_DIR}/fprettify.py")
  if(NOT TARGET fprettify)
  ExternalProject_Add(fprettify
          PREFIX ${FPRETTIFY_PREFIX}
          URL https://github.com/pseewald/fprettify/archive/refs/tags/v${_OQP_FPRETTIFY_VERSION}.tar.gz
          URL_HASH SHA256=052da19a9080a6641d3202e10572cf3d978e6bcc0e7db29c1eb8ba724e89adc7
          SOURCE_DIR ${FPRETTIFY_SOURCE_DIR}
          BINARY_DIR ""
          STAMP_DIR ""
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          INSTALL_COMMAND ""
          )
  endif()
  add_custom_target(format ${Python3_EXECUTABLE} ${FPRETTIFY_SOURCE_DIR}/fprettify.py --indent 2
                                                 --whitespace 1
                                                 --c-relations
                                                 --enable-replacements
                                                 --enable-decl
                                                 --whitespace-decl
                                                 --case 1 1 1 1
                                                 -r ${CMAKE_SOURCE_DIR}/tests ${CMAKE_SOURCE_DIR}/modules
                                                 -e int1.F90
                                                 -e int_gms_rotaxis.F90
                                                 -e grd1.F90
                                                 -e grd2.F90
                                                 -e opt.F90
                                                 -e dft_rhf_energy.f90
                                                 -e dft_urohf_energy.f90
                                                 -e rhf_energy.f90
                                                 -e rhf_JK.f90
                                                 -e dftgrid.F90
                                                 -e functionals.F90
                                                 -e guess.F90
                                                 -e constants.F90
                                                 -e svpleb.F90
                                                 -e polm_band.f90
                                                 -e hsandt.f90
                   )
  add_dependencies(format fprettify)
endif()



# ---------------------------------------------------------------------------
# ddX continuum-solvation backend (optional, ENABLE_DDX).
#
# ddX is LGPL-3.0 and stays an external, dynamically-linked dependency: its
# source is never vendored into OpenQP. When no prebuilt DDX_ROOT is supplied we
# build ddX v0.8.0 here with the same BLAS integer width as OpenQP.  ILP64
# builds use -fdefault-integer-8/-i8 and a small C-bridge patch that keeps
# ddX's int32 C ABI while letting the core be int64.
# ---------------------------------------------------------------------------
if(ENABLE_DDX AND OQP_DDX_AUTOBUILD)
    oqp_set_external_paths(DDX ddx
        "${CMAKE_BINARY_DIR}/external/ddx"
        "${CMAKE_BINARY_DIR}/external/ddx/src/ddx"
        "${CMAKE_BINARY_DIR}/external/ddx/build-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/ddx/stamp-${CMAKE_Fortran_COMPILER_ID}"
        "${CMAKE_BINARY_DIR}/external/install/ddx")
    set(DDX_INSTALL ${DDX_INSTALL_DIR} CACHE INTERNAL "DDX_INSTALL" FORCE)
    set(DDX_INCLUDE_DIR ${DDX_INSTALL}/include CACHE INTERNAL "DDX_INCLUDE_DIR")
    set(DDX_LIBRARY
        ${DDX_INSTALL}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ddx${CMAKE_SHARED_LIBRARY_SUFFIX}
        CACHE INTERNAL "DDX_LIBRARY" FORCE)

    # Match OpenTrustRegion's BLAS integer flags.
    set(ddx_integer_flags "")
    if(LINALG_LIB_INT64)
        if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang")
            set(ddx_integer_flags "-fdefault-integer-8")
        elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
            set(ddx_integer_flags "-i8")
        endif()
    endif()
    set(ddx_fortran_flags "${CMAKE_Fortran_FLAGS} ${ddx_integer_flags}")
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
        set(ddx_fortran_flags "${ddx_fortran_flags} -fallow-argument-mismatch")
    endif()

    set(DDX_CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${DDX_INSTALL}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBUILD_SHARED_LIBS=ON
        -DTESTS=OFF
        -DEXAMPLES=OFF
        -DWARNING_FLAGS=OFF
        -DBLA_SIZEOF_INTEGER=${BLA_SIZEOF_INTEGER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_Fortran_FLAGS=${ddx_fortran_flags}
        ${OQP_EXTERNAL_GENERATOR_ARGS}
    )
    # Forward CMAKE_PREFIX_PATH only when set, so a non-standard ILP64 OpenBLAS
    # location reaches the ddX subbuild. An empty -D would override the
    # environment CMAKE_PREFIX_PATH the ExternalProject otherwise inherits.
    if(CMAKE_PREFIX_PATH)
        oqp_external_cmake_list_arg(_ddx_prefix_path_arg CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
        list(APPEND DDX_CMAKE_ARGS "${_ddx_prefix_path_arg}")
    endif()
    # Forward the same BLAS/LAPACK choice OpenQP resolved (findLinearAlgebra()).
    if(_LINALG_LIB_TYPE STREQUAL NetLib)
        oqp_external_cmake_list_arg(_ddx_blas_arg BLAS_LIBRARIES ${LIBBLAS})
        oqp_external_cmake_list_arg(_ddx_lapack_arg LAPACK_LIBRARIES ${LIBLAPACK})
        list(APPEND DDX_CMAKE_ARGS
            "${_ddx_blas_arg}"
            "${_ddx_lapack_arg}")
    elseif(DEFINED BLAS_LIBRARIES AND DEFINED LAPACK_LIBRARIES)
        oqp_external_cmake_list_arg(_ddx_blas_arg BLAS_LIBRARIES ${BLAS_LIBRARIES})
        oqp_external_cmake_list_arg(_ddx_lapack_arg LAPACK_LIBRARIES ${LAPACK_LIBRARIES})
        list(APPEND DDX_CMAKE_ARGS
            "${_ddx_blas_arg}"
            "${_ddx_lapack_arg}")
    elseif(DEFINED BLA_VENDOR)
        list(APPEND DDX_CMAKE_ARGS -DBLA_VENDOR=${BLA_VENDOR})
        if(DEFINED BLA_PREFER_PKGCONFIG)
            list(APPEND DDX_CMAKE_ARGS -DBLA_PREFER_PKGCONFIG=${BLA_PREFER_PKGCONFIG})
        endif()
        if(DEFINED BLA_PKGCONFIG_BLAS)
            list(APPEND DDX_CMAKE_ARGS -DBLA_PKGCONFIG_BLAS=${BLA_PKGCONFIG_BLAS})
        endif()
        if(DEFINED BLA_PKGCONFIG_LAPACK)
            list(APPEND DDX_CMAKE_ARGS -DBLA_PKGCONFIG_LAPACK=${BLA_PKGCONFIG_LAPACK})
        endif()
    endif()

    oqp_reuse_or_build(ddx "${DDX_LIBRARY}" "${DDX_INCLUDE_DIR}")
    if(NOT TARGET ddx)
    ExternalProject_Add(ddx
        PREFIX ${DDX_PREFIX}
        URL https://github.com/ddsolvation/ddX/archive/refs/tags/v${_OQP_DDX_VERSION}.tar.gz
        URL_HASH SHA256=6f6d08c7699dd69d4b12dd0b1f8c519f5f3f84a3d23cc526dd1a4a60914969d6
        TLS_VERIFY TRUE
        LIST_SEPARATOR ${OQP_EXTERNAL_LIST_SEPARATOR}
        SOURCE_DIR ${DDX_SOURCE_DIR}
        PATCH_COMMAND patch -p1 -i ${CMAKE_SOURCE_DIR}/cmake/patches/ddx-v0.8.0-ilp64.patch
        BINARY_DIR ${DDX_BUILD_DIR}
        STAMP_DIR  ${DDX_STAMP_DIR}
        BUILD_BYPRODUCTS ${DDX_LIBRARY}
        UPDATE_COMMAND ""
        CMAKE_ARGS ${DDX_CMAKE_ARGS}
    )
    endif()
    if(_LINALG_LIB_TYPE STREQUAL NetLib)
        add_dependencies(ddx LAPACK)
    endif()

    # Imported target consumed by oqp and the smoke tests (see add_dependencies
    # on the `ddx` ExternalProject at each consumer). The include dir must exist
    # at configure time for INTERFACE_INCLUDE_DIRECTORIES.
    file(MAKE_DIRECTORY ${DDX_INCLUDE_DIR})
    add_library(DDX::ddx UNKNOWN IMPORTED GLOBAL)
    set_target_properties(DDX::ddx PROPERTIES
        IMPORTED_LOCATION ${DDX_LIBRARY}
        INTERFACE_INCLUDE_DIRECTORIES ${DDX_INCLUDE_DIR})
    message(STATUS "ddX will be built with ${BLA_SIZEOF_INTEGER}-byte BLAS integers from source: ${DDX_LIBRARY}")
endif()


set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)
