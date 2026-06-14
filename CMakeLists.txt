cmake_minimum_required(VERSION 3.25)

project(Open-Quantum-Package
  VERSION 1.0.0
  LANGUAGES C CXX Fortran)

set(CMAKE_Fortran_STANDARD 2003)

include(FortranCInterface)

option(ENABLE_MPI "Enable MPI support" OFF)
option(ENABLE_DDX "Enable optional ddX continuum-solvation backend" OFF)
option(OQP_REUSE_EXTERNALS
       "Reuse bundled external dependency builds across fresh OpenQP build directories"
       ON)
set(OQP_EXTERNALS_ROOT "" CACHE PATH
    "Optional override for the reusable bundled external dependency cache root.")
if(DEFINED OQP_EXTERNAL_ROOT AND NOT OQP_EXTERNALS_ROOT)
    set(OQP_EXTERNALS_ROOT "${OQP_EXTERNAL_ROOT}" CACHE PATH
        "Optional override for the reusable bundled external dependency cache root." FORCE)
    message(STATUS "OQP_EXTERNAL_ROOT is accepted as an alias for OQP_EXTERNALS_ROOT.")
endif()
mark_as_advanced(OQP_EXTERNALS_ROOT)

if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    message(STATUS "MPI support is enabled.")
endif()

option(BUILD_SHARED_LIBS "Enable shared build"               ON)

option(BUILD_TESTING     "Enable some simple tests"          ON)

option(USE_LIBINT        "Use Libint two-electron code"      ON)

set(LINALG_LIB auto CACHE STRING "Select BLAS & LAPACK Library")
# See BLA_VENDOR in FindBLAS documentation
set(linalg_lib_options none auto netlib
#    Generic
#    ACML ACML_MP ACML_GPU
#    Apple NAS
#    Arm Arm_mp Arm_ilp64 Arm_ilp64_mp
#    ATLAS
#    CXML DXML
#    EML EML_mt
    FLAME
#    FlexiBLAS
#    Fujitsu_SSL2 Fujitsu_SSL2BLAMP Fujitsu_SSL2SVE Fujitsu_SSL2BLAMPSVE
#    Goto
#    IBMESSL IBMESSL_SMP
#    Intel
#    Intel10_32
#    Intel10_64lp
#    Intel10_64lp_seq
    Intel10_64ilp
    Intel10_64ilp_seq
    Intel10_64_dyn
#    NVHPC
    OpenBLAS
#    PhiPACK
#    SCSL SCSL_mp
#    SGIMATH
#    SunPerf
    )
set_property(CACHE LINALG_LIB PROPERTY STRINGS ${linalg_lib_options})

if(NOT LINALG_LIB IN_LIST linalg_lib_options)
    message(STATUS "Error, LINALG_LIB='${LINALG_LIB}'")
    message(FATAL_ERROR "LINALG_LIB must be one of ${linalg_lib_options}")
endif()

option(LINALG_LIB_INT64  "Use 64-bit integer interface to BLAS/LAPACK" ON)
if(NOT LINALG_LIB_INT64)
    message(WARNING
      "LP64 BLAS/LAPACK is no longer supported after removal of defined/defind. "
      "Forcing LINALG_LIB_INT64=ON so stale CMake caches are repaired automatically."
    )
    set(LINALG_LIB_INT64 ON CACHE BOOL
        "Use 64-bit integer interface to BLAS/LAPACK" FORCE)
endif()

option(ENABLE_OPENMP     "Enable OpenMP parallelization"     OFF)

option(ENABLE_PYTHON     "Enable Python Interface"           ON)
option(ENABLE_MSAN       "Enable MemorySanitizer"            OFF) #not used now
option(ENABLE_ASAN       "Enable AddressSanitizer"           OFF)
option(ENABLE_LSAN       "Enable LeakSanitizer"              OFF)
option(ENABLE_UBSAN      "Enable UndefinedBehaviorSanitizer" OFF)
option(ENABLE_TSAN       "Enable ThreadSanitizer"            OFF)

option(ENABLE_Formatter  "Enable source-code formatters"     ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set(CMAKE_BUILD_TYPE Debug)
elseif(${CMAKE_BUILD_TYPE} MATCHES "Release")
  set(CMAKE_BUILD_TYPE Release)
elseif(${CMAKE_BUILD_TYPE} MATCHES "RelWithDebInfo")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
else()
  message(FATAL_ERROR "Unallowed build type!")
endif()

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_compile_options(
  $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-integer-8>
  $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
  )
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>)
  endif()
  if(ENABLE_ASAN)
    add_compile_options(
      -fsanitize=address
      -fsanitize-address-use-after-scope
      -fno-omit-frame-pointer
      -fsanitize=pointer-subtract
      )
    add_link_options(
      -fsanitize=address
      -fsanitize-address-use-after-scope
      -fno-omit-frame-pointer
      -fsanitize=pointer-subtract
      )
  endif()
  if(ENABLE_LSAN)
    add_compile_options(
      -fsanitize=leak
      )
    add_link_options(
      -fsanitize=leak
      )
  endif()
  if(ENABLE_UBSAN)
    add_compile_options(
      -fsanitize=undefined
      )
    add_link_options(
      -fsanitize=undefined
      )
  endif()
  if(ENABLE_TSAN)
    add_compile_options(
      -fsanitize=thread
      )
    add_link_options(
      -fsanitize=thread
      )
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang")
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fdefault-integer-8>)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-i8>)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-i8>)
endif()
set(BLA_SIZEOF_INTEGER 8)
set(OTR_SUFFIX "_64")
set(OTR_DEFS "-DUSE_ILP64")

# Build and link the external OpenTrustRegion (OpenTRAH) SCF solver. When OFF the
# external library is neither downloaded nor compiled, otr_interface.F90 is dropped
# from the build, and the TRAH dispatch falls back to the native solver for every
# trh_impl. Default ON preserves the validated gradient/MRSF reference paths.
option(ENABLE_OPENTRAH "Build and link the external OpenTrustRegion (OpenTRAH) solver" ON)
if(ENABLE_OPENMP)
  find_package(OpenMP 3.1 REQUIRED)
endif()

if(NOT BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_EXE_LINKER_FLAGS "-static")
  endif()
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(oqp_functions)

if(ENABLE_DDX)
  # ddX is LGPL-3.0: it stays an *external*, dynamically-linked dependency
  # (never vendored into the OpenQP tree). Two ways to provide it:
  #   * DDX_ROOT given -> use that prebuilt install (find_package).
  #   * otherwise      -> CMake builds ddX v0.8.0 from source as ILP64 to match
  #                       OpenQP's mandatory 64-bit-integer BLAS (see external/).
  # The ILP64 build is required: a stock LP64 ddX exports an unsuffixed dgemm_
  # that collides with OpenQP's ILP64 dgemm_ on Linux's flat namespace.
  if(DEFINED DDX_ROOT OR DEFINED ENV{DDX_ROOT})
    find_package(DDX REQUIRED)
    set(OQP_DDX_AUTOBUILD OFF)
    message(STATUS "ddX continuum-solvation backend (prebuilt): ${DDX_LIBRARY}")
  else()
    set(OQP_DDX_AUTOBUILD ON)
    message(STATUS "ddX continuum-solvation backend: building ILP64 from source (v0.8.0)")
  endif()
endif()

# Fix possible issues with DYLD_LIBRARY_PATH
if(APPLE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

if(${LINALG_LIB} STREQUAL none)
    message(FATAL_ERROR "LINALG_LIB=none is not supported: OpenQP (and the optional OpenTrustRegion) require BLAS/LAPACK. Use LINALG_LIB=auto, a concrete BLAS/LAPACK vendor, or netlib.")
endif()

# Resolve BLAS/LAPACK before configuring external projects. OpenTrustRegion
# needs the same BLAS integer-size decision as OpenQP; early resolution also
# lets a failed system BLAS search fall back to the bundled NetLib target before
# the external project graph is generated. Linking to oqp happens when
# findLinearAlgebra() is called again after the oqp target exists.
findLinearAlgebra()

add_subdirectory(external)

add_definitions(-DOQP)
add_subdirectory(source)

if(BUILD_TESTING)
  enable_testing()
  if(ENABLE_DDX)
    add_executable(oqp_ddx_link_smoke tests/ddx_link_smoke.c)
    target_link_libraries(oqp_ddx_link_smoke PRIVATE DDX::ddx)
    add_test(NAME oqp_ddx_link_smoke COMMAND oqp_ddx_link_smoke)

    add_executable(oqp_ddx_adapter_smoke
      tests/ddx_adapter_smoke.c
      source/solvent_ddx_adapter.c
    )
    target_compile_definitions(oqp_ddx_adapter_smoke PRIVATE OQP_ENABLE_DDX)
    target_include_directories(oqp_ddx_adapter_smoke PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/source)
    target_link_libraries(oqp_ddx_adapter_smoke PRIVATE DDX::ddx m)
    add_test(NAME oqp_ddx_adapter_smoke COMMAND oqp_ddx_adapter_smoke)

    # When ddX is built from source, the imported DDX::ddx target wraps a future
    # build byproduct; make the smoke executables wait for the ExternalProject.
    if(OQP_DDX_AUTOBUILD)
      add_dependencies(oqp_ddx_link_smoke ddx)
      add_dependencies(oqp_ddx_adapter_smoke ddx)
    endif()
  endif()
endif()

if(USE_LIBINT)
    FortranCInterface_VERIFY(CXX)
endif()

if(ENABLE_PYTHON)
    add_subdirectory(pyoqp)
endif()

install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/oqp.h
    DESTINATION include
)
