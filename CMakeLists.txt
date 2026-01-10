cmake_minimum_required(VERSION 3.25)

project(Open-Quantum-Package
  VERSION 1.0.0
  LANGUAGES C CXX Fortran)

set(CMAKE_Fortran_STANDARD 2003)

include(FortranCInterface)

option(ENABLE_MPI "Enable MPI support" OFF)

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
if(LINALG_LIB_INT64)
    set(BLA_SIZEOF_INTEGER 8)
    set(OTR_SUFFIX "_64")
    set(OTR_DEFS "-DUSE_ILP64")
else()
    set(BLA_SIZEOF_INTEGER 4)
    set(OTR_SUFFIX "_32")
    set(OTR_DEFS "")
endif()
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

# Fix possible issues with DYLD_LIBRARY_PATH
if(APPLE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

add_subdirectory(external)

add_definitions(-DOQP)
add_subdirectory(source)

if(BUILD_TESTING)
  enable_testing()
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
