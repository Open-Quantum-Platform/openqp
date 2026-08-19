cmake_minimum_required(VERSION 3.25)

# macOS: a bare configure picks Apple clang for C/C++ while Fortran must be
# GNU gfortran (Apple ships no Fortran compiler). Mixing the toolchains links
# BOTH OpenMP runtimes (clang's libomp + gfortran's libgomp) -- crashes and
# thread oversubscription -- and the validated macOS configuration is all-GNU
# (Homebrew gcc). When the user chose no compilers, auto-select the newest
# COMPLETE Homebrew GCC toolchain so a plain `pip install .` gets the
# supported compilers. The trio is picked as a matched set -- never a
# mixed-version gcc/gfortran pair.
#
# Two traps encoded here:
#  * this must run BEFORE project() (compilers are locked in there), and the
#    APPLE variable is only set BY project() -- pre-project() code must test
#    CMAKE_HOST_APPLE instead;
#  * pip's build frontend on macOS INJECTS CC=clang / CXX=clang++ into the
#    build environment. A bare unversioned default (clang/clang++/cc/c++) is
#    toolchain noise, not a user choice, so it does not suppress the
#    auto-selection. Real intent still always wins: any -DCMAKE_*_COMPILER,
#    any FC, or any non-default CC/CXX (e.g. CC=gcc-14 or an absolute path).
set(_oqp_env_cc  "$ENV{CC}")
set(_oqp_env_cxx "$ENV{CXX}")
if(CMAKE_HOST_APPLE AND _oqp_env_cc MATCHES "^(clang|cc)$")
  set(_oqp_env_cc "")
endif()
if(CMAKE_HOST_APPLE AND _oqp_env_cxx MATCHES "^(clang\\+\\+|c\\+\\+)$")
  set(_oqp_env_cxx "")
endif()
if(CMAKE_HOST_APPLE
   AND NOT DEFINED CMAKE_C_COMPILER
   AND NOT DEFINED CMAKE_CXX_COMPILER
   AND NOT DEFINED CMAKE_Fortran_COMPILER
   AND "${_oqp_env_cc}" STREQUAL ""
   AND "${_oqp_env_cxx}" STREQUAL ""
   AND "$ENV{FC}" STREQUAL "")
  foreach(_oqp_gcc_ver 16 15 14 13)
    find_program(_oqp_mac_gcc "gcc-${_oqp_gcc_ver}")
    find_program(_oqp_mac_gxx "g++-${_oqp_gcc_ver}")
    find_program(_oqp_mac_gfc "gfortran-${_oqp_gcc_ver}")
    if(_oqp_mac_gcc AND _oqp_mac_gxx AND _oqp_mac_gfc)
      set(CMAKE_C_COMPILER       "${_oqp_mac_gcc}")
      set(CMAKE_CXX_COMPILER     "${_oqp_mac_gxx}")
      set(CMAKE_Fortran_COMPILER "${_oqp_mac_gfc}")
      message(STATUS "macOS: auto-selected Homebrew GCC ${_oqp_gcc_ver} toolchain "
                     "(${_oqp_mac_gcc}); set CC/CXX/FC or -DCMAKE_*_COMPILER to override")
      break()
    endif()
    unset(_oqp_mac_gcc CACHE)
    unset(_oqp_mac_gxx CACHE)
    unset(_oqp_mac_gfc CACHE)
  endforeach()
  if(NOT CMAKE_C_COMPILER)
    message(WARNING
      "macOS: no complete Homebrew GCC toolchain (gcc-NN + g++-NN + gfortran-NN) "
      "was found; falling back to the system default compilers. A mixed "
      "Clang/gfortran build with OpenMP is unsupported (two OpenMP runtimes) -- "
      "install a full toolchain with `brew install gcc`.")
  endif()
  unset(_oqp_mac_gcc CACHE)
  unset(_oqp_mac_gxx CACHE)
  unset(_oqp_mac_gfc CACHE)
endif()
unset(_oqp_env_cc)
unset(_oqp_env_cxx)

project(Open-Quantum-Package
  VERSION 1.3.0
  LANGUAGES C CXX Fortran)

include(GNUInstallDirs)

set(CMAKE_Fortran_STANDARD 2003)

include(FortranCInterface)

option(ENABLE_MPI "Enable MPI support" OFF)
option(ENABLE_DDX "Enable optional ddX continuum-solvation backend" OFF)
option(OQP_REUSE_EXTERNALS
       "Reuse bundled external dependency builds across fresh OpenQP build directories"
       ON)
set(OQP_EXTERNALS_ROOT "" CACHE PATH
    "Optional override for the reusable bundled external dependency cache root.")
set(OQP_SOURCE_REVISION "" CACHE STRING
    "Verified 40-character Git revision recorded in packaged build metadata.")
if(DEFINED OQP_EXTERNAL_ROOT AND NOT OQP_EXTERNALS_ROOT)
    set(OQP_EXTERNALS_ROOT "${OQP_EXTERNAL_ROOT}" CACHE PATH
        "Optional override for the reusable bundled external dependency cache root." FORCE)
    message(STATUS "OQP_EXTERNAL_ROOT is accepted as an alias for OQP_EXTERNALS_ROOT.")
endif()
mark_as_advanced(OQP_EXTERNALS_ROOT OQP_SOURCE_REVISION)

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
    Apple
#    Generic
#    ACML ACML_MP ACML_GPU
#    NAS
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

# OpenQP is ILP64-only: 8-byte BLAS/LAPACK integers AND 8-byte default
# Fortran integers on EVERY platform -- one uniform integer model across all
# hardware. macOS uses Accelerate's modern $NEWLAPACK$ILP64 interface via
# symbol interposition (requires macOS >= 13.3; see findAccelerateILP64 in
# cmake/oqp_functions.cmake); Linux uses MKL / OpenBLAS ILP64.
# LP64 (4-byte BLAS integer) support has been REMOVED: the flag below is a
# tombstone so a stale CMake cache or an old build script cannot silently
# produce a mixed-width build.
if(DEFINED LINALG_LIB_INT64 AND NOT LINALG_LIB_INT64)
    message(FATAL_ERROR
      "LINALG_LIB_INT64=OFF: LP64 (4-byte BLAS integer) support has been "
      "removed -- OpenQP is ILP64-only on all platforms (macOS >= 13.3 uses "
      "Accelerate's \$NEWLAPACK\$ILP64 interface). Remove the flag and clear "
      "any stale CMake cache.")
endif()

# Escape hatch for the no-silent-NetLib policy in findLinearAlgebra(): when ON,
# a missing mandated/requested BLAS backend falls back to the bundled NetLib
# reference BLAS with a WARNING instead of failing the configuration.
option(OQP_ALLOW_REFERENCE_BLAS "Allow fallback to the slow bundled NetLib reference BLAS when the requested backend is missing" OFF)

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
endif()
# ILP64 everywhere (unconditional): 8-byte BLAS integers and 8-byte default
# Fortran integers. See the LINALG_LIB_INT64 tombstone above.
set(BLA_SIZEOF_INTEGER 8)
set(OTR_SUFFIX "_64")
set(OTR_DEFS "-DUSE_ILP64")
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang")
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fdefault-integer-8>)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-i8>)
    # Intel Fortran represents .TRUE. as -1 (all bits set) by default.  OpenQP
    # exports logical(c_bool) fields through the C API (include/oqp.h) and the
    # Python layer reads them as C _Bool, which only admits 0 or 1.
    # -fpscomp logicals selects the 1/0 representation that the Fortran
    # standard requires for C interoperability.
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-fpscomp logicals>")
    # ifx places automatic arrays and array temporaries on the stack, gfortran
    # on the heap.  liboqp runs inside a Python process with the default 8 MB
    # main-thread stack (and smaller OpenMP worker stacks), so large automatic
    # arrays (e.g. in the MRSF-TDDFT response driver) segfault with ifx unless
    # the user raises ulimit -s.  -heap-arrays gives ifx the gfortran
    # behaviour: arrays larger than 64 KB (or of run-time size) go on the heap.
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-heap-arrays 64>")
endif()

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

# Optional GPU acceleration (openqp-gpu). Default OFF; when ON, auto-downloads
# and links the library and builds the MRSF seam in its compile-time path.
include(openqp_gpu)

if(ENABLE_DDX)
  # ddX is LGPL-3.0: it stays an *external*, dynamically-linked dependency
  # (never vendored into the OpenQP tree). Two ways to provide it:
  #   * DDX_ROOT given -> use that prebuilt install (find_package).
  #   * otherwise      -> CMake builds ddX v0.8.0 from source to match OpenQP's
  #                       selected BLAS integer width (see external/).
  # For ILP64 builds, this also avoids a stock LP64 ddX exporting an unsuffixed
  # dgemm_ that collides with OpenQP's ILP64 dgemm_ on Linux's flat namespace.
  if(DEFINED DDX_ROOT OR DEFINED ENV{DDX_ROOT})
    find_package(DDX REQUIRED)
    set(OQP_DDX_AUTOBUILD OFF)
    message(STATUS "ddX continuum-solvation backend (prebuilt): ${DDX_LIBRARY}")
  else()
    set(OQP_DDX_AUTOBUILD ON)
    message(STATUS "ddX continuum-solvation backend: building ${BLA_SIZEOF_INTEGER}-byte BLAS integer variant from source (v0.8.0)")
  endif()
endif()

# Installed wheels must never inherit external build/cache directories from
# the link search path.  Each packaged shared dependency is either installed
# beside liboqp and found through a loader-relative RPATH, or linked into
# liboqp statically.  Keep this disabled globally as well as on the oqp target
# so helper targets cannot reintroduce machine-local paths into a wheel.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

if(${LINALG_LIB} STREQUAL none)
    message(FATAL_ERROR "LINALG_LIB=none is not supported: OpenQP (and the optional OpenTrustRegion) require BLAS/LAPACK. Use LINALG_LIB=auto, a concrete BLAS/LAPACK vendor, or netlib.")
endif()

# Resolve BLAS/LAPACK before configuring external projects. OpenTrustRegion
# needs the same BLAS integer-size decision as OpenQP; early resolution also
# lets a failed system BLAS search fall back to the bundled NetLib target before
# the external project graph is generated. Linking to oqp happens when
# findLinearAlgebra() is called again after the oqp target exists.
findLinearAlgebra()

# Prebuilt ddX (find_package(DDX) above) created the DDX::ddx imported target in
# cmake/FindDDX.cmake BEFORE findLinearAlgebra() resolved BLAS/LAPACK, so it
# could not record ddX's BLAS dependency. Now that the libraries are known,
# propagate them so the DDX_ROOT path links cleanly on a flat-namespace linker
# (Linux), exactly as the autobuilt path does via _ddx_interface_libs in
# external/CMakeLists.txt. Without this a prebuilt-ddX Linux build hits the same
# "undefined reference to drot_/dtrmm_/..." failure in the smoke tests and oqp.
if(ENABLE_DDX AND NOT OQP_DDX_AUTOBUILD AND TARGET DDX::ddx)
  if(_LINALG_LIB_TYPE STREQUAL "NetLib")
    set_property(TARGET DDX::ddx PROPERTY INTERFACE_LINK_LIBRARIES ${LIBLAPACK} ${LIBBLAS})
  elseif(DEFINED LAPACK_LIBRARIES OR DEFINED BLAS_LIBRARIES)
    set_property(TARGET DDX::ddx PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  endif()
endif()

add_subdirectory(external)

add_definitions(-DOQP)
add_subdirectory(source)

if(BUILD_TESTING)
  enable_testing()
  if(BUILD_SHARED_LIBS)
    add_test(NAME oqp_dftd4_dynamic_dependencies
      COMMAND ${CMAKE_COMMAND}
        "-DOQP_LIBRARY=$<TARGET_FILE:oqp>"
        "-DDFTD4_LIBRARY=${DFTD4_DFTD4_RUNTIME_SOURCE}"
        "-DMULTICHARGE_LIBRARY=${DFTD4_MULTICHARGE_RUNTIME_SOURCE}"
        "-DMCTC_LIBRARY=${DFTD4_MCTC_RUNTIME_SOURCE}"
        "-DDFTD4_LIBRARY_DIR=${DFTD4_INSTALL}/lib"
        -P "${CMAKE_SOURCE_DIR}/cmake/check_dftd4_dynamic.cmake")
  endif()
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

# Record the configured native-library layout for the Python runtime.  Wheels
# force this to "lib" in pyproject.toml, while native distro installs may use
# GNUInstallDirs values such as lib64 or lib/x86_64-linux-gnu.
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/openqp-runtime-libdir.txt"
     "${CMAKE_INSTALL_LIBDIR}\n")
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/openqp-runtime-libdir.txt"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/openqp"
)
