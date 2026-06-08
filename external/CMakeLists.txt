include(ExternalProject)

if(${USE_LIBINT})
set(LIBINT2_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/libint CACHE INTERNAL "LIBINT2_SOURCE_DIR")
set(LIBINT2_PREFIX ${CMAKE_BINARY_DIR}/external/libint CACHE INTERNAL "LIBINT2_PREFIX")

ExternalProject_Add(libint2
    PREFIX ${LIBINT2_PREFIX}
    URL https://qchemlab.knu.ac.kr/vmironov/oqp-libint/-/archive/v2.7.1.1-am4/oqp-libint-v2.7.1.1-am4.tar.gz
    TLS_VERIFY TRUE
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/libint
    BINARY_DIR ${LIBINT2_PREFIX}/build-${CMAKE_Fortran_COMPILER_ID}
    STAMP_DIR  ${LIBINT2_PREFIX}/stamp-${CMAKE_Fortran_COMPILER_ID}
    BUILD_BYPRODUCTS ${LIBINT2_PREFIX}/build-${CMAKE_Fortran_COMPILER_ID}/libint2.a
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
               -DCMAKE_GENERATOR:INTERNAL=${CMAKE_GENERATOR}
)
endif()

set(NLOPT_PREFIX ${CMAKE_BINARY_DIR}/external/nlopt CACHE INTERNAL "NLOPT_PREFIX")
ExternalProject_Add(nlopt
        PREFIX ${NLOPT_PREFIX}
        URL https://github.com/stevengj/nlopt/archive/v2.9.1.tar.gz
        URL_HASH SHA256=1e6c33f8cbdc4138d525f3326c231f14ed50d99345561e85285638c49b64ee93
        BINARY_DIR ${NLOPT_PREFIX}/build-${CMAKE_Fortran_COMPILER_ID}
        STAMP_DIR  ${NLOPT_PREFIX}/stamp-${CMAKE_Fortran_COMPILER_ID}
        BUILD_BYPRODUCTS ${NLOPT_PREFIX}/build-${CMAKE_Fortran_COMPILER_ID}/libnlopt.a
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
        )

ExternalProject_Add(libxc
        PREFIX ${CMAKE_SOURCE_DIR}/external/libxc
	URL https://gitlab.com/libxc/libxc/-/archive/7.0.0/libxc-7.0.0.tar.gz
        URL_HASH SHA256=8d4e343041c9cd869833822f57744872076ae709a613c118d70605539fb13a77
        BINARY_DIR ${CMAKE_BINARY_DIR}/external/libxc/src/libxc-build-${CMAKE_Fortran_COMPILER_ID}
        STAMP_DIR  ${CMAKE_BINARY_DIR}/external/libxc/src/libxc-stamp-${CMAKE_Fortran_COMPILER_ID}
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        CMAKE_ARGS -DDISABLE_VXC=OFF
                   -DDISABLE_FXC=OFF
                   -DDISABLE_KXC=OFF
                   -DENABLE_FORTRAN=ON
                   -DBUILD_TESTING=OFF
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                   -DCMAKE_POLICY_VERSION_MINIMUM=3.25
                   -DCMAKE_GENERATOR:INTERNAL=${CMAKE_GENERATOR}
        )

ExternalProject_Add(libtagarray
        PREFIX ${CMAKE_SOURCE_DIR}/external/libtagarray
        URL https://github.com/foxtran/libtagarray/archive/refs/tags/v0.0.6.tar.gz
        URL_HASH SHA256=9c844092f123276bb461e37e3451068c097b5a1fec9d65514ad3018baaeb2b32
#        GIT_REPOSITORY https://github.com/foxtran/libtagarray.git
#        GIT_TAG master
        BINARY_DIR ${CMAKE_BINARY_DIR}/external/libtagarray/src/libtagarray-build-${CMAKE_Fortran_COMPILER_ID}
        STAMP_DIR  ${CMAKE_BINARY_DIR}/external/libtagarray/src/libtagarray-stamp-${CMAKE_Fortran_COMPILER_ID}
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        CMAKE_ARGS -DENABLE_FORTRAN=ON
                   -DENABLE_TESTING=OFF
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                   -DCMAKE_GENERATOR:INTERNAL=${CMAKE_GENERATOR}
        )

ExternalProject_Add(libecpint
        PREFIX ${CMAKE_SOURCE_DIR}/external/libecpint
        URL https://github.com/robashaw/libecpint/archive/refs/tags/v1.0.7.tar.gz
        URL_HASH SHA256=e9c60fddb2614f113ab59ec620799d961db73979845e6e637c4a6fb72aee51cc
        BINARY_DIR ${CMAKE_BINARY_DIR}/external/libecpint/src/libecpint-build-${CMAKE_Fortran_COMPILER_ID}
        STAMP_DIR  ${CMAKE_BINARY_DIR}/external/libecpint/src/libecpint-stamp-${CMAKE_Fortran_COMPILER_ID}
	INSTALL_DIR ${CMAKE_BINARY_DIR}/external/libecpint/install
        UPDATE_COMMAND ""
        CMAKE_ARGS -DLIBECPINT_BUILD_DOCS=OFF
                   -DLIBECPINT_BUILD_TESTS=OFF
		   -DLIBECPINT_USE_PUGIXML=OFF
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
		   -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        INSTALL_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target install
        )
if(_LINALG_LIB_TYPE STREQUAL NetLib)
  set(BLA_LIB_SUFFIX "64")
  if(BUILD_SHARED_LIBS)
    set(la_fortran_flags "${CMAKE_Fortran_FLAGS} -fpic")
  else()
    set(la_fortran_flags ${CMAKE_Fortran_FLAGS})
  endif()
  set(LAPACK_PREFIX ${CMAKE_BINARY_DIR}/external/lapack CACHE INTERNAL "LAPACK_PREFIX")
  set(LIBLAPACK ${LAPACK_PREFIX}/src/lapack-build/lib/liblapack${BLA_LIB_SUFFIX}.a CACHE INTERNAL "LIBLAPACK")

  set(LIBBLAS   ${LAPACK_PREFIX}/src/lapack-build/lib/libblas${BLA_LIB_SUFFIX}.a CACHE INTERNAL "LIBBLAS")
  message(STATUS "BLAS_SUFFIX = ${BLA_LIB_SUFFIX}")


  ExternalProject_Add(LAPACK
        PREFIX ${CMAKE_SOURCE_DIR}/external/lapack
        URL https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz
        URL_HASH SHA256=328c1bea493a32cac5257d84157dc686cc3ab0b004e2bea22044e0a59f6f8a19
        BINARY_DIR ${CMAKE_BINARY_DIR}/external/lapack/src/lapack-build
        STAMP_DIR  ${CMAKE_BINARY_DIR}/external/lapack/src/lapack-stamp
        BUILD_BYPRODUCTS ${LAPACK_PREFIX}/src/lapack-build/lib/liblapack${BLA_LIB_SUFFIX}.a
                         ${LAPACK_PREFIX}/src/lapack-build/lib/libblas${BLA_LIB_SUFFIX}.a
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
                   -DCMAKE_GENERATOR:INTERNAL=${CMAKE_GENERATOR}
        )
endif()
if(ENABLE_OPENTRAH)
set(otr_integer_flags "")
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(otr_integer_flags "-fdefault-integer-8")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang")
    set(otr_integer_flags "-fdefault-integer-8")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
    set(otr_integer_flags "-i8")
endif()
set(otr_fortran_flags "${CMAKE_Fortran_FLAGS} ${OTR_DEFS} ${otr_integer_flags}")
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
    set(otr_fortran_flags "${otr_fortran_flags} -fallow-argument-mismatch")
endif()
set(OTR_CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/install/opentrustregion
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBLA_SIZEOF_INTEGER=${BLA_SIZEOF_INTEGER}
        -DINTEGER_SIZE=${BLA_SIZEOF_INTEGER}
        -DBUILD_SHARED_LIBS=OFF
        -DBUILD_TESTS=OFF
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_Fortran_FLAGS=${otr_fortran_flags}
        -DCMAKE_GENERATOR:INTERNAL=${CMAKE_GENERATOR}
        -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
)
if(_LINALG_LIB_TYPE STREQUAL NetLib)
    list(APPEND OTR_CMAKE_ARGS
        -DBLAS_LIBRARIES=${LIBBLAS}
        -DLAPACK_LIBRARIES=${LIBLAPACK}
    )
elseif(DEFINED BLA_VENDOR)
    list(APPEND OTR_CMAKE_ARGS -DBLA_VENDOR=${BLA_VENDOR})
endif()

ExternalProject_Add(libopentrustregion
    PREFIX ${CMAKE_SOURCE_DIR}/external/opentrustregion
    URL https://github.com/eriksen-lab/opentrustregion/archive/refs/tags/v2.0.0.tar.gz
    URL_HASH SHA256=682201cc9cc850a892d132e5e5fc11f6d87f8cd843fdbad88dafbf7277407f1c 
    UPDATE_COMMAND ""
    # Portability: the bundled OpenTrustRegion hardcodes -march=native/-xHost in its
    # own CMakeLists (CMAKE_Fortran_FLAGS_RELEASE), so on an AVX-512 build host it emits
    # AVX-512 into level_shifted_davidson which SIGILLs (exit 132) on AVX2-only compute
    # nodes -- only OTR routines crash, since main OpenQP uses no -march. Rewrite to a
    # portable AVX2 baseline so the binary runs on any AVX2+ node.
    PATCH_COMMAND /usr/bin/perl -0pi -e "s/-march=native/-march=x86-64-v3/g; s/-xHost/-axCORE-AVX2 -mno-avx512f/g" <SOURCE_DIR>/CMakeLists.txt
    CMAKE_ARGS ${OTR_CMAKE_ARGS}
    BINARY_DIR ${CMAKE_BINARY_DIR}/external/opentrustregion/src/libopentrustregion-build-${CMAKE_Fortran_COMPILER_ID}
    #    INSTALL_DIR ${CMAKE_BINARY_DIR}/external/install/opentrustregion
    #    STEP_TARGETS build install
)
if(_LINALG_LIB_TYPE STREQUAL NetLib)
  add_dependencies(libopentrustregion LAPACK)
endif()
endif() # ENABLE_OPENTRAH

set(OTR_LIB_NAME "${CMAKE_STATIC_LIBRARY_PREFIX}opentrustregion${OTR_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")
if(ENABLE_Formatter)
  find_package(Python3 REQUIRED)
  ExternalProject_Add(fprettify
          PREFIX ${CMAKE_SOURCE_DIR}/external/fprettify
          URL https://github.com/pseewald/fprettify/archive/refs/tags/v0.3.7.tar.gz
          URL_HASH SHA256=052da19a9080a6641d3202e10572cf3d978e6bcc0e7db29c1eb8ba724e89adc7
          BINARY_DIR ""
          STAMP_DIR ""
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          INSTALL_COMMAND ""
          )
  add_custom_target(format ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/external/fprettify/src/fprettify/fprettify.py --indent 2
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
endif()



# ---------------------------------------------------------------------------
# ddX continuum-solvation backend (optional, ENABLE_DDX).
#
# ddX is LGPL-3.0 and stays an external, dynamically-linked dependency: its
# source is never vendored into OpenQP. When no prebuilt DDX_ROOT is supplied we
# build ddX v0.8.0 here, as ILP64, so its BLAS integer width matches OpenQP's
# mandatory 64-bit-integer BLAS:
#   * -fdefault-integer-8 (same as OpenTrustRegion above), and
#   * a small C-bridge patch (cmake/patches/ddx-v0.8.0-ilp64.patch) that keeps
#     ddX's int32 C ABI while letting the core be int64.
# A stock LP64 ddX would export an unsuffixed dgemm_ that collides with OpenQP's
# ILP64 dgemm_ on Linux's flat symbol namespace.
# ---------------------------------------------------------------------------
if(ENABLE_DDX AND OQP_DDX_AUTOBUILD)
    set(DDX_PREFIX  ${CMAKE_BINARY_DIR}/external/ddx CACHE INTERNAL "DDX_PREFIX")
    set(DDX_INSTALL ${CMAKE_BINARY_DIR}/external/install/ddx CACHE INTERNAL "DDX_INSTALL")
    set(DDX_INCLUDE_DIR ${DDX_INSTALL}/include CACHE INTERNAL "DDX_INCLUDE_DIR")
    set(DDX_LIBRARY
        ${DDX_INSTALL}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ddx${CMAKE_SHARED_LIBRARY_SUFFIX}
        CACHE INTERNAL "DDX_LIBRARY")

    # Match OpenTrustRegion's ILP64 Fortran flags.
    set(ddx_integer_flags "")
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang")
        set(ddx_integer_flags "-fdefault-integer-8")
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
        set(ddx_integer_flags "-i8")
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
        -DCMAKE_GENERATOR:INTERNAL=${CMAKE_GENERATOR}
    )
    # Forward CMAKE_PREFIX_PATH only when set, so a non-standard ILP64 OpenBLAS
    # location reaches the ddX subbuild. An empty -D would override the
    # environment CMAKE_PREFIX_PATH the ExternalProject otherwise inherits.
    if(CMAKE_PREFIX_PATH)
        list(APPEND DDX_CMAKE_ARGS -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH})
    endif()
    # Forward the same BLAS/LAPACK choice OpenQP resolved (findLinearAlgebra()).
    if(_LINALG_LIB_TYPE STREQUAL NetLib)
        list(APPEND DDX_CMAKE_ARGS
            -DBLAS_LIBRARIES=${LIBBLAS}
            -DLAPACK_LIBRARIES=${LIBLAPACK})
    elseif(DEFINED BLA_VENDOR)
        list(APPEND DDX_CMAKE_ARGS -DBLA_VENDOR=${BLA_VENDOR})
    endif()

    ExternalProject_Add(ddx
        PREFIX ${DDX_PREFIX}
        URL https://github.com/ddsolvation/ddX/archive/refs/tags/v0.8.0.tar.gz
        URL_HASH SHA256=6f6d08c7699dd69d4b12dd0b1f8c519f5f3f84a3d23cc526dd1a4a60914969d6
        TLS_VERIFY TRUE
        PATCH_COMMAND patch -p1 -i ${CMAKE_SOURCE_DIR}/cmake/patches/ddx-v0.8.0-ilp64.patch
        BINARY_DIR ${DDX_PREFIX}/build-${CMAKE_Fortran_COMPILER_ID}
        STAMP_DIR  ${DDX_PREFIX}/stamp-${CMAKE_Fortran_COMPILER_ID}
        BUILD_BYPRODUCTS ${DDX_LIBRARY}
        UPDATE_COMMAND ""
        CMAKE_ARGS ${DDX_CMAKE_ARGS}
    )
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
    message(STATUS "ddX will be built ILP64 from source: ${DDX_LIBRARY}")
endif()


set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)
