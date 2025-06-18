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
        INSTALL_COMMAND ""
        CMAKE_ARGS -DBUILD_INDEX64=${LINALG_LIB_INT64}
                   -DBUILD_SINGLE=OFF
                   -DBUILD_DOUBLE=ON
                   -DBUILD_COMPLEX=OFF
                   -DBUILD_COMPLEX16=OFF
                   -DUSE_OPTIMIZED_BLAS=OFF
                   -DUSE_OPTIMIZED_LAPACK=OFF
                   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                   -DCMAKE_Fortran_FLAGS=${la_fortran_flags}
                   -DCMAKE_GENERATOR:INTERNAL=${CMAKE_GENERATOR}
        )
endif()

set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM 1)
