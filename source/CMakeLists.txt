set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

set(SOURCES)
file(GLOB SOURCES CONFIGURE_DEPENDS "*.F90" "*.c")

# otr_interface.F90 binds against the external OpenTrustRegion module; drop it when
# OpenTRAH is disabled (run_otr then dispatches everything to the native solver).
if(NOT ENABLE_OPENTRAH)
  list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/otr_interface.F90")
endif()

add_subdirectory(io)
add_subdirectory(mathlib)
add_subdirectory(modules)
add_subdirectory(integrals)
add_subdirectory(dftlib)

# Self-test bind(C) harnesses live under tests/fortran/ rather than the
# production modules/ tree, but must be linked into liboqp so the cffi-driven
# pytest drivers (tests/test_*_selftest.py) can call them. Listed explicitly,
# not globbed: tests/fortran/ also holds standalone test PROGRAMs that must NOT
# be compiled into the library.
list(APPEND SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/../tests/fortran/fock_deriv_selftest.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/../tests/fortran/hess1_selftest.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/../tests/fortran/grd2_hess_selftest.F90"
)

message(STATUS ${SOURCES})
message(STATUS "CMAKE_INSTALL_LIBDIR = ${CMAKE_INSTALL_LIBDIR}")
set(OTR_LIB_NAME "${CMAKE_STATIC_LIBRARY_PREFIX}opentrustregion${OTR_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")
add_library(oqp ${SOURCES})
add_dependencies(oqp libxc libtagarray libecpint)
target_link_libraries(oqp xcf03 xc tagarray_f tagarray ecpint Faddeeva)
if(ENABLE_OPENTRAH)
  add_dependencies(oqp libopentrustregion)
  target_link_libraries(oqp ${OTR_LIB_NAME})
  target_compile_definitions(oqp PRIVATE OQP_HAVE_OPENTRAH)
endif()
# dlsym() in mathlib/blas_thread_ctl.c
target_link_libraries(oqp ${CMAKE_DL_LIBS})
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_link_libraries(oqp stdc++)
endif()
target_link_directories(oqp PUBLIC ${LIBXC_LIBRARY_DIR})
target_include_directories(oqp PUBLIC ${LIBXC_INCLUDE_DIR})
target_link_directories(oqp PUBLIC ${TAGARRAY_LIBRARY_DIR})
target_include_directories(oqp PUBLIC ${TAGARRAY_INCLUDE_DIRS})
target_link_directories(oqp PUBLIC ${LIBECPINT_LIBRARY_DIR}
	                              ${FADDEEVA_LIBRARY_DIR})
target_include_directories(oqp PUBLIC ${LIBECPINT_INCLUDE_DIRS})
if(ENABLE_OPENTRAH)
  target_link_directories(oqp PUBLIC ${OPENTRUSTREGION_BUILD_DIR})
  target_include_directories(oqp PUBLIC ${OPENTRUSTREGION_BUILD_DIR})
endif()

set_target_properties(oqp PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(ENABLE_OPENMP)
  target_link_libraries(oqp OpenMP::OpenMP_Fortran OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()
if(ENABLE_DDX)
  target_compile_definitions(oqp PRIVATE OQP_ENABLE_DDX)
  if(OQP_DDX_AUTOBUILD)
    add_dependencies(oqp ddx)
  endif()
  target_link_libraries(oqp DDX::ddx)
  # solvent_pcm.F90 does `use ddx_harmonics`; gfortran needs ddX's Fortran
  # module dir (ddx_harmonics.mod) on the include path. The autobuild DDX::ddx
  # target does not always propagate it, so add it explicitly (harmless if the
  # imported target already carries it).
  target_include_directories(oqp PRIVATE ${DDX_INCLUDE_DIR})
  install(FILES ${DDX_LIBRARY} DESTINATION lib)
endif()
if(ENABLE_MPI)
    target_compile_definitions(oqp PRIVATE ENABLE_MPI)
    target_link_libraries(oqp MPI::MPI_Fortran)
endif()
set(WRAPPER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/wrapper/libecpint_wrapper.cpp)
target_sources(oqp PRIVATE ${WRAPPER_SOURCES})
if(USE_LIBINT)
  target_compile_definitions(oqp PRIVATE OQP_LIBINT_ENABLE)
  add_dependencies(oqp libint2)
  target_include_directories(oqp PUBLIC ${LIBINT2_SOURCE_DIR}/include)
  target_include_directories(oqp PUBLIC ${LIBINT2_BUILD_DIR}/include)
  target_include_directories(oqp PUBLIC ${LIBINT2_BUILD_DIR}/fortran)
  target_link_libraries(oqp ${LIBINT2_LIBRARY})
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    target_link_libraries(oqp stdc++)
  endif()
endif()

# BLAS/LAPACK
cleanBlasVars()
if(NOT ${LINALG_LIB} STREQUAL none)
  findLinearAlgebra()
endif()

target_compile_definitions(oqp PRIVATE OQP_BLAS_INT=${BLA_SIZEOF_INTEGER})
set_property(SOURCE ${SOURCES_mathlib_ilp64} APPEND PROPERTY COMPILE_DEFINITIONS OQP_BLAS_INT=${BLA_SIZEOF_INTEGER})

# NLopt include
add_dependencies(oqp nlopt)
target_include_directories(oqp PUBLIC ${NLOPT_BUILD_DIR})
target_link_libraries(oqp ${NLOPT_LIBRARY})

# DFT-D4 dispersion (native, see external/CMakeLists.txt). dftd4_interface.F90
# (globbed into SOURCES above) `use`s the dftd4/mctc-lib modules, so the module
# dirs must be on the include path and the ExternalProjects must build first.
add_dependencies(oqp dftd4)
set(DFTD4_MODTAG ${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION})
target_include_directories(oqp PUBLIC
        ${DFTD4_INSTALL}/include/mctc-lib/modules
        ${DFTD4_INSTALL}/include/multicharge/${DFTD4_MODTAG}
        ${DFTD4_INSTALL}/include/dftd4/${DFTD4_MODTAG})
# dftd4/multicharge/mctc-lib reference BLAS/LAPACK symbols (e.g. ssytri_) but are
# linked after the static liblapack/libblas archives. With a static NetLib build
# the linker has already passed those archives, so the referenced objects are
# never pulled in and liboqp.so ends up with an undefined symbol that only bites
# at load time (e.g. `import oqp`). On GNU ld, wrap the dispersion libs together
# with LAPACK/BLAS in a RESCAN group (--start-group/--end-group) so the
# cross-references resolve regardless of order. This is NOT used on Apple: ld64
# does not support --start-group, and it resolves cross-archive references by a
# global symbol scan anyway, so the plain link is correct there. (For shared
# BLAS/LAPACK, e.g. macOS Accelerate, LIBLAPACK is unset and symbols resolve
# dynamically, so the plain link is used too.)
if(LIBLAPACK AND NOT APPLE)
  target_link_libraries(oqp
          "$<LINK_GROUP:RESCAN,${DFTD4_INSTALL}/lib/libdftd4.a,${DFTD4_INSTALL}/lib/libmulticharge.a,${DFTD4_INSTALL}/lib/libmctc-lib.a,${LIBLAPACK},${LIBBLAS}>")
else()
  target_link_libraries(oqp
          ${DFTD4_INSTALL}/lib/libdftd4.a
          ${DFTD4_INSTALL}/lib/libmulticharge.a
          ${DFTD4_INSTALL}/lib/libmctc-lib.a)
endif()
# The dftd4 stack is built with default 4-byte integers (see external/CMakeLists.txt),
# but this file is compiled with OpenQP's -fdefault-integer-8. Tell the shim that
# mctc-lib's `new` expects a 4-byte `num` so it converts explicitly.
set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/dftd4_interface.F90
        APPEND PROPERTY COMPILE_DEFINITIONS OQP_D4_INT_KIND=4)

install(TARGETS oqp
EXPORT "${PROJECT_NAME}Targets"
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
