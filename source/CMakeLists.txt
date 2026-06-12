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

target_compile_definitions(oqp PRIVATE OQP_BLAS_INT=8)
set_property(SOURCE ${SOURCES_mathlib_ilp64} APPEND PROPERTY COMPILE_DEFINITIONS OQP_BLAS_INT=8)

# NLopt include
add_dependencies(oqp nlopt)
target_include_directories(oqp PUBLIC ${NLOPT_BUILD_DIR})
target_link_libraries(oqp ${NLOPT_LIBRARY})

install(TARGETS oqp
EXPORT "${PROJECT_NAME}Targets"
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
