cmake_minimum_required(VERSION 3.22)
project(opentrustregion Fortran C)
include(GNUInstallDirs)

# default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# set flags for build types
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wextra -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "-fpscomp logicals")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn all -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -xHost")
endif()

# try finding 32-bit BLAS/LAPACK first
set(BLA_SIZEOF_INTEGER 4 CACHE STRING "Size of BLAS/LAPACK integer type (4 for LP64, 8 for ILP64)")
find_package(BLAS)
find_package(LAPACK)

# check if 32-bit BLAS/LAPACK found
if (NOT BLAS_FOUND OR NOT LAPACK_FOUND)
    message(STATUS "BLAS/LAPACK with 32-bit integers not found, trying version with 64-bit integers...")

    # now try finding 64-bit BLAS/LAPACK
    set(BLA_SIZEOF_INTEGER 8 CACHE STRING "Size of BLAS/LAPACK integer type (4 for LP64, 8 for ILP64)" FORCE)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    if (BLAS_FOUND AND LAPACK_FOUND)
        message(STATUS "Using BLAS/LAPACK with 64-bit integers (ILP64) ")
        add_definitions(-DUSE_ILP64)
    else()
        message(FATAL_ERROR "Could not find a suitable BLAS/LAPACK library (neither with 32-bit or 64-bit integers).")
    endif()
else()
    message(STATUS "Using BLAS/LAPACK with 32-bit integers (LP64)")
endif()

# define source files
set(OPENTRUSTREGION_SOURCES src/opentrustregion.f90 src/c_interface.f90)

add_library(opentrustregion ${OPENTRUSTREGION_SOURCES})

# common settings (moved from foreach)
set_target_properties(opentrustregion PROPERTIES
    OUTPUT_NAME opentrustregion
    POSITION_INDEPENDENT_CODE ON
)

# link against BLAS and LAPACK
target_link_libraries(opentrustregion PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# set include directories
target_include_directories(opentrustregion PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# enable preprocessing
set_source_files_properties(${OPENTRUSTREGION_SOURCES} PROPERTIES Fortran_PREPROCESS ON)

# option to build testsuite
option(BUILD_TESTS "Build the testsuite" ON)

if (BUILD_TESTS)
    # create shared test library which is then called from Python for testing
    add_library(
        testsuite SHARED
        tests/opentrustregion_unit_tests.f90
        tests/c_interface_unit_tests.f90
        tests/opentrustregion_mock.f90
        tests/c_interface_mock.f90
        tests/opentrustregion_system_tests.f90
    )
  
    # link against BLAS and LAPACK and against OpenTrustRegion
    target_link_libraries(
        testsuite
        PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} opentrustregion
    )
  
    # get the absolute path to the test data directory and add it as a definition so the
    # Fortran code can access it
    set(TEST_DATA_PATH "${CMAKE_SOURCE_DIR}/tests/data")
    add_definitions(-DTEST_DATA_PATH=\"${TEST_DATA_PATH}\")
    target_compile_definitions(testsuite PRIVATE TEST_DATA_PATH="${TEST_DATA_PATH}")
    
    # enable preprocessing
    set_source_files_properties(tests/opentrustregion_system_tests.f90 PROPERTIES Fortran_PREPROCESS ON)
endif()

# set paths to find module files, header files and library in build directory
set(OpenTrustRegion_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/ ${CMAKE_BINARY_DIR})

set(OpenTrustRegion_LIBRARY $<TARGET_FILE:opentrustregion>)

# create the OpenTrustRegionConfig.cmake from the template
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/OpenTrustRegionConfig.cmake.in
  ${CMAKE_BINARY_DIR}/cmake/OpenTrustRegionConfig.cmake
  @ONLY
)

# export the target to make OpenTrustRegion::OpenTrustRegion usable
install(TARGETS opentrustregion
        EXPORT OpenTrustRegionTargets
        LIBRARY DESTINATION lib # shared libraries (.so, .dylib)
        ARCHIVE DESTINATION lib # static libraries (.a)
        INCLUDES DESTINATION include)
install(FILES include/opentrustregion.h DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/)

# export the target with namespace OpenTrustRegion::
install(EXPORT OpenTrustRegionTargets
        NAMESPACE OpenTrustRegion::
        DESTINATION lib/cmake/OpenTrustRegion)

# generate and install the config file
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/OpenTrustRegionConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/OpenTrustRegionConfig.cmake
    INSTALL_DESTINATION lib/cmake/OpenTrustRegion
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/OpenTrustRegionConfig.cmake
        DESTINATION lib/cmake/OpenTrustRegion)
