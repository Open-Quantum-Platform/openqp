set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

set(SOURCES)
file(GLOB SOURCES CONFIGURE_DEPENDS "*.F90" "*.c")

# otr_interface.F90 binds against the external OpenTrustRegion module; drop it when
# OpenTRAH is disabled (run_otr then dispatches everything to the native solver).
if(NOT ENABLE_OPENTRAH)
  list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/otr_interface.F90")
endif()

# dftd4_interface.F90 is the bind(C) shim over the external DFT-D4 library;
# drop it when the backend is disabled (see ENABLE_DFTD4 in the top-level
# CMakeLists.txt).
if(NOT ENABLE_DFTD4)
  list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/dftd4_interface.F90")
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
  "${CMAKE_CURRENT_SOURCE_DIR}/../tests/fortran/fci_symmetry_selftest.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/../tests/fortran/cholesky_eri_selftest.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/../tests/fortran/simplex_qp_selftest.F90"
)

message(STATUS ${SOURCES})
message(STATUS "CMAKE_INSTALL_LIBDIR = ${CMAKE_INSTALL_LIBDIR}")
set(OTR_LIB_NAME "${CMAKE_STATIC_LIBRARY_PREFIX}opentrustregion${OTR_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")
add_library(oqp ${SOURCES})
if(WIN32)
  # The Python runtime (pyoqp/oqp/runtime.py) loads "liboqp.<suffix>" from the
  # configured libdir on every platform.  On Windows CMake would otherwise
  # produce "oqp.dll" and install it to bin/.  Symbols are not exported from a
  # DLL by default; export them all so cffi's dlopen finds the C API.
  set_target_properties(oqp PROPERTIES
    PREFIX "lib"
    WINDOWS_EXPORT_ALL_SYMBOLS ON)
  set(_OQP_RUNTIME_DESTINATION ${CMAKE_INSTALL_LIBDIR})
else()
  set(_OQP_RUNTIME_DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
add_dependencies(oqp libxc libtagarray libecpint)
target_link_libraries(oqp xcf03 xc tagarray_f tagarray ecpint Faddeeva)
if(ENABLE_OPENTRAH)
  add_dependencies(oqp libopentrustregion)
  target_link_libraries(oqp ${OTR_LIB_NAME})
  target_compile_definitions(oqp PRIVATE OQP_HAVE_OPENTRAH)
endif()
# dlsym() in mathlib/blas_thread_ctl.c
target_link_libraries(oqp ${CMAKE_DL_LIBS})
if(WIN32)
  # routec_bridge/routec_sig reach their optional accelerator library through
  # the platform dynamic loader; on Windows that is kernel32, not libdl.
  target_compile_definitions(oqp PRIVATE OQP_DL_WIN32)
endif()
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_link_libraries(oqp stdc++)
endif()
# Optional GPU acceleration: no-op unless OPENQP_WITH_GPU=ON (see cmake/openqp_gpu.cmake).
openqp_gpu_attach(oqp)
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
  install(FILES ${DDX_LIBRARY} DESTINATION ${CMAKE_INSTALL_LIBDIR})
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

# DFT-D4 dispersion (native, see external/CMakeLists.txt). dftd4_interface.F90
# (globbed into SOURCES above) `use`s the dftd4/mctc-lib modules, so the module
# dirs must be on the include path and the ExternalProjects must build first.
if(ENABLE_DFTD4)
add_dependencies(oqp dftd4)
set(DFTD4_MODTAG ${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION})
target_include_directories(oqp PRIVATE
        ${DFTD4_INSTALL}/include/mctc-lib/modules
        ${DFTD4_INSTALL}/include/multicharge/${DFTD4_MODTAG}
        ${DFTD4_INSTALL}/include/dftd4/${DFTD4_MODTAG})
if(BUILD_SHARED_LIBS)
  # Link liboqp only to the shared libraries whose APIs the shim calls
  # directly. dftd4 records its canonical multicharge and mctc-lib SOVERSION
  # dependencies, so adding multicharge here would create a redundant edge.
  foreach(_d4_shared_library IN ITEMS
          "${DFTD4_DFTD4_LIB}"
          "${DFTD4_MULTICHARGE_LIB}"
          "${DFTD4_MCTC_LIB}")
    if(_d4_shared_library MATCHES "\\.a$")
      message(FATAL_ERROR
        "Shared OpenQP requires dynamically linked DFT-D4 libraries, found: "
        "${_d4_shared_library}")
    endif()
  endforeach()
  target_link_libraries(oqp
          "$<BUILD_INTERFACE:${DFTD4_DFTD4_LIB}>"
          "$<BUILD_INTERFACE:${DFTD4_MCTC_LIB}>")
else()
  # Static archives do not carry transitive link dependencies.  Keep all three
  # on the in-tree link line and export relocatable -L/-l requirements for
  # installed liboqp.a consumers.  Bundled NetLib is also static and was
  # linked earlier, before multicharge introduced its LAPACK references; put
  # the complete cycle in a GNU-ld rescan group so those late references are
  # resolved for both in-tree executables and downstream static consumers.
  if(LIBLAPACK AND NOT APPLE)
    target_link_libraries(oqp
            "$<BUILD_INTERFACE:$<LINK_GROUP:RESCAN,${DFTD4_DFTD4_LIB},${DFTD4_MULTICHARGE_LIB},${DFTD4_MCTC_LIB},${LIBLAPACK},${LIBBLAS}>>"
            "$<INSTALL_INTERFACE:$<LINK_GROUP:RESCAN,dftd4,multicharge,mctc-lib,lapack,blas>>")
  elseif(NOT APPLE AND (LAPACK_LIBRARIES OR BLAS_LIBRARIES))
    # findLinearAlgebra() linked the system archives before DFT-D4.  Repeat
    # them after the DFT-D4 rescan group so late multicharge references can
    # pull symbols from static OpenBLAS/MKL archives.  The repeated public
    # items also preserve this order for installed liboqp.a consumers.
    target_link_libraries(oqp
            "$<BUILD_INTERFACE:$<LINK_GROUP:RESCAN,${DFTD4_DFTD4_LIB},${DFTD4_MULTICHARGE_LIB},${DFTD4_MCTC_LIB}>>"
            "$<INSTALL_INTERFACE:$<LINK_GROUP:RESCAN,dftd4,multicharge,mctc-lib>>"
            ${LAPACK_LIBRARIES}
            ${BLAS_LIBRARIES})
  else()
    target_link_libraries(oqp
            "$<BUILD_INTERFACE:${DFTD4_DFTD4_LIB}>"
            "$<BUILD_INTERFACE:${DFTD4_MULTICHARGE_LIB}>"
            "$<BUILD_INTERFACE:${DFTD4_MCTC_LIB}>"
            "$<INSTALL_INTERFACE:dftd4>"
            "$<INSTALL_INTERFACE:multicharge>"
            "$<INSTALL_INTERFACE:mctc-lib>")
  endif()
  target_link_directories(oqp PUBLIC
          "$<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>")
endif()
endif()  # ENABLE_DFTD4
if(ENABLE_DFTD4 AND BUILD_SHARED_LIBS AND APPLE AND _LINALG_LIB_TYPE STREQUAL "Accelerate_ILP64")
  # The DFT-D4 stack is linked in separate ExternalProjects.  Fail the OpenQP
  # build if either BLAS-using dylib retained a classic LP64 Accelerate symbol
  # despite declaring WITH_ILP64=ON.
  add_custom_command(TARGET oqp POST_BUILD
    COMMAND ${CMAKE_COMMAND} "-DLIB=${DFTD4_MULTICHARGE_RUNTIME_SOURCE}"
            -P ${CMAKE_SOURCE_DIR}/cmake/check_accelerate_aliases.cmake
    COMMAND ${CMAKE_COMMAND} "-DLIB=${DFTD4_DFTD4_RUNTIME_SOURCE}"
            -P ${CMAKE_SOURCE_DIR}/cmake/check_accelerate_aliases.cmake
    COMMENT "Verifying DFT-D4 stack uses Accelerate ILP64 symbols")
endif()

# Both liboqp and every DFT-D4 dependency live in the wheel's oqp/lib
# directory.  A loader-relative RPATH makes that directory self-contained and
# lets users replace the three SOVERSION files without editing global loader
# paths.  Do not append the persistent build-cache path to installed binaries.
if(APPLE)
  set(_OQP_PACKAGE_LIB_RPATH "@loader_path")
else()
  set(_OQP_PACKAGE_LIB_RPATH "$ORIGIN")
endif()
if(BUILD_SHARED_LIBS AND ENABLE_DFTD4 AND NOT WIN32)
  set_target_properties(oqp PROPERTIES
          BUILD_RPATH "${DFTD4_INSTALL}/lib"
          INSTALL_RPATH "${_OQP_PACKAGE_LIB_RPATH}"
          INSTALL_RPATH_USE_LINK_PATH FALSE
          INSTALL_REMOVE_ENVIRONMENT_RPATH TRUE)
endif()
# The dftd4 stack is built with default 4-byte integers (see external/CMakeLists.txt),
# but this file is compiled with OpenQP's -fdefault-integer-8. Tell the shim that
# mctc-lib's `new` expects a 4-byte `num` so it converts explicitly.
if(ENABLE_DFTD4)
set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/dftd4_interface.F90
        APPEND PROPERTY COMPILE_DEFINITIONS OQP_D4_INT_KIND=4)
endif()

install(TARGETS oqp
EXPORT "${PROJECT_NAME}Targets"
  RUNTIME DESTINATION ${_OQP_RUNTIME_DESTINATION}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Wheels cannot portably preserve a shared-library symlink chain.  Install the
# fully-versioned upstream binaries as three ordinary files named for the
# runtime SONAME/install-name that liboqp and the DFT-D4 libraries reference.
# Replacing one of these files with an ABI-compatible build therefore changes
# the implementation used at the next process start.
if(NOT ENABLE_DFTD4)
  # nothing to install
elseif(BUILD_SHARED_LIBS)
  install(FILES "${DFTD4_MCTC_RUNTIME_SOURCE}"
          DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RENAME "${DFTD4_MCTC_RUNTIME_NAME}")
  install(FILES "${DFTD4_MULTICHARGE_RUNTIME_SOURCE}"
          DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RENAME "${DFTD4_MULTICHARGE_RUNTIME_NAME}")
  install(FILES "${DFTD4_DFTD4_RUNTIME_SOURCE}"
          DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RENAME "${DFTD4_DFTD4_RUNTIME_NAME}")
else()
  install(FILES
          "${DFTD4_MCTC_LIB}"
          "${DFTD4_MULTICHARGE_LIB}"
          "${DFTD4_DFTD4_LIB}"
          DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# CMake's generated macOS install_name_tool command can be defeated when the
# GNU Darwin linker has already injected @loader_path: its duplicate add makes
# the whole multi-edit command fail, leaving build-cache RPATHs in the wheel.
# Normalize the four final package files one edit at a time after all of them
# have been installed.  The configured script is intentionally macOS-only;
# Linux keeps CMake's native $ORIGIN install-RPATH handling.
if(BUILD_SHARED_LIBS AND APPLE)
  set(_OQP_MACOS_RPATH_SANITIZER
      "${CMAKE_CURRENT_BINARY_DIR}/sanitize_macos_package_rpaths.cmake")
  configure_file(
      "${CMAKE_SOURCE_DIR}/cmake/sanitize_macos_package_rpaths.cmake.in"
      "${_OQP_MACOS_RPATH_SANITIZER}"
      @ONLY)
  install(SCRIPT "${_OQP_MACOS_RPATH_SANITIZER}")
endif()

# Bundle the complete, patched source trees used to build the three libraries.
# This keeps wheel and Docker distributions self-contained for corresponding-
# source delivery instead of relying on an ephemeral BuildKit/cache directory
# or on an upstream URL remaining available. BUILD-INFO.json records resolved
# build settings; the checked-in unified patches record the two build-file-only
# modifications exactly.
if(ENABLE_DFTD4)
set(_OQP_DFTD4_SOURCE_DESTINATION
    "${CMAKE_INSTALL_DATADIR}/corresponding-source/dftd4-stack")
install(DIRECTORY "${DFTD4_MCTC_SOURCE_DIR}/"
        DESTINATION "${_OQP_DFTD4_SOURCE_DESTINATION}/${DFTD4_MCTC_SOURCE_NAME}")
install(DIRECTORY "${DFTD4_MC_SOURCE_DIR}/"
        DESTINATION "${_OQP_DFTD4_SOURCE_DESTINATION}/${DFTD4_MULTICHARGE_SOURCE_NAME}")
install(DIRECTORY "${DFTD4_D4_SOURCE_DIR}/"
        DESTINATION "${_OQP_DFTD4_SOURCE_DESTINATION}/${DFTD4_DFTD4_SOURCE_NAME}")
install(FILES
        "${CMAKE_SOURCE_DIR}/external/dftd4-corresponding-source/README.md"
        "${CMAKE_SOURCE_DIR}/external/dftd4-corresponding-source/apply-patch.cmake"
        "${CMAKE_SOURCE_DIR}/external/dftd4-corresponding-source/generate-build-info.cmake"
        "${DFTD4_BUILD_INFO_FILE}"
        DESTINATION "${_OQP_DFTD4_SOURCE_DESTINATION}")
install(FILES
        "${DFTD4_MCTC_PATCH}"
        "${DFTD4_DFTD4_PATCH}"
        DESTINATION "${_OQP_DFTD4_SOURCE_DESTINATION}/patches")
install(FILES "${CMAKE_SOURCE_DIR}/external/CMakeLists.txt"
        DESTINATION "${_OQP_DFTD4_SOURCE_DESTINATION}"
        RENAME "openqp-external-build.cmake")
endif()  # ENABLE_DFTD4
