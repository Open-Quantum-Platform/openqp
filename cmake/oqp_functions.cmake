function(add_oqp_executable name)
  set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

  set(SOURCES_${name})
  file(GLOB SOURCES_${name} CONFIGURE_DEPENDS
    "*.F90"
    "*.f90"
    "*.c"
  )

  add_executable(${name}.x ${SOURCES_${name}})
  target_link_libraries(${name}.x oqp)
  target_include_directories(${name}.x PUBLIC ${CMAKE_BINARY_DIR}/source)
  install(TARGETS ${name}.x
    EXPORT "${PROJECT_NAME}Targets"
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endfunction()

function(add_oqp_test name)
  set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

  set(SOURCES_${name})
  file(GLOB SOURCES_${name} CONFIGURE_DEPENDS
    "${name}.F90"
  )

  add_executable(${name}.x ${SOURCES_${name}})
  target_link_libraries(${name}.x oqp)
  target_link_libraries(${name}.x oqp_testing)
  target_include_directories(${name}.x PUBLIC ${CMAKE_BINARY_DIR}/source)
  target_include_directories(${name}.x PUBLIC ${CMAKE_BINARY_DIR}/tests/oqp_testing)
endfunction()

# BLAS/LAPACK
function (cleanBlasVars)
    get_cmake_property(_vars VARIABLES)
    foreach (_var ${_vars})
      STRING(REGEX MATCH "^BLAS.*" _res ${_var})
      if(_res)
        unset (${_var} CACHE)
      endif()
    endforeach()
endfunction()

macro(findOpenBLASConfig)
    find_package(OpenBLAS CONFIG QUIET)
    if(OpenBLAS_FOUND AND OpenBLAS_LIBRARIES)
        list(GET OpenBLAS_INCLUDE_DIRS 0 _openblas_include_dir)
        set(_openblas_config_h "${_openblas_include_dir}/openblas_config.h")
        if(EXISTS "${_openblas_config_h}")
            file(READ "${_openblas_config_h}" _openblas_config)
            if(LINALG_LIB_INT64 AND NOT _openblas_config MATCHES "#define[ \t]+OPENBLAS_USE64BITINT")
                message(FATAL_ERROR "OpenBLAS config at ${_openblas_config_h} is LP64, but OpenQP is configured for 8-byte BLAS integers.")
            elseif(NOT LINALG_LIB_INT64 AND _openblas_config MATCHES "#define[ \t]+OPENBLAS_USE64BITINT")
                message(FATAL_ERROR "OpenBLAS config at ${_openblas_config_h} is ILP64, but OpenQP is configured for 4-byte BLAS integers.")
            endif()
        endif()
        set(BLAS_FOUND TRUE)
        set(LAPACK_FOUND TRUE)
        set(BLAS_LIBRARIES ${OpenBLAS_LIBRARIES})
        set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
        set(BLAS_SIZEOF_INTEGER ${BLA_SIZEOF_INTEGER})
        set(LAPACK_SIZEOF_INTEGER ${BLA_SIZEOF_INTEGER})
        unset(_openblas_config)
        unset(_openblas_config_h)
        unset(_openblas_include_dir)
    endif()
endmacro()

macro(findBlasLapack)
    cleanBlasVars()
    find_package(BLAS)
    find_package(LAPACK)
    if(BLAS_FOUND AND DEFINED BLAS_SIZEOF_INTEGER AND NOT BLAS_SIZEOF_INTEGER EQUAL ${BLA_SIZEOF_INTEGER})
        message(FATAL_ERROR "Selected BLAS reports ${BLAS_SIZEOF_INTEGER}-byte integers; OpenQP configured for ${BLA_SIZEOF_INTEGER}-byte BLAS integers.")
    endif()
    if(LAPACK_FOUND AND DEFINED LAPACK_SIZEOF_INTEGER AND NOT LAPACK_SIZEOF_INTEGER EQUAL ${BLA_SIZEOF_INTEGER})
        message(FATAL_ERROR "Selected LAPACK reports ${LAPACK_SIZEOF_INTEGER}-byte integers; OpenQP configured for ${BLA_SIZEOF_INTEGER}-byte LAPACK integers.")
    endif()
endmacro()

macro(findLinearAlgebra)

    set(linalg_lib "${LINALG_LIB}")

    if(${linalg_lib} STREQUAL auto)
      findBlasLapack()
      if(NOT ${LAPACK_FOUND} OR NOT ${BLAS_FOUND})
          # Fall back to netlib if nothing have been found
          cleanBlasVars()
          set(linalg_lib netlib)
      endif()

    elseif(${linalg_lib} STREQUAL netlib)
        # do nothing

    else()
      if(${linalg_lib} STREQUAL OpenBLAS)
        set(BLA_PREFER_PKGCONFIG ON)
        if(LINALG_LIB_INT64)
          set(BLA_PKGCONFIG_BLAS openblas64)
          set(BLA_PKGCONFIG_LAPACK openblas64)
        else()
          set(BLA_PKGCONFIG_BLAS openblas)
          set(BLA_PKGCONFIG_LAPACK openblas)
        endif()
      endif()
      set(BLA_VENDOR ${linalg_lib})
      findBlasLapack()
      if(NOT ${LAPACK_FOUND} OR NOT ${BLAS_FOUND})
          if(${linalg_lib} STREQUAL OpenBLAS)
              cleanBlasVars()
              findOpenBLASConfig()
          endif()
      endif()
      if(NOT ${LAPACK_FOUND} OR NOT ${BLAS_FOUND})
          # If nothing was found, try to find any other
          # linear algebra library in system
          unset(BLA_VENDOR)
          set(linalg_lib auto)
          findBlasLapack()
          if(NOT ${LAPACK_FOUND} OR NOT ${BLAS_FOUND})
              # Fall back to netlib if nothing have been found
              cleanBlasVars()
              set(linalg_lib netlib)
          endif()
      endif()
    endif()

    if(NOT ${linalg_lib} STREQUAL netlib)
        # Fix MKL to use single dynamic library if BLAS selected automatically
        if(${linalg_lib} STREQUAL auto
            AND DEFINED BLAS_mkl_core_LIBRARY
            AND BUILD_SHARED_LIBS
            )

            set(linalg_lib Intel10_64_dyn)
            set(BLA_VENDOR Intel10_64_dyn)
            findBlasLapack()
            if(NOT ${LAPACK_FOUND} OR NOT ${BLAS_FOUND})
                set(linalg_lib netlib)
            endif()
        endif()
    endif()

    if(linalg_lib STREQUAL netlib)
      set(_LINALG_LIB_TYPE "NetLib" CACHE INTERNAL "_LIANLG_LIB_TYPE")
      if(TARGET oqp)
        add_dependencies(oqp LAPACK)
        target_link_libraries(oqp ${LIBLAPACK} ${LIBBLAS})
      endif()
    elseif(linalg_lib STREQUAL Intel10_64_dyn)
      set(_LINALG_LIB_TYPE "MKL_RT" CACHE INTERNAL "_LIANLG_LIB_TYPE")
      if(TARGET oqp)
        target_link_libraries(oqp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      endif()
      if(LINALG_LIB_INT64)
         set(_MKL_INTERFACE_LAYER "ILP64" CACHE INTERNAL "_MKL_INTERFACE_LAYER")
      else()
         set(_MKL_INTERFACE_LAYER "LP64" CACHE INTERNAL "_MKL_INTERFACE_LAYER")
      endif()
    else()
      set(_LINALG_LIB_TYPE "other" CACHE INTERNAL "_LIANLG_LIB_TYPE")
      if(TARGET oqp)
        target_link_libraries(oqp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      endif()
    endif()
    unset(linalg_lib)
endmacro()
