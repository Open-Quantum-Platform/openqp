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

macro(findBlasLapack)
    cleanBlasVars()
    find_package(BLAS)
    find_package(LAPACK)
endmacro()

macro(findLinearAlgebra)

    if(LINALG_LIB_INT64)
        set(BLA_SIZEOF_INTEGER 8)
    else()
        set(BLA_SIZEOF_INTEGER 4)
    endif()

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
      set(BLA_VENDOR ${linalg_lib})
      findBlasLapack()
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
      add_dependencies(oqp LAPACK)
      target_link_libraries(oqp ${LIBLAPACK} ${LIBBLAS})
    elseif(linalg_lib STREQUAL Intel10_64_dyn)
      set(_LINALG_LIB_TYPE "MKL_RT" CACHE INTERNAL "_LIANLG_LIB_TYPE")
      target_link_libraries(oqp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      if(LINALG_LIB_INT64)
         set(_MKL_INTERFACE_LAYER "ILP64" CACHE INTERNAL "_MKL_INTERFACE_LAYER")
      else()
         set(_MKL_INTERFACE_LAYER "LP64" CACHE INTERNAL "_MKL_INTERFACE_LAYER")
      endif()
    else()
      set(_LINALG_LIB_TYPE "other" CACHE INTERNAL "_LIANLG_LIB_TYPE")
      target_link_libraries(oqp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()
    unset(linalg_lib)
endmacro()
