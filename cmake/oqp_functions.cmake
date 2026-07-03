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
    set(_openblas_libraries)
    if(OpenBLAS_LIBRARY)
        list(APPEND _openblas_libraries "${OpenBLAS_LIBRARY}")
    endif()
    if(TARGET OpenBLAS::OpenBLAS)
        get_target_property(_openblas_imported_location OpenBLAS::OpenBLAS IMPORTED_LOCATION)
        if(NOT _openblas_imported_location OR _openblas_imported_location STREQUAL "_openblas_imported_location-NOTFOUND")
            get_target_property(_openblas_imported_configs OpenBLAS::OpenBLAS IMPORTED_CONFIGURATIONS)
            foreach(_openblas_imported_config ${_openblas_imported_configs})
                get_target_property(_openblas_imported_location OpenBLAS::OpenBLAS "IMPORTED_LOCATION_${_openblas_imported_config}")
                if(_openblas_imported_location AND NOT _openblas_imported_location STREQUAL "_openblas_imported_location-NOTFOUND")
                    break()
                endif()
            endforeach()
        endif()
        if(_openblas_imported_location AND NOT _openblas_imported_location STREQUAL "_openblas_imported_location-NOTFOUND")
            list(APPEND _openblas_libraries "${_openblas_imported_location}")
        endif()
        get_target_property(_openblas_interface_libraries OpenBLAS::OpenBLAS INTERFACE_LINK_LIBRARIES)
        if(_openblas_interface_libraries AND NOT _openblas_interface_libraries STREQUAL "_openblas_interface_libraries-NOTFOUND")
            list(APPEND _openblas_libraries ${_openblas_interface_libraries})
        endif()
    endif()
    foreach(_openblas_library ${OpenBLAS_LIBRARIES})
        if(NOT _openblas_library STREQUAL "OpenBLAS::OpenBLAS")
            list(APPEND _openblas_libraries "${_openblas_library}")
        endif()
    endforeach()
    if(_openblas_libraries)
        list(REMOVE_DUPLICATES _openblas_libraries)
    endif()
    if(OpenBLAS_FOUND AND _openblas_libraries)
        set(_openblas_include_dir)
        if(OpenBLAS_INCLUDE_DIRS)
            list(GET OpenBLAS_INCLUDE_DIRS 0 _openblas_include_dir)
        elseif(OpenBLAS_INCLUDE_DIR)
            set(_openblas_include_dir "${OpenBLAS_INCLUDE_DIR}")
        elseif(TARGET OpenBLAS::OpenBLAS)
            get_target_property(_openblas_include_dirs OpenBLAS::OpenBLAS INTERFACE_INCLUDE_DIRECTORIES)
            if(_openblas_include_dirs AND NOT _openblas_include_dirs STREQUAL "_openblas_include_dirs-NOTFOUND")
                list(GET _openblas_include_dirs 0 _openblas_include_dir)
            endif()
        endif()
        if(_openblas_include_dir)
            set(_openblas_config_h "${_openblas_include_dir}/openblas_config.h")
        endif()
        if(_openblas_config_h AND EXISTS "${_openblas_config_h}")
            file(READ "${_openblas_config_h}" _openblas_config)
            if(LINALG_LIB_INT64 AND NOT _openblas_config MATCHES "#define[ \t]+OPENBLAS_USE64BITINT")
                message(FATAL_ERROR "OpenBLAS config at ${_openblas_config_h} is LP64, but OpenQP is configured for 8-byte BLAS integers.")
            elseif(NOT LINALG_LIB_INT64 AND _openblas_config MATCHES "#define[ \t]+OPENBLAS_USE64BITINT")
                message(FATAL_ERROR "OpenBLAS config at ${_openblas_config_h} is ILP64, but OpenQP is configured for 4-byte BLAS integers.")
            endif()
        endif()
        set(BLAS_FOUND TRUE)
        set(LAPACK_FOUND TRUE)
        set(BLAS_LIBRARIES ${_openblas_libraries})
        set(LAPACK_LIBRARIES ${_openblas_libraries})
        set(BLAS_SIZEOF_INTEGER ${BLA_SIZEOF_INTEGER})
        set(LAPACK_SIZEOF_INTEGER ${BLA_SIZEOF_INTEGER})
    endif()
    unset(_openblas_config)
    unset(_openblas_config_h)
    unset(_openblas_imported_config)
    unset(_openblas_imported_configs)
    unset(_openblas_imported_location)
    unset(_openblas_include_dir)
    unset(_openblas_include_dirs)
    unset(_openblas_interface_libraries)
    unset(_openblas_libraries)
    unset(_openblas_library)
endmacro()

# Validate the ILP64/LP64 integer width of an OpenBLAS selected through
# pkg-config, reading openblas_config.h straight from the pkg-config include dir.
#
# Regression note: CMake's FindBLAS.cmake takes an early return out of its
# pkg-config branch (Modules/FindBLAS.cmake) the moment pkg_check_modules()
# reports the requested module (e.g. openblas64) as found -- it sets BLAS_FOUND
# and BLAS_LIBRARIES and return()s WITHOUT ever setting BLAS_SIZEOF_INTEGER.
# Because of that, the BLAS_SIZEOF_INTEGER consistency gate in findBlasLapack()
# is silently skipped on the pkg-config happy path, and findOpenBLASConfig() is
# only reached on find_package() FAILURE. A hand-rolled/vendored openblas64.pc
# that actually points at an LP64 build would therefore link cleanly (the symbol
# names match) yet every BLAS call would read an 8-byte dimension/leading-dim as
# two 4-byte ints => silent numerical corruption on ILP64 builds. So we must NOT
# depend on BLAS_SIZEOF_INTEGER here: re-query pkg-config ourselves and assert
# OPENBLAS_USE64BITINT is present iff LINALG_LIB_INT64 is ON.
macro(validateOpenBLASPkgConfigInt64 _pc_module)
    if(NOT "${_pc_module}" STREQUAL "")
        find_package(PkgConfig QUIET)
        if(PKG_CONFIG_FOUND)
            pkg_check_modules(_OQP_OPENBLAS QUIET "${_pc_module}")
            if(_OQP_OPENBLAS_FOUND)
                set(_oqp_openblas_config_h)
                foreach(_oqp_inc IN LISTS _OQP_OPENBLAS_INCLUDE_DIRS)
                    if(EXISTS "${_oqp_inc}/openblas_config.h")
                        set(_oqp_openblas_config_h "${_oqp_inc}/openblas_config.h")
                        break()
                    endif()
                endforeach()
                # Only assert when we can actually read the header; a missing
                # header means we cannot prove a mismatch, and failing hard there
                # would break otherwise-valid setups that ship the .pc but no dev
                # headers. This mirrors findOpenBLASConfig()'s guarded check.
                if(_oqp_openblas_config_h)
                    file(READ "${_oqp_openblas_config_h}" _oqp_openblas_config)
                    if(LINALG_LIB_INT64 AND NOT _oqp_openblas_config MATCHES "#define[ \t]+OPENBLAS_USE64BITINT")
                        message(FATAL_ERROR "OpenBLAS pkg-config module '${_pc_module}' resolves to an LP64 build (${_oqp_openblas_config_h} does not define OPENBLAS_USE64BITINT), but OpenQP is configured for 8-byte BLAS integers (LINALG_LIB_INT64=ON). Point pkg-config at a genuine ILP64 OpenBLAS (openblas64) or set LINALG_LIB_INT64=OFF.")
                    elseif(NOT LINALG_LIB_INT64 AND _oqp_openblas_config MATCHES "#define[ \t]+OPENBLAS_USE64BITINT")
                        message(FATAL_ERROR "OpenBLAS pkg-config module '${_pc_module}' resolves to an ILP64 build (${_oqp_openblas_config_h} defines OPENBLAS_USE64BITINT), but OpenQP is configured for 4-byte BLAS integers (LINALG_LIB_INT64=OFF). Point pkg-config at an LP64 OpenBLAS (openblas) or set LINALG_LIB_INT64=ON.")
                    endif()
                    unset(_oqp_openblas_config)
                endif()
                unset(_oqp_openblas_config_h)
                unset(_oqp_inc)
            endif()
        endif()
    endif()
endmacro()

# Handle "the required BLAS/LAPACK backend was not found".
# Policy: NEVER fall back to the bundled NetLib reference BLAS silently -- that
# is how unoptimized wheels shipped and how CI ended up testing a different
# BLAS than users ran. Configuration FAILS with instructions instead. The
# reference build remains available, but only as an explicit choice:
# -DLINALG_LIB=netlib, or -DOQP_ALLOW_REFERENCE_BLAS=ON to accept the fallback.
macro(oqpLinalgNotFound _backend _hint)
    if(OQP_ALLOW_REFERENCE_BLAS)
        message(WARNING
            "BLAS/LAPACK backend '${_backend}' was not found. Falling back to "
            "the bundled NetLib reference BLAS because OQP_ALLOW_REFERENCE_BLAS=ON. "
            "This build will be SLOW (typically 1.5-2x on BLAS-bound jobs).")
        cleanBlasVars()
        set(linalg_lib netlib)
    else()
        message(FATAL_ERROR
            "Required BLAS/LAPACK backend '${_backend}' was not found.\n"
            "${_hint}\n"
            "Alternatively: pick another backend explicitly with -DLINALG_LIB=<vendor>, "
            "or accept the slow bundled reference BLAS with -DLINALG_LIB=netlib "
            "(or -DOQP_ALLOW_REFERENCE_BLAS=ON).")
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

    # =================================================================
    # Deterministic native-BLAS platform policy.
    # 'auto' resolves to the mandated native backend per OS/arch instead of
    # probing whatever BLAS happens to be installed. Environment probing made
    # the selected BLAS differ between wheels, CI, and user machines (and
    # shipped NetLib reference BLAS in the PyPI Linux wheels); the mandate
    # makes every build path pick the same backend on the same hardware:
    #     macOS (arm64 and x86_64)  -> Apple Accelerate (LP64)
    #     Linux x86_64              -> Intel MKL (ILP64)
    #     Linux aarch64             -> OpenBLAS (ILP64)
    # Platforms outside this table keep the legacy environment probe.
    # There is no silent fallback to the bundled NetLib reference BLAS --
    # see oqpLinalgNotFound() for the explicit escape hatches.
    # =================================================================
    if(linalg_lib STREQUAL auto)
      set(_oqp_linalg_mandated FALSE)
      if(APPLE)
        set(linalg_lib Apple)
        set(_oqp_linalg_mandated TRUE)
      elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64)$")
        set(linalg_lib Intel10_64ilp)
        set(_oqp_linalg_mandated TRUE)
      elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
        set(linalg_lib OpenBLAS)
        set(_oqp_linalg_mandated TRUE)
      endif()
      if(_oqp_linalg_mandated)
        message(STATUS "LINALG_LIB=auto resolved to '${linalg_lib}' for "
                       "${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR} (native-BLAS platform policy)")
      endif()
      unset(_oqp_linalg_mandated)
    endif()

    if(linalg_lib STREQUAL auto)
      # Platform without a mandated backend: legacy environment probe,
      # but a missing BLAS is an explicit failure, not a silent NetLib build.
      findBlasLapack()
      if(NOT LAPACK_FOUND OR NOT BLAS_FOUND)
          oqpLinalgNotFound(auto "CMake's FindBLAS found no BLAS/LAPACK on this platform.")
      endif()

    elseif(linalg_lib STREQUAL netlib)
        # Explicit user choice of the bundled reference BLAS: do nothing.

    else()
      if(linalg_lib STREQUAL Apple AND LINALG_LIB_INT64)
        # Accelerate's classic (Fortran-symbol) interface is LP64-only.
        message(FATAL_ERROR
            "LINALG_LIB resolves to Apple Accelerate, which only provides an "
            "LP64 (4-byte integer) BLAS/LAPACK interface, but LINALG_LIB_INT64=ON "
            "requests 8-byte BLAS integers. Configure with -DLINALG_LIB_INT64=OFF "
            "(the macOS default) or pick an ILP64 backend explicitly, e.g. "
            "-DLINALG_LIB=OpenBLAS.")
      endif()
      if(linalg_lib STREQUAL OpenBLAS)
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
      if(BLAS_FOUND AND LAPACK_FOUND AND linalg_lib STREQUAL OpenBLAS)
          # FindBLAS's pkg-config branch returns early without setting
          # BLAS_SIZEOF_INTEGER, so the width gate inside findBlasLapack() never
          # fired for this configuration. Re-validate the ILP64/LP64 width from
          # openblas_config.h so a mislabeled openblas64.pc cannot deliver a
          # silently-truncating 4-byte BLAS (see validateOpenBLASPkgConfigInt64).
          validateOpenBLASPkgConfigInt64("${BLA_PKGCONFIG_BLAS}")
      endif()
      if(NOT LAPACK_FOUND OR NOT BLAS_FOUND)
          if(linalg_lib STREQUAL OpenBLAS)
              cleanBlasVars()
              findOpenBLASConfig()
          endif()
      endif()
      if(NOT LAPACK_FOUND OR NOT BLAS_FOUND)
          # No silent re-probe or NetLib fallback: fail with instructions.
          if(linalg_lib STREQUAL Apple)
            oqpLinalgNotFound(Apple "Apple Accelerate was not found; it ships with macOS -- check the Xcode Command Line Tools / SDK installation.")
          elseif(linalg_lib MATCHES "^Intel10_64ilp")
            oqpLinalgNotFound("${linalg_lib}" "Install Intel MKL (oneAPI) and load its environment (source setvars.sh, or 'module load imkl') so MKLROOT is set before configuring.")
          elseif(linalg_lib STREQUAL OpenBLAS)
            oqpLinalgNotFound(OpenBLAS "Install an ILP64 OpenBLAS so that 'pkg-config --exists openblas64' succeeds (e.g. libopenblas64-dev / openblas built with INTERFACE64=1), or an LP64 one when LINALG_LIB_INT64=OFF.")
          else()
            oqpLinalgNotFound("${linalg_lib}" "The requested BLA_VENDOR '${linalg_lib}' was not found by CMake's FindBLAS.")
          endif()
      endif()
      # Record the integer width as a *fact of the selected backend* on paths
      # where FindBLAS does not report BLAS_SIZEOF_INTEGER, so the width gate
      # can never be skipped (Accelerate classic = LP64 by definition; the
      # Intel10_64ilp* vendor names mandate ILP64 by definition).
      if(NOT linalg_lib STREQUAL netlib AND BLAS_FOUND)
        if(NOT DEFINED BLAS_SIZEOF_INTEGER)
          if(linalg_lib STREQUAL Apple)
            set(BLAS_SIZEOF_INTEGER 4)
          elseif(linalg_lib MATCHES "^Intel10_64ilp")
            set(BLAS_SIZEOF_INTEGER 8)
          endif()
        endif()
        if(DEFINED BLAS_SIZEOF_INTEGER AND NOT BLAS_SIZEOF_INTEGER EQUAL ${BLA_SIZEOF_INTEGER})
          message(FATAL_ERROR
              "Selected BLAS ('${linalg_lib}') has ${BLAS_SIZEOF_INTEGER}-byte "
              "integers; OpenQP is configured for ${BLA_SIZEOF_INTEGER}-byte "
              "BLAS integers (LINALG_LIB_INT64=${LINALG_LIB_INT64}).")
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
