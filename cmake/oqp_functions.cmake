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
            # OpenQP is ILP64-only: the OpenBLAS build must be 8-byte.
            if(NOT _openblas_config MATCHES "#define[ \t]+OPENBLAS_USE64BITINT")
                message(FATAL_ERROR "OpenBLAS config at ${_openblas_config_h} is LP64, but OpenQP requires 8-byte (ILP64) BLAS integers.")
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
# OPENBLAS_USE64BITINT is present (OpenQP is ILP64-only).
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
                    if(NOT _oqp_openblas_config MATCHES "#define[ \t]+OPENBLAS_USE64BITINT")
                        message(FATAL_ERROR "OpenBLAS pkg-config module '${_pc_module}' resolves to an LP64 build (${_oqp_openblas_config_h} does not define OPENBLAS_USE64BITINT), but OpenQP requires 8-byte (ILP64) BLAS integers. Point pkg-config at a genuine ILP64 OpenBLAS (openblas64).")
                    endif()
                    unset(_oqp_openblas_config)
                endif()
                unset(_oqp_openblas_config_h)
                unset(_oqp_inc)
            endif()
        endif()
    endif()
endmacro()

# =================================================================
# Accelerate ILP64 backend (macOS; OpenQP is ILP64-only).
#
# Accelerate's classic Fortran symbols (dgemm_, dgesvd_, ...) are LP64
# (4-byte integers) only. Since macOS 13.3 Accelerate ALSO ships a genuine
# ILP64 (8-byte integer) BLAS/LAPACK behind suffixed symbols
# (`dgemm$NEWLAPACK$ILP64`, ...). We interpose every classic name OpenQP
# references onto its ILP64 variant with an ld alias-list FILE, so macOS
# builds run 8-byte integers end-to-end exactly like MKL/OpenBLAS ILP64 on
# Linux -- uniform integer model on all platforms, no LP64 special case, and
# still zero bundled BLAS (Accelerate is a system framework).
#
# Two ld gotchas are baked in here (both verified the hard way):
#   * use an alias-list FILE (-Wl,-alias_list,<file>), never inline
#     -Wl,-alias,... -- the `$` in the symbol name gets shell-expanded to
#     nothing by the sh that drives the link, silently aliasing to `_dgemm`;
#   * the file format is "<real_symbol> <alias_name>" per line, i.e. the
#     ILP64 symbol FIRST and the classic name as its alias.
#
# The symbol list below was generated from the real liboqp:
#   nm -u liboqp.dylib | grep -E '^_[a-z][a-z0-9]+_$'
# with every candidate link-probed against Accelerate for a
# $NEWLAPACK$ILP64 variant (163/163 routable). If OpenQP gains a new
# BLAS/LAPACK call that is missing here, the POST_BUILD guard
# (cmake/check_accelerate_aliases.cmake) FAILS the build naming the symbol --
# a stale list cannot silently ship a 4-byte routine into an 8-byte build.
# =================================================================
set(OQP_ACCELERATE_ILP64_SYMS
    caxpy ccopy cdotc cdotu cgbmv cgemm cgemv cgerc cgeru chbmv chemm chemv
    cher cher2 cher2k cherk chpmv chpr chpr2 cscal csrot csscal cswap csymm
    csyr2k csyrk ctbmv ctbsv ctpmv ctpsv ctrmm ctrmv ctrsm ctrsv
    dasum daxpy dcopy ddot dgbmv dgemm dgemv dgeqrf dger dgesv dgesvd dgetri
    dgglse dlamch dnrm2 dorgqr dormqr drot drotm dsbmv dscal dsdot dspev
    dspevx dspmv dspr dspr2 dswap dsyev dsyevd dsymm dsymv dsyr dsyr2 dsyr2k
    dsyrk dsysv dsytrf dsytri dsytrs dtbmv dtbsv dtpmv dtpsv dtpttr dtrmm
    dtrmv dtrsm dtrsv dtrttp dzasum dznrm2
    icamax idamax isamax izamax
    sasum saxpy scasum scnrm2 scopy sdot sdsdot sgbmv sgemm sgemv sger snrm2
    srot srotm ssbmv sscal sspmv sspr sspr2 sswap ssymm ssymv ssyr ssyr2
    ssyr2k ssyrk ssytrf ssytri ssytrs stbmv stbsv stpmv stpsv strmm strmv
    strsm strsv
    xerbla
    zaxpy zcopy zdotc zdotu zdrot zdscal zgbmv zgemm zgemv zgerc zgeru zhbmv
    zheev zhemm zhemv zher zher2 zher2k zherk zhpmv zhpr zhpr2 zscal zswap
    zsymm zsyr2k zsyrk ztbmv ztbsv ztpmv ztpsv ztrmm ztrmv ztrsm ztrsv
    CACHE INTERNAL "BLAS/LAPACK routines interposed onto Accelerate ILP64")

function(oqp_accelerate_alias_file out_path)
    set(_txt "")
    foreach(_s IN LISTS OQP_ACCELERATE_ILP64_SYMS)
        # "$" is literal in CMake strings; the line format is
        # "<real symbol> <alias>": the ILP64 symbol is real, classic name aliases it.
        string(APPEND _txt "_${_s}$NEWLAPACK$ILP64 _${_s}_\n")
    endforeach()
    file(WRITE "${out_path}" "${_txt}")
endfunction()

macro(findAccelerateILP64)
    cleanBlasVars()
    # Probe: link a stub with the FULL alias list. This validates at configure
    # time that every symbol in OQP_ACCELERATE_ILP64_SYMS has a $NEWLAPACK$ILP64
    # variant in this SDK's Accelerate (needs macOS >= 13.3) -- a missing one
    # fails here with a clear story instead of at the liboqp link.
    oqp_accelerate_alias_file("${CMAKE_BINARY_DIR}/_oqp_accel_ilp64_probe.alias")
    file(WRITE "${CMAKE_BINARY_DIR}/_oqp_accel_ilp64_probe.f90" "program p\nend program\n")
    try_compile(_OQP_ACC_ILP64_OK "${CMAKE_BINARY_DIR}/_oqp_accel_ilp64_probe_dir"
        SOURCES "${CMAKE_BINARY_DIR}/_oqp_accel_ilp64_probe.f90"
        LINK_LIBRARIES "-framework Accelerate"
        LINK_OPTIONS "-Wl,-alias_list,${CMAKE_BINARY_DIR}/_oqp_accel_ilp64_probe.alias")
    if(NOT _OQP_ACC_ILP64_OK)
        message(FATAL_ERROR
            "Accelerate's ILP64 interface (\$NEWLAPACK\$ILP64 symbols) is not "
            "available from this SDK -- it requires macOS >= 13.3. Either update "
            "macOS/Xcode or pick another ILP64 backend explicitly, e.g. "
            "-DLINALG_LIB=OpenBLAS (OpenQP is ILP64-only; LP64 support was removed).")
    endif()
    oqp_accelerate_alias_file("${CMAKE_BINARY_DIR}/accelerate_ilp64.alias")
    set(OQP_ACC_ALIAS "${CMAKE_BINARY_DIR}/accelerate_ilp64.alias"
        CACHE INTERNAL "Accelerate ILP64 alias list")
    set(BLAS_FOUND   TRUE)
    set(LAPACK_FOUND TRUE)
    set(BLAS_LIBRARIES   "-framework Accelerate")
    set(LAPACK_LIBRARIES "-framework Accelerate")
    # ILP64: 8-byte integers on both sides, recorded as facts so the width
    # gate holds on this path too.
    set(BLAS_SIZEOF_INTEGER   8)
    set(LAPACK_SIZEOF_INTEGER 8)
    message(STATUS "BLAS/LAPACK: Apple Accelerate ILP64 (\$NEWLAPACK\$ILP64 interposition, "
                   "${CMAKE_BINARY_DIR}/accelerate_ilp64.alias)")
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
    # makes every build path pick the same backend on the same hardware,
    # ALL with the uniform 8-byte (ILP64) integer model by default:
    #     macOS (arm64 and x86_64)  -> Apple Accelerate ILP64
    #                                  ($NEWLAPACK$ILP64 interposition)
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

    elseif(linalg_lib STREQUAL Apple)
      # macOS: Accelerate's modern ILP64 interface via symbol interposition
      # (see findAccelerateILP64 above; OpenQP is ILP64-only).
      # FATALs internally if the SDK lacks $NEWLAPACK$ILP64 (macOS < 13.3).
      # OpenTrustRegion is a STATIC library absorbed into liboqp, so its
      # BLAS/LAPACK references resolve at the liboqp link where the alias list
      # applies -- the interposition covers it, and the POST_BUILD guard fails
      # the build naming any of its symbols missing from the alias list.
      # ddX however builds a SHARED libddx that links BLAS itself and is NOT
      # wired for the interposition yet: refuse it (OFF by default) rather
      # than let it bind the classic LP64 symbols and silently corrupt.
      if(ENABLE_DDX)
        message(FATAL_ERROR
            "ENABLE_DDX builds a shared libddx that links BLAS/LAPACK directly "
            "and is not yet wired for the Accelerate ILP64 alias interposition. "
            "Build it with an ILP64 library instead (-DLINALG_LIB=OpenBLAS) or "
            "disable ddX on macOS.")
      endif()
      findAccelerateILP64()

    else()
      if(linalg_lib STREQUAL OpenBLAS)
        # ILP64-only: always the openblas64 pkg-config module.
        set(BLA_PREFER_PKGCONFIG ON)
        set(BLA_PKGCONFIG_BLAS openblas64)
        set(BLA_PKGCONFIG_LAPACK openblas64)
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
            oqpLinalgNotFound(OpenBLAS "Install an ILP64 OpenBLAS so that 'pkg-config --exists openblas64' succeeds (e.g. libopenblas64-dev / openblas built with INTERFACE64=1),.")
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
              "integers; OpenQP requires ${BLA_SIZEOF_INTEGER}-byte (ILP64) "
              "BLAS integers.")
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
      set(_MKL_INTERFACE_LAYER "ILP64" CACHE INTERNAL "_MKL_INTERFACE_LAYER")
    elseif(linalg_lib STREQUAL Apple)
      set(_LINALG_LIB_TYPE "Accelerate_ILP64" CACHE INTERNAL "_LIANLG_LIB_TYPE")
      if(TARGET oqp)
        target_link_libraries(oqp ${BLAS_LIBRARIES})
        # Interpose every classic BLAS/LAPACK symbol onto its $NEWLAPACK$ILP64
        # variant (see findAccelerateILP64), then verify at POST_BUILD that no
        # classic LP64 symbol leaked through un-aliased -- a leak would run a
        # 4-byte-integer routine inside this 8-byte build.
        target_link_options(oqp PRIVATE "-Wl,-alias_list,${OQP_ACC_ALIAS}")
        add_custom_command(TARGET oqp POST_BUILD
            COMMAND ${CMAKE_COMMAND} -DLIB=$<TARGET_FILE:oqp>
                    -P ${CMAKE_SOURCE_DIR}/cmake/check_accelerate_aliases.cmake)
      endif()
    else()
      set(_LINALG_LIB_TYPE "other" CACHE INTERNAL "_LIANLG_LIB_TYPE")
      if(TARGET oqp)
        target_link_libraries(oqp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
      endif()
    endif()
    unset(linalg_lib)
endmacro()
