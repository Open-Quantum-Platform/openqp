cmake_minimum_required(VERSION 3.25)

foreach(_required_var IN ITEMS
        OQP_LIBRARY DFTD4_LIBRARY MULTICHARGE_LIBRARY MCTC_LIBRARY
        DFTD4_LIBRARY_DIR)
  if(NOT DEFINED ${_required_var} OR "${${_required_var}}" STREQUAL "")
    message(FATAL_ERROR "${_required_var} is required")
  endif()
endforeach()

foreach(_shared_library IN ITEMS
        "${OQP_LIBRARY}" "${DFTD4_LIBRARY}"
        "${MULTICHARGE_LIBRARY}" "${MCTC_LIBRARY}")
  if(NOT EXISTS "${_shared_library}")
    message(FATAL_ERROR "Expected shared library does not exist: ${_shared_library}")
  endif()
  if(_shared_library MATCHES "\\.a$")
    message(FATAL_ERROR "Static DFT-D4 artifact is not allowed: ${_shared_library}")
  endif()
endforeach()

if(APPLE)
  find_program(_OQP_DEPENDENCY_TOOL otool REQUIRED)
  set(_OQP_DEPENDENCY_ARGS -L)
  set(_OQP_RPATH_ARGS -l)
  set(_OQP_EXPECTED_LOCAL_RPATH "@loader_path")
elseif(UNIX)
  find_program(_OQP_DEPENDENCY_TOOL readelf REQUIRED)
  set(_OQP_DEPENDENCY_ARGS -d)
  set(_OQP_RPATH_ARGS -d)
  set(_OQP_EXPECTED_LOCAL_RPATH "$ORIGIN")
else()
  message(FATAL_ERROR "Dynamic dependency inspection supports Linux and macOS")
endif()

function(_oqp_inspect output_var library)
  execute_process(
    COMMAND "${_OQP_DEPENDENCY_TOOL}" ${_OQP_DEPENDENCY_ARGS} "${library}"
    RESULT_VARIABLE _status
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _error)
  if(NOT _status EQUAL 0)
    message(FATAL_ERROR
      "Could not inspect dependencies for ${library}: ${_error}")
  endif()
  set(${output_var} "${_output}" PARENT_SCOPE)
endfunction()

function(_oqp_inspect_rpath output_var library)
  execute_process(
    COMMAND "${_OQP_DEPENDENCY_TOOL}" ${_OQP_RPATH_ARGS} "${library}"
    RESULT_VARIABLE _status
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _error)
  if(NOT _status EQUAL 0)
    message(FATAL_ERROR "Could not inspect RPATH for ${library}: ${_error}")
  endif()
  set(${output_var} "${_output}" PARENT_SCOPE)
endfunction()

function(_oqp_require_text output_var needle context)
  string(FIND "${${output_var}}" "${needle}" _position)
  if(_position EQUAL -1)
    message(FATAL_ERROR
      "${context}: expected '${needle}' in dynamic metadata:\n${${output_var}}")
  endif()
endfunction()

_oqp_inspect(_oqp_deps "${OQP_LIBRARY}")
_oqp_require_text(_oqp_deps "libdftd4" "liboqp dynamic dependency")
_oqp_require_text(_oqp_deps "libmctc-lib" "liboqp dynamic dependency")

_oqp_inspect(_dftd4_deps "${DFTD4_LIBRARY}")
_oqp_require_text(_dftd4_deps "libmulticharge" "DFT-D4 dynamic dependency")
_oqp_require_text(_dftd4_deps "libmctc-lib" "DFT-D4 dynamic dependency")

_oqp_inspect(_multicharge_deps "${MULTICHARGE_LIBRARY}")
_oqp_require_text(_multicharge_deps "libmctc-lib" "multicharge dynamic dependency")

# The build-tree liboqp intentionally searches the external install prefix.
# The package install rewrites this to loader-relative form, which is checked
# against every repaired wheel by wheel_smoke_test.py.
_oqp_inspect_rpath(_oqp_rpath "${OQP_LIBRARY}")
_oqp_require_text(_oqp_rpath "${DFTD4_LIBRARY_DIR}" "liboqp build RPATH")

foreach(_dependency_library IN ITEMS
        "${DFTD4_LIBRARY}" "${MULTICHARGE_LIBRARY}" "${MCTC_LIBRARY}")
  _oqp_inspect_rpath(_dependency_rpath "${_dependency_library}")
  _oqp_require_text(
    _dependency_rpath "${_OQP_EXPECTED_LOCAL_RPATH}"
    "DFT-D4 stack install RPATH (${_dependency_library})")
endforeach()

message(STATUS "DFT-D4 shared dependency graph and loader-relative RPATH verified")
