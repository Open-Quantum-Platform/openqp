project(pyoqp)

# Reuse the parent's Python. Fail early if parent forgot to provide it.
if(NOT TARGET Python::Interpreter OR NOT DEFINED Python_EXECUTABLE)
  message(FATAL_ERROR
    "pyoqp expects the top-level CMake to call find_package(Python ...).")
endif()

# --- Detect Python’s extension suffix, e.g. ".cpython-311-darwin.so"
execute_process(
  COMMAND "${Python_EXECUTABLE}" "-c"
          "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '')"
  OUTPUT_VARIABLE _cffi_ext_suffix
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT _cffi_ext_suffix)
  set(_cffi_ext_suffix ".so")
endif()

# --- Decide abi3 alias suffix by platform (Python modules are .so on macOS)
if(WIN32)
  set(OQP_ABI3_SUFFIX ".pyd")
elseif(APPLE)
  set(OQP_ABI3_SUFFIX ".so")
else()
  set(OQP_ABI3_SUFFIX ".so")
endif()

# --- CFFI discovery using the chosen interpreter
execute_process(
  COMMAND "${Python_EXECUTABLE}" "-c"
          "import re, cffi; print(re.compile('/__init__.py.*').sub('', cffi.__file__))"
  RESULT_VARIABLE _cffi_status
  OUTPUT_VARIABLE _cffi_location
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT _cffi_status)
  set(CFFI "${_cffi_location}" CACHE STRING "Location of Python CFFI")
endif()

execute_process(
  COMMAND "${Python_EXECUTABLE}" "-c" "import cffi; print(cffi.__version__)"
  OUTPUT_VARIABLE _cffi_version
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# --- Build the CFFI extension (OUTPUT must match the real filename)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/_oqp${_cffi_ext_suffix}"
  COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/oqp_cffi_build.py"
          -I "${CMAKE_CURRENT_SOURCE_DIR}/../include"
          -L "${CMAKE_CURRENT_BINARY_DIR}/../source"
          -loqp
  DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/oqp_cffi_build.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/../include/oqp.h"
)

# --- Copy from the built filename to the abi3 alias once
add_custom_target(liboqp
  ALL
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_BINARY_DIR}/_oqp${_cffi_ext_suffix}"
          "${CMAKE_CURRENT_BINARY_DIR}/_oqp.abi3${OQP_ABI3_SUFFIX}"
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/_oqp${_cffi_ext_suffix}"
)

# Ensure the core library is built before the Python module copy
add_dependencies(liboqp oqp)

# --- Install the abi3 alias
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/_oqp.abi3${OQP_ABI3_SUFFIX}"
  DESTINATION lib
)

INSTALL(
    DIRECTORY ${CMAKE_SOURCE_DIR}/basis_sets/
    DESTINATION share/basis_sets
    FILES_MATCHING PATTERN "*.basis"
)
