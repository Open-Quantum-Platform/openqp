# Generate the installed DFT-D4 BUILD-INFO.json from resolved parent-CMake
# values.  Keep this helper dependency-free so it also serves as the executable
# specification for the manifest format.

function(_oqp_d4_json_quote out_var input_value)
  set(_value "${input_value}")
  string(REPLACE "\\" "\\\\" _value "${_value}")
  string(REPLACE "\"" "\\\"" _value "${_value}")
  string(REPLACE "\n" "\\n" _value "${_value}")
  string(REPLACE "\r" "\\r" _value "${_value}")
  string(REPLACE "\t" "\\t" _value "${_value}")
  set(${out_var} "\"${_value}\"" PARENT_SCOPE)
endfunction()

function(_oqp_d4_json_array out_var)
  set(_json "[")
  set(_separator "")
  foreach(_item IN LISTS ARGN)
    _oqp_d4_json_quote(_quoted "${_item}")
    string(APPEND _json "${_separator}${_quoted}")
    set(_separator ",")
  endforeach()
  string(APPEND _json "]")
  set(${out_var} "${_json}" PARENT_SCOPE)
endfunction()

function(_oqp_d4_sanitize_build_value out_var input_value)
  set(_value "${input_value}")
  foreach(_path_var IN ITEMS
      CMAKE_SOURCE_DIR CMAKE_BINARY_DIR OQP_EXTERNALS_ROOT _OQP_EXTERNALS_ROOT
      DFTD4_INSTALL
      DFTD4_MCTC_SOURCE_DIR DFTD4_MCTC_BUILD_DIR
      DFTD4_MC_SOURCE_DIR DFTD4_MC_BUILD_DIR
      DFTD4_D4_SOURCE_DIR DFTD4_D4_BUILD_DIR)
    if(DEFINED ${_path_var} AND NOT "${${_path_var}}" STREQUAL "")
      string(REPLACE "${${_path_var}}" "<${_path_var}>" _value "${_value}")
    endif()
  endforeach()
  # Frontends and CI systems can add paths outside CMAKE_BINARY_DIR (pip build
  # isolation, BuildKit cache mounts, or macOS per-user temporary folders).
  # Retain the flag itself while replacing only the volatile path token.
  string(REGEX REPLACE "(/private)?/tmp/[^ ;\t\r\n\"]+" "<transient-path>"
      _value "${_value}")
  string(REGEX REPLACE "/private/var/folders/[^ ;\t\r\n\"]+" "<transient-path>"
      _value "${_value}")
  string(REGEX REPLACE "/root/\\.cache/[^ ;\t\r\n\"]+" "<cache-path>"
      _value "${_value}")
  string(REGEX REPLACE "/home/[^/]+/\\.cache/[^ ;\t\r\n\"]+" "<cache-path>"
      _value "${_value}")
  string(REGEX REPLACE "/host-cache/[^ ;\t\r\n\"]+" "<cache-path>"
      _value "${_value}")
  string(REGEX REPLACE "/Users/[^/]+/Library/Caches/[^ ;\t\r\n\"]+"
      "<cache-path>" _value "${_value}")
  set(${out_var} "${_value}" PARENT_SCOPE)
endfunction()

function(_oqp_d4_library_labels out_var)
  set(_labels "")
  foreach(_library IN LISTS ARGN)
    if("${_library}" MATCHES "/")
      get_filename_component(_label "${_library}" NAME)
    else()
      set(_label "${_library}")
    endif()
    _oqp_d4_sanitize_build_value(_label "${_label}")
    list(APPEND _labels "${_label}")
  endforeach()
  set(${out_var} "${_labels}" PARENT_SCOPE)
endfunction()

function(oqp_generate_dftd4_build_info output_file)
  foreach(_required IN ITEMS
      PROJECT_VERSION CMAKE_GENERATOR CMAKE_SYSTEM_NAME CMAKE_SYSTEM_PROCESSOR
      CMAKE_C_COMPILER_ID CMAKE_C_COMPILER_VERSION
      CMAKE_CXX_COMPILER_ID CMAKE_CXX_COMPILER_VERSION
      CMAKE_Fortran_COMPILER_ID CMAKE_Fortran_COMPILER_VERSION
      ENABLE_OPENMP BUILD_SHARED_LIBS LINALG_LIB
      _OQP_MCTC_LIB_VERSION _OQP_MCTC_LIB_URL _OQP_MCTC_LIB_SHA256
      _OQP_MULTICHARGE_VERSION _OQP_MULTICHARGE_URL _OQP_MULTICHARGE_SHA256
      _OQP_DFTD4_VERSION _OQP_DFTD4_URL _OQP_DFTD4_SHA256
      DFTD4_MCTC_PATCH DFTD4_MCTC_PATCH_SHA256
      DFTD4_DFTD4_PATCH DFTD4_DFTD4_PATCH_SHA256
      DFTD4_MCTC_RUNTIME_NAME DFTD4_MULTICHARGE_RUNTIME_NAME
      DFTD4_DFTD4_RUNTIME_NAME BLA_SIZEOF_INTEGER _LINALG_LIB_TYPE)
    if(NOT DEFINED ${_required} OR "${${_required}}" STREQUAL "")
      message(FATAL_ERROR
        "Cannot generate DFT-D4 BUILD-INFO.json: ${_required} is unresolved")
    endif()
  endforeach()
  foreach(_sha_var IN ITEMS
      _OQP_MCTC_LIB_SHA256 _OQP_MULTICHARGE_SHA256 _OQP_DFTD4_SHA256
      DFTD4_MCTC_PATCH_SHA256 DFTD4_DFTD4_PATCH_SHA256)
    string(LENGTH "${${_sha_var}}" _sha_length)
    if(NOT "${${_sha_var}}" MATCHES "^[0-9a-fA-F]+$"
       OR NOT _sha_length EQUAL 64)
      message(FATAL_ERROR
        "Cannot generate DFT-D4 BUILD-INFO.json: ${_sha_var} is not SHA-256")
    endif()
  endforeach()

  if(_LINALG_LIB_TYPE STREQUAL "NetLib")
    set(_blas_libraries ${LIBBLAS})
    set(_lapack_libraries ${LIBLAPACK})
  else()
    set(_blas_libraries ${BLAS_LIBRARIES})
    set(_lapack_libraries ${LAPACK_LIBRARIES})
  endif()
  _oqp_d4_library_labels(_blas_labels ${_blas_libraries})
  _oqp_d4_library_labels(_lapack_labels ${_lapack_libraries})
  if(NOT _blas_labels OR NOT _lapack_labels)
    message(FATAL_ERROR
      "Cannot generate DFT-D4 BUILD-INFO.json: resolved BLAS/LAPACK libraries are empty")
  endif()
  _oqp_d4_json_array(_blas_json ${_blas_labels})
  _oqp_d4_json_array(_lapack_json ${_lapack_labels})

  foreach(_language IN ITEMS C CXX Fortran)
    get_filename_component(_compiler_name "${CMAKE_${_language}_COMPILER}" NAME)
    if(_compiler_name STREQUAL "")
      set(_compiler_name "unavailable")
    endif()
    _oqp_d4_json_quote(_compiler_${_language}_id
      "${CMAKE_${_language}_COMPILER_ID}")
    _oqp_d4_json_quote(_compiler_${_language}_version
      "${CMAKE_${_language}_COMPILER_VERSION}")
    _oqp_d4_json_quote(_compiler_${_language}_executable "${_compiler_name}")
  endforeach()

  _oqp_d4_sanitize_build_value(_c_flags "${CMAKE_C_FLAGS}")
  _oqp_d4_sanitize_build_value(_c_release_flags "${CMAKE_C_FLAGS_RELEASE}")
  _oqp_d4_sanitize_build_value(_fortran_flags "${CMAKE_Fortran_FLAGS}")
  _oqp_d4_sanitize_build_value(_fortran_release_flags
    "${CMAKE_Fortran_FLAGS_RELEASE}")
  foreach(_flag_var IN ITEMS
      _c_flags _c_release_flags _fortran_flags _fortran_release_flags)
    _oqp_d4_json_quote(${_flag_var}_json "${${_flag_var}}")
  endforeach()

  if(ENABLE_OPENMP)
    set(_openmp_json true)
  else()
    set(_openmp_json false)
  endif()

  if(BUILD_SHARED_LIBS)
    set(_build_shared_libs_json true)
  else()
    set(_build_shared_libs_json false)
  endif()

  set(_revision_json null)
  set(_dirty_json null)
  if(DEFINED OQP_SOURCE_REVISION AND NOT OQP_SOURCE_REVISION STREQUAL "")
    string(LENGTH "${OQP_SOURCE_REVISION}" _explicit_revision_length)
    if(NOT OQP_SOURCE_REVISION MATCHES "^[0-9a-fA-F]+$"
       OR NOT _explicit_revision_length EQUAL 40)
      message(FATAL_ERROR
        "OQP_SOURCE_REVISION must be a full 40-character Git SHA")
    endif()
    string(TOLOWER "${OQP_SOURCE_REVISION}" _explicit_revision)
    _oqp_d4_json_quote(_revision_json "${_explicit_revision}")
    set(_dirty_json false)
  else()
    find_package(Git QUIET)
  endif()
  if(_revision_json STREQUAL "null"
     AND Git_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" rev-parse --verify HEAD
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      RESULT_VARIABLE _revision_result
      OUTPUT_VARIABLE _revision
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(LENGTH "${_revision}" _revision_length)
    if(_revision_result EQUAL 0
       AND _revision MATCHES "^[0-9a-fA-F]+$"
       AND _revision_length EQUAL 40)
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" status --porcelain --untracked-files=normal
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        RESULT_VARIABLE _dirty_result
        OUTPUT_VARIABLE _dirty_output
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      if(_dirty_result EQUAL 0)
        if(_dirty_output STREQUAL "")
          _oqp_d4_json_quote(_revision_json "${_revision}")
          set(_dirty_json false)
        else()
          # HEAD alone does not exactly identify a dirty source tree. Keep the
          # revision null instead of presenting the base commit as the build.
          set(_dirty_json true)
        endif()
      endif()
    endif()
  endif()

  get_filename_component(_mctc_patch_name "${DFTD4_MCTC_PATCH}" NAME)
  get_filename_component(_dftd4_patch_name "${DFTD4_DFTD4_PATCH}" NAME)

  foreach(_string_var IN ITEMS
      CMAKE_VERSION CMAKE_GENERATOR CMAKE_SYSTEM_NAME CMAKE_SYSTEM_PROCESSOR
      PROJECT_VERSION LINALG_LIB _LINALG_LIB_TYPE
      _OQP_MCTC_LIB_VERSION _OQP_MCTC_LIB_URL _OQP_MCTC_LIB_SHA256
      _OQP_MULTICHARGE_VERSION _OQP_MULTICHARGE_URL _OQP_MULTICHARGE_SHA256
      _OQP_DFTD4_VERSION _OQP_DFTD4_URL _OQP_DFTD4_SHA256
      DFTD4_MCTC_PATCH_SHA256 DFTD4_DFTD4_PATCH_SHA256
      DFTD4_MCTC_RUNTIME_NAME DFTD4_MULTICHARGE_RUNTIME_NAME
      DFTD4_DFTD4_RUNTIME_NAME _mctc_patch_name _dftd4_patch_name)
    _oqp_d4_json_quote(${_string_var}_json "${${_string_var}}")
  endforeach()

  set(_json "{\n")
  string(APPEND _json
    "  \"schema\": \"org.open-quantum-platform.dftd4-build-info\",\n"
    "  \"schema_version\": 1,\n"
    "  \"openqp\": {\"version\": ${PROJECT_VERSION_json}, "
    "\"source_revision\": ${_revision_json}, "
    "\"source_tree_dirty\": ${_dirty_json}},\n"
    "  \"components\": [\n"
    "    {\"name\": \"mctc-lib\", \"version\": ${_OQP_MCTC_LIB_VERSION_json}, "
    "\"archive_url\": ${_OQP_MCTC_LIB_URL_json}, "
    "\"sha256\": ${_OQP_MCTC_LIB_SHA256_json}, \"license\": \"Apache-2.0\"},\n"
    "    {\"name\": \"multicharge\", \"version\": ${_OQP_MULTICHARGE_VERSION_json}, "
    "\"archive_url\": ${_OQP_MULTICHARGE_URL_json}, "
    "\"sha256\": ${_OQP_MULTICHARGE_SHA256_json}, \"license\": \"Apache-2.0\"},\n"
    "    {\"name\": \"dftd4\", \"version\": ${_OQP_DFTD4_VERSION_json}, "
    "\"archive_url\": ${_OQP_DFTD4_URL_json}, "
    "\"sha256\": ${_OQP_DFTD4_SHA256_json}, \"license\": \"LGPL-3.0-or-later\"}\n"
    "  ],\n"
    "  \"patches\": [\n"
    "    {\"component\": \"mctc-lib\", \"file\": ${_mctc_patch_name_json}, "
    "\"sha256\": ${DFTD4_MCTC_PATCH_SHA256_json}, \"strip\": 1},\n"
    "    {\"component\": \"dftd4\", \"file\": ${_dftd4_patch_name_json}, "
    "\"sha256\": ${DFTD4_DFTD4_PATCH_SHA256_json}, \"strip\": 1}\n"
    "  ],\n"
    "  \"build\": {\n"
    "    \"cmake_version\": ${CMAKE_VERSION_json},\n"
    "    \"generator\": ${CMAKE_GENERATOR_json},\n"
    "    \"system\": {\"name\": ${CMAKE_SYSTEM_NAME_json}, "
    "\"processor\": ${CMAKE_SYSTEM_PROCESSOR_json}},\n"
    "    \"build_type\": \"Release\",\n"
    "    \"compilers\": {\n"
    "      \"c\": {\"id\": ${_compiler_C_id}, \"version\": ${_compiler_C_version}, "
    "\"executable\": ${_compiler_C_executable}},\n"
    "      \"cxx\": {\"id\": ${_compiler_CXX_id}, \"version\": ${_compiler_CXX_version}, "
    "\"executable\": ${_compiler_CXX_executable}},\n"
    "      \"fortran\": {\"id\": ${_compiler_Fortran_id}, "
    "\"version\": ${_compiler_Fortran_version}, "
    "\"executable\": ${_compiler_Fortran_executable}}\n"
    "    },\n"
    "    \"forwarded_flags\": {\"c\": ${_c_flags_json}, "
    "\"c_release\": ${_c_release_flags_json}, "
    "\"fortran\": ${_fortran_flags_json}, "
    "\"fortran_release\": ${_fortran_release_flags_json}},\n"
    "    \"openmp\": ${_openmp_json},\n"
    "    \"blas\": {\"requested_provider\": ${LINALG_LIB_json}, "
    "\"resolved_provider\": ${_LINALG_LIB_TYPE_json}, "
    "\"resolved_blas_libraries\": ${_blas_json}, "
    "\"resolved_lapack_libraries\": ${_lapack_json}, "
    "\"integer_bytes\": ${BLA_SIZEOF_INTEGER}},\n"
    "    \"build_shared_libs\": ${_build_shared_libs_json}\n"
    "  },\n"
    "  \"canonical_runtime_names\": {"
    "\"mctc-lib\": ${DFTD4_MCTC_RUNTIME_NAME_json}, "
    "\"multicharge\": ${DFTD4_MULTICHARGE_RUNTIME_NAME_json}, "
    "\"dftd4\": ${DFTD4_DFTD4_RUNTIME_NAME_json}}\n"
    "}\n")

  foreach(_forbidden IN ITEMS
      "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${OQP_EXTERNALS_ROOT}"
      "${_OQP_EXTERNALS_ROOT}" "${DFTD4_INSTALL}")
    if(NOT "${_forbidden}" STREQUAL "")
      string(FIND "${_json}" "${_forbidden}" _path_position)
      if(NOT _path_position EQUAL -1)
        message(FATAL_ERROR
          "DFT-D4 BUILD-INFO.json contains transient path: ${_forbidden}")
      endif()
    endif()
  endforeach()

  string(JSON _schema ERROR_VARIABLE _json_error GET "${_json}" schema)
  if(_json_error)
    message(FATAL_ERROR "Generated invalid DFT-D4 BUILD-INFO.json: ${_json_error}")
  endif()
  get_filename_component(_output_dir "${output_file}" DIRECTORY)
  file(MAKE_DIRECTORY "${_output_dir}")
  file(WRITE "${output_file}" "${_json}")
endfunction()
