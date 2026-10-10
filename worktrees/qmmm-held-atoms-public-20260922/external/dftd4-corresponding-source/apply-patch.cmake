# Idempotently apply one checked-in corresponding-source patch.
#
# A reusable external cache can retain an already-patched source tree after its
# ExternalProject stamp directory is lost. In that state a plain `patch -N`
# fails. Distinguish a clean source, the exact already-applied transformation,
# and a genuinely incompatible/corrupt source without editing during probes.

foreach(_required IN ITEMS PATCH_EXECUTABLE SOURCE_DIR PATCH_FILE PATCH_STRIP)
  if(NOT DEFINED ${_required} OR "${${_required}}" STREQUAL "")
    message(FATAL_ERROR "DFT-D4 patch helper requires -D${_required}=...")
  endif()
endforeach()
if(NOT EXISTS "${PATCH_EXECUTABLE}")
  message(FATAL_ERROR "Configured patch executable does not exist: ${PATCH_EXECUTABLE}")
endif()
if(NOT IS_DIRECTORY "${SOURCE_DIR}")
  message(FATAL_ERROR "DFT-D4 patch source directory does not exist: ${SOURCE_DIR}")
endif()
if(NOT EXISTS "${PATCH_FILE}")
  message(FATAL_ERROR "DFT-D4 patch file does not exist: ${PATCH_FILE}")
endif()

set(_patch_common -F 0 -p${PATCH_STRIP} -i "${PATCH_FILE}")
execute_process(
  COMMAND "${PATCH_EXECUTABLE}" --dry-run -N ${_patch_common}
  WORKING_DIRECTORY "${SOURCE_DIR}"
  RESULT_VARIABLE _forward_result
  OUTPUT_VARIABLE _forward_stdout
  ERROR_VARIABLE _forward_stderr)
if(_forward_result EQUAL 0)
  execute_process(
    COMMAND "${PATCH_EXECUTABLE}" -N ${_patch_common}
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE _apply_result
    OUTPUT_VARIABLE _apply_stdout
    ERROR_VARIABLE _apply_stderr)
  if(NOT _apply_result EQUAL 0)
    message(FATAL_ERROR
      "DFT-D4 patch passed dry-run but failed to apply: ${PATCH_FILE}\n"
      "${_apply_stdout}\n${_apply_stderr}")
  endif()
  message(STATUS "Applied DFT-D4 corresponding-source patch: ${PATCH_FILE}")
  return()
endif()

execute_process(
  COMMAND "${PATCH_EXECUTABLE}" --dry-run -R ${_patch_common}
  WORKING_DIRECTORY "${SOURCE_DIR}"
  RESULT_VARIABLE _reverse_result
  OUTPUT_VARIABLE _reverse_stdout
  ERROR_VARIABLE _reverse_stderr)
if(_reverse_result EQUAL 0)
  message(STATUS
    "DFT-D4 corresponding-source patch is already applied: ${PATCH_FILE}")
  return()
endif()

message(FATAL_ERROR
  "DFT-D4 source matches neither side of patch: ${PATCH_FILE}\n"
  "Forward probe:\n${_forward_stdout}\n${_forward_stderr}\n"
  "Reverse probe:\n${_reverse_stdout}\n${_reverse_stderr}")
