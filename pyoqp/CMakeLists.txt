project(pyoqp)

find_package(Python3
    REQUIRED
    COMPONENTS
    Interpreter Development
)
#find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} EXACT REQUIRED)
#find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
execute_process(
    COMMAND
        ${Python3_EXECUTABLE} "-c" "import re, cffi; print(re.compile('/__init__.py.*').sub('',cffi.__file__))"
    RESULT_VARIABLE _cffi_status
    OUTPUT_VARIABLE _cffi_location
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT _cffi_status)
    set(CFFI ${_cffi_location} CACHE STRING "Location of Python CFFI")
endif()

execute_process(
    COMMAND
        ${Python3_EXECUTABLE} "-c" "import cffi; print(cffi.__version__)"
    OUTPUT_VARIABLE _cffi_version
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND
        ${Python3_EXECUTABLE} "-c" "from distutils.sysconfig import get_config_var; print(get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE _cffi_ext_suffix
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

#include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args(CFFI
#    FOUND_VAR CFFI_FOUND
#    REQUIRED_VARS CFFI
#    VERSION_VAR _cffi_version
#    )

add_custom_command(
    OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/_oqp${_cffi_ext_suffix}
    COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oqp_cffi_build.py -I ${CMAKE_CURRENT_SOURCE_DIR}/../include
                        -L ${CMAKE_CURRENT_BINARY_DIR}/../source
                        -loqp
    DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/oqp_cffi_build.py
        ${CMAKE_CURRENT_SOURCE_DIR}/../include/oqp.h
)

if (WIN32)
  set(OQP_ABI3_SUFFIX ".dll")
elseif(UNIX)
  set(OQP_ABI3_SUFFIX ".so")
endif()

add_custom_target(liboqp
    ALL
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different  ${CMAKE_CURRENT_BINARY_DIR}/_oqp${_cffi_ext_suffix}
                                               ${CMAKE_CURRENT_BINARY_DIR}/_oqp.abi3${OQP_ABI3_SUFFIX}
                                               #${CMAKE_CURRENT_BINARY_DIR}/_oqp.abi3${CMAKE_SHARED_LIBRARY_SUFFIX}
    DEPENDS
        ${CMAKE_CURRENT_BINARY_DIR}/_oqp${_cffi_ext_suffix}
)

add_dependencies(liboqp oqp)

get_filename_component(OQP_INSTALL ${CMAKE_INSTALL_PREFIX} ABSOLUTE)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/oqp.sh ${CMAKE_CURRENT_BINARY_DIR}/oqp.sh @ONLY)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/oqp_dl-find.sh ${CMAKE_CURRENT_BINARY_DIR}/oqp_dl-find.sh @ONLY)

install(
    #FILES ${CMAKE_CURRENT_BINARY_DIR}/_oqp.abi3${CMAKE_SHARED_LIBRARY_SUFFIX}
    FILES ${CMAKE_CURRENT_BINARY_DIR}/_oqp.abi3${OQP_ABI3_SUFFIX}
    #FILES ${CMAKE_CURRENT_BINARY_DIR}/_oqp${_cffi_ext_suffix}
    DESTINATION lib
)

#install(
#PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/oqp/oqp_runner.py  ${CMAKE_CURRENT_SOURCE_DIR}/oqp_nac.py
#    DESTINATION bin
#)

#install(
#    PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/oqp_runner.py ${CMAKE_CURRENT_SOURCE_DIR}/oqp_dl-find.py ${CMAKE_CURRENT_BINARY_DIR}/oqp.sh ${CMAKE_CURRENT_BINARY_DIR}/oqp_dl-find.sh
#    DESTINATION bin
#)

#install(
#    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
#    DESTINATION lib/oqp
#)

INSTALL(
    DIRECTORY ${CMAKE_SOURCE_DIR}/basis_sets/
    DESTINATION share/basis_sets
    FILES_MATCHING PATTERN "*.basis"
)
