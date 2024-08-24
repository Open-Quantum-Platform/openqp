set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR e2k)

set(CMAKE_SYSROOT /opt/mcst/lcc-1.26.04.e2k-v4.5.4/fs)

set(tools /opt/mcst/lcc-1.26.04.e2k-v4.5.4)
set(CMAKE_C_COMPILER       ${tools}/bin/lcc)
set(CMAKE_CXX_COMPILER     ${tools}/bin/l++)
set(CMAKE_Fortran_COMPILER ${tools}/bin/lfortran)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
