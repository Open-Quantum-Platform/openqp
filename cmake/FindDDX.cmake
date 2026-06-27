# Find the optional ddX continuum-solvation library.
#
# This module defines:
#   DDX_FOUND
#   DDX_INCLUDE_DIRS
#   DDX_LIBRARIES
#   DDX::ddx
#
# Hints:
#   -DDX_ROOT=/path/to/ddx/install/or/build
#   environment variable DDX_ROOT

find_path(DDX_INCLUDE_DIR
  NAMES ddx.h
  HINTS
    ${DDX_ROOT}
    $ENV{DDX_ROOT}
  PATH_SUFFIXES include src src/src
)

find_library(DDX_LIBRARY
  NAMES ddx libddx
  HINTS
    ${DDX_ROOT}
    $ENV{DDX_ROOT}
  PATH_SUFFIXES lib lib64 src build/src
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DDX
  REQUIRED_VARS DDX_INCLUDE_DIR DDX_LIBRARY
)

if(DDX_FOUND)
  set(DDX_INCLUDE_DIRS ${DDX_INCLUDE_DIR})
  set(DDX_LIBRARIES ${DDX_LIBRARY})

  if(NOT TARGET DDX::ddx)
    add_library(DDX::ddx UNKNOWN IMPORTED)
    set_target_properties(DDX::ddx PROPERTIES
      IMPORTED_LOCATION "${DDX_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${DDX_INCLUDE_DIR}"
    )
    # NOTE: a prebuilt shared libddx leaves its BLAS/LAPACK symbols undefined, so
    # on a flat-namespace linker (Linux) consumers must resolve them. INTERFACE_
    # LINK_LIBRARIES is attached by the top-level CMakeLists.txt right after
    # findLinearAlgebra() (this module runs before BLAS/LAPACK is resolved), so
    # the same libraries the autobuilt path uses are propagated to DDX::ddx.
  endif()
endif()

mark_as_advanced(DDX_INCLUDE_DIR DDX_LIBRARY)
