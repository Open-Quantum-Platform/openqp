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

# On Windows the linker consumes an import library while the DLL is the
# runtime artifact. Keep both paths so a prebuilt DDX_ROOT has the same target
# model as the source-built external dependency.
if(WIN32)
  unset(DDX_RUNTIME_LIBRARY CACHE)
  unset(DDX_RUNTIME_LIBRARY)
  find_file(DDX_RUNTIME_LIBRARY
    NAMES ddx.dll libddx.dll
    HINTS
      ${DDX_ROOT}
      $ENV{DDX_ROOT}
    PATH_SUFFIXES bin lib lib64 src build/src
  )
else()
  set(DDX_RUNTIME_LIBRARY "${DDX_LIBRARY}")
endif()

include(FindPackageHandleStandardArgs)
set(_DDX_REQUIRED_VARS DDX_INCLUDE_DIR DDX_LIBRARY)
if(WIN32)
  list(APPEND _DDX_REQUIRED_VARS DDX_RUNTIME_LIBRARY)
endif()
find_package_handle_standard_args(DDX
  REQUIRED_VARS ${_DDX_REQUIRED_VARS}
)

if(DDX_FOUND)
  set(DDX_INCLUDE_DIRS ${DDX_INCLUDE_DIR})
  set(DDX_LIBRARIES ${DDX_LIBRARY})

  if(NOT TARGET DDX::ddx)
    add_library(DDX::ddx UNKNOWN IMPORTED)
    set_target_properties(DDX::ddx PROPERTIES
      IMPORTED_LOCATION "${DDX_RUNTIME_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${DDX_INCLUDE_DIR}"
    )
    if(WIN32)
      set_target_properties(DDX::ddx PROPERTIES
        IMPORTED_IMPLIB "${DDX_LIBRARY}"
      )
    endif()
    # NOTE: a prebuilt shared libddx leaves its BLAS/LAPACK symbols undefined, so
    # on a flat-namespace linker (Linux) consumers must resolve them. INTERFACE_
    # LINK_LIBRARIES is attached by the top-level CMakeLists.txt right after
    # findLinearAlgebra() (this module runs before BLAS/LAPACK is resolved), so
    # the same libraries the autobuilt path uses are propagated to DDX::ddx.
  endif()
endif()

mark_as_advanced(DDX_INCLUDE_DIR DDX_LIBRARY DDX_RUNTIME_LIBRARY)
