set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB SOURCES_modules CONFIGURE_DEPENDS "*.F90" "*.f90")

# Compile-time gating of optional native feature modules (default ON via the
# top-level ENABLE_QMMM / ENABLE_NAMD options). When OFF, drop the module source
# so its symbols are absent from liboqp; the in-core call sites are #ifdef-guarded.
if(NOT ENABLE_QMMM)
  list(REMOVE_ITEM SOURCES_modules "${CMAKE_CURRENT_SOURCE_DIR}/qmmm.F90")
endif()
if(NOT ENABLE_NAMD)
  list(REMOVE_ITEM SOURCES_modules "${CMAKE_CURRENT_SOURCE_DIR}/namd.F90")
endif()

list(APPEND SOURCES ${SOURCES_modules})
set(SOURCES ${SOURCES} PARENT_SCOPE)
