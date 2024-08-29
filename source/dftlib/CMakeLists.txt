
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB SOURCES_mathlib CONFIGURE_DEPENDS "*.F90")

list(APPEND SOURCES ${SOURCES_mathlib})
set(SOURCES ${SOURCES} PARENT_SCOPE)

set(SOURCES_mathlib_ilp64 "${dir}/oqp_linalg.F90;${dir}/mathlib_types.F90")

set(SOURCES_mathlib_ilp64 ${SOURCES_mathlib_ilp64} PARENT_SCOPE)
