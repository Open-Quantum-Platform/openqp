set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

# Phaseless AFQMC. The kernels are vendored from Open-Quantum-Platform/AFQMC;
# afqmc_host.F90 supplies the host shims and afqmc_oqp.F90 is the OpenQP entry
# point. No separate executable and no file-based hand-off: everything is
# compiled straight into liboqp.
file(GLOB SOURCES_afqmc CONFIGURE_DEPENDS "*.F90")

list(APPEND SOURCES ${SOURCES_afqmc})
set(SOURCES ${SOURCES} PARENT_SCOPE)
