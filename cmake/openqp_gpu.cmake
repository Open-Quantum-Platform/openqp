# openqp_gpu.cmake — optional link to the SEPARATE openqp-gpu CUDA library.
#
# openqp-gpu is a standalone project with its own repository and build. OpenQP
# does NOT bundle, vendor, or download it -- so the public OpenQP tree has no
# dependency on (and no reference to) that repository. When OPENQP_WITH_GPU=ON,
# OpenQP links a library you have built/obtained separately:
#
#   * an installed package on CMAKE_PREFIX_PATH   (find_package(openqp_gpu)), or
#   * a local source checkout                     (-DOPENQP_GPU_SOURCE_DIR=<dir>).
#
# Default is OFF: stock OpenQP, unchanged. Even then the seams remain usable at
# runtime via dlopen (set the OQP_ROUTEC_* env vars to a prebuilt libopenqp_gpu),
# so the GPU backend is a fully separable, pluggable component.
#
# Add to OpenQP's top-level CMakeLists.txt:
#   include(${CMAKE_SOURCE_DIR}/path/to/openqp_gpu.cmake)
#   ...
#   add_library(oqp ...)
#   openqp_gpu_attach(oqp)
#
#   cmake -B build -DOPENQP_WITH_GPU=ON \
#         [-DCMAKE_PREFIX_PATH=<openqp-gpu install prefix>] \
#         [-DOPENQP_GPU_SOURCE_DIR=<openqp-gpu checkout>] \
#         [-DCMAKE_CUDA_ARCHITECTURES=80]

option(OPENQP_WITH_GPU "Link the separately built openqp-gpu CUDA library" OFF)
set(OPENQP_GPU_SOURCE_DIR "" CACHE PATH
    "Path to a local openqp-gpu source checkout to build in-tree (optional; otherwise find_package(openqp_gpu) locates an installed package)")

if(OPENQP_WITH_GPU)
  # Select the target GPU architecture BEFORE enabling CUDA, so CMake
  # initializes CMAKE_CUDA_ARCHITECTURES (and the GPU targets) from it.
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 80)   # A100; override on the cmake line
  endif()
  enable_language(CUDA)
  # Only the linkable library is needed; the standalone DF-tensor builder tool
  # is not built when openqp-gpu is added as a subproject.
  set(OPENQP_GPU_DFBUILD OFF CACHE BOOL "" FORCE)

  if(OPENQP_GPU_SOURCE_DIR)
    message(STATUS "openqp-gpu: building from local source ${OPENQP_GPU_SOURCE_DIR}")
    add_subdirectory(${OPENQP_GPU_SOURCE_DIR} openqp_gpu_build)
  else()
    find_package(openqp_gpu QUIET)
    if(NOT openqp_gpu_FOUND)
      message(FATAL_ERROR
        "OPENQP_WITH_GPU=ON but the openqp-gpu library was not found.\n"
        "openqp-gpu is a separate project (built and installed on its own); "
        "OpenQP does not download it. Build/install openqp-gpu, then either\n"
        "  - add its install prefix to CMAKE_PREFIX_PATH, or\n"
        "  - point -DOPENQP_GPU_SOURCE_DIR=<openqp-gpu checkout>.")
    endif()
    message(STATUS "openqp-gpu: using installed package (${openqp_gpu_DIR})")
  endif()
endif()

# openqp_gpu_attach(<target>): wire the GPU library into an OpenQP target.
# No-op when OPENQP_WITH_GPU is OFF, so it is always safe to call.
function(openqp_gpu_attach target)
  if(OPENQP_WITH_GPU)
    # OpenQP links `oqp` with the plain target_link_libraries() signature
    # (source/CMakeLists.txt); CMake forbids mixing plain and keyword
    # signatures on one target, so use the plain form here too.
    target_link_libraries(${target} openqp_gpu::openqp_gpu)
    target_compile_definitions(${target} PRIVATE OQP_GPU_LINKED)
    message(STATUS "openqp-gpu: linked into '${target}' (-DOQP_GPU_LINKED)")
  endif()
endfunction()
