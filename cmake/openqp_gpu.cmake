# openqp_gpu.cmake — drop-in CMake integration for the GPU acceleration library.
#
# Add to OpenQP's top-level CMakeLists.txt:
#
#   include(${CMAKE_SOURCE_DIR}/path/to/openqp_gpu.cmake)
#   ...
#   add_library(oqp ...)                 # the existing liboqp target
#   openqp_gpu_attach(oqp)               # <- link + enable the linked seam
#
# When OPENQP_WITH_GPU=ON the openqp-gpu library is auto-downloaded (FetchContent),
# built, and linked, and routec_sig.F90 is compiled in the linked path
# (-DOQP_GPU_LINKED) so the MRSF Davidson diverts to the GPU with no runtime
# OQP_ROUTEC_SIG dylib needed. When OFF (default) nothing is fetched and the
# seam keeps its runtime-dlopen behaviour.
#
#   cmake -B build -DOPENQP_WITH_GPU=ON [-DCMAKE_CUDA_ARCHITECTURES=80]
#
# The library is CUDA, so a CUDA toolkit (nvcc) must be present at OpenQP build
# time when the option is on. The repo is private; the fetch uses whatever git
# credentials are configured (ssh key or token). Override the source with
# -DOPENQP_GPU_REPO=<url> -DOPENQP_GPU_TAG=<ref>, or point at a local checkout /
# prebuilt install with -DOPENQP_GPU_SOURCE_DIR=<dir> to skip the download.

option(OPENQP_WITH_GPU "Link the auto-downloaded openqp-gpu CUDA library" OFF)
set(OPENQP_GPU_REPO "git@github.com:Open-Quantum-Platform/openqp-gpu.git"
    CACHE STRING "openqp-gpu git repository for auto-download")
set(OPENQP_GPU_TAG  "integrate-native-pieces"
    CACHE STRING "openqp-gpu git tag/branch to fetch")
set(OPENQP_GPU_SOURCE_DIR "" CACHE PATH
    "Use this local openqp-gpu checkout instead of downloading (optional)")

if(OPENQP_WITH_GPU)
  # The GPU library needs CUDA; enable it in this build tree.
  enable_language(CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 80)   # A100; override on the cmake line
  endif()

  # The complete GPU attachment links the SCF J/K, XC, MRSF sigma, and gradient
  # entry points, so build the full library. Only the standalone DF-tensor
  # builder tool is skipped -- it is not linked into OpenQP.
  set(OPENQP_GPU_DFBUILD OFF CACHE BOOL "" FORCE)

  if(OPENQP_GPU_SOURCE_DIR)
    message(STATUS "openqp-gpu: using local source ${OPENQP_GPU_SOURCE_DIR}")
    add_subdirectory(${OPENQP_GPU_SOURCE_DIR} openqp_gpu_build)
  else()
    include(FetchContent)
    message(STATUS "openqp-gpu: fetching ${OPENQP_GPU_REPO} @ ${OPENQP_GPU_TAG}")
    FetchContent_Declare(openqp_gpu
      GIT_REPOSITORY ${OPENQP_GPU_REPO}
      GIT_TAG        ${OPENQP_GPU_TAG}
      GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(openqp_gpu)
  endif()
endif()

# openqp_gpu_attach(<target>): wire the GPU library into an OpenQP target.
# No-op when OPENQP_WITH_GPU is OFF, so it is always safe to call.
function(openqp_gpu_attach target)
  if(OPENQP_WITH_GPU)
    target_link_libraries(${target} PRIVATE openqp_gpu::openqp_gpu)
    target_compile_definitions(${target} PRIVATE OQP_GPU_LINKED)
    message(STATUS "openqp-gpu: linked into '${target}' (routec_sig -DOQP_GPU_LINKED)")
  endif()
endfunction()
