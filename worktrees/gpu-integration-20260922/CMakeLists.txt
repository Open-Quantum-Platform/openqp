cmake_minimum_required(VERSION 3.18)
project(openqp_gpu LANGUAGES C CXX CUDA VERSION 0.1.0)

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80)   # A100; override with -DCMAKE_CUDA_ARCHITECTURES=
endif()

option(OPENQP_GPU_DFT     "Build the DFT (XC) path into the library"     ON)
option(OPENQP_GPU_GRAD    "Build the gradient path into the library"     ON)
option(OPENQP_GPU_DFBUILD "Build the native density-fitting tensor tool" ON)

find_package(CUDAToolkit REQUIRED)
include(GNUInstallDirs)
option(OPENQP_GPU_METC "Build the independent direct-integral METC backend" ON)
option(OPENQP_GPU_METC_ONLY "Build only METC (no DF/BLAS dependency)" OFF)
if(OPENQP_GPU_METC OR OPENQP_GPU_METC_ONLY)
  add_library(openqp_gpu_metc SHARED
    src/metc/gpu_metc_cuda.cu src/metc/gpu_workspace_runtime.c)
  target_compile_definitions(openqp_gpu_metc PRIVATE OQP_CUDA_ENABLE)
  target_include_directories(openqp_gpu_metc PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
  target_link_libraries(openqp_gpu_metc PRIVATE CUDA::cudart)
  set_target_properties(openqp_gpu_metc PROPERTIES C_STANDARD 99
    CUDA_STANDARD 17 POSITION_INDEPENDENT_CODE ON)
  install(TARGETS openqp_gpu_metc LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
if(OPENQP_GPU_METC_ONLY)
  install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  return()
endif()

# ---- core library: SCF driver (+ sigma, + XC, + gradient) ----------------
# scf.cu is the density-fitting HF/hybrid SCF + J/K engine (cuBLAS/cuSOLVER,
# no rotation code). sigma.cu is the MRSF-TDDFT sigma-session (v3 ABI for the
# OpenQP OQP_ROUTEC_SIG seam). xc.cu and grad.cu are self-contained CUDA
# translation units exporting their own C entry points (routec_vxc, gradient).
set(OQPGPU_SOURCES src/scf.cu src/sigma.cu src/df_store.cu)
if(OPENQP_GPU_GRAD)
  # rotation tables embedded as a byte array (compiled in its own -O0 data TU so
  # it never burdens grad.cu's heavy device compile); makes the .so file-free.
  list(APPEND OQPGPU_SOURCES src/routec_tables_embedded.cpp)
endif()
if(OPENQP_GPU_DFT)
  list(APPEND OQPGPU_SOURCES src/xc.cu)
endif()
if(OPENQP_GPU_GRAD)
  list(APPEND OQPGPU_SOURCES src/grad.cu)
endif()

add_library(openqp_gpu SHARED ${OQPGPU_SOURCES})
add_library(openqp_gpu::openqp_gpu ALIAS openqp_gpu)

# grad.cu carries ~13 MB of generated derivative-integral kernels; -O3 (even
# with --split-compile) hangs nvcc's optimizer for 30+ min on this single TU.
# -O1 compiles in ~7 min and is NOT a runtime bottleneck (warm-call kernel
# time is ~0.1 s on H2O; the heavy Gamma stage runs in cuBLAS regardless).
if(OPENQP_GPU_GRAD)
  # grad.cu carries the generated derivative-integral device helpers. Measured
  # facts on this TU (A100, chc4): BOTH cicc -O3 (frontend NVVM) and ptxas -O3
  # (SASS optimizer) fail to finish in 20 min on the giant grad3c_al/lab
  # kernels, even after the dead __global__ routec_grad3c_k_* kernels were
  # stripped (halving the count). And it does not matter: the pass2 kernels are
  # local-memory-bound (per-triple V/Pl/DR blocks spill to local for the high-
  # angular-momentum classes), so no ptxas opt level moves the dominant cost --
  # only a block-free, Gamma-folded kernel (generator rewrite) would. So build
  # at -O1: it compiles in a few minutes and is not the runtime bottleneck.
  set_source_files_properties(src/grad.cu PROPERTIES COMPILE_OPTIONS
    "$<$<COMPILE_LANGUAGE:CUDA>:-O1;-Xcicc;-O1;-Xptxas;-O1>")
endif()
set_target_properties(openqp_gpu PROPERTIES
  OUTPUT_NAME openqp_gpu
  POSITION_INDEPENDENT_CODE ON
  CUDA_SEPARABLE_COMPILATION ON)
target_compile_options(openqp_gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3>)
target_include_directories(openqp_gpu
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>)
  target_include_directories(openqp_gpu PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(openqp_gpu PRIVATE CUDA::cublas CUDA::cusolver)

# sigma.cu (the MRSF MO-frame GEMMs) and grad.cu use a host CBLAS, so BLAS is a
# core dependency of the library -- not gradient-only.
find_package(BLAS REQUIRED)
find_path(CBLAS_INCLUDE_DIR cblas.h
  HINTS ENV EBROOTOPENBLAS ENV OPENBLAS_ROOT ENV BLAS_ROOT ENV MKLROOT
  PATH_SUFFIXES include include/openblas)
if(NOT CBLAS_INCLUDE_DIR)
  message(FATAL_ERROR "openqp-gpu needs cblas.h (used by sigma.cu and grad.cu); "
    "set OPENBLAS_ROOT / MKLROOT or load an OpenBLAS module")
endif()
target_include_directories(openqp_gpu PRIVATE ${CBLAS_INCLUDE_DIR})
target_link_libraries(openqp_gpu PRIVATE ${BLAS_LIBRARIES})

# ---- native DF-tensor builder tool --------------------------------------
# Standalone executable: reads shell/aux text + the rotation tables, writes the
# density-fitting tensor B.  This is the ONE component that uses the rotation
# integral codegen (isolated here; the library above is rotation-free). Builds
# the 3-center integrals natively on the GPU, no external QC package needed.
#   openqp_gpu_build_df <routec_tables.bin> <shells.txt> <aux.txt> <B.bin> ...
if(OPENQP_GPU_DFBUILD)
  add_executable(openqp_gpu_build_df tools/dfbuild/build_df.cu)
  target_include_directories(openqp_gpu_build_df PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tools/dfbuild)
  # -Xcompiler -fopenmp: the native (pyscf-free) aux-metric loop is host OpenMP.
  target_compile_options(openqp_gpu_build_df PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 --diag-suppress=177 -Xcompiler -fopenmp>)
  # LAPACK (dpotrf/dtrtri) for the native aux-metric Cholesky inverse.
  find_package(BLAS REQUIRED)
  target_link_libraries(openqp_gpu_build_df PRIVATE openqp_gpu CUDA::cublas ${BLAS_LIBRARIES} -fopenmp)
  set_target_properties(openqp_gpu_build_df PROPERTIES BUILD_RPATH "\$ORIGIN")

  # Same engine as a SHARED library for python ctypes: exports oqpgpu_run(argc,
  # argv, mem) taking the orbital/aux shells + packed H/S in MEMORY (mem!=null),
  # so the whole geometry->energy path is one in-process call -- no builder
  # subprocess, no shells.txt/H.pk/operand files. -Bsymbolic binds this lib's
  # own copy of the integral machinery internally (the executable carries a
  # second copy; the two must never interpose).
  add_library(openqp_gpu_df SHARED tools/dfbuild/build_df.cu)
  target_include_directories(openqp_gpu_df PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/tools/dfbuild)
  target_compile_options(openqp_gpu_df PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 --diag-suppress=177 -Xcompiler -fopenmp>)
  target_link_libraries(openqp_gpu_df PRIVATE openqp_gpu CUDA::cublas ${BLAS_LIBRARIES} -fopenmp)
  target_link_options(openqp_gpu_df PRIVATE -Wl,-Bsymbolic)
  set_target_properties(openqp_gpu_df PROPERTIES BUILD_RPATH "\$ORIGIN")
endif()

include(GNUInstallDirs)
install(TARGETS openqp_gpu EXPORT openqp_gpuTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
if(OPENQP_GPU_DFBUILD)
  install(TARGETS openqp_gpu_build_df RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES data/routec_tables.bin DESTINATION ${CMAKE_INSTALL_DATADIR}/openqp_gpu)
endif()
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT openqp_gpuTargets NAMESPACE openqp_gpu::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/openqp_gpu)

# Package config so `find_package(openqp_gpu)` works after install -- this is
# how a separate project (e.g. OpenQP with OPENQP_WITH_GPU=ON) links the library
# without vendoring or downloading this repository.
include(CMakePackageConfigHelpers)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/openqp_gpuConfig.cmake
"include(CMakeFindDependencyMacro)\n\
find_dependency(CUDAToolkit)\n\
include(\"\${CMAKE_CURRENT_LIST_DIR}/openqp_gpuTargets.cmake\")\n")
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/openqp_gpuConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/openqp_gpuConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/openqp_gpuConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/openqp_gpu)
