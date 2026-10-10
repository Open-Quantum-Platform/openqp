// Portability shim: lets the math headers compile both with nvcc (GPU) and a
// plain host C++ compiler (CPU backend). Under a non-CUDA compiler the
// __host__/__device__ qualifiers are stripped and the CUDA runtime is omitted.
#pragma once
#ifdef __CUDACC__
  #include <cuda_runtime.h>
#else
  #ifndef __host__
  #define __host__
  #endif
  #ifndef __device__
  #define __device__
  #endif
#endif
