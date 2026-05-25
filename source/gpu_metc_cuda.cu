#include <cuda_runtime.h>

namespace {

__device__ __forceinline__ int idx4(int f, int m, int row, int col, int nf, int nmatrix, int nbf) {
  return f + nf * (m + nmatrix * (row + nbf * col));
}

__device__ __forceinline__ double atomic_add_double(double* address, double value) {
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
  unsigned long long int* address_as_ull = reinterpret_cast<unsigned long long int*>(address);
  unsigned long long int old = *address_as_ull;
  unsigned long long int assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(value + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
#else
  return atomicAdd(address, value);
#endif
}

__device__ __forceinline__ void add4(double* f3, int f, int m, int row, int col,
                                     int nf, int nmatrix, int nbf, double value) {
  atomic_add_double(&f3[idx4(f, m, row, col, nf, nmatrix, nbf)], value);
}

__device__ __forceinline__ double get4(const double* d3, int f, int m, int row, int col,
                                       int nf, int nmatrix, int nbf) {
  return d3[idx4(f, m, row, col, nf, nmatrix, nbf)];
}

__global__ void mrsf_metc_kernel(const int* ids, const double* ints, int ncur,
                                 double* f3, const double* d3, int nf, int nmatrix,
                                 int nbf, int cur_pass, double scale_exchange,
                                 double scale_coulomb) {
  int linear = blockIdx.x * blockDim.x + threadIdx.x;
  int total = ncur * nf * nmatrix;
  if (linear >= total) return;

  int m = linear % nmatrix;
  int tmp = linear / nmatrix;
  int f = tmp % nf;
  int n = tmp / nf;

  int i = ids[4 * n + 0] - 1;
  int j = ids[4 * n + 1] - 1;
  int k = ids[4 * n + 2] - 1;
  int l = ids[4 * n + 3] - 1;
  double val = ints[n];
  double xval = val * scale_exchange;
  double cval = val * scale_coulomb;

  if (cur_pass == 1) {
    if (m < 4) {
      add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
      add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
      add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));
      add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));
    }
    if (m < 7) {
      add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
  } else if (cur_pass == 2 && m == 6) {
    add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
    add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
    add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
    add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
  }
}

__global__ void umrsf_metc_kernel(const int* ids, const double* ints, int ncur,
                                  double* f3, const double* d3, int nf, int nmatrix,
                                  int nbf, int cur_pass, double scale_exchange,
                                  double scale_coulomb) {
  int linear = blockIdx.x * blockDim.x + threadIdx.x;
  int total = ncur * nf * nmatrix;
  if (linear >= total) return;

  int m = linear % nmatrix;
  int tmp = linear / nmatrix;
  int f = tmp % nf;
  int n = tmp / nf;

  int i = ids[4 * n + 0] - 1;
  int j = ids[4 * n + 1] - 1;
  int k = ids[4 * n + 2] - 1;
  int l = ids[4 * n + 3] - 1;
  double val = ints[n];
  double xval = val * scale_exchange;
  double cval = val * scale_coulomb;

  if (cur_pass == 1) {
    if (m < 8) {
      add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
      add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
      add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));
      add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));

      add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
    if (m == 8 || m == 9) {
      add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
      add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
      add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
      add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
      add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
      add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
      add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
      add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
    }
  }

  if ((cur_pass == 1 || cur_pass == 2) && m == 10) {
    add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
    add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
    add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
    add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
    add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
    add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
  }
}

}  // namespace

extern "C" int oqp_gpu_metc_contract(const int* ids, const double* ints, int ncur,
                                      double* f3, const double* d3, int nf,
                                      int nmatrix, int nbf, int cur_pass,
                                      double scale_exchange, double scale_coulomb,
                                      bool is_umrsf) {
  if (ncur <= 0 || nf <= 0 || nmatrix <= 0 || nbf <= 0) return 0;

  size_t ids_bytes = static_cast<size_t>(4) * ncur * sizeof(int);
  size_t ints_bytes = static_cast<size_t>(ncur) * sizeof(double);
  size_t tensor_count = static_cast<size_t>(nf) * nmatrix * nbf * nbf;
  size_t tensor_bytes = tensor_count * sizeof(double);

  int* d_ids = nullptr;
  double* d_ints = nullptr;
  double* d_f3 = nullptr;
  double* d_d3 = nullptr;

  cudaError_t err = cudaMalloc(&d_ids, ids_bytes);
  if (err != cudaSuccess) return static_cast<int>(err);
  err = cudaMalloc(&d_ints, ints_bytes);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMalloc(&d_f3, tensor_bytes);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMalloc(&d_d3, tensor_bytes);
  if (err != cudaSuccess) goto cleanup;

  err = cudaMemcpy(d_ids, ids, ids_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMemcpy(d_ints, ints, ints_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMemcpy(d_f3, f3, tensor_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) goto cleanup;
  err = cudaMemcpy(d_d3, d3, tensor_bytes, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) goto cleanup;

  {
    int threads = 256;
    int total = ncur * nf * nmatrix;
    int blocks = (total + threads - 1) / threads;
    if (is_umrsf) {
      umrsf_metc_kernel<<<blocks, threads>>>(d_ids, d_ints, ncur, d_f3, d_d3, nf,
                                             nmatrix, nbf, cur_pass, scale_exchange,
                                             scale_coulomb);
    } else {
      mrsf_metc_kernel<<<blocks, threads>>>(d_ids, d_ints, ncur, d_f3, d_d3, nf,
                                            nmatrix, nbf, cur_pass, scale_exchange,
                                            scale_coulomb);
    }
  }

  err = cudaGetLastError();
  if (err != cudaSuccess) goto cleanup;
  err = cudaDeviceSynchronize();
  if (err != cudaSuccess) goto cleanup;
  err = cudaMemcpy(f3, d_f3, tensor_bytes, cudaMemcpyDeviceToHost);

cleanup:
  cudaFree(d_ids);
  cudaFree(d_ints);
  cudaFree(d_f3);
  cudaFree(d_d3);
  return static_cast<int>(err);
}
