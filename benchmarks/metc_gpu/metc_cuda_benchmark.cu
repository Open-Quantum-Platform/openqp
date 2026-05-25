#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "../../source/gpu_metc_cuda.cu"

namespace {

inline std::size_t cpu_idx4(int f, int m, int row, int col, int nf, int nmatrix, int nbf) {
  return static_cast<std::size_t>(f) + static_cast<std::size_t>(nf) *
      (static_cast<std::size_t>(m) + static_cast<std::size_t>(nmatrix) *
      (static_cast<std::size_t>(row) + static_cast<std::size_t>(nbf) * static_cast<std::size_t>(col)));
}

inline void add4(std::vector<double>& f3, int f, int m, int row, int col,
                 int nf, int nmatrix, int nbf, double value) {
  f3[cpu_idx4(f, m, row, col, nf, nmatrix, nbf)] += value;
}

inline double get4(const std::vector<double>& d3, int f, int m, int row, int col,
                   int nf, int nmatrix, int nbf) {
  return d3[cpu_idx4(f, m, row, col, nf, nmatrix, nbf)];
}

void cpu_mrsf_contract(const std::vector<int>& ids, const std::vector<double>& ints, int ncur,
                       std::vector<double>& f3, const std::vector<double>& d3, int nf,
                       int nmatrix, int nbf, int cur_pass, double scale_exchange,
                       double scale_coulomb) {
  for (int n = 0; n < ncur; ++n) {
    int i = ids[4 * n + 0] - 1;
    int j = ids[4 * n + 1] - 1;
    int k = ids[4 * n + 2] - 1;
    int l = ids[4 * n + 3] - 1;
    double xval = ints[n] * scale_exchange;
    double cval = ints[n] * scale_coulomb;
    for (int f = 0; f < nf; ++f) {
      if (cur_pass == 1) {
        for (int m = 0; m < std::min(4, nmatrix); ++m) {
          add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
          add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
          add4(f3, f, m, i, j, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
          add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, i, j, nf, nmatrix, nbf));
          add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, k, l, nf, nmatrix, nbf));
          add4(f3, f, m, k, l, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));
          add4(f3, f, m, j, i, nf, nmatrix, nbf, cval * get4(d3, f, m, l, k, nf, nmatrix, nbf));
          add4(f3, f, m, l, k, nf, nmatrix, nbf, cval * get4(d3, f, m, j, i, nf, nmatrix, nbf));
        }
        for (int m = 0; m < std::min(7, nmatrix); ++m) {
          add4(f3, f, m, i, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, l, nf, nmatrix, nbf));
          add4(f3, f, m, k, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, j, nf, nmatrix, nbf));
          add4(f3, f, m, i, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, j, k, nf, nmatrix, nbf));
          add4(f3, f, m, l, i, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, j, nf, nmatrix, nbf));
          add4(f3, f, m, j, k, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, l, nf, nmatrix, nbf));
          add4(f3, f, m, k, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, l, i, nf, nmatrix, nbf));
          add4(f3, f, m, j, l, nf, nmatrix, nbf, -xval * get4(d3, f, m, i, k, nf, nmatrix, nbf));
          add4(f3, f, m, l, j, nf, nmatrix, nbf, -xval * get4(d3, f, m, k, i, nf, nmatrix, nbf));
        }
      } else if (cur_pass == 2 && nmatrix >= 7) {
        int m = 6;
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
  }
}

void cpu_umrsf_contract(const std::vector<int>& ids, const std::vector<double>& ints, int ncur,
                        std::vector<double>& f3, const std::vector<double>& d3, int nf,
                        int nmatrix, int nbf, int cur_pass, double scale_exchange,
                        double scale_coulomb) {
  for (int n = 0; n < ncur; ++n) {
    int i = ids[4 * n + 0] - 1;
    int j = ids[4 * n + 1] - 1;
    int k = ids[4 * n + 2] - 1;
    int l = ids[4 * n + 3] - 1;
    double xval = ints[n] * scale_exchange;
    double cval = ints[n] * scale_coulomb;
    for (int f = 0; f < nf; ++f) {
      if (cur_pass == 1) {
        for (int m = 0; m < std::min(8, nmatrix); ++m) {
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
        for (int m = 8; m < std::min(10, nmatrix); ++m) {
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
      if ((cur_pass == 1 || cur_pass == 2) && nmatrix >= 11) {
        int m = 10;
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
  }
}

struct CaseConfig {
  std::string mode;
  int nbf;
  int nf;
  int ncur;
  int pass;
  int seed;
};

double deterministic_uniform(unsigned long long& state, double lo, double hi) {
  state = state * 2862933555777941757ULL + 3037000493ULL;
  double unit = static_cast<double>((state >> 11) & 0x1fffffffffffffULL) / static_cast<double>(0x1fffffffffffffULL);
  return lo + (hi - lo) * unit;
}

int deterministic_int(unsigned long long& state, int lo, int hi) {
  double x = deterministic_uniform(state, 0.0, 1.0);
  int value = lo + static_cast<int>(x * static_cast<double>(hi - lo + 1));
  return std::min(hi, std::max(lo, value));
}

void run_case(const CaseConfig& c) {
  int nmatrix = (c.mode == "umrsf") ? 11 : 7;
  bool is_umrsf = (c.mode == "umrsf");
  unsigned long long rng = static_cast<unsigned long long>(c.seed) + 0x9e3779b97f4a7c15ULL;
  std::vector<int> ids(4 * c.ncur);
  std::vector<double> ints(c.ncur);
  std::size_t tensor_count = static_cast<std::size_t>(c.nf) * nmatrix * c.nbf * c.nbf;
  std::vector<double> d3(tensor_count), f_cpu(tensor_count, 0.0), f_gpu(tensor_count, 0.0);
  for (int n = 0; n < c.ncur; ++n) {
    ids[4*n+0] = deterministic_int(rng, 1, c.nbf);
    ids[4*n+1] = deterministic_int(rng, 1, c.nbf);
    ids[4*n+2] = deterministic_int(rng, 1, c.nbf);
    ids[4*n+3] = deterministic_int(rng, 1, c.nbf);
    ints[n] = deterministic_uniform(rng, -1.0, 1.0) * 0.1;
  }
  for (auto& x : d3) x = deterministic_uniform(rng, -1.0, 1.0) * 0.01;

  double sx = 0.5;
  double sc = 1.0;
  auto cpu0 = std::chrono::high_resolution_clock::now();
  if (is_umrsf) cpu_umrsf_contract(ids, ints, c.ncur, f_cpu, d3, c.nf, nmatrix, c.nbf, c.pass, sx, sc);
  else cpu_mrsf_contract(ids, ints, c.ncur, f_cpu, d3, c.nf, nmatrix, c.nbf, c.pass, sx, sc);
  auto cpu1 = std::chrono::high_resolution_clock::now();
  int ierr = oqp_gpu_metc_contract(ids.data(), ints.data(), c.ncur, f_gpu.data(), d3.data(),
                                   c.nf, nmatrix, c.nbf, c.pass, sx, sc, is_umrsf);
  cudaDeviceSynchronize();
  auto gpu1 = std::chrono::high_resolution_clock::now();

  double max_abs = 0.0, rms = 0.0, denom = 0.0;
  for (std::size_t i = 0; i < tensor_count; ++i) {
    double diff = std::abs(f_cpu[i] - f_gpu[i]);
    max_abs = std::max(max_abs, diff);
    rms += diff * diff;
    denom += f_cpu[i] * f_cpu[i];
  }
  rms = std::sqrt(rms / static_cast<double>(tensor_count));
  double rel = std::sqrt(rms * rms * tensor_count / std::max(denom, 1.0e-300));
  double cpu_ms = std::chrono::duration<double, std::milli>(cpu1 - cpu0).count();
  double gpu_ms = std::chrono::duration<double, std::milli>(gpu1 - cpu1).count();
  double speedup = cpu_ms / std::max(gpu_ms, 1.0e-12);
  std::cout << c.mode << ',' << c.nbf << ',' << c.nf << ',' << nmatrix << ',' << c.ncur << ','
            << c.pass << ',' << c.seed << ',' << ierr << ',' << cpu_ms << ',' << gpu_ms << ','
            << speedup << ',' << max_abs << ',' << rms << ',' << rel << '\n';
}

} // namespace

int main(int argc, char** argv) {
  std::cout << "mode,nbf,nf,nmatrix,ncur,pass,seed,ierr,cpu_ms,gpu_total_ms,speedup_total,max_abs_diff,rms_diff,rel_l2_diff\n";
  std::vector<CaseConfig> cases = {
    {"mrsf", 8, 1, 128, 1, 11}, {"mrsf", 16, 4, 1024, 1, 12},
    {"mrsf", 32, 8, 4096, 1, 13}, {"mrsf", 32, 16, 4096, 2, 14},
    {"mrsf", 64, 16, 8192, 1, 15}, {"umrsf", 16, 4, 1024, 1, 21},
    {"umrsf", 32, 8, 4096, 1, 22}, {"umrsf", 32, 16, 4096, 2, 23},
    {"umrsf", 64, 16, 8192, 1, 24}
  };
  for (const auto& c : cases) run_case(c);
  return 0;
}
