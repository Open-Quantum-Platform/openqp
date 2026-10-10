// M1: device radial integrals g_L^(n)(eta;R) = C_{L,n}(eta) R^L F_{L,n}(T), T=eta R^2.
// F_{L,n}(T)=1F1((n+L+1)/2; L+3/2; -T) is read from a host-tabulated grid
// (tools/radial_gen.py -> radial_table.bin); device uses 4-point (cubic) Lagrange
// interpolation in T plus a large-T asymptotic tail. No special functions or
// cancellation-prone series on device (mirrors the Boys-function device kernel).
#pragma once
#include "sph_eri/port.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstring>

namespace sph_eri {

struct RadialTable {
  int NT = 0, NMAX = 0;
  double TMAX = 0.0, h = 0.0;
  const double* d_F = nullptr;   // [pair*NT + i]
  const int*    d_off = nullptr; // [L*(NMAX+1)+n] -> pair*NT, or -1
  const double* d_K = nullptr;   // [L*(NMAX+1)+n] = sqrt(pi/2) Gamma(a)/(2^b Gamma(b))
                                 // a=(n+L+1)/2, b=L+1.5; the (L,n)-only prefactor
                                 // (precomputed to keep tgamma off the hot loop)
  const float*  d_Ff = nullptr;  // FP32 copies of d_F / d_K so radial_g<float> reads a
  const float*  d_Kf = nullptr;  // native FP32 table (no double loads/casts on the FP32 path)
};

// load radial_table.bin onto the device; returns a RadialTable with device ptrs.
inline RadialTable radial_load(const char* path) {
  FILE* f = std::fopen(path, "rb");
  if (!f) { std::fprintf(stderr, "radial_load: cannot open %s\n", path); std::exit(1); }
  int NT = 0, npairs = 0; double TMAX = 0;
  size_t rd = std::fread(&NT, 4, 1, f); rd += std::fread(&npairs, 4, 1, f);
  rd += std::fread(&TMAX, 8, 1, f); (void)rd;
  // read pairs, find NMAX
  std::vector<int> Ls(npairs), Ns(npairs);
  std::vector<double> F((size_t)npairs * NT);
  int NMAX = 0;
  for (int p = 0; p < npairs; ++p) {
    int L, n; if(std::fread(&L,4,1,f)!=1||std::fread(&n,4,1,f)!=1){std::fprintf(stderr,"radial_load: short read\n");std::exit(1);}
    if (std::fread(&F[(size_t)p*NT], 8, NT, f) != (size_t)NT){std::fprintf(stderr,"radial_load: short read\n");std::exit(1);}
    Ls[p]=L; Ns[p]=n; if (L>NMAX) NMAX=L; if (n>NMAX) NMAX=n;
  }
  std::fclose(f);
  std::vector<int> off((size_t)(NMAX+1)*(NMAX+1), -1);
  for (int p = 0; p < npairs; ++p) off[(size_t)Ls[p]*(NMAX+1)+Ns[p]] = p*NT;
  // (L,n)-only prefactor K = sqrt(pi/2) Gamma(a)/(2^b Gamma(b)); a=(n+L+1)/2, b=L+1.5
  std::vector<double> K((size_t)(NMAX+1)*(NMAX+1), 0.0);
  for (int L=0; L<=NMAX; ++L) for (int n=0; n<=NMAX; ++n) {
    double aa=(n+L+1)*0.5, bb=L+1.5;
    K[(size_t)L*(NMAX+1)+n] = std::sqrt(M_PI/2.0)*std::tgamma(aa)/(std::pow(2.0,bb)*std::tgamma(bb));
  }
  RadialTable t; t.NT=NT; t.NMAX=NMAX; t.TMAX=TMAX; t.h=TMAX/(NT-1);
  // FP32 copies of F and K (so radial_g<float> reads a native FP32 table)
  std::vector<float> Ff(F.size()); for(size_t i=0;i<F.size();++i) Ff[i]=(float)F[i];
  std::vector<float> Kf(K.size()); for(size_t i=0;i<K.size();++i) Kf[i]=(float)K[i];
#ifdef __CUDACC__
  double* dF; int* doff; double* dK; float* dFf; float* dKf;
  cudaMalloc(&dF, F.size()*sizeof(double)); cudaMemcpy(dF, F.data(), F.size()*sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc(&doff, off.size()*sizeof(int)); cudaMemcpy(doff, off.data(), off.size()*sizeof(int), cudaMemcpyHostToDevice);
  cudaMalloc(&dK, K.size()*sizeof(double)); cudaMemcpy(dK, K.data(), K.size()*sizeof(double), cudaMemcpyHostToDevice);
  cudaMalloc(&dFf, Ff.size()*sizeof(float)); cudaMemcpy(dFf, Ff.data(), Ff.size()*sizeof(float), cudaMemcpyHostToDevice);
  cudaMalloc(&dKf, Kf.size()*sizeof(float)); cudaMemcpy(dKf, Kf.data(), Kf.size()*sizeof(float), cudaMemcpyHostToDevice);
  t.d_F=dF; t.d_off=doff; t.d_K=dK; t.d_Ff=dFf; t.d_Kf=dKf;
#else                                   // CPU backend: keep tables in host memory
  double* hF=(double*)std::malloc(F.size()*sizeof(double)); std::memcpy(hF,F.data(),F.size()*sizeof(double));
  int* hoff=(int*)std::malloc(off.size()*sizeof(int)); std::memcpy(hoff,off.data(),off.size()*sizeof(int));
  double* hK=(double*)std::malloc(K.size()*sizeof(double)); std::memcpy(hK,K.data(),K.size()*sizeof(double));
  float* hFf=(float*)std::malloc(Ff.size()*sizeof(float)); std::memcpy(hFf,Ff.data(),Ff.size()*sizeof(float));
  float* hKf=(float*)std::malloc(Kf.size()*sizeof(float)); std::memcpy(hKf,Kf.data(),Kf.size()*sizeof(float));
  t.d_F=hF; t.d_off=hoff; t.d_K=hK; t.d_Ff=hFf; t.d_Kf=hKf;
#endif
  return t;
}

template <typename T>
__host__ __device__ inline T radial_g(const RadialTable& t, int L, int n, T eta, T R) {
  T Tv = eta * R * R;
  T a  = (T)(n + L + 1) * (T)0.5;
  T b  = (T)L + (T)1.5;
  T RL = (T)1; for(int kk=0;kk<L;++kk) RL*=R;       // R^L by integer power
  // pref = K(L,n) * (4 eta)^a * R^L ; K precomputed (no tgamma on the hot path)
  T pref;
  if (L<=t.NMAX && n<=t.NMAX) {
    T Kln = (sizeof(T)==4 && t.d_Kf) ? (T)t.d_Kf[L*(t.NMAX+1)+n] : (T)t.d_K[L*(t.NMAX+1)+n];
    pref = Kln * (T)pow(4.0*(double)eta,(double)a) * RL;
  }
  else                                               // outside table: full formula
    pref = (T)sqrt(M_PI/2.0)*(T)tgamma((double)a)*(T)pow(4.0*(double)eta,(double)a)
           / ((T)pow(2.0,(double)b)*(T)tgamma((double)b)) * RL;
  T F;
  if (Tv >= (T)t.TMAX) {                 // large-T asymptotic: Gamma(b)/Gamma(b-a) T^-a
    double ba = (double)b - (double)a;
    double coef = (ba <= 0.0 && fabs(ba-round(ba))<1e-12) ? 0.0 : tgamma((double)b)/tgamma(ba);
    F = (T)(coef * pow((double)Tv, -(double)a));
  } else {
    if (L>t.NMAX || n>t.NMAX) return (T)0;          // outside table coverage
    int off = t.d_off[L*(t.NMAX+1)+n];
    if (off < 0) return (T)0;                         // (L,n) not tabulated
    int i = (int)(Tv / t.h); if (i<1) i=1; if (i>t.NT-3) i=t.NT-3;
    T acc = 0;
    #pragma unroll
    for (int j=0;j<4;++j){
      int ij = i-1+j; T xj = (T)ij*(T)t.h;
      T term = (sizeof(T)==4 && t.d_Ff) ? (T)t.d_Ff[off+ij] : (T)t.d_F[off+ij];
      for (int k=0;k<4;++k) if(k!=j){ T xk=(T)(i-1+k)*(T)t.h; term *= (Tv-xk)/(xj-xk); }
      acc += term;
    }
    F = acc;
  }
  return pref * F;
}

} // namespace sph_eri
