// =============================================================================
// routec_oqp_jk2_gpu.cu — OpenQP Route-C GPU backend v2 (APIV2 / R3b).
//
// One dylib, two seams:
//   1. routec_fock_jk            — v1 SCF seam (packed lower-tri D in, packed
//                                  F2e = scale_coul*J - 0.5*scale_exch*K out),
//                                  superset of routec_oqp_jk_gpu.cu v1: also
//                                  accepts packed B files and nfocks >= 1.
//   2. routec_jkm_init/free/apply — APIV2 batched RAW J(M)/K(M) for arbitrary
//                                  NON-symmetric M (response bridge). All
//                                  consumer prefactors live at the Fortran
//                                  seam (APIV2.md), never here.
//   2b. routec_jkm_apply_lr      — APIV2.1 low-rank fast path: J/K of
//                                  M_o = sum_r u_r v_r^T from explicit factor
//                                  columns (MRSF d3 slots 1-4 rank-1, 5-6
//                                  rank-2). Drops the nbf^3 K term to
//                                  ~6*naux*nbf^2 per rank-1 output.
//
// Math (B = whitened 3c tensor, every aux slice B_P SYMMETRIC in (a,b)):
//   J(M)_ab = (ab|kl) M_kl = sum_P B_P[a,b] g_P,  g_P = sum_kl B_P[k,l] M_kl
//   K(M)_ab = (ak|bl) M_kl = sum_P (B_P · M · B_P^T)[a,b]
//           = sum_P (B_P · M · B_P)[a,b]                 (B_P = B_P^T)
//   NO symmetrization anywhere: K(M^T) = K(M)^T != K(M) is load-bearing
//   (the 0.31 Ha transposed-K control in APIV2.md). J needs no M_sym either:
//   g_P = B_P : M = B_P : M_sym exactly, because B_P is symmetric.
//
// Layout discipline (the #1 risk — derived twice at every GEMM below):
//   * device dB: C-order (naux, nbf, nbf): dB[P*nn + i*nbf + j] = B_P[i,j].
//     - col-major (nbf x nbf) view of dB+P*nn, ld=nbf: element (r,c) is at
//       c*nbf + r = B_P[c,r] = B_P[r,c] by slice symmetry  =>  view == B_P.
//     - col-major (nbf x naux*nbf) view "Bview" of dB, ld=nbf: column
//       s = P*nbf + l, row r is at s*nbf + r = P*nn + l*nbf + r
//       =>  Bview[r, (P,l)] = B_P[l,r]                       (NO symmetry used)
//     - col-major (nn x naux) view "Bflat" of dB, ld=nn: Bflat[z,P] =
//       dB[P*nn + z] = B_P[i,j] with z = i*nbf + j.
//   * caller x: Fortran X(nbf,nbf,nvec), X(a,b,v) = M_v[a,b], a fastest:
//     linear v*nn + b*nbf + a. The col-major (nbf x nbf) view of dX+v*nn,
//     ld=nbf, IS M_v exactly (no transpose — cuBLAS and Fortran agree).
//   * outputs dJ/dK mirror x's layout: buffer[v*nn + b*nbf + a] = J/K[a,b].
//
// B-file strategy (choice APIV2 leaves open): packed B files are unpacked to
// dense ON THE HOST at load time (verbatim the CPU-v2 loader loop from
// routec_oqp_jk.cpp), then uploaded once; device keeps DENSE B for the whole
// session. Rationale: every GEMM needs dense slices anyway, host-unpack
// reuses the already-validated CPU code path, and the one-time extra PCIe
// traffic (2x packed size, ~0.1 s at 2 GB over ~25 GB/s) is irrelevant
// against Davidson-iteration reuse. A device-side unpack kernel would save
// host RAM only and adds an untestable-on-this-Mac kernel.
//
// Build (chc4, A100):  module load CUDA/12.6.0 GCC/12.3.0   (build AND run)
//   nvcc -O3 -arch=sm_80 -std=c++17 -Xcompiler -fPIC -shared \
//        routec_oqp_jk2_gpu.cu -lcublas -o libroutec_oqp_gpu2.so
//
// Env:
//   OQP_ROUTEC_B      B file (v1 format; packed header variants per CPU v2)
//   OQP_JKM_PTRARRAY  =1: enable the pointer-array (P,v)-fold GEMM path
//   OQP_JKM_VCHUNK    K-intermediate chunk (vectors per chunk; implies fold)
//   OQP_JKM_LOOPV     =1: force vchunk=1 (kept for compatibility = default)
// Default is the per-vector strided-batched path (vchunk=1): measured on
// A100 (chc4 2026-06-12) it is uniformly >= the pointer-array fold and
// perfectly linear in nvec — 3.34/30.9/153/450 ms/vec at w8/w16/w24/w32,
// 75% of FP64-TC peak at w32; the fold pays host pointer-array builds and
// extra H2D per chunk and degrades when the chunk shrinks.
// =============================================================================
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <algorithm>
#include <vector>
#include <sys/stat.h>
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include "routec_df.h"   // shared resident B (routec_b_shared)

// bool-returning guards: callers translate `false` into info!=0 / rc!=0.
#define CUTRY(x) do{ cudaError_t e_=(x); if(e_!=cudaSuccess){ \
  fprintf(stderr,"[routec-gpu2] CUDA error %s at line %d\n", \
          cudaGetErrorString(e_),__LINE__); return false; } }while(0)
#define CBTRY(x) do{ cublasStatus_t s_=(x); if(s_!=CUBLAS_STATUS_SUCCESS){ \
  fprintf(stderr,"[routec-gpu2] cuBLAS error %d at line %d\n", \
          (int)s_,__LINE__); return false; } }while(0)

namespace {

// ------------------------------------------------------------ session state
double* dB = nullptr;            // dense (naux, nbf, nbf) C-order, device
bool    b_borrowed = false;      // dB is owned by df_store (shared) -- don't free
int g_naux = 0, g_nbf = 0;
cublasHandle_t hb = nullptr;
// CDF compacted-on-device state (set by load_B when OQP_CDF_ONDEV): B held as
// (naux x ncp) compacted; J/K tile-reconstruct dense aux-slices on the fly.
bool    g_cdf_on   = false;
const double* g_cdf_cB   = nullptr;   // fp64 arena (v2: naux x ncp; v3: n64 x ncp)
const int*    g_cdf_keep = nullptr;   // (ncp) lower-tri pair indices
int     g_cdf_ncp  = 0;
// CDF v3 (magnitude-tiered mixed precision): per-column arena dispatch.
int           g_cdf_lowprec = 0;
const void*   g_cdf_B32  = nullptr;   // float  arena (n32 x ncp)
const void*   g_cdf_B16  = nullptr;   // __half arena (n16 x ncp)
const void*   g_cdf_B8   = nullptr;   // int8   arena (n8  x ncp)
const unsigned char* g_cdf_tier = nullptr;  // (naux) per-column tier
const int*    g_cdf_slot  = nullptr;  // (naux) row within arena
const float*  g_cdf_scale = nullptr;  // (naux) per-column dequant scale
double* dBtile = nullptr; size_t capBtile = 0;   // per-tile dense scratch
static const int CDF_TILE = getenv("OQP_CDF_TILE") ? atoi(getenv("OQP_CDF_TILE")) : 256;  // aux per tile (bigger = fewer scatter/GEMM launches, more tile scratch)

// ---- FP32-K (KP step-2, the real K-speed lever) ----------------------------
// Optional single-precision exchange (K) path, gated by env OQP_JKM_FP32K=1 OR
// flipped at runtime by the adaptive SCF ramp (set_fp32k below). RI-K is
// ~95-98% of the DF J+K apply, so this is the speed lever. dBf is a float
// mirror of dB built once at load; the half-transform W=B_P·M and the Gram
// K=W·B^T run as cublasSgemmStridedBatched, K cast back to double on output.
// FP32-K is safe while the SCF density is rough (loose DIIS); the ramp switches
// to FP64-K once the DIIS commitator drops below a threshold.
bool   g_fp32k = false;          // load_B() reads OQP_JKM_FP32K; ramp toggles it
float* dBf  = nullptr;           // dense (naux, nbf, nbf) C-order float mirror
float* dWf  = nullptr;           // float K intermediate: wf_vchunk * naux * nn
float* dXkf = nullptr;           // float density slices for K (nn*nvec)
float* dKf  = nullptr;           // float K output (nn*nvec) before cast
int    wf_vchunk = 0;            // vectors per FP32 K chunk currently allocated
size_t capXkf = 0, capKf = 0;    // element capacities (floats)

// elementwise double<->float casts (FP32-K).
__global__ void cast_d2f(const double* in, float* out, long long n) {
    long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) out[i] = (float)in[i];
}
__global__ void cast_f2d(const float* in, double* out, long long n) {
    long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) out[i] = (double)in[i];
}

// workspaces (allocated max-once, grown on demand, freed by jkm_free)
double *dX = nullptr, *dJ = nullptr, *dK = nullptr, *dG = nullptr;
size_t capX = 0, capJ = 0, capK = 0, capG = 0;       // element capacities
double *dW = nullptr;            // K intermediate: w_vchunk * naux * nn
int w_vchunk = 0;                // vectors per K chunk currently allocated
size_t w_elems = 0;              // dW capacity in doubles (aux-tiled K needs only KT*nn)
// occ-K (#145): half-transform straight from the OCC factor Cocc (nbf x nocc)
// instead of the full density D. K(D) = 2 sum_P (B_P Cocc)(B_P Cocc)^T, cost
// O(naux*nocc*nbf^2) -- ~nbf/nocc (~5x) FEWER FLOPs than jkm_core's full-D K
// (W = B_P*D, O(naux*nbf^3)). dMocc holds the half-transform Mf^T (col-major
// nbf x naux*nocc): column (P,o) = (B_P * Cocc)[:,o]; K = dMocc * dMocc^T.
double *dMocc = nullptr; size_t capMocc = 0;     // FP64 occ half-transform
float  *dMoccf = nullptr; size_t capMoccf = 0;   // FP32-K twin
float  *dCof = nullptr;  size_t capCof = 0;      // FP32 Cocc factor
const double **dAarr = nullptr, **dBarr = nullptr;   // pointer-array GEMM
double **dCarr = nullptr;
int cap_ptr = 0;
double *ddp = nullptr, *dfp = nullptr;               // fock_jk packed buffers
size_t capP = 0, capF = 0;
double *dU = nullptr, *dV = nullptr;                 // lr factor columns
double *dYZ = nullptr, *dGlr = nullptr;              // lr Y/Z stack + J g's
size_t capU = 0, capV = 0, capYZ = 0, capGlr = 0;
std::vector<const double*> hA, hBp;                  // host staging for arrays
std::vector<double*> hC;

bool grow(double** p, size_t* cap, size_t need_elems) {
    if (need_elems <= *cap) return true;
    if (*p) cudaFree(*p);
    *p = nullptr; *cap = 0;
    CUTRY(cudaMalloc((void**)p, need_elems*sizeof(double)));
    *cap = need_elems;
    return true;
}

// float-buffer twin of grow() for the FP32-K workspaces.
bool growf(float** p, size_t* cap, size_t need_elems) {
    if (need_elems <= *cap) return true;
    if (*p) cudaFree(*p);
    *p = nullptr; *cap = 0;
    CUTRY(cudaMalloc((void**)p, need_elems*sizeof(float)));
    *cap = need_elems;
    return true;
}

// ---- multi-GPU (OQP_MULTI_GPU): B aux-split across two devices --------------
// dev0 holds B rows [0,n0), dev1 holds [n0,naux). Per cycle each device
// contracts its own slice (J partial nn*nvec, K partial nn) and dev1's partial
// is peer-copied to dev0 and added -- a few MB over NVLink, vs GBs of B reads
// saved. Device-1 GEMMs run on their own stream, CONCURRENT with device 0.
bool g_mg = false; int g_mgn0 = 0, g_mgn1 = 0;
double* dB1 = nullptr;                       // borrowed dev1 slice (df_store owns)
cublasHandle_t hb1 = nullptr; cudaStream_t s1 = nullptr;
double *dX1=nullptr, *dG1=nullptr, *dJ1=nullptr, *dCo1=nullptr,
       *dMocc1=nullptr, *dK1=nullptr, *dW1=nullptr, *dPeer=nullptr;
size_t capX1=0,capG1=0,capJ1=0,capCo1=0,capMocc1=0,capK1=0,capW1=0,capPeer=0;

// grow() on device 1 (temporarily switches the current device).
bool grow1(double** p, size_t* cap, size_t need_elems) {
    if (need_elems <= *cap) return true;
    CUTRY(cudaSetDevice(1));
    if (*p) cudaFree(*p);
    *p = nullptr; *cap = 0;
    if (cudaMalloc((void**)p, need_elems*sizeof(double)) != cudaSuccess) {
        cudaSetDevice(0); return false;
    }
    *cap = need_elems;
    CUTRY(cudaSetDevice(0));
    return true;
}

// build the float B mirror once (lazy: the SCF ramp can enable FP32-K after
// load_B already ran with the env unset). dB must already be resident.
bool ensure_Bf() {
    if (dBf) return true;
    if (!dB) return false;
    const long long nB = (long long)g_naux*g_nbf*g_nbf;
    CUTRY(cudaMalloc((void**)&dBf, (size_t)nB*sizeof(float)));
    cast_d2f<<<(unsigned)((nB+255)/256), 256>>>(dB, dBf, nB);
    CUTRY(cudaGetLastError());
    fprintf(stderr, "[routec-gpu2] FP32 B mirror resident (%.1f MB) — K=FP32\n",
            (size_t)nB*4.0/1048576.0);
    return true;
}

bool ensure_ptr_arrays(int need) {
    if (need <= cap_ptr) return true;
    if (dAarr) cudaFree((void*)dAarr);
    if (dBarr) cudaFree((void*)dBarr);
    if (dCarr) cudaFree((void*)dCarr);
    dAarr = dBarr = nullptr; dCarr = nullptr; cap_ptr = 0;
    CUTRY(cudaMalloc((void**)&dAarr, (size_t)need*sizeof(double*)));
    CUTRY(cudaMalloc((void**)&dBarr, (size_t)need*sizeof(double*)));
    CUTRY(cudaMalloc((void**)&dCarr, (size_t)need*sizeof(double*)));
    cap_ptr = need;
    return true;
}

// K intermediate W: ideally nvec*naux*nn doubles (full (P,v) fold), but that
// is B-sized PER VECTOR — chunk by available memory, halving on OOM.
bool ensure_w(int nvec) {
    const size_t per_v = (size_t)g_naux*g_nbf*g_nbf;       // doubles/vector
    int want = 1;                                  // strided per-vector default
    if (getenv("OQP_JKM_PTRARRAY")) want = nvec;   // opt into the (P,v) fold
    if (const char* ev = getenv("OQP_JKM_VCHUNK")) {
        int u = atoi(ev);
        if (u >= 1) want = std::min(nvec, u);
    }
    if (getenv("OQP_JKM_LOOPV")) want = 1;
    if (dW && w_vchunk >= want) return true;
    if (dW) { cudaFree(dW); dW = nullptr; w_vchunk = 0; w_elems = 0; }
    while (want >= 1) {
        if (cudaMalloc((void**)&dW, (size_t)want*per_v*sizeof(double))
            == cudaSuccess) {
            if (want > 1 && !ensure_ptr_arrays(g_naux*want)) {
                cudaFree(dW); dW = nullptr;   // ptr arrays failed: degrade
                want = 1;
                continue;
            }
            w_vchunk = want;
            w_elems = (size_t)want*per_v;
            return true;
        }
        dW = nullptr;
        want >>= 1;                                        // OOM: halve
    }
    fprintf(stderr, "[routec-gpu2] cannot allocate K workspace "
            "(naux=%d nbf=%d: %.1f MB/vector)\n",
            g_naux, g_nbf, per_v*8.0/1048576.0);
    return false;
}

// aux-tiled full-D K workspace: only KT*nbf^2 doubles (vs naux*nbf^2 -- the
// difference is ~16 GB at (H2O)32). Reuses a bigger already-resident dW.
bool ensure_w_ktile(int KT) {
    const size_t need = (size_t)KT*g_nbf*g_nbf;
    if (dW && w_elems >= need) return true;
    if (dW) { cudaFree(dW); dW = nullptr; w_vchunk = 0; w_elems = 0; }
    CUTRY(cudaMalloc((void**)&dW, need*sizeof(double)));
    w_elems = need; w_vchunk = 0;   // not sized for the (P,v)-fold path
    return true;
}

// FP32-K intermediate: float twin of ensure_w() (strided-batched only, so a
// per-vector float W is enough — half the FP64 footprint).
bool ensure_wf(int nvec) {
    const size_t per_v = (size_t)g_naux*g_nbf*g_nbf;       // floats/vector
    int want = 1;
    if (const char* ev = getenv("OQP_JKM_VCHUNK")) {
        int u = atoi(ev);
        if (u >= 1) want = std::min(nvec, u);
    }
    if (getenv("OQP_JKM_LOOPV")) want = 1;
    if (dWf && wf_vchunk >= want) return true;
    if (dWf) { cudaFree(dWf); dWf = nullptr; wf_vchunk = 0; }
    while (want >= 1) {
        if (cudaMalloc((void**)&dWf, (size_t)want*per_v*sizeof(float))
            == cudaSuccess) { wf_vchunk = want; return true; }
        dWf = nullptr;
        want >>= 1;
    }
    fprintf(stderr, "[routec-gpu2] cannot allocate FP32 K workspace "
            "(naux=%d nbf=%d: %.1f MB/vector)\n",
            g_naux, g_nbf, per_v*4.0/1048576.0);
    return false;
}

// ------------------------------------------------------------------ B loader
// EXACT mirror of the CPU-v2 loader in routec_oqp_jk.cpp (header variants +
// file-size validation), with the dense tensor uploaded to the device and
// the host copy released.
bool load_B() {
    if (g_cdf_on) return true;                    // compacted already resident
    if (dB) return true;
    // OQP_CDF_ONDEV: hold the COMPACTED B on device (naux x ncp) and tile J/K
    // from it, so the full dense B never resides on the GPU. Real device-memory
    // win (needs a CDF v2 file from ROUTEC_CDF_V2). HF/DFT path.
    if (getenv("OQP_CDF_ONDEV")) {
        CdfDev cdf;
        if (routec_cdf_dev(&cdf)!=0) {
            fprintf(stderr,"[routec-gpu2] OQP_CDF_ONDEV set but compacted B load failed "
                    "(need a ROUTEC_CDF_V2 file)\n"); return false; }
        g_naux=cdf.naux; g_nbf=cdf.nao; g_cdf_ncp=cdf.ncp; g_cdf_cB=cdf.B64; g_cdf_keep=cdf.keep;
        g_cdf_lowprec=cdf.lowprec; g_cdf_B32=cdf.B32; g_cdf_B16=cdf.B16; g_cdf_B8=cdf.B8;
        g_cdf_tier=cdf.tier; g_cdf_slot=cdf.slot; g_cdf_scale=cdf.scale;
        g_cdf_on=true;
        if (!hb) CBTRY(cublasCreate(&hb));
        fprintf(stderr,"[routec-gpu2] CDF on-device: naux=%d nbf=%d ncp=%d %s -> tiled J/K\n",
                cdf.naux,cdf.nao,cdf.ncp, cdf.lowprec?"(mixed fp64/fp32/fp16/int8)":"(fp64)");
        return true;
    }
    // Shared resident B: df_store loads OQP_ROUTEC_B once (dense/packed handled
    // there) and hands back a BORROWED device pointer, so the MRSF sigma
    // session reuses this same upload instead of reading the file again.
    int naux = 0, nbf = 0;
    const double* B = routec_b_shared(&naux, &nbf);
    if (!B) { fprintf(stderr, "[routec-gpu2] shared B load failed\n"); return false; }
    g_naux = naux; g_nbf = nbf;
    const size_t nn = (size_t)g_nbf*g_nbf;
    if ((long long)g_naux*g_nbf > (long long)INT_MAX) {  // K step-2 k-dim guard
        fprintf(stderr, "[routec-gpu2] naux*nbf exceeds INT_MAX\n");
        return false;
    }
    dB = const_cast<double*>(B);         // borrowed from df_store -- do NOT free
    b_borrowed = true;
    if (!hb) CBTRY(cublasCreate(&hb));
    fprintf(stderr, "[routec-gpu2] B resident (shared): naux=%d nbf=%d (%.1f MB dense)\n",
            g_naux, g_nbf, (size_t)g_naux*nn*8.0/1048576.0);
    // ---- multi-GPU: split B across two devices (aux dimension) --------------
    if (getenv("OQP_MULTI_GPU") && atoi(getenv("OQP_MULTI_GPU")) != 0) {
        int n0=0, n1=0, nb2=0; const double *B0=nullptr, *B1p=nullptr;
        if (routec_b_split_mg(&n0, &n1, &nb2, &B0, &B1p) == 1) {
            dB = const_cast<double*>(B0);          // dev0 half (df_store shrank it)
            dB1 = const_cast<double*>(B1p);        // dev1 half
            g_mgn0 = n0; g_mgn1 = n1; g_mg = true;
            CUTRY(cudaSetDevice(1));
            CBTRY(cublasCreate(&hb1));
            CUTRY(cudaStreamCreate(&s1));
            CBTRY(cublasSetStream(hb1, s1));
            CUTRY(cudaSetDevice(0));
            fprintf(stderr, "[routec-gpu2] multi-GPU J/K: dev0 %d + dev1 %d aux slices\n", n0, n1);
        } else {
            fprintf(stderr, "[routec-gpu2] OQP_MULTI_GPU set but split unavailable; single-GPU\n");
        }
    }
    // FP32-K: if requested via env, build the float B mirror now.
    if (const char* e = getenv("OQP_JKM_FP32K")) g_fp32k = (atoi(e) != 0);
    if (g_mg && g_fp32k) {
        fprintf(stderr, "[routec-gpu2] FP32-K disabled under multi-GPU (unsupported combo)\n");
        g_fp32k = false;
    }
    if (g_fp32k && !ensure_Bf()) return false;
    return true;
}

// ------------------------------------------------------------------- kernels
// unpack nslices packed lower-triangle matrices (OpenQP row-walk i>=j,
// t = i*(i+1)/2 + j) into full squares. Output is SYMMETRIC, so col-major
// (b*nbf+a) vs row-major orientation is moot — written col-major for clarity.
__global__ void unpack_tri_multi(const double* dp, double* X, int nbf,
                                 int nslices) {
    long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    const long long nn = (long long)nbf*nbf;
    const long long ntri = (long long)nbf*(nbf+1)/2;
    if (i >= nn*nslices) return;
    const int s = (int)(i / nn);
    const long long z = i % nn;                 // z = b*nbf + a (col-major)
    const int a = (int)(z % nbf), b = (int)(z / nbf);
    const int hi = a > b ? a : b, lo = a > b ? b : a;
    X[i] = dp[(long long)s*ntri + (long long)hi*(hi+1)/2 + lo];
}

// pack nslices: f[s][t(i,j)] = cj*J[a=i,b=j] - cx*K[i,j]. For symmetric D
// both J(D) and K(D) are symmetric, so reading the (i,j) col-major element
// (linear j*nbf+i) is orientation-safe. Row index from t via double sqrt
// (exact for nbf <= ~46000) + two safety loops.
__global__ void pack_tri_multi(const double* J, const double* K, double* f,
                               int nbf, int nslices, double cj, double cx) {
    long long t = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    const long long ntri = (long long)nbf*(nbf+1)/2;
    const long long nn = (long long)nbf*nbf;
    if (t >= ntri*nslices) return;
    const int s = (int)(t / ntri);
    const long long u = t % ntri;
    int i = (int)((sqrt(8.0*(double)u + 1.0) - 1.0)*0.5);
    while ((long long)(i+1)*(i+2)/2 <= u) ++i;
    while ((long long)i*(i+1)/2 > u) --i;
    const int j = (int)(u - (long long)i*(i+1)/2);
    const long long z = (long long)s*nn + (long long)j*nbf + i;
    f[t] = cj*J[z] - cx*K[z];
}

// ---- CDF compacted-on-device: reconstruct dense aux-slice tiles on the fly ---
// When OQP_CDF_ONDEV is set, load_B holds the COMPACTED B (naux x ncp) on device
// (via routec_cdf_dev) instead of the full dense (naux x nbf^2). J and K then
// loop over aux TILES, rebuilding a small dense B-tile per tile, so the full
// dense B never resides on the GPU -- the real device-memory win.
// (g_cdf_on / g_cdf_cB / g_cdf_keep / g_cdf_ncp / dBtile / CDF_TILE declared
// near dB above.)
__global__ void cdf_scatter_slices(const double* cB, const int* keep, int ncp,
                                   int nao, int nt, double* dBt) {
    const long long idx = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= (long long)nt*ncp) return;
    const int a = (int)(idx / ncp), c = (int)(idx % ncp);
    const long long t = keep[c];
    int i = (int)((sqrt(8.0*(double)t+1.0)-1.0)*0.5);   // unpack tri t=i*(i+1)/2+j
    while ((long long)(i+1)*(i+2)/2 <= t) ++i;
    while ((long long)i*(i+1)/2 > t) --i;
    const int j = (int)(t - (long long)i*(i+1)/2);
    const double v = cB[(long long)a*ncp + c];
    const long long base = (long long)a*nao*nao;
    dBt[base + (long long)i*nao + j] = v;
    dBt[base + (long long)j*nao + i] = v;
}
// CDF v3: reconstruct a dense B-tile from the MIXED-precision arenas, upcasting
// each aux column P=p0+a from its tier (fp64/fp16/int8) to fp64 -- so every
// downstream GEMM stays fp64 and UNCHANGED; only the resident storage shrank.
__global__ void cdf_scatter_slices_mp(const double* B64, const float* B32, const __half* B16,
                                      const signed char* B8, const unsigned char* tier,
                                      const int* slot, const float* scale,
                                      const int* keep, int ncp, int nao, int p0, int nt, double* dBt) {
    const long long idx = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= (long long)nt*ncp) return;
    const int a = (int)(idx / ncp), c = (int)(idx % ncp);
    const int P = p0 + a;
    const long long t = keep[c];
    int i = (int)((sqrt(8.0*(double)t+1.0)-1.0)*0.5);
    while ((long long)(i+1)*(i+2)/2 <= t) ++i;
    while ((long long)i*(i+1)/2 > t) --i;
    const int j = (int)(t - (long long)i*(i+1)/2);
    const int sl = slot[P];
    double v;
    switch (tier[P]) {                                       // 0 fp64 1 fp32 2 fp16 3 int8
        case 0:  v = B64[(long long)sl*ncp + c]; break;
        case 1:  v = (double)B32[(long long)sl*ncp + c]; break;
        case 2:  v = (double)__half2float(B16[(long long)sl*ncp + c]) * (double)scale[P]; break;
        default: v = (double)B8[(long long)sl*ncp + c] * (double)scale[P]; break;
    }
    const long long base = (long long)a*nao*nao;
    dBt[base + (long long)i*nao + j] = v;
    dBt[base + (long long)j*nao + i] = v;
}
// build the dense B-tile for aux [p0, p0+nt): zero then scatter the compacted rows
static bool cdf_tile(int p0, int nt, int nao) {
    const long long nn = (long long)nao*nao;
    if (!grow(&dBtile, &capBtile, (size_t)nt*nn)) return false;
    if (cudaMemset(dBtile, 0, (size_t)nt*nn*sizeof(double)) != cudaSuccess) return false;
    const long long tot = (long long)nt*g_cdf_ncp; const int TB = 256;
    if (g_cdf_lowprec)
        cdf_scatter_slices_mp<<<(unsigned)((tot+TB-1)/TB), TB>>>(
            g_cdf_cB, (const float*)g_cdf_B32, (const __half*)g_cdf_B16, (const signed char*)g_cdf_B8,
            g_cdf_tier, g_cdf_slot, g_cdf_scale,
            g_cdf_keep, g_cdf_ncp, nao, p0, nt, dBtile);
    else
        cdf_scatter_slices<<<(unsigned)((tot+TB-1)/TB), TB>>>(
            g_cdf_cB + (long long)p0*g_cdf_ncp, g_cdf_keep, g_cdf_ncp, nao, nt, dBtile);
    return cudaGetLastError() == cudaSuccess;
}

// ----------------------------------------------------------------- core math
// dXin : device, nvec col-major (nbf x nbf) slices, slice v = M_v exactly.
// dJout/dKout : device, same layout (either may be null). Pure device work;
// callers do all H2D/D2H. Returns false on any CUDA/cuBLAS failure.
bool jkm_core(const double* dXin, int nvec, double* dJout, double* dKout) {
    const int nbf = g_nbf, naux = g_naux;
    const long long nn = (long long)nbf*nbf;
    const double one = 1.0, zero = 0.0;

    if (dJout) {
        if (!grow(&dG, &capG, (size_t)naux*nvec)) return false;
      if (g_mg) {
        // ---- dual-GPU J: each device contracts its own aux slice ------------
        // dev1 (async on s1): G1 = B1^T X, J1 = B1 G1; dev0 concurrently the
        // same on its half; then J += peer(J1). Same algebra split by P.
        const long long nX = nn*nvec;
        if (!grow1(&dX1,&capX1,(size_t)nX)) return false;
        if (!grow1(&dG1,&capG1,(size_t)g_mgn1*nvec)) return false;
        if (!grow1(&dJ1,&capJ1,(size_t)nX)) return false;
        if (!grow(&dPeer,&capPeer,(size_t)nX)) return false;
        CUTRY(cudaMemcpyPeer(dX1,1,dXin,0,(size_t)nX*8));
        CUTRY(cudaSetDevice(1));
        CBTRY(cublasDgemm(hb1, CUBLAS_OP_T, CUBLAS_OP_N,
                          g_mgn1, nvec, (int)nn, &one,
                          dB1, (int)nn, dX1, (int)nn, &zero, dG1, g_mgn1));
        CBTRY(cublasDgemm(hb1, CUBLAS_OP_N, CUBLAS_OP_N,
                          (int)nn, nvec, g_mgn1, &one,
                          dB1, (int)nn, dG1, g_mgn1, &zero, dJ1, (int)nn));
        CUTRY(cudaSetDevice(0));
        CBTRY(cublasDgemm(hb, CUBLAS_OP_T, CUBLAS_OP_N,
                          g_mgn0, nvec, (int)nn, &one,
                          dB, (int)nn, dXin, (int)nn, &zero, dG, g_mgn0));
        CBTRY(cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_N,
                          (int)nn, nvec, g_mgn0, &one,
                          dB, (int)nn, dG, g_mgn0, &zero, dJout, (int)nn));
        CUTRY(cudaSetDevice(1)); CUTRY(cudaStreamSynchronize(s1)); CUTRY(cudaSetDevice(0));
        CUTRY(cudaMemcpyPeer(dPeer,0,dJ1,1,(size_t)nX*8));
        CBTRY(cublasDaxpy(hb,(int)nX,&one,dPeer,1,dJout,1));
      } else if (g_cdf_on) {
        // ---- tiled J from compacted B (no full dense B on device) ----------
        for (int p0=0;p0<naux;p0+=CDF_TILE){ const int nt=std::min(CDF_TILE,naux-p0);
            if(!cdf_tile(p0,nt,nbf)) return false;                    // [J1] g[p0..]=B_tile^T·X
            CBTRY(cublasDgemm(hb,CUBLAS_OP_T,CUBLAS_OP_N, nt,nvec,(int)nn,&one,
                dBtile,(int)nn, dXin,(int)nn, &zero, dG+p0, naux)); }
        for (int p0=0;p0<naux;p0+=CDF_TILE){ const int nt=std::min(CDF_TILE,naux-p0);
            if(!cdf_tile(p0,nt,nbf)) return false;                    // [J2] J += B_tile·g[p0..]
            const double bta=(p0==0)?0.0:1.0;
            CBTRY(cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N, (int)nn,nvec,nt,&one,
                dBtile,(int)nn, dG+p0,naux, &bta, dJout,(int)nn)); }
      } else {
        // ---- [J1]  G = Bflat^T · Xall      (one GEMM for ALL vectors)
        // Col-major algebra: Bflat = (nn x naux) view of dB (ld=nn),
        // Xall = (nn x nvec) view of dXin (ld=nn, columns ARE the x slices
        // because slice stride == nn == ld). op(A)=T: (naux x nn)·(nn x nvec)
        // -> G (naux x nvec), ldc=naux.
        // Element sum: G[P,v] = sum_z dB[P*nn+z] * dXin[v*nn+z]. With
        // z = i*nbf+j: dB term = B_P[i,j]; dXin linear v*nn + b*nbf + a
        // matches z with b=i, a=j, i.e. M_v[j,i]. So
        // G[P,v] = sum_ij B_P[i,j] M_v[j,i] = sum_ij B_P[j,i] M_v[j,i]
        //        = sum_kl B_P[k,l] M_v[k,l] = g_{P,v}.   (slice symmetry)
        CBTRY(cublasDgemm(hb, CUBLAS_OP_T, CUBLAS_OP_N,
                          naux, nvec, (int)nn, &one,
                          dB, (int)nn,            // A: (nn x naux), ld nn
                          dXin, (int)nn,          // B: (nn x nvec), ld nn
                          &zero, dG, naux));      // C: (naux x nvec)
        // ---- [J2]  Jall = Bflat · G
        // Col-major algebra: (nn x naux)·(naux x nvec) -> (nn x nvec), ld nn.
        // Element sum: dJout[v*nn+z] = sum_P dB[P*nn+z]*G[P,v]
        //   = sum_P B_P[i,j] g_{P,v} = J_v[i,j]  (z = i*nbf+j).
        // Fortran reads jm(a,b,v) at v*nn + b*nbf + a = z with i=b, j=a,
        // i.e. J_v[b,a] = J_v[a,b] by J symmetry. RAW J, no scaling.
        CBTRY(cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_N,
                          (int)nn, nvec, naux, &one,
                          dB, (int)nn, dG, naux,
                          &zero, dJout, (int)nn));
      }
    }

    if (dKout && g_mg) {
        // ---- dual-GPU full-D K (guess path): each device tiles its OWN aux
        // slice (same aux-tiled algebra as the default path below), dev1 into
        // dJ1 (reused as its K accumulator -- the J branch above has already
        // consumed it), then K += peer(K1).
        static int s_ktile_mg = getenv("OQP_JKM_KTILE") ? atoi(getenv("OQP_JKM_KTILE")) : 256;
        const int KT = std::max(1, std::min(s_ktile_mg > 0 ? s_ktile_mg : 256,
                                            std::max(g_mgn0, g_mgn1)));
        if (!ensure_w_ktile(KT)) return false;
        if (!grow1(&dW1,&capW1,(size_t)KT*nn)) return false;
        const long long nX = nn*nvec;
        if (!grow1(&dX1,&capX1,(size_t)nX)) return false;
        if (!grow1(&dJ1,&capJ1,(size_t)nX)) return false;
        if (!grow(&dPeer,&capPeer,(size_t)nX)) return false;
        CUTRY(cudaMemcpyPeer(dX1,1,dXin,0,(size_t)nX*8));
        CUTRY(cudaSetDevice(1));
        for (int v = 0; v < nvec; ++v) {
            for (int t0 = 0; t0 < g_mgn1; t0 += KT) {
                const int nt = std::min(KT, g_mgn1 - t0);
                CBTRY(cublasDgemmStridedBatched(hb1, CUBLAS_OP_N, CUBLAS_OP_N,
                    nbf, nbf, nbf, &one,
                    dB1 + (size_t)t0*nn, nbf, nn,
                    dX1 + (size_t)v*nn, nbf, 0LL,
                    &zero, dW1, nbf, nn, nt));
                const double bta = (t0 == 0) ? 0.0 : 1.0;
                CBTRY(cublasDgemm(hb1, CUBLAS_OP_N, CUBLAS_OP_T,
                    nbf, nbf, (int)((long long)nt*nbf), &one,
                    dW1, nbf, dB1 + (size_t)t0*nn, nbf,
                    &bta, dJ1 + (size_t)v*nn, nbf));
            }
        }
        CUTRY(cudaSetDevice(0));
        for (int v = 0; v < nvec; ++v) {
            for (int t0 = 0; t0 < g_mgn0; t0 += KT) {
                const int nt = std::min(KT, g_mgn0 - t0);
                CBTRY(cublasDgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N,
                    nbf, nbf, nbf, &one,
                    dB + (size_t)t0*nn, nbf, nn,
                    dXin + (size_t)v*nn, nbf, 0LL,
                    &zero, dW, nbf, nn, nt));
                const double bta = (t0 == 0) ? 0.0 : 1.0;
                CBTRY(cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T,
                    nbf, nbf, (int)((long long)nt*nbf), &one,
                    dW, nbf, dB + (size_t)t0*nn, nbf,
                    &bta, dKout + (size_t)v*nn, nbf));
            }
        }
        CUTRY(cudaSetDevice(1)); CUTRY(cudaStreamSynchronize(s1)); CUTRY(cudaSetDevice(0));
        CUTRY(cudaMemcpyPeer(dPeer,0,dJ1,1,(size_t)nX*8));
        CBTRY(cublasDaxpy(hb,(int)nX,&one,dPeer,1,dKout,1));
    } else if (dKout && g_cdf_on) {
        // ---- tiled full-D K from compacted B (guess path): K_v = Σ_P (B_P M_v) B_P^T
        static double* dWt=nullptr; static size_t capWt=0;
        if(!grow(&dWt,&capWt,(size_t)CDF_TILE*nn)) return false;
        for(int v=0; v<nvec; ++v){ double* Kv=dKout+(size_t)v*nn;
          for(int p0=0;p0<naux;p0+=CDF_TILE){ const int nt=std::min(CDF_TILE,naux-p0);
            if(!cdf_tile(p0,nt,nbf)) return false;
            CBTRY(cublasDgemmStridedBatched(hb,CUBLAS_OP_N,CUBLAS_OP_N, nbf,nbf,nbf,&one,
                dBtile,nbf,nn, dXin+(size_t)v*nn,nbf,0LL, &zero,dWt,nbf,nn, nt));   // W_a=B_a·M_v
            const double bta=(p0==0)?0.0:1.0;
            CBTRY(cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_T, nbf,nbf,(int)((long long)nt*nbf),
                &one, dWt,nbf, dBtile,nbf, &bta, Kv,nbf)); }                        // K_v+=ΣW·B^T
        }
    } else if (dKout && g_fp32k) {
        // ---- FP32-K: identical algebra, single precision (KP step-2) --------
        // RI-K is ~95-98% of the DF apply, so this is the real speed lever.
        // Per-vector strided batch only: [K1f] W=B_P·M and [K2f] K=W·B^T as
        // cublasSgemmStridedBatched against the float B mirror dBf; K cast back
        // to double on output. Safe while the SCF density is rough (loose DIIS);
        // the ramp flips back to FP64-K near convergence.
        const float onef = 1.0f, zerof = 0.0f;
        const long long nX = nn*nvec;
        if (!ensure_Bf()) return false;
        if (!ensure_wf(nvec)) return false;
        if (!growf(&dXkf, &capXkf, (size_t)nX)) return false;
        if (!growf(&dKf,  &capKf,  (size_t)nX)) return false;
        // TF32 tensor cores are the REAL A100 K lever (plain FP32 SIMT ties the
        // FP64-TC baseline): OQP_JKM_FP32K_TF32=1 routes the K Sgemms onto TF32
        // TC (~4.5-4.8x K-only vs FP64-K) at ~2e-4 rel error — loose-early-only,
        // exactly what the ramp gates. Default plain FP32 (1e-6 band).
        const bool tf32 = getenv("OQP_JKM_FP32K_TF32") &&
                          atoi(getenv("OQP_JKM_FP32K_TF32")) != 0;
        cublasMath_t prev_mm = CUBLAS_DEFAULT_MATH;
        if (tf32) { cublasGetMathMode(hb, &prev_mm);
                    cublasSetMathMode(hb, CUBLAS_TF32_TENSOR_OP_MATH); }
        cast_d2f<<<(unsigned)((nX+255)/256), 256>>>(dXin, dXkf, nX);
        CUTRY(cudaGetLastError());
        for (int v0 = 0; v0 < nvec; v0 += wf_vchunk) {
            const int nc = std::min(wf_vchunk, nvec - v0);
            for (int vi = 0; vi < nc; ++vi) {           // [K1f] per-vector batch
                CBTRY(cublasSgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N,
                    nbf, nbf, nbf, &onef,
                    dBf, nbf, nn,
                    dXkf + (size_t)(v0+vi)*nn, nbf, 0LL,
                    &zerof, dWf + (size_t)vi*(size_t)naux*nn, nbf, nn, naux));
            }
            CBTRY(cublasSgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_T,
                nbf, nbf, (int)((long long)naux*nbf), &onef,
                dWf, nbf, (long long)naux*nn,
                dBf, nbf, 0LL,
                &zerof, dKf + (size_t)v0*nn, nbf, nn, nc));
        }
        cast_f2d<<<(unsigned)((nX+255)/256), 256>>>(dKf, dKout, nX);
        CUTRY(cudaGetLastError());
        if (tf32) cublasSetMathMode(hb, prev_mm);   // restore FP64 GEMM math
    } else if (dKout) {
        // ---- aux-tiled full-D K (default): identical math, but W holds only
        // KT aux slices at a time, so the workspace is KT*nbf^2 instead of
        // naux*nbf^2 (~1.2 GB vs ~17 GB at (H2O)32). [K2] accumulates over the
        // tiles with beta=1; each chunk still contracts KT*nbf (>~100k) columns,
        // so the GEMMs stay at full cuBLAS efficiency. OQP_JKM_KTILE=0 restores
        // the monolithic workspace (and the opt-in PTRARRAY/VCHUNK folds use it).
        static int s_ktile = getenv("OQP_JKM_KTILE") ? atoi(getenv("OQP_JKM_KTILE")) : 256;
        if (s_ktile > 0 && !getenv("OQP_JKM_PTRARRAY") && !getenv("OQP_JKM_VCHUNK")) {
            const int KT = std::min(s_ktile, naux);
            if (!ensure_w_ktile(KT)) return false;
            for (int v = 0; v < nvec; ++v) {
                for (int t0 = 0; t0 < naux; t0 += KT) {
                    const int nt = std::min(KT, naux - t0);
                    CBTRY(cublasDgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N,
                        nbf, nbf, nbf, &one,
                        dB + (size_t)t0*nn, nbf, nn,        // A_P = B_P, P in tile
                        dXin + (size_t)v*nn, nbf, 0LL,      // B   = M_v (shared)
                        &zero, dW, nbf, nn, nt));           // W_P = B_P M_v
                    const double bta = (t0 == 0) ? 0.0 : 1.0;
                    CBTRY(cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T,
                        nbf, nbf, (int)((long long)nt*nbf), &one,
                        dW, nbf, dB + (size_t)t0*nn, nbf,
                        &bta, dKout + (size_t)v*nn, nbf));  // K_v += W · B_tile^T
                }
            }
            return true;
        }
        if (!ensure_w(nvec)) return false;
        for (int v0 = 0; v0 < nvec; v0 += w_vchunk) {
            const int nc = std::min(w_vchunk, nvec - v0);
            // ---- [K1]  W_{P,v} = B_P · M_v   for all (P, v-in-chunk)
            // Batch index b = vi*naux + P  (P fastest). Per batch, col-major
            // opN/opN nbf^3 GEMM:
            //   A_b = dB + P*nn, ld nbf: col-major view (r,c) at c*nbf+r =
            //         B_P[c,r] = B_P[r,c]  => A_b == B_P   (slice symmetry)
            //   B_b = dXin + v*nn, ld nbf: (k,l) at l*nbf+k = M_v[k,l]  exact
            //   C_b = dW + b*nn, ld nbf: stores W col-major:
            //         dW[b*nn + l*nbf + a] = W_{P,v}[a,l]
            // Element sum: C_b[a,l] = sum_k A_b[a,k] B_b[k,l]
            //   = sum_k B_P[a,k] M_v[k,l] = (B_P M_v)[a,l].          OK
            // C offsets are uniform in b (stride nn) but A repeats with
            // period naux and B with period 1/naux — NOT expressible as one
            // strided-batched call, hence pointer-array cublasDgemmBatched
            // when nc > 1 (the (P,v) fold), strided-batched when nc == 1.
            if (nc == 1) {
                CBTRY(cublasDgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N,
                    nbf, nbf, nbf, &one,
                    dB, nbf, nn,                       // A_P = B_P
                    dXin + (size_t)v0*nn, nbf, 0LL,    // B   = M_{v0} (shared)
                    &zero, dW, nbf, nn, naux));
            } else {
                const int nbatch = naux*nc;
                hA.resize(nbatch); hBp.resize(nbatch); hC.resize(nbatch);
                for (int vi = 0; vi < nc; ++vi) {
                    const double* xv = dXin + (size_t)(v0+vi)*nn;
                    const size_t boff = (size_t)vi*naux;
                    for (int P = 0; P < naux; ++P) {
                        hA [boff+P] = dB + (size_t)P*nn;
                        hBp[boff+P] = xv;
                        hC [boff+P] = dW + (boff+P)*nn;
                    }
                }
                CUTRY(cudaMemcpy((void*)dAarr, hA.data(),
                      (size_t)nbatch*sizeof(double*), cudaMemcpyHostToDevice));
                CUTRY(cudaMemcpy((void*)dBarr, hBp.data(),
                      (size_t)nbatch*sizeof(double*), cudaMemcpyHostToDevice));
                CUTRY(cudaMemcpy((void*)dCarr, hC.data(),
                      (size_t)nbatch*sizeof(double*), cudaMemcpyHostToDevice));
                CBTRY(cublasDgemmBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N,
                    nbf, nbf, nbf, &one,
                    dAarr, nbf, dBarr, nbf, &zero, dCarr, nbf, nbatch));
            }
            // ---- [K2]  K_v = W_v · Bview^T   (super-index s = (P,l))
            // Col-major algebra, strided-batched over vi in the chunk:
            //   A_vi = dW + vi*naux*nn, viewed (nbf x naux*nbf), ld nbf:
            //     A_vi[a, s=(P,l)] at vi*naux*nn + s*nbf + a
            //                       = dW[(vi*naux+P)*nn + l*nbf + a]
            //                       = W_{P,vi}[a,l]      (matches [K1] b-order)
            //   B    = dB, viewed Bview (nbf x naux*nbf), ld nbf, opT,
            //          stride 0 (shared):  Bview[r, s=(P,l)] at s*nbf + r
            //                       = dB[P*nn + l*nbf + r] = B_P[l,r]
            //   C_vi = dKout + (v0+vi)*nn, ld nbf, stride nn.
            // Element sum: C[a,c] = sum_s A[a,s]*Bview[c,s]
            //   = sum_{P,l} (B_P M_v)[a,l] * B_P[l,c]
            //   = sum_{P,k,l} B_P[a,k] M_v[k,l] B_P[l,c] = K(M_v)[a,c].
            // (No symmetry needed in this step; the definition lands exactly.)
            // Stored col-major: dKout[v*nn + c*nbf + a] = K_v[a,c] — precisely
            // the Fortran km(a,b,v). RAW K, no symmetrization, no scaling.
            CBTRY(cublasDgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_T,
                nbf, nbf, (int)((long long)naux*nbf), &one,
                dW, nbf, (long long)naux*nn,
                dB, nbf, 0LL,
                &zero, dKout + (size_t)v0*nn, nbf, nn, nc));
        }
    }
    return true;
}

// ---------------------------------------------------------- occ-K (#145) ----
// K(D) from the OCC factor Cocc (nbf x nocc, col-major) instead of the full D.
//   K(D)_ac = 2 sum_P sum_o (B_P Cocc)[a,o] (B_P Cocc)[c,o]
// Pass dCo = sqrt(2)*Cocc so K = sum_P (B_P dCo)(B_P dCo)^T exactly = K(D).
// Two GEMMs:
//   [K1occ] M_P = B_P * dCo   (nbf x nocc), strided-batched over P -> dMocc
//           stored col-major nbf x (naux*nocc): block P at +P*(nbf*nocc), ld nbf.
//   [K2occ] K = dMocc * dMocc^T   (nbf x nbf), one GEMM, contract dim naux*nocc.
// FP64 by default; OQP_JKM_FP32K=1 (or the SCF ramp via g_fp32k) routes both
// GEMMs through the float twin (optionally TF32 TC), K cast back to double.
// Cost O(naux*nocc*nbf^2): ~nbf/nocc fewer FLOPs than jkm_core's full-D K.
bool scf_K_occ(const double* dCo, int nocc, double* dKout) {
    const int nbf = g_nbf, naux = g_naux;
    const long long nn = (long long)nbf*nbf;
    const long long mlen = (long long)naux*nocc*nbf;   // dMocc element count
    const double one = 1.0, zero = 0.0;
    if (g_mg) {
        // ---- dual-GPU occ-RI-K: each device half-transforms its own aux
        // slice (M = B_g Cocc, K_g = M M^T), then K += peer(K1). Per-cycle
        // exchange: Cocc down (nbf*nocc) + K partial back (nn) -- a few MB.
        if (!grow1(&dCo1,&capCo1,(size_t)nbf*nocc)) return false;
        if (!grow1(&dMocc1,&capMocc1,(size_t)g_mgn1*nocc*nbf)) return false;
        if (!grow1(&dK1,&capK1,(size_t)nn)) return false;
        if (!grow(&dMocc,&capMocc,(size_t)g_mgn0*nocc*nbf)) return false;
        if (!grow(&dPeer,&capPeer,(size_t)nn)) return false;
        CUTRY(cudaMemcpyPeer(dCo1,1,dCo,0,(size_t)nbf*nocc*8));
        CUTRY(cudaSetDevice(1));
        CBTRY(cublasDgemmStridedBatched(hb1, CUBLAS_OP_N, CUBLAS_OP_N,
            nbf, nocc, nbf, &one,
            dB1, nbf, nn, dCo1, nbf, 0LL,
            &zero, dMocc1, nbf, (long long)nbf*nocc, g_mgn1));
        CBTRY(cublasDgemm(hb1, CUBLAS_OP_N, CUBLAS_OP_T,
            nbf, nbf, (int)((long long)g_mgn1*nocc), &one,
            dMocc1, nbf, dMocc1, nbf, &zero, dK1, nbf));
        CUTRY(cudaSetDevice(0));
        CBTRY(cublasDgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N,
            nbf, nocc, nbf, &one,
            dB, nbf, nn, dCo, nbf, 0LL,
            &zero, dMocc, nbf, (long long)nbf*nocc, g_mgn0));
        CBTRY(cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T,
            nbf, nbf, (int)((long long)g_mgn0*nocc), &one,
            dMocc, nbf, dMocc, nbf, &zero, dKout, nbf));
        CUTRY(cudaSetDevice(1)); CUTRY(cudaStreamSynchronize(s1)); CUTRY(cudaSetDevice(0));
        CUTRY(cudaMemcpyPeer(dPeer,0,dK1,1,(size_t)nn*8));
        CBTRY(cublasDaxpy(hb,(int)nn,&one,dPeer,1,dKout,1));
        return true;
    }
    if (g_cdf_on) {
        // ---- tiled occ-RI-K from compacted B (no full dense B on device) ----
        // K = Σ_tiles (B_tile·Cocc)(B_tile·Cocc)^T, accumulated over aux tiles.
        if (!grow(&dMocc, &capMocc, (size_t)CDF_TILE*nocc*nbf)) return false;
        for (int p0=0;p0<naux;p0+=CDF_TILE){ const int nt=std::min(CDF_TILE,naux-p0);
            if(!cdf_tile(p0,nt,nbf)) return false;
            CBTRY(cublasDgemmStridedBatched(hb,CUBLAS_OP_N,CUBLAS_OP_N, nbf,nocc,nbf,&one,
                dBtile,nbf,nn, dCo,nbf,0LL, &zero,dMocc,nbf,(long long)nbf*nocc, nt));  // [K1]
            const double bta=(p0==0)?0.0:1.0;
            CBTRY(cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_T, nbf,nbf,(int)((long long)nt*nocc),
                &one, dMocc,nbf, dMocc,nbf, &bta, dKout,nbf)); }                        // [K2]
        return true;
    }
    (void)mlen;
    if (g_fp32k) {
        const float onef = 1.0f, zerof = 0.0f;
        if (!ensure_Bf()) return false;
        if (!growf(&dMoccf, &capMoccf, (size_t)mlen)) return false;
        if (!growf(&dCof,   &capCof,   (size_t)nbf*nocc)) return false;
        if (!growf(&dKf,    &capKf,    (size_t)nn)) return false;   // float K out
        cast_d2f<<<(unsigned)(((long long)nbf*nocc+255)/256),256>>>(dCo, dCof,
            (long long)nbf*nocc);
        CUTRY(cudaGetLastError());
        const bool tf32 = getenv("OQP_JKM_FP32K_TF32") &&
                          atoi(getenv("OQP_JKM_FP32K_TF32")) != 0;
        cublasMath_t prev_mm = CUBLAS_DEFAULT_MATH;
        if (tf32) { cublasGetMathMode(hb, &prev_mm);
                    cublasSetMathMode(hb, CUBLAS_TF32_TENSOR_OP_MATH); }
        // [K1occ] M_P = B_P * Cocc  (batched, shared Cocc, stride 0)
        CBTRY(cublasSgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N,
            nbf, nocc, nbf, &onef,
            dBf, nbf, nn,                          // A_P = B_P
            dCof, nbf, 0LL,                        // B   = Cocc (shared)
            &zerof, dMoccf, nbf, (long long)nbf*nocc, naux));
        // [K2occ] K = M * M^T  (contract naux*nocc)
        CBTRY(cublasSgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T,
            nbf, nbf, (int)((long long)naux*nocc), &onef,
            dMoccf, nbf, dMoccf, nbf, &zerof, dKf, nbf));
        if (tf32) cublasSetMathMode(hb, prev_mm);
        cast_f2d<<<(unsigned)((nn+255)/256),256>>>(dKf, dKout, nn);
        CUTRY(cudaGetLastError());
        return true;
    }
    if (!grow(&dMocc, &capMocc, (size_t)mlen)) return false;
    // [K1occ] M_P = B_P * Cocc   (strided-batched over P; Cocc shared, stride 0)
    CBTRY(cublasDgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N,
        nbf, nocc, nbf, &one,
        dB, nbf, nn,                               // A_P = B_P
        dCo, nbf, 0LL,                             // B   = Cocc (shared)
        &zero, dMocc, nbf, (long long)nbf*nocc, naux));
    // [K2occ] K = M * M^T   (nbf x nbf, contraction over naux*nocc)
    CBTRY(cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T,
        nbf, nbf, (int)((long long)naux*nocc), &one,
        dMocc, nbf, dMocc, nbf, &zero, dKout, nbf));
    return true;
}

}  // namespace

// =============================================================================
extern "C" {

// ------------------------------------------------------------ v1 SCF seam ---
// d, f: packed lower-triangle (ntri, nfocks), Fortran column-major (each
// slice contiguous). d holds plain elements of the occupation-scaled density.
// nfocks=1 (R):    F = scale_coul*J(D) - 0.5*scale_exch*K(D).
// nfocks=2 (ROHF): F_s = scale_coul*J(D1+D2) - scale_exch*K(D_s) — full x,
//   total-density J (the v1 per-slot/half-x bug reproduced here was fixed
//   2026-06-12, mirroring the CPU dylib fix; ~23 Ha wrong in-loop before).
// info: 0 ok, nonzero => caller falls back to the native path.
// v2 superset of v1: any nfocks (batched through the core), packed B accepted.
void routec_fock_jk(const double* d, double* f, const int* nbf_,
                    const int* nfocks_, const double* scale_exch,
                    const double* scale_coul, int* info) {
    *info = 1;
    if (!load_B()) return;
    const int nbf = *nbf_, nfocks = *nfocks_;
    if (nbf != g_nbf) {
        fprintf(stderr, "[routec-gpu2] fock_jk nbf mismatch %d vs %d\n",
                nbf, g_nbf);
        return;
    }
    if (nfocks < 1) return;
    const size_t ntri = (size_t)nbf*(nbf+1)/2, nn = (size_t)nbf*nbf;
    if (!grow(&ddp, &capP, ntri*nfocks) || !grow(&dfp, &capF, ntri*nfocks) ||
        !grow(&dX, &capX, nn*nfocks)    || !grow(&dJ, &capJ, nn*nfocks) ||
        !grow(&dK, &capK, nn*nfocks)) return;
    if (cudaMemcpy(ddp, d, ntri*nfocks*sizeof(double),
                   cudaMemcpyHostToDevice) != cudaSuccess) return;
    {
        const long long tot = (long long)nn*nfocks;
        unpack_tri_multi<<<(unsigned)((tot+255)/256), 256>>>(ddp, dX, nbf,
                                                             nfocks);
        if (cudaGetLastError() != cudaSuccess) return;
    }
    if (!jkm_core(dX, nfocks, dJ, dK)) return;
    // OpenQP urohf convention (nfocks=2): F_s = c*J(D1+D2) - x*K(D_s) with
    // FULL x (not 0.5x). J is linear, so J(D1+D2) = J0+J1 in place.
    if (nfocks == 2) {
        const double one = 1.0;
        if (cublasDaxpy(hb, (int)nn, &one, dJ+nn, 1, dJ, 1)
            != CUBLAS_STATUS_SUCCESS) return;
        if (cudaMemcpy(dJ+nn, dJ, nn*sizeof(double),
                       cudaMemcpyDeviceToDevice) != cudaSuccess) return;
    }
    const double cx = (nfocks == 1) ? 0.5*(*scale_exch) : (*scale_exch);
    {
        const long long tot = (long long)ntri*nfocks;
        pack_tri_multi<<<(unsigned)((tot+255)/256), 256>>>(dJ, dK, dfp, nbf,
            nfocks, *scale_coul, cx);
        if (cudaGetLastError() != cudaSuccess) return;
    }
    if (cudaMemcpy(f, dfp, ntri*nfocks*sizeof(double),
                   cudaMemcpyDeviceToHost) != cudaSuccess) return;
    *info = 0;
}

// ----------------------------------------------------------- APIV2 (R3b) ---
int routec_jkm_init(const int* nbf) {
    if (!load_B()) return 1;
    if (*nbf != g_nbf) {
        fprintf(stderr, "[routec-gpu2] jkm_init nbf mismatch: %d vs B file %d\n",
                *nbf, g_nbf);
        return 2;
    }
    return 0;
}

void routec_jkm_free(void) {
    // full teardown; a later call re-loads lazily (mirrors CPU v2 semantics).
    // dB is borrowed from df_store (shared with sigma) -- release the pointer
    // but do NOT free the buffer; df_store owns it (routec_b_free).
    if (dB)  { if (!b_borrowed) cudaFree(dB); dB = nullptr; b_borrowed = false; }
    if (dBf) { cudaFree(dBf); dBf = nullptr; }
    if (dWf) { cudaFree(dWf); dWf = nullptr; wf_vchunk = 0; }
    if (dXkf){ cudaFree(dXkf); dXkf = nullptr; capXkf = 0; }
    if (dKf) { cudaFree(dKf); dKf = nullptr; capKf = 0; }
    g_fp32k = false;
    if (dX)  { cudaFree(dX);  dX  = nullptr; capX = 0; }
    if (dJ)  { cudaFree(dJ);  dJ  = nullptr; capJ = 0; }
    if (dK)  { cudaFree(dK);  dK  = nullptr; capK = 0; }
    if (dG)  { cudaFree(dG);  dG  = nullptr; capG = 0; }
    if (dW)  { cudaFree(dW);  dW  = nullptr; w_vchunk = 0; w_elems = 0; }
    if (dMocc) { cudaFree(dMocc); dMocc = nullptr; capMocc = 0; }
    if (dMoccf){ cudaFree(dMoccf); dMoccf = nullptr; capMoccf = 0; }
    if (dCof)  { cudaFree(dCof);  dCof = nullptr; capCof = 0; }
    if (ddp) { cudaFree(ddp); ddp = nullptr; capP = 0; }
    if (dfp) { cudaFree(dfp); dfp = nullptr; capF = 0; }
    if (dU)  { cudaFree(dU);  dU  = nullptr; capU = 0; }
    if (dV)  { cudaFree(dV);  dV  = nullptr; capV = 0; }
    if (dYZ) { cudaFree(dYZ); dYZ = nullptr; capYZ = 0; }
    if (dGlr){ cudaFree(dGlr); dGlr = nullptr; capGlr = 0; }
    if (dAarr) { cudaFree((void*)dAarr); dAarr = nullptr; }
    if (dBarr) { cudaFree((void*)dBarr); dBarr = nullptr; }
    if (dCarr) { cudaFree((void*)dCarr); dCarr = nullptr; }
    cap_ptr = 0;
    if (hb) { cublasDestroy(hb); hb = nullptr; }
    g_naux = 0; g_nbf = 0;
}

// x  : X(nbf,nbf,nvec) Fortran order; X(a,b,v) = M^v_ab (full square, raw,
//      generally NON-symmetric). jm/km: J(M_v)/K(M_v), same layout; either
//      may be NULL to skip. RAW output — no prefactors, no symmetrization.
void routec_jkm_apply(const double* x, const int* nbf_, const int* nvec_,
                      double* jm, double* km, int* info) {
    *info = 1;
    if (!load_B()) return;
    const int nbf = *nbf_, nvec = *nvec_;
    if (nbf != g_nbf) {
        fprintf(stderr, "[routec-gpu2] jkm_apply nbf mismatch: %d vs %d\n",
                nbf, g_nbf);
        return;
    }
    if (nvec <= 0 || (!jm && !km)) { *info = (nvec >= 0) ? 0 : 1; return; }
    const size_t nn = (size_t)nbf*nbf;
    if (!grow(&dX, &capX, nn*nvec)) return;
    if (jm && !grow(&dJ, &capJ, nn*nvec)) return;
    if (km && !grow(&dK, &capK, nn*nvec)) return;
    if (cudaMemcpy(dX, x, nn*nvec*sizeof(double),
                   cudaMemcpyHostToDevice) != cudaSuccess) return;
    if (!jkm_core(dX, nvec, jm ? dJ : nullptr, km ? dK : nullptr)) return;
    if (jm && cudaMemcpy(jm, dJ, nn*nvec*sizeof(double),
                         cudaMemcpyDeviceToHost) != cudaSuccess) return;
    if (km && cudaMemcpy(km, dK, nn*nvec*sizeof(double),
                         cudaMemcpyDeviceToHost) != cudaSuccess) return;
    *info = 0;
}

// --------------------------------------------------------- APIV2.1 (R3d) ---
// Low-rank apply: nout outputs, output o is M_o = sum_r u_r v_r^T with
// rank_o = ranks[o] (>= 0; designed for <= 2 — MRSF slots 1-6). u, v:
// Fortran (nbf, nrt) column-major, nrt = sum(ranks); column r of u pairs
// with column r of v, grouped by output in order; caller folds signs into
// the columns. jm/km: (nbf,nbf,nout) Fortran, RAW, either may be NULL.
//
// Math (B_P symmetric):  K(u v^T)_ab = sum_P (B_P u)_a (B_P v)_b ;
//                        J(u v^T)    : g_P = u^T B_P v, J = sum_P g_P B_P.
// No nbf^3 term anywhere: ~6*naux*nbf^2 FLOPs per rank-1 output vs
// 4*naux*nbf^3 dense.
//
// Layout discipline (derived twice, as everywhere in this file):
//   * factor staging dUbuf: device (nbf x 2*nc) col-major = [U-cols|V-cols]
//     of the current chunk (nc rank terms); column t at t*nbf (U side),
//     (nc+t)*nbf (V side) — vectors, no transpose ambiguity.
//   * [LR1] one GEMM per chunk, ONE pass over B for ALL factor columns:
//     A = dB viewed col-major (nbf x naux*nbf), ld nbf ("Bcols"):
//       Bcols[r, s=(P,l)] at s*nbf + r = dB[P*nn + l*nbf + r] = B_P[l,r].
//     op(A)=T: (naux*nbf x nbf).  C = dYZ (naux*nbf x 2*nc), ld naux*nbf:
//       C[s,c] = sum_j Bcols[j,s] * W[j,c] = sum_j B_P[l,j] w_c[j]
//              = (B_P w_c)_l                                (NO symmetry used)
//     Column c of dYZ viewed col-major (nbf x naux), ld nbf:
//       Yc[l,P] at c*naux*nbf + P*nbf + l                   (matches C[s,c])
//   * [LR2] K_o accumulation, per rank term t of output o:
//     opN/opT GEMM (nbf x nbf x naux), A = Y_t (nbf x naux, ld nbf),
//     B = Z_t (same view, V side), beta = 0 first term else 1:
//       C[a,b] = sum_P Y_t[a,P] * Z_t[b,P] = sum_P (B_P u)_a (B_P v)_b
//              = K(u_t v_t^T)[a,b]
//     stored col-major at o*nn + b*nbf + a — precisely Fortran km(a,b,o).
//     RAW: swapping u<->v transposes K (the lr tripwire pins this).
//   * [LR3] J g-vectors: per (o, t) DGEMV opT on Y_t (nbf x naux, ld nbf)
//     with x = raw factor column v_t (dUbuf V side), beta = 1 into the
//     pre-zeroed column o of dGlr (naux x no, col-major, ld naux):
//       g[P] += sum_l Y_t[l,P] * v_t[l] = (B_P u_t)^T v_t = u_t^T B_P v_t
//     (slice symmetry used once, exactly as in the math line above).
//   * [LR4] Jall = Bflat * G, ONE GEMM straight into the jm slab:
//     A = dB viewed col-major (nn x naux), ld nn: Bflat[z,P] = dB[P*nn+z]
//       = B_P[i,j] with z = i*nbf + j.  C = dJout + o0*nn, ld nn:
//       C[z,o] = sum_P B_P[i,j] g_{P,o} = J_o[i,j].
//     Fortran reads jm(a,b,o) at z = b*nbf + a => J_o[b,a] = J_o[a,b] by
//     J symmetry — orientation-safe, same argument as dense [J2].
// Chunking: outputs are grouped so dYZ (2*nc*naux*nbf doubles) stays under
// OQP_JKM_LR_CAP_MB (default 2048 on the A100); each chunk is one extra
// pass over B. FLOPs/chunk: 4*naux*nbf^2*nc [LR1] + 2*naux*nbf^2 per rank
// term [LR2] + 2*naux*nbf^2 per J output [LR4].
void routec_jkm_apply_lr(const double* u, const double* v, const int* nbf_,
                         const int* nout_, const int* ranks,
                         double* jm, double* km, int* info) {
    *info = 1;
    if (!load_B()) return;
    const int nbf = *nbf_, nout = *nout_;
    if (nbf != g_nbf) {
        fprintf(stderr, "[routec-gpu2] jkm_apply_lr nbf mismatch: %d vs %d\n",
                nbf, g_nbf);
        return;
    }
    if (nout <= 0 || (!jm && !km)) { *info = (nout >= 0) ? 0 : 1; return; }
    long long nrt = 0;
    for (int o = 0; o < nout; ++o) {
        if (ranks[o] < 0) {
            fprintf(stderr, "[routec-gpu2] jkm_apply_lr bad rank %d at %d\n",
                    ranks[o], o);
            return;
        }
        nrt += ranks[o];
    }
    const size_t nn = (size_t)nbf*nbf;
    const size_t row = (size_t)g_naux*nbf;          // one Y/Z column (doubles)
    const double one = 1.0, zero = 0.0;

    // chunk outputs so dYZ (2*nc columns of naux*nbf doubles) fits the cap
    size_t cap_mb = 2048;
    if (const char* e = getenv("OQP_JKM_LR_CAP_MB")) {
        long x = atol(e);
        if (x > 0) cap_mb = (size_t)x;
    }
    size_t max_rt = (cap_mb*1048576) / (2*row*sizeof(double));
    if (max_rt < 2) max_rt = 2;                     // >= one rank-2 output

    if (jm && !grow(&dJ, &capJ, nn*nout)) return;
    if (km && !grow(&dK, &capK, nn*nout)) return;

    int o0 = 0;
    long long r0 = 0;                               // chunk start (output, term)
    while (o0 < nout) {
        int o1 = o0, nc = 0;                        // [o0,o1): outputs in chunk
        while (o1 < nout && (nc == 0 || nc + ranks[o1] <= (int)max_rt)) {
            nc += ranks[o1];
            ++o1;
        }
        const int no = o1 - o0;
        if (nc > 0) {
            // stage [U|V] chunk columns and run [LR1]
            if (!grow(&dU, &capU, 2*(size_t)nc*nbf)) return;
            if (!grow(&dYZ, &capYZ, 2*(size_t)nc*row)) return;
            if (cudaMemcpy(dU, u + (size_t)r0*nbf,
                           (size_t)nc*nbf*sizeof(double),
                           cudaMemcpyHostToDevice) != cudaSuccess) return;
            if (cudaMemcpy(dU + (size_t)nc*nbf, v + (size_t)r0*nbf,
                           (size_t)nc*nbf*sizeof(double),
                           cudaMemcpyHostToDevice) != cudaSuccess) return;
            if (cublasDgemm(hb, CUBLAS_OP_T, CUBLAS_OP_N,
                            (int)row, 2*nc, nbf, &one,
                            dB, nbf,                // A: Bcols (nbf x naux*nbf)
                            dU, nbf,                // B: [U|V]  (nbf x 2*nc)
                            &zero, dYZ, (int)row)   // C: (naux*nbf x 2*nc)
                != CUBLAS_STATUS_SUCCESS) return;
        }
        if (jm) {
            if (!grow(&dGlr, &capGlr, (size_t)g_naux*no)) return;
            if (cudaMemset(dGlr, 0, (size_t)g_naux*no*sizeof(double))
                != cudaSuccess) return;
        }
        int t = 0;                                  // rank-term index in chunk
        for (int o = o0; o < o1; ++o) {
            if (km && ranks[o] == 0 &&
                cudaMemset(dK + (size_t)o*nn, 0, nn*sizeof(double))
                    != cudaSuccess) return;
            for (int q = 0; q < ranks[o]; ++q, ++t) {
                const double* Yt = dYZ + (size_t)t*row;
                const double* Zt = dYZ + (size_t)(nc + t)*row;
                if (km &&   // [LR2] K_o += Y_t Z_t^T over P
                    cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T,
                                nbf, nbf, g_naux, &one, Yt, nbf, Zt, nbf,
                                (q == 0) ? &zero : &one,
                                dK + (size_t)o*nn, nbf)
                        != CUBLAS_STATUS_SUCCESS) return;
                if (jm &&   // [LR3] g_{:,o} += Y_t^T v_t
                    cublasDgemv(hb, CUBLAS_OP_T, nbf, g_naux, &one, Yt, nbf,
                                dU + (size_t)(nc + t)*nbf, 1, &one,
                                dGlr + (size_t)(o - o0)*g_naux, 1)
                        != CUBLAS_STATUS_SUCCESS) return;
            }
        }
        if (jm &&           // [LR4] jm slab = Bflat * G (one GEMM)
            cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_N,
                        (int)nn, no, g_naux, &one, dB, (int)nn,
                        dGlr, g_naux, &zero, dJ + (size_t)o0*nn, (int)nn)
                != CUBLAS_STATUS_SUCCESS) return;
        o0 = o1; r0 += nc;
    }
    if (jm && cudaMemcpy(jm, dJ, nn*nout*sizeof(double),
                         cudaMemcpyDeviceToHost) != cudaSuccess) return;
    if (km && cudaMemcpy(km, dK, nn*nout*sizeof(double),
                         cudaMemcpyDeviceToHost) != cudaSuccess) return;
    *info = 0;
}


// ===========================================================================
// routec_scf_solve — A1 Path-G device-resident DF-RHF SCF solver.
//
// Mirrors sessions/20260617_a1cu/pathG_scf_scaffold.py::scf_loop_pathG
// bit-for-physics. OpenQP owns the guess, Hcore/S build, nuclear repulsion
// and final output; it calls this ONCE and gets back the converged density
// (and MO coeffs/energies). Inside the loop there is ZERO host<->device
// traffic in the J/K hot path: H, S uploaded once; X = S^{-1/2} built once
// (cuSOLVER eigh); B already device-resident (load_B). Per cycle, all on
// device: eigh(X^T F X) -> C -> D = 2 Cocc Cocc^T -> RAW J/K via jkm_core ->
// F = H + sc*J - 0.5*se*K -> E -> CDIIS (FDS-SDF residual on device; the
// tiny (<=9) B-matrix solve on host, 81 doubles, negligible & off the hot
// path). The only D2H per cycle are the scalar E and the small CDIIS dot
// products; the density/MO results come back ONCE at the end.
//
// C interface (Fortran seam mirrors routec_fock_jk gating):
//   h, s   : packed lower-tri (ntri), Fortran col-major, OpenQP row-walk
//            i>=j, t = i*(i+1)/2 + j  (same packing as routec_fock_jk d).
//   nbf    : basis dim (must match the resident B).
//   nocc   : doubly-occupied count (nelec/2).
//   enuc   : nuclear repulsion (added to electronic E here).
//   sc, se : Coulomb / exchange scales (RHF: sc=se=1).
//   conv_e : |dE| tolerance; conv_d : ||FDS-SDF|| tolerance.
//   maxit  : cycle cap.
//   e_out  : converged total energy (scalar).
//   d_out  : packed lower-tri converged density (plain elements, occ-2).
//   c_out  : MO coefficients (nbf*nbf, col-major C(a,i)) — may be NULL.
//   eps_out: MO energies (nbf)                            — may be NULL.
//   ncyc_out: cycles used.  info: 0 ok, nonzero => caller falls back.
// ===========================================================================

// unpack ONE packed lower-tri (OpenQP i>=j, t=i*(i+1)/2+j) into a full
// symmetric square, col-major (b*nbf+a). Symmetric => orientation-safe.
__global__ void scf_unpack_tri(const double* dp, double* X, int nbf) {
    long long z = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    const long long nn = (long long)nbf*nbf;
    if (z >= nn) return;
    const int a = (int)(z % nbf), b = (int)(z / nbf);
    const int hi = a > b ? a : b, lo = a > b ? b : a;
    X[z] = dp[(long long)hi*(hi+1)/2 + lo];
}

// pack a full symmetric square (col-major) back to packed lower-tri.
__global__ void scf_pack_tri(const double* M, double* dp, int nbf) {
    long long u = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    const long long ntri = (long long)nbf*(nbf+1)/2;
    if (u >= ntri) return;
    int i = (int)((sqrt(8.0*(double)u + 1.0) - 1.0)*0.5);
    while ((long long)(i+1)*(i+2)/2 <= u) ++i;
    while ((long long)i*(i+1)/2 > u) --i;
    const int j = (int)(u - (long long)i*(i+1)/2);
    dp[u] = M[(long long)j*nbf + i];   // symmetric => (i,j) == (j,i)
}

// col-major column scaling: out[:,c] = in[:,c] * d[c]  (in = U eigenvectors,
// d = s^{-1/2}; out = U diag(s^{-1/2})).
__global__ void scf_scale_cols(const double* in, const double* d, double* out,
                               int nbf) {
    long long z = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    const long long nn = (long long)nbf*nbf;
    if (z >= nn) return;
    const int c = (int)(z / nbf);
    out[z] = in[z] * d[c];
}

// Generalized Wolfsberg-Helmholz (GWH / extended-Hückel) initial Fock, built on
// device from Hcore (dense square) and S (dense square) alone — no external
// guess, no pyscf. This replaces the bare core guess (F0 = H), which on these
// water clusters wastes ~9 SCF cycles thrashing (the first few CDIIS steps
// amplify the core-guess garbage — energy explodes from -1010 to -633 Ha before
// recovering). GWH lands the first density near the SAD/Hückel basin, matching
// native OpenQP's Hückel guess and gpu4pyscf's minao guess in cycle count.
//   F0[i,i] = H[i,i]
//   F0[i,j] = 0.5 * cx * S[i,j] * (H[i,i] + H[j,j])   (i != j),  cx = 1.75
// dHd, dS : dense symmetric (col-major) Hcore / overlap.  out : dense F0.
__global__ void scf_build_gwh(const double* dHd, const double* dS, double* out,
                              int nbf, double cx) {
    long long z = (long long)blockIdx.x*blockDim.x + threadIdx.x;
    const long long nn = (long long)nbf*nbf;
    if (z >= nn) return;
    const int a = (int)(z % nbf), b = (int)(z / nbf);
    const double haa = dHd[(long long)a*nbf + a];
    const double hbb = dHd[(long long)b*nbf + b];
    out[z] = (a == b) ? haa
                      : 0.5 * cx * dS[z] * (haa + hbb);
}

// add a scalar to the diagonal of a col-major n x n matrix (level shift: F += sigma*I).
__global__ void k_add_diag(double* M, double s, int n) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) M[(long)i*n + i] += s;
}

// Exchange-correlation potential (xc.cu). When OQP_OWNXC_DIR is set,
// routec_scf_solve runs RKS/hybrid DFT: each iteration adds Vxc to the Fock and
// Exc to the energy. Caller passes se = hybrid exact-exchange fraction (0 = pure
// GGA, 0.2 = B3LYP, 0.5 = BHHLYP); the functional (OQP_OWNXC_FUNC) must match.
void routec_vxc(const double* d, double* fxc, const int* nbf,
                const int* nfocks, double* eexc, double* totele, int* info);

void routec_scf_solve(const double* h, const double* s, const int* nbf_,
                      const int* nocc_, const double* enuc_,
                      const double* sc_, const double* se_,
                      const double* conv_e_, const double* conv_d_,
                      const int* maxit_,
                      double* e_out, double* d_out,
                      double* c_out, double* eps_out,
                      int* ncyc_out, int* info) {
    *info = 1;
    if (!load_B()) return;
    const int nbf = *nbf_, nocc = *nocc_;
    if (nbf != g_nbf) {
        fprintf(stderr, "[routec-gpu2] scf_solve nbf mismatch: %d vs %d\n",
                nbf, g_nbf);
        return;
    }
    const double sc = *sc_, se = *se_;
    const double conv_e = *conv_e_, conv_d = *conv_d_;
    const int maxit = *maxit_;
    const size_t nn = (size_t)nbf*nbf, ntri = (size_t)nbf*(nbf+1)/2;
    // DFT (RKS/hybrid): active when the XC grid dir is set. RHF path (unset) is
    // untouched and bit-identical.  Host buffers for the per-iteration Vxc call.
    const bool dft_on = (getenv("OQP_OWNXC_DIR") != nullptr);
    // The XC works in the CARTESIAN AO basis; the SCF here is in whatever basis
    // H/S/B use (spherical is well-conditioned).  If OQP_OWNXC_C2S is set (a
    // {ncart,nbf} matrix c2s with AO_sph_i = sum_a c2s[a,i] AO_cart_a), the DFT
    // block maps D_sph->D_cart = c2s D c2s^T for the XC and Vxc_cart->Vxc_sph =
    // c2s^T Vxc c2s back.  Without it, ncart == nbf (pure-cartesian path).
    int ncart = nbf;
    std::vector<double> hC2S;
    std::vector<double> hPk, hVxc;
    if (dft_on) {
        if (const char* cf = getenv("OQP_OWNXC_C2S")) {
            FILE* f = fopen(cf, "rb");
            if (f) { long long a=0,b=0; size_t r=fread(&a,8,1,f); r=fread(&b,8,1,f); (void)r;
                     ncart=(int)a; hC2S.resize((size_t)ncart*nbf);
                     r=fread(hC2S.data(),8,hC2S.size(),f); (void)r; fclose(f); }
        }
        size_t nct=(size_t)ncart*(ncart+1)/2;
        hPk.resize(nct); hVxc.resize(nct);
        fprintf(stderr, "[routec-gpu2] DFT mode: XC from %s func=%s exx=%.3f (nbf=%d ncart=%d)\n",
                getenv("OQP_OWNXC_DIR"), getenv("OQP_OWNXC_FUNC")?getenv("OQP_OWNXC_FUNC"):"?", se, nbf, ncart); }
    const size_t ncc = (size_t)ncart*ncart, nctri = (size_t)ncart*(ncart+1)/2;
    const double one = 1.0, zero = 0.0, mone = -1.0;
    cudaError_t ce; cublasStatus_t cs; cusolverStatus_t ss;

    cusolverDnHandle_t hs = nullptr;
    if (cusolverDnCreate(&hs) != CUSOLVER_STATUS_SUCCESS) return;

    // device buffers (local; freed at exit). dB/dJ/dK reuse the session pool.
    double *dH=nullptr,*dS=nullptr,*dX=nullptr,*dF=nullptr,*dFp=nullptr,
           *dC=nullptr,*dCp=nullptr,*dD=nullptr,*dCo=nullptr,*dJ1=nullptr,
           *dK1=nullptr,*dT=nullptr,*dT2=nullptr,*dpk=nullptr,*dW=nullptr,
           *deps=nullptr,*dFprev=nullptr;
    int *dInfo=nullptr; double *dwork=nullptr; int lwork=0;
    // CDIIS history (device residual/Fock stacks), tiny B-matrix on host.
    const int mdiis = 8;
    double *dEhist=nullptr, *dFhist=nullptr;
    // device-resident CDIIS (#145, fix #2): gather the nd active error vectors
    // into a contiguous (nn x nd) buffer, form the whole B-matrix with ONE GEMM
    // Bm = Egather^T Egather, copy nd^2 doubles to host ONCE. Replaces the
    // nd^2 synchronous cublasDdot D2H round-trips of the original host-CDIIS.
    double *dEgather=nullptr, *dBmat=nullptr;
    // g4p-matched convergence (#153): orbital-rotation gradient g = 2 Cvir^T F Cocc
    // computed on the PRE-DIIS Fock (mirrors gpu4pyscf's norm_gorb). dFCo = F Cocc
    // (nbf x nocc), dGrad = Cvir^T (F Cocc) (nvir x nocc). Only allocated/used when
    // OQP_SCF_G4PCONV=1; the default (legacy) path is untouched / bit-identical.
    double *dFCo=nullptr, *dGrad=nullptr;
    // level-shift: previous cycle's occupied orthonormal MO coeffs (nbf x nocc),
    // used to build the projector P' for the virtual shift sigma*(I - P').
    double *dCpo_prev=nullptr;
    // DFT c2s bridge, fully on device: dC2S holds host-row-major c2s (ncart x nbf),
    // which cuBLAS (column-major) sees as C2S^T (nbf x ncart). dDc/dVc/dTc are the
    // cartesian-frame density/potential/temp, dpkc the packed cartesian staging.
    double *dC2S=nullptr, *dDc=nullptr, *dVc=nullptr, *dTc=nullptr, *dpkc=nullptr;

    #define ALLOC(p,nelem) do{ if(cudaMalloc((void**)&(p),(nelem)*sizeof(double))!=cudaSuccess){ goto cleanup; } }while(0)
    ALLOC(dH,nn); ALLOC(dS,nn); ALLOC(dX,nn); ALLOC(dF,nn); ALLOC(dFp,nn);
    ALLOC(dC,nn); ALLOC(dCp,nn); ALLOC(dD,nn); ALLOC(dCo,(size_t)nbf*nocc);
    ALLOC(dJ1,nn); ALLOC(dK1,nn); ALLOC(dT,nn); ALLOC(dT2,nn); ALLOC(dpk,ntri);
    ALLOC(dW,nbf); ALLOC(deps,nbf); ALLOC(dFprev,nn);
    ALLOC(dEhist,nn*mdiis); ALLOC(dFhist,nn*mdiis);
    ALLOC(dEgather,nn*mdiis); ALLOC(dBmat,(size_t)mdiis*mdiis);
    // #153 g4p-matched convergence buffers (always allocated; cheap: nbf*nocc each).
    ALLOC(dFCo,(size_t)nbf*nocc);
    ALLOC(dGrad,(size_t)(nbf-nocc)*nocc > 0 ? (size_t)(nbf-nocc)*nocc : 1);
    ALLOC(dCpo_prev,(size_t)nbf*nocc);
    if (dft_on) {
        ALLOC(dDc,ncc); ALLOC(dVc,ncc); ALLOC(dTc,ncc); ALLOC(dpkc,nctri);
        if (!hC2S.empty()) {
            ALLOC(dC2S,(size_t)ncart*nbf);
            if (cudaMemcpy(dC2S,hC2S.data(),(size_t)ncart*nbf*sizeof(double),
                           cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
        }
    }
    if (cudaMalloc((void**)&dInfo, sizeof(int)) != cudaSuccess) goto cleanup;

    // upload packed H, S; unpack to dense symmetric squares.
    if (cudaMemcpy(dpk, h, ntri*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
    scf_unpack_tri<<<(unsigned)((nn+255)/256),256>>>(dpk, dH, nbf);
    if (cudaGetLastError()!=cudaSuccess) goto cleanup;
    if (cudaMemcpy(dpk, s, ntri*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
    scf_unpack_tri<<<(unsigned)((nn+255)/256),256>>>(dpk, dS, nbf);
    if (cudaGetLastError()!=cudaSuccess) goto cleanup;

    // workspace for cusolverDnDsyevd (reuse for all eighs: same nbf).
    if (cusolverDnDsyevd_bufferSize(hs, CUSOLVER_EIG_MODE_VECTOR,
            CUBLAS_FILL_MODE_LOWER, nbf, dX, nbf, dW, &lwork)
            != CUSOLVER_STATUS_SUCCESS) goto cleanup;
    if (cudaMalloc((void**)&dwork, (size_t)lwork*sizeof(double))!=cudaSuccess) goto cleanup;

    // ---- X = S^{-1/2} once: eigh(S) -> S = U diag(s) U^T; X = U s^{-1/2} U^T.
    // copy S into dX (syevd overwrites with eigenvectors), eigenvalues -> dW.
    if (cudaMemcpy(dX, dS, nn*sizeof(double), cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
    if (cusolverDnDsyevd(hs, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
            nbf, dX, nbf, dW, dwork, lwork, dInfo)!=CUSOLVER_STATUS_SUCCESS) goto cleanup;
    {
        // host: s^{-1/2}. dW holds ascending eigenvalues (nbf).
        std::vector<double> sev(nbf);
        if (cudaMemcpy(sev.data(), dW, nbf*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup;
        for (int i=0;i<nbf;++i) sev[i] = (sev[i]>1e-12)? 1.0/sqrt(sev[i]) : 0.0;
        if (cudaMemcpy(dW, sev.data(), nbf*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
    }
    // dX currently = U (col-major, columns are eigenvectors). Form T = U * diag(s^{-1/2})
    // by column scaling, then X = T * U^T.  dX col-major: col i scaled by sev[i].
    scf_scale_cols<<<(unsigned)((nn+255)/256),256>>>(dX, dW, dT, nbf); // dT = U diag
    if (cudaGetLastError()!=cudaSuccess) goto cleanup;
    // X = dT * dX^T  (U diag U^T), symmetric.  store into dX via dT2 then copy.
    cs = cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T, nbf, nbf, nbf, &one,
                     dT, nbf, dX, nbf, &zero, dT2, nbf);
    if (cs!=CUBLAS_STATUS_SUCCESS) goto cleanup;
    if (cudaMemcpy(dX, dT2, nn*sizeof(double), cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;

    // ---- SCF loop ----
    // Initial guess (priority):
    //   1. OQP_SCF_GUESS_D=<file>: a packed lower-tri (ntri doubles) guess
    //      density — e.g. OpenQP's own Hückel guess DM_A. Build F0 = H + J - 0.5K
    //      from it (ONE jkm_core pass). This is the in-engine way to match native
    //      OpenQP's Hückel start (no pyscf), and lands the fewest cycles.
    //   2. else GWH (extended-Hückel) Fock from H,S on device — dependency-free,
    //      far better than the bare core guess (F0=H) which wastes ~9 cycles.
    //   3. OQP_SCF_COREGUESS=1 forces the old core guess (A/B test).
    {
        bool coreguess = false;
        if (const char* e = getenv("OQP_SCF_COREGUESS")) coreguess = (atoi(e) != 0);
        const char* gpath = getenv("OQP_SCF_GUESS_D");
        if (coreguess) {
            if (cudaMemcpy(dF, dH, nn*sizeof(double), cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        } else if (gpath && gpath[0]) {
            FILE* gf = fopen(gpath, "rb");
            if (!gf) { fprintf(stderr,"[routec-gpu2] cannot open guess D %s\n", gpath); goto cleanup; }
            std::vector<double> gpk(ntri);
            size_t grd = fread(gpk.data(), sizeof(double), ntri, gf); fclose(gf);
            if (grd != ntri) { fprintf(stderr,"[routec-gpu2] guess D short read %zu/%zu\n",grd,ntri); goto cleanup; }
            if (cudaMemcpy(dpk, gpk.data(), ntri*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
            scf_unpack_tri<<<(unsigned)((nn+255)/256),256>>>(dpk, dD, nbf);
            if (cudaGetLastError()!=cudaSuccess) goto cleanup;
            // pure GGA/LDA (se=0) needs no exchange: skip K so its 17 GB full-D
            // workspace (naux*nbf^2) never allocates — the whole reason a pure
            // functional fits far larger systems than a hybrid on the same card.
            if (!jkm_core(dD, 1, dJ1, se!=0.0 ? dK1 : nullptr)) goto cleanup;  // J(D0)[,K(D0) if hybrid]
            if (cudaMemcpy(dF, dH, nn*sizeof(double), cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            { double a1=sc;      cs=cublasDaxpy(hb,(int)nn,&a1,dJ1,1,dF,1); if(cs)goto cleanup; }
            if (se!=0.0) { double a2=-0.5*se; cs=cublasDaxpy(hb,(int)nn,&a2,dK1,1,dF,1); if(cs)goto cleanup; }
            fprintf(stderr,"[routec-gpu2] SCF guess: density from %s\n", gpath);
            // the guess's full-D K allocated the naux*nbf^2 workspace dW (~17 GB at
            // (H2O)32) which the per-cycle occ-K never touches -- free it so it does
            // not sit resident for the whole SCF. (If a full-D K is ever needed again,
            // ensure_w simply reallocates.)
            if (dW) { cudaFree(dW); dW = nullptr; w_vchunk = 0; w_elems = 0; }
            // DFT: the guess Fock must ALSO carry Vxc(D0) -- without it the first
            // diagonalization sees a pure-J potential (for se=0 literally H+J),
            // which cost several extra SCF cycles (BLYP (H2O)16: 23 vs 18).
            if (dft_on) {
                if (dC2S) {
                    cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,ncart,nbf,&one,
                                   dD,nbf,dC2S,nbf,&zero,dTc,nbf); if(cs)goto cleanup;
                    cs=cublasDgemm(hb,CUBLAS_OP_T,CUBLAS_OP_N,ncart,ncart,nbf,&one,
                                   dC2S,nbf,dTc,nbf,&zero,dDc,ncart); if(cs)goto cleanup;
                } else {
                    if (cudaMemcpy(dDc,dD,ncc*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
                }
                scf_pack_tri<<<(unsigned)((nctri+255)/256),256>>>(dDc, dpkc, ncart);
                if (cudaGetLastError()!=cudaSuccess) goto cleanup;
                if (cudaMemcpy(hPk.data(), dpkc, nctri*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup;
                double exc0=0.0, tot0=0.0; int xinfo0=1, one_f0=1;
                routec_vxc(hPk.data(), hVxc.data(), &ncart, &one_f0, &exc0, &tot0, &xinfo0);
                if (xinfo0!=0){ fprintf(stderr,"[routec-gpu2] guess routec_vxc failed (info=%d)\n",xinfo0); goto cleanup; }
                if (cudaMemcpy(dpkc, hVxc.data(), nctri*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
                scf_unpack_tri<<<(unsigned)((ncc+255)/256),256>>>(dpkc, dVc, ncart);
                if (cudaGetLastError()!=cudaSuccess) goto cleanup;
                if (dC2S) {
                    cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,ncart,ncart,&one,
                                   dC2S,nbf,dVc,ncart,&zero,dTc,nbf); if(cs)goto cleanup;
                    cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_T,nbf,nbf,ncart,&one,
                                   dTc,nbf,dC2S,nbf,&zero,dT2,nbf); if(cs)goto cleanup;
                    cs=cublasDaxpy(hb,(int)nn,&one,dT2,1,dF,1); if(cs)goto cleanup;
                } else {
                    cs=cublasDaxpy(hb,(int)nn,&one,dVc,1,dF,1); if(cs)goto cleanup;
                }
            }
        } else {
            double gwh_cx = 1.75;
            if (const char* e = getenv("OQP_SCF_GWH_CX")) { double v=atof(e); if(v>0) gwh_cx=v; }
            scf_build_gwh<<<(unsigned)((nn+255)/256),256>>>(dH, dS, dF, nbf, gwh_cx);
            if (cudaGetLastError()!=cudaSuccess) goto cleanup;
        }
    }
    {
    double E=0.0, E_old=0.0; int it; int ndiis=0;
    std::vector<double> Bm; // host CDIIS B-matrix scratch
    // ---- adaptive FP32-K ramp ------------------------------------------------
    // Use FP32-K while the DIIS commutator (rnorm) is loose; switch to FP64-K
    // once it drops below OQP_SCF_FP32K_RAMP (default 0 = OFF; set e.g. 1e-3).
    // The decision for cycle `it` is made from the PREVIOUS cycle's rnorm
    // (the density fed to J/K this cycle). Once we drop below the threshold we
    // latch to FP64 for the rest (no flip-flopping). Cycle 0 has no rnorm yet
    // -> FP32-K (the initial density is the roughest of the run).
    double ramp = 0.0;            // 0 disables the ramp (pure FP64-K)
    if (const char* e = getenv("OQP_SCF_FP32K_RAMP")) ramp = atof(e);
    const bool ramp_on = (ramp > 0.0);
    if (ramp_on && !ensure_Bf()) goto cleanup;
    double prev_rnorm = 1e300;    // forces FP32-K on cycle 0
    bool   latched_fp64 = false;  // once true, stay FP64-K
    int    n_fp32 = 0, n_fp64 = 0;
    // ---- early-cycle DIIS conditioning ---------------------------------------
    // From a one-electron / GWH guess the first few error vectors are huge and a
    // naive CDIIS extrapolation AMPLIFIES them (the energy explodes, e.g. -1718
    // -> -1311 -> -1728 on w24), costing many cycles. Two guards, both bit-safe
    // at convergence (they only act while rnorm is large, vanish near the fixed
    // point so the converged density/energy are untouched):
    //   diis_on : extrapolate only once rnorm < this (linear DIIS model valid).
    //   damp    : while rnorm >= diis_on, damp F toward the previous Fock
    //             F <- (1-damp)*F + damp*F_prev  (Roothaan damping) to settle
    //             the guess instead of letting DIIS diverge.
    // DEFAULTS: diis_on huge (always extrapolate) and damp 0 (off) — i.e. plain
    // CDIIS. MEASURED: with the GWH guess, plain CDIIS already converges cleanly
    // (w16 14 cyc, w24 20 cyc at the 1e-5 gate); the rnorm-gated damping HURT
    // (w16 +1, and w24 DIVERGED when its large-rnorm guess never crossed the
    // gate) — so it is opt-in only, for pathological one-electron-guess cases.
    double diis_on = 1e300, damp = 0.0;
    if (const char* e = getenv("OQP_SCF_DIIS_ON")) diis_on = atof(e);
    if (const char* e = getenv("OQP_SCF_DAMP"))    damp    = atof(e);
    // occ-K (#145): the ~nbf/nocc K-FLOP lever; default ON, OQP_SCF_OCCK=0 disables.
    bool use_occk = true;
    if (const char* e = getenv("OQP_SCF_OCCK")) use_occk = (atoi(e) != 0);
    // device-resident CDIIS B-matrix (#145 fix #2): default ON, OQP_SCF_DEVDIIS=0 disables.
    bool use_devdiis = true;
    if (const char* e = getenv("OQP_SCF_DEVDIIS")) use_devdiis = (atoi(e) != 0);
    // #153 g4p-matched stopping rule. NOW THE DEFAULT (the fairness fix): the gate is
    // EXACTLY gpu4pyscf's:  |dE| < conv_e (its conv_tol)  AND  norm_gorb < conv_grad
    // (its conv_tol_grad, norm_gorb = ||2 Cvir^T F Cocc||, the MO occ-vir orbital-
    // rotation gradient L2 norm on the PRE-DIIS Fock). conv_grad from OQP_SCF_CONVG
    // (default 1e-5, matching the g4p harness; pyscf's own default is sqrt(conv_tol)).
    // We also drop the legacy it>3 floor here to mirror g4p (which checks from cycle 1);
    // the water clusters never satisfy the gate that early anyway.
    // OQP_SCF_G4PCONV=0 reverts to the legacy gate (|dE|<conv_e AND rnorm<conv_d on the
    // raw-AO ||FDS-SDF||) for A/B testing — kept for reproducibility, not the default.
    // (#153b: cycle count GROWING with size [14->18 GWH] vs g4p's flat-10 was NOT the
    //  stop rule — it was the GUESS. GWH/extended-Hückel starts ~36 Ha above the basin
    //  and thrashes the first ~3 cycles; a minao/SAD density guess fed via
    //  OQP_SCF_GUESS_D [an offline pyscf-free input artifact, like H/S/B] starts ~1.7 Ha
    //  above and lands at 10/10/11 cycles, matching g4p, with bit-identical converged E.
    //  The pruned-orthonormal-CDIIS lever [diagnosis RANK 1] is a NO-OP on cc-pVDZ water:
    //  min overlap eig ~0.011, nothing below the 1e-6 prune, and WITHOUT pruning the
    //  canonical and symmetric-Löwdin error vectors give the IDENTICAL CDIIS B-matrix
    //  [both Tr(g_i S^-1 g_j S^-1)]. So the guess is the whole lever, not the metric.)
    bool g4pconv = true;
    if (const char* e = getenv("OQP_SCF_G4PCONV")) g4pconv = (atoi(e) != 0);
    double conv_grad = 1e-5;
    if (const char* e = getenv("OQP_SCF_CONVG")) { double v=atof(e); if(v>0) conv_grad=v; }
    // CDIIS restart (limit-cycle guard). Large systems can drive the SCF into a
    // DIIS limit cycle: the residual falls to a good point, then the extrapolation
    // overshoots (ill-conditioned B-matrix -> huge coefficients) and rnorm jumps
    // back up and stalls (e.g. (H2O)32 BLYP: rnorm 0.59 @ it8 -> 2.0 -> stuck).
    // When rnorm jumps by > restart_fac vs the previous cycle while not yet
    // converged, flush the DIIS ring so the next step is a clean, un-extrapolated
    // full Fock (the standard recovery), then rebuild the subspace. No-op on a
    // smooth run (rnorm keeps falling -> never fires). Default on; OQP_SCF_NORESTART=1
    // reverts to plain CDIIS for A/B.
    // DIIS restart is OPT-IN (OQP_SCF_RESTART=1): default off, because the 1.5x
    // jump trigger also fires on the benign rnorm wiggles of an easy run (it cost
    // (H2O)16 ~3 extra cycles) and it does NOT actually converge the hard near-
    // degenerate case anyway (level shift does). Kept for A/B / experimentation.
    bool do_restart = getenv("OQP_SCF_RESTART") && atoi(getenv("OQP_SCF_RESTART"));
    double restart_fac = getenv("OQP_SCF_RESTART_FAC") ? atof(getenv("OQP_SCF_RESTART_FAC")) : 1.5;
    double rnorm_prev_cyc = 1e300;
    int n_restart = 0;
    // Level shift (virtual-orbital shift sigma; cures the occ-virt oscillation that
    // a plain CDIIS limit-cycles on for large systems). F' <- F' + sigma*(I - P'_prev),
    // P'_prev = the previous cycle's occupied orthonormal projector, so occupied
    // orbitals are unshifted and virtuals rise by sigma -> the density update is
    // damped without changing the converged solution (the shift term vanishes on the
    // stationary occupied space). Tapered off once rnorm < lshift_off so the final
    // MOs/eigenvalues are unshifted. Opt-in via OQP_SCF_LSHIFT=sigma (Ha); when on,
    // it supersedes the flush-restart. Default off (bit-identical to before).
    double lshift = getenv("OQP_SCF_LSHIFT") ? atof(getenv("OQP_SCF_LSHIFT")) : 0.0;
    double lshift_off = getenv("OQP_SCF_LSHIFT_OFF") ? atof(getenv("OQP_SCF_LSHIFT_OFF")) : 1e-8;
    double lshift_rref = getenv("OQP_SCF_LSHIFT_RREF") ? atof(getenv("OQP_SCF_LSHIFT_RREF")) : 1.0;
    for (it=0; it<maxit; ++it) {
        // decide K precision for THIS cycle from the previous rnorm
        if (ramp_on) {
            if (!latched_fp64 && prev_rnorm < ramp) latched_fp64 = true;
            g_fp32k = !latched_fp64;
        }
        // 3. Fp = X^T F X ; eigh(Fp) -> Cp ; C = X Cp
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dX,nbf,dF,nbf,&zero,dT,nbf); if(cs)goto cleanup;
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dT,nbf,dX,nbf,&zero,dFp,nbf); if(cs)goto cleanup;
        // level shift on the orthonormal Fock: Fp += sigma*(I - P'_prev). Applied
        // while still far from convergence (rnorm > lshift_off) and only from it>0
        // (need the previous cycle's occupied projector). Occupied space unshifted
        // -> converged density/energy unchanged; virtuals raised -> occ-virt mixing
        // per step is damped, breaking the CDIIS limit cycle.
        if (lshift > 0.0 && it > 0 && rnorm_prev_cyc > lshift_off) {
            // DYNAMIC shift: sigma scales with the residual, sig = lshift*min(1,rnorm/rref).
            // Far out (rnorm>=rref) the full shift breaks the oscillation; as the SCF
            // converges sig -> 0, so the residual floor it induces (~sig*dP ~ lshift*rnorm^2)
            // vanishes faster than rnorm -> reaches the tight gate. A CONSTANT shift instead
            // floors rnorm (shift on) or lets the oscillation return (shift off).
            double sig = lshift * fmin(1.0, rnorm_prev_cyc / lshift_rref);
            int TB=128, nb1=(nbf+TB-1)/TB;
            k_add_diag<<<nb1,TB>>>(dFp, sig, nbf);             // Fp += sig*I
            if (cudaGetLastError()!=cudaSuccess) goto cleanup;
            const double msig = -sig;                          // Fp -= sig*P'_prev
            cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_T,nbf,nbf,nocc,&msig,
                           dCpo_prev,nbf,dCpo_prev,nbf,&one,dFp,nbf); if(cs)goto cleanup;
        }
        if (cudaMemcpy(dCp,dFp,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        ss=cusolverDnDsyevd(hs,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,nbf,dCp,nbf,deps,dwork,lwork,dInfo);
        if(ss!=CUSOLVER_STATUS_SUCCESS) goto cleanup;
        // stash this cycle's occupied orthonormal MOs (first nocc cols, col-major
        // contiguous) for next cycle's level-shift projector P'.
        if (lshift > 0.0 &&
            cudaMemcpy(dCpo_prev,dCp,(size_t)nbf*nocc*sizeof(double),
                       cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        // C = X Cp
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dX,nbf,dCp,nbf,&zero,dC,nbf); if(cs)goto cleanup;
        // 4. Cocc = C[:, :nocc] ; D = 2 Cocc Cocc^T
        if (cudaMemcpy(dCo,dC,(size_t)nbf*nocc*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        {
            const double two=2.0;
            cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_T,nbf,nbf,nocc,&two,dCo,nbf,dCo,nbf,&zero,dD,nbf); if(cs)goto cleanup;
        }
        // 5. J,K.  occ-K (#145, default ON): K(D) from the OCC factor sqrt(2)*Cocc
        //    (cost O(naux*nocc*nbf^2), ~nbf/nocc fewer FLOPs) + J(D) from jkm_core
        //    (J only, O(naux*nbf^2)).  OQP_SCF_OCCK=0 falls back to the full-D
        //    jkm_core K (A/B test).  D is symmetric.
        if (se == 0.0) {
            if (!jkm_core(dD, 1, dJ1, nullptr)) goto cleanup;   // pure GGA/LDA: J only, no K (se=0)
        } else if (use_occk) {
            if (!jkm_core(dD, 1, dJ1, nullptr)) goto cleanup;   // J(D) only
            // scale dCo -> sqrt(2)*Cocc so scf_K_occ yields K(2 Cocc Cocc^T)=K(D)
            { const double rt2=1.4142135623730951;
              cs=cublasDscal(hb,(int)((long long)nbf*nocc),&rt2,dCo,1); if(cs)goto cleanup; }
            if (!scf_K_occ(dCo, nocc, dK1)) goto cleanup;       // K(D) via occ path
        } else {
            if (!jkm_core(dD, 1, dJ1, dK1)) goto cleanup;       // full-D J & K
        }
        // 1. F = H + sc*J - 0.5*se*K   (RHF: -0.5*K; pure DFT: se=0 -> no K term)
        if (cudaMemcpy(dF,dH,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        { double a1=sc; cs=cublasDaxpy(hb,(int)nn,&a1,dJ1,1,dF,1); if(cs)goto cleanup; }
        if (se != 0.0) { double a2=-0.5*se; cs=cublasDaxpy(hb,(int)nn,&a2,dK1,1,dF,1); if(cs)goto cleanup; }
        // 6. E = 0.5 * sum(D*(H+F)) + Enuc   (device dot)
        if (cudaMemcpy(dT,dH,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        cs=cublasDaxpy(hb,(int)nn,&one,dF,1,dT,1); if(cs)goto cleanup; // T = H+F
        { double ehf; cs=cublasDdot(hb,(int)nn,dD,1,dT,1,&ehf); if(cs)goto cleanup;
          E = 0.5*ehf + *enuc_; }
        // DFT: F += Vxc, E += Exc.  Energy above used the pre-XC Fock (H + J +
        // se-scaled K), so adding the functional Exc gives the DFT total; adding
        // Vxc to dF here feeds the next diagonalization and the DIIS residual.
        if (dft_on) {
            // All matrix work on device; only the small PACKED cartesian density
            // and Vxc cross the PCIe bus (routec_vxc's host interface).
            // dC2S buffer, viewed column-major by cuBLAS, is A = C2S^T (nbf x ncart).
            if (dC2S) {
                // Tc = D_sph * A          (nbf x ncart)  = D * C2S^T
                cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,ncart,nbf,&one,
                               dD,nbf,dC2S,nbf,&zero,dTc,nbf); if(cs)goto cleanup;
                // Dc = A^T * Tc           (ncart x ncart) = C2S * D * C2S^T
                cs=cublasDgemm(hb,CUBLAS_OP_T,CUBLAS_OP_N,ncart,ncart,nbf,&one,
                               dC2S,nbf,dTc,nbf,&zero,dDc,ncart); if(cs)goto cleanup;
            } else {
                if (cudaMemcpy(dDc,dD,ncc*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            }
            scf_pack_tri<<<(unsigned)((nctri+255)/256),256>>>(dDc, dpkc, ncart);
            if (cudaGetLastError()!=cudaSuccess) goto cleanup;
            if (cudaMemcpy(hPk.data(), dpkc, nctri*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup;
            double exc=0.0, totele=0.0; int xinfo=1, one_f=1;
            routec_vxc(hPk.data(), hVxc.data(), &ncart, &one_f, &exc, &totele, &xinfo);
            if (xinfo!=0){ fprintf(stderr,"[routec-gpu2] routec_vxc failed (info=%d)\n",xinfo); goto cleanup; }
            if (cudaMemcpy(dpkc, hVxc.data(), nctri*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
            scf_unpack_tri<<<(unsigned)((ncc+255)/256),256>>>(dpkc, dVc, ncart);
            if (cudaGetLastError()!=cudaSuccess) goto cleanup;
            if (dC2S) {
                // Tc = A * Vc             (nbf x ncart)  = C2S^T * Vxc_cart
                cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,ncart,ncart,&one,
                               dC2S,nbf,dVc,ncart,&zero,dTc,nbf); if(cs)goto cleanup;
                // Vs = Tc * A^T           (nbf x nbf)    = C2S^T * Vxc_cart * C2S
                cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_T,nbf,nbf,ncart,&one,
                               dTc,nbf,dC2S,nbf,&zero,dT2,nbf); if(cs)goto cleanup;
                cs=cublasDaxpy(hb,(int)nn,&one,dT2,1,dF,1); if(cs)goto cleanup;  // F += Vxc_sph
            } else {
                cs=cublasDaxpy(hb,(int)nn,&one,dVc,1,dF,1); if(cs)goto cleanup;
            }
            E += exc;
        }
        // 2. CDIIS residual.  Raw AO commutator err = F D S - S D F:
        //    dT = D S ; err = F(DS) - (F(DS))^T via two GEMMs + a geam.
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dD,nbf,dS,nbf,&zero,dT,nbf); if(cs)goto cleanup; // T=D S
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dF,nbf,dT,nbf,&zero,dT2,nbf); if(cs)goto cleanup; // T2=F D S
        cs=cublasDgeam(hb,CUBLAS_OP_N,CUBLAS_OP_T,nbf,nbf,&one,dT2,nbf,&mone,dT2,nbf,dT,nbf); if(cs)goto cleanup; // dT=err_AO
        // Pulay CDIIS canonically uses the error in the ORTHONORMAL (Löwdin)
        //   basis e = X^T(FDS-SDF)X, X = S^{-1/2}. MEASURED on these clusters it
        //   is WORSE by ~2 cycles (24->20 vs 24->18 with the raw-AO commutator):
        //   the cc-pVDZ water cluster S has near-linear-dependent modes, and
        //   s^{-1/2} amplifies their noise in the error, hurting the DIIS metric.
        //   So we DEFAULT to the raw-AO commutator (better conditioned here) and
        //   keep the orthonormalized path behind OQP_SCF_ORTH_ERR=1 for systems
        //   where S is well conditioned.
        bool raw_err = true;
        if (const char* e = getenv("OQP_SCF_ORTH_ERR")) raw_err = (atoi(e) == 0);
        if (const char* e = getenv("OQP_SCF_RAW_ERR")) raw_err = (atoi(e) != 0);
        if (!raw_err) {
            cs=cublasDgemm(hb,CUBLAS_OP_T,CUBLAS_OP_N,nbf,nbf,nbf,&one,dX,nbf,dT,nbf,&zero,dT2,nbf); if(cs)goto cleanup; // T2 = X^T err
            cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dT2,nbf,dX,nbf,&zero,dT,nbf); if(cs)goto cleanup;  // dT = (X^T err) X
        }
        double rnorm; cs=cublasDnrm2(hb,(int)nn,dT,1,&rnorm); if(cs)goto cleanup;
        // #153 g4p-matched convergence: orbital-rotation gradient norm
        //   norm_gorb = || 2 * Cvir^T F Cocc ||_2   (exactly gpu4pyscf get_grad),
        // computed on the PRE-DIIS Fock dF and the CURRENT MO coeffs dC (col-major;
        // Cocc = first nocc cols, Cvir = cols [nocc..nbf)). Only used when
        // OQP_SCF_G4PCONV=1; computed unconditionally (cheap) so the trace can show
        // both metrics. dT/dT2 are free here (commutator done); use dFCo/dGrad.
        double norm_gorb = -1.0;
        {
            const int nvir = nbf - nocc;
            if (nvir > 0) {
                // dFCo = F * Cocc   (nbf x nocc)
                cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nocc,nbf,&one,
                               dF,nbf,dC,nbf,&zero,dFCo,nbf); if(cs)goto cleanup;
                // dGrad = Cvir^T * dFCo  (nvir x nocc), Cvir = dC + nocc*nbf
                cs=cublasDgemm(hb,CUBLAS_OP_T,CUBLAS_OP_N,nvir,nocc,nbf,&one,
                               dC+(size_t)nocc*nbf,nbf,dFCo,nbf,&zero,dGrad,nvir); if(cs)goto cleanup;
                double gn; cs=cublasDnrm2(hb,nvir*nocc,dGrad,1,&gn); if(cs)goto cleanup;
                norm_gorb = 2.0*gn;
            } else {
                norm_gorb = 0.0;
            }
        }
        if (ramp_on) {
            (g_fp32k ? n_fp32 : n_fp64)++;
            prev_rnorm = rnorm;
        }
        static bool s_trace = getenv("OQP_SCF_TRACE") && atoi(getenv("OQP_SCF_TRACE"));
        if (ramp_on || s_trace)
            fprintf(stderr, "[routec-gpu2][scf] it=%2d  K=%s  rnorm=%.3e  "
                    "|g|=%.3e  dE=%.3e  E=%.10f\n", it, g_fp32k ? "FP32" : "FP64",
                    rnorm, norm_gorb, E-E_old, E);
        // CDIIS overshoot -> restart: flush the ring so the push below starts a
        // fresh subspace (nd becomes 1 -> the extrapolation is skipped -> dF stays
        // the plain full Fock: a clean Roothaan recovery step).
        if (do_restart && lshift == 0.0 && ndiis > 1 && rnorm > restart_fac * rnorm_prev_cyc
                && rnorm > conv_grad) {
            ndiis = 0; ++n_restart;
            if (s_trace) fprintf(stderr, "[routec-gpu2][scf]   DIIS restart @ it=%d"
                    " (rnorm %.3e > %.1fx prev %.3e)\n", it, rnorm, restart_fac,
                    rnorm_prev_cyc);
        }
        rnorm_prev_cyc = rnorm;
        // push into history ring.  Optionally drop the cycle-0 error (the GWH
        // guess error is large and can pollute the DIIS B-matrix) via
        // OQP_SCF_DIIS_START=1 (start DIIS history at cycle 1).
        static int s_diis_start = getenv("OQP_SCF_DIIS_START") ? atoi(getenv("OQP_SCF_DIIS_START")) : 0;
        if (it >= s_diis_start) {
            int slot = ndiis % mdiis;
            if (cudaMemcpy(dEhist+(size_t)slot*nn,dT,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            if (cudaMemcpy(dFhist+(size_t)slot*nn,dF,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            ndiis++;
        }
        int nd = ndiis<mdiis ? ndiis : mdiis;
        const bool do_diis = (rnorm < diis_on);
        if (nd>1 && do_diis) {
            // B-matrix: Bm[i,j] = <e_i, e_j>.  map ring slots (oldest..newest) -> 0..nd-1
            Bm.assign((size_t)(nd+1)*(nd+1), -1.0); Bm[(size_t)(nd+1)*(nd+1)-1]=0.0;
            int start = ndiis<=mdiis ? 0 : ndiis%mdiis;
            if (use_devdiis) {
                // ---- device-resident B-matrix (#145 fix #2) ------------------
                // gather the nd active error vectors (ring -> contiguous), then
                // ONE GEMM Bmat = Egather^T Egather, then ONE D2H of nd^2 doubles.
                for (int i=0;i<nd;++i){
                    int si=(start+i)%mdiis;
                    if (cudaMemcpy(dEgather+(size_t)i*nn, dEhist+(size_t)si*nn,
                        nn*sizeof(double), cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
                }
                cs=cublasDgemm(hb,CUBLAS_OP_T,CUBLAS_OP_N,nd,nd,(int)nn,&one,
                               dEgather,(int)nn,dEgather,(int)nn,&zero,dBmat,nd); if(cs)goto cleanup;
                std::vector<double> bb((size_t)nd*nd);
                if (cudaMemcpy(bb.data(), dBmat, (size_t)nd*nd*sizeof(double),
                    cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup;   // ONE D2H
                for (int i=0;i<nd;++i) for(int j=0;j<nd;++j)
                    Bm[(size_t)i*(nd+1)+j]=bb[(size_t)j*nd+i];   // dBmat col-major -> Bm
            } else {
                // ---- original host-CDIIS (nd^2 synchronous dots), A/B fallback
                for (int i=0;i<nd;++i){
                    int si=(start+i)%mdiis;
                    for(int j=0;j<nd;++j){
                        int sj=(start+j)%mdiis; double dij;
                        cs=cublasDdot(hb,(int)nn,dEhist+(size_t)si*nn,1,dEhist+(size_t)sj*nn,1,&dij); if(cs)goto cleanup;
                        Bm[(size_t)i*(nd+1)+j]=dij;
                    }
                }
            }
            // solve Bm c = rhs (rhs = [0..0,-1]) via Gaussian elim on host.
            std::vector<double> A=Bm; std::vector<double> rhs(nd+1,0.0); rhs[nd]=-1.0;
            int N=nd+1; bool ok=true;
            for(int col=0;col<N&&ok;++col){
                int piv=col; double best=fabs(A[(size_t)col*N+col]);
                for(int r=col+1;r<N;++r){double v=fabs(A[(size_t)r*N+col]); if(v>best){best=v;piv=r;}}
                if(best<1e-300){ok=false;break;}
                if(piv!=col){ for(int k=0;k<N;++k) std::swap(A[(size_t)col*N+k],A[(size_t)piv*N+k]); std::swap(rhs[col],rhs[piv]); }
                double d=A[(size_t)col*N+col];
                for(int r=0;r<N;++r){ if(r==col)continue; double f=A[(size_t)r*N+col]/d;
                    for(int k=0;k<N;++k) A[(size_t)r*N+k]-=f*A[(size_t)col*N+k]; rhs[r]-=f*rhs[col]; }
            }
            if(ok){
                for(int i=0;i<N;++i) rhs[i]/=A[(size_t)i*N+i];
                // extrapolate F = sum_i c_i F_i
                if(cudaMemset(dF,0,nn*sizeof(double))!=cudaSuccess) goto cleanup;
                for(int i=0;i<nd;++i){ int si=(start+i)%mdiis; double ci=rhs[i];
                    cs=cublasDaxpy(hb,(int)nn,&ci,dFhist+(size_t)si*nn,1,dF,1); if(cs)goto cleanup; }
            }
        } else if (!do_diis && it>0 && damp>0.0) {
            // early, large-error cycle: damp toward the previous Fock instead of
            // a (divergent) DIIS extrapolation.  F <- (1-damp)*F + damp*F_prev
            double a=1.0-damp, b=damp;
            cs=cublasDscal(hb,(int)nn,&a,dF,1); if(cs)goto cleanup;
            cs=cublasDaxpy(hb,(int)nn,&b,dFprev,1,dF,1); if(cs)goto cleanup;
        }
        // remember the (conditioned) Fock used this cycle for next-cycle damping
        if (cudaMemcpy(dFprev,dF,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        if (g4pconv) {
            // #153 g4p-EXACT gate: |dE| < conv_e AND ||2 Cvir^T F Cocc|| < conv_grad.
            // dE at cycle 0 is meaningless (E_old=0) so require it>0, mirroring g4p
            // (which compares e_tot to last_hf_e starting cycle 1).
            if (it>0 && fabs(E-E_old)<conv_e && norm_gorb<conv_grad) { it++; break; }
        } else {
            if (it>3 && fabs(E-E_old)<conv_e && rnorm<conv_d) { it++; break; }
        }
        E_old=E;
    }
    *ncyc_out = it;
    *e_out = E;
    if (ramp_on) {
        fprintf(stderr, "[routec-gpu2][ramp] DONE: %d cycles FP32-K, %d cycles "
                "FP64-K (threshold %.1e)\n", n_fp32, n_fp64, ramp);
        g_fp32k = false;   // reset ramp state (env default re-read on next load)
    }
    }
    // ---- return converged D (packed), and optional MO C / eps, ONCE ----
    scf_pack_tri<<<(unsigned)((ntri+255)/256),256>>>(dD, dpk, nbf);
    if (cudaGetLastError()!=cudaSuccess) goto cleanup;
    if (cudaMemcpy(d_out, dpk, ntri*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup;
    if (c_out)   { if (cudaMemcpy(c_out,   dC,   nn*sizeof(double),  cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup; }
    if (eps_out) { if (cudaMemcpy(eps_out, deps, nbf*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup; }
    *info = 0;

cleanup:
    #undef ALLOC
    if(dH)cudaFree(dH); if(dS)cudaFree(dS); if(dX)cudaFree(dX); if(dF)cudaFree(dF);
    if(dFp)cudaFree(dFp); if(dC)cudaFree(dC); if(dCp)cudaFree(dCp); if(dD)cudaFree(dD);
    if(dCo)cudaFree(dCo); if(dJ1)cudaFree(dJ1); if(dK1)cudaFree(dK1); if(dT)cudaFree(dT);
    if(dT2)cudaFree(dT2); if(dpk)cudaFree(dpk); if(dW)cudaFree(dW); if(deps)cudaFree(deps);
    if(dFprev)cudaFree(dFprev);
    if(dEhist)cudaFree(dEhist); if(dFhist)cudaFree(dFhist);
    if(dEgather)cudaFree(dEgather); if(dBmat)cudaFree(dBmat);
    if(dFCo)cudaFree(dFCo); if(dGrad)cudaFree(dGrad); if(dCpo_prev)cudaFree(dCpo_prev);
    if(dC2S)cudaFree(dC2S); if(dDc)cudaFree(dDc); if(dVc)cudaFree(dVc);
    if(dTc)cudaFree(dTc); if(dpkc)cudaFree(dpkc);
    if(dInfo)cudaFree(dInfo); if(dwork)cudaFree(dwork);
    if(hs)cusolverDnDestroy(hs);
    (void)ce;
}

// =============================================================================
// routec_scf_solve_uhf — open-shell UHF-HF device-resident DF-SCF solver.
// Mirrors routec_scf_solve (RHF) but carries two spin channels:
//   Da = Ca_occ Ca_occ^T  (nalpha cols, NO factor 2)
//   Db = Cb_occ Cb_occ^T  (nbeta  cols)
//   Jtot = J(Da) + J(Db)
//   Fa = H + sc*Jtot - se*K(Da)
//   Fb = H + sc*Jtot - se*K(Db)
//   E  = 0.5*( Tr[Da (H+Fa)] + Tr[Db (H+Fb)] ) + Enuc
// J(Da),J(Db),K(Da),K(Db) come from ONE batched jkm_core call (nvec=2 over the
// stacked [Da,Db]) — reuses the validated RHF J/K machinery verbatim.
// CDIIS: combined per-spin error e = [Fa Da S - S Da Fa ; Fb Db S - S Db Fb]
// (2*nn-length residual) and combined Fock history [Fa;Fb]; standard UHF DIIS.
// Convergence gate (mirrors RHF): |dE| < conv_e AND ||e_combined|| < conv_d.
// Returns packed Da,Db (+ optional MO Ca,Cb, eps_a,eps_b) to host ONCE.
// =============================================================================
void routec_scf_solve_uhf(const double* h, const double* s, const int* nbf_,
                          const int* nalpha_, const int* nbeta_,
                          const double* enuc_,
                          const double* sc_, const double* se_,
                          const double* conv_e_, const double* conv_d_,
                          const int* maxit_,
                          const double* da0, const double* db0,
                          double* e_out, double* da_out, double* db_out,
                          double* ca_out, double* cb_out,
                          double* epsa_out, double* epsb_out,
                          int* ncyc_out, int* info) {
    *info = 1;
    if (!load_B()) return;
    const int nbf = *nbf_, na = *nalpha_, nb = *nbeta_;
    if (nbf != g_nbf) {
        fprintf(stderr, "[routec-gpu2] scf_solve_uhf nbf mismatch: %d vs %d\n",
                nbf, g_nbf);
        return;
    }
    const double sc = *sc_, se = *se_;
    const double conv_e = *conv_e_, conv_d = *conv_d_;
    const int maxit = *maxit_;
    const int nocc_max = na > nb ? na : nb;
    const size_t nn = (size_t)nbf*nbf, ntri = (size_t)nbf*(nbf+1)/2;
    const double one = 1.0, zero = 0.0, mone = -1.0;
    cublasStatus_t cs; cusolverStatus_t ss;

    cusolverDnHandle_t hs = nullptr;
    if (cusolverDnCreate(&hs) != CUSOLVER_STATUS_SUCCESS) return;

    // device buffers. dDab is the stacked [Da|Db] (2 nn-slices) fed to jkm_core;
    // dJab / dKab receive J(Da),J(Db) / K(Da),K(Db).
    double *dH=nullptr,*dS=nullptr,*dX=nullptr,
           *dFa=nullptr,*dFb=nullptr,*dFp=nullptr,
           *dCa=nullptr,*dCb=nullptr,*dCp=nullptr,
           *dDab=nullptr,*dCo=nullptr,
           *dJab=nullptr,*dKab=nullptr,*dJtot=nullptr,
           *dT=nullptr,*dT2=nullptr,*dpk=nullptr,*dW=nullptr,
           *depsa=nullptr,*depsb=nullptr;
    int *dInfo=nullptr; double *dwork=nullptr; int lwork=0;
    const int mdiis = 8;
    // combined error/Fock history: each entry is 2*nn (alpha then beta).
    double *dEhist=nullptr, *dFhist=nullptr;

    #define ALLOC(p,nelem) do{ if(cudaMalloc((void**)&(p),(nelem)*sizeof(double))!=cudaSuccess){ goto cleanup; } }while(0)
    ALLOC(dH,nn); ALLOC(dS,nn); ALLOC(dX,nn);
    ALLOC(dFa,nn); ALLOC(dFb,nn); ALLOC(dFp,nn);
    ALLOC(dCa,nn); ALLOC(dCb,nn); ALLOC(dCp,nn);
    ALLOC(dDab,2*nn); ALLOC(dCo,(size_t)nbf*nocc_max);
    ALLOC(dJab,2*nn); ALLOC(dKab,2*nn); ALLOC(dJtot,nn);
    ALLOC(dT,nn); ALLOC(dT2,nn); ALLOC(dpk,ntri);
    ALLOC(dW,nbf); ALLOC(depsa,nbf); ALLOC(depsb,nbf);
    ALLOC(dEhist,2*nn*mdiis); ALLOC(dFhist,2*nn*mdiis);
    if (cudaMalloc((void**)&dInfo, sizeof(int)) != cudaSuccess) goto cleanup;

    // upload packed H, S; unpack to dense symmetric squares.
    if (cudaMemcpy(dpk, h, ntri*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
    scf_unpack_tri<<<(unsigned)((nn+255)/256),256>>>(dpk, dH, nbf);
    if (cudaGetLastError()!=cudaSuccess) goto cleanup;
    if (cudaMemcpy(dpk, s, ntri*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
    scf_unpack_tri<<<(unsigned)((nn+255)/256),256>>>(dpk, dS, nbf);
    if (cudaGetLastError()!=cudaSuccess) goto cleanup;

    // workspace for cusolverDnDsyevd (reuse for all eighs: same nbf).
    if (cusolverDnDsyevd_bufferSize(hs, CUSOLVER_EIG_MODE_VECTOR,
            CUBLAS_FILL_MODE_LOWER, nbf, dX, nbf, dW, &lwork)
            != CUSOLVER_STATUS_SUCCESS) goto cleanup;
    if (cudaMalloc((void**)&dwork, (size_t)lwork*sizeof(double))!=cudaSuccess) goto cleanup;

    // ---- X = S^{-1/2} once (identical to RHF path).
    if (cudaMemcpy(dX, dS, nn*sizeof(double), cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
    if (cusolverDnDsyevd(hs, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
            nbf, dX, nbf, dW, dwork, lwork, dInfo)!=CUSOLVER_STATUS_SUCCESS) goto cleanup;
    {
        std::vector<double> sev(nbf);
        if (cudaMemcpy(sev.data(), dW, nbf*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup;
        for (int i=0;i<nbf;++i) sev[i] = (sev[i]>1e-12)? 1.0/sqrt(sev[i]) : 0.0;
        if (cudaMemcpy(dW, sev.data(), nbf*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
    }
    scf_scale_cols<<<(unsigned)((nn+255)/256),256>>>(dX, dW, dT, nbf);
    if (cudaGetLastError()!=cudaSuccess) goto cleanup;
    cs = cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_T, nbf, nbf, nbf, &one,
                     dT, nbf, dX, nbf, &zero, dT2, nbf);
    if (cs!=CUBLAS_STATUS_SUCCESS) goto cleanup;
    if (cudaMemcpy(dX, dT2, nn*sizeof(double), cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;

    // ---- SCF loop ----
    // Initial guess. If da0/db0 (packed lower-tri density guess) are supplied,
    // build Fa/Fb from them (one J/K pass) so the loop starts at the SAME state
    // as a reference SCF — this isolates loop-math correctness from guess-driven
    // multiple-UHF-solution landing. Otherwise use the core guess Fa=Fb=H.
    if (da0 && db0) {
        if (cudaMemcpy(dpk, da0, ntri*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
        scf_unpack_tri<<<(unsigned)((nn+255)/256),256>>>(dpk, dDab, nbf);
        if (cudaGetLastError()!=cudaSuccess) goto cleanup;
        if (cudaMemcpy(dpk, db0, ntri*sizeof(double), cudaMemcpyHostToDevice)!=cudaSuccess) goto cleanup;
        scf_unpack_tri<<<(unsigned)((nn+255)/256),256>>>(dpk, dDab+nn, nbf);
        if (cudaGetLastError()!=cudaSuccess) goto cleanup;
        if (!jkm_core(dDab, 2, dJab, dKab)) goto cleanup;
        if (cudaMemcpy(dJtot,dJab,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        cs=cublasDaxpy(hb,(int)nn,&one,dJab+nn,1,dJtot,1); if(cs)goto cleanup;
        if (cudaMemcpy(dFa,dH,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        { double a1=sc;  cs=cublasDaxpy(hb,(int)nn,&a1,dJtot,1,dFa,1); if(cs)goto cleanup; }
        { double a2=-se; cs=cublasDaxpy(hb,(int)nn,&a2,dKab,1,dFa,1);  if(cs)goto cleanup; }
        if (cudaMemcpy(dFb,dH,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        { double a1=sc;  cs=cublasDaxpy(hb,(int)nn,&a1,dJtot,1,dFb,1); if(cs)goto cleanup; }
        { double a2=-se; cs=cublasDaxpy(hb,(int)nn,&a2,dKab+nn,1,dFb,1); if(cs)goto cleanup; }
    } else {
        if (cudaMemcpy(dFa, dH, nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        if (cudaMemcpy(dFb, dH, nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
    }
    {
    double E=0.0, E_old=0.0; int it; int ndiis=0;
    std::vector<double> Bm;
    for (it=0; it<maxit; ++it) {
        // ---- alpha: Fp = X^T Fa X ; eigh -> Cp ; Ca = X Cp ; Da = Ca_occ Ca_occ^T
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dX,nbf,dFa,nbf,&zero,dT,nbf); if(cs)goto cleanup;
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dT,nbf,dX,nbf,&zero,dFp,nbf); if(cs)goto cleanup;
        if (cudaMemcpy(dCp,dFp,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        ss=cusolverDnDsyevd(hs,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,nbf,dCp,nbf,depsa,dwork,lwork,dInfo);
        if(ss!=CUSOLVER_STATUS_SUCCESS) goto cleanup;
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dX,nbf,dCp,nbf,&zero,dCa,nbf); if(cs)goto cleanup;
        if (na>0) {
            if (cudaMemcpy(dCo,dCa,(size_t)nbf*na*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_T,nbf,nbf,na,&one,dCo,nbf,dCo,nbf,&zero,dDab,nbf); if(cs)goto cleanup;
        } else { if(cudaMemset(dDab,0,nn*sizeof(double))!=cudaSuccess) goto cleanup; }
        // ---- beta:  Fp = X^T Fb X ; eigh -> Cp ; Cb = X Cp ; Db = Cb_occ Cb_occ^T
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dX,nbf,dFb,nbf,&zero,dT,nbf); if(cs)goto cleanup;
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dT,nbf,dX,nbf,&zero,dFp,nbf); if(cs)goto cleanup;
        if (cudaMemcpy(dCp,dFp,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        ss=cusolverDnDsyevd(hs,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,nbf,dCp,nbf,depsb,dwork,lwork,dInfo);
        if(ss!=CUSOLVER_STATUS_SUCCESS) goto cleanup;
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dX,nbf,dCp,nbf,&zero,dCb,nbf); if(cs)goto cleanup;
        if (nb>0) {
            if (cudaMemcpy(dCo,dCb,(size_t)nbf*nb*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_T,nbf,nbf,nb,&one,dCo,nbf,dCo,nbf,&zero,dDab+nn,nbf); if(cs)goto cleanup;
        } else { if(cudaMemset(dDab+nn,0,nn*sizeof(double))!=cudaSuccess) goto cleanup; }

        // ---- J,K for BOTH spins in one batched core call (nvec=2 over [Da|Db]).
        //   dJab[0..nn)=J(Da), dJab[nn..2nn)=J(Db); dKab likewise K(Da),K(Db).
        if (!jkm_core(dDab, 2, dJab, dKab)) goto cleanup;
        // Jtot = J(Da) + J(Db)
        if (cudaMemcpy(dJtot,dJab,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        cs=cublasDaxpy(hb,(int)nn,&one,dJab+nn,1,dJtot,1); if(cs)goto cleanup;
        // Fa = H + sc*Jtot - se*K(Da)
        if (cudaMemcpy(dFa,dH,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        { double a1=sc;    cs=cublasDaxpy(hb,(int)nn,&a1,dJtot,1,dFa,1); if(cs)goto cleanup; }
        { double a2=-se;   cs=cublasDaxpy(hb,(int)nn,&a2,dKab,1,dFa,1);  if(cs)goto cleanup; }
        // Fb = H + sc*Jtot - se*K(Db)
        if (cudaMemcpy(dFb,dH,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
        { double a1=sc;    cs=cublasDaxpy(hb,(int)nn,&a1,dJtot,1,dFb,1); if(cs)goto cleanup; }
        { double a2=-se;   cs=cublasDaxpy(hb,(int)nn,&a2,dKab+nn,1,dFb,1); if(cs)goto cleanup; }

        // ---- E = 0.5*( Tr[Da(H+Fa)] + Tr[Db(H+Fb)] ) + Enuc
        {
            double ea, eb;
            if (cudaMemcpy(dT,dH,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            cs=cublasDaxpy(hb,(int)nn,&one,dFa,1,dT,1); if(cs)goto cleanup;          // T = H+Fa
            cs=cublasDdot(hb,(int)nn,dDab,1,dT,1,&ea);  if(cs)goto cleanup;          // Tr[Da(H+Fa)]
            if (cudaMemcpy(dT,dH,nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            cs=cublasDaxpy(hb,(int)nn,&one,dFb,1,dT,1); if(cs)goto cleanup;          // T = H+Fb
            cs=cublasDdot(hb,(int)nn,dDab+nn,1,dT,1,&eb); if(cs)goto cleanup;        // Tr[Db(H+Fb)]
            E = 0.5*(ea+eb) + *enuc_;
        }

        // ---- combined CDIIS residual: e = [Fa Da S - S Da Fa ; Fb Db S - S Db Fb]
        // stored into the combined history slot (2*nn). Reuse dEhist slot as scratch.
        int slot = ndiis % mdiis;
        double* eslot = dEhist + (size_t)slot*2*nn;
        // alpha block
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dDab,nbf,dS,nbf,&zero,dT,nbf); if(cs)goto cleanup;  // T=Da S
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dFa,nbf,dT,nbf,&zero,dT2,nbf); if(cs)goto cleanup;  // T2=Fa Da S
        cs=cublasDgeam(hb,CUBLAS_OP_N,CUBLAS_OP_T,nbf,nbf,&one,dT2,nbf,&mone,dT2,nbf,eslot,nbf); if(cs)goto cleanup;   // e_a=T2-T2^T
        // beta block
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dDab+nn,nbf,dS,nbf,&zero,dT,nbf); if(cs)goto cleanup; // T=Db S
        cs=cublasDgemm(hb,CUBLAS_OP_N,CUBLAS_OP_N,nbf,nbf,nbf,&one,dFb,nbf,dT,nbf,&zero,dT2,nbf); if(cs)goto cleanup;    // T2=Fb Db S
        cs=cublasDgeam(hb,CUBLAS_OP_N,CUBLAS_OP_T,nbf,nbf,&one,dT2,nbf,&mone,dT2,nbf,eslot+nn,nbf); if(cs)goto cleanup;  // e_b
        double rnorm; cs=cublasDnrm2(hb,(int)(2*nn),eslot,1,&rnorm); if(cs)goto cleanup;
        // push combined Fock [Fa;Fb] into history slot
        {
            double* fslot = dFhist + (size_t)slot*2*nn;
            if (cudaMemcpy(fslot,    dFa, nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            if (cudaMemcpy(fslot+nn, dFb, nn*sizeof(double),cudaMemcpyDeviceToDevice)!=cudaSuccess) goto cleanup;
            ndiis++;
        }
        int nd = ndiis<mdiis ? ndiis : mdiis;
        if (nd>1) {
            Bm.assign((size_t)(nd+1)*(nd+1), -1.0); Bm[(size_t)(nd+1)*(nd+1)-1]=0.0;
            int start = ndiis<=mdiis ? 0 : ndiis%mdiis;
            for (int i=0;i<nd;++i){
                int si=(start+i)%mdiis;
                for(int j=0;j<nd;++j){
                    int sj=(start+j)%mdiis; double dij;
                    cs=cublasDdot(hb,(int)(2*nn),dEhist+(size_t)si*2*nn,1,dEhist+(size_t)sj*2*nn,1,&dij); if(cs)goto cleanup;
                    Bm[(size_t)i*(nd+1)+j]=dij;
                }
            }
            std::vector<double> A=Bm; std::vector<double> rhs(nd+1,0.0); rhs[nd]=-1.0;
            int N=nd+1; bool ok=true;
            for(int col=0;col<N&&ok;++col){
                int piv=col; double best=fabs(A[(size_t)col*N+col]);
                for(int r=col+1;r<N;++r){double v=fabs(A[(size_t)r*N+col]); if(v>best){best=v;piv=r;}}
                if(best<1e-300){ok=false;break;}
                if(piv!=col){ for(int k=0;k<N;++k) std::swap(A[(size_t)col*N+k],A[(size_t)piv*N+k]); std::swap(rhs[col],rhs[piv]); }
                double d=A[(size_t)col*N+col];
                for(int r=0;r<N;++r){ if(r==col)continue; double f=A[(size_t)r*N+col]/d;
                    for(int k=0;k<N;++k) A[(size_t)r*N+k]-=f*A[(size_t)col*N+k]; rhs[r]-=f*rhs[col]; }
            }
            if(ok){
                for(int i=0;i<N;++i) rhs[i]/=A[(size_t)i*N+i];
                // extrapolate Fa,Fb = sum_i c_i [Fa_i;Fb_i]
                if(cudaMemset(dFa,0,nn*sizeof(double))!=cudaSuccess) goto cleanup;
                if(cudaMemset(dFb,0,nn*sizeof(double))!=cudaSuccess) goto cleanup;
                for(int i=0;i<nd;++i){ int si=(start+i)%mdiis; double ci=rhs[i];
                    double* fslot = dFhist + (size_t)si*2*nn;
                    cs=cublasDaxpy(hb,(int)nn,&ci,fslot,1,   dFa,1); if(cs)goto cleanup;
                    cs=cublasDaxpy(hb,(int)nn,&ci,fslot+nn,1,dFb,1); if(cs)goto cleanup;
                }
            }
        }
        if (it>3 && fabs(E-E_old)<conv_e && rnorm<conv_d) { it++; break; }
        E_old=E;
    }
    *ncyc_out = it;
    *e_out = E;
    }
    // ---- return converged Da,Db (packed), + optional MO Ca,Cb / eps, ONCE ----
    scf_pack_tri<<<(unsigned)((ntri+255)/256),256>>>(dDab, dpk, nbf);
    if (cudaGetLastError()!=cudaSuccess) goto cleanup;
    if (cudaMemcpy(da_out, dpk, ntri*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup;
    scf_pack_tri<<<(unsigned)((ntri+255)/256),256>>>(dDab+nn, dpk, nbf);
    if (cudaGetLastError()!=cudaSuccess) goto cleanup;
    if (cudaMemcpy(db_out, dpk, ntri*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup;
    if (ca_out)   { if (cudaMemcpy(ca_out,   dCa,   nn*sizeof(double),  cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup; }
    if (cb_out)   { if (cudaMemcpy(cb_out,   dCb,   nn*sizeof(double),  cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup; }
    if (epsa_out) { if (cudaMemcpy(epsa_out, depsa, nbf*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup; }
    if (epsb_out) { if (cudaMemcpy(epsb_out, depsb, nbf*sizeof(double), cudaMemcpyDeviceToHost)!=cudaSuccess) goto cleanup; }
    *info = 0;

cleanup:
    #undef ALLOC
    if(dH)cudaFree(dH); if(dS)cudaFree(dS); if(dX)cudaFree(dX);
    if(dFa)cudaFree(dFa); if(dFb)cudaFree(dFb); if(dFp)cudaFree(dFp);
    if(dCa)cudaFree(dCa); if(dCb)cudaFree(dCb); if(dCp)cudaFree(dCp);
    if(dDab)cudaFree(dDab); if(dCo)cudaFree(dCo);
    if(dJab)cudaFree(dJab); if(dKab)cudaFree(dKab); if(dJtot)cudaFree(dJtot);
    if(dT)cudaFree(dT); if(dT2)cudaFree(dT2); if(dpk)cudaFree(dpk);
    if(dW)cudaFree(dW); if(depsa)cudaFree(depsa); if(depsb)cudaFree(depsb);
    if(dEhist)cudaFree(dEhist); if(dFhist)cudaFree(dFhist);
    if(dInfo)cudaFree(dInfo); if(dwork)cudaFree(dwork);
    if(hs)cusolverDnDestroy(hs);
}

}  // extern "C"
