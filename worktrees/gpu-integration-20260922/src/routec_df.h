// routec_df.h -- device-resident DF-tensor store shared across the energy
// engines. df_store.cu defines these; scf.cu (HF/DFT J/K) and sigma.cu (MRSF
// sigma) call routec_b_shared() so the B tensor is loaded+uploaded ONCE and
// shared, instead of each fopen'ing OQP_ROUTEC_B separately.
#pragma once

extern "C" {
// Borrowed device pointer to the single resident dense B (naux, nbf, nbf),
// row-major, symmetric per aux slice. Owned by df_store -- do NOT cudaFree it.
// Loaded once from OQP_ROUTEC_B (keyed by file identity), nullptr on error.
const double* routec_b_shared(int* naux, int* nbf);
// Multi-GPU aux split of the resident B (OQP_MULTI_GPU): device 0 keeps aux
// [0,n0), device 1 gets [n0,naux); dev0 residence SHRINKS to its half (the
// per-GPU B footprint halves). Performs the split on first call. Returns 1
// with outputs filled when split, 0 when unavailable (1 GPU / no P2P / no B).
// After a split, routec_b_shared() refuses (consumers must be split-aware).
int routec_b_split_mg(int* n0, int* n1, int* nbf,
                      const double** B0, const double** B1);
// Force-release the shared B (e.g. the seam on a geometry change).
void routec_b_free(void);
// In-process handoff: ADOPT an already-device-resident FULL dense B
// (naux, nbf, nbf). df_store takes ownership (cudaFree on release);
// routec_b_shared() then returns it without touching OQP_ROUTEC_B.
void routec_b_adopt(double* dB_device, int naux, int nbf);

// Compacted-on-device B (real device-memory win): load a CDF v2/v3 file and keep
// the compacted (naux x ncp) tensor + keep_pairs map resident WITHOUT expanding
// to dense. Consumers reconstruct dense aux-slices in small tiles. 0 on success.
//
// v2 = pure fp64 (lowprec=0: B64 is the whole naux x ncp fp64, tier/slot/scale
// null). v3 = magnitude-tiered mixed precision (lowprec=1): each aux column P
// lives in exactly one arena (B64 fp64 / B32 float / B16 __half / B8 int8) at row
// slot[P], with per-column tier[P] (0/1/2/3) and dequant scale[P] (fp64/fp32: 1,
// fp16: pow2, int8: absmax/127). fp32 is the safe workhorse tier (rel err ~6e-8,
// 2x); fp16/int8 compress the small-magnitude tail further. All pointers are
// BORROWED device pointers owned by df_store -- do NOT free.
struct CdfDev {
  int naux, nao, ncp;
  const int*  keep;             // [ncp] original lower-tri pair indices (device)
  int lowprec;                  // 0 = v2 (fp64 only), 1 = v3 (mixed precision)
  const double* B64;            // [n64 x ncp] fp64  arena (v2: n64=naux = all of B)
  const void*   B32;            // [n32 x ncp] float arena (v3 only, else null)
  const void*   B16;            // [n16 x ncp] __half arena (v3 only, else null)
  const void*   B8;             // [n8  x ncp] int8  arena (v3 only, else null)
  const unsigned char* tier;    // [naux] 0=fp64 1=fp32 2=fp16 3=int8 (v3, else null)
  const int*    slot;           // [naux] row of the column within its arena (v3)
  const float*  scale;          // [naux] per-column dequant scale (v3)
  int n64, n32, n16, n8;
};
int routec_cdf_dev(struct CdfDev* out);

// Generic handle-based store (foundation for tagarray / integral-direct).
int           routec_df_put(const double* host, int naux, int ncol);
const double* routec_df_dptr(int handle, int* naux, int* ncol);
int           routec_df_to_host(int handle, double* out);
void          routec_df_free(int handle);
}
