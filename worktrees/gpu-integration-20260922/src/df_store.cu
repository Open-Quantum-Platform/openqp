// df_store.cu — device-resident DF-tensor store, keyed by an integer handle.
//
// The foundation for passing the density-fitting tensor (and, later, densities
// and shells) through OpenQP's tagarray instead of disk files: OpenQP puts a
// tensor once, the library holds it RESIDENT on GPU RAM, and every seam
// (SCF J/K, DFT Vxc, MRSF sigma, gradient) shares that one device copy with no
// re-upload and no per-run file I/O. `routec_df_to_host` mirrors it back so the
// tag stays visible from the Python layer.
//
// Thread-safe; handles are opaque ints (>=1). Freeing a handle releases the
// device memory. A geometry change is signalled by freeing + re-putting.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <mutex>
#include <sys/stat.h>
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include "routec_df.h"

namespace {
struct Entry { double* d = nullptr; long naux = 0, ncol = 0; };
std::map<int, Entry> g_store;
std::mutex           g_mtx;
int                  g_next = 1;

// A dedicated SINGLE resident B, shared by the two energy engines (SCF J/K in
// scf.cu and the MRSF sigma session in sigma.cu) so an MRSF run loads+uploads
// the ~GB density-fitting tensor ONCE, not once per consumer. Keyed by the
// source-file identity (path|mtime|size) so a geometry change auto-invalidates.
double*     g_B      = nullptr;
long        g_B_naux = 0, g_B_nbf = 0;
std::string g_B_key;
std::mutex  g_B_mtx;
// multi-GPU aux split: after routec_b_split_mg, g_B holds only aux [0,n0) on
// device 0 and g_B1 holds aux [n0,naux) on device 1. n0==0 means "not split".
double*     g_B1     = nullptr;
long        g_B_n0   = 0;
}  // namespace

extern "C" {

// Upload an (naux x ncol) row-major tensor to the device and register it.
// Returns a handle >= 1, or -1 on error.
int routec_df_put(const double* host, int naux, int ncol) {
  if (!host || naux <= 0 || ncol <= 0) return -1;
  const size_t bytes = (size_t)naux * ncol * sizeof(double);
  double* d = nullptr;
  if (cudaMalloc(&d, bytes) != cudaSuccess) {
    fprintf(stderr, "[df_store] cudaMalloc %.1f MB failed\n", bytes / 1048576.0);
    return -1;
  }
  if (cudaMemcpy(d, host, bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
    cudaFree(d); return -1;
  }
  std::lock_guard<std::mutex> lk(g_mtx);
  int h = g_next++;
  g_store[h] = {d, (long)naux, (long)ncol};
  fprintf(stderr, "[df_store] handle %d: %ld x %ld resident (%.1f MB)\n",
          h, (long)naux, (long)ncol, bytes / 1048576.0);
  return h;
}

// Device pointer for a handle (nullptr if unknown). Optionally returns dims.
const double* routec_df_dptr(int handle, int* naux, int* ncol) {
  std::lock_guard<std::mutex> lk(g_mtx);
  auto it = g_store.find(handle);
  if (it == g_store.end()) return nullptr;
  if (naux) *naux = (int)it->second.naux;
  if (ncol) *ncol = (int)it->second.ncol;
  return it->second.d;
}

// Mirror the resident tensor back to a host buffer (naux*ncol doubles).
// 0 on success, nonzero otherwise. Used to keep the tagarray/Python copy live.
int routec_df_to_host(int handle, double* out) {
  std::lock_guard<std::mutex> lk(g_mtx);
  auto it = g_store.find(handle);
  if (it == g_store.end() || !out) return 1;
  const size_t bytes = (size_t)it->second.naux * it->second.ncol * sizeof(double);
  return cudaMemcpy(out, it->second.d, bytes, cudaMemcpyDeviceToHost) == cudaSuccess ? 0 : 2;
}

// Release the device memory for a handle.
void routec_df_free(int handle) {
  std::lock_guard<std::mutex> lk(g_mtx);
  auto it = g_store.find(handle);
  if (it == g_store.end()) return;
  if (it->second.d) cudaFree(it->second.d);
  g_store.erase(it);
}

// ---- shared resident B (the SCF J/K + MRSF sigma choke point) --------------
// Load the density-fitting tensor from OQP_ROUTEC_B ONCE and hold it resident,
// so both energy engines SHARE one device copy. Returns a BORROWED device
// pointer to dense B, laid out (naux, nbf, nbf) row-major, symmetric per aux
// slice -- the store owns it; consumers must NOT cudaFree it. Keyed by the
// file identity so a changed geometry (new/rewritten B) auto-frees the old one
// and reloads. Returns nullptr on error. (Later, integral-direct mode will
// populate this same buffer/handle from an in-process 3c rebuild instead of a
// file, so the consumers never change.)
void routec_b_adopt(double* dB_device, int naux, int nbf) {
  std::lock_guard<std::mutex> lk(g_B_mtx);
  if (g_B && g_B != dB_device) cudaFree(g_B);
  g_B = dB_device; g_B_naux = naux; g_B_nbf = nbf; g_B_key = "ADOPTED";
  fprintf(stderr, "[df_store] B adopted in-process (naux=%d nbf=%d, %.1f MB)\n",
          naux, nbf, (double)naux*nbf*nbf*8/1048576.0);
}

// Split the RESIDENT B across two GPUs along the aux dimension: device 0 keeps
// aux [0, n0), device 1 gets aux [n0, naux) (contiguous slices of the
// (naux, nbf, nbf) layout, so this is two flat copies). Device 0's residence is
// then SHRUNK to its half -- the per-GPU footprint of B halves, which is what
// lets bigger systems fit. Requires bidirectional P2P (NVLink on chc4).
// Returns 1 with outputs filled if split (now or already), 0 if unavailable.
// After a split routec_b_shared() refuses to hand out the half tensor.
int routec_b_split_mg(int* n0_out, int* n1_out, int* nbf_out,
                      const double** B0, const double** B1) {
  std::lock_guard<std::mutex> lk(g_B_mtx);
  if (!g_B) return 0;
  const long nn = g_B_nbf * g_B_nbf;
  if (g_B_n0 == 0) {
    int ndev = 0;
    if (cudaGetDeviceCount(&ndev) != cudaSuccess || ndev < 2) return 0;
    int can01 = 0, can10 = 0;
    cudaDeviceCanAccessPeer(&can01, 0, 1);
    cudaDeviceCanAccessPeer(&can10, 1, 0);
    if (!can01 || !can10) {
      fprintf(stderr, "[df_store] mg: no P2P between GPU0/1, staying single-GPU\n");
      return 0;
    }
    cudaSetDevice(0); cudaDeviceEnablePeerAccess(1, 0); cudaGetLastError();
    cudaSetDevice(1); cudaDeviceEnablePeerAccess(0, 0); cudaGetLastError();
    const long n0 = g_B_naux / 2, n1 = g_B_naux - n0;
    double* b1 = nullptr;                             // dev1 half
    if (cudaMalloc(&b1, (size_t)n1*nn*8) != cudaSuccess) {
      cudaSetDevice(0);
      fprintf(stderr, "[df_store] mg: dev1 alloc failed, staying single-GPU\n");
      return 0;
    }
    if (cudaMemcpyPeer(b1, 1, g_B + (size_t)n0*nn, 0, (size_t)n1*nn*8) != cudaSuccess) {
      cudaFree(b1); cudaSetDevice(0); return 0;
    }
    cudaSetDevice(0);
    double* b0 = nullptr;                             // shrink dev0 to its half
    if (cudaMalloc(&b0, (size_t)n0*nn*8) == cudaSuccess &&
        cudaMemcpy(b0, g_B, (size_t)n0*nn*8, cudaMemcpyDeviceToDevice) == cudaSuccess) {
      cudaFree(g_B); g_B = b0;
    } else if (b0) { cudaFree(b0); }   // keep full B on dev0 (no shrink); split still valid
    g_B1 = b1; g_B_n0 = n0;
    fprintf(stderr, "[df_store] multi-GPU B split: dev0 %ld + dev1 %ld aux slices (%.1f MB each)\n",
            n0, n1, (double)n0*nn*8/1048576.0);
  }
  if (n0_out) *n0_out = (int)g_B_n0;
  if (n1_out) *n1_out = (int)(g_B_naux - g_B_n0);
  if (nbf_out) *nbf_out = (int)g_B_nbf;
  if (B0) *B0 = g_B;
  if (B1) *B1 = g_B1;
  return 1;
}

const double* routec_b_shared(int* naux_out, int* nbf_out) {
  std::lock_guard<std::mutex> lk(g_B_mtx);
  if (g_B && g_B_n0 > 0) {
    fprintf(stderr, "[df_store] B is multi-GPU split; this consumer needs routec_b_split_mg\n");
    return nullptr;
  }
  if (g_B && g_B_key == "ADOPTED") {              // in-process adopted tensor
    if (naux_out) *naux_out = (int)g_B_naux;
    if (nbf_out)  *nbf_out  = (int)g_B_nbf;
    return g_B;
  }
  const char* path = getenv("OQP_ROUTEC_B");
  if (!path) { fprintf(stderr, "[df_store] OQP_ROUTEC_B unset\n"); return nullptr; }
  struct stat st;
  if (stat(path, &st) != 0) { fprintf(stderr, "[df_store] cannot stat %s\n", path); return nullptr; }
  std::string key = std::string(path) + "|" + std::to_string((long long)st.st_mtime)
                  + "|" + std::to_string((long long)st.st_size);
  if (g_B && g_B_key == key) {                    // resident + same geometry
    if (naux_out) *naux_out = (int)g_B_naux;
    if (nbf_out)  *nbf_out  = (int)g_B_nbf;
    fprintf(stderr, "[df_store] B shared: reuse resident (naux=%ld nbf=%ld)\n",
            g_B_naux, g_B_nbf);
    return g_B;
  }
  if (g_B) { cudaFree(g_B); g_B = nullptr; g_B_key.clear(); }   // geometry changed
  FILE* f = fopen(path, "rb");
  if (!f) { fprintf(stderr, "[df_store] cannot open %s\n", path); return nullptr; }
  const size_t fbytes = (size_t)st.st_size;
  int hdr[2];
  if (fread(hdr, sizeof(int), 2, f) != 2) { fclose(f); return nullptr; }
  // ---- CDF v2 compacted: {magic, naux, nao, ncp} | keep[ncp] | (naux x ncp) --
  // Read the compacted B and scatter it back to a dense (naux, nao, nao) buffer
  // so the J/K/sigma consumers are unchanged (Stage ii). The compaction is
  // lossless (dropped columns were exact zeros), so energies must be identical.
  if (hdr[0] == 0x32464443 /*"CDF2"*/) {
    const int naux = hdr[1]; int nao = 0, ncp = 0;
    if (fread(&nao,4,1,f)!=1 || fread(&ncp,4,1,f)!=1) { fclose(f); return nullptr; }
    std::vector<int> keep(ncp);
    if (ncp>0 && fread(keep.data(),4,ncp,f)!=(size_t)ncp) { fclose(f); return nullptr; }
    std::vector<double> cB((size_t)naux*ncp);
    if (cB.size() && fread(cB.data(),8,cB.size(),f)!=cB.size()) { fclose(f); return nullptr; }
    fclose(f);
    const size_t nn = (size_t)nao * nao;
    std::vector<double> dense((size_t)naux*nn, 0.0);
    for (int c=0;c<ncp;++c) {
      const long long t = keep[c];                     // lower-tri index i*(i+1)/2 + j
      int i = (int)((std::sqrt(8.0*(double)t+1.0)-1.0)/2.0);
      while ((long long)(i+1)*(i+2)/2 <= t) ++i;        // correct fp rounding
      while ((long long)i*(i+1)/2 > t) --i;
      const int j = (int)(t - (long long)i*(i+1)/2);
      for (int P=0;P<naux;++P) {
        const double v = cB[(size_t)P*ncp+c];
        dense[(size_t)P*nn + (size_t)i*nao + j] = v;
        dense[(size_t)P*nn + (size_t)j*nao + i] = v;
      }
    }
    double* d=nullptr; const size_t bytes=(size_t)naux*nn*sizeof(double);
    if (cudaMalloc(&d,bytes)!=cudaSuccess) {
      fprintf(stderr,"[df_store] cudaMalloc %.1f MB (CDF v2) failed\n",bytes/1048576.0);
      return nullptr; }
    if (cudaMemcpy(d,dense.data(),bytes,cudaMemcpyHostToDevice)!=cudaSuccess) {
      cudaFree(d); return nullptr; }
    g_B=d; g_B_naux=naux; g_B_nbf=nao; g_B_key=key;
    if (naux_out) *naux_out=naux;
    if (nbf_out)  *nbf_out=nao;
    fprintf(stderr,"[df_store] B shared (CDF v2): naux=%d nbf=%d ncp=%d/%zu -> dense %.1f MB\n",
            naux,nao,ncp,nn?(size_t)nao*(nao+1)/2:0,bytes/1048576.0);
    return g_B;
  }
  bool packed = false; int naux = 0, nbf = 0;
  if      (hdr[0] > 0 && hdr[1] > 0) { naux =  hdr[0]; nbf =  hdr[1]; }
  else if (hdr[0] < 0 && hdr[1] > 0) { packed = true; nbf = -hdr[0]; naux =  hdr[1]; }
  else if (hdr[0] > 0 && hdr[1] < 0) { packed = true; naux =  hdr[0]; nbf = -hdr[1]; }
  else { fprintf(stderr, "[df_store] bad B header (%d,%d) in %s\n", hdr[0], hdr[1], path);
         fclose(f); return nullptr; }
  const size_t nn    = (size_t)nbf * nbf;
  const size_t npair = (size_t)nbf * (nbf + 1) / 2;
  const size_t nread = packed ? (size_t)naux * npair : (size_t)naux * nn;
  if (fbytes != 2 * sizeof(int) + nread * sizeof(double)) {
    fprintf(stderr, "[df_store] B size mismatch in %s\n", path); fclose(f); return nullptr; }
  std::vector<double> raw(nread);
  if (fread(raw.data(), sizeof(double), nread, f) != nread) { fclose(f); return nullptr; }
  fclose(f);
  std::vector<double> dense;
  const double* src = raw.data();
  if (packed) {                                   // host unpack tri rows -> dense
    dense.assign((size_t)naux * nn, 0.0);
    for (int P = 0; P < naux; ++P) {
      const double* rp = raw.data() + (size_t)P * npair;
      double* bp = dense.data() + (size_t)P * nn;
      size_t t = 0;
      for (int i = 0; i < nbf; ++i)
        for (int j = 0; j <= i; ++j, ++t) { bp[(size_t)i*nbf+j] = rp[t]; bp[(size_t)j*nbf+i] = rp[t]; }
    }
    src = dense.data();
  }
  double* d = nullptr;
  const size_t bytes = (size_t)naux * nn * sizeof(double);
  if (cudaMalloc(&d, bytes) != cudaSuccess) {
    fprintf(stderr, "[df_store] cudaMalloc %.1f MB for shared B failed\n", bytes / 1048576.0);
    return nullptr; }
  if (cudaMemcpy(d, src, bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
    cudaFree(d); return nullptr; }
  g_B = d; g_B_naux = naux; g_B_nbf = nbf; g_B_key = key;
  if (naux_out) *naux_out = naux;
  if (nbf_out)  *nbf_out  = nbf;
  fprintf(stderr, "[df_store] B shared: loaded naux=%d nbf=%d (%.1f MB dense, %s) from %s\n",
          naux, nbf, bytes / 1048576.0, packed ? "packed" : "dense", path);
  return g_B;
}

// Release the shared resident B (the seam calls this on a geometry change if it
// wants to force a reload before the file's mtime/size would show it).
void routec_b_free(void) {
  std::lock_guard<std::mutex> lk(g_B_mtx);
  if (g_B) cudaFree(g_B);
  if (g_B1) { int cur=0; cudaGetDevice(&cur); cudaSetDevice(1); cudaFree(g_B1); cudaSetDevice(cur); }
  g_B = nullptr; g_B1 = nullptr; g_B_n0 = 0; g_B_key.clear(); g_B_naux = g_B_nbf = 0;
}

// ---- compacted-on-device B (the real device-memory win) --------------------
// Load a CDF v2/v3 file and keep the COMPACTED tensor resident on the GPU plus
// the keep_pairs map -- WITHOUT scattering to dense. Consumers (tiled occ-RI-K /
// J / sigma) reconstruct dense aux-slices in small tiles on the fly, so the full
// dense B (naux x nbf^2) never exists on device. v2 = fp64 compacted (footprint
// drops by lower-tri 2x x compaction ncp/npair); v3 = magnitude-tiered mixed
// precision (fp64/fp16/int8), the compression multiplier on top. Returns 0 ok.
static double*        g_cB    = nullptr;  // fp64 arena (v2: all naux x ncp; v3: n64 x ncp)
static int*           g_keep  = nullptr;  // (ncp) original lower-tri pair indices
static void*          g_cB32  = nullptr;  // float  arena (n32 x ncp), v3
static void*          g_cB16  = nullptr;  // __half arena (n16 x ncp), v3
static void*          g_cB8   = nullptr;  // int8   arena (n8  x ncp), v3
static unsigned char* g_tier  = nullptr;  // (naux) per-column tier, v3
static int*           g_slot  = nullptr;  // (naux) row within arena, v3
static float*         g_scale = nullptr;  // (naux) per-column dequant scale, v3
static long g_cdf_naux=0, g_cdf_nao=0, g_cdf_ncp=0;
static int  g_cdf_lowprec=0, g_cdf_n64=0, g_cdf_n32=0, g_cdf_n16=0, g_cdf_n8=0;
static std::string g_cdf_key;

static void cdf_free_all() {
  if (g_cB)   { cudaFree(g_cB);   g_cB=nullptr; }
  if (g_keep) { cudaFree(g_keep); g_keep=nullptr; }
  if (g_cB32) { cudaFree(g_cB32); g_cB32=nullptr; }
  if (g_cB16) { cudaFree(g_cB16); g_cB16=nullptr; }
  if (g_cB8)  { cudaFree(g_cB8);  g_cB8=nullptr; }
  if (g_tier) { cudaFree(g_tier); g_tier=nullptr; }
  if (g_slot) { cudaFree(g_slot); g_slot=nullptr; }
  if (g_scale){ cudaFree(g_scale);g_scale=nullptr; }
}
static void cdf_fill_out(struct CdfDev* out) {
  out->naux=(int)g_cdf_naux; out->nao=(int)g_cdf_nao; out->ncp=(int)g_cdf_ncp;
  out->keep=g_keep; out->lowprec=g_cdf_lowprec;
  out->B64=g_cB; out->B32=g_cB32; out->B16=g_cB16; out->B8=g_cB8;
  out->tier=g_tier; out->slot=g_slot; out->scale=g_scale;
  out->n64=g_cdf_n64; out->n32=g_cdf_n32; out->n16=g_cdf_n16; out->n8=g_cdf_n8;
}

int routec_cdf_dev(struct CdfDev* out) {
  std::lock_guard<std::mutex> lk(g_B_mtx);
  const char* path = getenv("OQP_ROUTEC_B");
  if (!path) { fprintf(stderr, "[df_store] OQP_ROUTEC_B unset\n"); return 1; }
  struct stat st;
  if (stat(path, &st) != 0) { fprintf(stderr, "[df_store] cannot stat %s\n", path); return 1; }
  std::string key = std::string(path) + "|" + std::to_string((long long)st.st_mtime)
                  + "|" + std::to_string((long long)st.st_size);
  if (!g_cdf_key.empty() && g_cdf_key == key) { cdf_fill_out(out); return 0; }  // resident + same geom
  cdf_free_all(); g_cdf_key.clear();
  g_cdf_lowprec=0; g_cdf_n64=g_cdf_n16=g_cdf_n8=0;
  FILE* f = fopen(path, "rb");
  if (!f) { fprintf(stderr, "[df_store] cannot open %s\n", path); return 1; }
  int hdr[2];
  if (fread(hdr,sizeof(int),2,f)!=2) { fclose(f); return 1; }
  const bool is_v2 = (hdr[0]==0x32464443), is_v3 = (hdr[0]==0x33464443);
  if (!is_v2 && !is_v3) {
    // ---- DENSE file: compact on load --------------------------------------
    // Plain dense B {int naux, int nbf} + fp64 (naux x nbf x nbf), e.g. a
    // pyscf-built tensor with no exact-zero rows. OQP_CDF_LOAD_TAU (absolute
    // |B| threshold; 0 keeps every pair = lossless) drops the lower-tri pair
    // ROWS whose max |B| over aux is below tau, streaming the file in aux
    // slabs so the dense tensor never resides in host OR device memory.
    const char* te = getenv("OQP_CDF_LOAD_TAU");
    if (!te) {
      fprintf(stderr, "[df_store] routec_cdf_dev: dense B needs OQP_CDF_LOAD_TAU "
              "(or use a CDF v2/v3 file)\n");
      fclose(f); return 1; }
    const double tau = atof(te);
    const int naux_d = hdr[0], nbf = hdr[1];
    if (naux_d <= 0 || nbf <= 0 || naux_d > 1000000 || nbf > 100000) {
      fprintf(stderr, "[df_store] dense B header looks wrong (%d, %d)\n", naux_d, nbf);
      fclose(f); return 1; }
    const long long nn = (long long)nbf*nbf, ntri = (long long)nbf*(nbf+1)/2;
    const int SLAB = 16;
    std::vector<double> slab((size_t)SLAB*nn);
    // pass 1: per lower-tri pair, max |B| over all aux
    std::vector<double> colmax(ntri, 0.0);
    for (long long p0=0; p0<naux_d; p0+=SLAB) {
      const int nt = (int)std::min<long long>(SLAB, naux_d-p0);
      if (fread(slab.data(),8,(size_t)nt*nn,f)!=(size_t)nt*nn) { fclose(f); return 1; }
      for (int a=0; a<nt; ++a) {
        const double* Bp = slab.data()+(size_t)a*nn;
        long long t=0;
        for (int i=0;i<nbf;++i)
          for (int j=0;j<=i;++j,++t) {
            const double v = fabs(Bp[(long long)i*nbf+j]);
            if (v > colmax[t]) colmax[t] = v;
          }
      }
    }
    std::vector<int> keepv; keepv.reserve((size_t)ntri);
    for (long long t=0;t<ntri;++t) if (colmax[t] >= tau) keepv.push_back((int)t);
    const int ncp_d = (int)keepv.size();
    if (ncp_d == 0) { fprintf(stderr,"[df_store] LOAD_TAU=%.1e drops ALL pairs\n",tau); fclose(f); return 1; }
    std::vector<int> ki(ncp_d), kj(ncp_d);            // unpack tri once on host
    for (int c=0;c<ncp_d;++c) {
      long long t=keepv[c];
      int i=(int)((sqrt(8.0*(double)t+1.0)-1.0)*0.5);
      while ((long long)(i+1)*(i+2)/2 <= t) ++i;
      while ((long long)i*(i+1)/2 > t) --i;
      ki[c]=i; kj[c]=(int)(t-(long long)i*(i+1)/2);
    }
    const size_t cbytes = (size_t)naux_d*ncp_d*8;
    if (cudaMalloc(&g_keep,(size_t)ncp_d*4)!=cudaSuccess ||
        cudaMalloc(&g_cB, cbytes)!=cudaSuccess) {
      fprintf(stderr,"[df_store] cudaMalloc compacted B (%.1f MB) failed\n",cbytes/1048576.0);
      cdf_free_all(); fclose(f); return 1; }
    cudaMemcpy(g_keep,keepv.data(),(size_t)ncp_d*4,cudaMemcpyHostToDevice);
    // pass 2: gather kept pairs per aux row, upload slab-wise
    if (fseek(f, 8, SEEK_SET)!=0) { cdf_free_all(); fclose(f); return 1; }
    std::vector<double> rows((size_t)SLAB*ncp_d);
    for (long long p0=0; p0<naux_d; p0+=SLAB) {
      const int nt = (int)std::min<long long>(SLAB, naux_d-p0);
      if (fread(slab.data(),8,(size_t)nt*nn,f)!=(size_t)nt*nn) { cdf_free_all(); fclose(f); return 1; }
      for (int a=0; a<nt; ++a) {
        const double* Bp = slab.data()+(size_t)a*nn;
        double* r = rows.data()+(size_t)a*ncp_d;
        for (int c=0;c<ncp_d;++c) r[c] = Bp[(long long)ki[c]*nbf + kj[c]];
      }
      cudaMemcpy(g_cB+(size_t)p0*ncp_d, rows.data(), (size_t)nt*ncp_d*8,
                 cudaMemcpyHostToDevice);
    }
    fclose(f);
    g_cdf_naux=naux_d; g_cdf_nao=nbf; g_cdf_ncp=ncp_d; g_cdf_key=key;
    g_cdf_lowprec=0; g_cdf_n64=naux_d; g_cdf_n32=g_cdf_n16=g_cdf_n8=0;
    fprintf(stderr,"[df_store] B on device CDF (dense compacted on load): naux=%d nbf=%d "
            "ncp=%d/%lld (%.2f kept, tau=%.1e) -> %.1f MB (vs %.1f dense)\n",
            naux_d,nbf,ncp_d,ntri,(double)ncp_d/(double)ntri,tau,
            cbytes/1048576.0,(double)naux_d*nn*8/1048576.0);
    cdf_fill_out(out); return 0;
  }
  const int naux = hdr[1]; int nao=0, ncp=0;
  if (fread(&nao,4,1,f)!=1 || fread(&ncp,4,1,f)!=1) { fclose(f); return 1; }
  std::vector<int> keep(ncp);
  if (ncp>0 && fread(keep.data(),4,ncp,f)!=(size_t)ncp) { fclose(f); return 1; }

  if (is_v3) {                                    // ---- mixed-precision arenas ----
    std::vector<unsigned char> tier(naux); std::vector<int> slot(naux);
    std::vector<float> scale(naux); int n64=0,n32=0,n16=0,n8=0;
    if (fread(tier.data(),1,naux,f)!=(size_t)naux ||
        fread(slot.data(),4,naux,f)!=(size_t)naux ||
        fread(scale.data(),4,naux,f)!=(size_t)naux ||
        fread(&n64,4,1,f)!=1 || fread(&n32,4,1,f)!=1 ||
        fread(&n16,4,1,f)!=1 || fread(&n8,4,1,f)!=1) { fclose(f); return 1; }
    std::vector<double> B64((size_t)n64*ncp);
    std::vector<float>  B32((size_t)n32*ncp);
    std::vector<unsigned short> B16((size_t)n16*ncp);   // __half bits, 2 bytes
    std::vector<signed char> B8((size_t)n8*ncp);
    if ((B64.size()&&fread(B64.data(),8,B64.size(),f)!=B64.size()) ||
        (B32.size()&&fread(B32.data(),4,B32.size(),f)!=B32.size()) ||
        (B16.size()&&fread(B16.data(),2,B16.size(),f)!=B16.size()) ||
        (B8.size() &&fread(B8.data(), 1,B8.size(), f)!=B8.size())) { fclose(f); return 1; }
    fclose(f);
    const size_t b64=(size_t)n64*ncp*8, b32=(size_t)n32*ncp*4, b16=(size_t)n16*ncp*2, b8=(size_t)n8*ncp*1;
    bool ok = cudaMalloc(&g_keep,(size_t)ncp*4)==cudaSuccess
           && cudaMalloc(&g_tier,(size_t)naux*1)==cudaSuccess
           && cudaMalloc(&g_slot,(size_t)naux*4)==cudaSuccess
           && cudaMalloc(&g_scale,(size_t)naux*4)==cudaSuccess
           && (b64==0 || cudaMalloc(&g_cB,b64)==cudaSuccess)
           && (b32==0 || cudaMalloc(&g_cB32,b32)==cudaSuccess)
           && (b16==0 || cudaMalloc(&g_cB16,b16)==cudaSuccess)
           && (b8 ==0 || cudaMalloc(&g_cB8, b8 )==cudaSuccess);
    if (!ok) { fprintf(stderr,"[df_store] cudaMalloc v3 arenas failed\n"); cdf_free_all(); return 1; }
    cudaMemcpy(g_keep,keep.data(),(size_t)ncp*4,cudaMemcpyHostToDevice);
    cudaMemcpy(g_tier,tier.data(),(size_t)naux*1,cudaMemcpyHostToDevice);
    cudaMemcpy(g_slot,slot.data(),(size_t)naux*4,cudaMemcpyHostToDevice);
    cudaMemcpy(g_scale,scale.data(),(size_t)naux*4,cudaMemcpyHostToDevice);
    if (b64) cudaMemcpy(g_cB, B64.data(),b64,cudaMemcpyHostToDevice);
    if (b32) cudaMemcpy(g_cB32,B32.data(),b32,cudaMemcpyHostToDevice);
    if (b16) cudaMemcpy(g_cB16,B16.data(),b16,cudaMemcpyHostToDevice);
    if (b8)  cudaMemcpy(g_cB8, B8.data(), b8, cudaMemcpyHostToDevice);
    g_cdf_naux=naux; g_cdf_nao=nao; g_cdf_ncp=ncp; g_cdf_key=key;
    g_cdf_lowprec=1; g_cdf_n64=n64; g_cdf_n32=n32; g_cdf_n16=n16; g_cdf_n8=n8;
    const double mb=(b64+b32+b16+b8+(size_t)ncp*4+(size_t)naux*9)/1048576.0;
    fprintf(stderr,"[df_store] B on device CDF-v3 (mixed): naux=%d nao=%d ncp=%d "
            "tiers[fp64/fp32/fp16/int8]=%d/%d/%d/%d -> %.1f MB (vs %.1f fp64-compacted, %.1f dense)\n",
            naux,nao,ncp,n64,n32,n16,n8,mb,(double)naux*ncp*8/1048576.0,
            (double)naux*nao*nao*8/1048576.0);
    cdf_fill_out(out); return 0;
  }

  // ---- v2: fp64 compacted (with optional host accuracy probes) ----
  std::vector<double> cB((size_t)naux*ncp);
  if (cB.size() && fread(cB.data(),8,cB.size(),f)!=cB.size()) { fclose(f); return 1; }
  fclose(f);
  // ACCURACY PROBE: OQP_CDF_FP32 round-trips every element through fp32 (quantize)
  // to measure the ENERGY impact of uniform fp32 storage. Device stays fp64.
  if (getenv("OQP_CDF_FP32")) {
    for (size_t i=0;i<cB.size();++i) cB[i] = (double)(float)cB[i];
    fprintf(stderr,"[df_store] CDF fp32-quantize probe (accuracy test)\n");
  }
  // ACCURACY PROBE: OQP_CDF_LOWPREC round-trips each aux column through a
  // MAGNITUDE-TIERED low-precision store (fp64 / fp16 / int8, per-column block
  // scale) IN PLACE, to measure the energy impact of the precision tiers before
  // committing to a v3 build. Device stays fp64 here (accuracy, not memory). Same
  // tiering rule the v3 builder bakes in: tier by column absmax s_P vs global max
  // S (s>=TH1*S -> fp64; TH2*S<=s<TH1*S -> fp16 pow2 scale; s<TH2*S -> int8).
  if (getenv("OQP_CDF_LOWPREC")) {
    // 4-tier decision, identical to build_df's v3 emit. fp32 is the safe
    // workhorse (2x); fp16/int8 only for the small-magnitude tail.
    const double R64 = getenv("OQP_CDF_TH64") ? atof(getenv("OQP_CDF_TH64")) : 0.5;
    const double R32 = getenv("OQP_CDF_TH32") ? atof(getenv("OQP_CDF_TH32")) : 1e-3;
    const double R16 = getenv("OQP_CDF_TH16") ? atof(getenv("OQP_CDF_TH16")) : 1e-5;
    const double A64 = getenv("OQP_CDF_A64") ? atof(getenv("OQP_CDF_A64")) : -1.0;
    const double A32 = getenv("OQP_CDF_A32") ? atof(getenv("OQP_CDF_A32")) : -1.0;
    const double A16 = getenv("OQP_CDF_A16") ? atof(getenv("OQP_CDF_A16")) : -1.0;
    const bool NOINT8 = getenv("OQP_CDF_NOINT8")!=nullptr, NOFP16 = getenv("OQP_CDF_NOFP16")!=nullptr;
    std::vector<double> smax(naux, 0.0);
    double S = 0.0;
    for (int P=0;P<naux;++P) {
      double m=0.0; const double* col=&cB[(size_t)P*ncp];
      for (int c=0;c<ncp;++c){ double a=fabs(col[c]); if(a>m) m=a; }
      smax[P]=m; if(m>S) S=m;
    }
    const double b64=(A64>0.0)?A64:R64*S, b32=(A32>0.0)?A32:R32*S, b16=(A16>0.0)?A16:R16*S;
    long n64=0,n32=0,n16=0,n8=0;
    for (int P=0;P<naux;++P) {
      const double s=smax[P]; double* col=&cB[(size_t)P*ncp];
      if (s==0.0 || s>=b64) { ++n64; continue; }             // fp64: leave exact
      if (s>=b32 || NOFP16) {                                // fp32
        for (int c=0;c<ncp;++c) col[c] = (double)(float)col[c];
        ++n32;
      } else if (NOINT8 || s>=b16) {                         // fp16, pow2 block scale
        const double scale = exp2(round(log2(s)));
        for (int c=0;c<ncp;++c)
          col[c] = (double)__half2float(__float2half((float)(col[c]/scale)))*scale;
        ++n16;
      } else {                                               // int8 symmetric
        const double scale8 = s/127.0, inv=1.0/scale8;
        for (int c=0;c<ncp;++c){ double q=nearbyint(col[c]*inv);
          if(q>127.0)q=127.0; else if(q<-127.0)q=-127.0; col[c]=q*scale8; }
        ++n8;
      }
    }
    const double bpc = (8.0*n64+4.0*n32+2.0*n16+1.0*n8)/(double)naux;
    fprintf(stderr,"[df_store routec] CDF LOWPREC probe: fp64=%ld fp32=%ld fp16=%ld int8=%ld / %d "
            "(b64=%.1e b32=%.1e b16=%.1e S=%.1e) -> %.2f B/col = %.2fx vs fp64\n",
            n64,n32,n16,n8,naux,b64,b32,b16,S,bpc,8.0/bpc);
  }
  const size_t bBytes=(size_t)naux*ncp*8, kBytes=(size_t)ncp*4;
  if (cudaMalloc(&g_cB,bBytes)!=cudaSuccess ||
      cudaMalloc(&g_keep,kBytes)!=cudaSuccess) {
    fprintf(stderr,"[df_store] cudaMalloc compacted (%.1f MB) failed\n",bBytes/1048576.0);
    cdf_free_all(); return 1; }
  cudaMemcpy(g_cB,cB.data(),bBytes,cudaMemcpyHostToDevice);
  cudaMemcpy(g_keep,keep.data(),kBytes,cudaMemcpyHostToDevice);
  g_cdf_naux=naux; g_cdf_nao=nao; g_cdf_ncp=ncp; g_cdf_key=key;
  g_cdf_lowprec=0; g_cdf_n64=naux; g_cdf_n16=0; g_cdf_n8=0;
  fprintf(stderr,"[df_store] B on device COMPACTED: naux=%d nao=%d ncp=%d (%.1f MB vs %.1f dense = %.1fx)\n",
          naux,nao,ncp,bBytes/1048576.0,(double)naux*nao*nao*8/1048576.0,
          ncp?((double)naux*nao*nao)/((double)naux*ncp):1.0);
  cdf_fill_out(out); return 0;
}

}  // extern "C"
