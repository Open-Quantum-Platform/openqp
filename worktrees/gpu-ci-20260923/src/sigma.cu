// sigma.cu — MRSF-TDDFT sigma-session (v3 ABI) on the GPU.
//
// Exports routec_sig_init / routec_sig_set_scale / routec_sig_iter /
// routec_sig_free — the session the OpenQP-side OQP_ROUTEC_SIG seam
// (routec_sig.F90) dlopens.  Given MO trial vectors it returns
// (A-B).X sigma vectors in the MO frame; the Davidson solver is untouched.
//
// Math contract: the validated CPU stub (tools/sigval/routec_sig_cpustub.cpp,
// itself G-sigma3-gated against native OpenQP MRSF to 1.7e-11 Ha).  The MO-frame
// steps (6a density factors, 6c mntoia + esum) are kept VERBATIM from the stub
// on the host; only the flop-dominant DF J/K middle runs on the device.
//
// Device J/K exploits that all 7 MRSF trial densities are LOW-RANK,
// D_m = U_m V_m^T with r <= noccb+8:
//   W_P = B_P U   Y_P = B_P V          (batched GEMM over P, B_P symmetric)
//   K   = sum_P W_P Y_P^T              (one GEMM, contraction naux*r)
//   g_P = <W_P, V>                     (small reduction kernel)
//   J   = sum_P g_P B_P                (GEMV over the naux axis)
// so a sigma application costs O(naux*nbf^2*(noccb+14)) instead of
// O(naux*nbf^3*7).  B comes from OQP_ROUTEC_B (dense or packed), resident
// on the device for the whole session.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <sys/stat.h>
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <cublas_v2.h>
#include <cblas.h>

#include "routec_df.h"   // shared resident B (routec_b_shared)

namespace {

const double ISQ2 = 1.0 / std::sqrt(2.0);

#define SIG_CUCHK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
  fprintf(stderr, "[sig-gpu] CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return false; } } while (0)
#define SIG_CBCHK(x) do { cublasStatus_t s_ = (x); if (s_ != CUBLAS_STATUS_SUCCESS) { \
  fprintf(stderr, "[sig-gpu] cuBLAS error %d at %s:%d\n", (int)s_, __FILE__, __LINE__); return false; } } while (0)

// ---- session state --------------------------------------------------------
struct Sess {
    bool   active = false;
    int    nbf = 0, nocca = 0, noccb = 0, kind = 1, naux = 0, rmax = 0;
    double scale = 0.5;
    std::vector<double> Ca, Cb, Fa, Fb;   // (nbf,nbf) row-major, a[i*nbf+j]
    cublasHandle_t hb = nullptr;
    double *dB = nullptr;                  // (naux, nbf, nbf) slices, each symmetric
    bool    b_borrowed = false;            // dB shared from df_store -- don't free
    double *dU = nullptr, *dV = nullptr;   // (nbf, rmax) col-major staging
    double *dW = nullptr, *dY = nullptr;   // (nbf, rmax) x naux batch blocks
    double *dJ = nullptr, *dK = nullptr;   // (nbf, nbf)
    double *dg = nullptr;                  // (naux)
    double *dF = nullptr;                  // (nbf, nbf) combined result
};
Sess S;

// ---- CDF compacted-on-device for the sigma product (mirrors scf.cu) --------
static bool   g_scdf_on   = false;
static const double* g_scdf_cB   = nullptr;   // fp64 arena (v2: naux x ncp; v3: n64 x ncp)
static const int*    g_scdf_keep = nullptr;
static int    g_scdf_ncp  = 0;
static int           g_scdf_lowprec = 0;      // CDF v3 mixed precision
static const void*   g_scdf_B32  = nullptr;   // float arena
static const void*   g_scdf_B16  = nullptr;   // __half arena
static const void*   g_scdf_B8   = nullptr;   // int8 arena
static const unsigned char* g_scdf_tier = nullptr;
static const int*    g_scdf_slot  = nullptr;
static const float*  g_scdf_scale = nullptr;
static double* dBtile_s = nullptr; static size_t capBtile_s = 0;
static const int SCDF_TILE = 256;
__global__ void scdf_scatter(const double* cB,const int* keep,int ncp,int nao,int nt,double* dBt){
    const long long idx=(long long)blockIdx.x*blockDim.x+threadIdx.x;
    if(idx>=(long long)nt*ncp) return;
    const int a=(int)(idx/ncp), c=(int)(idx%ncp);
    const long long t=keep[c];
    int i=(int)((sqrt(8.0*(double)t+1.0)-1.0)*0.5);
    while((long long)(i+1)*(i+2)/2<=t)++i; while((long long)i*(i+1)/2>t)--i;
    const int j=(int)(t-(long long)i*(i+1)/2);
    const double v=cB[(long long)a*ncp+c];
    const long long base=(long long)a*nao*nao;
    dBt[base+(long long)i*nao+j]=v; dBt[base+(long long)j*nao+i]=v;
}
// CDF v3: mixed-precision reconstruction, upcasting each aux column to fp64.
__global__ void scdf_scatter_mp(const double* B64,const float* B32,const __half* B16,
                                const signed char* B8,const unsigned char* tier,
                                const int* slot,const float* scale,
                                const int* keep,int ncp,int nao,int p0,int nt,double* dBt){
    const long long idx=(long long)blockIdx.x*blockDim.x+threadIdx.x;
    if(idx>=(long long)nt*ncp) return;
    const int a=(int)(idx/ncp), c=(int)(idx%ncp);
    const int P=p0+a;
    const long long t=keep[c];
    int i=(int)((sqrt(8.0*(double)t+1.0)-1.0)*0.5);
    while((long long)(i+1)*(i+2)/2<=t)++i; while((long long)i*(i+1)/2>t)--i;
    const int j=(int)(t-(long long)i*(i+1)/2);
    const int sl=slot[P]; double v;
    switch(tier[P]){                                         // 0 fp64 1 fp32 2 fp16 3 int8
        case 0:  v=B64[(long long)sl*ncp+c]; break;
        case 1:  v=(double)B32[(long long)sl*ncp+c]; break;
        case 2:  v=(double)__half2float(B16[(long long)sl*ncp+c])*(double)scale[P]; break;
        default: v=(double)B8[(long long)sl*ncp+c]*(double)scale[P]; break;
    }
    const long long base=(long long)a*nao*nao;
    dBt[base+(long long)i*nao+j]=v; dBt[base+(long long)j*nao+i]=v;
}
static bool scdf_tile(int p0,int nt,int nao){
    const long long nnl=(long long)nao*nao;
    if((size_t)nt*nnl>capBtile_s){ if(dBtile_s)cudaFree(dBtile_s);
        if(cudaMalloc(&dBtile_s,(size_t)nt*nnl*8)!=cudaSuccess) return false;
        capBtile_s=(size_t)nt*nnl; }
    if(cudaMemset(dBtile_s,0,(size_t)nt*nnl*8)!=cudaSuccess) return false;
    const long long tot=(long long)nt*g_scdf_ncp; const int TB=256;
    if(g_scdf_lowprec)
        scdf_scatter_mp<<<(unsigned)((tot+TB-1)/TB),TB>>>(
            g_scdf_cB,(const float*)g_scdf_B32,(const __half*)g_scdf_B16,(const signed char*)g_scdf_B8,
            g_scdf_tier,g_scdf_slot,g_scdf_scale,
            g_scdf_keep,g_scdf_ncp,nao,p0,nt,dBtile_s);
    else
        scdf_scatter<<<(unsigned)((tot+TB-1)/TB),TB>>>(
            g_scdf_cB+(long long)p0*g_scdf_ncp,g_scdf_keep,g_scdf_ncp,nao,nt,dBtile_s);
    return cudaGetLastError()==cudaSuccess;
}

// ---- B loader: use the shared resident B (loaded once, shared with SCF) ----
bool load_B_to_device() {
    if (getenv("OQP_CDF_ONDEV")) {                 // compacted on device + tiled sigma
        CdfDev cdf;
        if (routec_cdf_dev(&cdf)!=0) {
            fprintf(stderr,"[sig-gpu] CDF compacted load failed\n"); return false; }
        if (cdf.nao != S.nbf) { fprintf(stderr,"[sig-gpu] CDF nbf %d != %d\n",cdf.nao,S.nbf); return false; }
        S.naux=cdf.naux; g_scdf_on=true; g_scdf_cB=cdf.B64; g_scdf_keep=cdf.keep; g_scdf_ncp=cdf.ncp;
        g_scdf_lowprec=cdf.lowprec; g_scdf_B32=cdf.B32; g_scdf_B16=cdf.B16; g_scdf_B8=cdf.B8;
        g_scdf_tier=cdf.tier; g_scdf_slot=cdf.slot; g_scdf_scale=cdf.scale;
        fprintf(stderr,"[sig-gpu] B COMPACTED on device (CDF): naux=%d nbf=%d ncp=%d %s (tiled sigma)\n",
                cdf.naux,cdf.nao,cdf.ncp, cdf.lowprec?"(mixed fp64/fp32/fp16/int8)":"(fp64)");
        return true;
    }
    int naux = 0, bnbf = 0;
    const double* B = routec_b_shared(&naux, &bnbf);
    if (!B) { fprintf(stderr, "[sig-gpu] shared B load failed\n"); return false; }
    if (bnbf != S.nbf) {
        fprintf(stderr, "[sig-gpu] B nbf %d != session nbf %d\n", bnbf, S.nbf);
        return false;
    }
    S.naux = naux;
    S.dB = const_cast<double*>(B);        // borrowed from df_store -- do NOT free
    S.b_borrowed = true;
    fprintf(stderr, "[sig-gpu] B on device (shared): naux=%d nbf=%d (%.1f MB)\n",
            naux, bnbf, (size_t)naux * bnbf * bnbf * 8.0 / 1048576.0);
    return true;
}

void store_fortran(std::vector<double>& dst, const double* src, int n) {
    dst.assign((size_t)n * n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            dst[(size_t)i * n + j] = src[(size_t)j * n + i];  // F col-major -> row-major
}

void x_from_col(const double* col, std::vector<double>& X) {
    const int n = S.nbf, nocca = S.nocca, noccb = S.noccb;
    std::fill(X.begin(), X.end(), 0.0);
    size_t k = 0;
    for (int j = noccb; j < n; ++j)
        for (int i = 0; i < nocca; ++i)
            X[(size_t)i * n + j] = col[k++];
}

// ---- 6a: low-rank factors of the 7 mrsfcbc densities ----------------------
// Same math as the stub's dense mrsfcbc, kept as explicit U_m V_m^T factors.
// U/V columns are nbf-vectors appended into col-major (nbf, r) buffers.
struct LRFac { int r = 0; std::vector<double> U, V; };  // col-major nbf x r

void push_col(LRFac& f, const std::vector<double>& u, const std::vector<double>& v,
              double su, int n) {
    f.U.resize((size_t)n * (f.r + 1));
    f.V.resize((size_t)n * (f.r + 1));
    for (int i = 0; i < n; ++i) {
        f.U[(size_t)n * f.r + i] = su * u[i];
        f.V[(size_t)n * f.r + i] = v[i];
    }
    f.r++;
}

void mrsfcbc_lr(const std::vector<double>& X, LRFac* F) {
    const int n = S.nbf, nocca = S.nocca, noccb = S.noccb, mrst = S.kind;
    const int l1 = nocca - 2, l2 = nocca - 1;
    for (int m = 0; m < 7; ++m) { F[m].r = 0; F[m].U.clear(); F[m].V.clear(); }
    auto col = [&](const std::vector<double>& C, int c) {
        std::vector<double> v(n);
        for (int i = 0; i < n; ++i) v[i] = C[(size_t)i * n + c];
        return v;
    };
    // t vectors — verbatim stub loops
    std::vector<double> t_o2(n, 0.0), t_o1(n, 0.0), t_c1(n, 0.0), t_c2(n, 0.0);
    for (int a = 0; a < n; ++a) {
        double s2 = 0, s1 = 0;
        for (int j = nocca; j < n; ++j) {
            s2 += S.Cb[(size_t)a * n + j] * X[(size_t)l2 * n + j];
            s1 += S.Cb[(size_t)a * n + j] * X[(size_t)l1 * n + j];
        }
        t_o2[a] = s2; t_o1[a] = s1;
    }
    for (int a = 0; a < n; ++a) {
        double c1 = 0, c2 = 0;
        for (int i = 0; i < noccb; ++i) {
            c1 += S.Ca[(size_t)a * n + i] * X[(size_t)i * n + l1];
            c2 += S.Ca[(size_t)a * n + i] * X[(size_t)i * n + l2];
        }
        t_c1[a] = c1; t_c2[a] = c2;
    }
    std::vector<double> Ca_l1 = col(S.Ca, l1), Ca_l2 = col(S.Ca, l2);
    std::vector<double> Cb_l1 = col(S.Cb, l1), Cb_l2 = col(S.Cb, l2);
    // D0 = bo2v = Ca_l2 (x) t_o2 ; D1 = bo1v ; D2 = bco1 ; D3 = bco2
    push_col(F[0], Ca_l2, t_o2, 1.0, n);
    push_col(F[1], Ca_l1, t_o1, 1.0, n);
    push_col(F[2], t_c1, Cb_l1, 1.0, n);
    push_col(F[3], t_c2, Cb_l2, 1.0, n);
    // D4 = o21v = t_o2 (x) Ca_l1 - t_o1 (x) Ca_l2
    push_col(F[4], t_o2, Ca_l1, 1.0, n);
    push_col(F[4], t_o1, Ca_l2, -1.0, n);
    // D5 = co12 = Cb_l2 (x) t_c1 - Cb_l1 (x) t_c2
    push_col(F[5], Cb_l2, t_c1, 1.0, n);
    push_col(F[5], Cb_l1, t_c2, -1.0, n);
    // D6 = ball = D0+D1+D2+D3 + Ca[:, :noccb] tmp^T + kind terms
    push_col(F[6], Ca_l2, t_o2, 1.0, n);
    push_col(F[6], Ca_l1, t_o1, 1.0, n);
    push_col(F[6], t_c1, Cb_l1, 1.0, n);
    push_col(F[6], t_c2, Cb_l2, 1.0, n);
    {   // tmp[a,i] = sum_{j>=nocca} Cb[a,j] X[i,j]  (verbatim), one column per i
        std::vector<double> ui(n), vi(n);
        for (int i = 0; i < noccb; ++i) {
            for (int a = 0; a < n; ++a) {
                ui[a] = S.Ca[(size_t)a * n + i];
                double s = 0;
                for (int j = nocca; j < n; ++j) s += S.Cb[(size_t)a * n + j] * X[(size_t)i * n + j];
                vi[a] = s;
            }
            push_col(F[6], ui, vi, 1.0, n);
        }
    }
    const double xll = X[(size_t)l1 * n + l1];
    if (mrst == 1) {
        const double x21 = X[(size_t)l2 * n + l1], x12 = X[(size_t)l1 * n + l2];
        push_col(F[6], Ca_l2, Cb_l1, x21, n);
        push_col(F[6], Ca_l1, Cb_l2, x12, n);
        push_col(F[6], Ca_l1, Cb_l1,  xll * ISQ2, n);
        push_col(F[6], Ca_l2, Cb_l2, -xll * ISQ2, n);
    } else {
        push_col(F[6], Ca_l1, Cb_l1, xll * ISQ2, n);
        push_col(F[6], Ca_l2, Cb_l2, xll * ISQ2, n);
    }
}

// ---- 6c part 1: mrsfmntoia — VERBATIM from the CPU stub --------------------
void mrsfmntoia(const std::vector<double>* F, std::vector<double>& outX) {
    const int n = S.nbf, nocca = S.nocca, noccb = S.noccb, mrst = S.kind;
    const int l1 = nocca - 2, l2 = nocca - 1;
    const double* ado2v = F[0].data(); const double* ado1v = F[1].data();
    const double* adco1 = F[2].data(); const double* adco2 = F[3].data();
    const double* ao21v = F[4].data(); const double* aco12 = F[5].data();
    const double* agdlr = F[6].data();
    // scr = Ca^T agdlr Cb — the two nbf^3 products of the MO-frame step,
    // as BLAS GEMMs (same math as the stub's loops)
    std::vector<double> AG((size_t)n * n, 0.0);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n,
                1.0, S.Ca.data(), n, agdlr, n, 0.0, AG.data(), n);
    std::vector<double>& scr = outX;
    scr.assign((size_t)n * n, 0.0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n,
                1.0, AG.data(), n, S.Cb.data(), n, 0.0, scr.data(), n);
    std::vector<double> wrk = scr;
    auto colCa = [&](int c, int i) { return S.Ca[(size_t)i * n + c]; };
    auto colCb = [&](int c, int i) { return S.Cb[(size_t)i * n + c]; };
    {
        std::vector<double> tmp(n, 0.0);
        for (int a = 0; a < n; ++a) {
            double s = 0;
            for (int b = 0; b < n; ++b) s += ado1v[(size_t)a * n + b] * colCb(l2, b)
                                           + aco12[(size_t)a * n + b] * colCb(l1, b);
            tmp[a] = s;
        }
        for (int p = 0; p < nocca - 2; ++p) {
            double s = 0;
            for (int a = 0; a < n; ++a) s += colCa(p, a) * tmp[a];
            wrk[(size_t)p * n + l2] += s;
        }
    }
    {
        std::vector<double> tmp(n, 0.0);
        for (int a = 0; a < n; ++a) {
            double s = 0;
            for (int b = 0; b < n; ++b) s += ado2v[(size_t)a * n + b] * colCb(l1, b)
                                           - aco12[(size_t)a * n + b] * colCb(l2, b);
            tmp[a] = s;
        }
        for (int p = 0; p < nocca - 2; ++p) {
            double s = 0;
            for (int a = 0; a < n; ++a) s += colCa(p, a) * tmp[a];
            wrk[(size_t)p * n + l1] += s;
        }
    }
    {
        std::vector<double> tmp(n, 0.0);
        for (int a = 0; a < n; ++a) {
            double s = 0;
            for (int b = 0; b < n; ++b) s += adco2[(size_t)b * n + a] * colCa(l1, b)
                                           + ao21v[(size_t)b * n + a] * colCa(l2, b);
            tmp[a] = s;
        }
        for (int q = nocca; q < n; ++q) {
            double s = 0;
            for (int a = 0; a < n; ++a) s += colCb(q, a) * tmp[a];
            wrk[(size_t)l1 * n + q] += s;
        }
    }
    {
        std::vector<double> tmp(n, 0.0);
        for (int a = 0; a < n; ++a) {
            double s = 0;
            for (int b = 0; b < n; ++b) s += adco1[(size_t)b * n + a] * colCa(l2, b)
                                           - ao21v[(size_t)b * n + a] * colCa(l1, b);
            tmp[a] = s;
        }
        for (int q = nocca; q < n; ++q) {
            double s = 0;
            for (int a = 0; a < n; ++a) s += colCb(q, a) * tmp[a];
            wrk[(size_t)l2 * n + q] += s;
        }
    }
    if (mrst == 1) {
        wrk[(size_t)l1 * n + l1] = (scr[(size_t)l1 * n + l1] - scr[(size_t)l2 * n + l2]) * ISQ2;
        wrk[(size_t)l2 * n + l2] = 0.0;
    } else {
        wrk[(size_t)l1 * n + l1] = (scr[(size_t)l1 * n + l1] + scr[(size_t)l2 * n + l2]) * ISQ2;
        wrk[(size_t)l2 * n + l1] = wrk[(size_t)l1 * n + l2] = wrk[(size_t)l2 * n + l2] = 0.0;
    }
    outX.assign((size_t)nocca * (n - noccb), 0.0);
    size_t k = 0;
    for (int j = noccb; j < n; ++j)
        for (int i = 0; i < nocca; ++i)
            outX[k++] = wrk[(size_t)i * n + j];
}

// ---- 6c part 2: mrsfesum — VERBATIM from the CPU stub ----------------------
void mrsfesum(const std::vector<double>& Xin, std::vector<double>& sigflat) {
    const int n = S.nbf, nocca = S.nocca, noccb = S.noccb, mrst = S.kind;
    const int l1 = nocca - 2, l2 = nocca - 1;
    const double* fij = S.Fa.data();
    const double* fab = S.Fb.data();
    std::vector<double> scr = Xin;
    scr[(size_t)l1 * n + l1] = 0.0;
    scr[(size_t)l2 * n + l2] = 0.0;
    std::vector<double> wrk1((size_t)n * n, 0.0);
    for (int i = 0; i < nocca; ++i)
        for (int b = noccb; b < n; ++b) {
            double s = 0;
            for (int c = noccb; c < n; ++c) s += scr[(size_t)i * n + c] * fab[(size_t)b * n + c];
            double s2 = 0;
            for (int k = 0; k < nocca; ++k) s2 += fij[(size_t)i * n + k] * scr[(size_t)k * n + b];
            wrk1[(size_t)i * n + b] = s - s2;
        }
    const double xlr = Xin[(size_t)l1 * n + l1];
    double dumn = 0.0;
    if (mrst == 1) {
        for (int j = noccb; j < n; ++j) {
            wrk1[(size_t)l1 * n + j] += fab[(size_t)j * n + l1] * xlr * ISQ2;
            wrk1[(size_t)l2 * n + j] -= fab[(size_t)j * n + l2] * xlr * ISQ2;
        }
        for (int i = 0; i < nocca; ++i) {
            wrk1[(size_t)i * n + l1] -= fij[(size_t)i * n + l1] * xlr * ISQ2;
            wrk1[(size_t)i * n + l2] += fij[(size_t)i * n + l2] * xlr * ISQ2;
        }
        for (int k = 0; k < nocca; ++k)
            dumn += -fij[(size_t)l1 * n + k] * scr[(size_t)k * n + l1]
                    +  fij[(size_t)l2 * n + k] * scr[(size_t)k * n + l2];
        for (int b = noccb; b < n; ++b)
            dumn +=  fab[(size_t)l1 * n + b] * scr[(size_t)l1 * n + b]
                    - fab[(size_t)l2 * n + b] * scr[(size_t)l2 * n + b];
    } else {
        for (int j = noccb; j < n; ++j) {
            wrk1[(size_t)l1 * n + j] += fab[(size_t)j * n + l1] * xlr * ISQ2;
            wrk1[(size_t)l2 * n + j] += fab[(size_t)j * n + l2] * xlr * ISQ2;
        }
        for (int i = 0; i < nocca; ++i) {
            wrk1[(size_t)i * n + l1] -= fij[(size_t)i * n + l1] * xlr * ISQ2;
            wrk1[(size_t)i * n + l2] -= fij[(size_t)i * n + l2] * xlr * ISQ2;
        }
        for (int k = 0; k < nocca; ++k)
            dumn += -fij[(size_t)l1 * n + k] * scr[(size_t)k * n + l1]
                    -  fij[(size_t)l2 * n + k] * scr[(size_t)k * n + l2];
        for (int b = noccb; b < n; ++b)
            dumn +=  fab[(size_t)l1 * n + b] * scr[(size_t)l1 * n + b]
                    + fab[(size_t)l2 * n + b] * scr[(size_t)l2 * n + b];
    }
    wrk1[(size_t)l1 * n + l1] = dumn * ISQ2
        + xlr * (fab[(size_t)l1 * n + l1] + fab[(size_t)l2 * n + l2]
                 - fij[(size_t)l1 * n + l1] - fij[(size_t)l2 * n + l2]) * 0.5;
    if (mrst == 1) {
        wrk1[(size_t)l2 * n + l2] = 0.0;
    } else {
        wrk1[(size_t)l2 * n + l1] = wrk1[(size_t)l1 * n + l2] = wrk1[(size_t)l2 * n + l2] = 0.0;
    }
    size_t k = 0;
    for (int j = noccb; j < n; ++j)
        for (int i = 0; i < nocca; ++i)
            sigflat[k++] += wrk1[(size_t)i * n + j];
}

void free_device() {
    // S.dB is borrowed from df_store (shared with SCF) -- release the pointer
    // but do NOT free the buffer; df_store owns it.
    if (S.dB && !S.b_borrowed) cudaFree(S.dB);
    S.dB = nullptr; S.b_borrowed = false;
    if (S.dU) cudaFree(S.dU);
    if (S.dV) cudaFree(S.dV);
    if (S.dW) cudaFree(S.dW);
    if (S.dY) cudaFree(S.dY);
    if (S.dJ) cudaFree(S.dJ);
    if (S.dK) cudaFree(S.dK);
    if (S.dg) cudaFree(S.dg);
    if (S.dF) cudaFree(S.dF);
    if (S.hb) cublasDestroy(S.hb);
    S = Sess();
}

}  // namespace

// ---- device kernels --------------------------------------------------------
__global__ void sigk_gdot(long len, const double* W, const double* V, double* g) {
    // g[P] = <W_P, V>, one block per aux slice P
    extern __shared__ double sh[];
    const double* w = W + (size_t)blockIdx.x * len;
    double s = 0.0;
    for (long t = threadIdx.x; t < len; t += blockDim.x) s += w[t] * V[t];
    sh[threadIdx.x] = s;
    __syncthreads();
    for (int k = blockDim.x / 2; k; k >>= 1) {
        if (threadIdx.x < k) sh[threadIdx.x] += sh[threadIdx.x + k];
        __syncthreads();
    }
    if (threadIdx.x == 0) g[blockIdx.x] = sh[0];
}

__global__ void sigk_fcomb(long nn, double cou, double scale,
                           const double* J, const double* K, double* F) {
    long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < nn) F[t] = (cou != 0.0 ? cou * J[t] : 0.0) - scale * K[t];
}

// =========================================================== v3 ABI ========
extern "C" {

int routec_sig_init(const int* nbf_, const double* mo_a, const double* mo_b,
                    const double* fmo_a, const double* fmo_b,
                    const int* nocca_, const int* noccb_, const int* kind_) {
    if (S.active) free_device();
    const int nbf = *nbf_;
    S.nbf = nbf; S.nocca = *nocca_; S.noccb = *noccb_; S.kind = *kind_;
    S.scale = 0.5;
    store_fortran(S.Ca, mo_a, nbf);
    store_fortran(S.Cb, mo_b, nbf);
    store_fortran(S.Fa, fmo_a, nbf);
    store_fortran(S.Fb, fmo_b, nbf);
    if (!load_B_to_device()) { free_device(); return 1; }
    if (cublasCreate(&S.hb) != CUBLAS_STATUS_SUCCESS) { free_device(); return 3; }
    S.rmax = S.noccb + 16;
    const size_t nn = (size_t)nbf * nbf;
    bool ok = true;
    ok = ok && cudaMalloc(&S.dU, (size_t)nbf * S.rmax * sizeof(double)) == cudaSuccess;
    ok = ok && cudaMalloc(&S.dV, (size_t)nbf * S.rmax * sizeof(double)) == cudaSuccess;
    ok = ok && cudaMalloc(&S.dW, (size_t)S.naux * nbf * S.rmax * sizeof(double)) == cudaSuccess;
    ok = ok && cudaMalloc(&S.dY, (size_t)S.naux * nbf * S.rmax * sizeof(double)) == cudaSuccess;
    ok = ok && cudaMalloc(&S.dJ, nn * sizeof(double)) == cudaSuccess;
    ok = ok && cudaMalloc(&S.dK, nn * sizeof(double)) == cudaSuccess;
    ok = ok && cudaMalloc(&S.dg, (size_t)S.naux * sizeof(double)) == cudaSuccess;
    ok = ok && cudaMalloc(&S.dF, nn * sizeof(double)) == cudaSuccess;
    if (!ok) { fprintf(stderr, "[sig-gpu] device alloc failed\n"); free_device(); return 4; }
    S.active = true;
    fprintf(stderr, "[sig-gpu] session: nbf=%d nocca=%d noccb=%d kind=%d ntrial=%d naux=%d\n",
            nbf, S.nocca, S.noccb, S.kind, S.nocca * (nbf - S.noccb), S.naux);
    return 0;
}

void routec_sig_set_scale(const double* s) { S.scale = *s; }

void routec_sig_iter(const double* bvec_mo, const int* nv_new_,
                     double* sigma_mo, int* info) {
    *info = 1;
    if (!S.active) { fprintf(stderr, "[sig-gpu] not init\n"); return; }
    const int n = S.nbf, nocca = S.nocca, noccb = S.noccb, naux = S.naux;
    const int nv = *nv_new_;
    const int ntrial = nocca * (n - noccb);
    const double scale = S.scale;
    const int mrst = S.kind;
    const size_t nn = (size_t)n * n;
    const double one = 1.0, zero = 0.0;

    std::vector<double> X(nn);
    LRFac lr[7];
    std::vector<double> Fh[7], Fcol(nn);
    for (int m = 0; m < 7; ++m) Fh[m].assign(nn, 0.0);
    std::vector<double> sflat;

    for (int c = 0; c < nv; ++c) {
        const double* col = bvec_mo + (size_t)c * ntrial;
        x_from_col(col, X);
        mrsfcbc_lr(X, lr);
        for (int m = 0; m < 7; ++m) {
            const int r = lr[m].r;
            if (r > S.rmax) { fprintf(stderr, "[sig-gpu] r=%d > rmax=%d\n", r, S.rmax); return; }
            const double cou = (m < 4) ? scale : 0.0;
            // stage factors, W_P = B_P U, Y_P = B_P V (col-major batches)
            if (cudaMemcpy(S.dU, lr[m].U.data(), (size_t)n * r * sizeof(double),
                           cudaMemcpyHostToDevice) != cudaSuccess) return;
            if (cudaMemcpy(S.dV, lr[m].V.data(), (size_t)n * r * sizeof(double),
                           cudaMemcpyHostToDevice) != cudaSuccess) return;
          if (g_scdf_on) {
            // ---- tiled sigma from compacted B (no full dense B on device) ----
            if (cudaMemset(S.dK, 0, nn*sizeof(double)) != cudaSuccess) return;
            for (int p0=0;p0<naux;p0+=SCDF_TILE){ const int nt=(naux-p0<SCDF_TILE)?(naux-p0):SCDF_TILE;
                if(!scdf_tile(p0,nt,n)) return;
                if (cublasDgemmStridedBatched(S.hb,CUBLAS_OP_N,CUBLAS_OP_N, n,r,n,&one,
                    dBtile_s,n,nn, S.dU,n,0, &zero, S.dW,n,(long long)n*r, nt)!=CUBLAS_STATUS_SUCCESS) return;
                if (cublasDgemmStridedBatched(S.hb,CUBLAS_OP_N,CUBLAS_OP_N, n,r,n,&one,
                    dBtile_s,n,nn, S.dV,n,0, &zero, S.dY,n,(long long)n*r, nt)!=CUBLAS_STATUS_SUCCESS) return;
                const double bta=(p0==0)?0.0:1.0;
                if (cublasDgemm(S.hb,CUBLAS_OP_N,CUBLAS_OP_T, n,n, nt*r, &one,
                    S.dW,n, S.dY,n, &bta, S.dK,n)!=CUBLAS_STATUS_SUCCESS) return;    // K += Σ W_P Y_P^T
                if (cou != 0.0)                                                     // g[p0..] = <W_P,V>
                    sigk_gdot<<<nt,256,256*sizeof(double)>>>((long)n*r, S.dW, S.dV, S.dg+p0);
            }
            if (cou != 0.0) {                                                        // J = Σ_P g_P B_P (tiled)
                if (cudaMemset(S.dJ, 0, nn*sizeof(double)) != cudaSuccess) return;
                for (int p0=0;p0<naux;p0+=SCDF_TILE){ const int nt=(naux-p0<SCDF_TILE)?(naux-p0):SCDF_TILE;
                    if(!scdf_tile(p0,nt,n)) return;
                    const double bta=(p0==0)?0.0:1.0;
                    if (cublasDgemv(S.hb,CUBLAS_OP_N, (int)nn, nt, &one, dBtile_s,(int)nn,
                        S.dg+p0, 1, &bta, S.dJ, 1)!=CUBLAS_STATUS_SUCCESS) return;
                }
            }
          } else {
            if (cublasDgemmStridedBatched(S.hb, CUBLAS_OP_N, CUBLAS_OP_N,
                    n, r, n, &one, S.dB, n, nn, S.dU, n, 0, &zero,
                    S.dW, n, (long long)n * r, naux) != CUBLAS_STATUS_SUCCESS) return;
            if (cublasDgemmStridedBatched(S.hb, CUBLAS_OP_N, CUBLAS_OP_N,
                    n, r, n, &one, S.dB, n, nn, S.dV, n, 0, &zero,
                    S.dY, n, (long long)n * r, naux) != CUBLAS_STATUS_SUCCESS) return;
            // K = sum_P W_P Y_P^T  (flat contraction over naux*r)
            if (cublasDgemm(S.hb, CUBLAS_OP_N, CUBLAS_OP_T, n, n, naux * r,
                    &one, S.dW, n, S.dY, n, &zero, S.dK, n) != CUBLAS_STATUS_SUCCESS) return;
            if (cou != 0.0) {
                // g_P = <W_P, V> ; J = sum_P g_P B_P
                sigk_gdot<<<naux, 256, 256 * sizeof(double)>>>((long)n * r, S.dW, S.dV, S.dg);
                if (cublasDgemv(S.hb, CUBLAS_OP_N, (int)nn, naux, &one, S.dB, (int)nn,
                        S.dg, 1, &zero, S.dJ, 1) != CUBLAS_STATUS_SUCCESS) return;
            }
          }
            sigk_fcomb<<<(unsigned)((nn + 255) / 256), 256>>>((long)nn, cou, scale, S.dJ, S.dK, S.dF);
            if (cudaMemcpy(Fcol.data(), S.dF, nn * sizeof(double),
                           cudaMemcpyDeviceToHost) != cudaSuccess) return;
            // device buffers are col-major; host 6c code is row-major
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    Fh[m][(size_t)i * n + j] = Fcol[(size_t)j * n + i];
        }
        if (mrst == 3)
            for (int m = 0; m < 6; ++m)
                for (size_t t = 0; t < nn; ++t) Fh[m][t] = -Fh[m][t];
        mrsfmntoia(Fh, sflat);
        mrsfesum(X, sflat);
        std::memcpy(sigma_mo + (size_t)c * ntrial, sflat.data(), ntrial * sizeof(double));
    }
    *info = 0;
}

void routec_sig_free(void) {
    free_device();
}

}  // extern "C"
