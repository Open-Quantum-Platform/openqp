// routec_sig_cpustub.cpp — CPU validation stub for the R3 sigma-session v3 ABI.
//
// Implements routec_sig_init/set_scale/iter/free with the EXACT Stage-0 MRSF
// sigma pipeline (verbatim port of sig_ref.py build_sigma_stage0 — itself the
// validated == device referee, G-sigma0 rel 2e-15). DF-B J/K from the same
// B file the GPU engine reads (OQP_ROUTEC_B). No CUDA, no GPU.
//
// Purpose: prove the Fortran seam plumbing in tdhf_mrsf_energy.F90 — the
// bvec_mo <-> sigma_mo (ntrial, nv) col-major layout, the fa/fb MO-Fock
// handoff, the column slicing, the kind/scale handoff, and Davidson subspace
// consistency — WITHOUT involving the GPU engine. When the driver is pointed
// at this stub via OQP_ROUTEC_SIG, the amo it produces must match a standalone
// Stage-0 computation on the same (bvec, Ca, Cb, fa, fb, B, scale, kind) to
// machine precision (G-sigma1 plumbing gate).
//
// Build:
//   g++ -O2 -fPIC -shared -fopenmp routec_sig_cpustub.cpp -o libroutec_sig_cpustub.so
//
// All scalars passed by pointer (matches the .cu ABI and the Fortran bind(C)
// interface in routec_sig.F90).

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <sys/stat.h>

namespace {

const double ISQ2 = 1.0 / std::sqrt(2.0);

// ---- session state -------------------------------------------------------
struct Sig {
    bool   active = false;
    int    nbf = 0, nocca = 0, noccb = 0, kind = 1;
    double scale = 0.5;
    std::vector<double> Ca, Cb, Fa, Fb;   // (nbf,nbf) row-major here (a[i*nbf+j])
};
Sig g;

// ---- B (DF aux slices), dense (naux, nbf, nbf), B_P symmetric ------------
int g_naux = 0, g_Bnbf = 0;
std::vector<double> g_B;   // naux*nbf*nbf, row-major per slice

bool load_B() {
    if (!g_B.empty()) return true;
    const char* path = getenv("OQP_ROUTEC_B");
    if (!path) { fprintf(stderr, "[sig-cpustub] OQP_ROUTEC_B unset\n"); return false; }
    FILE* f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "[sig-cpustub] cannot open %s\n", path); return false; }
    struct stat st;
    if (fstat(fileno(f), &st) != 0) { fclose(f); return false; }
    const size_t fbytes = (size_t)st.st_size;
    int hdr[2];
    if (fread(hdr, sizeof(int), 2, f) != 2) { fclose(f); return false; }
    bool packed = false;
    if (hdr[0] > 0 && hdr[1] > 0) { g_naux = hdr[0]; g_Bnbf = hdr[1]; }
    else if (hdr[0] < 0 && hdr[1] > 0) { packed = true; g_Bnbf = -hdr[0]; g_naux = hdr[1]; }
    else if (hdr[0] > 0 && hdr[1] < 0) { packed = true; g_naux = hdr[0]; g_Bnbf = -hdr[1]; }
    else { fprintf(stderr, "[sig-cpustub] bad B header (%d,%d)\n", hdr[0], hdr[1]); fclose(f); return false; }
    const size_t nn = (size_t)g_Bnbf * g_Bnbf;
    const size_t npair = (size_t)g_Bnbf * (g_Bnbf + 1) / 2;
    const size_t nread = packed ? (size_t)g_naux * npair : (size_t)g_naux * nn;
    if (fbytes != 2 * sizeof(int) + nread * sizeof(double)) {
        fprintf(stderr, "[sig-cpustub] B size mismatch\n"); fclose(f); return false;
    }
    std::vector<double> raw(nread);
    if (fread(raw.data(), sizeof(double), nread, f) != nread) { fclose(f); return false; }
    fclose(f);
    g_B.assign((size_t)g_naux * nn, 0.0);
    if (packed) {  // unpack lower/upper-tri rows -> dense symmetric
        for (int p = 0; p < g_naux; ++p) {
            const double* s = raw.data() + (size_t)p * npair;
            double* d = g_B.data() + (size_t)p * nn;
            size_t k = 0;
            for (int i = 0; i < g_Bnbf; ++i)
                for (int j = 0; j <= i; ++j) {
                    d[(size_t)i * g_Bnbf + j] = s[k];
                    d[(size_t)j * g_Bnbf + i] = s[k];
                    ++k;
                }
        }
    } else {
        std::memcpy(g_B.data(), raw.data(), nread * sizeof(double));
    }
    fprintf(stderr, "[sig-cpustub] B loaded: naux=%d nbf=%d (%s)\n",
            g_naux, g_Bnbf, packed ? "packed" : "dense");
    return true;
}

// ---- small dense helpers (all matrices nbf*nbf row-major, a[i*n+j]) ------
inline double Cael(int i, int j) { return g.Ca[(size_t)i * g.nbf + j]; }  // Ca[i,j]
inline double Cbel(int i, int j) { return g.Cb[(size_t)i * g.nbf + j]; }

// RAW J(M), K(M) from DF-B, M arbitrary (nbf,nbf). Matches jk_dfb in sig_ref.py:
//   g_P = B_P : M ;  J = sum_P g_P B_P ;  K = sum_P B_P M B_P.
void jk_dfb(const double* M, double* J, double* K) {
    const int n = g.nbf;
    const size_t nn = (size_t)n * n;
    std::fill(J, J + nn, 0.0);
    std::fill(K, K + nn, 0.0);
    std::vector<double> BM((size_t)n * n);
    for (int p = 0; p < g_naux; ++p) {
        const double* Bp = g_B.data() + (size_t)p * nn;
        // gP = Bp : M
        double gp = 0.0;
        for (size_t t = 0; t < nn; ++t) gp += Bp[t] * M[t];
        // J += gP * Bp
        for (size_t t = 0; t < nn; ++t) J[t] += gp * Bp[t];
        // BM = Bp * M ; K += BM * Bp
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                double s = 0.0;
                for (int k = 0; k < n; ++k) s += Bp[(size_t)i * n + k] * M[(size_t)k * n + j];
                BM[(size_t)i * n + j] = s;
            }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                double s = 0.0;
                for (int k = 0; k < n; ++k) s += BM[(size_t)i * n + k] * Bp[(size_t)k * n + j];
                K[(size_t)i * n + j] += s;
            }
    }
}

// outer-product accumulate: A[i,j] += s * u[i] * v[j]
inline void outer_acc(std::vector<double>& A, const std::vector<double>& u,
                      const std::vector<double>& v, double s, int n) {
    for (int i = 0; i < n; ++i) {
        double su = s * u[i];
        double* Ai = A.data() + (size_t)i * n;
        for (int j = 0; j < n; ++j) Ai[j] += su * v[j];
    }
}

// X(nbf,nbf) from a flat trial column (iatogen): active block [0:nocca, noccb:nbf],
// col-major i-fast over (i=0..nocca-1, j=noccb..nbf-1).
void x_from_col(const double* col, std::vector<double>& X) {
    const int n = g.nbf, nocca = g.nocca, noccb = g.noccb;
    std::fill(X.begin(), X.end(), 0.0);
    size_t k = 0;
    for (int j = noccb; j < n; ++j)
        for (int i = 0; i < nocca; ++i)
            X[(size_t)i * n + j] = col[k++];
}

// 6a: mrsfcbc -> 7 AO densities D[m] (nbf,nbf). Verbatim sig_ref.py mrsfcbc.
void mrsfcbc(const std::vector<double>& X, std::vector<double>* D) {
    const int n = g.nbf, nocca = g.nocca, noccb = g.noccb, mrst = g.kind;
    const int l1 = nocca - 2, l2 = nocca - 1;
    const size_t nn = (size_t)n * n;
    for (int m = 0; m < 7; ++m) std::fill(D[m].begin(), D[m].end(), 0.0);
    auto& bo2v = D[0]; auto& bo1v = D[1]; auto& bco1 = D[2]; auto& bco2 = D[3];
    auto& o21v = D[4]; auto& co12 = D[5]; auto& ball = D[6];
    // column views of Ca/Cb: Ca[:,c] = Ca[i*n+c]
    auto col = [&](const std::vector<double>& C, int c) {
        std::vector<double> v(n);
        for (int i = 0; i < n; ++i) v[i] = C[(size_t)i * n + c];
        return v;
    };
    // t_o2 = Cb[:,nocca:] @ X[l2, nocca:]
    std::vector<double> t_o2(n, 0.0), t_o1(n, 0.0), t_c1(n, 0.0), t_c2(n, 0.0);
    for (int a = 0; a < n; ++a) {
        double s2 = 0, s1 = 0;
        for (int j = nocca; j < n; ++j) {
            s2 += g.Cb[(size_t)a * n + j] * X[(size_t)l2 * n + j];
            s1 += g.Cb[(size_t)a * n + j] * X[(size_t)l1 * n + j];
        }
        t_o2[a] = s2; t_o1[a] = s1;
    }
    // t_c1 = Ca[:, :noccb] @ X[:noccb, l1] ; t_c2 = Ca[:, :noccb] @ X[:noccb, l2]
    for (int a = 0; a < n; ++a) {
        double c1 = 0, c2 = 0;
        for (int i = 0; i < noccb; ++i) {
            c1 += g.Ca[(size_t)a * n + i] * X[(size_t)i * n + l1];
            c2 += g.Ca[(size_t)a * n + i] * X[(size_t)i * n + l2];
        }
        t_c1[a] = c1; t_c2[a] = c2;
    }
    std::vector<double> Ca_l1 = col(g.Ca, l1), Ca_l2 = col(g.Ca, l2);
    std::vector<double> Cb_l1 = col(g.Cb, l1), Cb_l2 = col(g.Cb, l2);
    outer_acc(bo2v, Ca_l2, t_o2, 1.0, n);
    outer_acc(bo1v, Ca_l1, t_o1, 1.0, n);
    outer_acc(bco1, t_c1, Cb_l1, 1.0, n);
    outer_acc(bco2, t_c2, Cb_l2, 1.0, n);
    outer_acc(o21v, t_o2, Ca_l1, 1.0, n);
    outer_acc(o21v, t_o1, Ca_l2, -1.0, n);
    outer_acc(co12, Cb_l2, t_c1, 1.0, n);
    outer_acc(co12, Cb_l1, t_c2, -1.0, n);
    for (size_t t = 0; t < nn; ++t) ball[t] = bo2v[t] + bo1v[t] + bco1[t] + bco2[t];
    // ball += Ca[:, :noccb] @ (Cb[:, nocca:] @ X[:noccb, nocca:].T).T
    //   tmp(nbf, noccb): tmp[a,i] = sum_{j>=nocca} Cb[a,j] * X[i,j]
    std::vector<double> tmp((size_t)n * noccb, 0.0);
    for (int a = 0; a < n; ++a)
        for (int i = 0; i < noccb; ++i) {
            double s = 0;
            for (int j = nocca; j < n; ++j) s += g.Cb[(size_t)a * n + j] * X[(size_t)i * n + j];
            tmp[(size_t)a * noccb + i] = s;
        }
    //   ball[a,b] += sum_i Ca[a,i] * tmp[b,i]    (= Ca_o @ tmp.T)
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
            double s = 0;
            for (int i = 0; i < noccb; ++i) s += g.Ca[(size_t)a * n + i] * tmp[(size_t)b * noccb + i];
            ball[(size_t)a * n + b] += s;
        }
    const double xll = X[(size_t)l1 * n + l1];
    if (mrst == 1) {
        const double x21 = X[(size_t)l2 * n + l1], x12 = X[(size_t)l1 * n + l2];
        outer_acc(ball, Ca_l2, Cb_l1, x21, n);
        outer_acc(ball, Ca_l1, Cb_l2, x12, n);
        outer_acc(ball, Ca_l1, Cb_l1,  xll * ISQ2, n);
        outer_acc(ball, Ca_l2, Cb_l2, -xll * ISQ2, n);
    } else { // triplet
        outer_acc(ball, Ca_l1, Cb_l1, xll * ISQ2, n);
        outer_acc(ball, Ca_l2, Cb_l2, xll * ISQ2, n);
    }
}

// 6c part 1: mrsfmntoia(F[7]) -> sigma block, returns flat (ntrial). Verbatim.
void mrsfmntoia(const std::vector<double>* F, std::vector<double>& outX) {
    const int n = g.nbf, nocca = g.nocca, noccb = g.noccb, mrst = g.kind;
    const int l1 = nocca - 2, l2 = nocca - 1;
    const double* ado2v = F[0].data(); const double* ado1v = F[1].data();
    const double* adco1 = F[2].data(); const double* adco2 = F[3].data();
    const double* ao21v = F[4].data(); const double* aco12 = F[5].data();
    const double* agdlr = F[6].data();
    // scr = Ca.T @ agdlr @ Cb  (nbf,nbf): scr[p,q] = sum_{a,b} Ca[a,p] agdlr[a,b] Cb[b,q]
    std::vector<double> AG((size_t)n * n, 0.0);   // AG = Ca.T @ agdlr : AG[p,b]
    for (int p = 0; p < n; ++p)
        for (int b = 0; b < n; ++b) {
            double s = 0;
            for (int a = 0; a < n; ++a) s += g.Ca[(size_t)a * n + p] * agdlr[(size_t)a * n + b];
            AG[(size_t)p * n + b] = s;
        }
    std::vector<double>& scr = outX;  // reuse; scr[p,q]
    scr.assign((size_t)n * n, 0.0);
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) {
            double s = 0;
            for (int b = 0; b < n; ++b) s += AG[(size_t)p * n + b] * g.Cb[(size_t)b * n + q];
            scr[(size_t)p * n + q] = s;
        }
    std::vector<double> wrk = scr;  // copy
    // column views
    auto colCa = [&](int c, int i) { return g.Ca[(size_t)i * n + c]; };
    auto colCb = [&](int c, int i) { return g.Cb[(size_t)i * n + c]; };
    // tmp = ado1v @ Cb[:,l2] + aco12 @ Cb[:,l1]; wrk[:nocca-2,l2] += Ca[:,:nocca-2].T @ tmp
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
    // tmp = adco2.T @ Ca[:,l1] + ao21v.T @ Ca[:,l2]; wrk[l1, nocca:] += Cb[:,nocca:].T @ tmp
    {
        std::vector<double> tmp(n, 0.0);
        for (int a = 0; a < n; ++a) {  // (.T) -> sum over first index
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
    // flatten wrk[:nocca, noccb:] in F-order (i fast)
    outX.assign((size_t)nocca * (n - noccb), 0.0);
    size_t k = 0;
    for (int j = noccb; j < n; ++j)
        for (int i = 0; i < nocca; ++i)
            outX[k++] = wrk[(size_t)i * n + j];
}

// 6c part 2: mrsfesum(X, fij=Fa, fab=Fb) -> add to flat sigma. Verbatim.
void mrsfesum(const std::vector<double>& Xin, std::vector<double>& sigflat) {
    const int n = g.nbf, nocca = g.nocca, noccb = g.noccb, mrst = g.kind;
    const int l1 = nocca - 2, l2 = nocca - 1;
    const double* fij = g.Fa.data();   // MO Fock alpha (occ-occ used)
    const double* fab = g.Fb.data();   // MO Fock beta  (virt-virt used)
    std::vector<double> scr = Xin;     // copy
    scr[(size_t)l1 * n + l1] = 0.0;
    scr[(size_t)l2 * n + l2] = 0.0;
    std::vector<double> wrk1((size_t)n * n, 0.0);
    // wrk1[:nocca, noccb:] = scr[:nocca, noccb:] @ fab[noccb:, noccb:].T
    //                      - fij[:nocca, :nocca] @ scr[:nocca, noccb:]
    for (int i = 0; i < nocca; ++i)
        for (int b = noccb; b < n; ++b) {
            double s = 0;
            for (int c = noccb; c < n; ++c) s += scr[(size_t)i * n + c] * fab[(size_t)b * n + c]; // fab[b,c] = fab[noccb:,noccb:].T
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

}  // namespace

// =========================================================== v3 ABI ========
extern "C" {

// NOTE: mo_a/mo_b/fmo_a/fmo_b arrive as Fortran col-major (nbf,nbf). We store
// them transposed into row-major so that a[i*nbf+j] == A(i,j) (1-based i,j-1).
// Fortran A(i,j) at flat (j-1)*nbf+(i-1); we set ours[i*nbf+j] = inA[j*nbf+i].
static void store_fortran(std::vector<double>& dst, const double* src, int n) {
    dst.assign((size_t)n * n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            dst[(size_t)i * n + j] = src[(size_t)j * n + i];  // transpose F->row-major
}

int routec_sig_init(const int* nbf_, const double* mo_a, const double* mo_b,
                    const double* fmo_a, const double* fmo_b,
                    const int* nocca_, const int* noccb_, const int* kind_) {
    if (!load_B()) return 1;
    const int nbf = *nbf_;
    if (nbf != g_Bnbf) {
        fprintf(stderr, "[sig-cpustub] init nbf mismatch: %d vs B %d\n", nbf, g_Bnbf);
        return 2;
    }
    g.nbf = nbf; g.nocca = *nocca_; g.noccb = *noccb_; g.kind = *kind_;
    g.scale = 0.5;
    store_fortran(g.Ca, mo_a, nbf);
    store_fortran(g.Cb, mo_b, nbf);
    store_fortran(g.Fa, fmo_a, nbf);
    store_fortran(g.Fb, fmo_b, nbf);
    g.active = true;
    fprintf(stderr, "[sig-cpustub] session: nbf=%d nocca=%d noccb=%d kind=%d ntrial=%d\n",
            nbf, g.nocca, g.noccb, g.kind, g.nocca * (nbf - g.noccb));
    return 0;
}

void routec_sig_set_scale(const double* s) { g.scale = *s; }

void routec_sig_iter(const double* bvec_mo, const int* nv_new_,
                     double* sigma_mo, int* info) {
    *info = 1;
    if (!g.active) { fprintf(stderr, "[sig-cpustub] not init\n"); return; }
    const int n = g.nbf, nocca = g.nocca, noccb = g.noccb;
    const int nv = *nv_new_;
    const int ntrial = nocca * (n - noccb);
    const double scale = g.scale;
    const int mrst = g.kind;
    std::vector<double> X((size_t)n * n);
    std::vector<double> D[7], F[7];
    for (int m = 0; m < 7; ++m) { D[m].assign((size_t)n * n, 0.0); F[m].assign((size_t)n * n, 0.0); }
    std::vector<double> J((size_t)n * n), K((size_t)n * n);
    std::vector<double> sflat;
    for (int c = 0; c < nv; ++c) {
        const double* col = bvec_mo + (size_t)c * ntrial;
        x_from_col(col, X);
        mrsfcbc(X, D);
        for (int m = 0; m < 7; ++m) {
            jk_dfb(D[m].data(), J.data(), K.data());
            const double cou = (m < 4) ? scale : 0.0;
            for (size_t t = 0; t < (size_t)n * n; ++t)
                F[m][t] = cou * J[t] - scale * K[t];
        }
        if (mrst == 3)
            for (int m = 0; m < 6; ++m)
                for (size_t t = 0; t < (size_t)n * n; ++t) F[m][t] = -F[m][t];
        mrsfmntoia(F, sflat);            // writes sflat (ntrial)
        mrsfesum(X, sflat);              // adds esum term
        std::memcpy(sigma_mo + (size_t)c * ntrial, sflat.data(), ntrial * sizeof(double));
    }
    *info = 0;
}

void routec_sig_free(void) {
    g = Sig();
}

}  // extern "C"
