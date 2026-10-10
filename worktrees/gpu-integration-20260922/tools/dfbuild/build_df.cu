// Fork B.2b: matched 3-center head-to-head vs gpu4pyscf get_int3c2e.
// Same workload as calib/gpu4pyscf_int3c2e_20260611.py: (H2O)8 cc-pVDZ /
// def2-universal-jkfit, full (nao,nao,naux) FP64 tensor.
//
// Pipeline: orbital shells -> canonical pairs (l_hi>=l_lo) with per-primitive-
// pair channel tensors (coefficients folded); aux shells -> dummy-s pairs
// (exponent-0 partner). Triplets (pair, aux) Schwarz-screened at 1e-13
// (conservative: keep more than gpu4pyscf does), grouped by
// (pairclass, auxclass, kprim), launched chunked through the generated
// scalar+block kernels, scattered into the dense host tensor, written out.
// Timing: sum of GPU event times over all launches (kernel-only analog) and
// wall time of the whole build (setup analog).
//
//   build: nvcc -O3 -arch=sm_80 -std=c++17 -I include -diag-suppress 177 \
//            -x cu test/bench_routec_int3c.cu -o bench_routec_int3c
//   run:   ./bench_routec_int3c tools/routec_tables.bin shells_w8.txt aux_w8.txt out_int3c.bin
#define main routec_ref_main_disabled
#include "test_routec.cpp"
#undef main

#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <algorithm>
#include "routec_gen.inc"
#include "routec_gen_cuda.cuh"

#include <cublas_v2.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>

// LAPACK for the native (pyscf-free) auxiliary-metric Cholesky inverse.
extern "C" {
  void dpotrf_(const char*,const int*,double*,const int*,int*);
  void dtrtri_(const char*,const char*,const int*,double*,const int*,int*);
}

// in-process solve handoff (linked from libopenqp_gpu)
extern "C" {
  void routec_b_adopt(double* dB_device, int naux, int nbf);
  void routec_scf_solve(const double* h, const double* s, const int* nbf,
                        const int* nocc, const double* enuc,
                        const double* sc, const double* se,
                        const double* conv_e, const double* conv_d,
                        const int* maxit,
                        double* e_out, double* d_out,
                        double* c_out, double* eps_out,
                        int* ncyc_out, int* info);
}

// unpack the whitened dB (col-major npair x naux, dB[t + npair*P], lower-tri
// pair t) into a FULL dense (naux, nao, nao) tensor for adoption by df_store.
__global__ void unpack_b_full_kernel(double* __restrict__ full,
                                     const double* __restrict__ B,
                                     const int* __restrict__ iu,
                                     const int* __restrict__ ju,
                                     long npair, int naux, int nao){
  long t = (long)blockIdx.x*blockDim.x + threadIdx.x;
  if(t >= npair*(long)naux) return;
  long pr = t % npair; int P = (int)(t / npair);
  int i = iu[pr], j = ju[pr];
  double v = B[pr + npair*(long)P];
  full[(size_t)P*nao*nao + (size_t)i*nao + j] = v;
  full[(size_t)P*nao*nao + (size_t)j*nao + i] = v;
}

// ---- native pyscf-free aux Coulomb metric V_PQ=(P|Q) via the SAME routec
// kernel that builds the 3-center aux index (each aux shell paired with a
// dummy s, exactly as the aux-Schwarz bound at line ~268).  Because V and the
// 3-center M come from one kernel in one convention, the fit is exact with the
// aux transform set to identity (Cag=I, aperm=identity, ascale=1); the aux
// index is summed in B.B^T so its convention cancels.  Returns Linv=inv(chol(V))
// (lower), row-major (naux,naux) -- same object prep's pyscf path produced.
template<class ShellT>
static void native_aux_linv(const RTC& T,int mode,const std::vector<ShellT>& ax,
                            int naux,std::vector<double>& Linv){
  std::vector<double> V((size_t)naux*naux,0.0);
  const int nsh=(int)ax.size();
  #pragma omp parallel for schedule(dynamic)
  for(int p=0;p<nsh;++p){
    const ShellT& P=ax[p];
    std::vector<double> blk;
    for(int q=0;q<nsh;++q){
      const ShellT& Q=ax[q];
      const int npr=P.ns, nqr=Q.ns;
      std::vector<double> acc((size_t)npr*nqr,0.0);
      for(int pa=0;pa<P.np;++pa)
        for(int qc=0;qc<Q.np;++qc){
          routec_block(T,mode,P.l,P.X,P.e[pa],0,P.X,0.0,
                              Q.l,Q.X,Q.e[qc],0,Q.X,0.0,blk,false);
          const double w=P.c[pa]*Q.c[qc];
          for(size_t z=0;z<acc.size();++z) acc[z]+=w*blk[z];
        }
      for(int mi=0;mi<npr;++mi)
        for(int mj=0;mj<nqr;++mj)
          V[(size_t)(P.ao0+mi)*naux+(Q.ao0+mj)]=acc[(size_t)mi*nqr+mj];
    }
  }
  // V is row-major symmetric.  LAPACK is column-major, so "U" (col-major upper)
  // operates on the row-major LOWER triangle: dpotrf_("U") yields L (L L^T = V)
  // in the row-major lower triangle; dtrtri_("U") inverts it to Linv = L^{-1}.
  // Zero the strict row-major upper so the object matches pyscf's lower Linv.
  int n=naux,info=0;
  dpotrf_("U",&n,V.data(),&n,&info);
  if(info!=0){ fprintf(stderr,"native metric dpotrf info=%d\n",info); exit(1); }
  dtrtri_("U","N",&n,V.data(),&n,&info);
  if(info!=0){ fprintf(stderr,"native metric dtrtri info=%d\n",info); exit(1); }
  for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) V[(size_t)i*n+j]=0.0;
  Linv.assign(V.begin(),V.end());
}

// orbital perm/scale + lower-tri pack, fused.  dMpk[t,P] = s_i s_j * dM[(pi*nao+pj),P]
// where (i,j)=lower-tri index t (i>=j), pi=perm[i], pj=perm[j],
// s_i=dscale[i], s_j=dscale[j].  One thread per (t,P): coalesced along P.
__global__ void pack_kernel(double* __restrict__ Mpk, const double* __restrict__ M,
                            const int* __restrict__ perm, const double* __restrict__ ds,
                            const int* __restrict__ iu, const int* __restrict__ ju,
                            int npair, int naux, int nao){
  long t = (long)blockIdx.x*blockDim.x + threadIdx.x;
  if(t >= (long)npair*naux) return;
  int pidx = (int)(t/naux); int P = (int)(t%naux);
  int i=iu[pidx], j=ju[pidx];
  int pi=perm[i], pj=perm[j];
  double s=ds[i]*ds[j];
  Mpk[(long)pidx*naux+P] = s * M[((long)pi*nao+pj)*naux+P];
}
static void extern_pack(double* Mpk,const double* M,const int* perm,const double* ds,
                        const int* iu,const int* ju,int npair,int naux,int nao){
  long tot=(long)npair*naux;
  pack_kernel<<<(unsigned)((tot+255)/256),256>>>(Mpk,M,perm,ds,iu,ju,npair,naux,nao);
}

// pack FP64 dense T3 (nao,nao,naux) -> FP32 packed (npair,naux) with orbital
// perm/scale, BEFORE the aux-conv (commutes; halves the aux-conv GEMM size).
__global__ void pack_t3_f32_kernel(float* __restrict__ T3pk, const double* __restrict__ T3,
                            const int* __restrict__ perm, const double* __restrict__ ds,
                            const int* __restrict__ iu, const int* __restrict__ ju,
                            int npair, int naux, int nao){
  long t = (long)blockIdx.x*blockDim.x + threadIdx.x;
  if(t >= (long)npair*naux) return;
  int pidx = (int)(t/naux); int P = (int)(t%naux);
  int i=iu[pidx], j=ju[pidx];
  int pi=perm[i], pj=perm[j];
  double s=ds[i]*ds[j];
  T3pk[(long)pidx*naux+P] = (float)(s * T3[((long)pi*nao+pj)*naux+P]);
}
static void extern_pack_t3_f32(float* T3pk,const double* T3,const int* perm,const double* ds,
                        const int* iu,const int* ju,int npair,int naux,int nao){
  long tot=(long)npair*naux;
  pack_t3_f32_kernel<<<(unsigned)((tot+255)/256),256>>>(T3pk,T3,perm,ds,iu,ju,npair,naux,nao);
}
__global__ void cast_f32_f64_kernel(double* o,const float* in,long n){
  long t=(long)blockIdx.x*blockDim.x+threadIdx.x; if(t<n) o[t]=(double)in[t]; }
static void extern_cast_f32_f64(double* o,const float* in,long n){
  cast_f32_f64_kernel<<<(unsigned)((n+255)/256),256>>>(o,in,n); }

// FULL convention pack in one kernel: orbital perm/scale (axes i,j) + aux perm/scale
// (axis P) + lower-tri pack, reading dense T3 directly. Eliminates the dense
// block-diagonal aux-conv GEMM (which is ~99% zeros). Mathematically exact: Cag and
// Co are both blockdiag inv(M[l]) = scalar*permutation per shell (verified err 0).
//   Mpk[t,P] = ds[i]*ds[j]*ascale[P] * T3[perm[i], perm[j], aperm[P]]
__global__ void pack_full_kernel(double* __restrict__ Mpk, const double* __restrict__ T3,
                            const int* __restrict__ perm, const double* __restrict__ ds,
                            const int* __restrict__ aperm, const double* __restrict__ as,
                            const int* __restrict__ iu, const int* __restrict__ ju,
                            int npair, int naux, int nao){
  long t = (long)blockIdx.x*blockDim.x + threadIdx.x;
  if(t >= (long)npair*naux) return;
  int pidx = (int)(t/naux); int P = (int)(t%naux);
  int i=iu[pidx], j=ju[pidx];
  int pi=perm[i], pj=perm[j], pp=aperm[P];
  double s=ds[i]*ds[j]*as[P];
  Mpk[(long)pidx*naux+P] = s * T3[((long)pi*nao+pj)*naux+pp];
}
static void extern_pack_full(double* Mpk,const double* T3,const int* perm,const double* ds,
                        const int* aperm,const double* as,const int* iu,const int* ju,
                        int npair,int naux,int nao){
  long tot=(long)npair*naux;
  pack_full_kernel<<<(unsigned)((tot+255)/256),256>>>(Mpk,T3,perm,ds,aperm,as,iu,ju,npair,naux,nao);
}
// FP32 output variant of the full-convention pack (for the TF32/FP32 whiten path).
__global__ void pack_full_f32_kernel(float* __restrict__ Mpk, const double* __restrict__ T3,
                            const int* __restrict__ perm, const double* __restrict__ ds,
                            const int* __restrict__ aperm, const double* __restrict__ as,
                            const int* __restrict__ iu, const int* __restrict__ ju,
                            int npair, int naux, int nao){
  long t = (long)blockIdx.x*blockDim.x + threadIdx.x;
  if(t >= (long)npair*naux) return;
  int pidx = (int)(t/naux); int P = (int)(t%naux);
  int i=iu[pidx], j=ju[pidx];
  int pi=perm[i], pj=perm[j], pp=aperm[P];
  double s=ds[i]*ds[j]*as[P];
  Mpk[(long)pidx*naux+P] = (float)(s * T3[((long)pi*nao+pj)*naux+pp]);
}
static void extern_pack_full_f32(float* Mpk,const double* T3,const int* perm,const double* ds,
                        const int* aperm,const double* as,const int* iu,const int* ju,
                        int npair,int naux,int nao){
  long tot=(long)npair*naux;
  pack_full_f32_kernel<<<(unsigned)((tot+255)/256),256>>>(Mpk,T3,perm,ds,aperm,as,iu,ju,npair,naux,nao);
}

#define CUCHECK(x) do{ cudaError_t e=(x); if(e!=cudaSuccess){ \
  fprintf(stderr,"CUDA error %s at %s:%d\n",cudaGetErrorString(e),__FILE__,__LINE__); exit(1);} }while(0)

static void to_chm(const PairF& P, double* FT, double scale){
  for(int ic=0;ic<P.ncomp;++ic)
    for(int ch=0;ch<P.nchtot;++ch)
      FT[(size_t)ch*P.ncomp+ic]=scale*P.F[(size_t)ic*P.nchtot+ch];
}

struct Shell { int l,np; double X[3]; vector<double> e,c; int ao0,ns; };
static vector<Shell> read_shells(const char* path,int& nao,bool cart=false){
  vector<Shell> sh;
  FILE* f=fopen(path,"r"); if(!f){fprintf(stderr,"cannot open %s\n",path);exit(1);}
  int ns; if(fscanf(f,"%d",&ns)!=1) exit(1);
  int ao=0;
  for(int s=0;s<ns;++s){ Shell S2;
    if(fscanf(f,"%d %lf %lf %lf %d",&S2.l,&S2.X[0],&S2.X[1],&S2.X[2],&S2.np)!=5) exit(1);
    S2.e.resize(S2.np); S2.c.resize(S2.np);
    for(int k=0;k<S2.np;++k) if(fscanf(f,"%lf %lf",&S2.e[k],&S2.c[k])!=2) exit(1);
    S2.ao0=ao; S2.ns=cart?((S2.l+1)*(S2.l+2))/2:2*S2.l+1; ao+=S2.ns; sh.push_back(S2); }
  fclose(f); nao=ao;
  return sh;
}


template<typename TR>
__global__ void scatter_t3(int n, int nsA, int nsB, int nsP, int nao, int naux,
                           const int* __restrict__ aoA, const int* __restrict__ aoB,
                           const int* __restrict__ aoP, const TR* __restrict__ blk,
                           double* __restrict__ T3){
  long t = (long)blockIdx.x*blockDim.x + threadIdx.x;
  long per = (long)nsA*nsB*nsP;
  if(t >= (long)n*per) return;
  int i = (int)(t/per); long r = t%per;
  int ma = (int)(r/(nsB*nsP)); int rb = (int)(r%(nsB*nsP));
  int mb = rb/nsP, mp = rb%nsP;
  double v = (double)blk[t];
  long a = aoA[i]+ma, b = aoB[i]+mb, pp = aoP[i]+mp;
  T3[(a*nao+b)*(long)naux+pp] = v;
  T3[(b*nao+a)*(long)naux+pp] = v;
}

// DIRECT-PACKED scatter: fold the convention transform (inverse perms + scales)
// into the scatter and write straight into the packed (npair, naux) whiten input
// Mpk[t(i,j)*naux + P] -- the dense (nao,nao,naux) T3 (2.3 GB at (H2O)16,
// 18.5 GB at (H2O)32) never exists, and the pack pass disappears. Same-shell
// diagonal blocks write the same (t,P) twice with the identical value (benign).
template<typename TR>
__global__ void scatter_t3_packed(int n, int nsA, int nsB, int nsP, int naux,
                           const int* __restrict__ aoA, const int* __restrict__ aoB,
                           const int* __restrict__ aoP, const TR* __restrict__ blk,
                           const int* __restrict__ iperm, const double* __restrict__ ds,
                           const int* __restrict__ iaperm, const double* __restrict__ as,
                           double* __restrict__ Mpk){
  long t = (long)blockIdx.x*blockDim.x + threadIdx.x;
  long per = (long)nsA*nsB*nsP;
  if(t >= (long)n*per) return;
  int i = (int)(t/per); long r = t%per;
  int ma = (int)(r/(nsB*nsP)); int rb = (int)(r%(nsB*nsP));
  int mb = rb/nsP, mp = rb%nsP;
  int a = aoA[i]+ma, b = aoB[i]+mb, pr = aoP[i]+mp;
  int io = iperm[a], jo = iperm[b], Po = iaperm[pr];
  int hi = io>jo?io:jo, lo = io>jo?jo:io;
  long tt = (long)hi*(hi+1)/2 + lo;
  Mpk[tt*(long)naux+Po] = ds[io]*ds[jo]*as[Po]*(double)blk[t];
}

// ===================== library entry (python ctypes) ========================
// The GPU engine as a LIBRARY: python passes the big convention-critical arrays
// (orbital/aux shells, packed H/S) in memory and gets the solve back; the small
// async artifacts (guess density, DFT grid) keep their existing file+poll
// channels so the caller can overlap their construction with the GPU B-build.
// Paths/knobs still arrive argv/env-style (python sets os.environ, then calls).
// With mem==nullptr this is exactly the old file-driven builder (main below).
struct OqpGpuShells {
  int nsh;                       // number of shells
  const int* l;                  // [nsh] angular momentum
  const double *x,*y,*z;         // [nsh] centers (Bohr)
  const int* np;                 // [nsh] primitives per shell
  const double *ex,*cc;          // concatenated primitives (coefs pre-normalized)
};
struct OqpGpuMemIn {
  OqpGpuShells orb, aux;
  const double* Hpk;             // packed lower-tri core Hamiltonian (ntri)
  const double* Spk;             // packed lower-tri overlap (ntri)
};
static vector<Shell> shells_from_mem(const OqpGpuShells& m,int& nao,bool cart){
  vector<Shell> sh; sh.reserve(m.nsh); int ao=0; long ip=0;
  for(int s=0;s<m.nsh;++s){ Shell S2;
    S2.l=m.l[s]; S2.X[0]=m.x[s]; S2.X[1]=m.y[s]; S2.X[2]=m.z[s]; S2.np=m.np[s];
    S2.e.assign(m.ex+ip,m.ex+ip+S2.np); S2.c.assign(m.cc+ip,m.cc+ip+S2.np); ip+=S2.np;
    S2.ao0=ao; S2.ns=cart?((S2.l+1)*(S2.l+2))/2:2*S2.l+1; ao+=S2.ns; sh.push_back(S2); }
  nao=ao; return sh;
}
extern "C" int oqpgpu_run(int argc,char** argv,const OqpGpuMemIn* mem);
int main(int argc,char** argv){ return oqpgpu_run(argc,argv,nullptr); }

extern "C" int oqpgpu_run(int argc,char** argv,const OqpGpuMemIn* mem){
  if(argc<5){fprintf(stderr,"usage: %s tables.bin shells.txt aux.txt out.bin\n",argv[0]);return 1;}
  double twA=now_s();
  RTC T;
  if(!load_tables(argv[1],T)){fprintf(stderr,"cannot load tables\n");return 1;}
  cudaFree(0); printf("  [t] tables+ctx %.2f\n",now_s()-twA);
  int mode=pin_rotc_mode(T);
  if(mode!=0){fprintf(stderr,"FATAL: pinned mode 0 required\n");return 1;}
  double TAU = argc>5 ? atof(argv[5]) : 1e-13;   // cascade knob (default: drop-only)
  // Distance-screening threshold (the rot Schwarz bound qS_pair*qS_aux*F0(eta*R^2)):
  // drop (pair,aux) blocks whose whitened magnitude is provably <= TAU_S, BEFORE
  // the global whiten (so the fit stays self-consistent on the surviving support).
  // Separate from TAU (which also drives the fp32/lcut banding of survivors).
  // Default = TAU (no extra screening); set ROUTEC_TAU_S>0 to screen for big
  // systems. Error ~ TAU_S, so pick it below the DF error (~1e-5 Ha).
  double TAU_S = getenv("ROUTEC_TAU_S") ? atof(getenv("ROUTEC_TAU_S")) : TAU;
  if(TAU_S < TAU) TAU_S = TAU;
  const double E32=3e-6, EMILD=0.05, EDEEP=0.15; // worst-case near-field error rules
  double tw0=now_s();

  const bool CART = getenv("ROUTEC_CART")!=nullptr;   // OpenQP 6d/10f bra
  int nao,naux;
  vector<Shell> sh = mem ? shells_from_mem(mem->orb,nao,CART) : read_shells(argv[2],nao,CART);
  vector<Shell> ax = mem ? shells_from_mem(mem->aux,naux,false) : read_shells(argv[3],naux,false);
  printf("int3c: %zu orbital shells nao=%d; %zu aux shells naux=%d\n",
         sh.size(),nao,ax.size(),naux);

  auto findM=[&](int la,int lb,int lc,int ld)->const RoutecMixedClass*{
    if(CART){
      for(int i=0;i<routec_naux_cart;++i){ const auto& M=routec_aux_cart[i];
        if(M.la==la&&M.lb==lb&&M.lc==lc&&M.ld==ld) return &M; }
      return nullptr; }
    for(int i=0;i<routec_nmixed;++i){ const auto& M=routec_mixed[i];
      if(M.la==la&&M.lb==lb&&M.lc==lc&&M.ld==ld) return &M; }
    for(int i=0;i<routec_naux;++i){ const auto& M=routec_aux[i];
      if(M.la==la&&M.lb==lb&&M.lc==lc&&M.ld==ld) return &M; }
    return nullptr; };

  // -------- orbital canonical pairs: pools per class, Schwarz
  struct CPair { int si,sj,cls,ppOff,npp; double P[3],pmin,qS; };
  std::map<int,int> clsmap; vector<std::array<int,2>> clsl;
  vector<vector<double>> poolFT, poolPC;
  auto clsid=[&](int lhi,int llo){ return lhi*8+llo; };
  auto getcls=[&](int lhi,int llo)->int{
    int cid=clsid(lhi,llo);
    if(!clsmap.count(cid)){ clsmap[cid]=(int)clsl.size();
      clsl.push_back({lhi,llo}); poolFT.push_back({}); poolPC.push_back({}); }
    return clsmap[cid]; };
  const int NSH=(int)sh.size();
  vector<CPair> pr;
  // The pair-table build (build_pair/to_chm) and the contracted Schwarz below
  // were SERIAL host loops scaling with npair -- 11+ s at (H2O)16 while the GPU
  // integral pass itself is 0.2 s. Two-pass restructure: pass 1 (cheap, serial)
  // assigns classes and pool offsets and pre-sizes the pools; pass 2 fills the
  // slots in PARALLEL (identical layout and content, bit-for-bit).
  vector<size_t> prFtOff, prFtSz;
  for(int i=0;i<NSH;++i)for(int j=0;j<=i;++j){
    int si=i,sj=j;
    if(sh[sj].l>sh[si].l) std::swap(si,sj);
    CPair P2; P2.si=si; P2.sj=sj;
    const Shell &A=sh[si],&B=sh[sj];
    P2.cls=getcls(A.l,B.l);
    const RoutecMixedClass* M=nullptr;
    if(CART){ for(int t=0;t<routec_naux_cart;++t){ const auto& Mc=routec_aux_cart[t];
        if(Mc.la==A.l&&Mc.lb==B.l){ M=&Mc; break; } } }
    else M=findM(A.l,B.l,A.l,B.l);
    size_t ftsz=(size_t)M->nchb*M->ncompb;
    auto& FT=poolFT[P2.cls]; auto& PC=poolPC[P2.cls];
    P2.ppOff=(int)(PC.size()/4); P2.npp=A.np*B.np;
    prFtOff.push_back(FT.size()); prFtSz.push_back(ftsz);
    FT.resize(FT.size()+(size_t)P2.npp*ftsz);
    PC.resize(PC.size()+(size_t)P2.npp*4);
    for(int x=0;x<3;++x) P2.P[x]=0.5*(A.X[x]+B.X[x]);
    pr.push_back(P2);
  }
  #pragma omp parallel for schedule(dynamic)
  for(int k2=0;k2<(int)pr.size();++k2){
    CPair& P2=pr[k2];
    const Shell &A=sh[P2.si],&B=sh[P2.sj];
    auto& FT=poolFT[P2.cls]; auto& PC=poolPC[P2.cls];
    const size_t ftsz=prFtSz[k2], fo=prFtOff[k2];
    double pmin=1e300; int pp=0;
    for(int pa=0;pa<A.np;++pa)for(int pb=0;pb<B.np;++pb,++pp){
      PairF PF; build_pair(T,A.l,B.l,A.X,B.X,A.e[pa],B.e[pb],PF,false,!CART);
      to_chm(PF,&FT[fo+(size_t)pp*ftsz],A.c[pa]*B.c[pb]);
      double p=A.e[pa]+B.e[pb]; pmin=std::min(pmin,p);
      double* pc=&PC[((size_t)P2.ppOff+pp)*4];
      for(int x=0;x<3;++x) pc[x]=(A.e[pa]*A.X[x]+B.e[pb]*B.X[x])/p;
      pc[3]=p;
    }
    P2.pmin=pmin;
  }
  // contracted Schwarz (parallel; np^4 routec_block calls per shell pair)
  #pragma omp parallel for schedule(dynamic)
  for(int k2=0;k2<(int)pr.size();++k2){
    CPair& P2=pr[k2];
    const Shell &A=sh[P2.si],&B=sh[P2.sj];
    int nb=(2*A.l+1)*(2*B.l+1);     // routec_block output is spherical
    vector<double> acc((size_t)nb*nb,0.0), blk;
    for(int pa=0;pa<A.np;++pa)for(int pb=0;pb<B.np;++pb)
      for(int pc=0;pc<A.np;++pc)for(int pd=0;pd<B.np;++pd){
        routec_block(T,mode,A.l,A.X,A.e[pa],B.l,B.X,B.e[pb],
                            A.l,A.X,A.e[pc],B.l,B.X,B.e[pd],blk,false);
        double w=A.c[pa]*B.c[pb]*A.c[pc]*B.c[pd];
        for(size_t z=0;z<acc.size();++z) acc[z]+=w*blk[z]; }
    double mx=0.0; for(int z=0;z<nb;++z) mx=std::max(mx,std::fabs(acc[(size_t)z*nb+z]));
    P2.qS=std::sqrt(mx);
  }

  // -------- aux shells as dummy-s pairs: pools per aux class, Schwarz
  struct APair { int s,cls,ppOff,npp; double qS, cmin; };
  std::map<int,int> aclsmap; vector<int> aclsl;
  vector<vector<double>> apoolFT, apoolPC;
  auto getacls=[&](int l)->int{
    if(!aclsmap.count(l)){ aclsmap[l]=(int)aclsl.size();
      aclsl.push_back(l); apoolFT.push_back({}); apoolPC.push_back({}); }
    return aclsmap[l]; };
  vector<APair> apr;
  for(int s=0;s<(int)ax.size();++s){
    const Shell& P=ax[s];
    APair A2; A2.s=s; A2.cls=getacls(P.l);
    int nchp=0; for(int k=0;k<=P.l;++k) nchp+=T.deg[k].nch;   // pair (l_P,0): K=l_P
    size_t ftsz=(size_t)nchp*(2*P.l+1);
    auto& FT=apoolFT[A2.cls]; auto& PC=apoolPC[A2.cls];
    A2.ppOff=(int)(PC.size()/4); A2.npp=P.np;
    for(int pa=0;pa<P.np;++pa){
      PairF PF; build_pair(T,P.l,0,P.X,P.X,P.e[pa],0.0,PF,false);
      size_t o=FT.size(); FT.resize(o+ftsz);
      to_chm(PF,&FT[o],P.c[pa]);
      for(int x=0;x<3;++x) PC.push_back(P.X[x]);
      PC.push_back(P.e[pa]);
    }
    // aux Schwarz: contracted (P|P)
    int nb=P.ns;
    vector<double> acc((size_t)nb*nb,0.0), blk;
    for(int pa=0;pa<P.np;++pa)for(int pc=0;pc<P.np;++pc){
      routec_block(T,mode,P.l,P.X,P.e[pa],0,P.X,0.0,
                          P.l,P.X,P.e[pc],0,P.X,0.0,blk,false);
      double w=P.c[pa]*P.c[pc];
      for(size_t z=0;z<acc.size();++z) acc[z]+=w*blk[z]; }
    double mx=0.0; for(int z=0;z<nb;++z) mx=std::max(mx,std::fabs(acc[(size_t)z*nb+z]));
    A2.qS=std::sqrt(mx);
    A2.cmin=1e300; for(int pa=0;pa<P.np;++pa) A2.cmin=std::min(A2.cmin,P.e[pa]);
    apr.push_back(A2);
  }
  double tw1=now_s();
  printf("  setup (pair tensors + Schwarz, CPU): %.2f s\n",tw1-tw0);

  // -------- device pools
  double _tpu0=now_s(); size_t _pby=0;
  for(auto&v:poolFT)_pby+=v.size()*8; for(auto&v:poolPC)_pby+=v.size()*8;
  for(auto&v:apoolFT)_pby+=v.size()*8; for(auto&v:apoolPC)_pby+=v.size()*8;
  vector<double*> dFT(clsl.size()), dPC(clsl.size());
  for(size_t c=0;c<clsl.size();++c){
    CUCHECK(cudaMalloc(&dFT[c],poolFT[c].size()*8));
    CUCHECK(cudaMemcpy(dFT[c],poolFT[c].data(),poolFT[c].size()*8,cudaMemcpyHostToDevice));
    CUCHECK(cudaMalloc(&dPC[c],poolPC[c].size()*8));
    CUCHECK(cudaMemcpy(dPC[c],poolPC[c].data(),poolPC[c].size()*8,cudaMemcpyHostToDevice));
  }
  vector<double*> dAFT(aclsl.size()), dAPC(aclsl.size());
  for(size_t c=0;c<aclsl.size();++c){
    CUCHECK(cudaMalloc(&dAFT[c],apoolFT[c].size()*8));
    CUCHECK(cudaMemcpy(dAFT[c],apoolFT[c].data(),apoolFT[c].size()*8,cudaMemcpyHostToDevice));
    CUCHECK(cudaMalloc(&dAPC[c],apoolPC[c].size()*8));
    CUCHECK(cudaMemcpy(dAPC[c],apoolPC[c].data(),apoolPC[c].size()*8,cudaMemcpyHostToDevice));
  }
  vector<float*> dFT32(clsl.size()), dAFT32(aclsl.size());
  { vector<float> tmp;
    for(size_t c=0;c<clsl.size();++c){
      tmp.assign(poolFT[c].begin(),poolFT[c].end());
      CUCHECK(cudaMalloc(&dFT32[c],tmp.size()*4));
      CUCHECK(cudaMemcpy(dFT32[c],tmp.data(),tmp.size()*4,cudaMemcpyHostToDevice)); }
    for(size_t c=0;c<aclsl.size();++c){
      tmp.assign(apoolFT[c].begin(),apoolFT[c].end());
      CUCHECK(cudaMalloc(&dAFT32[c],tmp.size()*4));
      CUCHECK(cudaMemcpy(dAFT32[c],tmp.data(),tmp.size()*4,cudaMemcpyHostToDevice)); }
  }
  printf("  [t] pools upload %.2f (%.0f MB fp64 + fp32 mirror)\n", now_s()-_tpu0, _pby/1e6);

  // -------- triplet groups
  struct Key{int bc,ac,kp,band;   // band: 0=fp64 1=fp32 2=mild-lcut+32 3=deep-lcut+32
    bool operator<(const Key&o)const{return std::tie(bc,ac,kp,band)<std::tie(o.bc,o.ac,o.kp,o.band);} };
  std::map<Key,vector<int2>> grp;
  // pair-row survival: a bra shell-pair with ZERO surviving aux is an exact zero
  // row of B (survives whitening) -> compactible storage. Track it to size the
  // real O(N^2) storage win (vs the triplet drop which is only a BUILD saving).
  std::vector<char> pairhit(pr.size(),0);
  long nb_[4]={0,0,0,0}; long ndrop=0;
  // Parallel triplet grouping (was a SERIAL double loop + std::map inserts = the
  // single biggest CPU cost, ~1.2s @ (H2O)16). Each thread buckets its (b,a)
  // slice into a thread-local map; merge afterwards. The 3c kernel writes each
  // triplet to its OWN T3 block, so triplet order within a group does not affect
  // the result -> the built B is bit-identical regardless of the parallel order.
  static const bool ALL32 = getenv("ROUTEC_ALL32")!=nullptr;
  static const bool FORCE64= getenv("ROUTEC_FORCE64")!=nullptr;
  {
    int NT=omp_get_max_threads();
    std::vector<std::map<Key,vector<int2>>> tg(NT);
    long tdrop=0, tnb0=0,tnb1=0,tnb2=0,tnb3=0;
    #pragma omp parallel reduction(+:tdrop,tnb0,tnb1,tnb2,tnb3)
    {
      auto& lg=tg[omp_get_thread_num()];
      #pragma omp for schedule(dynamic,32)
      for(int b=0;b<(int)pr.size();++b){
        char hit=0;
        for(int a=0;a<(int)apr.size();++a){
          const Shell& P=ax[apr[a].s];
          double R=0; for(int x=0;x<3;++x){double d=pr[b].P[x]-P.X[x]; R+=d*d;}
          R=std::sqrt(R);
          const double cmin=apr[a].cmin;
          double eta=pr[b].pmin*cmin/(pr[b].pmin+cmin);
          double F0[1]; boys(0,eta*R*R,F0);
          double bound=pr[b].qS*apr[a].qS*F0[0];
          if(bound<=TAU_S){ ++tdrop; continue; }
          hit=1;
          int band;
          if(FORCE64)               band=0;
          else if(ALL32)            band=1;
          else if(EDEEP*bound<=TAU) band=3;
          else if(EMILD*bound<=TAU) band=2;
          else if(E32*bound<=TAU)   band=1;
          else                      band=0;
          if(CART)                  band=0;
          if(band==0)tnb0++; else if(band==1)tnb1++; else if(band==2)tnb2++; else tnb3++;
          lg[{pr[b].cls,apr[a].cls,pr[b].npp*apr[a].npp,band}].push_back({b,a});
        }
        if(hit) pairhit[b]=1;
      }
    }
    ndrop=tdrop; nb_[0]=tnb0; nb_[1]=tnb1; nb_[2]=tnb2; nb_[3]=tnb3;
    for(auto& lg:tg) for(auto& kv:lg){ auto& dst=grp[kv.first];
      dst.insert(dst.end(),kv.second.begin(),kv.second.end()); }
  }
  long nkeep=nb_[0]+nb_[1]+nb_[2]+nb_[3];
  const long ntot=(long)pr.size()*(long)apr.size();
  const long nkeep_dbg=nb_[0]+nb_[1]+nb_[2]+nb_[3];
  long pzero=0; for(char h:pairhit) if(!h) ++pzero;       // fully-screened pair rows
  printf("  triplets (tau=%.0e tau_s=%.0e): fp64 %ld  fp32 %ld  mild+32 %ld  deep+32 %ld  drop %ld"
         "  (keep %ld/%ld = %.4f)  pair_rows: zero %ld/%ld sig=%.4f\n",
         TAU,TAU_S,nb_[0],nb_[1],nb_[2],nb_[3],ndrop,nkeep_dbg,ntot,
         ntot? (double)nkeep_dbg/(double)ntot : 1.0,
         pzero,(long)pr.size(), pr.size()? 1.0-(double)pzero/(double)pr.size() : 1.0);
  (void)nkeep;
  printf("  [t] post-grouping %.2f\n",now_s()-twA);


  // ============================================================
  // A2-FP32 device-resident B-build: 3c(FP32/cascade) -> on-device T3
  //   -> aux-conv GEMM -> orbital perm/scale + lower-tri PACK kernel
  //   -> Coulomb whiten GEMM (frozen Linv) -> B  (naux,npair), all on device.
  // No host T3, no file round-trip. Timed warm with min/median/max over reps.
  // Extra args: argv[6]=prep_dir argv[7]=tag(e.g. w16) argv[8]=nrep argv[9]=Bout(optional)
  // ============================================================
  const char* prep = argc>6 ? argv[6] : nullptr;
  const char* tag  = argc>7 ? argv[7] : "w";
  int NREP = argc>8 ? atoi(argv[8]) : 11;
  const char* Bout = argc>9 ? argv[9] : nullptr;
  if(!prep){ fprintf(stderr,"need prep dir as argv[6]\n"); return 1; }
  long npair=(long)nao*(nao+1)/2;

  // ---- load frozen transform operands (computed once by Python prep) ----
  auto loadbin=[&](const char* nm, void* dst, size_t bytes){
    char path[1024]; snprintf(path,sizeof(path),"%s/%s_%s.bin",prep,nm,tag);
    FILE* f=fopen(path,"rb"); if(!f){fprintf(stderr,"cannot open %s\n",path);exit(1);}
    size_t got=fread(dst,1,bytes,f); if(got!=bytes){fprintf(stderr,"short read %s %zu/%zu\n",path,got,bytes);exit(1);}
    fclose(f); };
  vector<double> hCag((size_t)naux*naux), hLinv((size_t)naux*naux), hDscale(nao), hAscale(naux);
  vector<int> hPerm(nao), hAperm(naux), hIu(npair), hJu(npair);
  if(!mem){
    loadbin("cag",hCag.data(),hCag.size()*8);
    loadbin("linv",hLinv.data(),hLinv.size()*8);
    loadbin("dscale",hDscale.data(),hDscale.size()*8);
    loadbin("perm",hPerm.data(),hPerm.size()*4);
    loadbin("ascale",hAscale.data(),hAscale.size()*8);
    loadbin("aperm",hAperm.data(),hAperm.size()*4);
    loadbin("iu",hIu.data(),hIu.size()*4);
    loadbin("ju",hJu.data(),hJu.size()*4);
  } else {
    // LIBRARY MODE: every transform operand is derived here -- no prep files.
    // Requires the cartesian frame + native metric (the pyscf-free production
    // path); the spherical file-driven variant keeps the loadbin route above.
    if(!CART){ fprintf(stderr,"oqpgpu_run: memory mode requires ROUTEC_CART=1\n"); return 1; }
    if(!getenv("ROUTEC_NATIVE_METRIC")){ fprintf(stderr,"oqpgpu_run: memory mode requires ROUTEC_NATIVE_METRIC=1\n"); return 1; }
    { long t2=0; for(long i2=0;i2<nao;++i2) for(long j2=0;j2<=i2;++j2,++t2){ hIu[t2]=(int)i2; hJu[t2]=(int)j2; } }
    // OpenQP cartesian component order + dscale (mirror of the python prep's
    // cart_transform: perm maps OpenQP component -> routec component index,
    // dscale = sqrt((2l-1)!!/((2a-1)!!(2b-1)!!(2c-1)!!)).
    static const int OQ[5][15][3]={
      {{0,0,0}},
      {{1,0,0},{0,1,0},{0,0,1}},
      {{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}},
      {{3,0,0},{0,3,0},{0,0,3},{2,1,0},{2,0,1},{1,2,0},{0,2,1},{1,0,2},{0,1,2},{1,1,1}},
      {{4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},{1,3,0},{0,3,1},{1,0,3},{0,1,3},
       {2,2,0},{2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}}};
    auto df2=[](int n){ double r=1.0; while(n>1){ r*=n; n-=2; } return r; };
    { int o=0;
      for(const Shell& S2: sh){ int l=S2.l, ncmp=(l+1)*(l+2)/2;
        if(l>4){ fprintf(stderr,"oqpgpu_run: l>4 unsupported\n"); return 1; }
        int rc[15][3]; int k=0;
        for(int axp=l;axp>=0;--axp) for(int ay=l-axp;ay>=0;--ay){ rc[k][0]=axp; rc[k][1]=ay; rc[k][2]=l-axp-ay; ++k; }
        for(int b2=0;b2<ncmp;++b2){ const int* t=OQ[l][b2];
          int pm=-1; for(int k2=0;k2<ncmp;++k2) if(rc[k2][0]==t[0]&&rc[k2][1]==t[1]&&rc[k2][2]==t[2]){pm=k2;break;}
          hPerm[o+b2]=o+pm;
          hDscale[o+b2]=sqrt(df2(2*l-1)/(df2(2*t[0]-1)*df2(2*t[1]-1)*df2(2*t[2]-1)));
        }
        o+=ncmp; } }
    // aux transform is identity in the native-metric path (Linv computed below)
    for(int P=0;P<naux;++P){ hAperm[P]=P; hAscale[P]=1.0; }
    std::fill(hCag.begin(),hCag.end(),0.0);
    for(int P=0;P<naux;++P) hCag[(size_t)P*naux+P]=1.0;
  }

  // ---- native (pyscf-free) auxiliary metric: recompute Linv from routec_block
  // and set the aux transform to identity (the fit index convention cancels in
  // B.B^T).  Keeps the orbital perm/dscale + iu/ju from the loaded prep so this
  // gate isolates the metric replacement.
  if(getenv("ROUTEC_NATIVE_METRIC")){
    double tm0=now_s();
    native_aux_linv(T,mode,ax,naux,hLinv);
    for(int P=0;P<naux;++P){ hAperm[P]=P; hAscale[P]=1.0; }
    std::fill(hCag.begin(),hCag.end(),0.0);
    for(int P=0;P<naux;++P) hCag[(size_t)P*naux+P]=1.0;
    printf("  native aux metric (routec_block + LAPACK chol/inv): %.2f s\n",now_s()-tm0);
  }
  double *dCag,*dLinv,*dDscale,*dAscale,*dM,*dMpk,*dB; int *dPerm,*dAperm,*dIu,*dJu;
  CUCHECK(cudaMalloc(&dCag,(size_t)naux*naux*8));
  CUCHECK(cudaMemcpy(dCag,hCag.data(),(size_t)naux*naux*8,cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dLinv,(size_t)naux*naux*8));
  CUCHECK(cudaMemcpy(dLinv,hLinv.data(),(size_t)naux*naux*8,cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dDscale,(size_t)nao*8));
  CUCHECK(cudaMemcpy(dDscale,hDscale.data(),(size_t)nao*8,cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dAscale,(size_t)naux*8));
  CUCHECK(cudaMemcpy(dAscale,hAscale.data(),(size_t)naux*8,cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dPerm,(size_t)nao*4));
  CUCHECK(cudaMemcpy(dPerm,hPerm.data(),(size_t)nao*4,cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dAperm,(size_t)naux*4));
  CUCHECK(cudaMemcpy(dAperm,hAperm.data(),(size_t)naux*4,cudaMemcpyHostToDevice));
  // inverse permutations for the direct-packed scatter (routec index -> OpenQP index)
  vector<int> hIperm(nao), hIaperm(naux);
  for(int i2=0;i2<nao;++i2)  hIperm[hPerm[i2]]=i2;
  for(int P2=0;P2<naux;++P2) hIaperm[hAperm[P2]]=P2;
  int *dIperm=nullptr,*dIaperm=nullptr;
  CUCHECK(cudaMalloc(&dIperm,(size_t)nao*4));
  CUCHECK(cudaMemcpy(dIperm,hIperm.data(),(size_t)nao*4,cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dIaperm,(size_t)naux*4));
  CUCHECK(cudaMemcpy(dIaperm,hIaperm.data(),(size_t)naux*4,cudaMemcpyHostToDevice));
  // PRODUCTION default: scatter straight into the packed whiten input; the dense
  // T3 is only needed by the bench (NREP>0) and the ASM/WHITEN variant paths.
  const bool PACKED = (NREP==0 && !getenv("ASM_FP32") && !getenv("ASM_FULL64")
                       && !getenv("WHITEN_TF32") && !getenv("WHITEN_FP32")
                       && !getenv("ROUTEC_DENSE_T3"));
  CUCHECK(cudaMalloc(&dIu,(size_t)npair*4));
  CUCHECK(cudaMemcpy(dIu,hIu.data(),(size_t)npair*4,cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dJu,(size_t)npair*4));
  CUCHECK(cudaMemcpy(dJu,hJu.data(),(size_t)npair*4,cudaMemcpyHostToDevice));
  dM=nullptr;   // aux-conv (nao,nao,naux) buffer -- 8*nao^2*naux bytes; ONLY the
                // ASM_FULL64 cross-check path uses it, so it is allocated lazily
                // there.  Skipping it here frees ~nao^2*naux*8 B (7.8 GB at 600
                // bf / 2712 aux), which is what lets large systems fit in 40 GB.
  CUCHECK(cudaMalloc(&dMpk,(size_t)npair*naux*8));    // packed (npair,naux)
  CUCHECK(cudaMalloc(&dB,(size_t)naux*npair*8));      // cderi (naux,npair)

  // ---- device T3 + kernel scratch ----
  double* dT3=nullptr;
  printf("  [t] pre-dT3 %.2f\n",now_s()-twA);
  if(!PACKED) CUCHECK(cudaMalloc(&dT3,(size_t)nao*nao*naux*8));
  const long PB=1u<<20; long CHUNKMAX=131072;
  double *dRs,*dTc,*dOut; float *dRs32,*dTc32,*dOut32; int2* dQall=nullptr;
  CUCHECK(cudaMalloc(&dRs,(size_t)PB*164*8));
  CUCHECK(cudaMalloc(&dTc,(size_t)PB*195*8));
  CUCHECK(cudaMalloc(&dOut,(size_t)CHUNKMAX*625*8));
  CUCHECK(cudaMalloc(&dRs32,(size_t)PB*164*4));
  CUCHECK(cudaMalloc(&dTc32,(size_t)PB*195*4));
  CUCHECK(cudaMalloc(&dOut32,(size_t)CHUNKMAX*625*4));
  auto findLV=[&](int la,int lb,int lc,int want)->const RoutecAuxLcut*{
    for(int i=0;i<routec_naux_lcut;++i){ const auto& L2=routec_aux_lcut[i];
      if(L2.la==la&&L2.lb==lb&&L2.lc==lc&&L2.lcut==want) return &L2; }
    return nullptr; };

  // ---- PRE-STAGE the full launch plan into device-resident buffers (host work
  //      done ONCE, off the timed path). Each "launch" = one chunk of one group. ----
  struct Launch { const RoutecMixedClass* M; const RoutecAuxLcut* LV; bool f32; int bc,ac,kp,n;
                  long qoff; int aoff; int nsA,nsB,nsP; };
  vector<Launch> plan;
  vector<int2> allQ; vector<int> allAOa,allAOb,allAOp;
  for(auto& kv: grp){
    const Key& K=kv.first; const auto& ql=kv.second;
    auto [lhi,llo]=clsl[K.bc]; int lp=aclsl[K.ac];
    const RoutecMixedClass* M=findM(lhi,llo,lp,0);
    if(!M){fprintf(stderr,"no kernel (%d%d|%d0)\n",lhi,llo,lp);return 1;}
    int Lt=lhi+llo+lp;
    const RoutecAuxLcut* LV=nullptr;
    if(K.band==3) LV=findLV(lhi,llo,lp,Lt-4);
    if(K.band==3&&!LV) LV=findLV(lhi,llo,lp,Lt-2);
    if(K.band==2) LV=findLV(lhi,llo,lp,Lt-2);
    bool f32 = (K.band>=1);
    long CH=std::min(CHUNKMAX,std::max(1L,PB/K.kp));
    for(size_t o=0;o<ql.size();o+=CH){
      int n=(int)std::min((size_t)CH,ql.size()-o);
      Launch L; L.M=M; L.LV=LV; L.f32=f32; L.bc=K.bc; L.ac=K.ac; L.kp=K.kp; L.n=n;
      L.qoff=(long)allQ.size(); L.aoff=(int)allAOa.size();
      int nsA=0,nsB=0,nsP=0;
      for(int i=0;i<n;++i){
        const CPair& B2=pr[ql[o+i].x]; const APair& A2=apr[ql[o+i].y];
        int z=0;
        for(int a2=0;a2<B2.npp;++a2)for(int b2=0;b2<A2.npp;++b2,++z)
          allQ.push_back({B2.ppOff+a2,A2.ppOff+b2});
        const Shell &SA=sh[B2.si],&SB=sh[B2.sj],&SP=ax[A2.s];
        allAOa.push_back(SA.ao0); allAOb.push_back(SB.ao0); allAOp.push_back(SP.ao0);
        nsA=SA.ns; nsB=SB.ns; nsP=SP.ns;
      }
      L.nsA=nsA; L.nsB=nsB; L.nsP=nsP;
      plan.push_back(L);
    }
  }
  int2* dQp; int *dAa,*dAb,*dAp;
  CUCHECK(cudaMalloc(&dQp,allQ.size()*sizeof(int2)));
  CUCHECK(cudaMemcpy(dQp,allQ.data(),allQ.size()*sizeof(int2),cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dAa,allAOa.size()*4)); CUCHECK(cudaMemcpy(dAa,allAOa.data(),allAOa.size()*4,cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dAb,allAOb.size()*4)); CUCHECK(cudaMemcpy(dAb,allAOb.data(),allAOb.size()*4,cudaMemcpyHostToDevice));
  CUCHECK(cudaMalloc(&dAp,allAOp.size()*4)); CUCHECK(cudaMemcpy(dAp,allAOp.data(),allAOp.size()*4,cudaMemcpyHostToDevice));
  printf("  plan: %zu launches, %zu quartet-prims staged on device\n",plan.size(),allQ.size());

  cublasHandle_t cub; cublasCreate(&cub);
  const double one=1.0, zero=0.0;

  // ---- the device-resident B-build (everything on device, no host transfer) ----
  auto run_3c = [&](){
    if(PACKED) CUCHECK(cudaMemsetAsync(dMpk,0,(size_t)npair*naux*8));
    else       CUCHECK(cudaMemsetAsync(dT3,0,(size_t)nao*nao*naux*8));
    for(auto& L: plan){
      long np=(long)L.n*L.kp;
      const int2* Q=dQp+L.qoff;
      if(!L.f32){
        L.M->scalar<<<(unsigned)((np+255)/256),256>>>(np,Q,dPC[L.bc],dAPC[L.ac],dRs,dTc);
        L.M->block<<<(unsigned)L.n,L.M->blockT,L.M->smemn*8>>>(L.n,L.kp,dRs,dTc,dFT[L.bc],dAFT[L.ac],Q,dOut);
        long tot=(long)L.n*L.nsA*L.nsB*L.nsP;
        if(PACKED) scatter_t3_packed<double><<<(unsigned)((tot+255)/256),256>>>(L.n,L.nsA,L.nsB,L.nsP,(int)naux,dAa+L.aoff,dAb+L.aoff,dAp+L.aoff,dOut,dIperm,dDscale,dIaperm,dAscale,dMpk);
        else scatter_t3<double><<<(unsigned)((tot+255)/256),256>>>(L.n,L.nsA,L.nsB,L.nsP,nao,naux,dAa+L.aoff,dAb+L.aoff,dAp+L.aoff,dOut,dT3);
      } else if(L.LV){
        L.LV->scalar32<<<(unsigned)((np+255)/256),256>>>(np,Q,dPC[L.bc],dAPC[L.ac],dRs32,dTc32);
        L.LV->block32<<<(unsigned)L.n,L.LV->blockT,L.LV->smemn*4>>>(L.n,L.kp,dRs32,dTc32,dFT32[L.bc],dAFT32[L.ac],Q,dOut32);
        long tot=(long)L.n*L.nsA*L.nsB*L.nsP;
        if(PACKED) scatter_t3_packed<float><<<(unsigned)((tot+255)/256),256>>>(L.n,L.nsA,L.nsB,L.nsP,(int)naux,dAa+L.aoff,dAb+L.aoff,dAp+L.aoff,dOut32,dIperm,dDscale,dIaperm,dAscale,dMpk);
        else scatter_t3<float><<<(unsigned)((tot+255)/256),256>>>(L.n,L.nsA,L.nsB,L.nsP,nao,naux,dAa+L.aoff,dAb+L.aoff,dAp+L.aoff,dOut32,dT3);
      } else {
        L.M->scalar32<<<(unsigned)((np+255)/256),256>>>(np,Q,dPC[L.bc],dAPC[L.ac],dRs32,dTc32);
        L.M->block32<<<(unsigned)L.n,L.M->blockT,L.M->smemn*4>>>(L.n,L.kp,dRs32,dTc32,dFT32[L.bc],dAFT32[L.ac],Q,dOut32);
        long tot=(long)L.n*L.nsA*L.nsB*L.nsP;
        if(PACKED) scatter_t3_packed<float><<<(unsigned)((tot+255)/256),256>>>(L.n,L.nsA,L.nsB,L.nsP,(int)naux,dAa+L.aoff,dAb+L.aoff,dAp+L.aoff,dOut32,dIperm,dDscale,dIaperm,dAscale,dMpk);
        else scatter_t3<float><<<(unsigned)((tot+255)/256),256>>>(L.n,L.nsA,L.nsB,L.nsP,nao,naux,dAa+L.aoff,dAb+L.aoff,dAp+L.aoff,dOut32,dT3);
      }
    }
  };
  // aux-conv: dM(nao*nao,naux) = dT3(nao*nao,naux) @ Cag^T.  Row-major X@Cag^T:
  // treat as col-major: compute (Cag @ T3rows^T)^T via one dgemm with swapped args.
  // We want M[r,P] = sum_c T3[r,c] Cag[P,c].  In col-major cublas with leading dims:
  //   C[P + r*naux] = sum_c Cag[P + c*?]... -> use: op(A)=Cag (naux x naux, but it is
  //   stored row-major => as col-major it's Cag^T). Simplest: M^T (naux x R) = Cag * T3^T.
  // We instead build M as (naux x R) col-major == M[P,r] row-major (R=nao*nao):
  //   Mt[P,r] = sum_c Cag[P,c] T3[r,c].  col-major: Mt = Cagcm^T? Handle below.
  long R=(long)nao*nao;
  // FP32 operands for the (benign, well-conditioned) aux-conv GEMM. The Cag rotation
  // is blockdiag inv(M[l]) with O(1) condition number, so FP32 there is lossless at
  // the fit-class level; the ill-conditioned Coulomb whiten stays FP64.
  // Assembly precision modes:
  //   default (PACK-FIRST FP64): pack orbital perm/scale + lower-tri FIRST so the
  //     aux-conv GEMM runs on npair rows (half of nao*nao); aux-conv & whiten FP64.
  //     -> preserves the cascade-bounded J/K (the near-field FP64 from the 3c band
  //        is NOT rounded). This is the validated, a-priori-bounded path.
  //   ASM_FP32 : aux-conv in FP32 (uniform ~1e-7 floor on ALL pairs, breaks the
  //     cascade bound; still < the 1e-6 SCF gate but reported separately).
  //   ASM_FULL64: original full nao*nao FP64 aux-conv (slowest), for cross-check.
  const bool ASM_FP32  = getenv("ASM_FP32")!=nullptr;
  const bool ASM_FULL64= getenv("ASM_FULL64")!=nullptr;
  // WHITEN_TF32: run the (otherwise FP64-bound) Coulomb whiten GEMM on TF32 tensor
  // cores (10-bit mantissa). The whiten matrix Linv=inv(chol(int2c2e)) is ill-
  // conditioned, so this is a precision risk -> validated empirically vs the gate.
  // WHITEN_FP32: plain FP32 (24-bit) SGEMM whiten -- exact-ish, no tensor cores.
  const bool WHITEN_TF32 = getenv("WHITEN_TF32")!=nullptr;
  const bool WHITEN_FP32 = getenv("WHITEN_FP32")!=nullptr;
  double *dT3pk=nullptr;
  if(!PACKED) CUCHECK(cudaMalloc(&dT3pk,(size_t)npair*naux*8));  // packed (npair,naux) FP64 (ASM/bench paths)
  // FP32 mirrors only when something can use them (the bench, or the FP32/TF32
  // assembly/whiten variants) -- production (nrep=0, default path) skips ~1.2 GB
  // of allocations + uploads it never touches.
  float *dCag32=nullptr,*dT3pk32=nullptr,*dMpk32=nullptr; const float onef=1.f, zerof=0.f;
  float *dLinv32=nullptr,*dB32=nullptr;
  if(NREP>0 || ASM_FP32 || WHITEN_TF32 || WHITEN_FP32){
    CUCHECK(cudaMalloc(&dCag32,(size_t)naux*naux*4));
    { vector<float> t(hCag.begin(),hCag.end()); CUCHECK(cudaMemcpy(dCag32,t.data(),t.size()*4,cudaMemcpyHostToDevice)); }
    CUCHECK(cudaMalloc(&dT3pk32,(size_t)npair*naux*4));
    CUCHECK(cudaMalloc(&dMpk32,(size_t)npair*naux*4));
    // FP32 Linv + FP32 B for the TF32/FP32 whiten paths.
    CUCHECK(cudaMalloc(&dLinv32,(size_t)naux*naux*4));
    { vector<float> t(hLinv.begin(),hLinv.end()); CUCHECK(cudaMemcpy(dLinv32,t.data(),t.size()*4,cudaMemcpyHostToDevice)); }
    CUCHECK(cudaMalloc(&dB32,(size_t)naux*npair*4));
  }
  auto run_assembly = [&](){
    if(ASM_FULL64){
      if(!dM) CUCHECK(cudaMalloc(&dM,(size_t)nao*nao*naux*8));
      cublasDgemm(cub,CUBLAS_OP_T,CUBLAS_OP_N,(int)naux,(int)R,(int)naux,
                  &one,dCag,(int)naux,dT3,(int)naux,&zero,dM,(int)naux);
      extern_pack(dMpk,dM,dPerm,dDscale,dIu,dJu,(int)npair,(int)naux,nao);
      cublasDgemm(cub,CUBLAS_OP_T,CUBLAS_OP_N,(int)npair,(int)naux,(int)naux,
                  &one,dMpk,(int)naux,dLinv,(int)naux,&zero,dB,(int)npair);
      return;
    }
    if(ASM_FP32){
      extern_pack_t3_f32(dT3pk32,dT3,dPerm,dDscale,dIu,dJu,(int)npair,(int)naux,nao);
      cublasSgemm(cub,CUBLAS_OP_T,CUBLAS_OP_N,(int)naux,(int)npair,(int)naux,
                  &onef,dCag32,(int)naux,dT3pk32,(int)naux,&zerof,dMpk32,(int)naux);
      extern_cast_f32_f64(dMpk,dMpk32,(long)npair*naux);
      cublasDgemm(cub,CUBLAS_OP_T,CUBLAS_OP_N,(int)npair,(int)naux,(int)naux,
                  &one,dMpk,(int)naux,dLinv,(int)naux,&zero,dB,(int)npair);
      return;
    }
    // DEFAULT: GEMM-FREE aux-conv. Cag is exactly blockdiag scalar*permutation
    // (verified err 0), so fold the aux perm/scale into the pack with the orbital
    // perm/scale + lower-tri, reading dense T3 in one kernel ->
    //   Mpk[t,P] = ds[i]ds[j]ascale[P] * T3[perm[i],perm[j],aperm[P]].
    // Only the (ill-conditioned, FP64) Coulomb whiten GEMM remains. J/K stays at
    // the cascade floor: no extra rounding from the convention transform.
    if(WHITEN_TF32 || WHITEN_FP32){
      // pack the full convention straight to FP32, then whiten on FP32/TF32 cores.
      extern_pack_full_f32(dMpk32,dT3,dPerm,dDscale,dAperm,dAscale,dIu,dJu,(int)npair,(int)naux,nao);
      cublasSetMathMode(cub, WHITEN_TF32 ? CUBLAS_TF32_TENSOR_OP_MATH
                                         : CUBLAS_PEDANTIC_MATH);
      cublasSgemm(cub,CUBLAS_OP_T,CUBLAS_OP_N,(int)npair,(int)naux,(int)naux,
                  &onef,dMpk32,(int)naux,dLinv32,(int)naux,&zerof,dB32,(int)npair);
      cublasSetMathMode(cub, CUBLAS_DEFAULT_MATH);
      extern_cast_f32_f64(dB,dB32,(long)naux*npair);
      return;
    }
    if(!PACKED)   // packed mode: the scatter already produced Mpk (transform folded)
      extern_pack_full(dMpk,dT3,dPerm,dDscale,dAperm,dAscale,dIu,dJu,(int)npair,(int)naux,nao);
    cublasDgemm(cub,CUBLAS_OP_T,CUBLAS_OP_N,(int)npair,(int)naux,(int)naux,
                &one,dMpk,(int)naux,dLinv,(int)naux,&zero,dB,(int)npair);
  };

  // NREP=0 -> PRODUCTION mode: skip the whole benchmark block (3x pipeline
  // warmup + clock-warmer GEMMs + timed reps ~= 4 s of pure redundancy when the
  // caller only wants the tensor) and build B exactly once below.
  if(NREP>0){
  cudaEvent_t e0,e1; cudaEventCreate(&e0); cudaEventCreate(&e1);
  // GEMM warmer to pull SM to boost clock (mirrors the python warm())
  { double* wa; double* wb; double* wc; int W=4096;
    CUCHECK(cudaMalloc(&wa,(size_t)W*W*8)); CUCHECK(cudaMalloc(&wb,(size_t)W*W*8)); CUCHECK(cudaMalloc(&wc,(size_t)W*W*8));
    CUCHECK(cudaMemset(wa,1,(size_t)W*W*8)); CUCHECK(cudaMemset(wb,1,(size_t)W*W*8));
    auto warm=[&](int it){ for(int k=0;k<it;k++) cublasDgemm(cub,CUBLAS_OP_N,CUBLAS_OP_N,W,W,W,&one,wa,W,wb,W,&zero,wc,W); CUCHECK(cudaDeviceSynchronize()); };

    // ---- warmup the actual pipeline (and JIT) ----
    for(int r=0;r<3;r++){ run_3c(); run_assembly(); }
    CUCHECK(cudaDeviceSynchronize());

    // ---- timed reps: full device-resident B-build, warm before each ----
    vector<double> tk(NREP), tasm(NREP), tfull(NREP);
    for(int r=0;r<NREP;r++){
      warm(40);
      cudaEventRecord(e0); run_3c(); cudaEventRecord(e1); CUCHECK(cudaEventSynchronize(e1));
      float ms; cudaEventElapsedTime(&ms,e0,e1); tk[r]=ms;
      cudaEventRecord(e0); run_assembly(); cudaEventRecord(e1); CUCHECK(cudaEventSynchronize(e1));
      cudaEventElapsedTime(&ms,e0,e1); tasm[r]=ms;
      tfull[r]=tk[r]+tasm[r];
    }
    // also time full as a single fused event (no inter-stage sync) for honesty
    vector<double> tfused(NREP);
    for(int r=0;r<NREP;r++){
      warm(40);
      cudaEventRecord(e0); run_3c(); run_assembly(); cudaEventRecord(e1); CUCHECK(cudaEventSynchronize(e1));
      float ms; cudaEventElapsedTime(&ms,e0,e1); tfused[r]=ms;
    }
    auto stat=[&](vector<double> v,const char* nm){
      std::sort(v.begin(),v.end());
      printf("  %s: min=%.1f median=%.1f max=%.1f ms (N=%d)\n",nm,v.front(),v[v.size()/2],v.back(),(int)v.size()); };
    printf("\n=== A2-FP32 device-resident B-build %s (single process, warm, N=%d) ===\n",tag,NREP);
    stat(tk,"3c kernel+scatter (FP32/cascade->device T3)");
    stat(tasm,"assembly+whiten (aux-conv GEMM + pack + Linv GEMM)");
    stat(tfull,"B-build (3c + assembly), stage-summed");
    stat(tfused,"B-build (3c + assembly), fused single event");
    double kmed=tk[NREP/2]; std::sort(tk.begin(),tk.end());
    std::sort(tfused.begin(),tfused.end());
    printf("@@FP32@@ tag=%s k3c_min=%.1f k3c_med=%.1f asm_med=%.1f full_med=%.1f full_min=%.1f\n",
           tag,tk.front(),tk[tk.size()/2],tasm[tasm.size()/2],tfused[tfused.size()/2],tfused.front());
  }
  }  // NREP>0 bench block

  // ---- emit B for validation (one extra build, off the timed path) ----
  if(Bout){
    run_3c(); run_assembly(); CUCHECK(cudaDeviceSynchronize());
    vector<double> hB((size_t)naux*npair);
    CUCHECK(cudaMemcpy(hB.data(),dB,hB.size()*8,cudaMemcpyDeviceToHost));
    // ---- CDF v2 (compacted): drop exactly-zero pair-row columns of B --------
    // hB is aux-major (hB[P*npair + t]). A column t that is zero over ALL aux is
    // a fully-screened pair row (its M row was dropped, and 0*Linv = 0), so it
    // is LOSSLESSLY droppable -- the O(N^2) storage lever. Emit keep_pairs (the
    // surviving original AO-pair indices) + the compacted (naux x ncp) payload.
    if(getenv("ROUTEC_CDF_V2")){
      std::vector<int> keep; keep.reserve(npair);
      for(int t=0;t<npair;++t){
        bool nz=false;
        for(int P=0;P<naux;++P) if(hB[(size_t)P*npair+t]!=0.0){ nz=true; break; }
        if(nz) keep.push_back(t);
      }
      const int ncp=(int)keep.size();
      std::vector<double> cB((size_t)naux*ncp);
      for(int P=0;P<naux;++P)
        for(int c=0;c<ncp;++c) cB[(size_t)P*ncp+c]=hB[(size_t)P*npair+keep[c]];
      // ---- CDF v3: magnitude-tiered mixed precision (the compression multiplier) --
      // Tier each aux column P (its ncp surviving pair-rows) by column absmax s_P
      // vs the global max S: s>=TH1*S -> fp64 (kept exact; the few large-mass
      // columns); TH2*S<=s<TH1*S -> fp16 (pow2 block scale); s<TH2*S -> int8
      // (symmetric). Arenas stay in natural aux order + per-column (tier,slot,scale)
      // so the consumer upcasts per column on reconstruction. Same tiering rule as
      // the df_store OQP_CDF_LOWPREC accuracy probe.
      if(getenv("ROUTEC_CDF_LOWPREC")){
        // 4-tier boundaries by column absmax s vs global max S (RATIO knobs, or
        // ABSOLUTE overrides A64/A32/A16). fp32 (rel err ~6e-8) is the safe
        // workhorse -> 2x; fp16/int8 compress only the small-magnitude tail.
        const double R64=getenv("OQP_CDF_TH64")?atof(getenv("OQP_CDF_TH64")):0.5;   // s>=R64*S -> fp64
        const double R32=getenv("OQP_CDF_TH32")?atof(getenv("OQP_CDF_TH32")):1e-3;  // s>=R32*S -> fp32
        const double R16=getenv("OQP_CDF_TH16")?atof(getenv("OQP_CDF_TH16")):1e-5;  // s>=R16*S -> fp16
        const double A64=getenv("OQP_CDF_A64")?atof(getenv("OQP_CDF_A64")):-1.0;
        const double A32=getenv("OQP_CDF_A32")?atof(getenv("OQP_CDF_A32")):-1.0;
        const double A16=getenv("OQP_CDF_A16")?atof(getenv("OQP_CDF_A16")):-1.0;
        const bool NOINT8=getenv("OQP_CDF_NOINT8")!=nullptr, NOFP16=getenv("OQP_CDF_NOFP16")!=nullptr;
        std::vector<double> smax(naux,0.0); double S=0.0;
        for(int P=0;P<naux;++P){ double m=0.0; const double* col=&cB[(size_t)P*ncp];
          for(int c=0;c<ncp;++c){ double a=fabs(col[c]); if(a>m)m=a; } smax[P]=m; if(m>S)S=m; }
        const double b64=(A64>0.0)?A64:R64*S, b32=(A32>0.0)?A32:R32*S, b16=(A16>0.0)?A16:R16*S;
        std::vector<unsigned char> tier(naux); std::vector<int> slot(naux);
        std::vector<float> scale(naux,1.0f); int n64=0,n32=0,n16=0,n8=0;
        for(int P=0;P<naux;++P){ double s=smax[P];
          if(s==0.0||s>=b64){ tier[P]=0; slot[P]=n64++; }              // fp64
          else if(s>=b32){    tier[P]=1; slot[P]=n32++; }              // fp32
          else if(NOFP16){    tier[P]=1; slot[P]=n32++; }              // (fp16 disabled -> fp32)
          else if(NOINT8||s>=b16){ tier[P]=2; slot[P]=n16++; scale[P]=(float)exp2(round(log2(s))); } // fp16
          else { tier[P]=3; slot[P]=n8++; scale[P]=(float)(s/127.0); } }// int8
        std::vector<double> B64((size_t)n64*ncp);
        std::vector<float>   B32((size_t)n32*ncp);
        std::vector<__half>  B16((size_t)n16*ncp);
        std::vector<signed char> B8((size_t)n8*ncp);
        for(int P=0;P<naux;++P){ const double* col=&cB[(size_t)P*ncp]; int sl=slot[P];
          if(tier[P]==0){ double* d=&B64[(size_t)sl*ncp]; for(int c=0;c<ncp;++c) d[c]=col[c]; }
          else if(tier[P]==1){ float* d=&B32[(size_t)sl*ncp]; for(int c=0;c<ncp;++c) d[c]=(float)col[c]; }
          else if(tier[P]==2){ __half* d=&B16[(size_t)sl*ncp]; double sc=scale[P];
            for(int c=0;c<ncp;++c) d[c]=__float2half((float)(col[c]/sc)); }
          else { signed char* d=&B8[(size_t)sl*ncp]; double inv=1.0/scale[P];
            for(int c=0;c<ncp;++c){ double q=nearbyint(col[c]*inv);
              if(q>127.0)q=127.0; else if(q<-127.0)q=-127.0; d[c]=(signed char)q; } } }
        FILE* f=fopen(Bout,"wb");
        const int magic=0x33464443 /*"CDF3"*/, n4=(int)naux, nao4=(int)nao, ncp4=ncp;
        fwrite(&magic,4,1,f); fwrite(&n4,4,1,f); fwrite(&nao4,4,1,f); fwrite(&ncp4,4,1,f);
        fwrite(keep.data(),4,ncp,f);
        fwrite(tier.data(),1,naux,f); fwrite(slot.data(),4,naux,f); fwrite(scale.data(),4,naux,f);
        fwrite(&n64,4,1,f); fwrite(&n32,4,1,f); fwrite(&n16,4,1,f); fwrite(&n8,4,1,f);
        fwrite(B64.data(),8,B64.size(),f);
        fwrite(B32.data(),4,B32.size(),f);
        fwrite(B16.data(),2,B16.size(),f);
        fwrite(B8.data(),1,B8.size(),f);
        fclose(f);
        const double mb=(4.0*ncp+naux*9.0+8.0*B64.size()+4.0*B32.size()+2.0*B16.size()+1.0*B8.size())/1e6;
        printf("  wrote CDF-v3 %s naux=%d nao=%d ncp=%d tiers[fp64/fp32/fp16/int8]=%d/%d/%d/%d "
               "(b64=%.1e b32=%.1e b16=%.1e) %.0f MB (v2 %.0f, dense %.0f)\n",
               Bout,naux,(int)nao,ncp,n64,n32,n16,n8,b64,b32,b16,mb,
               (4.0*ncp+8.0*(double)naux*ncp)/1e6, hB.size()*8/1e6);
        cudaDeviceSynchronize(); return 0;
      }
      FILE* f=fopen(Bout,"wb");
      const int magic=0x32464443 /*"CDF2"*/, n4=(int)naux, nao4=(int)nao, ncp4=ncp;
      fwrite(&magic,4,1,f); fwrite(&n4,4,1,f); fwrite(&nao4,4,1,f); fwrite(&ncp4,4,1,f);
      fwrite(keep.data(),4,ncp,f);
      fwrite(cB.data(),8,cB.size(),f); fclose(f);
      printf("  wrote CDF-v2 %s (naux=%d nao=%d ncp=%d/%d=%.4f, %.0f MB vs %.0f dense)\n",
             Bout,naux,(int)nao,ncp,npair,npair?(double)ncp/npair:1.0,
             (4.0*ncp+cB.size()*8)/1e6, hB.size()*8/1e6);
      cudaDeviceSynchronize(); return 0;
    }
    FILE* f=fopen(Bout,"wb"); int n4=(int)naux, p4=(int)npair;
    // ROUTEC_B_SIGMA: write the packed header the MRSF sigma-session / SCF
    // loader expects -- (naux, -nao), lower-tri i>=j slices -- instead of the
    // raw (naux, npair).  The payload (P-major, npair lower-tri) is identical.
    if(getenv("ROUTEC_B_SIGMA")){ int nn=-(int)nao; fwrite(&n4,4,1,f); fwrite(&nn,4,1,f); }
    else { fwrite(&n4,4,1,f); fwrite(&p4,4,1,f); }
    fwrite(hB.data(),8,hB.size(),f); fclose(f);
    printf("  [t] pre-write %.2f\n",now_s()-twA);
    printf("  wrote B %s (naux=%d %s=%d, %.0f MB)\n",Bout,n4,
           getenv("ROUTEC_B_SIGMA")?"nao":"npair",
           getenv("ROUTEC_B_SIGMA")?(int)nao:p4,hB.size()*8/1e6);
  }
  // ---- IN-PROCESS SOLVE (single process, B never touches disk) --------------
  // ROUTEC_SOLVE_HS=<prefix>: read <prefix>H.pk / <prefix>S.pk (packed lower-tri
  // fp64) and run routec_scf_solve right here on the just-built device tensor
  // (adopted by df_store; XC/guess through the usual envs).
  if(const char* hs = getenv("ROUTEC_SOLVE_HS")){
    // DFT: the grid generator runs concurrently with the build -- wait for its
    // (atomically renamed) grid.bin before the XC init inside the solve.
    if(const char* xd = getenv("OQP_OWNXC_DIR")){
      std::string gp = std::string(xd) + "/grid.bin"; struct stat gst;
      for(int w=0; w<1200 && stat(gp.c_str(),&gst)!=0; ++w) usleep(100000);
    }
    struct timespec _tb0,_tb1,_ts0,_ts1; clock_gettime(CLOCK_MONOTONIC,&_tb0);
    if(!Bout){ run_3c(); run_assembly(); CUCHECK(cudaDeviceSynchronize()); }
    // Free the whiten INPUT before allocating the full tensor -- it is dead
    // after run_assembly, and at (H2O)32 the difference decides OOM: peak drops
    // from Mpk 9.3 + Bpacked 9.3 + Bfull 18.5 = 37 GB to 27.8 GB.
    cudaFree(dMpk); dMpk=nullptr;
    double* dBfull=nullptr;
    CUCHECK(cudaMalloc(&dBfull,(size_t)naux*nao*nao*8));
    { long tot=(long)npair*naux;
      unpack_b_full_kernel<<<(unsigned)((tot+255)/256),256>>>(dBfull,dB,dIu,dJu,(long)npair,(int)naux,nao);
      CUCHECK(cudaGetLastError()); CUCHECK(cudaDeviceSynchronize()); }
    routec_b_adopt(dBfull,(int)naux,nao);
    cudaFree(dB); dB=nullptr;    // packed B dead once the full tensor is adopted
    clock_gettime(CLOCK_MONOTONIC,&_tb1);
    const long ntri=(long)nao*(nao+1)/2;
    vector<double> hpk, spk;
    const double *hp_=nullptr, *sp_=nullptr;
    if(mem){ hp_=mem->Hpk; sp_=mem->Spk; }
    else {
      hpk.resize(ntri); spk.resize(ntri);
      std::string hp=std::string(hs)+"H.pk", sp=std::string(hs)+"S.pk";
      FILE* fh=fopen(hp.c_str(),"rb"); FILE* fs=fopen(sp.c_str(),"rb");
      if(!fh||!fs){ fprintf(stderr,"solve: cannot open %s / %s\n",hp.c_str(),sp.c_str()); return 1; }
      if(fread(hpk.data(),8,ntri,fh)!=(size_t)ntri||fread(spk.data(),8,ntri,fs)!=(size_t)ntri){
        fprintf(stderr,"solve: short H/S read\n"); return 1; }
      fclose(fh); fclose(fs);
      hp_=hpk.data(); sp_=spk.data();
    }
    // The guess density may be produced by the (concurrent) prep process AFTER
    // the builder launched -- e.g. the native minao/SAD guess, whose analytic
    // overlap build overlaps this B-build. Wait for its atomically-renamed file
    // before routec_scf_solve reads it (no-op when the prep wrote it upfront).
    if(const char* gp = getenv("OQP_SCF_GUESS_D")){
      struct stat gs;
      for(int w=0; w<1200 && stat(gp,&gs)!=0; ++w) usleep(100000);
    }
    int nocc=getenv("ROUTEC_SOLVE_NOCC")?atoi(getenv("ROUTEC_SOLVE_NOCC")):0;
    double enuc=getenv("ROUTEC_SOLVE_ENUC")?atof(getenv("ROUTEC_SOLVE_ENUC")):0.0;
    double sc=1.0, se=getenv("ROUTEC_SOLVE_SE")?atof(getenv("ROUTEC_SOLVE_SE")):1.0;
    double ce=1e-8, cd=1e-7; int nb2=nao;
    int maxit=getenv("ROUTEC_SOLVE_MAXIT")?atoi(getenv("ROUTEC_SOLVE_MAXIT")):100;
    double E=0.0; int ncyc=0, info=1;
    vector<double> dout(ntri), cmo((size_t)nao*nao), eps(nao);
    clock_gettime(CLOCK_MONOTONIC,&_ts0);
    routec_scf_solve(hp_,sp_,&nb2,&nocc,&enuc,&sc,&se,&ce,&cd,&maxit,
                     &E,dout.data(),cmo.data(),eps.data(),&ncyc,&info);
    clock_gettime(CLOCK_MONOTONIC,&_ts1);
    double _bb=(_tb1.tv_sec-_tb0.tv_sec)+(_tb1.tv_nsec-_tb0.tv_nsec)*1e-9;
    double _sv=(_ts1.tv_sec-_ts0.tv_sec)+(_ts1.tv_nsec-_ts0.tv_nsec)*1e-9;
    printf("[inproc-solve %s] E=%.8f cyc=%d info=%d\n",tag,E,ncyc,info);
    printf("[inproc-timing %s] Bbuild=%.2fs solve=%.2fs (%.0f ms/cyc)\n",tag,_bb,_sv,ncyc>0?_sv*1000.0/ncyc:0.0);
  }
  cudaDeviceSynchronize();
  return 0;
}
