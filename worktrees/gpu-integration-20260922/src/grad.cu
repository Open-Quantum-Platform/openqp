// OpenQP G-4/G-5 DF 2e-GRADIENT supplier -- GPU (Gamma-contracted) variant.
// Same ABI as the CPU dylib (test/routec_oqp_grad.cpp):
//   routec_grad2(d, xyz, de, nbf, natm, hfscale, coulscale, info)
//   routec_grad2_mrsf(d, p, spc, xyz, de, nbf, natm, hfscale, hfscale2,
//                     coulscale, spcscale, mrst, info)
//
// Host side (input parsing, pair tensors, pass-1 energy integrals, the ONE
// Cholesky, Gamma/gamma assembly) is VERBATIM the CPU dylib.  Only pass 2
// -- the contraction of the derivative integrals d(ab|P)/dx, d(P|Q)/dx with
// Gamma/gamma -- moves to the GPU: the production Gamma-contracted path of
// GRAD3_NOTES.md S5.  One thread per primitive triple computes the emitted
// per-triple derivative pieces (V, P_loc, D_R; routec_grad_gen_cuda.cuh
// __device__ kernels) and immediately dots them against the gathered
// Gamma/gamma block of its shell triple; the 9 per-atom gradient scalars are
// accumulated through a shared-memory per-block atom buffer + global
// atomicAdd.  No derivative blocks are ever materialized in global memory.
//
// 2c metric term: same kernels via class (l_P,0|l_Q) with FdT = 0 (zero
// block at the head of the FdT pool); atA=atB=atom(P), atC=atom(Q) makes the
// (dA+dB)->P, dC->Q mapping automatic.
//
// Environment (in addition to the CPU dylib's):
//   OQP_ROUTEC_GRAD_GPU=0   force the CPU pass-2 (debug fallback)
//   OQP_ROUTEC_GRAD_BS      block size (default 128)
//
// Build (chc4):
//   module load CUDA/12.6.0 GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0
//   nvcc -O2 -arch=sm_80 -std=c++17 -Xcompiler -fPIC -shared \
//        routec_oqp_grad_gpu.cu -lopenblas -o libroutec_oqp_grad_gpu.so
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <array>
#include <map>
#include <string>
#include <atomic>
#include <thread>
#include <chrono>
#include <mutex>

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include "routec_tables_embedded.h"  // self-contained rotation tables (no file)

// ---- device Gamma/gamma helpers (see the gam_dev block in grad2e) ----------
// permute the (a,c,P) tensor: out[a][c][P] = in[c][a][P]
__global__ void gradgpu_perm102(const double* in, double* out,
                                int nao, int naux, long long tot){
  long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
  if(i>=tot) return;
  int P = (int)(i % naux);
  long long ac = i / naux;
  int c = (int)(ac % nao), a = (int)(ac / nao);
  out[((long long)a*nao + c)*naux + P] = in[((long long)c*nao + a)*naux + P];
}
// Gamma[a][b][P] = cs*D[a][b]*c[P] - 0.5*hs*Gamma[a][b][P]  (in place)
__global__ void gradgpu_gamma_elem(double* G, const double* D, const double* c,
                                   int nao, int naux, double hs, double cs,
                                   long long tot){
  long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
  if(i>=tot) return;
  int P = (int)(i % naux);
  long long ab = i / naux;
  int b = (int)(ab % nao), a = (int)(ab / nao);
  G[i] = cs*D[(long long)a*nao+b]*c[P] - 0.5*hs*G[i];
}
// MRSF J-term rank-1: Gam[t][P] += w*(A[t]*cB[P] + B[t]*cA[P]), t=(a*nao+b)
__global__ void gradgpu_mrsf_rank1(double* Gam, const double* A, const double* B,
                                   const double* cA, const double* cB,
                                   double w, int naux, long long tot){
  long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
  if(i>=tot) return;
  int P = (int)(i % naux);
  long long t = i / naux;
  Gam[i] += w*(A[t]*cB[P] + B[t]*cA[P]);
}
// MRSF K-term gam: gam[p][q] -= 0.5*w*(T[p][q] + T[q][p])
__global__ void gradgpu_mrsf_gamsym(double* gam, const double* T,
                                    double w, int naux, long long tot){
  long long i = (long long)blockIdx.x*blockDim.x + threadIdx.x;
  if(i>=tot) return;
  int q = (int)(i % naux);
  int p = (int)(i / naux);
  gam[i] -= 0.5*w*(T[i] + T[(long long)q*naux+p]);
}

#include <cblas.h>
extern "C" {
void dpotrf_(const char*, const int*, double*, const int*, int*);
void dpotrs_(const char*, const int*, const int*, const double*, const int*,
             double*, const int*, int*);
}
typedef int lapack_int_t;

#include "routec_grad_gen.inc"        // emitted CPU kernels + extern-C drivers
#include "routec_grad_gen_cuda.cuh"   // emitted __device__ kernels

using std::vector;

namespace oqpgrad {

static double now_s(){using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();}

// ------------------------------------------------------------- monomials
static vector<std::array<int,3>> monos_of(int k){
  vector<std::array<int,3>> v;
  for(int lx=k; lx>=0; --lx) for(int ly=k-lx; ly>=0; --ly)
    v.push_back({lx,ly,k-lx-ly});
  return v;
}
static int ncart(int l){ return (l+1)*(l+2)/2; }

// ------------------------------------- RTC tables (same loader as G2 C++)
struct DegTab {
  int nmono=0, nch=0, nch0=0;
  vector<std::array<int,3>> ch;
  vector<double> C;               // nch x nmono row-major
  vector<int> ch0;
  vector<double> B0;
};
struct RTC {
  int kmax=0, kpair=0;
  vector<DegTab> deg;
};
// core parser: reads an RTC2 stream from an already-open FILE* (file or, via
// fmemopen, the embedded byte array).  Does NOT close f.
static bool load_tables_FILE(FILE* f, RTC& T){
  char magic[4]; int km,kp;
  if(fread(magic,1,4,f)!=4||memcmp(magic,"RTC2",4)) return false;
  if(fread(&km,4,1,f)!=1||fread(&kp,4,1,f)!=1) return false;
  T.kmax=km; T.kpair=kp;
  T.deg.resize(km+1);
  for(int k=0;k<=km;++k){ DegTab& D=T.deg[k];
    if(fread(&D.nmono,4,1,f)!=1||fread(&D.nch,4,1,f)!=1||
       fread(&D.nch0,4,1,f)!=1) return false;
    D.ch.resize(D.nch);
    for(int j=0;j<D.nch;++j){int b[3];
      if(fread(b,4,3,f)!=3) return false;
      D.ch[j]={b[0],b[1],b[2]};}
    D.C.resize((size_t)D.nch*D.nmono);
    if(fread(D.C.data(),8,D.C.size(),f)!=D.C.size()) return false;
    D.ch0.resize(D.nch0);
    if(D.nch0&&fread(D.ch0.data(),4,D.nch0,f)!=(size_t)D.nch0) return false;
    D.B0.resize((size_t)D.nmono*D.nch0);
    if(D.B0.size()&&fread(D.B0.data(),8,D.B0.size(),f)!=D.B0.size()) return false;
  }
  return true;
}
static bool load_tables(const char* path, RTC& T){
  FILE* f=fopen(path,"rb"); if(!f) return false;
  bool ok=load_tables_FILE(f,T); fclose(f); return ok;
}
// same parser, over the embedded byte array (no external file needed).
static bool load_tables_mem(const unsigned char* buf, size_t len, RTC& T){
  FILE* f=fmemopen((void*)buf,len,"rb"); if(!f) return false;
  bool ok=load_tables_FILE(f,T); fclose(f); return ok;
}

// ------------------------------------------- Hermite E full table (one axis)
struct ETab {
  int I=0, J=0, W=0;
  vector<double> v;
  const double* at(int i,int j) const { return &v[((size_t)i*J+j)*W]; }
};
static void e_table(int imax,int jmax,double Qx,double a,double b,ETab& E){
  const double p=a+b, mu=(b==0.0?0.0:a*b/p), inv2p=0.5/p;
  E.I=imax+1; E.J=jmax+1; E.W=imax+jmax+1;
  E.v.assign((size_t)E.I*E.J*E.W,0.0);
  auto at=[&](int ii,int jj)->double*{ return &E.v[((size_t)ii*E.J+jj)*E.W]; };
  at(0,0)[0]=std::exp(-mu*Qx*Qx);
  for(int ii=1;ii<E.I;++ii){ double* cur=at(ii,0); double* lo=at(ii-1,0);
    for(int t=0;t<=ii;++t){ double v=0.0;
      if(t-1>=0) v += inv2p*lo[t-1];
      v += (-(mu*Qx)/a)*lo[t];
      if(t+1<=ii-1) v += (t+1)*lo[t+1];
      cur[t]=v; } }
  for(int jj=1;jj<E.J;++jj) for(int ii=0;ii<E.I;++ii){
    double* cur=at(ii,jj); double* lo=at(ii,jj-1);
    for(int t=0;t<=ii+jj;++t){ double v=0.0;
      if(t-1>=0) v += inv2p*lo[t-1];
      v += ((mu*Qx)/b)*lo[t];
      if(t+1<=ii+jj-1) v += (t+1)*lo[t+1];
      cur[t]=v; } }
}

// -------------------------------------------------- pair channel tensors
static void pair_T(const RTC& T,int la,int lb,const double* A,const double* B,
                   double a,double b,int K,int nchK,double* FT /*zeroed*/){
  ETab E[3];
  for(int ax=0;ax<3;++ax) e_table(la,lb,A[ax]-B[ax],a,b,E[ax]);
  auto ca=monos_of(la), cb=monos_of(lb);
  const int ncomp=(int)(ca.size()*cb.size());
  vector<int> off(K+1); int tot=0;
  for(int k=0;k<=K;++k){ off[k]=tot; tot+=T.deg[k].nch; }
  vector<double> e;
  int ic=0;
  for(auto& ia: ca) for(auto& ib: cb){
    const double* Ex=E[0].at(ia[0],ib[0]); const int lx=ia[0]+ib[0]+1;
    const double* Ey=E[1].at(ia[1],ib[1]); const int ly=ia[1]+ib[1]+1;
    const double* Ez=E[2].at(ia[2],ib[2]); const int lz=ia[2]+ib[2]+1;
    const int kmaxc=(lx-1)+(ly-1)+(lz-1);
    for(int k=0;k<=std::min(K,kmaxc);++k){
      const DegTab& D=T.deg[k];
      auto ml=monos_of(k);
      e.assign(D.nmono,0.0); bool nz=false;
      for(int m=0;m<D.nmono;++m){ auto& mo=ml[m];
        if(mo[0]<lx&&mo[1]<ly&&mo[2]<lz){
          double v=Ex[mo[0]]*Ey[mo[1]]*Ez[mo[2]];
          if(v!=0.0){e[m]=v;nz=true;} } }
      if(!nz) continue;
      for(int chi=0;chi<D.nch;++chi){ double s=0.0;
        const double* Cr=&D.C[(size_t)chi*D.nmono];
        for(int m=0;m<D.nmono;++m) s+=Cr[m]*e[m];
        FT[(size_t)(off[k]+chi)*ncomp+ic]=s; }
    }
    ++ic;
  }
  (void)nchK;
}

static void pair_dT(const RTC& T,int la,int lb,const double* A,const double* B,
                    double a,double b,int nch1,double* FdT /*zeroed*/){
  const double p=a+b;
  const int K1=la+lb;
  ETab E[3];
  for(int ax=0;ax<3;++ax) e_table(la+1,lb,A[ax]-B[ax],a,b,E[ax]);
  auto ca=monos_of(la), cb=monos_of(lb);
  const int ncomp=(int)(ca.size()*cb.size());
  vector<int> off(K1+2); int tot=0;
  for(int k=0;k<=K1+1;++k){ off[k]=tot; tot+=T.deg[k].nch; }
  auto dE=[&](int ax,int i,int j,double* out){
    const double* Ep=E[ax].at(i+1,j); const int lp_=i+1+j+1;
    const double* Em=(i>=1)?E[ax].at(i-1,j):nullptr; const int lm=(i>=1)?(i-1+j+1):0;
    const double* E0=E[ax].at(i,j); const int l0=i+j+1;
    for(int t=0;t<=i+j+1;++t){
      double v=2.0*a*((t<lp_)?Ep[t]:0.0);
      if(Em&&t<lm) v-=i*Em[t];
      if(t>=1&&t<=l0) v-=(a/p)*E0[t-1];
      out[t]=v; }
  };
  vector<double> dbuf(K1+3), e;
  for(int w=0;w<3;++w){
    double* Fw=FdT+(size_t)w*nch1*ncomp;
    int ic=0;
    for(auto& ia: ca) for(auto& ib: cb){
      const double* Et[3]; int len[3];
      for(int ax=0;ax<3;++ax){
        if(ax==w){ dE(ax,ia[ax],ib[ax],dbuf.data());
          Et[ax]=dbuf.data(); len[ax]=ia[ax]+ib[ax]+2; }
        else { Et[ax]=E[ax].at(ia[ax],ib[ax]); len[ax]=ia[ax]+ib[ax]+1; }
      }
      const int kmaxc=(len[0]-1)+(len[1]-1)+(len[2]-1);
      for(int k=0;k<=std::min(K1+1,kmaxc);++k){
        const DegTab& D=T.deg[k];
        auto ml=monos_of(k);
        e.assign(D.nmono,0.0); bool nz=false;
        for(int m=0;m<D.nmono;++m){ auto& mo=ml[m];
          if(mo[0]<len[0]&&mo[1]<len[1]&&mo[2]<len[2]){
            double v=Et[0][mo[0]]*Et[1][mo[1]]*Et[2][mo[2]];
            if(v!=0.0){e[m]=v;nz=true;} } }
        if(!nz) continue;
        for(int chi=0;chi<D.nch;++chi){ double s=0.0;
          const double* Cr=&D.C[(size_t)chi*D.nmono];
          for(int m=0;m<D.nmono;++m) s+=Cr[m]*e[m];
          Fw[(size_t)(off[k]+chi)*ncomp+ic]=s; }
      }
      ++ic;
    }
  }
}

// ----------------------------------------------------------- input parsing
struct Sub { int l=0, atom=0, np=0, ao0=0, ns=0; vector<double> e, c; };

struct Ctx {
  bool loaded=false, ok=false;
  RTC T;
  int natm=0, nbf=0, naux=0;
  vector<Sub> sub, asub;
  vector<int> qmap;
  vector<double> smap;
  struct Dims { int nchb,nchb1,nchk,ncb,nck; };
  std::map<int,Dims> dims;
  struct PairT { vector<double> FT, FdT; int npp=0; };
  vector<PairT> pcache;
  vector<std::array<int,2>> pairs;
  std::map<std::pair<int,long long>,vector<double>> ketcache;
  int nthreads=1;
  bool verbose=false;
  double last_xyz_hash=0.0/0.0;
};
static Ctx g;
static std::mutex g_mtx;
// In-memory deck handed in via routec_grad_set_deck() (the tagarray path): when
// non-empty it replaces the OQP_ROUTEC_GRAD_INP file entirely. Set once, before
// the first gradient call.
static std::string g_deck_buf;

static bool read_subs(FILE* f,const char* tag,int& count,vector<Sub>& out,int& nao){
  char buf[128];
  if(fscanf(f,"%127s %d",buf,&count)!=2||strcmp(buf,tag)!=0) return false;
  out.resize(count);
  nao=0;
  for(int s=0;s<count;++s){ Sub& S=out[s];
    if(fscanf(f,"%d %d %d",&S.l,&S.atom,&S.np)!=3) return false;
    S.e.resize(S.np); S.c.resize(S.np);
    for(int k2=0;k2<S.np;++k2)
      if(fscanf(f,"%lf %lf",&S.e[k2],&S.c[k2])!=2) return false;
    S.ao0=nao; S.ns=ncart(S.l); nao+=S.ns;
  }
  return true;
}

static const Ctx::Dims& class_dims(int la,int lb,int lp){
  int key=(la*8+lb)*8+lp;
  auto it=g.dims.find(key);
  if(it!=g.dims.end()) return it->second;
  int d[5];
  if(routec_grad3c_dims(la,lb,lp,d)!=0){
    fprintf(stderr,"[routec-grad-gpu] FATAL: no emitted class (%d%d|%d)\n",la,lb,lp);
    abort();
  }
  Ctx::Dims dd{d[0],d[1],d[2],d[3],d[4]};
  return g.dims.emplace(key,dd).first->second;
}

static const vector<double>& ket_T(int lp,double q);   // fwd

static bool load_input(){
  if(g.loaded) return g.ok;
  g.loaded=true; g.ok=false;
  // Rotation tables: prefer an explicit OQP_ROUTEC_TABLES file (lets a user
  // supply a newer table), else the copy embedded in the library -- so the
  // shipped .so is self-contained and needs no external routec_tables.bin.
  const char* tabp=getenv("OQP_ROUTEC_TABLES");
  bool tok = tabp ? load_tables(tabp,g.T) : false;
  if(!tok) tok = load_tables_mem(routec_tables_embedded,
                                 routec_tables_embedded_len, g.T);
  if(!tok){
    fprintf(stderr,"[routec-grad-gpu] cannot load tables (env=%s, embedded=%zu B)\n",
            tabp?tabp:"(unset)", routec_tables_embedded_len); return false; }
  // Deck (shells/aux/AO map): prefer an in-memory deck handed in via
  // routec_grad_set_deck() -- the tagarray path, no file at all -- otherwise
  // fall back to the OQP_ROUTEC_GRAD_INP file. Parsed once, at first call.
  FILE* f=nullptr; const char* src=nullptr;
  if(!g_deck_buf.empty()){
    f=fmemopen((void*)g_deck_buf.data(),g_deck_buf.size(),"r"); src="<deck-tag>";
  }else{
    const char* inp=getenv("OQP_ROUTEC_GRAD_INP");
    if(!inp){ fprintf(stderr,"[routec-grad-gpu] no deck: call routec_grad_set_deck() "
                      "or set OQP_ROUTEC_GRAD_INP\n"); return false; }
    f=fopen(inp,"r"); src=inp;
  }
  if(!f){ fprintf(stderr,"[routec-grad-gpu] cannot open deck %s\n",src); return false; }
  char buf[128], ver[32];
  int nbf_hdr=0, naux_hdr=0;
  if(fscanf(f,"%127s %31s",buf,ver)!=2||strcmp(buf,"routec_grad_inp")!=0){
    fprintf(stderr,"[routec-grad-gpu] bad header in deck %s\n",src); fclose(f); return false; }
  if(fscanf(f,"%127s %d",buf,&g.natm)!=2||strcmp(buf,"natm")) {fclose(f);return false;}
  if(fscanf(f,"%127s %d",buf,&nbf_hdr)!=2||strcmp(buf,"nbf"))  {fclose(f);return false;}
  if(fscanf(f,"%127s %d",buf,&naux_hdr)!=2||strcmp(buf,"naux")){fclose(f);return false;}
  int ns=0,nas=0;
  if(!read_subs(f,"nsub",ns,g.sub,g.nbf)){fclose(f);return false;}
  if(!read_subs(f,"nauxsub",nas,g.asub,g.naux)){fclose(f);return false;}
  if(g.nbf!=nbf_hdr||g.naux!=naux_hdr){
    fprintf(stderr,"[routec-grad-gpu] header/subshell count mismatch (%d/%d vs %d/%d)\n",
            nbf_hdr,naux_hdr,g.nbf,g.naux); fclose(f); return false; }
  if(fscanf(f,"%127s",buf)!=1||strcmp(buf,"map")){fclose(f);return false;}
  g.qmap.resize(g.nbf); g.smap.resize(g.nbf);
  for(int i=0;i<g.nbf;++i)
    if(fscanf(f,"%d %lf",&g.qmap[i],&g.smap[i])!=2){fclose(f);return false;}
  fclose(f);
  for(auto& S:g.sub) if(S.l>2){
    fprintf(stderr,"[routec-grad-gpu] orbital l=%d unsupported (classes are s,p,d)\n",S.l);
    return false; }
  for(auto& S:g.asub) if(S.l>4){
    fprintf(stderr,"[routec-grad-gpu] aux l=%d unsupported (classes go to g)\n",S.l);
    return false; }
  for(int i=0;i<(int)g.sub.size();++i)
    for(int j=0;j<=i;++j){
      if(g.sub[i].l>=g.sub[j].l) g.pairs.push_back({i,j});
      else                       g.pairs.push_back({j,i});
    }
  g.pcache.resize(g.pairs.size());
  const char* th=getenv("OQP_ROUTEC_GRAD_THREADS");
  g.nthreads = th ? std::max(1,atoi(th))
                  : std::max(1u,std::thread::hardware_concurrency());
  g.verbose = getenv("OQP_ROUTEC_GRAD_VERBOSE")!=nullptr;
  for(auto& pr:g.pairs)
    for(auto& sp:g.asub)
      (void)class_dims(g.sub[pr[0]].l,g.sub[pr[1]].l,sp.l);
  for(auto& sp:g.asub)
    for(auto& sq:g.asub)
      (void)class_dims(sp.l,0,sq.l);
  for(auto& sp:g.asub)
    for(int k2=0;k2<sp.np;++k2)
      (void)ket_T(sp.l,sp.e[k2]);
  fprintf(stderr,"[routec-grad-gpu] input loaded: natm=%d nbf=%d naux=%d "
          "(%zu orbital / %zu aux subshells), %d threads\n",
          g.natm,g.nbf,g.naux,g.sub.size(),g.asub.size(),g.nthreads);
  g.ok=true;
  return true;
}

static const vector<double>& ket_T(int lp,double q){
  long long bits; memcpy(&bits,&q,8);
  auto key=std::make_pair(lp,bits);
  auto it=g.ketcache.find(key);
  if(it!=g.ketcache.end()) return it->second;
  int nch=0; for(int k=0;k<=lp;++k) nch+=g.T.deg[k].nch;
  vector<double> Fk((size_t)nch*ncart(lp),0.0);
  double Z[3]={0,0,0};
  pair_T(g.T,lp,0,Z,Z,q,0.0,lp,nch,Fk.data());
  return g.ketcache.emplace(key,std::move(Fk)).first->second;
}

static const Ctx::PairT& pair_get(int ip,const double* coords){
  Ctx::PairT& P=g.pcache[ip];
  if(P.npp) return P;
  const Sub& sa=g.sub[g.pairs[ip][0]];
  const Sub& sb=g.sub[g.pairs[ip][1]];
  const double* A=coords+3*sa.atom;
  const double* B=coords+3*sb.atom;
  const Ctx::Dims& d=class_dims(sa.l,sb.l,0);
  const int ncb=d.ncb;
  const size_t sbsz=(size_t)d.nchb*ncb, sdsz=3*(size_t)d.nchb1*ncb;
  P.npp=sa.np*sb.np;
  P.FT.assign((size_t)P.npp*sbsz,0.0);
  P.FdT.assign((size_t)P.npp*sdsz,0.0);
  int r=0;
  for(int pa=0;pa<sa.np;++pa)
    for(int pb=0;pb<sb.np;++pb,++r){
      pair_T (g.T,sa.l,sb.l,A,B,sa.e[pa],sb.e[pb],sa.l+sb.l,d.nchb,
              P.FT.data()+(size_t)r*sbsz);
      pair_dT(g.T,sa.l,sb.l,A,B,sa.e[pa],sb.e[pb],d.nchb1,
              P.FdT.data()+(size_t)r*sdsz);
    }
  return P;
}

// one contracted 3c block (CPU; pass 1 and the CPU pass-2 fallback)
struct BlockWS {
  vector<double> FT, FdT, FkT, geo, V, dA, dB, dC;
};
static void block3c(int ipair,const Sub& sa,const Sub& sb,const Sub& sp,
                    const double* coords,BlockWS& w,bool /*unused*/){
  const Ctx::Dims& d=class_dims(sa.l,sb.l,sp.l);
  const int ncb=d.ncb, nck=d.nck;
  const Ctx::PairT& P=pair_get(ipair,coords);
  const size_t sbsz=(size_t)d.nchb*ncb, sdsz=3*(size_t)d.nchb1*ncb,
               sksz=(size_t)d.nchk*nck;
  const int npp=P.npp, npq=sp.np, np=npp*npq;
  w.FT.resize((size_t)np*sbsz); w.FdT.resize((size_t)np*sdsz);
  w.FkT.resize((size_t)np*sksz); w.geo.resize((size_t)np*8);
  const double* A=coords+3*sa.atom;
  const double* B=coords+3*sb.atom;
  const double* C=coords+3*sp.atom;
  int ip=0;
  for(int r=0;r<npp;++r){
    const int pa=r/sb.np, pb=r%sb.np;
    const double a=sa.e[pa], b=sb.e[pb], p=a+b;
    const double cc=sa.c[pa]*sb.c[pb];
    double Pc[3]={(a*A[0]+b*B[0])/p,(a*A[1]+b*B[1])/p,(a*A[2]+b*B[2])/p};
    for(int kq=0;kq<npq;++kq,++ip){
      const double q=sp.e[kq];
      memcpy(&w.FT [ip*sbsz],&P.FT [(size_t)r*sbsz],sbsz*8);
      memcpy(&w.FdT[ip*sdsz],&P.FdT[(size_t)r*sdsz],sdsz*8);
      const vector<double>& Fk=ket_T(sp.l,q);
      memcpy(&w.FkT[ip*sksz],Fk.data(),sksz*8);
      double* gg=&w.geo[(size_t)ip*8];
      gg[0]=Pc[0]-C[0]; gg[1]=Pc[1]-C[1]; gg[2]=Pc[2]-C[2];
      gg[3]=p*q/(p+q);
      gg[4]=2.0*std::pow(M_PI,2.5)/(p*q*std::sqrt(p+q));
      gg[5]=a/p; gg[6]=b/p;
      gg[7]=cc*sp.c[kq];
    }
  }
  const int nbk=ncb*nck;
  w.V.resize(nbk); w.dA.resize(3*nbk); w.dB.resize(3*nbk); w.dC.resize(3*nbk);
  int rr=routec_grad3c_block(sa.l,sb.l,sp.l,0,np,w.FT.data(),w.FdT.data(),
                             w.FkT.data(),w.geo.data(),w.V.data(),
                             w.dA.data(),w.dB.data(),w.dC.data());
  if(rr!=0){ fprintf(stderr,"[routec-grad-gpu] missing class (%d%d|%d)\n",
                     sa.l,sb.l,sp.l); abort(); }
}

static void block2c(const Sub& sp,const Sub& sq,const double* coords,BlockWS& w){
  const Ctx::Dims& d=class_dims(sp.l,0,sq.l);
  const int ncb=d.ncb, nck=d.nck;
  const size_t sbsz=(size_t)d.nchb*ncb, sdsz=3*(size_t)d.nchb1*ncb,
               sksz=(size_t)d.nchk*nck;
  const int np=sp.np*sq.np;
  w.FT.resize((size_t)np*sbsz);
  w.FdT.assign((size_t)np*sdsz,0.0);
  w.FkT.resize((size_t)np*sksz); w.geo.resize((size_t)np*8);
  const double* Cp=coords+3*sp.atom;
  const double* Cq=coords+3*sq.atom;
  int ip=0;
  for(int pa=0;pa<sp.np;++pa){
    const double a=sp.e[pa];
    const vector<double>& Fb=ket_T(sp.l,a);
    for(int pb=0;pb<sq.np;++pb,++ip){
      const double b=sq.e[pb];
      memcpy(&w.FT [ip*sbsz],Fb.data(),sbsz*8);
      const vector<double>& Fk=ket_T(sq.l,b);
      memcpy(&w.FkT[ip*sksz],Fk.data(),sksz*8);
      double* gg=&w.geo[(size_t)ip*8];
      gg[0]=Cp[0]-Cq[0]; gg[1]=Cp[1]-Cq[1]; gg[2]=Cp[2]-Cq[2];
      gg[3]=a*b/(a+b);
      gg[4]=2.0*std::pow(M_PI,2.5)/(a*b*std::sqrt(a+b));
      gg[5]=1.0; gg[6]=0.0;
      gg[7]=sp.c[pa]*sq.c[pb];
    }
  }
  const int nbk=ncb*nck;
  w.V.resize(nbk); w.dA.resize(3*nbk); w.dB.resize(3*nbk); w.dC.resize(3*nbk);
  int rr=routec_grad3c_block(sp.l,0,sq.l,0,np,w.FT.data(),w.FdT.data(),
                             w.FkT.data(),w.geo.data(),w.V.data(),
                             w.dA.data(),w.dB.data(),w.dC.data());
  if(rr!=0){ fprintf(stderr,"[routec-grad-gpu] missing 2c class (%d0|%d)\n",
                     sp.l,sq.l); abort(); }
}

// parallel for over [0, n)
template<class F> static void par_for(int n,int nthreads,F&& body){
  nthreads=std::min(nthreads,std::max(1,n));
  if(nthreads<=1){ for(int i=0;i<n;++i) body(i,0); return; }
  std::atomic<int> next(0);
  vector<std::thread> ths;
  for(int t=0;t<nthreads;++t)
    ths.emplace_back([&,t](){
      for(;;){ int i=next.fetch_add(1); if(i>=n) break; body(i,t); } });
  for(auto& th:ths) th.join();
}

// pass 1: M = (ab|P), V2 = (P|Q)  (host; identical to the CPU dylib)
static void build_MV(const double* xyz,vector<double>& M,vector<double>& V2){
  const int nao=g.nbf, naux=g.naux;
  const size_t nn=(size_t)nao*nao;
  double h=0.0; for(int i=0;i<3*g.natm;++i) h+=xyz[i]*(i+0.5);
  if(!(h==g.last_xyz_hash)){
    for(auto& P:g.pcache){ P.npp=0; P.FT.clear(); P.FdT.clear(); }
    g.last_xyz_hash=h;
  }
  M.assign(nn*naux,0.0); V2.assign((size_t)naux*naux,0.0);
  const int npair=(int)g.pairs.size(), nasub=(int)g.asub.size();

  // ---- aux 2-center metric V2 first (independent of M): needed for the aux
  // Schwarz factor when the optional screen is on; numerically identical to the
  // old post-M order (disjoint writes). ----
  {
    vector<std::array<int,2>> apairs;
    for(int i=0;i<nasub;++i) for(int j=0;j<=i;++j) apairs.push_back({i,j});
    par_for((int)apairs.size(),g.nthreads,[&](int k2,int){
      static thread_local BlockWS w;
      const Sub& sp=g.asub[apairs[k2][0]];
      const Sub& sq=g.asub[apairs[k2][1]];
      block2c(sp,sq,xyz,w);
      const int ncb=ncart(sp.l), nck=ncart(sq.l);
      for(int caa=0;caa<ncb;++caa)
        for(int cqq=0;cqq<nck;++cqq){
          double v=w.V[(size_t)caa*nck+cqq];
          V2[(size_t)(sp.ao0+caa)*naux+(sq.ao0+cqq)]=v;
          V2[(size_t)(sq.ao0+cqq)*naux+(sp.ao0+caa)]=v;
        }
    });
  }

  // ---- optional rot-distance screen (ROUTEC_GRAD_TAU_S, default 0 = OFF =
  // bit-identical). Mirrors build_df's test bound = qS_pair * qS_aux * F0(eta R^2)
  // <= tau_s over (pair,aux), so the gradient's M/V is the derivative of the SAME
  // screened energy tensor AND build_MV skips negligible (pair,aux) blocks. The
  // aux Schwarz sqrt((P|P)) is exact (V2 diagonal); the pair factor is the max
  // contracted primitive overlap (a reasonable, empirically-validated bound).
  // NOTE: the gradient is ~1000x MORE sensitive to the screen than the energy
  // (derivatives amplify the dropped blocks), so use a TIGHTER tau here than for
  // the energy: (H2O)16 measured dG_max = 9.6e-5 at 1e-7 but 2.3e-7 at 1e-8, so
  // tau ~= 1e-8 keeps the gradient ~100x under the DF-gradient error at ~11%
  // build-time saving. ----
  static const double GTAU = getenv("ROUTEC_GRAD_TAU_S")?atof(getenv("ROUTEC_GRAD_TAU_S")):0.0;
  vector<double> qSaux, Ppair, pminp, qSpair;
  auto F0=[](double t)->double{ return t<1e-12?1.0:0.5*std::sqrt(M_PI/t)*std::erf(std::sqrt(t)); };
  if(GTAU>0.0){
    qSaux.assign(nasub,0.0);
    for(int s=0;s<nasub;++s){ const Sub& sp=g.asub[s]; double m=0;
      for(int c=0;c<ncart(sp.l);++c){ double d=V2[(size_t)(sp.ao0+c)*naux+(sp.ao0+c)]; if(d>m)m=d; }
      qSaux[s]=std::sqrt(m>0?m:0.0); }
    Ppair.assign(3*npair,0.0); pminp.assign(npair,0.0); qSpair.assign(npair,0.0);
    par_for(npair,g.nthreads,[&](int ip,int){
      const Sub& sa=g.sub[g.pairs[ip][0]]; const Sub& sb=g.sub[g.pairs[ip][1]];
      const double* A=&xyz[3*sa.atom]; const double* B=&xyz[3*sb.atom];
      double al=sa.e[0]; for(double e:sa.e) if(e<al)al=e;
      double be=sb.e[0]; for(double e:sb.e) if(e<be)be=e;
      double p=al+be; double r2=0; for(int x=0;x<3;++x){double d=A[x]-B[x];r2+=d*d;}
      for(int x=0;x<3;++x) Ppair[3*ip+x]=(al*A[x]+be*B[x])/p;
      pminp[ip]=p;
      double q=0; for(int i=0;i<sa.np;++i)for(int j=0;j<sb.np;++j){
        double a=sa.e[i],b=sb.e[j],pp=a+b,mu=a*b/pp;
        double K=std::fabs(sa.c[i]*sb.c[j])*std::pow(M_PI/pp,1.5)*std::exp(-mu*r2); if(K>q)q=K; }
      qSpair[ip]=q;
    });
  }

  par_for(npair,g.nthreads,[&](int ip,int){ pair_get(ip,xyz); });
  par_for(npair,g.nthreads,[&](int ip,int){
    static thread_local BlockWS w;
    const Sub& sa=g.sub[g.pairs[ip][0]];
    const Sub& sb=g.sub[g.pairs[ip][1]];
    for(int s=0;s<nasub;++s){
      const Sub& sp=g.asub[s];
      if(GTAU>0.0){
        double r2=0; for(int x=0;x<3;++x){double d=Ppair[3*ip+x]-xyz[3*sp.atom+x]; r2+=d*d;}
        double cmin=sp.e[0]; for(double e:sp.e) if(e<cmin)cmin=e;
        double eta=pminp[ip]*cmin/(pminp[ip]+cmin);
        if(qSpair[ip]*qSaux[s]*F0(eta*r2) <= GTAU) continue;   // rot-distance screen
      }
      block3c(ip,sa,sb,sp,xyz,w,false);
      const int nca=ncart(sa.l), ncbb=ncart(sb.l), nck=ncart(sp.l);
      for(int caa=0;caa<nca;++caa)
        for(int cbb=0;cbb<ncbb;++cbb){
          const double* src=&w.V[(size_t)(caa*ncbb+cbb)*nck];
          double* d1=&M[((size_t)(sa.ao0+caa)*nao+(sb.ao0+cbb))*naux+sp.ao0];
          for(int p2=0;p2<nck;++p2) d1[p2]=src[p2];
          if(sa.ao0!=sb.ao0){
            double* d2=&M[((size_t)(sb.ao0+cbb)*nao+(sa.ao0+caa))*naux+sp.ao0];
            for(int p2=0;p2<nck;++p2) d2[p2]=src[p2];
          }
        }
    }
  });
}

// ===================================================================== GPU
// pass-2 Gamma-contracted device path.

#define CUCHK(x) do { cudaError_t e_=(x); if(e_!=cudaSuccess){ \
  fprintf(stderr,"[routec-grad-gpu] CUDA error %s at %s:%d\n", \
          cudaGetErrorString(e_),__FILE__,__LINE__); return false; } } while(0)

struct GTask {
  long long ppb;       // bra prim base (index into BraPrim[])
  long long akb;       // ket prim base (index into KetPrim[])
  long long gam;       // offset into the gathered Gamma/gamma pool
  long long tfirst;    // first triple id of this task within its class
  int npp, npq;
  int atA, atB, atC;
  int pad;
  double Cx, Cy, Cz;   // ket centre
  double wpair;
};
struct BraPrim { double Pcx,Pcy,Pcz,p,ap,bp,cc; long long ft, fd; };
struct KetPrim { double q, ck; long long fk; };

#define DEF_GKERN(SUF, NCB, NCK) \
__global__ void __launch_bounds__(128) gk_##SUF(long ntrip, \
    const GTask* __restrict__ tasks, int ntask, \
    const double* __restrict__ FT, const double* __restrict__ FdT, \
    const double* __restrict__ FkT, const double* __restrict__ GAM, \
    const BraPrim* __restrict__ bra, const KetPrim* __restrict__ ket, \
    double* de, int natm) { \
  extern __shared__ double sde[]; \
  for (int i = threadIdx.x; i < 3*natm; i += blockDim.x) sde[i] = 0.0; \
  __syncthreads(); \
  const long qid = (long)blockIdx.x*blockDim.x + threadIdx.x; \
  if (qid < ntrip) { \
    int lo = 0, hi = ntask - 1; \
    while (lo < hi) { int mid = (lo + hi + 1) >> 1; \
      if (tasks[mid].tfirst <= qid) lo = mid; else hi = mid - 1; } \
    const GTask tk = tasks[lo]; \
    const long lidx = qid - tk.tfirst; \
    const BraPrim B = bra[tk.ppb + lidx / tk.npq]; \
    const KetPrim K = ket[tk.akb + lidx % tk.npq]; \
    const double Rx = B.Pcx - tk.Cx, Ry = B.Pcy - tk.Cy, Rz = B.Pcz - tk.Cz; \
    const double eta = B.p*K.q/(B.p+K.q); \
    const double pref = 34.986836655249725/(B.p*K.q*sqrt(B.p+K.q)); \
    double V[NCB*NCK], Pl[3*NCB*NCK], DR[3*NCB*NCK]; \
    if (Rx*Rx+Ry*Ry+Rz*Rz < 1e-16) \
      grad3c_lab_##SUF(FT+B.ft, FdT+B.fd, FkT+K.fk, Rx,Ry,Rz, eta, pref, V, Pl, DR); \
    else \
      grad3c_al_##SUF(FT+B.ft, FdT+B.fd, FkT+K.fk, Rx,Ry,Rz, eta, pref, V, Pl, DR); \
    const double cf = B.cc * K.ck * tk.wpair; \
    const double* Gm = GAM + tk.gam; \
    for (int w = 0; w < 3; ++w) { \
      double sA=0.0, sB=0.0, sC=0.0; \
      for (int i = 0; i < NCB*NCK; ++i) { \
        const double gv = Gm[i]; \
        const double pl = Pl[w*NCB*NCK+i], dr = DR[w*NCB*NCK+i]; \
        sA += (pl + B.ap*dr)*gv; \
        sB += (-pl + B.bp*dr)*gv; \
        sC -= dr*gv; } \
      atomicAdd(&sde[3*tk.atA+w], cf*sA); \
      atomicAdd(&sde[3*tk.atB+w], cf*sB); \
      atomicAdd(&sde[3*tk.atC+w], cf*sC); \
    } \
  } \
  __syncthreads(); \
  for (int i = threadIdx.x; i < 3*natm; i += blockDim.x) \
    if (sde[i] != 0.0) atomicAdd(&de[i], sde[i]); \
}

DEF_GKERN(00_0, 1, 1)
DEF_GKERN(00_1, 1, 3)
DEF_GKERN(00_2, 1, 6)
DEF_GKERN(00_3, 1, 10)
DEF_GKERN(00_4, 1, 15)
DEF_GKERN(10_0, 3, 1)
DEF_GKERN(10_1, 3, 3)
DEF_GKERN(10_2, 3, 6)
DEF_GKERN(10_3, 3, 10)
DEF_GKERN(10_4, 3, 15)
DEF_GKERN(11_0, 9, 1)
DEF_GKERN(11_1, 9, 3)
DEF_GKERN(11_2, 9, 6)
DEF_GKERN(11_3, 9, 10)
DEF_GKERN(11_4, 9, 15)
DEF_GKERN(20_0, 6, 1)
DEF_GKERN(20_1, 6, 3)
DEF_GKERN(20_2, 6, 6)
DEF_GKERN(20_3, 6, 10)
DEF_GKERN(20_4, 6, 15)
DEF_GKERN(21_0, 18, 1)
DEF_GKERN(21_1, 18, 3)
DEF_GKERN(21_2, 18, 6)
DEF_GKERN(21_3, 18, 10)
DEF_GKERN(21_4, 18, 15)
DEF_GKERN(22_0, 36, 1)
DEF_GKERN(22_1, 36, 3)
DEF_GKERN(22_2, 36, 6)
DEF_GKERN(22_3, 36, 10)
DEF_GKERN(22_4, 36, 15)
DEF_GKERN(30_0, 10, 1)
DEF_GKERN(30_1, 10, 3)
DEF_GKERN(30_2, 10, 6)
DEF_GKERN(30_3, 10, 10)
DEF_GKERN(30_4, 10, 15)
DEF_GKERN(40_0, 15, 1)
DEF_GKERN(40_1, 15, 3)
DEF_GKERN(40_2, 15, 6)
DEF_GKERN(40_3, 15, 10)
DEF_GKERN(40_4, 15, 15)

typedef void (*gk_launch_t)(long, const GTask*, int, const double*,
    const double*, const double*, const double*, const BraPrim*,
    const KetPrim*, double*, int, int, size_t);
#define DEF_GLNCH(SUF) \
static void gl_##SUF(long nt, const GTask* tk, int ntask, const double* ft, \
    const double* fdt, const double* fkt, const double* gm, \
    const BraPrim* bp_, const KetPrim* kp_, double* de, int natm, \
    int bs, size_t shmem) { \
  gk_##SUF<<<(unsigned)((nt + bs - 1) / bs), bs, shmem>>>( \
      nt, tk, ntask, ft, fdt, fkt, gm, bp_, kp_, de, natm); \
}
DEF_GLNCH(00_0) DEF_GLNCH(00_1) DEF_GLNCH(00_2) DEF_GLNCH(00_3) DEF_GLNCH(00_4)
DEF_GLNCH(10_0) DEF_GLNCH(10_1) DEF_GLNCH(10_2) DEF_GLNCH(10_3) DEF_GLNCH(10_4)
DEF_GLNCH(11_0) DEF_GLNCH(11_1) DEF_GLNCH(11_2) DEF_GLNCH(11_3) DEF_GLNCH(11_4)
DEF_GLNCH(20_0) DEF_GLNCH(20_1) DEF_GLNCH(20_2) DEF_GLNCH(20_3) DEF_GLNCH(20_4)
DEF_GLNCH(21_0) DEF_GLNCH(21_1) DEF_GLNCH(21_2) DEF_GLNCH(21_3) DEF_GLNCH(21_4)
DEF_GLNCH(22_0) DEF_GLNCH(22_1) DEF_GLNCH(22_2) DEF_GLNCH(22_3) DEF_GLNCH(22_4)
DEF_GLNCH(30_0) DEF_GLNCH(30_1) DEF_GLNCH(30_2) DEF_GLNCH(30_3) DEF_GLNCH(30_4)
DEF_GLNCH(40_0) DEF_GLNCH(40_1) DEF_GLNCH(40_2) DEF_GLNCH(40_3) DEF_GLNCH(40_4)

struct GLEntry { int la, lb, lp; gk_launch_t fn; };
static const GLEntry g_glaunch[] = {
  {0,0,0,gl_00_0},{0,0,1,gl_00_1},{0,0,2,gl_00_2},{0,0,3,gl_00_3},{0,0,4,gl_00_4},
  {1,0,0,gl_10_0},{1,0,1,gl_10_1},{1,0,2,gl_10_2},{1,0,3,gl_10_3},{1,0,4,gl_10_4},
  {1,1,0,gl_11_0},{1,1,1,gl_11_1},{1,1,2,gl_11_2},{1,1,3,gl_11_3},{1,1,4,gl_11_4},
  {2,0,0,gl_20_0},{2,0,1,gl_20_1},{2,0,2,gl_20_2},{2,0,3,gl_20_3},{2,0,4,gl_20_4},
  {2,1,0,gl_21_0},{2,1,1,gl_21_1},{2,1,2,gl_21_2},{2,1,3,gl_21_3},{2,1,4,gl_21_4},
  {2,2,0,gl_22_0},{2,2,1,gl_22_1},{2,2,2,gl_22_2},{2,2,3,gl_22_3},{2,2,4,gl_22_4},
  {3,0,0,gl_30_0},{3,0,1,gl_30_1},{3,0,2,gl_30_2},{3,0,3,gl_30_3},{3,0,4,gl_30_4},
  {4,0,0,gl_40_0},{4,0,1,gl_40_1},{4,0,2,gl_40_2},{4,0,3,gl_40_3},{4,0,4,gl_40_4},
};

// gather metadata kept host-side, parallel to the task lists
struct GGather { int src; int a0,b0,p0; int nca,ncb,nck; long long gam; };

struct GpuCtx {
  bool tried=false, ready=false;
  int bs=128;
  // per-class task groups
  struct CG { int la,lb,lp,ncb,nck; long ntrip=0;
              vector<GTask> tasks; GTask* d_tasks=nullptr;
              gk_launch_t launch=nullptr; };
  vector<CG> cg;
  vector<GGather> gath;          // all tasks, gather order
  // bra prim staging (rebuilt per geometry)
  vector<BraPrim> bra; BraPrim* d_bra=nullptr;
  vector<long long> pair_ppb;    // BraPrim base per 3c pair
  vector<long long> aux_ppb;     // BraPrim base per aux subshell (2c bra)
  vector<KetPrim> ket; KetPrim* d_ket=nullptr;
  vector<long long> aux_kb;      // KetPrim base per aux subshell
  // pools
  vector<long long> pair_ft, pair_fd;   // pool offsets per 3c pair
  std::map<std::pair<int,long long>,long long> ketoff; // (l,expbits)->ft pool
  long long ftpool_sz=0, fdpool_sz=0, zfd_sz=0;
  double *d_ft=nullptr, *d_fd=nullptr;
  vector<double> h_ft, h_fd;
  long long gampool_sz=0; double* d_gam=nullptr; vector<double> h_gam;
  double* d_de=nullptr;
  double t_build=0, t_gather=0, t_upload=0, t_kern=0;
  double up_gb=0;
};
static GpuCtx gpu;

static GpuCtx::CG& gpu_cg(int la,int lb,int lp){
  for(auto& c:gpu.cg) if(c.la==la&&c.lb==lb&&c.lp==lp) return c;
  gpu.cg.push_back({});
  GpuCtx::CG& c=gpu.cg.back();
  c.la=la; c.lb=lb; c.lp=lp;
  const Ctx::Dims& d=class_dims(la,lb,lp);
  c.ncb=d.ncb; c.nck=d.nck;
  for(auto& e:g_glaunch)
    if(e.la==la&&e.lb==lb&&e.lp==lp){ c.launch=e.fn; break; }
  if(!c.launch){ fprintf(stderr,"[routec-grad-gpu] no launcher (%d%d|%d)\n",
                         la,lb,lp); abort(); }
  return c;
}

// one-time skeleton build: task lists, ket prim table, pool offsets
static bool gpu_init(){
  if(gpu.tried) return gpu.ready;
  gpu.tried=true;
  const char* dis=getenv("OQP_ROUTEC_GRAD_GPU");
  if(dis&&atoi(dis)==0){ fprintf(stderr,"[routec-grad-gpu] GPU disabled by env\n");
    return false; }
  int ndev=0;
  if(cudaGetDeviceCount(&ndev)!=cudaSuccess||ndev<1){
    fprintf(stderr,"[routec-grad-gpu] no CUDA device -- CPU fallback\n");
    return false; }
  const char* bse=getenv("OQP_ROUTEC_GRAD_BS");
  if(bse) gpu.bs=std::max(32,std::min(128,atoi(bse)));
  double t0=now_s();
  const int nasub=(int)g.asub.size(), npair=(int)g.pairs.size();

  // ---- FT/FdT pool layout.  FT pool head = all ket tensors (used as 2c
  // bra AND as 3c/2c ket); then per-3c-pair blocks.  FdT pool head = one
  // zero block (2c tasks); then per-3c-pair blocks.
  long long ftoff=0;
  for(auto& kv:g.ketcache){ gpu.ketoff[kv.first]=ftoff; ftoff+=(long long)kv.second.size(); }
  const long long ket_region=ftoff;
  long long zmax=0;
  gpu.pair_ft.resize(npair); gpu.pair_fd.resize(npair);
  for(auto& sq:g.asub){
    const Ctx::Dims& d=class_dims(sq.l,0,0);
    zmax=std::max(zmax,(long long)3*d.nchb1*d.ncb);
  }
  long long fdoff=zmax;   // zero block at head
  for(int ip=0;ip<npair;++ip){
    const Sub& sa=g.sub[g.pairs[ip][0]];
    const Sub& sb=g.sub[g.pairs[ip][1]];
    const Ctx::Dims& d=class_dims(sa.l,sb.l,0);
    const long long sbsz=(long long)d.nchb*d.ncb, sdsz=3ll*d.nchb1*d.ncb;
    gpu.pair_ft[ip]=ftoff; gpu.pair_fd[ip]=fdoff;
    ftoff+=(long long)sa.np*sb.np*sbsz;
    fdoff+=(long long)sa.np*sb.np*sdsz;
  }
  gpu.ftpool_sz=ftoff; gpu.fdpool_sz=fdoff; gpu.zfd_sz=zmax;

  // ---- ket prim table (geometry independent except via nothing)
  gpu.aux_kb.resize(nasub);
  for(int s=0;s<nasub;++s){
    const Sub& sp=g.asub[s];
    gpu.aux_kb[s]=(long long)gpu.ket.size();
    for(int k2=0;k2<sp.np;++k2){
      long long bits; memcpy(&bits,&sp.e[k2],8);
      KetPrim kp; kp.q=sp.e[k2]; kp.ck=sp.c[k2];
      kp.fk=gpu.ketoff.at(std::make_pair(sp.l,bits));
      gpu.ket.push_back(kp);
    }
  }
  // ---- bra prim bases: 3c pairs then 2c (aux subshells)
  long long nbra=0;
  gpu.pair_ppb.resize(npair); gpu.aux_ppb.resize(nasub);
  for(int ip=0;ip<npair;++ip){
    gpu.pair_ppb[ip]=nbra;
    nbra+=(long long)g.sub[g.pairs[ip][0]].np*g.sub[g.pairs[ip][1]].np;
  }
  for(int s=0;s<nasub;++s){ gpu.aux_ppb[s]=nbra; nbra+=g.asub[s].np; }
  gpu.bra.resize(nbra);

  // ---- tasks: 3c in CPU pass-2 order (pairs x aux subshells), then 2c
  long long gam_off=0;
  for(int ip=0;ip<npair;++ip){
    const Sub& sa=g.sub[g.pairs[ip][0]];
    const Sub& sb=g.sub[g.pairs[ip][1]];
    const double wp=(g.pairs[ip][0]==g.pairs[ip][1])?1.0:2.0;
    for(int s=0;s<nasub;++s){
      const Sub& sp=g.asub[s];
      GpuCtx::CG& c=gpu_cg(sa.l,sb.l,sp.l);
      GTask t; memset(&t,0,sizeof(t));
      t.ppb=gpu.pair_ppb[ip]; t.npp=sa.np*sb.np;
      t.akb=gpu.aux_kb[s];    t.npq=sp.np;
      t.gam=gam_off; t.tfirst=c.ntrip;
      t.atA=sa.atom; t.atB=sb.atom; t.atC=sp.atom;
      t.wpair=wp;     // Cx..Cz filled per geometry
      c.tasks.push_back(t);
      c.ntrip+=(long long)t.npp*t.npq;
      gpu.gath.push_back({0, sa.ao0, sb.ao0, sp.ao0,
                          ncart(sa.l), ncart(sb.l), ncart(sp.l), gam_off});
      gam_off+=(long long)ncart(sa.l)*ncart(sb.l)*ncart(sp.l);
    }
  }
  for(int i=0;i<nasub;++i)
    for(int j=0;j<=i;++j){
      const Sub& sp=g.asub[i];
      const Sub& sq=g.asub[j];
      const double wp=(i==j)?1.0:2.0;
      GpuCtx::CG& c=gpu_cg(sp.l,0,sq.l);
      GTask t; memset(&t,0,sizeof(t));
      t.ppb=gpu.aux_ppb[i]; t.npp=sp.np;
      t.akb=gpu.aux_kb[j];  t.npq=sq.np;
      t.gam=gam_off; t.tfirst=c.ntrip;
      t.atA=sp.atom; t.atB=sp.atom; t.atC=sq.atom;
      t.wpair=wp;
      c.tasks.push_back(t);
      c.ntrip+=(long long)t.npp*t.npq;
      gpu.gath.push_back({1, sp.ao0, 0, sq.ao0,
                          ncart(sp.l), 1, ncart(sq.l), gam_off});
      gam_off+=(long long)ncart(sp.l)*ncart(sq.l);
    }
  gpu.gampool_sz=gam_off;

  // ---- device allocations (sizes are geometry-independent)
  auto mb=[](long long x){ return (double)x*8.0/1048576.0; };
  long long ntrip_tot=0, ntask_tot=0;
  for(auto& c:gpu.cg){ ntrip_tot+=c.ntrip; ntask_tot+=(long long)c.tasks.size(); }
  fprintf(stderr,"[routec-grad-gpu] GPU init: %lld tasks, %lld prim triples, "
          "pools FT %.1f MB FdT %.1f MB Gam %.1f MB bra %lld ket %zu\n",
          ntask_tot,ntrip_tot,mb(gpu.ftpool_sz),mb(gpu.fdpool_sz),
          mb(gpu.gampool_sz),nbra,gpu.ket.size());
  CUCHK(cudaMalloc(&gpu.d_ft, gpu.ftpool_sz*8));
  CUCHK(cudaMalloc(&gpu.d_fd, gpu.fdpool_sz*8));
  CUCHK(cudaMalloc(&gpu.d_gam, gpu.gampool_sz*8));
  CUCHK(cudaMalloc(&gpu.d_bra, gpu.bra.size()*sizeof(BraPrim)));
  CUCHK(cudaMalloc(&gpu.d_ket, gpu.ket.size()*sizeof(KetPrim)));
  CUCHK(cudaMalloc(&gpu.d_de, 3*g.natm*8));
  for(auto& c:gpu.cg){
    CUCHK(cudaMalloc(&c.d_tasks, c.tasks.size()*sizeof(GTask)));
    CUCHK(cudaMemcpy(c.d_tasks, c.tasks.data(), c.tasks.size()*sizeof(GTask),
                     cudaMemcpyHostToDevice));
  }
  CUCHK(cudaMemcpy(gpu.d_ket, gpu.ket.data(), gpu.ket.size()*sizeof(KetPrim),
                   cudaMemcpyHostToDevice));
  gpu.h_ft.resize(gpu.ftpool_sz);
  gpu.h_fd.resize(gpu.fdpool_sz);
  gpu.h_gam.resize(gpu.gampool_sz);
  // ket region of the FT pool is geometry-independent: fill staging now
  std::fill(gpu.h_fd.begin(),gpu.h_fd.begin()+gpu.zfd_sz,0.0);
  for(auto& kv:g.ketcache){
    long long off=gpu.ketoff.at(kv.first);
    memcpy(gpu.h_ft.data()+off,kv.second.data(),kv.second.size()*8);
  }
  gpu.t_build=now_s()-t0;
  gpu.ready=true;
  return true;
}

// per-geometry: rebuild BraPrim + FT/FdT staging and upload
static bool gpu_geometry(const double* xyz){
  const int nasub=(int)g.asub.size(), npair=(int)g.pairs.size();
  double t0=now_s();
  par_for(npair,g.nthreads,[&](int ip,int){
    const Sub& sa=g.sub[g.pairs[ip][0]];
    const Sub& sb=g.sub[g.pairs[ip][1]];
    const Ctx::Dims& d=class_dims(sa.l,sb.l,0);
    const long long sbsz=(long long)d.nchb*d.ncb, sdsz=3ll*d.nchb1*d.ncb;
    const Ctx::PairT& P=pair_get(ip,xyz);
    memcpy(gpu.h_ft.data()+gpu.pair_ft[ip],P.FT.data(),P.FT.size()*8);
    memcpy(gpu.h_fd.data()+gpu.pair_fd[ip],P.FdT.data(),P.FdT.size()*8);
    const double* A=xyz+3*sa.atom;
    const double* B=xyz+3*sb.atom;
    BraPrim* bp_=gpu.bra.data()+gpu.pair_ppb[ip];
    int r=0;
    for(int pa=0;pa<sa.np;++pa)
      for(int pb=0;pb<sb.np;++pb,++r){
        const double a=sa.e[pa], b=sb.e[pb], p=a+b;
        BraPrim& q=bp_[r];
        q.Pcx=(a*A[0]+b*B[0])/p; q.Pcy=(a*A[1]+b*B[1])/p;
        q.Pcz=(a*A[2]+b*B[2])/p;
        q.p=p; q.ap=a/p; q.bp=b/p; q.cc=sa.c[pa]*sb.c[pb];
        q.ft=gpu.pair_ft[ip]+r*sbsz; q.fd=gpu.pair_fd[ip]+r*sdsz;
      }
  });
  for(int s=0;s<nasub;++s){
    const Sub& sp=g.asub[s];
    const double* C=xyz+3*sp.atom;
    BraPrim* bp_=gpu.bra.data()+gpu.aux_ppb[s];
    for(int pa=0;pa<sp.np;++pa){
      long long bits; memcpy(&bits,&sp.e[pa],8);
      BraPrim& q=bp_[pa];
      q.Pcx=C[0]; q.Pcy=C[1]; q.Pcz=C[2];
      q.p=sp.e[pa]; q.ap=1.0; q.bp=0.0; q.cc=sp.c[pa];
      q.ft=gpu.ketoff.at(std::make_pair(sp.l,bits)); q.fd=0;
    }
  }
  // fill task ket centres (depend on geometry through atC)
  for(auto& c:gpu.cg){
    bool any=false;
    for(auto& t:c.tasks){
      const double* C=xyz+3*t.atC;
      if(t.Cx!=C[0]||t.Cy!=C[1]||t.Cz!=C[2]){ t.Cx=C[0]; t.Cy=C[1]; t.Cz=C[2];
        any=true; }
    }
    if(any)
      if(cudaMemcpy(c.d_tasks,c.tasks.data(),c.tasks.size()*sizeof(GTask),
                    cudaMemcpyHostToDevice)!=cudaSuccess) return false;
  }
  double t1=now_s();
  CUCHK(cudaMemcpy(gpu.d_ft,gpu.h_ft.data(),gpu.ftpool_sz*8,
                   cudaMemcpyHostToDevice));
  CUCHK(cudaMemcpy(gpu.d_fd,gpu.h_fd.data(),gpu.fdpool_sz*8,
                   cudaMemcpyHostToDevice));
  CUCHK(cudaMemcpy(gpu.d_bra,gpu.bra.data(),gpu.bra.size()*sizeof(BraPrim),
                   cudaMemcpyHostToDevice));
  double t2=now_s();
  gpu.t_build+=t1-t0; gpu.t_upload+=t2-t1;
  gpu.up_gb+=(gpu.ftpool_sz+gpu.fdpool_sz)*8e-9
             +gpu.bra.size()*sizeof(BraPrim)*1e-9;
  return true;
}

// the GPU pass 2
static bool contract_pass2_gpu(const double* xyz,const double* Gam,
                               const double* gam,double* de_out){
  if(!gpu_init()) return false;
  gpu.t_gather=0; gpu.t_upload=0; gpu.t_kern=0; gpu.up_gb=0; gpu.t_build=0;
  if(!gpu_geometry(xyz)) return false;
  const int nao=g.nbf, naux=g.naux, natm=g.natm;
  // gather Gamma/gamma blocks (host, parallel)
  double t0=now_s();
  const int ngt=(int)gpu.gath.size();
  par_for(ngt,g.nthreads,[&](int it,int){
    const GGather& G2=gpu.gath[it];
    double* dst=gpu.h_gam.data()+G2.gam;
    if(G2.src==0){
      for(int ca=0;ca<G2.nca;++ca)
        for(int cb=0;cb<G2.ncb;++cb){
          const double* src=&Gam[((size_t)(G2.a0+ca)*nao+(G2.b0+cb))*naux+G2.p0];
          for(int p2=0;p2<G2.nck;++p2) *dst++=src[p2];
        }
    } else {
      for(int ca=0;ca<G2.nca;++ca){
        const double* src=&gam[(size_t)(G2.a0+ca)*naux+G2.p0];
        for(int p2=0;p2<G2.nck;++p2) *dst++=src[p2];
      }
    }
  });
  double t1=now_s();
  if(cudaMemcpy(gpu.d_gam,gpu.h_gam.data(),gpu.gampool_sz*8,
                cudaMemcpyHostToDevice)!=cudaSuccess) return false;
  CUCHK(cudaMemset(gpu.d_de,0,3*natm*8));
  double t2=now_s();
  const size_t shmem=(size_t)3*natm*8;
  for(auto& c:gpu.cg){
    if(!c.ntrip) continue;
    c.launch(c.ntrip,c.d_tasks,(int)c.tasks.size(),gpu.d_ft,gpu.d_fd,
             gpu.d_ft /*ket tensors live at the FT pool head*/,
             gpu.d_gam,gpu.d_bra,gpu.d_ket,gpu.d_de,natm,gpu.bs,shmem);
  }
  cudaError_t le=cudaDeviceSynchronize();
  if(le!=cudaSuccess){
    fprintf(stderr,"[routec-grad-gpu] kernel failure: %s\n",
            cudaGetErrorString(le));
    return false;
  }
  double t3=now_s();
  CUCHK(cudaMemcpy(de_out,gpu.d_de,3*natm*8,cudaMemcpyDeviceToHost));
  gpu.t_gather=t1-t0; gpu.t_upload+= t2-t1; gpu.t_kern=t3-t2;
  gpu.up_gb+=gpu.gampool_sz*8e-9;
  return true;
}

// ===================================================== pass-2 entry point
// CPU fallback = the validated CPU dylib contraction, kept verbatim.
static void contract_pass2_cpu(const double* xyz,const double* Gam,
                               const double* gam,double* de_out){
  const int nao=g.nbf, naux=g.naux, natm=g.natm;
  const int npair=(int)g.pairs.size(), nasub=(int)g.asub.size();
  vector<vector<double>> departs(g.nthreads,vector<double>(3*(size_t)natm,0.0));
  par_for(npair,g.nthreads,[&](int ip,int t){
    static thread_local BlockWS w;
    double* dep=departs[t].data();
    const Sub& sa=g.sub[g.pairs[ip][0]];
    const Sub& sb=g.sub[g.pairs[ip][1]];
    const double wpair=(g.pairs[ip][0]==g.pairs[ip][1])?1.0:2.0;
    const int nca=ncart(sa.l), ncbb=ncart(sb.l);
    for(int s=0;s<nasub;++s){
      const Sub& sp=g.asub[s];
      block3c(ip,sa,sb,sp,xyz,w,false);
      const int nck=ncart(sp.l);
      const int nbk=nca*ncbb*nck;
      for(int w2=0;w2<3;++w2){
        double sA=0.0,sB=0.0,sC=0.0;
        const double *pa2=&w.dA[(size_t)w2*nbk], *pb2=&w.dB[(size_t)w2*nbk],
                     *pc2=&w.dC[(size_t)w2*nbk];
        int idx=0;
        for(int caa=0;caa<nca;++caa)
          for(int cbb=0;cbb<ncbb;++cbb){
            const double* gm=&Gam[((size_t)(sa.ao0+caa)*nao+(sb.ao0+cbb))*naux+sp.ao0];
            for(int p2=0;p2<nck;++p2,++idx){
              const double gv=gm[p2];
              sA+=pa2[idx]*gv; sB+=pb2[idx]*gv; sC+=pc2[idx]*gv;
            }
          }
        dep[3*sa.atom+w2]+=wpair*sA;
        dep[3*sb.atom+w2]+=wpair*sB;
        dep[3*sp.atom+w2]+=wpair*sC;
      }
    }
  });
  {
    vector<std::array<int,2>> apairs;
    for(int i=0;i<nasub;++i) for(int j=0;j<=i;++j) apairs.push_back({i,j});
    par_for((int)apairs.size(),g.nthreads,[&](int k2,int t){
      static thread_local BlockWS w;
      double* dep=departs[t].data();
      const Sub& sp=g.asub[apairs[k2][0]];
      const Sub& sq=g.asub[apairs[k2][1]];
      const double wpair=(apairs[k2][0]==apairs[k2][1])?1.0:2.0;
      block2c(sp,sq,xyz,w);
      const int ncb=ncart(sp.l), nck=ncart(sq.l), nbk=ncb*nck;
      for(int w2=0;w2<3;++w2){
        double sP=0.0,sQ=0.0;
        const double *pa2=&w.dA[(size_t)w2*nbk], *pb2=&w.dB[(size_t)w2*nbk],
                     *pc2=&w.dC[(size_t)w2*nbk];
        int idx=0;
        for(int caa=0;caa<ncb;++caa){
          const double* gm=&gam[(size_t)(sp.ao0+caa)*naux+sq.ao0];
          for(int q2=0;q2<nck;++q2,++idx){
            sP+=(pa2[idx]+pb2[idx])*gm[q2];
            sQ+=pc2[idx]*gm[q2];
          }
        }
        dep[3*sp.atom+w2]+=wpair*sP;
        dep[3*sq.atom+w2]+=wpair*sQ;
      }
    });
  }
  for(int i=0;i<3*natm;++i){
    double s=0.0;
    for(int t=0;t<g.nthreads;++t) s+=departs[t][i];
    de_out[i]=s;
  }
}

static void contract_pass2(const double* xyz,const double* Gam,
                           const double* gam,double* de_out){
  double t0=now_s();
  if(contract_pass2_gpu(xyz,Gam,gam,de_out)){
    if(g.verbose)
      fprintf(stderr,"[routec-grad-gpu] pass2 GPU: build %.3f s | gather "
              "%.3f s | upload %.3f s (%.2f GB) | kernels %.3f s | total "
              "%.3f s\n",gpu.t_build,gpu.t_gather,gpu.t_upload,gpu.up_gb,
              gpu.t_kern,now_s()-t0);
    return;
  }
  fprintf(stderr,"[routec-grad-gpu] pass2 falling back to CPU\n");
  contract_pass2_cpu(xyz,Gam,gam,de_out);
}

// ============================================== dE2e/dx assembly (HOST,
// verbatim from the CPU dylib except the pass-2 call)
static int grad2e(const double* dpk,const double* xyz,double* de_out,
                  int nbf,int natm,double hs,double cs){
  if(nbf!=g.nbf||natm!=g.natm){
    fprintf(stderr,"[routec-grad-gpu] dim mismatch: seam nbf=%d natm=%d vs input "
            "%d/%d\n",nbf,natm,g.nbf,g.natm);
    return 1;
  }
  const int nao=g.nbf, naux=g.naux;
  const size_t nn=(size_t)nao*nao;
  double t0=now_s();

  vector<double> D(nn,0.0);
  { vector<double> Do(nn);
    size_t t=0;
    for(int i=0;i<nao;++i)
      for(int j=0;j<=i;++j,++t){ Do[(size_t)i*nao+j]=dpk[t];
                                 Do[(size_t)j*nao+i]=dpk[t]; }
    for(int i=0;i<nao;++i)
      for(int j=0;j<nao;++j)
        D[(size_t)g.qmap[i]*nao+g.qmap[j]]=g.smap[i]*g.smap[j]*Do[(size_t)i*nao+j];
  }

  double t1=now_s();
  vector<double> M, V2;
  build_MV(xyz,M,V2);
  double t2=now_s();

  lapack_int_t n_l=naux, info_l=0, one_l=1;
  vector<double> Vch(V2);
  dpotrf_("L",&n_l,Vch.data(),&n_l,&info_l);
  if(info_l!=0){
    fprintf(stderr,"[routec-grad-gpu] dpotrf failed (info=%d)\n",(int)info_l);
    return 2;
  }
  vector<double> dvec(naux), c(naux);
  cblas_dgemv(CblasRowMajor,CblasTrans,(int)nn,naux,1.0,M.data(),naux,
              D.data(),1,0.0,dvec.data(),1);
  c=dvec;
  dpotrs_("L",&n_l,&one_l,Vch.data(),&n_l,c.data(),&n_l,&info_l);
  vector<double> G(M);
  { lapack_int_t nrhs=(lapack_int_t)nn;
    dpotrs_("L",&n_l,&nrhs,Vch.data(),&n_l,G.data(),&n_l,&info_l);
    if(info_l!=0){ fprintf(stderr,"[routec-grad-gpu] dpotrs failed\n"); return 2; }
  }
  double t3=now_s();

  // ---- Gamma/gamma on DEVICE (default; OQP_ROUTEC_GRAD_GAMMA_GPU=0 for the
  // host fallback).  The three big contractions (G := D.G per-slab batched;
  // gam = sum over the permuted product; Gamma = D.G + rank-1) are ~2 nao^3
  // naux flops -- decisive on the GPU at production sizes.
  bool gam_dev = true;
  if(const char* e=getenv("OQP_ROUTEC_GRAD_GAMMA_GPU")) gam_dev = atoi(e)!=0;
  vector<double> gam((size_t)naux*naux,0.0);
  bool gam_done=false;
  if(gam_dev){
    const size_t nnau=(size_t)nao*nao*naux, dnn=(size_t)nao*nao;
    cublasHandle_t bh=nullptr; cublasStatus_t bs=cublasCreate(&bh);
    double *dD=nullptr,*dG=nullptr,*dX=nullptr,*dc=nullptr,*dgm=nullptr;
    auto ok=[&](cudaError_t r){return r==cudaSuccess;};
    if(bs==CUBLAS_STATUS_SUCCESS
       && ok(cudaMalloc(&dD,dnn*8)) && ok(cudaMalloc(&dG,nnau*8))
       && ok(cudaMalloc(&dX,nnau*8)) && ok(cudaMalloc(&dc,(size_t)naux*8))
       && ok(cudaMalloc(&dgm,(size_t)naux*naux*8))
       && ok(cudaMemcpy(dD,D.data(),dnn*8,cudaMemcpyHostToDevice))
       && ok(cudaMemcpy(dG,G.data(),nnau*8,cudaMemcpyHostToDevice))
       && ok(cudaMemcpy(dc,c.data(),(size_t)naux*8,cudaMemcpyHostToDevice))){
      const double one=1.0,zero=0.0;
      const long long sAB=(long long)nao*naux;
      // (1) per-slab G_a := D * G_a   (row-major) == col-major C_a = G_a^c * D
      bs=cublasDgemmStridedBatched(bh,CUBLAS_OP_N,CUBLAS_OP_N,naux,nao,nao,
            &one,dG,naux,sAB,dD,nao,0,&zero,dX,naux,sAB,nao);
      if(bs==CUBLAS_STATUS_SUCCESS){
        std::swap(dG,dX);                        // dG = updated G
        // (2) permute (a,c,P)<-(c,a,P) then gam = P^T . G over the (a,c) index
        { long long tot=(long long)nnau; int TB=256; long long nb=(tot+TB-1)/TB;
          gradgpu_perm102<<<(unsigned)nb,TB>>>(dG,dX,nao,naux,tot);
        }
        if(cudaGetLastError()==cudaSuccess){
          // row-major gam(naux x naux) = Xflat^T (N x naux)^T . Gflat (N x naux)
          bs=cublasDgemm(bh,CUBLAS_OP_N,CUBLAS_OP_T,naux,naux,(int)dnn,
                &one,dG,naux,dX,naux,&zero,dgm,naux);
          if(bs==CUBLAS_STATUS_SUCCESS
             && ok(cudaMemcpy(gam.data(),dgm,(size_t)naux*naux*8,cudaMemcpyDeviceToHost))){
            // (3) Gamma = cs*D(x)c - 0.5*hs * (D . G)   (big GEMM + elementwise)
            bs=cublasDgemm(bh,CUBLAS_OP_N,CUBLAS_OP_N,(int)sAB,nao,nao,
                  &one,dG,(int)sAB,dD,nao,&zero,dX,(int)sAB);
            if(bs==CUBLAS_STATUS_SUCCESS){
              long long tot=(long long)nnau; int TB=256; long long nb=(tot+TB-1)/TB;
              gradgpu_gamma_elem<<<(unsigned)nb,TB>>>(dX,dD,dc,nao,naux,hs,cs,tot);
              if(cudaGetLastError()==cudaSuccess
                 && ok(cudaMemcpy(M.data(),dX,nnau*8,cudaMemcpyDeviceToHost))){
                for(size_t i=0;i<gam.size();++i) gam[i]*=0.25*hs;
                for(int p2=0;p2<naux;++p2)
                  for(int q2=0;q2<naux;++q2)
                    gam[(size_t)p2*naux+q2]-=0.5*cs*c[p2]*c[q2];
                gam_done=true;
              }
            }
          }
        }
      }
    }
    if(dD)cudaFree(dD); if(dG)cudaFree(dG); if(dX)cudaFree(dX);
    if(dc)cudaFree(dc); if(dgm)cudaFree(dgm);
    if(bh)cublasDestroy(bh);
    if(!gam_done) fprintf(stderr,"[routec-grad-gpu] device Gamma failed; host fallback\n");
  }
  if(!gam_done){
  {
    vector<vector<double>> tmps(g.nthreads,vector<double>((size_t)nao*naux));
    par_for(nao,g.nthreads,[&](int a2,int t){
      double* tmp=tmps[t].data();
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nao,naux,nao,
                  1.0,D.data(),nao,G.data()+(size_t)a2*nao*naux,naux,
                  0.0,tmp,naux);
      memcpy(G.data()+(size_t)a2*nao*naux,tmp,(size_t)nao*naux*8);
    });
  }
  }
  double E_J=0.0;
  for(int p2=0;p2<naux;++p2) E_J+=dvec[p2]*c[p2];
  E_J*=0.5;
  double E_K=0.0;
  if(g.verbose && !gam_done){
    vector<double> MD((size_t)nao*naux);
    for(int c2=0;c2<nao;++c2){
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nao,naux,nao,
                  1.0,D.data(),nao,M.data()+(size_t)c2*nao*naux,naux,
                  0.0,MD.data(),naux);
      for(int a2=0;a2<nao;++a2)
        E_K+=cblas_ddot(naux,G.data()+((size_t)a2*nao+c2)*naux,1,
                        MD.data()+(size_t)a2*naux,1);
    }
    E_K*=-0.25;
    fprintf(stderr,"[routec-grad-gpu] E_J=%.12f E_K=%.12f E2e(cs,hs)=%.12f\n",
            E_J,E_K,cs*E_J+hs*E_K);
  }
  if(!gam_done){
  {
    vector<vector<double>> strips(g.nthreads,vector<double>((size_t)nao*naux));
    vector<vector<double>> parts(g.nthreads,vector<double>((size_t)naux*naux,0.0));
    par_for(nao,g.nthreads,[&](int a2,int t){
      double* strip=strips[t].data();
      for(int c2=0;c2<nao;++c2)
        memcpy(strip+(size_t)c2*naux,
               G.data()+((size_t)c2*nao+a2)*naux,(size_t)naux*8);
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,naux,naux,nao,
                  1.0,G.data()+(size_t)a2*nao*naux,naux,strip,naux,
                  1.0,parts[t].data(),naux);
    });
    for(int t=0;t<g.nthreads;++t)
      for(size_t i=0;i<gam.size();++i) gam[i]+=parts[t][i];
    for(size_t i=0;i<gam.size();++i) gam[i]*=0.25*hs;
    for(int p2=0;p2<naux;++p2)
      for(int q2=0;q2<naux;++q2)
        gam[(size_t)p2*naux+q2]-=0.5*cs*c[p2]*c[q2];
  }
  {
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nao,(int)((size_t)nao*naux),
                nao,1.0,D.data(),nao,G.data(),(int)((size_t)nao*naux),
                0.0,M.data(),(int)((size_t)nao*naux));
    par_for(nao,g.nthreads,[&](int a2,int){
      for(int b2=0;b2<nao;++b2){
        double dab=cs*D[(size_t)a2*nao+b2];
        double* row=&M[((size_t)a2*nao+b2)*naux];
        for(int p2=0;p2<naux;++p2) row[p2]=dab*c[p2]-0.5*hs*row[p2];
      }
    });
  }
  }
  vector<double>().swap(G);
  const vector<double>& Gam=M;
  double t4=now_s();

  contract_pass2(xyz,Gam.data(),gam.data(),de_out);
  double t5=now_s();
  if(g.verbose)
    fprintf(stderr,"[routec-grad-gpu] timings: pairs %.3f s | pass1(int) %.3f s | "
            "chol+G %.3f s | Gamma/gamma %.3f s | pass2(deriv) %.3f s | "
            "total %.3f s\n",t1-t0,t2-t1,t3-t2,t4-t3,t5-t4,t5-t0);
  return 0;
}

// ===================================================================== G-5
struct GTerm { const double* A; const double* B; double w; };

static int grad2e_mrsf(const vector<vector<double>>& dens,
                       const GTerm* Jt,int nJ,const GTerm* Kt,int nK,
                       const double* xyz,double* de_out){
  const int nao=g.nbf, naux=g.naux;
  const size_t nn=(size_t)nao*nao, gsz=nn*(size_t)naux;
  (void)dens;
  double t0=now_s();
  vector<double> M, V2;
  build_MV(xyz,M,V2);
  double t1=now_s();

  lapack_int_t n_l=naux, info_l=0, one_l=1;
  vector<double> Vch(V2);
  dpotrf_("L",&n_l,Vch.data(),&n_l,&info_l);
  if(info_l!=0){
    fprintf(stderr,"[routec-grad-gpu] mrsf dpotrf failed (info=%d)\n",(int)info_l);
    return 2;
  }
  vector<double> Gam(gsz,0.0), gam((size_t)naux*naux,0.0);
  double E2e=0.0;

  // G := V^{-1} M is needed by the K terms on either path, and the device
  // path derives the J projections from it directly: G^T a = V^{-1} M^T a,
  // so the per-term triangular solves disappear.
  vector<double> G(M);
  { lapack_int_t nrhs=(lapack_int_t)nn;
    dpotrs_("L",&n_l,&nrhs,Vch.data(),&n_l,G.data(),&n_l,&info_l);
    if(info_l!=0){ fprintf(stderr,"[routec-grad-gpu] mrsf dpotrs failed\n"); return 2; }
  }

  // ---- MRSF Gamma/gamma on DEVICE (default; OQP_ROUTEC_GRAD_GAMMA_GPU=0
  // for the host fallback).  Each K term is 3 GEMM-shaped contractions of
  // O(nao^3 naux); with up to 9 K terms this is the dominant gradient cost.
  bool mdev = true;
  if(const char* e=getenv("OQP_ROUTEC_GRAD_GAMMA_GPU")) mdev = atoi(e)!=0;
  bool mdev_done=false;
  double t2=now_s();
  if(mdev){
    cublasHandle_t bh=nullptr; cublasStatus_t bs=cublasCreate(&bh);
    double *dG=nullptr,*dGam=nullptr,*dT1=nullptr,*dT2=nullptr,*dU=nullptr,
           *dA=nullptr,*dB=nullptr,*dcA=nullptr,*dcB=nullptr,*dT=nullptr,*dgm=nullptr;
    auto ok=[&](cudaError_t r){return r==cudaSuccess;};
    const double one=1.0,zero=0.0;
    const long long sAB=(long long)nao*naux;
    const int ncol=(int)((size_t)nao*naux);
    const int TB=256;
    bool okall = bs==CUBLAS_STATUS_SUCCESS
      && ok(cudaMalloc(&dG,gsz*8)) && ok(cudaMalloc(&dGam,gsz*8))
      && ok(cudaMalloc(&dT1,gsz*8)) && ok(cudaMalloc(&dT2,gsz*8))
      && ok(cudaMalloc(&dU,gsz*8))
      && ok(cudaMalloc(&dA,nn*8)) && ok(cudaMalloc(&dB,nn*8))
      && ok(cudaMalloc(&dcA,(size_t)naux*8)) && ok(cudaMalloc(&dcB,(size_t)naux*8))
      && ok(cudaMalloc(&dT,(size_t)naux*naux*8)) && ok(cudaMalloc(&dgm,(size_t)naux*naux*8))
      && ok(cudaMemcpy(dG,G.data(),gsz*8,cudaMemcpyHostToDevice))
      && ok(cudaMemset(dGam,0,gsz*8)) && ok(cudaMemset(dgm,0,(size_t)naux*naux*8));
    // J terms: solved projections via G^T, rank-1 Gamma on device; the tiny
    // naux^2 gam/E2e pieces committed on host only after full success.
    vector<vector<double>> cAs(nJ,vector<double>(naux)), cBs(nJ,vector<double>(naux));
    for(int it=0;okall&&it<nJ;++it){
      const double* A=Jt[it].A; const double* B=Jt[it].B; const double w=Jt[it].w;
      okall = ok(cudaMemcpy(dA,A,nn*8,cudaMemcpyHostToDevice))
           && ok(cudaMemcpy(dB,B,nn*8,cudaMemcpyHostToDevice))
           && cublasDgemv(bh,CUBLAS_OP_N,naux,(int)nn,&one,dG,naux,dA,1,&zero,dcA,1)
                ==CUBLAS_STATUS_SUCCESS
           && cublasDgemv(bh,CUBLAS_OP_N,naux,(int)nn,&one,dG,naux,dB,1,&zero,dcB,1)
                ==CUBLAS_STATUS_SUCCESS;
      if(okall){
        long long tot=(long long)gsz, nb=(tot+TB-1)/TB;
        gradgpu_mrsf_rank1<<<(unsigned)nb,TB>>>(dGam,dA,dB,dcA,dcB,w,naux,tot);
        okall = cudaGetLastError()==cudaSuccess
             && ok(cudaMemcpy(cAs[it].data(),dcA,(size_t)naux*8,cudaMemcpyDeviceToHost))
             && ok(cudaMemcpy(cBs[it].data(),dcB,(size_t)naux*8,cudaMemcpyDeviceToHost));
      }
    }
    t2=now_s();
    // K terms
    for(int it=0;okall&&it<nK;++it){
      const double* A=Kt[it].A; const double* B=Kt[it].B; const double w=Kt[it].w;
      auto is_sym=[&](const double* X){
        for(int i=0;i<nao;++i)
          for(int j=0;j<i;++j)
            if(std::fabs(X[(size_t)i*nao+j]-X[(size_t)j*nao+i])>1e-14*
               (1.0+std::fabs(X[(size_t)i*nao+j]))) return false;
        return true;
      };
      const bool symA=is_sym(A), symB=is_sym(B);
      okall = ok(cudaMemcpy(dA,A,nn*8,cudaMemcpyHostToDevice))
           && ok(cudaMemcpy(dB,B,nn*8,cudaMemcpyHostToDevice))
           // tmp1[c,a,P] = sum_b B[a,b] G[c,b,P]  (per-slab batched GEMM)
           && cublasDgemmStridedBatched(bh,CUBLAS_OP_N,CUBLAS_OP_N,naux,nao,nao,
                &one,dG,naux,sAB,dB,nao,0,&zero,dT1,naux,sAB,nao)
                ==CUBLAS_STATUS_SUCCESS;
      if(okall&&!symB)   // tmp2: same with B transposed
        okall = cublasDgemmStridedBatched(bh,CUBLAS_OP_N,CUBLAS_OP_T,naux,nao,nao,
                  &one,dG,naux,sAB,dB,nao,0,&zero,dT2,naux,sAB,nao)
                  ==CUBLAS_STATUS_SUCCESS;
      if(okall){
        double* t2p = symB ? dT1 : dT2;
        if(symA&&symB){
          const double a2w=2.0*w;
          okall = cublasDgemm(bh,CUBLAS_OP_N,CUBLAS_OP_N,ncol,nao,nao,
                    &a2w,dT1,ncol,dA,nao,&one,dGam,ncol)==CUBLAS_STATUS_SUCCESS;
        } else {
          okall = cublasDgemm(bh,CUBLAS_OP_N,CUBLAS_OP_N,ncol,nao,nao,
                    &w,dT1,ncol,dA,nao,&one,dGam,ncol)==CUBLAS_STATUS_SUCCESS
               && cublasDgemm(bh,CUBLAS_OP_N,CUBLAS_OP_T,ncol,nao,nao,
                    &w,t2p,ncol,dA,nao,&one,dGam,ncol)==CUBLAS_STATUS_SUCCESS;
        }
      }
      if(okall){
        // U[a,b,P] = sum_k A[k,a] G[k,b,P];  T[P,Q] = sum_t U[t,P] tmp1[t,Q]
        okall = cublasDgemm(bh,CUBLAS_OP_N,CUBLAS_OP_T,ncol,nao,nao,
                  &one,dG,ncol,dA,nao,&zero,dU,ncol)==CUBLAS_STATUS_SUCCESS
             && cublasDgemm(bh,CUBLAS_OP_N,CUBLAS_OP_T,naux,naux,(int)nn,
                  &one,dT1,naux,dU,naux,&zero,dT,naux)==CUBLAS_STATUS_SUCCESS;
      }
      if(okall){
        long long tot=(long long)naux*naux, nb=(tot+TB-1)/TB;
        gradgpu_mrsf_gamsym<<<(unsigned)nb,TB>>>(dgm,dT,w,naux,tot);
        okall = cudaGetLastError()==cudaSuccess;
      }
    }
    if(okall){
      vector<double> gamK((size_t)naux*naux);
      okall = ok(cudaMemcpy(Gam.data(),dGam,gsz*8,cudaMemcpyDeviceToHost))
           && ok(cudaMemcpy(gamK.data(),dgm,(size_t)naux*naux*8,cudaMemcpyDeviceToHost));
      if(okall){
        vector<double> Vc(naux);
        for(int it=0;it<nJ;++it){
          const double w=Jt[it].w;
          const double* cA=cAs[it].data(); const double* cB=cBs[it].data();
          for(int p=0;p<naux;++p)
            for(int q=0;q<naux;++q)
              gam[(size_t)p*naux+q]-=0.5*w*(cA[p]*cB[q]+cB[p]*cA[q]);
          // E2e = w * cA^T V cB  (cA,cB are the solved projections)
          cblas_dgemv(CblasRowMajor,CblasNoTrans,naux,naux,1.0,V2.data(),naux,
                      cB,1,0.0,Vc.data(),1);
          double aB=0.0;
          for(int p=0;p<naux;++p) aB+=cA[p]*Vc[p];
          E2e+=w*aB;
        }
        for(size_t i=0;i<gam.size();++i) gam[i]+=gamK[i];
        mdev_done=true;
      }
    }
    if(dG)cudaFree(dG); if(dGam)cudaFree(dGam); if(dT1)cudaFree(dT1);
    if(dT2)cudaFree(dT2); if(dU)cudaFree(dU); if(dA)cudaFree(dA);
    if(dB)cudaFree(dB); if(dcA)cudaFree(dcA); if(dcB)cudaFree(dcB);
    if(dT)cudaFree(dT); if(dgm)cudaFree(dgm);
    if(bh)cublasDestroy(bh);
    if(!mdev_done){
      fprintf(stderr,"[routec-grad-gpu] mrsf device Gamma failed; host fallback\n");
      std::fill(Gam.begin(),Gam.end(),0.0);
      std::fill(gam.begin(),gam.end(),0.0);
      E2e=0.0;
    }
  }

  if(!mdev_done){
  for(int it=0;it<nJ;++it){
    const double* A=Jt[it].A; const double* B=Jt[it].B; const double w=Jt[it].w;
    vector<double> cA(naux), cB(naux);
    cblas_dgemv(CblasRowMajor,CblasTrans,(int)nn,naux,1.0,M.data(),naux,
                A,1,0.0,cA.data(),1);
    cblas_dgemv(CblasRowMajor,CblasTrans,(int)nn,naux,1.0,M.data(),naux,
                B,1,0.0,cB.data(),1);
    double aB=0.0;
    vector<double> aA(cA);
    dpotrs_("L",&n_l,&one_l,Vch.data(),&n_l,cA.data(),&n_l,&info_l);
    dpotrs_("L",&n_l,&one_l,Vch.data(),&n_l,cB.data(),&n_l,&info_l);
    for(int p=0;p<naux;++p) aB+=aA[p]*cB[p];
    E2e+=w*aB;
    par_for(nao,g.nthreads,[&](int a2,int){
      for(int b2=0;b2<nao;++b2){
        double* row=&Gam[((size_t)a2*nao+b2)*naux];
        const double av=w*A[(size_t)a2*nao+b2], bv=w*B[(size_t)a2*nao+b2];
        for(int p=0;p<naux;++p) row[p]+=av*cB[p]+bv*cA[p];
      }
    });
    for(int p=0;p<naux;++p)
      for(int q=0;q<naux;++q)
        gam[(size_t)p*naux+q]-=0.5*w*(cA[p]*cB[q]+cB[p]*cA[q]);
  }
  t2=now_s();

  {
    vector<double> tmp1(gsz), tmp2(gsz), U(gsz);
    vector<vector<double>> parts(g.nthreads,vector<double>((size_t)naux*naux,0.0));
    for(int it=0;it<nK;++it){
      const double* A=Kt[it].A; const double* B=Kt[it].B; const double w=Kt[it].w;
      auto is_sym=[&](const double* X){
        for(int i=0;i<nao;++i)
          for(int j=0;j<i;++j)
            if(std::fabs(X[(size_t)i*nao+j]-X[(size_t)j*nao+i])>1e-14*
               (1.0+std::fabs(X[(size_t)i*nao+j]))) return false;
        return true;
      };
      const bool symA=is_sym(A), symB=is_sym(B);
      par_for(nao,g.nthreads,[&](int c2,int){
        const double* Gc=G.data()+(size_t)c2*nao*naux;
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nao,naux,nao,
                    1.0,B,nao,Gc,naux,0.0,tmp1.data()+(size_t)c2*nao*naux,naux);
        if(!symB)
          cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,nao,naux,nao,
                      1.0,B,nao,Gc,naux,0.0,tmp2.data()+(size_t)c2*nao*naux,naux);
      });
      const double* t2p = symB ? tmp1.data() : tmp2.data();
      const int ncol=(int)((size_t)nao*naux);
      if(symA&&symB){
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nao,ncol,nao,
                    2.0*w,A,nao,tmp1.data(),ncol,1.0,Gam.data(),ncol);
      } else {
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nao,ncol,nao,
                    w,A,nao,tmp1.data(),ncol,1.0,Gam.data(),ncol);
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,nao,ncol,nao,
                    w,A,nao,t2p,ncol,1.0,Gam.data(),ncol);
      }
      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,nao,ncol,nao,
                  1.0,A,nao,G.data(),ncol,0.0,U.data(),ncol);
      std::fill(parts.begin(),parts.end(),
                vector<double>((size_t)naux*naux,0.0));
      par_for(nao,g.nthreads,[&](int k2,int t){
        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,naux,naux,nao,
                    1.0,U.data()+(size_t)k2*nao*naux,naux,
                    tmp1.data()+(size_t)k2*nao*naux,naux,
                    1.0,parts[t].data(),naux);
      });
      vector<double> T((size_t)naux*naux,0.0);
      for(int t=0;t<g.nthreads;++t)
        for(size_t i=0;i<T.size();++i) T[i]+=parts[t][i];
      for(int p=0;p<naux;++p)
        for(int q=0;q<naux;++q)
          gam[(size_t)p*naux+q]-=0.5*w*(T[(size_t)p*naux+q]
                                        +T[(size_t)q*naux+p]);
    }
  }
  }
  double t3=now_s();

  par_for(nao,g.nthreads,[&](int a2,int){
    for(int b2=0;b2<a2;++b2){
      double* r1=&Gam[((size_t)a2*nao+b2)*naux];
      double* r2=&Gam[((size_t)b2*nao+a2)*naux];
      for(int p=0;p<naux;++p){
        const double s=0.5*(r1[p]+r2[p]);
        r1[p]=s; r2[p]=s;
      }
    }
  });
  for(int p=0;p<naux;++p)
    for(int q=0;q<p;++q){
      const double s=0.5*(gam[(size_t)p*naux+q]+gam[(size_t)q*naux+p]);
      gam[(size_t)p*naux+q]=s; gam[(size_t)q*naux+p]=s;
    }

  contract_pass2(xyz,Gam.data(),gam.data(),de_out);
  double t4=now_s();
  if(g.verbose)
    fprintf(stderr,"[routec-grad-gpu] mrsf timings: pass1 %.3f s | chol+J %.3f s "
            "| K-terms %.3f s | pass2 %.3f s | total %.3f s | E2e(J) %.12f\n",
            t1-t0,t2-t1,t3-t2,t4-t3,t4-t0,E2e);
  return 0;
}

static void map_oqp_full(const double* a_f,double* out){
  const int nao=g.nbf;
  for(int i=0;i<nao;++i)
    for(int j=0;j<nao;++j)
      out[(size_t)g.qmap[i]*nao+g.qmap[j]] =
          g.smap[i]*g.smap[j]*a_f[(size_t)j*nao+i];
}

} // namespace oqpgrad

extern "C" {

// Hand the gradient deck (shells/aux/AO map, the routec_grad_inp text) to the
// library in memory -- the tagarray path, so no OQP_ROUTEC_GRAD_INP file is
// needed. Call once before the first gradient; len<=0 or null clears it (revert
// to the file). This is what the OpenQP seam calls with the deck it pulls from
// the tagarray tag OpenQP populates from get_basis().
void routec_grad_set_deck(const char* buf, int len){
  if(buf && len>0) oqpgrad::g_deck_buf.assign(buf,(size_t)len);
  else             oqpgrad::g_deck_buf.clear();
}

void routec_grad2(const double* d, const double* xyz, double* de,
                  const int* nbf, const int* natm, const double* hfscale,
                  const double* coulscale, int* info){
  using namespace oqpgrad;
  *info=1;
  std::lock_guard<std::mutex> lk(g_mtx);
  if(!load_input()) return;
  double tt0=now_s();
  int rc=grad2e(d,xyz,de,*nbf,*natm,*hfscale,*coulscale);
  if(rc==0){
    *info=0;
    fprintf(stderr,"[routec-grad-gpu] dE2e/dx done in %.3f s (hfscale=%.3f)\n",
            now_s()-tt0,*hfscale);
  }
}

void routec_grad2_mrsf(const double* d, const double* p, const double* spc,
                       const double* xyz, double* de, const int* nbf,
                       const int* natm, const double* hfscale,
                       const double* hfscale2, const double* coulscale,
                       const double* spcscale, const int* mrst, int* info){
  using namespace oqpgrad;
  *info=1;
  std::lock_guard<std::mutex> lk(g_mtx);
  if(!load_input()) return;
  if(*nbf!=g.nbf||*natm!=g.natm){
    fprintf(stderr,"[routec-grad-gpu] mrsf dim mismatch: seam nbf=%d natm=%d vs "
            "input %d/%d\n",*nbf,*natm,g.nbf,g.natm);
    return;
  }
  if(*mrst!=1&&*mrst!=3){
    fprintf(stderr,"[routec-grad-gpu] mrsf: unsupported mrst=%d\n",*mrst);
    return;
  }
  double tt0=now_s();
  const int nao=g.nbf;
  const size_t nn=(size_t)nao*nao;

  vector<vector<double>> dn;
  dn.reserve(20);
  for(int s=0;s<2;++s){
    dn.emplace_back(nn);
    map_oqp_full(d+(size_t)s*nn,dn.back().data());
  }
  for(int s=0;s<2;++s){
    dn.emplace_back(nn);
    map_oqp_full(p+(size_t)s*nn,dn.back().data());
  }
  for(int m=0;m<7;++m){
    dn.emplace_back(nn,0.0);
    double* out=dn.back().data();
    for(int i=0;i<nao;++i)
      for(int j=0;j<nao;++j)
        out[(size_t)g.qmap[i]*nao+g.qmap[j]] =
            g.smap[i]*g.smap[j]*spc[m+7*((size_t)i+(size_t)nao*j)];
  }
  auto comb=[&](int ia,int ib,double wa,double wb){
    dn.emplace_back(nn);
    double* o=dn.back().data();
    const double *A=dn[ia].data(), *B=dn[ib].data();
    for(size_t k=0;k<nn;++k) o[k]=wa*A[k]+wb*B[k];
    return (int)dn.size()-1;
  };
  auto transp=[&](int ia){
    dn.emplace_back(nn);
    double* o=dn.back().data();
    const double* A=dn[ia].data();
    for(int i=0;i<nao;++i)
      for(int j=0;j<nao;++j) o[(size_t)i*nao+j]=A[(size_t)j*nao+i];
    return (int)dn.size()-1;
  };
  auto symm=[&](int ia){
    dn.emplace_back(nn);
    double* o=dn.back().data();
    const double* A=dn[ia].data();
    for(int i=0;i<nao;++i)
      for(int j=0;j<nao;++j)
        o[(size_t)i*nao+j]=A[(size_t)i*nao+j]+A[(size_t)j*nao+i];
    return (int)dn.size()-1;
  };
  const int iD1=comb(0,1,1.0,1.0), iD2=comb(0,1,1.0,-1.0);
  const int iP1=comb(2,3,1.0,1.0), iP2=comb(2,3,1.0,-1.0);
  const int iball=10, iballT=transp(10);
  const int ico12=9, ico12T=transp(9);
  const int io21v=8, io21vT=transp(8);
  const int isbco1=symm(6), isbco2=symm(7), isbo2v=symm(4), isbo1v=symm(5);

  const double xc=*hfscale, xc2=*hfscale2, cs=*coulscale;
  const double sgnk=(*mrst==3)?-1.0:1.0;
  const double q1=sgnk*spcscale[0], q2=sgnk*spcscale[1], q3=sgnk*spcscale[2];

  std::vector<GTerm> Jt, Kt;
  auto Dp=[&](int i){ return dn[i].data(); };
  Jt.push_back({Dp(iD1),Dp(iD1),cs/2});
  Jt.push_back({Dp(iP1),Dp(iD1),cs/2});
  Jt.push_back({Dp(iD1),Dp(iP1),cs/2});
  if(q3!=0.0){
    Jt.push_back({Dp(isbco1),Dp(isbo2v),q3/2});
    Jt.push_back({Dp(isbco2),Dp(isbo1v),q3/2});
  }
  (void)iballT;
  Kt.push_back({Dp(iD1),Dp(iD1),-xc/4});
  Kt.push_back({Dp(iP1),Dp(iD1),-xc/2});
  Kt.push_back({Dp(iD2),Dp(iD2),-xc/4});
  Kt.push_back({Dp(iP2),Dp(iD2),-xc/2});
  Kt.push_back({Dp(iball),Dp(iball),-xc2});
  if(q1!=0.0) Kt.push_back({Dp(ico12),Dp(ico12T),q1});
  if(q2!=0.0) Kt.push_back({Dp(io21v),Dp(io21vT),q2});
  if(q3!=0.0){
    Kt.push_back({Dp(6),Dp(4),-2*q3});
    Kt.push_back({Dp(7),Dp(5),-2*q3});
  }

  int rc=grad2e_mrsf(dn,Jt.data(),(int)Jt.size(),Kt.data(),(int)Kt.size(),
                     xyz,de);
  if(rc==0){
    *info=0;
    fprintf(stderr,"[routec-grad-gpu] MRSF dE2e/dx done in %.3f s "
            "(xc=%.3f xc2=%.3f mrst=%d, %d J + %d K terms)\n",
            now_s()-tt0,xc,xc2,*mrst,(int)Jt.size(),(int)Kt.size());
  }
}

}  // extern "C"
