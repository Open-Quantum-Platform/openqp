// OpenQP G-4/G-5 DF 2e-GRADIENT supplier (GRADIENTS.md gates G-4, G-5):
// external bind(C) symbols called from the patched hf_2e_grad /
// mrsf_2e_grad (openqp feat/routec-dfjk-bridge, $OQP_ROUTEC_GRAD_LIB seam).
//
//   routec_grad2(d, xyz, de, nbf, natm, hfscale, coulscale, info)
//     d    : packed lower-tri TOTAL density (occupation-2), OpenQP AO frame
//     xyz  : (3, natm) atom coordinates, Bohr (column per atom)
//     de   : (3, natm) OUT, dE_2e/dx = cs*dE_J + hs*dE_K of the DF energy
//     info : 0 ok; nonzero -> caller falls back to the native grd2 path
//
//   routec_grad2_mrsf(d, p, spc, xyz, de, nbf, natm, hfscale, hfscale2,
//                     coulscale, spcscale, mrst, info)        (G-5)
//     d/p  : (nbf,nbf,2) RAW alpha/beta SCF and relaxed difference
//            densities (full square, Fortran order); spc: (7,nbf,nbf)
//            MRSF response densities (slot 7 = ball); spcscale: (3)
//            spc_coco/ovov/coov; mrst: 1|3. Generalized DF Gamma/gamma
//            assembly over the phase-A(ii)-verified bilinear term list
//            (see routec_grad2_mrsf below).
//
// dE_2e/dx = sum_abP Gamma_abP d(ab|P)/dx + sum_PQ gamma_PQ d(P|Q)/dx with
// the G-3-validated intermediates (test/grad_dfrhf_gate.py):
//   c = V^{-1} d,  G = V^{-1} M   (ONE Cholesky solve; the cart-aux metric
//   has cond(V) ~ 8e11 -- NEVER form inv(V) or the explicit A = M'V^{-1}M)
//   Gamma = cs * D (x) c  - hs/2 * D G_P D
//   gamma = -cs/2 * c (x) c + hs/4 * Tr[G_P D G_Q D]
// Derivative 3c/2c integrals come from the EMITTED G-3 kernels
// (test/routec_grad_gen.inc, classes {ss,ps,pp,ds,dp,ds,dd,fs,gs} x aux
// l=0..4) -- no pyscf at runtime. Pair tensors F / F' are built HERE in C++
// (the GRAD_DERIV.md S2 E'-recursion port; layouts match the kernel ABI:
// FT [nchb][ncb], FdT [3][nchb1][ncb], FkT [nchk][nck], channel-major).
//
// Conventions: all integral work happens in the pyscf-cart frame of the
// exported molecule (the validated G-3 frame). The OpenQP density is mapped
// in by the q/s congruence (D_py[q_i,q_j] = s_i s_j D_oqp[i,j], densities
// contravariant -- stage2_openqp_dfjk.py); the output dE/dx is a derivative
// w.r.t. nuclear coordinates and is frame-invariant, so no back-map.
//
// Environment:
//   OQP_ROUTEC_GRAD_INP     molecule/basis/map file (written by the G-4
//                           harness; format below)  [required]
//   OQP_ROUTEC_TABLES       routec_tables.bin (RTC2)  [default repo path]
//   OQP_ROUTEC_GRAD_THREADS worker threads  [default: hw concurrency]
//   OQP_ROUTEC_GRAD_VERBOSE stage timings + energy diagnostics to stderr
//
// Input file format (text):
//   routec_grad_inp v1
//   natm <n>   nbf <nbf>   naux <naux>
//   nsub <n>                 # orbital subshells, pyscf-cart AO order
//   <l> <atom> <nprim>       # per subshell, then nprim lines "exp coef"
//   ...                      # coef = raw * gto_norm(l,exp) * f(l)  (libcint
//                            #        common factor; written by the harness)
//   nauxsub <n>              # aux subshells (cart!), same layout
//   ...
//   map                      # nbf lines: q[i] s[i]  (OQP AO i -> pyscf AO
//                            #   q[i], scale s[i] = nrm[q[i]])
//
// Build (macOS):
//   clang++ -O1 -std=c++17 -shared -fPIC -framework Accelerate \
//       test/routec_oqp_grad.cpp -o libroutec_oqp_grad.dylib
//   (NOT g++-15: GCC cannot compile Accelerate's NEON headers.)
// Build (linux): g++ -O2 -std=c++17 -shared -fPIC test/routec_oqp_grad.cpp \
//       -lopenblas -o libroutec_oqp_grad.so
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <array>
#include <map>
#include <atomic>
#include <thread>
#include <chrono>
#include <mutex>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
typedef __CLPK_integer lapack_int_t;
#else
#include <cblas.h>
extern "C" {
void dpotrf_(const char*, const int*, double*, const int*, int*);
void dpotrs_(const char*, const int*, const int*, const double*, const int*,
             double*, const int*, int*);
}
typedef int lapack_int_t;
#endif

#include "routec_grad_gen.inc"   // emitted G-3 kernels + extern-C drivers

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
static bool load_tables(const char* path, RTC& T){
  FILE* f=fopen(path,"rb"); if(!f) return false;
  char magic[4]; int km,kp;
  if(fread(magic,1,4,f)!=4||memcmp(magic,"RTC2",4)) {fclose(f);return false;}
  if(fread(&km,4,1,f)!=1||fread(&kp,4,1,f)!=1){fclose(f);return false;}
  T.kmax=km; T.kpair=kp;
  T.deg.resize(km+1);
  for(int k=0;k<=km;++k){ DegTab& D=T.deg[k];
    if(fread(&D.nmono,4,1,f)!=1||fread(&D.nch,4,1,f)!=1||
       fread(&D.nch0,4,1,f)!=1){fclose(f);return false;}
    D.ch.resize(D.nch);
    for(int j=0;j<D.nch;++j){int b[3];
      if(fread(b,4,3,f)!=3){fclose(f);return false;}
      D.ch[j]={b[0],b[1],b[2]};}
    D.C.resize((size_t)D.nch*D.nmono);
    if(fread(D.C.data(),8,D.C.size(),f)!=D.C.size()){fclose(f);return false;}
    D.ch0.resize(D.nch0);
    if(D.nch0&&fread(D.ch0.data(),4,D.nch0,f)!=(size_t)D.nch0){fclose(f);return false;}
    D.B0.resize((size_t)D.nmono*D.nch0);
    if(D.B0.size()&&fread(D.B0.data(),8,D.B0.size(),f)!=D.B0.size()){fclose(f);return false;}
  }
  // rest of the file (g0 COO, rotation refs, c2s) is not needed here: the
  // emitted kernels carry their own coupling literals.
  fclose(f); return true;
}

// ------------------------------------------- Hermite E full table (one axis)
// E^{ij}_t per the Python cartesian_eri.hermite_E recursion (exact match,
// identical to test_routec.cpp e_axis but keeps the whole (i,j) table).
struct ETab {
  int I=0, J=0, W=0;              // I=imax+1, J=jmax+1, W=imax+jmax+1
  vector<double> v;               // [(i*J+j)*W + t]
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
// FT layout [ch][comp] (channel-major TRANSPOSED -- the kernel ABI), comp
// = ia_idx*ncart(lb) + ib_idx in monos_of order; degrees 0..K concatenated
// in original channel order. Port of grad3c_prototype.pair_tensor.
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

// dF^{ab}/dA_w, w=x,y,z; degrees 0..K1+1; layout [w][ch][comp] (kernel ABI).
// dE^{ij}_t/dA_w = 2a E^{i+1,j}_t - i E^{i-1,j}_t - (a/p) E^{ij}_{t-1}
// (GRAD_DERIV.md S2; dF/dB = -dF/dA handled in the emitted driver).
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
  // dE rows cache per (axis, i, j): t = 0..i+j+1
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
  vector<int> qmap;            // OQP AO i -> pyscf AO q[i]
  vector<double> smap;         // s_i
  // per-class dims cache
  struct Dims { int nchb,nchb1,nchk,ncb,nck; };
  std::map<int,Dims> dims;
  // pair tensor cache (geometry-dependent: rebuilt per gradient call)
  // key = pair index in the (i>=j swapped) pair list
  struct PairT { vector<double> FT, FdT; int npp=0; };
  vector<PairT> pcache;
  vector<std::array<int,2>> pairs;        // (ia_sub, ib_sub) with l_a >= l_b
  std::map<std::pair<int,long long>,vector<double>> ketcache; // (l, exp bits)
  int nthreads=1;
  bool verbose=false;
  double last_xyz_hash=0.0/0.0;
};
static Ctx g;
static std::mutex g_mtx;

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
    fprintf(stderr,"[routec-grad] FATAL: no emitted class (%d%d|%d)\n",la,lb,lp);
    abort();
  }
  Ctx::Dims dd{d[0],d[1],d[2],d[3],d[4]};
  return g.dims.emplace(key,dd).first->second;
}

static const vector<double>& ket_T(int lp,double q);   // fwd

static bool load_input(){
  if(g.loaded) return g.ok;
  g.loaded=true; g.ok=false;
  const char* tabp=getenv("OQP_ROUTEC_TABLES");
  if(!tabp) tabp="/Volumes/External_Storage/claude/libintRot/tools/routec_tables.bin";
  if(!load_tables(tabp,g.T)){
    fprintf(stderr,"[routec-grad] cannot load tables %s\n",tabp); return false; }
  const char* inp=getenv("OQP_ROUTEC_GRAD_INP");
  if(!inp){ fprintf(stderr,"[routec-grad] OQP_ROUTEC_GRAD_INP unset\n"); return false; }
  FILE* f=fopen(inp,"r");
  if(!f){ fprintf(stderr,"[routec-grad] cannot open %s\n",inp); return false; }
  char buf[128], ver[32];
  int nbf_hdr=0, naux_hdr=0;
  if(fscanf(f,"%127s %31s",buf,ver)!=2||strcmp(buf,"routec_grad_inp")!=0){
    fprintf(stderr,"[routec-grad] bad header in %s\n",inp); fclose(f); return false; }
  if(fscanf(f,"%127s %d",buf,&g.natm)!=2||strcmp(buf,"natm")) {fclose(f);return false;}
  if(fscanf(f,"%127s %d",buf,&nbf_hdr)!=2||strcmp(buf,"nbf"))  {fclose(f);return false;}
  if(fscanf(f,"%127s %d",buf,&naux_hdr)!=2||strcmp(buf,"naux")){fclose(f);return false;}
  int ns=0,nas=0;
  if(!read_subs(f,"nsub",ns,g.sub,g.nbf)){fclose(f);return false;}
  if(!read_subs(f,"nauxsub",nas,g.asub,g.naux)){fclose(f);return false;}
  if(g.nbf!=nbf_hdr||g.naux!=naux_hdr){
    fprintf(stderr,"[routec-grad] header/subshell count mismatch (%d/%d vs %d/%d)\n",
            nbf_hdr,naux_hdr,g.nbf,g.naux); fclose(f); return false; }
  if(fscanf(f,"%127s",buf)!=1||strcmp(buf,"map")){fclose(f);return false;}
  g.qmap.resize(g.nbf); g.smap.resize(g.nbf);
  for(int i=0;i<g.nbf;++i)
    if(fscanf(f,"%d %lf",&g.qmap[i],&g.smap[i])!=2){fclose(f);return false;}
  fclose(f);
  for(auto& S:g.sub) if(S.l>2){
    fprintf(stderr,"[routec-grad] orbital l=%d unsupported (classes are s,p,d)\n",S.l);
    return false; }
  for(auto& S:g.asub) if(S.l>4){
    fprintf(stderr,"[routec-grad] aux l=%d unsupported (classes go to g)\n",S.l);
    return false; }
  // unordered subshell pair list, (sa,sb) ordered so l_a >= l_b
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
  // prebuild every class-dims entry and ket tensor the worker threads will
  // touch (both caches are std::map -- populate single-threaded HERE, then
  // they are read-only in the parallel sections)
  for(auto& pr:g.pairs)
    for(auto& sp:g.asub)
      (void)class_dims(g.sub[pr[0]].l,g.sub[pr[1]].l,sp.l);
  for(auto& sp:g.asub)
    for(auto& sq:g.asub)
      (void)class_dims(sp.l,0,sq.l);
  for(auto& sp:g.asub)
    for(int k2=0;k2<sp.np;++k2)
      (void)ket_T(sp.l,sp.e[k2]);
  fprintf(stderr,"[routec-grad] input loaded: natm=%d nbf=%d naux=%d "
          "(%zu orbital / %zu aux subshells), %d threads\n",
          g.natm,g.nbf,g.naux,g.sub.size(),g.asub.size(),g.nthreads);
  g.ok=true;
  return true;
}

// ket tensor (lp,0) at one centre: geometry-independent, cached per (l,exp).
// NOTE: the cache is fully populated in load_input() (single-threaded);
// lookups from worker threads are read-only afterwards.
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

// build (or fetch) the bra pair tensors for pair index ip at the current
// geometry; coords = (3,natm)
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

// one contracted 3c block (sa sb | sp): fills V[ncb*nck], dA/dB/dC[3*ncb*nck]
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
  if(rr!=0){ fprintf(stderr,"[routec-grad] missing class (%d%d|%d)\n",
                     sa.l,sb.l,sp.l); abort(); }
}

// one contracted 2c block (P|Q) via class (lp,0|lq) with FdT = 0:
// d/dC_P = dA+dB (= +D_R), d/dC_Q = dC (= -D_R)
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
  if(rr!=0){ fprintf(stderr,"[routec-grad] missing 2c class (%d0|%d)\n",
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

// pass 1: build M = (ab|P) (a,b,P C-order) and V2 = (P|Q); also refreshes
// the per-geometry pair cache.
static void build_MV(const double* xyz,vector<double>& M,vector<double>& V2){
  const int nao=g.nbf, naux=g.naux;
  const size_t nn=(size_t)nao*nao;
  // geometry changed -> drop the per-geometry pair cache
  double h=0.0; for(int i=0;i<3*g.natm;++i) h+=xyz[i]*(i+0.5);
  if(!(h==g.last_xyz_hash)){
    for(auto& P:g.pcache){ P.npp=0; P.FT.clear(); P.FdT.clear(); }
    g.last_xyz_hash=h;
  }
  M.assign(nn*naux,0.0); V2.assign((size_t)naux*naux,0.0);
  const int npair=(int)g.pairs.size(), nasub=(int)g.asub.size();
  // pair tensors first (parallel, lazily built inside too)
  par_for(npair,g.nthreads,[&](int ip,int){ pair_get(ip,xyz); });
  par_for(npair,g.nthreads,[&](int ip,int){
    static thread_local BlockWS w;
    const Sub& sa=g.sub[g.pairs[ip][0]];
    const Sub& sb=g.sub[g.pairs[ip][1]];
    for(int s=0;s<nasub;++s){
      const Sub& sp=g.asub[s];
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
  // 2c metric (small; serial in aux pairs but parallel over rows)
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
}

// pass 2: contract the derivative blocks with Gamma (nn x naux; MUST be
// symmetric in (a,b)) and gamma (naux x naux; symmetric).
static void contract_pass2(const double* xyz,const double* Gam,
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

// the full dE2e/dx assembly; returns 0 on success
static int grad2e(const double* dpk,const double* xyz,double* de_out,
                  int nbf,int natm,double hs,double cs){
  if(nbf!=g.nbf||natm!=g.natm){
    fprintf(stderr,"[routec-grad] dim mismatch: seam nbf=%d natm=%d vs input "
            "%d/%d\n",nbf,natm,g.nbf,g.natm);
    return 1;
  }
  const int nao=g.nbf, naux=g.naux;
  const size_t nn=(size_t)nao*nao;
  double t0=now_s();

  // density: packed OQP -> dense OQP -> pyscf frame (contravariant map)
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

  // ---- pass 1: M = (ab|P) (a,b,P C-order) and V2 = (P|Q)
  double t1=now_s();
  vector<double> M, V2;
  build_MV(xyz,M,V2);
  double t2=now_s();

  // ---- DF intermediates: ONE Cholesky factorization of the cart-aux metric
  lapack_int_t n_l=naux, info_l=0, one_l=1;
  vector<double> Vch(V2);
  dpotrf_("L",&n_l,Vch.data(),&n_l,&info_l);
  if(info_l!=0){
    fprintf(stderr,"[routec-grad] dpotrf failed (info=%d)\n",(int)info_l);
    return 2;
  }
  // dvec_P = sum_ab M_abP D_ab ; c = V^{-1} dvec
  vector<double> dvec(naux), c(naux);
  cblas_dgemv(CblasRowMajor,CblasTrans,(int)nn,naux,1.0,M.data(),naux,
              D.data(),1,0.0,dvec.data(),1);
  c=dvec;
  dpotrs_("L",&n_l,&one_l,Vch.data(),&n_l,c.data(),&n_l,&info_l);
  // G = V^{-1} M : flat M (nn, naux) row-major == (naux, nn) col-major
  vector<double> G(M);
  { lapack_int_t nrhs=(lapack_int_t)nn;
    dpotrs_("L",&n_l,&nrhs,Vch.data(),&n_l,G.data(),&n_l,&info_l);
    if(info_l!=0){ fprintf(stderr,"[routec-grad] dpotrs failed\n"); return 2; }
  }
  double t3=now_s();

  // GD_acP = sum_b G_abP D_bc (per-a GEMM, in place via temp)
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
  // diagnostics: E_J, E_K of the DF energy at this density
  double E_J=0.0;
  for(int p2=0;p2<naux;++p2) E_J+=dvec[p2]*c[p2];
  E_J*=0.5;
  double E_K=0.0;
  if(g.verbose){
    // E_K = -1/4 sum_acP GD_acP * MD_caP, MD_caP = sum_b M_cbP D_ba
    vector<double> MD((size_t)nao*naux);
    for(int c2=0;c2<nao;++c2){
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nao,naux,nao,
                  1.0,D.data(),nao,M.data()+(size_t)c2*nao*naux,naux,
                  0.0,MD.data(),naux);
      // MD[a][P] here = sum_b D_ab M_cbP = MD_caP (D symmetric)
      for(int a2=0;a2<nao;++a2)
        E_K+=cblas_ddot(naux,G.data()+((size_t)a2*nao+c2)*naux,1,
                        MD.data()+(size_t)a2*naux,1);
    }
    E_K*=-0.25;
    fprintf(stderr,"[routec-grad] E_J=%.12f E_K=%.12f E2e(cs,hs)=%.12f\n",
            E_J,E_K,cs*E_J+hs*E_K);
  }
  // gamma_PQ = -cs/2 c_P c_Q + hs/4 sum_ac GD_acP GD_caQ
  vector<double> gam((size_t)naux*naux,0.0);
  {
    vector<vector<double>> strips(g.nthreads,vector<double>((size_t)nao*naux));
    vector<vector<double>> parts(g.nthreads,vector<double>((size_t)naux*naux,0.0));
    par_for(nao,g.nthreads,[&](int a2,int t){
      double* strip=strips[t].data();           // strip[c][Q] = GD_caQ
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
  // Gamma_abP = cs D_ab c_P - hs/2 (D GD)_abP   (overwrite M)
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
  vector<double>().swap(G);   // free
  const vector<double>& Gam=M;
  double t4=now_s();

  // ---- pass 2: derivative blocks contracted with Gamma / gamma
  contract_pass2(xyz,Gam.data(),gam.data(),de_out);
  double t5=now_s();
  if(g.verbose)
    fprintf(stderr,"[routec-grad] timings: pairs %.3f s | pass1(int) %.3f s | "
            "chol+G %.3f s | Gamma/gamma %.3f s | pass2(deriv) %.3f s | "
            "total %.3f s\n",t1-t0,t2-t1,t3-t2,t4-t3,t5-t4,t5-t0);
  return 0;
}

// ===================================================================== G-5
// MRSF-TDDFT 2e gradient: generalized DF Gamma/gamma assembly over bilinear
// term lists (verified numpy reference: sessions/20260612_grad_g5/
// g5_probeA2_grad2pdm.py; term list from the (1/8)*sum' decomposition of
// grd2_mrsf_compute_data_t%get_density, exact-vs-native 2e term 2e-9):
//   E_J(A,B;w) = w A_ij (ij|kl) B_kl :
//     Gamma += w (A (x) cB + B (x) cA),  gamma -= w/2 (cA cB^T + cB cA^T)
//   E_K(A,B;w) = w A_ik (ij|kl) B_jl :
//     Gamma += w (A G_P B^T + A^T G_P B)
//     gamma -= w/2 (T + T^T),  T_PQ = G_ijP A_ik G_klQ B_jl
//   (cX = V^{-1}(M:X), G = V^{-1}M; one Cholesky as in grad2e)
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
    fprintf(stderr,"[routec-grad] mrsf dpotrf failed (info=%d)\n",(int)info_l);
    return 2;
  }
  vector<double> Gam(gsz,0.0), gam((size_t)naux*naux,0.0);
  double E2e=0.0;

  // ---- J terms
  for(int it=0;it<nJ;++it){
    const double* A=Jt[it].A; const double* B=Jt[it].B; const double w=Jt[it].w;
    vector<double> cA(naux), cB(naux);
    cblas_dgemv(CblasRowMajor,CblasTrans,(int)nn,naux,1.0,M.data(),naux,
                A,1,0.0,cA.data(),1);
    cblas_dgemv(CblasRowMajor,CblasTrans,(int)nn,naux,1.0,M.data(),naux,
                B,1,0.0,cB.data(),1);
    double aB=0.0;
    vector<double> aA(cA);                 // raw vectors before solve
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

  // ---- G = V^{-1} M
  vector<double> G(M);
  { lapack_int_t nrhs=(lapack_int_t)nn;
    dpotrs_("L",&n_l,&nrhs,Vch.data(),&n_l,G.data(),&n_l,&info_l);
    if(info_l!=0){ fprintf(stderr,"[routec-grad] mrsf dpotrs failed\n"); return 2; }
  }
  double t2=now_s();

  // ---- K terms
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
      // tmp1[c][b][P] = sum_d B[b][d] G[c][d][P]    (per-c GEMM)
      // tmp2[c][b][P] = sum_d B[d][b] G[c][d][P]    (== tmp1 if B symmetric)
      par_for(nao,g.nthreads,[&](int c2,int){
        const double* Gc=G.data()+(size_t)c2*nao*naux;
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nao,naux,nao,
                    1.0,B,nao,Gc,naux,0.0,tmp1.data()+(size_t)c2*nao*naux,naux);
        if(!symB)
          cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,nao,naux,nao,
                      1.0,B,nao,Gc,naux,0.0,tmp2.data()+(size_t)c2*nao*naux,naux);
      });
      const double* t2p = symB ? tmp1.data() : tmp2.data();
      // AGBt[a][b][P] = sum_c A[a][c] tmp1[c][b][P];  Gam += w*AGBt
      // AGB [a][b][P] = sum_c A[c][a] tmp2[c][b][P];  Gam += w*AGB
      // U   [k][j][P] = sum_i A[i][k] G[i][j][P]      (for T)
      const int ncol=(int)((size_t)nao*naux);
      if(symA&&symB){      // AGBt == AGB: one GEMM with doubled weight
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
      // E_K = w sum_abP M_abP (A G_P B^T)_ab = w sum: M : (A*tmp1)
      // (reuse U for T first, E via dot of M with the Gam increment is
      //  awkward after accumulation -- compute E directly:)
      // T_PQ = sum_{kj} U[k][j][P] tmp1[k][j][Q]
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
      // E_K = w sum_P Tr[M_P^T A G_P B^T]; (A G B^T)[a][b][P] = A*tmp1
      // = w * dot(M, AGBt). AGBt was folded into Gam; recompute cheaply:
      double e=0.0;
      for(int a2=0;a2<nao;++a2){
        // row a of AGBt = sum_c A[a][c] tmp1[c][:,:]
        // accumulate e += sum_bP M[a][b][P] * AGBt[a][b][P]
        // do it via temporary row buffer
        ;
      }
      (void)e;  // E_K diagnostic omitted (Gam/gam carry the physics)
    }
  }
  double t3=now_s();

  // symmetrize Gamma in (a,b) and gamma in (P,Q): the pass-2 contraction
  // visits unordered pairs with weight 2 and therefore requires the
  // symmetric parts (exact, since d(ab|P) and (P|Q)' are pair-symmetric).
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
    fprintf(stderr,"[routec-grad] mrsf timings: pass1 %.3f s | chol+J %.3f s "
            "| K-terms %.3f s | pass2 %.3f s | total %.3f s | E2e(J) %.12f\n",
            t1-t0,t2-t1,t3-t2,t4-t3,t4-t0,E2e);
  return 0;
}

// map a full-square OQP-frame matrix (Fortran column-major slice) to the
// pyscf frame, C row-major: A_py[q_i][q_j] = s_i s_j A_oqp(i,j)
static void map_oqp_full(const double* a_f,double* out){
  const int nao=g.nbf;
  for(int i=0;i<nao;++i)
    for(int j=0;j<nao;++j)
      out[(size_t)g.qmap[i]*nao+g.qmap[j]] =
          g.smap[i]*g.smap[j]*a_f[(size_t)j*nao+i];
}

} // namespace oqpgrad

extern "C" {

// the Fortran seam (all args by reference, bind(C) Fortran-side)
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
    fprintf(stderr,"[routec-grad] dE2e/dx done in %.3f s (hfscale=%.3f)\n",
            now_s()-tt0,*hfscale);
  }
}

// G-5 MRSF seam (tdhf_mrsf_gradient.F90 -> routec_bridge):
//   d, p   : (nbf,nbf,2) RAW alpha/beta SCF and relaxed difference
//            densities (full square, Fortran order, BEFORE the native
//            init() sum/diff transform)
//   spc    : (7,nbf,nbf) MRSF response densities (slot 7 = ball)
//   spcscale: (3) spc_coco, spc_ovov, spc_coov
//   mrst   : 1 singlet / 3 triplet (sign of the spin-pair terms)
// Term list = the G-5 phase A(ii) verified decomposition (q_i =
// sgnk*spcscale_i, c = coulscale, xc = hfscale, xc2 = hfscale2):
//   J: (D1,D1,c/2) (P1,D1,c/2) (D1,P1,c/2)
//      (sym bco1, sym bo2v, q3/2) (sym bco2, sym bo1v, q3/2)
//   K: s=1,2: (Ds,Ds,-xc/4) (Ps,Ds,-xc/4) (Ds,Ps,-xc/4)
//      (ball,ball,-xc2/2) (ball^T,ball^T,-xc2/2)
//      (co12,co12^T,q1) (o21v,o21v^T,q2)
//      (bco1,bo2v,-2 q3) (bco2,bo1v,-2 q3)
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
    fprintf(stderr,"[routec-grad] mrsf dim mismatch: seam nbf=%d natm=%d vs "
            "input %d/%d\n",*nbf,*natm,g.nbf,g.natm);
    return;
  }
  if(*mrst!=1&&*mrst!=3){
    fprintf(stderr,"[routec-grad] mrsf: unsupported mrst=%d\n",*mrst);
    return;
  }
  double tt0=now_s();
  const int nao=g.nbf;
  const size_t nn=(size_t)nao*nao;

  // map all inputs to the pyscf frame
  vector<vector<double>> dn;   // 0:Da 1:Db 2:Pa 3:Pb 4..10: spc slots 1..7
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
  // combinations: D1, D2, P1, P2, transposes, sym(bco/bo2v)
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
  // spc slot indices in dn: bo2v=4 bo1v=5 bco1=6 bco2=7 o21v=8 co12=9 ball=10
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
  // merged pairs (exact because Gamma/gamma are symmetrized before the
  // pass-2 contraction): (P,D)+(D,P) -> one term with doubled weight;
  // (ball,ball)+(ball^T,ball^T) -> (ball,ball) with doubled weight.
  (void)iballT;
  Kt.push_back({Dp(iD1),Dp(iD1),-xc/4});
  Kt.push_back({Dp(iP1),Dp(iD1),-xc/2});
  Kt.push_back({Dp(iD2),Dp(iD2),-xc/4});
  Kt.push_back({Dp(iP2),Dp(iD2),-xc/2});
  Kt.push_back({Dp(iball),Dp(iball),-xc2});
  if(q1!=0.0) Kt.push_back({Dp(ico12),Dp(ico12T),q1});
  if(q2!=0.0) Kt.push_back({Dp(io21v),Dp(io21vT),q2});
  if(q3!=0.0){
    Kt.push_back({Dp(6),Dp(4),-2*q3});   // (bco1, bo2v)
    Kt.push_back({Dp(7),Dp(5),-2*q3});   // (bco2, bo1v)
  }

  int rc=grad2e_mrsf(dn,Jt.data(),(int)Jt.size(),Kt.data(),(int)Kt.size(),
                     xyz,de);
  if(rc==0){
    *info=0;
    fprintf(stderr,"[routec-grad] MRSF dE2e/dx done in %.3f s "
            "(xc=%.3f xc2=%.3f mrst=%d, %d J + %d K terms)\n",
            now_s()-tt0,xc,xc2,*mrst,(int)Jt.size(),(int)Kt.size());
  }
}

}  // extern "C"
