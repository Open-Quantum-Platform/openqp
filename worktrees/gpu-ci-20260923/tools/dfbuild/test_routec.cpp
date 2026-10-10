// Phase 2 / gate G2 of PLAN_ROUTE_C.md: Route-C factorized assembly on CPU.
//
// Self-contained except (i) the exported algebra tables (routec_tables.bin,
// from Spherical_ERI/tools/export_routec_tables.py — the G0/G1-verified
// conventions) and (ii) rotc_order from libintRot for the per-lambda channel
// rotations (convention pinned at startup against Python reference matrices
// embedded in the table file).
//
// Validation: anchored MD reference (4 embedded ground-truth integrals from
// the Python eri_primitive), then Route-C vs MD per class.
// Benchmark: from-scratch and pair-cached q/s + exact MAC counts per stage.
//
//   build: g++ -O3 -std=c++17 -I include test/test_routec.cpp -o test_routec
//   run:   ./test_routec tools/routec_tables.bin [l=2] [K=1] [Nq=2000]
#include "sph_eri/rotc_cuda.cuh"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <array>
#include <chrono>
#include <random>
#include <map>

using std::vector;
static double now_s(){using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();}

// ---------------------------------------------------------------- op counters
struct Ops { uint64_t pair=0, rot=0, hbuild=0, tbuild=0, gemm=0;
  void zero(){pair=rot=hbuild=tbuild=gemm=0;}
  uint64_t total() const {return pair+rot+hbuild+tbuild+gemm;} };
static Ops g_ops;

// per-stage wall-clock profile of the quartet-only path (ROUTEC_PROF=1)
struct Prof { double frame=0, rotb=0, rotk=0, rt=0, h0=0, grp=0; bool on=false;
  void zero(){frame=rotb=rotk=rt=h0=grp=0;} } static g_prof;

// ---------------------------------------------------------------- monomials
static vector<std::array<int,3>> monos_of(int k){
  vector<std::array<int,3>> v;
  for(int lx=k; lx>=0; --lx) for(int ly=k-lx; ly>=0; --ly) v.push_back({lx,ly,k-lx-ly});
  return v;
}
static int nmono(int k){ return (k+1)*(k+2)/2; }

// ---------------------------------------------------------------- Boys
static void boys(int nmax, double T, double* F){
  if(T < 1e-13){ for(int n=0;n<=nmax;++n) F[n]=1.0/(2*n+1); return; }
  if(T > 40.0){
    F[0] = 0.5*std::sqrt(M_PI/T);
    for(int n=1;n<=nmax;++n) F[n] = F[n-1]*(2*n-1)/(2.0*T);
    return;
  }
  // top by series, downward recursion
  double e = std::exp(-T);
  double term = 1.0/(2*nmax+1), s = term;
  for(int i=1;i<500;++i){ term *= 2.0*T/(2*nmax+2*i+1); s += term;
    if(term < 1e-17*s) break; }
  F[nmax] = s*e;
  for(int n=nmax;n>0;--n) F[n-1] = (2.0*T*F[n]+e)/(2*n-1);
}

// ------------------------------------------------- Hermite E (one axis, DP)
// E^{ij}_t per the Python cartesian_eri.hermite_E recursion (exact match).
static void e_axis(int i,int j,double Qx,double a,double b,double* out /*len i+j+1*/){
  const double p=a+b, mu=a*b/p, inv2p=0.5/p;
  // DP over (ii,jj): table of arrays length (ii+jj+1)
  static thread_local vector<double> buf; // [(ii)*(j+1)+(jj)] -> offset table
  int I=i+1, J=j+1, W=i+j+1;
  buf.assign((size_t)I*J*W, 0.0);
  auto at=[&](int ii,int jj)->double*{ return &buf[((size_t)ii*J+jj)*W]; };
  at(0,0)[0]=std::exp(-mu*Qx*Qx);
  for(int ii=1;ii<I;++ii){ double* cur=at(ii,0); double* lo=at(ii-1,0);
    for(int t=0;t<=ii;++t){ double v=0.0;
      if(t-1>=0) v += inv2p*lo[t-1];
      v += (-(mu*Qx)/a)*lo[t];
      if(t+1<=ii-1) v += (t+1)*lo[t+1];
      cur[t]=v; } }
  for(int jj=1;jj<J;++jj) for(int ii=0;ii<I;++ii){
    double* cur=at(ii,jj); double* lo=at(ii,jj-1);
    for(int t=0;t<=ii+jj;++t){ double v=0.0;
      if(t-1>=0) v += inv2p*lo[t-1];
      v += ((mu*Qx)/b)*lo[t];
      if(t+1<=ii+jj-1) v += (t+1)*lo[t+1];
      cur[t]=v; } }
  std::memcpy(out, at(i,j), sizeof(double)*W);
}

// ------------------------------------------------- erf-RSH attenuation seam
// Operator erf(w*r12)/r12 (long-range part of CAM range separation): EXACT
// substitution F_n(eta R^2) -> s^{n+1/2} F_n(s*eta R^2), s = w^2/(eta+w^2)
// [standard erf-attenuated Boys rescale; rho == our eta]. Equivalently the
// SAME R-tensor with eta -> s*eta and one overall sqrt(s) on the F_n seeds:
// seed_n = (-2 eta)^n F_n = sqrt(s) * (-2 eta_eff)^n F_n(eta_eff R^2),
// eta_eff = s*eta. The spatial recursion and ALL angular/rotation/coupling
// structure are untouched -- the attenuation enters ONLY here.
// Enabled at runtime via env ROUTEC_OMEGA (Bohr^-1); unset/0 = full Coulomb.
static const double g_rsh_omega2 = [](){
  const char* e = std::getenv("ROUTEC_OMEGA");
  if(!e) return -1.0;
  const double w = std::atof(e);
  return w > 0.0 ? w*w : -1.0; }();

// ---------------------------------------------------------------- R-tensor
// S(T) = d_R^T F0(eta R^2), dense over |T| <= L. layout idx = (tx*(L+1)+ty)*(L+1)+tz
struct RT {
  int L; vector<double> v; // only |T|<=L valid
  double at(int tx,int ty,int tz) const { return v[((size_t)tx*(L+1)+ty)*(L+1)+tz]; }
};
static void r_tensor(double eta, const double* R, int L, RT& out, bool count=false){
  out.L=L; int L1=L+1;
  double seed_scale=1.0;             // erf-RSH (see above); 1.0 when off
  if(g_rsh_omega2>0.0){
    const double s_=g_rsh_omega2/(eta+g_rsh_omega2);
    eta*=s_; seed_scale=std::sqrt(s_);
  }
  double F[16]; boys(L, eta*(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]), F);
  // levels n: flat workspace reused across calls (L<=12 -> 13^4 doubles max)
  static thread_local vector<double> lev;
  size_t cube=(size_t)L1*L1*L1;
  lev.assign((size_t)L1*cube, 0.0);
  double m2e=seed_scale;
  for(int n=0;n<L1;++n){ lev[(size_t)n*cube] = m2e*F[n]; m2e *= -2.0*eta; }
  const bool aligned = (R[0]==0.0 && R[1]==0.0);
  auto idx=[&](int a,int b,int c){return ((size_t)a*L1+b)*L1+c;};
  for(int tot=1; tot<=L; ++tot)
    for(int tx=tot; tx>=0; --tx) for(int ty=tot-tx; ty>=0; --ty){ int tz=tot-tx-ty;
      // aligned frame: entries with odd tx or ty vanish identically
      if(aligned && ((tx|ty)&1)) continue;
      for(int n=0;n<=L-tot;++n){ double val=0.0;
        double* up=&lev[(size_t)(n+1)*cube];
        if(tx>0){ if(!aligned) val = R[0]*up[idx(tx-1,ty,tz)];
                  if(tx>1) val += (tx-1)*up[idx(tx-2,ty,tz)]; }
        else if(ty>0){ if(!aligned) val = R[1]*up[idx(tx,ty-1,tz)];
                  if(ty>1) val += (ty-1)*up[idx(tx,ty-2,tz)]; }
        else { val = R[2]*up[idx(tx,ty,tz-1)];
                  if(tz>1) val += (tz-1)*up[idx(tx,ty,tz-2)]; }
        lev[(size_t)n*cube+idx(tx,ty,tz)] = val;
        if(count) g_ops.hbuild += 2; } }
  out.v.assign(lev.begin(), lev.begin()+cube);
}

// ------------------------------------------------------------- MD reference
static double md_reference(const int* lA,const double* A,double a,
                           const int* lB,const double* B,double b,
                           const int* lC,const double* C,double c,
                           const int* lD,const double* D,double d){
  double p=a+b,q=c+d;
  double P[3],Q[3],PQ[3];
  for(int x=0;x<3;++x){ P[x]=(a*A[x]+b*B[x])/p; Q[x]=(c*C[x]+d*D[x])/q; PQ[x]=P[x]-Q[x]; }
  double eta=p*q/(p+q);
  int Ltot=lA[0]+lA[1]+lA[2]+lB[0]+lB[1]+lB[2]+lC[0]+lC[1]+lC[2]+lD[0]+lD[1]+lD[2];
  RT S; r_tensor(eta,PQ,Ltot,S);
  vector<double> Ebx(lA[0]+lB[0]+1),Eby(lA[1]+lB[1]+1),Ebz(lA[2]+lB[2]+1);
  vector<double> Ekx(lC[0]+lD[0]+1),Eky(lC[1]+lD[1]+1),Ekz(lC[2]+lD[2]+1);
  e_axis(lA[0],lB[0],A[0]-B[0],a,b,Ebx.data());
  e_axis(lA[1],lB[1],A[1]-B[1],a,b,Eby.data());
  e_axis(lA[2],lB[2],A[2]-B[2],a,b,Ebz.data());
  e_axis(lC[0],lD[0],C[0]-D[0],c,d,Ekx.data());
  e_axis(lC[1],lD[1],C[1]-D[1],c,d,Eky.data());
  e_axis(lC[2],lD[2],C[2]-D[2],c,d,Ekz.data());
  double pref=2.0*std::pow(M_PI,2.5)/(p*q*std::sqrt(p+q));
  double tot=0.0;
  for(size_t t=0;t<Ebx.size();++t) for(size_t u=0;u<Eby.size();++u) for(size_t v=0;v<Ebz.size();++v){
    double eb=Ebx[t]*Eby[u]*Ebz[v]; if(eb==0.0) continue;
    for(size_t tp=0;tp<Ekx.size();++tp) for(size_t up=0;up<Eky.size();++up) for(size_t vp=0;vp<Ekz.size();++vp){
      double ek=Ekx[tp]*Eky[up]*Ekz[vp]; if(ek==0.0) continue;
      double sgn=((tp+up+vp)&1)?-1.0:1.0;
      tot += eb*ek*sgn*S.at(t+tp,u+up,v+vp); } }
  return pref*tot;
}

// ---------------------------------------------------------------- tables
struct DegTab {
  int nmono=0, nch=0, nch0=0;
  vector<std::array<int,3>> ch;   // (n,lam,mu)
  vector<double> C;               // nch x nmono row-major
  vector<int> ch0;                // indices of mu==0 channels
  vector<double> B0;              // nmono x nch0 row-major
};
struct Coo { int i1,i2,j0; double v; };
struct RTC {
  int kmax=0, kpair=0;
  vector<DegTab> deg;
  vector<vector<Coo>> g0;         // [(k1*(kpair+1))+k2]
  // rotation reference: per test U: U[9] + Rl per lam<=6
  struct RotRef { double U[9]; vector<vector<double>> Rl; };
  vector<RotRef> rotref;
  // c2s per l: ncart x nsph row-major
  vector<vector<double>> c2s; vector<int> c2s_nc, c2s_ns;
};
static bool load_tables(const char* path, RTC& T){
  FILE* f=fopen(path,"rb"); if(!f) return false;
  char magic[4]; int km,kp;
  if(fread(magic,1,4,f)!=4||memcmp(magic,"RTC2",4)) {fclose(f);return false;}
  fread(&km,4,1,f); fread(&kp,4,1,f); T.kmax=km; T.kpair=kp;
  T.deg.resize(km+1);
  for(int k=0;k<=km;++k){ DegTab& D=T.deg[k];
    fread(&D.nmono,4,1,f); fread(&D.nch,4,1,f); fread(&D.nch0,4,1,f);
    D.ch.resize(D.nch);
    for(int j=0;j<D.nch;++j){int b[3];fread(b,4,3,f);D.ch[j]={b[0],b[1],b[2]};}
    D.C.resize((size_t)D.nch*D.nmono); fread(D.C.data(),8,D.C.size(),f);
    D.ch0.resize(D.nch0); fread(D.ch0.data(),4,D.nch0,f);
    D.B0.resize((size_t)D.nmono*D.nch0); fread(D.B0.data(),8,D.B0.size(),f);
  }
  int npairs; fread(&npairs,4,1,f);
  T.g0.assign((size_t)(kp+1)*(kp+1), {});
  for(int p=0;p<npairs;++p){ int k1,k2,nnz; fread(&k1,4,1,f);fread(&k2,4,1,f);fread(&nnz,4,1,f);
    auto& v=T.g0[(size_t)k1*(kp+1)+k2]; v.resize(nnz);
    for(int i=0;i<nnz;++i){int b[3];fread(b,4,3,f);double val;fread(&val,8,1,f);
      v[i]={b[0],b[1],b[2],val}; } }
  int nrot; fread(&nrot,4,1,f); T.rotref.resize(nrot);
  for(int r=0;r<nrot;++r){ fread(T.rotref[r].U,8,9,f);
    T.rotref[r].Rl.resize(7);
    for(int lam=0;lam<=6;++lam){ T.rotref[r].Rl[lam].resize((size_t)(2*lam+1)*(2*lam+1));
      fread(T.rotref[r].Rl[lam].data(),8,T.rotref[r].Rl[lam].size(),f);} }
  int nl; fread(&nl,4,1,f);
  T.c2s.resize(nl); T.c2s_nc.resize(nl); T.c2s_ns.resize(nl);
  for(int l=0;l<nl;++l){ int nc,ns; fread(&nc,4,1,f); fread(&ns,4,1,f);
    T.c2s_nc[l]=nc; T.c2s_ns[l]=ns;
    T.c2s[l].resize((size_t)nc*ns); fread(T.c2s[l].data(),8,T.c2s[l].size(),f); }
  fclose(f); return true;
}

// ------------------------------------------------- rotc convention pinning
// Python G1 pinned: F_tilde = F @ Rl_py(U)  (columns of F are channels).
// Find how to get Rl_py(U) from sph_eri::rotc_order: input U or U^T, output
// direct or transposed.  Returns mode 0..3; -1 if none matches.
static int pin_rotc_mode(const RTC& T){
  for(int mode=0;mode<4;++mode){
    double worst=0.0;
    for(const auto& rr : T.rotref){
      double R3[9];
      if(mode&1){ for(int i=0;i<3;++i)for(int j=0;j<3;++j) R3[i*3+j]=rr.U[j*3+i]; }
      else      { memcpy(R3,rr.U,sizeof R3); }
      for(int lam=1;lam<=6;++lam){
        int n=2*lam+1; vector<double> Rl((size_t)n*n);
        sph_eri::rotc_order<double>(R3,lam,Rl.data());
        for(int i=0;i<n;++i)for(int j=0;j<n;++j){
          double ref=rr.Rl[lam][(size_t)i*n+j];
          double got=(mode&2)?Rl[(size_t)j*n+i]:Rl[(size_t)i*n+j];
          worst=std::max(worst,std::fabs(got-ref)); } } }
    if(worst<1e-10){ printf("  rotc convention pinned: mode=%d (input %s, output %s), err %.2e\n",
        mode,(mode&1)?"U^T":"U",(mode&2)?"transposed":"direct",worst); return mode; }
  }
  return -1;
}
static void get_Rl(const double* U,int lam,int mode,double* Rl){
  double R3[9];
  if(mode&1){ for(int i=0;i<3;++i)for(int j=0;j<3;++j) R3[i*3+j]=U[j*3+i]; }
  else memcpy(R3,U,9*sizeof(double));
  int n=2*lam+1;
  if(mode&2){ vector<double> tmp((size_t)n*n);
    sph_eri::rotc_order<double>(R3,lam,tmp.data());
    for(int i=0;i<n;++i)for(int j=0;j<n;++j) Rl[(size_t)i*n+j]=tmp[(size_t)j*n+i];
  } else sph_eri::rotc_order<double>(R3,lam,Rl);
}

// ------------------------------------------------- pair channel tensors
struct PairF {
  int la, lb, K;                // K = la+lb
  int ncomp=0, nchtot=0;
  vector<int> off;              // per degree k: column offset
  vector<double> F;             // ncomp x nchtot row-major
};
static void build_pair(const RTC& T,int la,int lb,const double* A,const double* B,
                       double a,double b,PairF& P,bool count,bool sph=true){
  P.la=la; P.lb=lb; P.K=la+lb;
  auto ca=monos_of(la), cb=monos_of(lb);
  int ncart=(int)(ca.size()*cb.size());
  P.ncomp=ncart;
  P.off.assign(P.K+1,0);
  int tot=0; for(int k=0;k<=P.K;++k){ P.off[k]=tot; tot+=T.deg[k].nch; }
  P.nchtot=tot;
  P.F.assign((size_t)ncart*tot,0.0);
  vector<double> Ex(2*P.K+2),Ey(2*P.K+2),Ez(2*P.K+2), e;
  int ic=0;
  for(auto& ia: ca) for(auto& ib: cb){
    e_axis(ia[0],ib[0],A[0]-B[0],a,b,Ex.data());
    e_axis(ia[1],ib[1],A[1]-B[1],a,b,Ey.data());
    e_axis(ia[2],ib[2],A[2]-B[2],a,b,Ez.data());
    int kmaxc=ia[0]+ib[0]+ia[1]+ib[1]+ia[2]+ib[2];
    for(int k=0;k<=std::min(P.K,kmaxc);++k){
      const DegTab& D=T.deg[k];
      auto mlist=monos_of(k);
      e.assign(D.nmono,0.0); bool nz=false;
      for(int m=0;m<D.nmono;++m){ auto& mo=mlist[m];
        if(mo[0]<=ia[0]+ib[0]&&mo[1]<=ia[1]+ib[1]&&mo[2]<=ia[2]+ib[2]){
          double v=Ex[mo[0]]*Ey[mo[1]]*Ez[mo[2]];
          if(v!=0.0){e[m]=v;nz=true;} } }
      if(!nz) continue;
      double* Fr=&P.F[(size_t)ic*P.nchtot+P.off[k]];
      for(int chi=0;chi<D.nch;++chi){ double s2=0.0;
        const double* Crow=&D.C[(size_t)chi*D.nmono];
        for(int m=0;m<D.nmono;++m) s2+=Crow[m]*e[m];
        Fr[chi]=s2; }
      if(count) g_ops.pair += (uint64_t)D.nch*D.nmono;
    }
    ++ic;
  }
  if(!sph) return;
  // spherical fold of the component axis (pair-level, cacheable):
  // F_sph[(ma,mb),ch] = sum_{ia,ib} c2sA[ia,ma] c2sB[ib,mb] F[(ia,ib),ch]
  int nca=(int)ca.size(), ncb=(int)cb.size();
  int nsa=T.c2s_ns[la], nsb=T.c2s_ns[lb];
  const double* CA=T.c2s[la].data();
  const double* CB=T.c2s[lb].data();
  static thread_local vector<double> half, outF;
  half.assign((size_t)nsa*ncb*tot,0.0);
  for(int ma=0;ma<nsa;++ma)
    for(int iaq=0;iaq<nca;++iaq){ double cma=CA[(size_t)iaq*nsa+ma];
      if(cma==0.0) continue;
      for(int ibq=0;ibq<ncb;++ibq){
        const double* src=&P.F[((size_t)iaq*ncb+ibq)*tot];
        double* dst=&half[((size_t)ma*ncb+ibq)*tot];
        for(int ch=0;ch<tot;++ch) dst[ch]+=cma*src[ch];
        if(count) g_ops.pair += (uint64_t)tot; } }
  outF.assign((size_t)nsa*nsb*tot,0.0);
  for(int mb=0;mb<nsb;++mb)
    for(int ibq=0;ibq<ncb;++ibq){ double cmb=CB[(size_t)ibq*nsb+mb];
      if(cmb==0.0) continue;
      for(int ma=0;ma<nsa;++ma){
        const double* src=&half[((size_t)ma*ncb+ibq)*tot];
        double* dst=&outF[((size_t)ma*nsb+mb)*tot];
        for(int ch=0;ch<tot;++ch) dst[ch]+=cmb*src[ch];
        if(count) g_ops.pair += (uint64_t)tot; } }
  P.ncomp=nsa*nsb;
  P.F.assign(outF.begin(), outF.begin()+(size_t)P.ncomp*tot);
}
// rotation matrices for one quartet frame, built once and shared by both pairs
struct QuartetRot { int lmax=0; double Rl[9][17*17]; };
static void build_qrot(const double* U,int lmax,int mode,QuartetRot& QR){
  QR.lmax=lmax;
  for(int lam=1;lam<=lmax;++lam) get_Rl(U,lam,mode,QR.Rl[lam]);
}
// rotate channel (mu) blocks: F_tilde = F @ Rl per (k, lambda-block)
static void rotate_pair(const RTC& T,PairF& P,const QuartetRot& QR,bool count){
  double tmp[17];
  for(int k=0;k<=P.K;++k){ const DegTab& D=T.deg[k];
    int pos=0;
    while(pos<D.nch){ int lam=D.ch[pos][1]; int n=2*lam+1;
      if(lam>0){
        const double* Rl=QR.Rl[lam];
        for(int icp=0;icp<P.ncomp;++icp){
          double* Fb=&P.F[(size_t)icp*P.nchtot+P.off[k]+pos];
          for(int jj=0;jj<n;++jj){ double s2=0.0;
            for(int ii=0;ii<n;++ii) s2+=Fb[ii]*Rl[(size_t)ii*n+jj];
            tmp[jj]=s2; }
          for(int jj=0;jj<n;++jj) Fb[jj]=tmp[jj];
        }
        if(count) g_ops.rot += (uint64_t)P.ncomp*n*n;
      }
      pos+=n; }
  }
}

// ------------------------------------------------- per-class coupling plan
// mu-sorted contiguous layout: channel columns permuted so each mu-group is a
// contiguous range; gamma0 entries pre-resolved to dense (row,col) targets with
// the ket sign folded in.  Built once per (K1,K2), cached.
struct ClassPlan {
  int K1=0,K2=0,nchb=0,nchk=0;
  vector<int> permB, permK;             // permuted pos -> original column
  vector<int> gposB, gposK;             // original column -> permuted pos
  struct Grp{int ob,nb,ok,nk;};
  vector<Grp> grps;
  struct Ent{int row,col,j0; double v;};
  vector<vector<Ent>> ents;             // per grp
  vector<int> h0off; int h0tot=0;
};
static void build_side(const RTC& T,int K,vector<int>& off,vector<int>& mu_of,int& nch){
  off.assign(K+1,0); int tot=0;
  for(int k=0;k<=K;++k){ off[k]=tot; tot+=T.deg[k].nch; }
  nch=tot; mu_of.assign(tot,0);
  for(int k=0;k<=K;++k) for(int j=0;j<T.deg[k].nch;++j)
    mu_of[off[k]+j]=T.deg[k].ch[j][2];
}
static const ClassPlan& get_plan(const RTC& T,int K1,int K2){
  // thread_local: build_df now runs the pair build + Schwarz under OpenMP; a
  // shared map with unguarded concurrent insert races (intermittent SIGSEGV).
  // Plans are small (a few classes), so per-thread copies are cheap.
  static thread_local std::map<int,ClassPlan> cache;
  int key=K1*64+K2;
  auto it=cache.find(key); if(it!=cache.end()) return it->second;
  ClassPlan P; P.K1=K1; P.K2=K2;
  vector<int> off1,off2,mu1,mu2;
  build_side(T,K1,off1,mu1,P.nchb);
  build_side(T,K2,off2,mu2,P.nchk);
  // permutations: stable order by (mu, original pos)
  auto mk_perm=[&](const vector<int>& mu,int n,vector<int>& perm,
                   vector<int>& gpos,vector<int>& gof,vector<int>& gn){
    perm.clear(); gpos.assign(n,-1); gof.assign(25,0); gn.assign(25,0);
    for(int m=-12;m<=12;++m){ gof[m+12]=(int)perm.size();
      for(int p=0;p<n;++p) if(mu[p]==m){ gpos[p]=(int)perm.size(); perm.push_back(p); }
      gn[m+12]=(int)perm.size()-gof[m+12]; } };
  vector<int> gofB,gnB,gofK,gnK;
  vector<int>& gposB=P.gposB; vector<int>& gposK=P.gposK;
  mk_perm(mu1,P.nchb,P.permB,gposB,gofB,gnB);
  mk_perm(mu2,P.nchk,P.permK,gposK,gofK,gnK);
  // groups present on both sides
  vector<int> grp_of_mu(25,-1);
  for(int mi=0;mi<25;++mi) if(gnB[mi]>0&&gnK[mi]>0){
    grp_of_mu[mi]=(int)P.grps.size();
    P.grps.push_back({gofB[mi],gnB[mi],gofK[mi],gnK[mi]}); }
  P.ents.resize(P.grps.size());
  // H0 flat offsets
  P.h0off.assign(K1+K2+1,0); int h=0;
  for(int k=0;k<=K1+K2;++k){ P.h0off[k]=h; h+=T.deg[k].nch0; }
  P.h0tot=h;
  // resolve coupling entries
  for(int k1=0;k1<=K1;++k1) for(int k2=0;k2<=K2;++k2){
    const auto& coo=T.g0[(size_t)k1*(T.kpair+1)+k2];
    double sgn=(k2&1)?-1.0:1.0;
    for(const auto& c0 : coo){
      int gb=off1[k1]+c0.i1, gk=off2[k2]+c0.i2;
      int m=mu1[gb];
      if(mu2[gk]!=m) continue;        // strict rule (exporter enforces)
      int gi=grp_of_mu[m+12]; if(gi<0) continue;
      P.ents[gi].push_back({gposB[gb]-P.grps[gi].ob, gposK[gk]-P.grps[gi].ok,
                            P.h0off[k1+k2]+c0.j0, sgn*c0.v}); } }
  return cache.emplace(key,std::move(P)).first->second;
}
static void permute_cols(const PairF& P,const vector<int>& perm,vector<double>& out){
  out.assign((size_t)P.ncomp*P.nchtot,0.0);
  for(int c=0;c<P.ncomp;++c){
    const double* src=&P.F[(size_t)c*P.nchtot];
    double* dst=&out[(size_t)c*P.nchtot];
    for(int p=0;p<P.nchtot;++p) dst[p]=src[perm[p]]; }
}
// fused rotate+permute: read the (cached, unrotated) pair tensor, write the
// quartet-frame mu-sorted tensor.  No pair copy, one pass.
static void rotate_permute(const RTC& T,const PairF& P,const QuartetRot& QR,
                           const vector<int>& gpos,vector<double>& out,bool count){
  out.assign((size_t)P.ncomp*P.nchtot,0.0);
  double tmp[17];
  for(int icp=0;icp<P.ncomp;++icp){
    const double* src=&P.F[(size_t)icp*P.nchtot];
    double* dst=&out[(size_t)icp*P.nchtot];
    for(int k=0;k<=P.K;++k){ const DegTab& D=T.deg[k];
      int pos=0;
      while(pos<D.nch){ int lam=D.ch[pos][1]; int n=2*lam+1;
        const double* Fb=src+P.off[k]+pos;
        if(lam==0){ dst[gpos[P.off[k]+pos]]=Fb[0]; }
        else {
          const double* Rl=QR.Rl[lam];
          for(int jj=0;jj<n;++jj){ double s2=0.0;
            for(int ii=0;ii<n;++ii) s2+=Fb[ii]*Rl[(size_t)ii*n+jj];
            tmp[jj]=s2; }
          for(int jj=0;jj<n;++jj) dst[gpos[P.off[k]+pos+jj]]=tmp[jj];
        }
        pos+=n; } }
    if(count){ for(int k=0;k<=P.K;++k){ const DegTab& D=T.deg[k];
      int pos=0; while(pos<D.nch){ int lam=D.ch[pos][1];
        if(lam>0) g_ops.rot += (uint64_t)(2*lam+1)*(2*lam+1);
        pos+=2*lam+1; } } }
  }
}

// ------------------------------------------------- planned coupling
static void couple_planned(const RTC& T,const ClassPlan& CP,
    const double* Fbp,int ncompb,const double* Fkp,int ncompk,
    double eta,double Rn,double pref,double* out,bool count){
  int L=CP.K1+CP.K2;
  double Rz[3]={0.0,0.0,Rn};
  double tp0=g_prof.on?now_s():0.0;
  static thread_local RT S;
  r_tensor(eta,Rz,L,S,count);
  double tp1=g_prof.on?now_s():0.0;
  static thread_local vector<double> H0;
  H0.assign(CP.h0tot,0.0);
  for(int k=0;k<=L;++k){ const DegTab& D=T.deg[k];
    auto mlist=monos_of(k);
    for(int j0=0;j0<D.nch0;++j0){ double s2=0.0;
      for(int m=0;m<D.nmono;++m){ auto& mo=mlist[m];
        if((mo[0]|mo[1])&1) continue;
        s2+=D.B0[(size_t)m*D.nch0+j0]*S.at(mo[0],mo[1],mo[2]); }
      H0[CP.h0off[k]+j0]=s2; }
    if(count) g_ops.hbuild += (uint64_t)D.nch0*((D.nmono+3)/4+1);
  }
  double tp2=g_prof.on?now_s():0.0;
  if(g_prof.on){ g_prof.rt+=tp1-tp0; g_prof.h0+=tp2-tp1; }
  static thread_local vector<double> Tmu,W;
  for(size_t gi=0;gi<CP.grps.size();++gi){
    const auto& G=CP.grps[gi];
    const auto& E=CP.ents[gi];
    if(E.empty()) continue;
    Tmu.assign((size_t)G.nb*G.nk,0.0);
    for(const auto& e : E) Tmu[(size_t)e.row*G.nk+e.col]+=e.v*H0[e.j0];
    if(count) g_ops.tbuild += E.size();
    W.assign((size_t)ncompb*G.nk,0.0);
    for(int ic=0;ic<ncompb;++ic){
      const double* Fr=Fbp+(size_t)ic*CP.nchb+G.ob;
      double* Wr=&W[(size_t)ic*G.nk];
      for(int i=0;i<G.nb;++i){ double f=Fr[i];
        if(f==0.0) continue;
        const double* Tr=&Tmu[(size_t)i*G.nk];
        for(int j=0;j<G.nk;++j) Wr[j]+=f*Tr[j]; } }
    if(count) g_ops.gemm += (uint64_t)ncompb*G.nb*G.nk;
    for(int ic=0;ic<ncompb;++ic){
      const double* Wr=&W[(size_t)ic*G.nk];
      double* Or=out+(size_t)ic*ncompk;
      for(int jc=0;jc<ncompk;++jc){
        const double* Fr=Fkp+(size_t)jc*CP.nchk+G.ok;
        double s2=0.0;
        for(int j=0;j<G.nk;++j) s2+=Wr[j]*Fr[j];
        Or[jc]+=s2; } }
    if(count) g_ops.gemm += (uint64_t)ncompb*ncompk*G.nk;
  }
  for(int i=0;i<ncompb*ncompk;++i) out[i]*=pref;
  if(g_prof.on) g_prof.grp+=now_s()-tp2;
}

// ------------------------------------------------- Route-C quartet block
// out: ncomp_b x ncomp_k (caller-sized). Pair tensors must be PRE-ROTATED
// copies for this quartet (aligned path).  [legacy unplanned path kept for
// reference below]
static void routec_couple(const RTC& T,const PairF& Pb,const PairF& Pk,
                          double eta,double Rn,double pref,double* out,bool count){
  int K1=Pb.K, K2=Pk.K, L=K1+K2;
  // aligned H0 per total degree
  double Rz[3]={0.0,0.0,Rn};
  static thread_local RT S;
  r_tensor(eta,Rz,L,S,count);
  static thread_local vector<vector<double>> H0;
  if((int)H0.size()<L+1) H0.resize(L+1);
  for(int k=0;k<=L;++k){ const DegTab& D=T.deg[k];
    H0[k].assign(D.nch0,0.0);
    auto mlist=monos_of(k);
    for(int j0=0;j0<D.nch0;++j0){ double s2=0.0;
      for(int m=0;m<D.nmono;++m){ auto& mo=mlist[m];
        if((mo[0]|mo[1])&1) continue;            // aligned S: odd x/y vanish
        s2+=D.B0[(size_t)m*D.nch0+j0]*S.at(mo[0],mo[1],mo[2]); }
      H0[k][j0]=s2; }
    if(count) g_ops.hbuild += (uint64_t)D.nch0*((D.nmono+3)/4+1);
  }
  // T blocks (full bra-channels x ket-channels), mu-paired sparse by construction
  static thread_local vector<double> Tm;
  Tm.assign((size_t)Pb.nchtot*Pk.nchtot,0.0);
  for(int k1=0;k1<=K1;++k1) for(int k2=0;k2<=K2;++k2){
    const auto& coo=T.g0[(size_t)k1*(T.kpair+1)+k2];
    double sgn=(k2&1)?-1.0:1.0;     // ket (-1)^{|u|}
    const auto& H=H0[k1+k2];
    for(const auto& c0 : coo)
      Tm[(size_t)(Pb.off[k1]+c0.i1)*Pk.nchtot + (Pk.off[k2]+c0.i2)] += sgn*c0.v*H[c0.j0];
    if(count) g_ops.tbuild += coo.size();
  }
  // mu-blocked contraction: group channels by mu (signed)
  // build index lists once per call (cheap; could be cached per class)
  auto groups=[&](const PairF& P){
    // mu in [-Lmax..Lmax] -> list of channel column indices
    vector<vector<int>> g(2*12+1);
    for(int k=0;k<=P.K;++k){ const DegTab& D=T.deg[k];
      for(int j=0;j<D.nch;++j) g[D.ch[j][2]+12].push_back(P.off[k]+j); }
    return g; };
  auto gb=groups(Pb), gk=groups(Pk);
  static thread_local vector<double> Tsub, W;
  for(int mi=0;mi<25;++mi){
    const auto& Ib=gb[mi]; const auto& Jk=gk[mi];
    if(Ib.empty()||Jk.empty()) continue;
    int nb=(int)Ib.size(), nk=(int)Jk.size();
    Tsub.assign((size_t)nb*nk,0.0);
    bool any=false;
    for(int i=0;i<nb;++i) for(int j=0;j<nk;++j){
      double v=Tm[(size_t)Ib[i]*Pk.nchtot+Jk[j]];
      Tsub[(size_t)i*nk+j]=v; if(v!=0.0) any=true; }
    if(!any) continue;
    // W = Fb[:,Ib] * Tsub   (ncompb x nk)
    W.assign((size_t)Pb.ncomp*nk,0.0);
    for(int icp=0;icp<Pb.ncomp;++icp){
      const double* Fr=&Pb.F[(size_t)icp*Pb.nchtot];
      double* Wr=&W[(size_t)icp*nk];
      for(int i=0;i<nb;++i){ double f=Fr[Ib[i]];
        if(f==0.0) continue;
        const double* Tr=&Tsub[(size_t)i*nk];
        for(int j=0;j<nk;++j) Wr[j]+=f*Tr[j]; } }
    if(count) g_ops.gemm += (uint64_t)Pb.ncomp*nb*nk;
    // out += W * Fk[:,Jk]^T
    for(int icp=0;icp<Pb.ncomp;++icp){
      const double* Wr=&W[(size_t)icp*nk];
      double* Or=&out[(size_t)icp*Pk.ncomp];
      for(int jcp=0;jcp<Pk.ncomp;++jcp){
        const double* Fr=&Pk.F[(size_t)jcp*Pk.nchtot];
        double s2=0.0;
        for(int j=0;j<nk;++j) s2+=Wr[j]*Fr[Jk[j]];
        Or[jcp]+=s2; } }
    if(count) g_ops.gemm += (uint64_t)Pb.ncomp*Pk.ncomp*nk;
  }
  for(int i=0;i<Pb.ncomp*Pk.ncomp;++i) out[i]*=pref;
}

// full from-scratch block for one primitive quartet
static void routec_block(const RTC& T,int mode,
                         int la,const double* A,double a,int lb,const double* B,double b,
                         int lc,const double* C,double c,int ld,const double* D,double d,
                         vector<double>& out,bool count){
  double p=a+b,q=c+d;
  double P[3],Q[3],R[3];
  for(int x=0;x<3;++x){P[x]=(a*A[x]+b*B[x])/p;Q[x]=(c*C[x]+d*D[x])/q;R[x]=P[x]-Q[x];}
  double eta=p*q/(p+q);
  double Rn=std::sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
  double pref=2.0*std::pow(M_PI,2.5)/(p*q*std::sqrt(p+q));
  PairF Pb,Pk;
  build_pair(T,la,lb,A,B,a,b,Pb,count);
  build_pair(T,lc,ld,C,D,c,d,Pk,count);
  // frame: z = Rhat; R->0: any orthonormal frame is exact (V continuous,
  // frame-dependence O(Rn)) — previously NaN, now identity frame
  double z[3];
  if(Rn>1e-12){ z[0]=R[0]/Rn; z[1]=R[1]/Rn; z[2]=R[2]/Rn; }
  else { z[0]=0.0; z[1]=0.0; z[2]=1.0; }
  double aa[3]; if(std::fabs(z[0])<0.9){aa[0]=1;aa[1]=0;aa[2]=0;} else {aa[0]=0;aa[1]=1;aa[2]=0;}
  double dot=aa[0]*z[0]+aa[1]*z[1]+aa[2]*z[2];
  double x[3]={aa[0]-dot*z[0],aa[1]-dot*z[1],aa[2]-dot*z[2]};
  double xn=std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  for(int i=0;i<3;++i)x[i]/=xn;
  double y[3]={z[1]*x[2]-z[2]*x[1], z[2]*x[0]-z[0]*x[2], z[0]*x[1]-z[1]*x[0]};
  double U[9]={x[0],y[0],z[0], x[1],y[1],z[1], x[2],y[2],z[2]}; // columns x,y,z
  static thread_local QuartetRot QR;
  build_qrot(U,std::max(Pb.K,Pk.K),mode,QR);
  const ClassPlan& CP=get_plan(T,Pb.K,Pk.K);
  static thread_local vector<double> Fbp,Fkp;
  rotate_permute(T,Pb,QR,CP.gposB,Fbp,count);
  rotate_permute(T,Pk,QR,CP.gposK,Fkp,count);
  out.assign((size_t)Pb.ncomp*Pk.ncomp,0.0);
  couple_planned(T,CP,Fbp.data(),Pb.ncomp,Fkp.data(),Pk.ncomp,eta,Rn,pref,out.data(),count);
}

// ---------------------------------------------------------------- main
int main(int argc,char** argv){
  if(argc<2){fprintf(stderr,"usage: %s routec_tables.bin [l=2] [K=1] [Nq=2000]\n",argv[0]);return 1;}
  RTC T;
  if(!load_tables(argv[1],T)){fprintf(stderr,"cannot load %s\n",argv[1]);return 1;}
  int lcls=argc>2?atoi(argv[2]):2;
  int K=argc>3?atoi(argv[3]):1;
  int Nq=argc>4?atoi(argv[4]):2000;
  printf("Route-C C++ (Phase 2): tables kmax=%d kpair=%d\n",T.kmax,T.kpair);

  int mode=pin_rotc_mode(T);
  if(mode<0){fprintf(stderr,"FATAL: rotc convention does not match Python reference\n");return 1;}

  // ---- anchored MD validation
  double A[3]={0.1,-0.2,0.3},B[3]={0.9,0.5,-0.4},Cc[3]={-0.3,1.1,0.8},D[3]={0.7,-0.6,1.2};
  double a=1.1,b=0.7,c=0.9,d=1.3;
  struct Anc{int lA[3],lB[3],lC[3],lD[3];double ref;};
  Anc anc[4]={
    {{0,0,0},{0,0,0},{0,0,0},{0,0,0},0.19037665793372102},
    {{1,0,0},{0,1,0},{0,0,1},{1,0,0},0.0012685809054877757},
    {{1,1,0},{0,0,2},{1,0,1},{0,2,0},0.0008790037393669891},
    {{1,1,1},{2,1,0},{0,2,1},{3,0,0},-4.283400638256463e-06}};
  double worst=0.0;
  for(auto& q:anc){
    double v=md_reference(q.lA,A,a,q.lB,B,b,q.lC,Cc,c,q.lD,D,d);
    worst=std::max(worst,std::fabs(v-q.ref)/std::max(1e-30,std::fabs(q.ref)));}
  printf("  MD reference vs Python anchors: max rel err %.2e  %s\n",worst,worst<1e-12?"PASS":"FAIL");
  if(worst>=1e-12) return 1;

  // ---- Route-C (spherical output) vs c2s^4-transformed MD validation
  for(int lv=0; lv<=std::min(3,T.kpair/2); ++lv){
    vector<double> blk;
    routec_block(T,mode,lv,A,a,lv,B,b,lv,Cc,c,lv,D,d,blk,false);
    auto ca=monos_of(lv);
    int nc=(int)ca.size(), ns=2*lv+1;
    vector<double> cart((size_t)nc*nc*nc*nc);
    { size_t i=0;
      for(auto& i1:ca)for(auto& i2:ca)for(auto& i3:ca)for(auto& i4:ca)
        cart[i++]=md_reference(i1.data(),A,a,i2.data(),B,b,i3.data(),Cc,c,i4.data(),D,d); }
    // transform each of the 4 indices cart->spherical with c2s (nc x ns)
    const double* Cm=T.c2s[lv].data();
    vector<double> tmp;
    int dims[4]={nc,nc,nc,nc};
    for(int ax=0;ax<4;++ax){
      int d0=1; for(int x2=0;x2<ax;++x2) d0*=dims[x2];
      int d2=1; for(int x2=ax+1;x2<4;++x2) d2*=dims[x2];
      tmp.assign((size_t)d0*ns*d2,0.0);
      for(int o=0;o<d0;++o)
        for(int icq=0;icq<nc;++icq)
          for(int m=0;m<ns;++m){ double cm=Cm[(size_t)icq*ns+m];
            if(cm==0.0) continue;
            const double* src=&cart[((size_t)o*nc+icq)*d2];
            double* dst=&tmp[((size_t)o*ns+m)*d2];
            for(int q2=0;q2<d2;++q2) dst[q2]+=cm*src[q2]; }
      cart.swap(tmp); dims[ax]=ns;
    }
    double w=0.0, scale=0.0;
    for(double v:cart) scale=std::max(scale,std::fabs(v));
    for(size_t i=0;i<cart.size();++i)
      w=std::max(w,std::fabs(blk[i]-cart[i])/std::max(scale,1e-30));
    printf("  Route-C(sph) vs MD(c2s^4), (l=%d)^4 full block: max rel err %.2e  %s\n",
           lv,w,w<1e-11?"PASS":"FAIL");
    if(w>=1e-11) return 1;
  }

  // ---- benchmark: from-scratch and pair-cached
  std::mt19937_64 rng(12345);
  std::uniform_real_distribution<double> jit(-0.05,0.05);
  auto bench=[&](bool cached){
    // K primitives per shell, unit coefficients
    vector<double> exps(K); for(int i2=0;i2<K;++i2) exps[i2]=0.6+0.5*i2;
    double t0=now_s(); int reps=0; double tmin=0.4;
    uint64_t ops_total=0;
    g_ops.zero();
    vector<double> blk, acc;
    while(now_s()-t0<tmin){
      double Aj[3],Bj[3],Cj[3],Dj[3];
      for(int x=0;x<3;++x){Aj[x]=A[x]+jit(rng);Bj[x]=B[x]+jit(rng);
                           Cj[x]=Cc[x]+jit(rng)+1.5;Dj[x]=D[x]+jit(rng)+1.5;}
      bool count = (reps==0);
      if(!cached){
        // full from-scratch: K^4 primitive quartets, everything inside
        acc.clear();
        for(int pa=0;pa<K;++pa)for(int pb=0;pb<K;++pb)for(int pc=0;pc<K;++pc)for(int pd=0;pd<K;++pd){
          routec_block(T,mode,lcls,Aj,exps[pa],lcls,Bj,exps[pb],lcls,Cj,exps[pc],lcls,Dj,exps[pd],blk,count);
          if(acc.empty()) acc.assign(blk.size(),0.0);
          for(size_t i2=0;i2<blk.size();++i2) acc[i2]+=blk[i2]; }
      } else {
        // pair-cached: build+rotate... rotation is per-quartet (frame depends on
        // primitive P,Q) so only the UNROTATED pair build is cacheable.
        // measure: per-quartet = rotate(copy) + H + T + couple.
        static vector<PairF> braP, ketP; braP.clear(); ketP.clear();
        for(int pa=0;pa<K;++pa)for(int pb=0;pb<K;++pb){ PairF P0;
          build_pair(T,lcls,lcls,Aj,Bj,exps[pa],exps[pb],P0,false); braP.push_back(P0);}
        for(int pc=0;pc<K;++pc)for(int pd=0;pd<K;++pd){ PairF P0;
          build_pair(T,lcls,lcls,Cj,Dj,exps[pc],exps[pd],P0,false); ketP.push_back(P0);}
        acc.clear();
        for(int pa=0;pa<K;++pa)for(int pb=0;pb<K;++pb)for(int pc=0;pc<K;++pc)for(int pd=0;pd<K;++pd){
          double aa2=exps[pa],bb=exps[pb],cc2=exps[pc],dd=exps[pd];
          double p=aa2+bb,q=cc2+dd;
          double P[3],Q[3],R[3];
          for(int x=0;x<3;++x){P[x]=(aa2*Aj[x]+bb*Bj[x])/p;Q[x]=(cc2*Cj[x]+dd*Dj[x])/q;R[x]=P[x]-Q[x];}
          double eta=p*q/(p+q), Rn=std::sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
          double pref=2.0*std::pow(M_PI,2.5)/(p*q*std::sqrt(p+q));
          double z[3]={R[0]/Rn,R[1]/Rn,R[2]/Rn};
          double aa3[3]; if(std::fabs(z[0])<0.9){aa3[0]=1;aa3[1]=0;aa3[2]=0;} else {aa3[0]=0;aa3[1]=1;aa3[2]=0;}
          double dot=aa3[0]*z[0]+aa3[1]*z[1]+aa3[2]*z[2];
          double x[3]={aa3[0]-dot*z[0],aa3[1]-dot*z[1],aa3[2]-dot*z[2]};
          double xn=std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
          for(int i2=0;i2<3;++i2)x[i2]/=xn;
          double y[3]={z[1]*x[2]-z[2]*x[1], z[2]*x[0]-z[0]*x[2], z[0]*x[1]-z[1]*x[0]};
          double U[9]={x[0],y[0],z[0], x[1],y[1],z[1], x[2],y[2],z[2]};
          const PairF& Pb=braP[(size_t)pa*K+pb]; const PairF& Pk=ketP[(size_t)pc*K+pd];
          static thread_local QuartetRot QR;
          build_qrot(U,std::max(Pb.K,Pk.K),mode,QR);
          const ClassPlan& CP=get_plan(T,Pb.K,Pk.K);
          static thread_local vector<double> Fbp,Fkp;
          rotate_permute(T,Pb,QR,CP.gposB,Fbp,count);
          rotate_permute(T,Pk,QR,CP.gposK,Fkp,count);
          if(acc.empty()) acc.assign((size_t)Pb.ncomp*Pk.ncomp,0.0);
          blk.assign((size_t)Pb.ncomp*Pk.ncomp,0.0);
          couple_planned(T,CP,Fbp.data(),Pb.ncomp,Fkp.data(),Pk.ncomp,eta,Rn,pref,blk.data(),count);
          for(size_t i2=0;i2<blk.size();++i2) acc[i2]+=blk[i2]; }
      }
      if(reps==0) ops_total=g_ops.total();
      ++reps;
    }
    double el=now_s()-t0;
    double per_contracted = el/reps;            // one contracted quartet (K^4 prims)
    double per_prim = per_contracted/std::pow((double)K,4);
    printf("  %-12s l=%d K=%d: %8.2f us/contracted-quartet  %8.2f us/prim-quartet  "
           "%9.0f prim-q/s\n",
           cached?"pair-cached":"from-scratch",lcls,K,1e6*per_contracted,1e6*per_prim,1.0/per_prim);
    if(!cached){
      double kk=std::pow((double)K,4);
      printf("    MACs/prim-quartet: pair %.0f  rot %.0f  H %.0f  T %.0f  gemm %.0f  total %.0f\n",
        g_ops.pair/kk,g_ops.rot/kk,g_ops.hbuild/kk,g_ops.tbuild/kk,g_ops.gemm/kk,(double)ops_total/kk);
    }
  };
  printf("benchmark (single thread; pin externally with taskset):\n");
  bench(false);
  bench(true);
  // ---- quartet-only: pair tensors prebuilt OUTSIDE the timer (production
  // semantics; matches the libcint baseline which runs with CINTOpt pair data)
  {
    const int NPOOL=64;
    double exp_a=0.9,exp_b=1.4,exp_c=0.7,exp_d=1.1;
    PairF bra; build_pair(T,lcls,lcls,A,B,exp_a,exp_b,bra,false);
    vector<PairF> kets(NPOOL);
    vector<std::array<double,6>> ketgeo(NPOOL);
    for(int i2=0;i2<NPOOL;++i2){
      double Cj[3],Dj[3];
      for(int x=0;x<3;++x){Cj[x]=Cc[x]+jit(rng)*8+1.5;Dj[x]=D[x]+jit(rng)*8+1.5;}
      build_pair(T,lcls,lcls,Cj,Dj,exp_c,exp_d,kets[i2],false);
      ketgeo[i2]={Cj[0],Cj[1],Cj[2],Dj[0],Dj[1],Dj[2]};
    }
    double p=exp_a+exp_b,q=exp_c+exp_d;
    double P[3]; for(int x=0;x<3;++x)P[x]=(exp_a*A[x]+exp_b*B[x])/p;
    double eta=p*q/(p+q);
    double pref=2.0*std::pow(M_PI,2.5)/(p*q*std::sqrt(p+q));
    const ClassPlan& CP=get_plan(T,bra.K,bra.K);
    static thread_local vector<double> Fbp,Fkp; vector<double> blk;
    g_prof.on = (getenv("ROUTEC_PROF")!=nullptr); g_prof.zero();
    double t0=now_s(); long iq=0; g_ops.zero(); uint64_t ops1=0;
    while(now_s()-t0<0.6){
      const PairF& ket=kets[iq%NPOOL];
      const auto& gg2=ketgeo[iq%NPOOL];
      double tq0=g_prof.on?now_s():0.0;
      double Q[3]={(exp_c*gg2[0]+exp_d*gg2[3])/q,(exp_c*gg2[1]+exp_d*gg2[4])/q,
                   (exp_c*gg2[2]+exp_d*gg2[5])/q};
      double R[3]={P[0]-Q[0],P[1]-Q[1],P[2]-Q[2]};
      double Rn=std::sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
      double z[3]={R[0]/Rn,R[1]/Rn,R[2]/Rn};
      double aa3[3]; if(std::fabs(z[0])<0.9){aa3[0]=1;aa3[1]=0;aa3[2]=0;} else {aa3[0]=0;aa3[1]=1;aa3[2]=0;}
      double dot=aa3[0]*z[0]+aa3[1]*z[1]+aa3[2]*z[2];
      double x[3]={aa3[0]-dot*z[0],aa3[1]-dot*z[1],aa3[2]-dot*z[2]};
      double xn=std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      for(int i2=0;i2<3;++i2)x[i2]/=xn;
      double y[3]={z[1]*x[2]-z[2]*x[1], z[2]*x[0]-z[0]*x[2], z[0]*x[1]-z[1]*x[0]};
      double U[9]={x[0],y[0],z[0], x[1],y[1],z[1], x[2],y[2],z[2]};
      static thread_local QuartetRot QR;
      build_qrot(U,bra.K,mode,QR);
      double tq1=g_prof.on?now_s():0.0;
      bool count=(iq==0);
      rotate_permute(T,bra,QR,CP.gposB,Fbp,count);
      double tq2=g_prof.on?now_s():0.0;
      rotate_permute(T,ket,QR,CP.gposK,Fkp,count);
      double tq3=g_prof.on?now_s():0.0;
      if(g_prof.on){ g_prof.frame+=tq1-tq0; g_prof.rotb+=tq2-tq1; g_prof.rotk+=tq3-tq2; }
      blk.assign((size_t)bra.ncomp*ket.ncomp,0.0);
      couple_planned(T,CP,Fbp.data(),bra.ncomp,Fkp.data(),ket.ncomp,eta,Rn,pref,blk.data(),count);
      if(iq==0) ops1=g_ops.total();
      ++iq;
    }
    double el=now_s()-t0;
    printf("  %-12s l=%d    : %8.2f us/quartet  %9.0f q/s   (MACs/q: %llu)\n",
           "quartet-only",lcls,1e6*el/iq,iq/el,(unsigned long long)ops1);
    if(g_prof.on){
      double tot=g_prof.frame+g_prof.rotb+g_prof.rotk+g_prof.rt+g_prof.h0+g_prof.grp;
      printf("    stage profile (us/q, %% of instrumented):\n");
      auto pr=[&](const char* nm,double v){
        printf("      %-10s %7.3f  %5.1f%%\n",nm,1e6*v/iq,100.0*v/tot); };
      pr("frame+qrot",g_prof.frame); pr("rot bra",g_prof.rotb); pr("rot ket",g_prof.rotk);
      pr("r_tensor",g_prof.rt); pr("H0 fold",g_prof.h0); pr("mu-groups",g_prof.grp);
      printf("      %-10s %7.3f  (instrumented sum vs %0.3f total)\n","SUM",
             1e6*tot/iq,1e6*el/iq);
    }
  }
  // MD reference timing for in-binary sanity
  {
    auto ca=monos_of(lcls); int nc=(int)ca.size();
    double t0=now_s(); int reps=0;
    while(now_s()-t0<0.4){
      double Aj[3],Bj[3],Cj[3],Dj[3];
      for(int x=0;x<3;++x){Aj[x]=A[x]+jit(rng);Bj[x]=B[x]+jit(rng);Cj[x]=Cc[x]+jit(rng)+1.5;Dj[x]=D[x]+jit(rng)+1.5;}
      double s2=0.0;
      for(auto& i1:ca)for(auto& i2:ca)for(auto& i3:ca)for(auto& i4:ca)
        s2+=md_reference(i1.data(),Aj,a,i2.data(),Bj,b,i3.data(),Cj,c,i4.data(),Dj,d);
      (void)s2; ++reps; }
    double el=now_s()-t0;
    printf("  plain-MD ref l=%d K=1: %8.2f us/quartet (%9.0f q/s) [unoptimized reference]\n",
      lcls,1e6*el/reps,reps/el);
  }
  (void)Nq;
  return 0;
}
