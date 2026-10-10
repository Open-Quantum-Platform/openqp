// M3: device assembly of the spherical-resolution ERI (Route A, in-frame):
//   (ab|cd) = pref * Re sum_{t,t'} E_ab[t] E_cd[t'] i^{|t|-|t'|}
//                       sum_{l,m} i^l A^{t+t'}_{lm} g_l^{|t+t'|}(eta,R) Y_lm(R^)
// ties M1 (radial), M2 (Hermite E), Y_lm, and the Gaunt table. Validated vs MD.
#pragma once
#include "sph_eri/port.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include "sph_eri/radial_cuda.cuh"
#include "sph_eri/hermite_cuda.cuh"
#include "sph_eri/ylm_cuda.cuh"

namespace sph_eri {

struct GauntTable {
  int maxdeg = 0;
  const int* d_start = nullptr;  // encode(T) -> first entry index, or -1
  const int* d_count = nullptr;  // encode(T) -> number of (l,m) entries
  const int* d_l = nullptr;
  const int* d_m = nullptr;
  const double* d_A = nullptr;
};

__host__ __device__ inline int gaunt_encode(int tx,int ty,int tz,int md){ return (tx*(md+1)+ty)*(md+1)+tz; }

inline GauntTable gaunt_load(const char* path) {
  FILE* f=std::fopen(path,"rb"); if(!f){std::fprintf(stderr,"gaunt_load: open %s\n",path);std::exit(1);}
  int md=0, ne=0; if(std::fread(&md,4,1,f)!=1||std::fread(&ne,4,1,f)!=1){std::exit(1);}
  struct E{int tx,ty,tz,l,m;double A;};
  std::vector<E> ent(ne);
  for(int i=0;i<ne;++i){ int t[5]; double A; if(std::fread(t,4,5,f)!=5||std::fread(&A,8,1,f)!=1){std::exit(1);} ent[i]={t[0],t[1],t[2],t[3],t[4],A}; }
  std::fclose(f);
  // group by encode(T)
  std::sort(ent.begin(),ent.end(),[&](const E&a,const E&b){return gaunt_encode(a.tx,a.ty,a.tz,md)<gaunt_encode(b.tx,b.ty,b.tz,md);});
  int ncell=(md+1)*(md+1)*(md+1);
  std::vector<int> start(ncell,-1), count(ncell,0), Ls(ne), Ms(ne); std::vector<double> As(ne);
  for(int i=0;i<ne;++i){ int key=gaunt_encode(ent[i].tx,ent[i].ty,ent[i].tz,md); if(start[key]<0)start[key]=i; count[key]++; Ls[i]=ent[i].l; Ms[i]=ent[i].m; As[i]=ent[i].A; }
  GauntTable g; g.maxdeg=md;
#ifdef __CUDACC__
  int *ds,*dc,*dl,*dm; double* dA;
  cudaMalloc(&ds,ncell*4);cudaMemcpy(ds,start.data(),ncell*4,cudaMemcpyHostToDevice);
  cudaMalloc(&dc,ncell*4);cudaMemcpy(dc,count.data(),ncell*4,cudaMemcpyHostToDevice);
  cudaMalloc(&dl,ne*4);cudaMemcpy(dl,Ls.data(),ne*4,cudaMemcpyHostToDevice);
  cudaMalloc(&dm,ne*4);cudaMemcpy(dm,Ms.data(),ne*4,cudaMemcpyHostToDevice);
  cudaMalloc(&dA,ne*8);cudaMemcpy(dA,As.data(),ne*8,cudaMemcpyHostToDevice);
  g.d_start=ds;g.d_count=dc;g.d_l=dl;g.d_m=dm;g.d_A=dA;
#else                                   // CPU backend: host-resident table
  int* ds=(int*)std::malloc(ncell*4); std::memcpy(ds,start.data(),ncell*4);
  int* dc=(int*)std::malloc(ncell*4); std::memcpy(dc,count.data(),ncell*4);
  int* dl=(int*)std::malloc(ne*4); std::memcpy(dl,Ls.data(),ne*4);
  int* dm=(int*)std::malloc(ne*4); std::memcpy(dm,Ms.data(),ne*4);
  double* dA=(double*)std::malloc(ne*8); std::memcpy(dA,As.data(),ne*8);
  g.d_start=ds;g.d_count=dc;g.d_l=dl;g.d_m=dm;g.d_A=dA;
#endif
  return g;
}

// i^k as (re,im); k may be negative
__host__ __device__ inline void ipow(int k, double& re, double& im){ int r=((k%4)+4)%4; re=(r==0)?1:(r==2)?-1:0; im=(r==1)?1:(r==3)?-1:0; }

// One primitive Cartesian-component ERI with all four shells of momentum (lx,ly,lz)=L*.
template <typename T>
__host__ __device__ inline T eri_assemble(int Ax,int Ay,int Az, int Bx,int By,int Bz,
                                 int Cx,int Cy,int Cz, int Dx,int Dy,int Dz,
                                 const T* A,const T a,const T* B,const T b,
                                 const T* C,const T c,const T* D,const T d,
                                 const RadialTable& rt, const GauntTable& gt,
                                 const T* gtab=nullptr, const T* Ytab=nullptr, int Lmax=0,
                                 bool canon=false) {  // canon: R is along z (only m=0 survives) -> skip m!=0
  // NOTE: gtab/Ytab are read in assembly precision T -- pass FP32 tables with T=float
  // to run the dominant angular sum (below) entirely in FP32 (FP64 accumulate kept).
  T p=a+b, q=c+d, eta=p*q/(p+q);
  T P[3],Q[3],Rv[3];
  for(int i=0;i<3;++i){ P[i]=(a*A[i]+b*B[i])/p; Q[i]=(c*C[i]+d*D[i])/q; Rv[i]=P[i]-Q[i]; }
  T R=(T)sqrt((double)(Rv[0]*Rv[0]+Rv[1]*Rv[1]+Rv[2]*Rv[2]));
  T r=fmax(R,(T)1e-300);
  T thR=(T)acos(fmax(-1.0,fmin(1.0,(double)(Rv[2]/r))));
  T phR=(T)atan2((double)Rv[1],(double)Rv[0]);
  // Hermite per axis (E_ab uses Q_ax=A-B; E_cd uses C-D)
  T ex[14],ey[14],ez[14], fx[14],fy[14],fz[14];
  hermite_E_axis<T>(Ax,Bx,A[0]-B[0],a,b,ex); hermite_E_axis<T>(Ay,By,A[1]-B[1],a,b,ey); hermite_E_axis<T>(Az,Bz,A[2]-B[2],a,b,ez);
  hermite_E_axis<T>(Cx,Dx,C[0]-D[0],c,d,fx); hermite_E_axis<T>(Cy,Dy,C[1]-D[1],c,d,fy); hermite_E_axis<T>(Cz,Dz,C[2]-D[2],c,d,fz);
  int axb=Ax+Bx,ayb=Ay+By,azb=Az+Bz, cxd=Cx+Dx,cyd=Cy+Dy,czd=Cz+Dz;
  double totre=0, totim=0;
  for(int tx=0;tx<=axb;++tx)for(int ty=0;ty<=ayb;++ty)for(int tz=0;tz<=azb;++tz){
    T Eab=ex[tx]*ey[ty]*ez[tz]; if(Eab==(T)0) continue; int dt=tx+ty+tz;
    for(int ux=0;ux<=cxd;++ux)for(int uy=0;uy<=cyd;++uy)for(int uz=0;uz<=czd;++uz){
      T Ecd=fx[ux]*fy[uy]*fz[uz]; if(Ecd==(T)0) continue; int du=ux+uy+uz;
      int Tx=tx+ux,Ty=ty+uy,Tz=tz+uz, nT=Tx+Ty+Tz;
      // phase i^{dt-du}
      double pr,pi; ipow(dt-du,pr,pi);
      double w=(double)(Eab*Ecd);
      double cre=w*pr, cim=w*pi;  // c-coefficient for this (t,t')
      // angular sum S(T) = sum_{l,m} i^l A g_l Y
      int key=gaunt_encode(Tx,Ty,Tz,gt.maxdeg); if(Tx>gt.maxdeg||Ty>gt.maxdeg||Tz>gt.maxdeg) continue;
      int s=gt.d_start[key], cnt=gt.d_count[key];
      for(int e=0;e<cnt;++e){
        int l=gt.d_l[s+e], m=gt.d_m[s+e]; double Almv=gt.d_A[s+e];
        if(canon && m!=0) continue;     // canonical frame: Y_lm(z)=0 for m!=0
        // precomputed tables (per primitive quartet) when provided, else inline
        T g = gtab ? gtab[l*(Lmax+1)+nT]            : radial_g<T>(rt,l,nT,eta,R);
        T Y = Ytab ? Ytab[l*(2*Lmax+1)+(m+Lmax)]    : real_sph_harm<T>(l,m,thR,phR);
        double ilr,ili; ipow(l,ilr,ili);
        T s0=(T)Almv*g*Y;                          // A g Y in assembly precision T
        double sre=(double)s0*ilr, sim=(double)s0*ili; // i^l A g Y (FP64 accumulate)
        // total += c * S  (complex multiply)
        totre += cre*sre - cim*sim;
        totim += cre*sim + cim*sre;
      }
    }
  }
  double pref=(2.0/M_PI)*pow(M_PI/(double)p,1.5)*pow(M_PI/(double)q,1.5);
  return (T)(pref*totre);
}

// Like eri_contract_RT but reads PRECOMPUTED per-axis Hermite E (ex..fz, already
// built once per primitive quartet and shared across all Cartesian components) --
// removing the dominant per-component Hermite recomputation. axb=Ax+Bx etc bound
// the sums; ex[t] valid for t=0..axb, fx[u] for u=0..cxd, etc.
template <typename T>
__host__ __device__ inline T eri_contract_RT_pre(
    int axb,int ayb,int azb,int cxd,int cyd,int czd,
    const T* ex,const T* ey,const T* ez,const T* fx,const T* fy,const T* fz,
    const T* RTre, const T* RTim, int Lmax, T p, T q) {
  int L1=Lmax+1; double totre=0;
  for(int tx=0;tx<=axb;++tx)for(int ty=0;ty<=ayb;++ty)for(int tz=0;tz<=azb;++tz){
    T Eab=ex[tx]*ey[ty]*ez[tz]; if(Eab==(T)0) continue; int dt=tx+ty+tz;
    for(int ux=0;ux<=cxd;++ux)for(int uy=0;uy<=cyd;++uy)for(int uz=0;uz<=czd;++uz){
      T Ecd=fx[ux]*fy[uy]*fz[uz]; if(Ecd==(T)0) continue; int du=ux+uy+uz;
      int Ti=((tx+ux)*L1+(ty+uy))*L1+(tz+uz);
      double pr,pi; ipow(dt-du,pr,pi);
      double w=(double)(Eab*Ecd);
      totre += (w*pr)*(double)RTre[Ti] - (w*pi)*(double)RTim[Ti];
    }
  }
  double pref=(2.0/M_PI)*pow(M_PI/(double)p,1.5)*pow(M_PI/(double)q,1.5);
  return (T)(pref*totre);
}

// MD-hoist: the angular sum S(T)=sum_{l,m} i^l A^T_lm g_l^(|T|) Y_lm depends ONLY on the
// Hermite index T=(Tx,Ty,Tz) and the geometry/exponents, NOT on the Cartesian component.
// So precompute it once per primitive quartet (RTre/RTim, indexed Ti=(Tx*(Lmax+1)+Ty)*(Lmax+1)+Tz),
// then each component is just the cheap Hermite-E contraction below -- removing the
// per-component Gaunt/g/Y recomputation. Works for CONTRACTED shells (unlike Route-B).
template <typename T>
__host__ __device__ inline T eri_contract_RT(
    int Ax,int Ay,int Az,int Bx,int By,int Bz,int Cx,int Cy,int Cz,int Dx,int Dy,int Dz,
    const T* A,const T a,const T* B,const T b,const T* C,const T c,const T* D,const T d,
    const T* RTre, const T* RTim, int Lmax, T p, T q) {
  T ex[14],ey[14],ez[14], fx[14],fy[14],fz[14];
  hermite_E_axis<T>(Ax,Bx,A[0]-B[0],a,b,ex); hermite_E_axis<T>(Ay,By,A[1]-B[1],a,b,ey); hermite_E_axis<T>(Az,Bz,A[2]-B[2],a,b,ez);
  hermite_E_axis<T>(Cx,Dx,C[0]-D[0],c,d,fx); hermite_E_axis<T>(Cy,Dy,C[1]-D[1],c,d,fy); hermite_E_axis<T>(Cz,Dz,C[2]-D[2],c,d,fz);
  int axb=Ax+Bx,ayb=Ay+By,azb=Az+Bz, cxd=Cx+Dx,cyd=Cy+Dy,czd=Cz+Dz, L1=Lmax+1;
  double totre=0;
#ifdef CONTRACT_NOTU
  return (T)(ex[0]*fx[0]);   // ablation: skip the (t,u) loop, keep Hermite cost only
#endif
  for(int tx=0;tx<=axb;++tx)for(int ty=0;ty<=ayb;++ty)for(int tz=0;tz<=azb;++tz){
    T Eab=ex[tx]*ey[ty]*ez[tz]; if(Eab==(T)0) continue; int dt=tx+ty+tz;
    for(int ux=0;ux<=cxd;++ux)for(int uy=0;uy<=cyd;++uy)for(int uz=0;uz<=czd;++uz){
      T Ecd=fx[ux]*fy[uy]*fz[uz]; if(Ecd==(T)0) continue; int du=ux+uy+uz;
      int Ti=((tx+ux)*L1+(ty+uy))*L1+(tz+uz);
      double pr,pi; ipow(dt-du,pr,pi);
      double w=(double)(Eab*Ecd);
      totre += (w*pr)*(double)RTre[Ti] - (w*pi)*(double)RTim[Ti];
    }
  }
  double pref=(2.0/M_PI)*pow(M_PI/(double)p,1.5)*pow(M_PI/(double)q,1.5);
  return (T)(pref*totre);
}

// Precompute the geometry-only g_l^(nT)(eta,R) and Y_lm(theta,phi) for ONE
// primitive quartet, for all l,nT,m in [0,Lmax] (Lmax = la+lb+lc+ld).  These are
// identical for every Cartesian component of the quartet, so computing them once
// here and passing them to eri_assemble removes the per-component recomputation
// of the special functions (the J/K hot-loop cost).  Layout:
//   gtab[l*(Lmax+1)+nT],  Ytab[l*(2*Lmax+1)+(m+Lmax)].
template <typename T>
__host__ __device__ inline void precompute_gY(
    const T* A,const T a,const T* B,const T b,const T* C,const T c,const T* D,const T d,
    const RadialTable& rt, int Lmax, double* gtab, double* Ytab){
  T p=a+b, q=c+d, eta=p*q/(p+q);
  T P[3],Q[3],Rv[3];
  for(int i=0;i<3;++i){ P[i]=(a*A[i]+b*B[i])/p; Q[i]=(c*C[i]+d*D[i])/q; Rv[i]=P[i]-Q[i]; }
  T R=(T)sqrt((double)(Rv[0]*Rv[0]+Rv[1]*Rv[1]+Rv[2]*Rv[2]));
  T r=fmax(R,(T)1e-300);
  T thR=(T)acos(fmax(-1.0,fmin(1.0,(double)(Rv[2]/r))));
  T phR=(T)atan2((double)Rv[1],(double)Rv[0]);
  for(int l=0;l<=Lmax;++l){
    for(int nT=0;nT<=Lmax;++nT) gtab[l*(Lmax+1)+nT]=(double)radial_g<T>(rt,l,nT,eta,R);
    for(int m=-l;m<=l;++m)       Ytab[l*(2*Lmax+1)+(m+Lmax)]=(double)real_sph_harm<T>(l,m,thR,phR);
  }
}

// Per shell-PAIR Hermite E tables (the reusable "cloud-shape" data).  For a pair
// (la,lb) with displacement disp=A-B and exponents a,b, fills 3 axis tables
//   E[axis][(i*(lb+1)+j)*(la+lb+1)+t]  for i=0..la, j=0..lb, t=0..i+j.
// Computed once per primitive pair and reused across every quartet containing it.
template <typename T>
__host__ __device__ inline void precompute_hermite(int la,int lb,const T* disp,T a,T b,
                                                   T* Ex,T* Ey,T* Ez){
  int W=la+lb+1;
  for(int i=0;i<=la;++i)for(int j=0;j<=lb;++j){
    T tx[14],ty[14],tz[14];
    hermite_E_axis<T>(i,j,disp[0],a,b,tx); hermite_E_axis<T>(i,j,disp[1],a,b,ty); hermite_E_axis<T>(i,j,disp[2],a,b,tz);
    int base=(i*(lb+1)+j)*W;
    for(int t=0;t<=i+j;++t){ Ex[base+t]=tx[t]; Ey[base+t]=ty[t]; Ez[base+t]=tz[t]; }
  }
}

// Fully-precomputed Cartesian-component ERI: NO special functions, NO Hermite
// recursion at call time -- reads the per-pair Hermite tables (Ex..Fz), the
// per-quartet radial/Y_lm tables (gtab,Ytab), and the global Gaunt table.  This
// is the shell-pair-driven hot kernel: everything expensive was precomputed.
template <typename T>
__host__ __device__ inline T eri_assemble_pre(
    int Ax,int Ay,int Az,int Bx,int By,int Bz,int Cx,int Cy,int Cz,int Dx,int Dy,int Dz,
    int la,int lb,int lc,int ld, double p, double q,
    const T* Ex,const T* Ey,const T* Ez, const T* Fx,const T* Fy,const T* Fz,
    const GauntTable& gt, const double* gtab,const double* Ytab,int Lmax){
  int Wb=la+lb+1, Wk=lc+ld+1;
  const T* ex=&Ex[(Ax*(lb+1)+Bx)*Wb]; const T* ey=&Ey[(Ay*(lb+1)+By)*Wb]; const T* ez=&Ez[(Az*(lb+1)+Bz)*Wb];
  const T* fx=&Fx[(Cx*(ld+1)+Dx)*Wk]; const T* fy=&Fy[(Cy*(ld+1)+Dy)*Wk]; const T* fz=&Fz[(Cz*(ld+1)+Dz)*Wk];
  int axb=Ax+Bx,ayb=Ay+By,azb=Az+Bz, cxd=Cx+Dx,cyd=Cy+Dy,czd=Cz+Dz;
  double totre=0;
  for(int tx=0;tx<=axb;++tx)for(int ty=0;ty<=ayb;++ty)for(int tz=0;tz<=azb;++tz){
    T Eab=ex[tx]*ey[ty]*ez[tz]; if(Eab==(T)0) continue; int dt=tx+ty+tz;
    for(int ux=0;ux<=cxd;++ux)for(int uy=0;uy<=cyd;++uy)for(int uz=0;uz<=czd;++uz){
      T Ecd=fx[ux]*fy[uy]*fz[uz]; if(Ecd==(T)0) continue; int du=ux+uy+uz;
      int Tx=tx+ux,Ty=ty+uy,Tz=tz+uz, nT=Tx+Ty+Tz;
      double pr,pi; ipow(dt-du,pr,pi);
      double w=(double)(Eab*Ecd), cre=w*pr, cim=w*pi;
      if(Tx>gt.maxdeg||Ty>gt.maxdeg||Tz>gt.maxdeg) continue;
      int key=gaunt_encode(Tx,Ty,Tz,gt.maxdeg), s=gt.d_start[key], cnt=gt.d_count[key];
      for(int e=0;e<cnt;++e){
        int l=gt.d_l[s+e], m=gt.d_m[s+e]; double Almv=gt.d_A[s+e];
        double g=gtab[l*(Lmax+1)+nT], Y=Ytab[l*(2*Lmax+1)+(m+Lmax)];
        double ilr,ili; ipow(l,ilr,ili);
        double sre=Almv*g*Y*ilr, sim=Almv*g*Y*ili;
        totre += cre*sre - cim*sim;
      }
    }
  }
  double pref=(2.0/M_PI)*pow(M_PI/p,1.5)*pow(M_PI/q,1.5);
  return (T)(pref*totre);
}

// Contracted Cartesian-component ERI: sum over primitives with contraction coeffs.
template <typename T>
__host__ __device__ inline T eri_contracted(int lx,int ly,int lz,
    const T* A,const T* B,const T* C,const T* D,
    const T* ea,const T* ca,int Ka, const T* eb,const T* cb,int Kb,
    const T* ec,const T* cc,int Kc, const T* ed,const T* cd,int Kd,
    const RadialTable& rt, const GauntTable& gt) {
  T s=(T)0;
  for(int i=0;i<Ka;++i)for(int j=0;j<Kb;++j)for(int k=0;k<Kc;++k)for(int l=0;l<Kd;++l)
    s += ca[i]*cb[j]*cc[k]*cd[l]*eri_assemble<T>(lx,ly,lz,lx,ly,lz,lx,ly,lz,lx,ly,lz,
            A,ea[i],B,eb[j],C,ec[k],D,ed[l], rt,gt);
  return s;
}

} // namespace sph_eri
