// Real spherical-harmonic rotation matrices R^l by direct recursion (Choi/
// Ivanic-Ruedenberg), __host__ __device__.  Given a 3x3 Cartesian rotation R
// (row-major), build R^1..R^L, each (2l+1)x(2l+1) indexed [m+l, n+l], with
//   Y_lm(R n) = sum_{m'} R^l_{m m'} Y_lm'(n).
// R^1 is R reordered to the (y,z,x) = (m=-1,0,+1) basis.  This is the rotation
// that carries all orientation dependence in the spherical-resolution ERI.
#pragma once
#include "sph_eri/port.h"
#include <cmath>

namespace sph_eri {

template <typename T> __host__ __device__ inline T rotc_r1(const T* R1, int p, int q){
  return R1[(p+1)*3 + (q+1)];                 // p,q in {-1,0,1}
}
template <typename T> __host__ __device__ inline T rotc_rm(const T* Rlm1, int l, int p, int q){
  int lm1=l-1; if(p<-lm1||p>lm1||q<-lm1||q>lm1) return (T)0;
  int s=2*lm1+1; return Rlm1[(p+lm1)*s + (q+lm1)];
}
template <typename T> __host__ __device__ inline T rotc_P(int i,int a,int b,int l,const T* R1,const T* Rlm1){
  if(b==l)  return rotc_r1<T>(R1,i,1)*rotc_rm<T>(Rlm1,l,a,l-1) - rotc_r1<T>(R1,i,-1)*rotc_rm<T>(Rlm1,l,a,-(l-1));
  if(b==-l) return rotc_r1<T>(R1,i,1)*rotc_rm<T>(Rlm1,l,a,-(l-1)) + rotc_r1<T>(R1,i,-1)*rotc_rm<T>(Rlm1,l,a,l-1);
  return rotc_r1<T>(R1,i,0)*rotc_rm<T>(Rlm1,l,a,b);
}
// Build R^l (size (2l+1)^2) from R^1 and R^{l-1} into Rl (row-major [m+l,n+l]).
template <typename T> __host__ __device__ inline void rotc_next(const T* R1,const T* Rlm1,int l,T* Rl){
  int sz=2*l+1;
  for(int m=-l;m<=l;++m)for(int n=-l;n<=l;++n){
    int d0=(m==0)?1:0;
    T denom = (abs(n)<l) ? (T)((l+n)*(l-n)) : (T)((2*l)*(2*l-1));
    T u = (T)sqrt((double)((l+m)*(l-m))/(double)denom);
    T v = (T)0.5*(T)sqrt((double)((1+d0)*(l+abs(m)-1)*(l+abs(m)))/(double)denom)*(T)(1-2*d0);
    T w = (T)(-0.5)*(T)sqrt((double)((l-abs(m)-1)*(l-abs(m)))/(double)denom);
    T U = rotc_P<T>(0,m,n,l,R1,Rlm1);
    T V,W;
    if(m==0){ V = rotc_P<T>(1,1,n,l,R1,Rlm1) + rotc_P<T>(-1,-1,n,l,R1,Rlm1); W=(T)0; }
    else if(m>0){
      int dm1=(m==1)?1:0;
      V = rotc_P<T>(1,m-1,n,l,R1,Rlm1)*(T)sqrt((double)(1+dm1)) - rotc_P<T>(-1,-m+1,n,l,R1,Rlm1)*(T)(1-dm1);
      W = rotc_P<T>(1,m+1,n,l,R1,Rlm1) + rotc_P<T>(-1,-m-1,n,l,R1,Rlm1);
    } else {
      int dm1=(m==-1)?1:0;
      V = rotc_P<T>(1,m+1,n,l,R1,Rlm1)*(T)(1-dm1) + rotc_P<T>(-1,-m-1,n,l,R1,Rlm1)*(T)sqrt((double)(1+dm1));
      W = rotc_P<T>(1,m-1,n,l,R1,Rlm1) - rotc_P<T>(-1,-m+1,n,l,R1,Rlm1);
    }
    Rl[(m+l)*sz+(n+l)] = u*U + v*V + w*W;
  }
}
// R1 from the 3x3 Cartesian R (row-major), reordered to (y,z,x): perm=(1,2,0).
template <typename T> __host__ __device__ inline void rotc_R1(const T* R3,T* R1){
  const int p[3]={1,2,0};
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) R1[i*3+j]=R3[p[i]*3+p[j]];
}
// Build R^l for a single order l (<=8) into Rl[(2l+1)^2]; needs R3[9].
template <typename T> __host__ __device__ inline void rotc_order(const T* R3,int l,T* Rl){
  if(l==0){ Rl[0]=(T)1; return; }
  T R1buf[9]; rotc_R1<T>(R3, R1buf);
  if(l==1){ for(int i=0;i<9;++i)Rl[i]=R1buf[i]; return; }
  T a[289], b[289];                 // scratch up to l=8 -> (17)^2
  for(int i=0;i<9;++i) a[i]=R1buf[i];   // prev = R^1
  T *prev=a,*cur=b;
  for(int ll=2; ll<=l; ++ll){ rotc_next<T>(R1buf,prev,ll,cur); T* t=prev; prev=cur; cur=t; }
  int sz=2*l+1; for(int i=0;i<sz*sz;++i) Rl[i]=prev[i];
}

} // namespace sph_eri
