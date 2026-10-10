// =============================================================================
// routec_ownxc_gga_vxc.cu  --  OpenQP Route-C OWN GPU semilocal-XC backend (rung-2).
//
// Extends the rung-1 CUDA LDA backend (Slater + VWN5) to GGA + the flagship
// hybrids.  Hand-coded functionals, NO libxc link (no symbol clash with
// liboqp's own libxc).
//
// Exports ONE bind(C) symbol consumed by OpenQP's routec_bridge:
//   void routec_vxc(const double* d, double* fxc, const int* nbf,
//                   const int* nfocks, double* eexc, double* totele, int* info);
//
// Identity target: bit-for-bit at the SAME grid as native OpenQP, component by
// component.  Mirrors native dft_gridint_energy.F90 update() / compDRhoAO /
// basis_tools.F90 compAOvg exactly.
//
// FUNCTIONALS (OQP_OWNXC_FUNC):
//   SLATER  : XC_LDA_X                                         (rung-1)
//   SVWN    : Slater + VWN5                                    (rung-1)
//   BLYP    : B88 exchange + LYP correlation                  (GGA, rung-2)
//   B3LYP   : 0.08 Slater + 0.72 B88 + 0.19 VWN_RPA + 0.81 LYP (hybrid semilocal)
//   BHHLYP  : 0.5 B88 + LYP                                    (hybrid semilocal)
//   (the HF-exchange fraction of B3LYP/BHHLYP is carried by fock_jk, NOT here)
//
// SPIN CONVENTION (closed-shell RHF, hasBeta=false), mirrors native:
//   rho(1)=rho_a=0.5*rho_tot, rho(2)=rho_b=rho_a.
//   drhoa(j) = sum_u aoG1[u,j]*X[u]   (X = D_tot @ phi == native moVA; NO factor 2)
//   drhob = drhoa.
//   sigma_aa = drhoa.drhoa, sigma_ab = drhoa.drhob, sigma_bb = drhob.drhob (all equal).
//   Functionals evaluated POLARIZED (libxc convention): inputs (rho_a,rho_b,
//   sigma_aa,sigma_ab,sigma_bb); outputs eps (per electron of TOTAL),
//   d1dr(ra)=dE/drho_a, d1ds(ga)=dE/dsigma_aa, d1ds(gc)=dE/dsigma_ab.
//
// CONTRACTION (native update(), alpha branch only since hasBeta=false):
//   Z[u] = 0.5*w*d1dr(ra)*phi[u]                                 (LDA part)
//        + w*( c.aoG1 ),  c = 2*d1ds(ga)*drhoa + d1ds(gc)*drhob  (GGA part)
//   Vfull = phi*Z^T + Z*phi^T ; fold bfnrm; pack lower-tri OQP frame.
//   (weight w folded into d1dr/d1ds here, as native scalexc folds it before update.)
// =============================================================================
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#define CK(x) do{ cudaError_t e=(x); if(e!=cudaSuccess){ \
  fprintf(stderr,"[ownxc] CUDA %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
  return false; } }while(0)
#define CKv(x) do{ cudaError_t e=(x); if(e!=cudaSuccess){ \
  fprintf(stderr,"[ownxc] CUDA %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); } }while(0)

// ---- functional selection ---------------------------------------------------
enum { FUNC_SLATER=0, FUNC_SVWN=1, FUNC_BLYP=2, FUNC_B3LYP=3, FUNC_BHHLYP=4 };
__host__ __device__ __forceinline__ bool func_is_gga(int f){ return f>=FUNC_BLYP; }

// ---- persistent device state ------------------------------------------------
struct Ctx {
  bool ready=false;
  int  func=FUNC_SLATER;
  int  nbf=0, nshell=0, nprim=0, npts=0, blk=0;  // blk = phi/gphi sizing (npts when cached, bwork when tiled)
  int  bwork=0;           // working-block size: X/Z/rho/... are [nbf x bwork] ALWAYS (even when phi is cached,
                          // the loop walks the cached phi in bwork slices -> no [nbf x npts] work buffers)
  bool phi_cached=false;  // small systems: collocate the full grid ONCE (fast); large: re-collocate per block per iter
  bool uks_alloc=false;   // beta-channel buffers are allocated lazily on the first UKS call (RKS never pays)
  double *d_gx=nullptr, *d_gy=nullptr, *d_gz=nullptr, *d_gw=nullptr;
  int    *d_sh_am=nullptr, *d_sh_g0=nullptr, *d_sh_nc=nullptr, *d_sh_ao=nullptr;
  double *d_sh_cx=nullptr, *d_sh_cy=nullptr, *d_sh_cz=nullptr, *d_sh_md2=nullptr;
  double *d_ex=nullptr, *d_cc=nullptr, *d_pmd2=nullptr, *d_bfnrm=nullptr;
  int    *d_cartx=nullptr, *d_carty=nullptr, *d_cartz=nullptr;
  int    maxang=0, maxcart=0;
  // working buffers
  double *d_phi=nullptr;     // [nbf*npts]   AO values
  double *d_gphi=nullptr;    // [nbf*npts*3] AO gradients (x,y,z planes)
  double *d_Dfold=nullptr;   // [nbf*nbf]
  double *d_X=nullptr;       // [nbf*npts]   Dfold@phi (== native moVA)
  double *d_rho=nullptr;     // [npts]       total density
  double *d_drho=nullptr;    // [npts*3]     ALPHA density gradient (x,y,z planes)
  double *d_vrho=nullptr;    // [npts]       w*d1dr(ra)
  double *d_vsa=nullptr;     // [npts]       w*d1ds(ga)
  double *d_vsc=nullptr;     // [npts]       w*d1ds(gc)
  double *d_eps=nullptr;     // [npts]       eps_xc per electron (total)
  double *d_Z=nullptr;       // [nbf*npts]
  double *d_Vfull=nullptr;   // [nbf*nbf]
  // ---- UKS (open-shell, nfocks==2) extra channels (alpha reuses RKS bufs:
  //      d_vrho=vra_w, d_vsa=vsaa_w, d_vsc=vsab_w; d_rho/d_drho = alpha) ----
  double *d_Dfold_b=nullptr; // [nbf*nbf]   beta folded density
  double *d_Xb=nullptr;      // [nbf*npts]  Dfold_b @ phi
  double *d_rho_b=nullptr;   // [npts]      beta density
  double *d_drho_b=nullptr;  // [npts*3]    beta density gradient
  double *d_vrb=nullptr;     // [npts]      w*dE/drho_b
  double *d_vsbb=nullptr;    // [npts]      w*dE/dsigma_bb
  double *d_Zb=nullptr;      // [nbf*npts]
  double *d_Vfull_b=nullptr; // [nbf*nbf]
  double *d_excd=nullptr;    // [npts]      XC energy density (UKS reduction)
  double *d_eexc=nullptr, *d_totele=nullptr;
  cublasHandle_t cub=nullptr;
  // ---- active-AO SPARSE XC (OQP_OWNXC_SPARSE): per grid block only the shells
  //      that reach it are collocated/contracted, turning the two dense
  //      nbf^2*bp GEMMs into nact^2*bp (O(N^3)->O(N^2)). ----
  bool sparse=false;
  int  sblk=0, nblock=0, nact_max=0;
  std::vector<int> h_bsh_off, h_ao_off;   // host CSR offsets (per-block sizes)
  int* d_bsh_gid=nullptr;    // [sum nsh_b]  global shell index, per block
  int* d_bsh_row=nullptr;    // [sum nsh_b]  compact starting AO row of that shell
  int* d_ao_idx=nullptr;     // [sum nact]   global AO index, per block (gather/scatter map)
  double *d_phic=nullptr, *d_gphic=nullptr; // [nact_max x sblk] (+*3)
  double *d_Xc=nullptr, *d_Zc=nullptr;      // [nact_max x sblk]
  double *d_Dc=nullptr, *d_Vc=nullptr;      // [nact_max x nact_max]
  double *d_Xcb=nullptr, *d_Zcb=nullptr;    // UKS beta compact [nact_max x sblk]
  double *d_Dcb=nullptr, *d_Vcb=nullptr;    // UKS beta compact [nact_max x nact_max]
  // multi-GPU XC (OQP_MULTI_GPU): one Ctx per device, each owning the grid
  // slice [npts*rank/nranks, npts*(rank+1)/nranks). XC is a sum over points,
  // so the split is exact; partials are combined on the HOST (no peer copies).
  int dev=0, rank=0, nranks=1;
} g, g_xc1;

// =============================================================================
// collocation: AO values phi AND gradients gphi (3 planes), matching
// basis_tools.F90 compAOv (values) + compAOvg (gradients).
//   gphi layout: gphi[u + nbf*p + nbf*npts*j], j=0(x),1(y),2(z).
// vexp1 = sum cc*exp(-ex r2);  vexp2 = sum 2*ex*cc*exp(-ex r2).
// general l: aogx = (-dr[ix+1]*vexp2 + vexp1*ix*dr[ix-1]) * dr[iy]*dr[iz], etc.
// =============================================================================
__global__ void k_collocate_g(int npts, int nshell, int nbf, int gga,
    const double* __restrict__ gx, const double* __restrict__ gy, const double* __restrict__ gz,
    const int* sh_am, const int* sh_g0, const int* sh_nc, const int* sh_ao,
    const double* sh_cx, const double* sh_cy, const double* sh_cz, const double* sh_md2,
    const double* ex, const double* cc, const double* pmd2,
    const int* cartx, const int* carty, const int* cartz, int maxcart,
    double* __restrict__ phi, double* __restrict__ gphi)
{
  int p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= npts) return;
  long npL = npts, nbL = nbf;
  double px = gx[p], py = gy[p], pz = gz[p];
  for (int ish = 0; ish < nshell; ish++) {
    int am  = sh_am[ish];
    int g0  = sh_g0[ish];
    int nc  = sh_nc[ish];
    int ao0 = sh_ao[ish];
    double dx = px - sh_cx[ish];
    double dy = py - sh_cy[ish];
    double dz = pz - sh_cz[ish];
    double r2 = dx*dx + dy*dy + dz*dz;
    int ncart = (am+1)*(am+2)/2;
    if (r2 > sh_md2[ish]) {
      for (int ic = 0; ic < ncart; ic++) {
        long o = (ao0-1+ic) + nbL*p;
        phi[o] = 0.0;
        if (gga) { gphi[o]=0.0; gphi[o+nbL*npL]=0.0; gphi[o+2*nbL*npL]=0.0; }
      }
      continue;
    }
    double vexp1 = 0.0, vexp2 = 0.0;
    for (int k = 0; k < nc; k++) {
      int ig = g0-1 + k;
      if (r2 > pmd2[ig]) continue;
      double v = exp(-ex[ig]*r2) * cc[ig];
      vexp1 += v;
      vexp2 += 2.0*ex[ig]*v;
    }
    if (am == 0) {
      long o = (ao0-1) + nbL*p;
      phi[o] = vexp1;
      if (gga) { gphi[o]=-vexp2*dx; gphi[o+nbL*npL]=-vexp2*dy; gphi[o+2*nbL*npL]=-vexp2*dz; }
    } else if (am == 1) {
      long o0=(ao0-1)+nbL*p, o1=o0+1, o2=o0+2;
      phi[o0]=vexp1*dx; phi[o1]=vexp1*dy; phi[o2]=vexp1*dz;
      if (gga) {
        gphi[o0]            = vexp1 - vexp2*dx*dx;
        gphi[o1]            =       - vexp2*dx*dy;
        gphi[o2]            =       - vexp2*dx*dz;
        gphi[o0+nbL*npL]    =       - vexp2*dx*dy;
        gphi[o1+nbL*npL]    = vexp1 - vexp2*dy*dy;
        gphi[o2+nbL*npL]    =       - vexp2*dy*dz;
        gphi[o0+2*nbL*npL]  =       - vexp2*dx*dz;
        gphi[o1+2*nbL*npL]  =       - vexp2*dy*dz;
        gphi[o2+2*nbL*npL]  = vexp1 - vexp2*dz*dz;
      }
    } else {
      // dr powers with dr[-1]=0, dr[0]=1, dr[i]=comp^i (index shifted by 1 -> use [i+1])
      // store drx[i] for i in 0..am+1 mapped to native dr1(i,:) with i=-1..am+1
      double drx[12], dry[12], drz[12];   // index 0 == native dr1(-1)=0
      drx[0]=dry[0]=drz[0]=0.0;           // dr1(-1)
      drx[1]=dry[1]=drz[1]=1.0;           // dr1(0)
      drx[2]=dx; dry[2]=dy; drz[2]=dz;    // dr1(1)
      for (int i = 2; i <= am+1; i++) {   // dr1(i) = dr1(i-1)*dr1(1)
        drx[i+1]=drx[i]*dx; dry[i+1]=dry[i]*dy; drz[i+1]=drz[i]*dz;
      }
      for (int ic = 0; ic < ncart; ic++) {
        int ix = cartx[ic + maxcart*am];
        int iy = carty[ic + maxcart*am];
        int iz = cartz[ic + maxcart*am];
        // native dr1(k) maps to local index (k+1): dr1(ix)->[ix+1], dr1(ix-1)->[ix], dr1(ix+1)->[ix+2]
        double X = drx[ix+1], Y = dry[iy+1], Z = drz[iz+1];
        long o = (ao0-1+ic) + nbL*p;
        phi[o] = vexp1 * X*Y*Z;
        if (gga) {
          double xm = ix*drx[ix];          // ix*dr1(ix-1)
          double ym = iy*dry[iy];
          double zm = iz*drz[iz];
          double xp = drx[ix+2]*vexp2;      // dr1(ix+1)*vexp2
          double yp = dry[iy+2]*vexp2;
          double zp = drz[iz+2]*vexp2;
          gphi[o]           = (-xp + vexp1*xm)*Y*Z;
          gphi[o+nbL*npL]   = (-yp + vexp1*ym)*X*Z;
          gphi[o+2*nbL*npL] = (-zp + vexp1*zm)*X*Y;
        }
      }
    }
  }
}

// rho_tot(p) = sum_u phi[u,p]*X[u,p]  (X = Dfold @ phi). closed shell rho_a=0.5*rho_tot.
// One WARP per grid point, reducing over the AO index u in stride-32 chunks so
// the phi/X reads are COALESCED (lanes 0..31 hit consecutive u = u0..u0+31 for a
// fixed p). Was: one thread per p looping u -> stride-nbf uncoalesced.
__global__ void k_rho(int npts, int nbf, const double* phi, const double* X, double* rho) {
  const int p = (int)(((long)blockIdx.x*blockDim.x + threadIdx.x) >> 5);
  const int lane = threadIdx.x & 31;
  if (p >= npts) return;
  const long nbL=nbf;
  double s = 0.0;
  for (int u = lane; u < nbf; u += 32) s += phi[u + nbL*p] * X[u + nbL*p];
  for (int o = 16; o > 0; o >>= 1) s += __shfl_down_sync(0xffffffffu, s, o);
  if (lane == 0) rho[p] = s;
}

// ALPHA density gradient (native compDRhoAO closed-shell, NO factor 2):
//   drhoa(j) = sum_u gphi[u,j]*X[u].   stored drho[p + npts*j], j=0,1,2.
__global__ void k_drho(int npts, int nbf, const double* gphi, long gstr,
                       const double* X, double* drho) {
  // gstr = gphi COMPONENT stride in elements (nbf*npts_of_the_gphi_buffer): when the
  // full grid is cached it is nbf*npts_total while this call sees only a bp-wide slice.
  const int p = (int)(((long)blockIdx.x*blockDim.x + threadIdx.x) >> 5);
  const int lane = threadIdx.x & 31;
  if (p >= npts) return;
  const long nbL=nbf, npL=npts;
  double sx=0.0, sy=0.0, sz=0.0;
  for (int u = lane; u < nbf; u += 32) {
    double xu = X[u + nbL*p];
    sx += gphi[u + nbL*p]         * xu;
    sy += gphi[u + nbL*p + gstr]  * xu;
    sz += gphi[u + nbL*p + 2*gstr]* xu;
  }
  for (int o = 16; o > 0; o >>= 1) {
    sx += __shfl_down_sync(0xffffffffu, sx, o);
    sy += __shfl_down_sync(0xffffffffu, sy, o);
    sz += __shfl_down_sync(0xffffffffu, sz, o);
  }
  if (lane == 0) { drho[p] = sx; drho[p + npL] = sy; drho[p + 2*npL] = sz; }
}

// ============================ FUNCTIONALS ====================================
// All hand-coded closed forms. Spin-polarized libxc convention.

// ---- VWN5 paramagnetic eps + dEps/drs (rung-1, retained for SVWN) ----------
__device__ __forceinline__ double vwn_eps_and_v(double rho_tot, double &vc,
        double A1, double b1, double c1, double x01) {
  const double pi = 3.14159265358979323846;
  double Q  = sqrt(4.0*c1 - b1*b1);
  double rs = cbrt(3.0/(4.0*pi*rho_tot));
  double x  = sqrt(rs);
  double X  = rs + b1*x + c1;
  double X0 = x01*x01 + b1*x01 + c1;
  double atanarg = Q/(2.0*x + b1);
  double f1 = 2.0*b1/Q;
  double f2 = b1*x01/X0;
  double f3 = 2.0*(2.0*x01 + b1)/Q;
  double eps = A1*( log(rs/X) + (f1 - f2*f3)*atan(atanarg)
                  - f2*log((x - x01)*(x - x01)/X) );
  double dX_drs = 1.0 + b1/(2.0*x);
  double d_log_rsX = (1.0/rs) - dX_drs/X;
  double d_atan = -(Q/((2.0*x+b1)*(2.0*x+b1) + Q*Q)) * (1.0/x);
  double d_term3 = -f2*( 2.0*(1.0/(x - x01))*(1.0/(2.0*x)) - dX_drs/X );
  double deps_drs = A1*( d_log_rsX + (f1 - f2*f3)*d_atan + d_term3 );
  vc = eps - (rs/3.0)*deps_drs;   // d(rho*eps)/drho = eps - (rs/3) deps/drs
  return eps;
}

// ---- spin-polarized VWN (parameterized: works for VWN5 and VWN_RPA) --------
// vwn_aux2 = the VWN auxiliary G(rs; A,b,c,x0) returning eps and d eps/d rs.
__device__ __forceinline__ void vwn_aux2(double rs, double A, double b, double c, double x0,
                                         double &eps, double &deps_drs) {
  double Q = sqrt(4.0*c - b*b);
  double x = sqrt(rs);
  double X = rs + b*x + c;
  double X0 = x0*x0 + b*x0 + c;
  double f1 = 2.0*b/Q, f2 = b*x0/X0, f3 = 2.0*(2.0*x0 + b)/Q;
  eps = A*( log(rs/X) + (f1 - f2*f3)*atan(Q/(2.0*x + b)) - f2*log((x - x0)*(x - x0)/X) );
  double dX = 1.0 + b/(2.0*x);
  double d_atan = -(Q/((2.0*x+b)*(2.0*x+b) + Q*Q)) * (1.0/x);
  deps_drs = A*( (1.0/rs - dX/X) + (f1 - f2*f3)*d_atan
               - f2*( 2.0*(1.0/(x - x0))*(1.0/(2.0*x)) - dX/X ) );
}
// eps_c(rs,zeta) with the VWN5 spin interpolation, given the three fits
// (paramagnetic / ferromagnetic / spin-stiffness). Returns eps_c + v_c,a + v_c,b.
__device__ __forceinline__ void vwn_spin_param(double ra, double rb,
        double Ap,double bp,double cp,double x0p,
        double Af,double bf,double cf,double x0f,
        double Aa,double ba,double ca,double x0a,
        double &eps_c, double &vca, double &vcb) {
  const double pi = 3.14159265358979323846;
  double rho = ra + rb;
  if (rho < 1e-30) { eps_c=0.0; vca=0.0; vcb=0.0; return; }
  double zeta = (ra - rb)/rho;
  if (zeta >  1.0) zeta =  1.0;
  if (zeta < -1.0) zeta = -1.0;
  double rs = cbrt(3.0/(4.0*pi*rho));
  double epsP,dP, epsF,dF, alc,dal;
  vwn_aux2(rs, Ap,bp,cp,x0p, epsP, dP);
  vwn_aux2(rs, Af,bf,cf,x0f, epsF, dF);
  vwn_aux2(rs, Aa,ba,ca,x0a, alc, dal);
  double two13 = cbrt(2.0), fden = 2.0*two13 - 2.0;
  double opz = 1.0+zeta, omz = 1.0-zeta, o13 = cbrt(opz), m13 = cbrt(omz);
  double fz  = (opz*o13 + omz*m13 - 2.0)/fden;
  double fpz = (4.0/3.0)*(o13 - m13)/fden;
  const double fpp0 = 4.0/(9.0*(two13 - 1.0));
  double z3 = zeta*zeta*zeta, z4 = z3*zeta, beta = fz/fpp0, dbeta = fpz/fpp0;
  eps_c = epsP + alc*beta*(1.0 - z4) + (epsF - epsP)*fz*z4;
  double deps_drs = dP + dal*beta*(1.0 - z4) + (dF - dP)*fz*z4;
  double deps_dz  = alc*( dbeta*(1.0 - z4) + beta*(-4.0*z3) )
                  + (epsF - epsP)*( fpz*z4 + fz*4.0*z3 );
  double drs_drho = -rs/(3.0*rho);
  double dz_dra =  (1.0 - zeta)/rho;
  double dz_drb = -(1.0 + zeta)/rho;
  vca = eps_c + rho*(deps_drs*drs_drho + deps_dz*dz_dra);
  vcb = eps_c + rho*(deps_drs*drs_drho + deps_dz*dz_drb);
}

// VWN "functional I" spin interpolation (used by VWN1-RPA in B3LYP/B3LYPV1R):
//   eps_c(rs,zeta) = eps_P + (eps_F - eps_P)*f(zeta)   -- NO spin-stiffness alpha_c.
__device__ __forceinline__ void vwn1_spin(double ra, double rb,
        double Ap,double bp,double cp,double x0p,
        double Af,double bf,double cf,double x0f,
        double &eps_c, double &vca, double &vcb) {
  const double pi = 3.14159265358979323846;
  double rho = ra + rb;
  if (rho < 1e-30) { eps_c=0.0; vca=0.0; vcb=0.0; return; }
  double zeta = (ra - rb)/rho;
  if (zeta >  1.0) zeta =  1.0;
  if (zeta < -1.0) zeta = -1.0;
  double rs = cbrt(3.0/(4.0*pi*rho));
  double epsP,dP, epsF,dF;
  vwn_aux2(rs, Ap,bp,cp,x0p, epsP, dP);
  vwn_aux2(rs, Af,bf,cf,x0f, epsF, dF);
  double two13 = cbrt(2.0), fden = 2.0*two13 - 2.0;
  double opz = 1.0+zeta, omz = 1.0-zeta, o13 = cbrt(opz), m13 = cbrt(omz);
  double fz  = (opz*o13 + omz*m13 - 2.0)/fden;
  double fpz = (4.0/3.0)*(o13 - m13)/fden;
  eps_c = epsP + (epsF - epsP)*fz;
  double deps_drs = dP + (dF - dP)*fz;
  double deps_dz  = (epsF - epsP)*fpz;
  double drs_drho = -rs/(3.0*rho);
  double dz_dra =  (1.0 - zeta)/rho;
  double dz_drb = -(1.0 + zeta)/rho;
  vca = eps_c + rho*(deps_drs*drs_drho + deps_dz*dz_dra);
  vcb = eps_c + rho*(deps_drs*drs_drho + deps_dz*dz_drb);
}

// ---- B88 exchange (spin-polarized), libxc gga_x_b88 ------------------------
// Per spin: E_x[rho_s, sigma_ss] with the spin-scaling Ex = sum_s ex(2 rho_s? )...
// We use the standard form on a single spin density n and its |grad n|^2 = s2:
//   Slater per-spin LDA exchange:  e_lda = -Cx2 * n^{4/3},  Cx2=(3/4)(6/pi)^{1/3}/?
// libxc uses: f(x) = 1 + beta/Cx_b88 * x^2/(1+6 beta x asinh(x)),  x=|grad n|/n^{4/3}
// with Ex = -Cx * n^{4/3} * f(x), Cx = (3/2)(3/(4pi))^{1/3}  (per-spin LDA-X const).
// Returns: e_dens = Ex (energy density for this spin channel), and
//   dn = dEx/dn , ds2 = dEx/dsigma_ss.
__device__ __forceinline__ void b88_x_spin(double n, double s2,
        double &e_dens, double &dn, double &ds2) {
  // per-spin LDA exchange constant Cx = (3/2)*(3/(4 pi))^(1/3)
  const double pi = 3.14159265358979323846;
  const double Cx = 1.5 * cbrt(3.0/(4.0*pi));   // = -coeff of n^{4/3} (sign added below)
  const double beta = 0.0042;
  if (n < 1e-30) { e_dens=0.0; dn=0.0; ds2=0.0; return; }
  double n13 = cbrt(n);
  double n43 = n13*n;                 // n^{4/3}
  double gr  = sqrt(s2);              // |grad n|
  double x   = gr / n43;              // reduced gradient
  double asx = asinh(x);              // asinh(x)=log(x+sqrt(1+x^2))
  double denom = 1.0 + 6.0*beta*x*asx;
  double Fx = 1.0 + (beta/Cx) * (x*x)/denom;   // enhancement over LDA (Cx defined +)
  // Ex = -Cx * n^{4/3} * Fx
  e_dens = -Cx * n43 * Fx;
  // derivatives. Let g(x) = (beta/Cx) x^2 / (1+6 beta x asinh x).
  // dEx/dn = -Cx*(4/3)n^{1/3}*Fx  - Cx*n^{4/3} * g'(x) * dx/dn
  // dEx/dsigma = -Cx*n^{4/3} * g'(x) * dx/dsigma
  // x = gr * n^{-4/3} => dx/dn = -(4/3) x / n ; dx/dgr = 1/n^{4/3} ;
  //   dx/dsigma = dx/dgr * dgr/dsigma = (1/n^{4/3}) * 1/(2 gr) = x/(2 sigma)
  double sqrt1x2 = sqrt(1.0 + x*x);
  double dasx = 1.0/sqrt1x2;                 // d/dx asinh x
  double dden = 6.0*beta*(asx + x*dasx);     // d denom / dx
  double num  = (beta/Cx)*x*x;
  double dg   = (beta/Cx)*(2.0*x*denom - x*x*dden)/(denom*denom); // g'(x)
  double dFx_dx = dg;                        // Fx = 1 + g(x)
  double dx_dn  = -(4.0/3.0)*x/n;
  double dEx_dn_lda = -Cx*(4.0/3.0)*n13*Fx;
  dn  = dEx_dn_lda + (-Cx*n43)*dFx_dx*dx_dn;
  double dx_dsig = (s2>1e-50) ? x/(2.0*s2) : 0.0;
  ds2 = (-Cx*n43)*dFx_dx*dx_dsig;
}

// ---- LYP correlation (spin-polarized closed-shell-capable), libxc gga_c_lyp -
// Standard Lee-Yang-Parr (Miehlich-Savin-Stoll-Preuss form).
// Inputs: ra, rb (spin densities), sigma_aa, sigma_ab, sigma_bb.
// Outputs: E_dens (correlation energy density, total), and partials
//   dra=dE/dra, drb=dE/drb, dsaa=dE/dsigma_aa, dsab=dE/dsigma_ab, dsbb=dE/dsigma_bb.
// Constants:
__device__ __forceinline__ void lyp_c(double ra, double rb,
        double saa, double sab, double sbb,
        double &E, double &dra, double &drb,
        double &dsaa, double &dsab, double &dsbb) {
  const double A = 0.04918, B = 0.132, C = 0.2533, Dd = 0.349;
  const double pi = 3.14159265358979323846;
  double rho = ra + rb;
  E=0; dra=0; drb=0; dsaa=0; dsab=0; dsbb=0;
  if (rho < 1e-30) return;
  double r13 = cbrt(rho);
  double rm13 = 1.0/r13;
  double den = 1.0 + Dd*rm13;          // 1 + d rho^{-1/3}
  double omega = exp(-C*rm13)/den * pow(rho, -11.0/3.0);
  double delta = C*rm13 + Dd*rm13/den; // delta = c rho^{-1/3} + d rho^{-1/3}/(1+d rho^{-1/3})
  double CF = 0.3 * pow(3.0*pi*pi, 2.0/3.0);   // (3/10)(3 pi^2)^{2/3}
  double rara = ra*rb;
  // sigma totals
  double gaa = saa, gbb = sbb, gtot = saa + 2.0*sab + sbb;  // |grad rho|^2
  // The Miehlich et al. closed expression for E_c:
  //  E = -A*4/den*ra*rb/rho
  //      -A*B*omega*[ ra*rb*( 2^{11/3} CF (ra^{8/3}+rb^{8/3})
  //          + (47/18 - 7/18 delta) gtot
  //          - (5/2 - 1/18 delta)(gaa+gbb)
  //          - (delta-11)/9 *(ra/rho*gaa + rb/rho*gbb) )
  //        - 2/3 rho^2 gtot + (2/3 rho^2 - ra^2) gbb + (2/3 rho^2 - rb^2) gaa ]
  double t1 = -A*4.0/den * rara/rho;
  double f83a = pow(ra, 8.0/3.0);
  double f83b = pow(rb, 8.0/3.0);
  double c113 = pow(2.0, 11.0/3.0);
  double term_kin = c113*CF*(f83a + f83b);
  double k1 = (47.0/18.0 - 7.0/18.0*delta);
  double k2 = (2.5 - 1.0/18.0*delta);
  double k3 = (delta - 11.0)/9.0;
  double inbr = rara*( term_kin
                     + k1*gtot
                     - k2*(gaa+gbb)
                     - k3*(ra/rho*gaa + rb/rho*gbb) )
              - 2.0/3.0*rho*rho*gtot
              + (2.0/3.0*rho*rho - ra*ra)*gbb
              + (2.0/3.0*rho*rho - rb*rb)*gaa;
  double t2 = -A*B*omega*inbr;
  E = t1 + t2;

  // ---- numerical-safe analytic derivatives via finite differences on the
  //      closed form would break "no extra error"; instead derive analytically.
  // We compute partials w.r.t the 5 inputs analytically.
  // Precompute derivatives of scalar helpers wrt rho (then chain to ra,rb).
  // d(rm13)/drho = -1/3 rho^{-4/3}
  double drm13 = -1.0/3.0*pow(rho, -4.0/3.0);
  // den = 1 + Dd*rm13 ; dden/drho = Dd*drm13
  double dden = Dd*drm13;
  // delta = C*rm13 + Dd*rm13/den
  //   ddelta/drho = C*drm13 + Dd*(drm13*den - rm13*dden)/den^2
  double ddelta = C*drm13 + Dd*(drm13*den - rm13*dden)/(den*den);
  // omega = exp(-C rm13) * rho^{-11/3} / den
  //   d ln omega/drho = -C*drm13 - 11/3/rho - dden/den
  double dln_omega = -C*drm13 - (11.0/3.0)/rho - dden/den;
  double domega = omega*dln_omega;

  // ---- sigma partials (linear in sigmas) -------------------------------------
  // gtot = gaa + 2 gab + gbb ; so d/dsaa contributions:
  //   from inbr: rara*(k1 - k2 - k3*ra/rho) [gaa coeff via gtot,(gaa+gbb),and ra term]
  //              wait separate gaa, gbb, gtot, sab carefully.
  // Collect coefficient of each sigma in inbr:
  //   coeff_gaa = rara*(k1 - k2 - k3*ra/rho) - 2/3 rho^2 + (2/3 rho^2 - rb*rb)
  //   coeff_gbb = rara*(k1 - k2 - k3*rb/rho) - 2/3 rho^2 + (2/3 rho^2 - ra*ra)
  //   coeff_gab = rara*(2*k1)               - 2*2/3 rho^2     (gtot has 2*sab)
  double coeff_gaa = rara*(k1 - k2 - k3*ra/rho) - 2.0/3.0*rho*rho + (2.0/3.0*rho*rho - rb*rb);
  double coeff_gbb = rara*(k1 - k2 - k3*rb/rho) - 2.0/3.0*rho*rho + (2.0/3.0*rho*rho - ra*ra);
  double coeff_gab = rara*(2.0*k1) - 2.0*(2.0/3.0*rho*rho);
  // E partials wrt sigma = -A*B*omega * coeff
  dsaa = -A*B*omega*coeff_gaa;
  dsbb = -A*B*omega*coeff_gbb;
  dsab = -A*B*omega*coeff_gab;

  // ---- density partials ------------------------------------------------------
  // t1 = -4A * ra*rb/(rho*den)
  //   d t1/dra = -4A * [ (rb*rho*den - ra*rb*(den + rho*dden)) / (rho*den)^2 ]
  // Use product/quotient. Let P = ra*rb, Q = rho*den.
  double P = rara, Qd = rho*den;
  double dP_dra = rb, dP_drb = ra;
  double dQ_dra = den + rho*dden;   // drho/dra=1
  double dQ_drb = den + rho*dden;
  double dt1_dra = -4.0*A*(dP_dra*Qd - P*dQ_dra)/(Qd*Qd);
  double dt1_drb = -4.0*A*(dP_drb*Qd - P*dQ_drb)/(Qd*Qd);

  // t2 = -A*B*omega*inbr.  d t2 = -A*B*(domega*inbr + omega*dinbr)
  // need d inbr/dra and d inbr/drb (sigmas held fixed).
  // inbr = rara*S + T, where
  //   S = term_kin + k1*gtot - k2*(gaa+gbb) - k3*(ra/rho*gaa + rb/rho*gbb)
  //   T = -2/3 rho^2 gtot + (2/3 rho^2 - ra^2) gbb + (2/3 rho^2 - rb^2) gaa
  double S = term_kin + k1*gtot - k2*(gaa+gbb) - k3*(ra/rho*gaa + rb/rho*gbb);
  // d rara/dra = rb ; d rara/drb = ra
  // dS/dra: d term_kin/dra = c113*CF*(8/3 ra^{5/3})
  double dtk_dra = c113*CF*(8.0/3.0)*pow(ra,5.0/3.0);
  double dtk_drb = c113*CF*(8.0/3.0)*pow(rb,5.0/3.0);
  // k1,k2 depend on delta(rho): dk1=-7/18 ddelta ; dk2=-1/18 ddelta ; dk3=ddelta/9
  double dk1 = -7.0/18.0*ddelta;
  double dk2 = -1.0/18.0*ddelta;
  double dk3 = (1.0/9.0)*ddelta;
  // ra/rho derivative wrt ra = (rho-ra)/rho^2 = rb/rho^2 ; wrt rb = -ra/rho^2
  double draOrho_dra = rb/(rho*rho);
  double draOrho_drb = -ra/(rho*rho);
  double drbOrho_dra = -rb/(rho*rho);
  double drbOrho_drb = ra/(rho*rho);
  double raOrho = ra/rho, rbOrho = rb/rho;
  double dS_dra = dtk_dra + dk1*gtot - dk2*(gaa+gbb)
                - ( dk3*(raOrho*gaa + rbOrho*gbb)
                  + k3*(draOrho_dra*gaa + drbOrho_dra*gbb) );
  double dS_drb = dtk_drb + dk1*gtot - dk2*(gaa+gbb)
                - ( dk3*(raOrho*gaa + rbOrho*gbb)
                  + k3*(draOrho_drb*gaa + drbOrho_drb*gbb) );
  // T partials: T = -2/3 rho^2 gtot + (2/3 rho^2 - ra^2) gbb + (2/3 rho^2 - rb^2) gaa
  // drho/dra=1
  double dT_dra = -2.0/3.0*2.0*rho*gtot + (2.0/3.0*2.0*rho - 2.0*ra)*gbb + (2.0/3.0*2.0*rho)*gaa;
  double dT_drb = -2.0/3.0*2.0*rho*gtot + (2.0/3.0*2.0*rho)*gbb + (2.0/3.0*2.0*rho - 2.0*rb)*gaa;
  double dinbr_dra = rb*S + rara*dS_dra + dT_dra;
  double dinbr_drb = ra*S + rara*dS_drb + dT_drb;
  double dt2_dra = -A*B*(domega*inbr + omega*dinbr_dra);
  double dt2_drb = -A*B*(domega*inbr + omega*dinbr_drb);
  dra = dt1_dra + dt2_dra;
  drb = dt1_drb + dt2_drb;
}

// ---- main per-point functional kernel --------------------------------------
// Produces eps_out (per-electron eps_xc, for Exc=sum w rho eps) and the
// weight-folded potentials vrho=w*d1dr(ra), vsa=w*d1ds(ga), vsc=w*d1ds(gc).
__global__ void k_func_gga(int npts, int func,
        const double* rho_tot, const double* drho, const double* w,
        double* eps_out, double* vrho, double* vsa, double* vsc) {
  int p = blockIdx.x*blockDim.x + threadIdx.x;
  if (p >= npts) return;
  long npL=npts;
  const double pi = 3.14159265358979323846;
  double rt = rho_tot[p];
  double eps_total = 0.0;     // per-electron eps_xc of TOTAL density
  double d_ra = 0.0;          // dE/drho_a
  double d_ga = 0.0;          // dE/dsigma_aa
  double d_gc = 0.0;          // dE/dsigma_ab
  if (rt > 1e-30) {
    double ra = 0.5*rt, rb = ra;
    // closed-shell alpha gradient components
    double dxa=0.0, dya=0.0, dza=0.0;
    if (func_is_gga(func)) { dxa=drho[p]; dya=drho[p+npL]; dza=drho[p+2*npL]; }
    double saa = dxa*dxa + dya*dya + dza*dza;
    double sab = saa;   // drhoa.drhob, all equal closed shell
    double sbb = saa;

    double Etot = 0.0;  // energy density (per volume), total

    // ----- exchange -----
    if (func == FUNC_SLATER) {
      // Slater spin-unpol total: handled via per-spin sum below for consistency.
      // Ex_total = sum_s -Cx*(rho_s)^{4/3}*2^{1/3}? Use the well-known total form:
      const double Cx = 0.75 * pow(3.0/pi, 1.0/3.0);
      double r13 = cbrt(rt);
      double eps_x = -Cx*r13;
      eps_total += eps_x;
      d_ra += (4.0/3.0)*eps_x;     // dEx/drho_a (closed shell == dEx/drho_tot)
    } else {
      // B88 exchange per spin (sum a,b). For closed shell symmetric.
      double cxS=0.0, cxB=1.0;  // exchange-coefficient for the mix
      if (func == FUNC_B3LYP)  { cxS=0.08; cxB=0.72; }
      else if (func==FUNC_BHHLYP){ cxS=0.0; cxB=0.5; }
      else /*BLYP*/            { cxS=0.0; cxB=1.0; }
      // Slater piece (B3LYP carries 0.08 LDA exchange)
      if (cxS != 0.0) {
        const double Cx = 0.75 * pow(3.0/pi, 1.0/3.0);
        double r13 = cbrt(rt);
        double eps_x = -Cx*r13;
        Etot += cxS * rt * eps_x;
        d_ra += cxS * (4.0/3.0)*eps_x;
      }
      // B88 piece per spin
      double ea, dna, dsa;  b88_x_spin(ra, saa, ea, dna, dsa);
      double eb, dnb, dsb;  b88_x_spin(rb, sbb, eb, dnb, dsb);
      Etot += cxB*(ea + eb);
      d_ra += cxB*dna;                 // dEx/drho_a (alpha channel)
      d_ga += cxB*dsa;                 // dEx/dsigma_aa
      // sigma_ab does not enter B88 -> d_gc += 0
    }

    // ----- correlation -----
    if (func == FUNC_SVWN) {
      double vc; double eps_c = vwn_eps_and_v(rt, vc, 0.0310907,3.72744,12.9352,-0.10498);
      eps_total += eps_c;
      d_ra += vc;
    } else if (func == FUNC_B3LYP) {
      // 0.19 VWN_RPA + 0.81 LYP
      double vc; double eps_c = vwn_eps_and_v(rt, vc, 0.0310907,13.0720,42.7198,-0.409286);
      Etot += 0.19 * rt * eps_c;
      d_ra += 0.19 * vc;
      double E,dra2,drb2,dsaa,dsab,dsbb;
      lyp_c(ra,rb,saa,sab,sbb, E,dra2,drb2,dsaa,dsab,dsbb);
      Etot += 0.81*E;
      d_ra += 0.81*dra2;
      d_ga += 0.81*dsaa;
      d_gc += 0.81*dsab;
    } else if (func == FUNC_BLYP || func == FUNC_BHHLYP) {
      double E,dra2,drb2,dsaa,dsab,dsbb;
      lyp_c(ra,rb,saa,sab,sbb, E,dra2,drb2,dsaa,dsab,dsbb);
      Etot += E;
      d_ra += dra2;
      d_ga += dsaa;
      d_gc += dsab;
    }

    // fold the energy-density pieces (Etot) into per-electron eps and add the
    // explicit eps_total (Slater/VWN per-electron) contributions.
    eps_total += Etot / rt;
  }
  eps_out[p] = eps_total;     // Exc = sum w*rho*eps_out
  vrho[p] = w[p]*d_ra;        // weight-folded d1dr(ra)
  vsa[p]  = w[p]*d_ga;        // weight-folded d1ds(ga)
  vsc[p]  = w[p]*d_gc;        // weight-folded d1ds(gc)
}

// ---- UKS per-point functional: true ra,rb + per-channel gradients ----------
// Emits XC energy density excd and the weight-folded potentials:
//   vra_w=w*dE/drho_a, vrb_w=w*dE/drho_b, vsaa_w=w*dE/dsigma_aa,
//   vsbb_w=w*dE/dsigma_bb, vsab_w=w*dE/dsigma_ab. B88/LYP already spin-polarized;
//   Slater per-spin (2^1/3 form); VWN5 / VWN_RPA via vwn_spin_param.
__global__ void k_func_gga_uks(int npts, int func,
        const double* rho_a, const double* rho_b,
        const double* drho_a, const double* drho_b, const double* w,
        double* excd, double* vra_w, double* vrb_w,
        double* vsaa_w, double* vsbb_w, double* vsab_w) {
  int p = blockIdx.x*blockDim.x + threadIdx.x;
  if (p >= npts) return;
  long npL = npts;
  const double pi = 3.14159265358979323846;
  double ra = rho_a[p], rb = rho_b[p];
  double E = 0.0, d_ra = 0.0, d_rb = 0.0, d_saa = 0.0, d_sbb = 0.0, d_sab = 0.0;
  bool gga = func_is_gga(func);
  double ax=0,ay=0,az=0, bx=0,by=0,bz=0;
  if (gga) { ax=drho_a[p]; ay=drho_a[p+npL]; az=drho_a[p+2*npL];
             bx=drho_b[p]; by=drho_b[p+npL]; bz=drho_b[p+2*npL]; }
  double saa = ax*ax+ay*ay+az*az;
  double sbb = bx*bx+by*by+bz*bz;
  double sab = ax*bx+ay*by+az*bz;
  // ----- exchange -----
  if (func == FUNC_SLATER) {
    const double Cx = 0.75*pow(3.0/pi,1.0/3.0); const double two13 = cbrt(2.0);
    if (ra>1e-30){ E += -Cx*two13*pow(ra,4.0/3.0); d_ra += (4.0/3.0)*(-Cx*two13)*cbrt(ra); }
    if (rb>1e-30){ E += -Cx*two13*pow(rb,4.0/3.0); d_rb += (4.0/3.0)*(-Cx*two13)*cbrt(rb); }
  } else {
    double cxS=0.0, cxB=1.0;
    if (func==FUNC_B3LYP){ cxS=0.08; cxB=0.72; }
    else if (func==FUNC_BHHLYP){ cxS=0.0; cxB=0.5; }
    else { cxS=0.0; cxB=1.0; }
    if (cxS != 0.0) {
      const double Cx = 0.75*pow(3.0/pi,1.0/3.0); const double two13 = cbrt(2.0);
      if (ra>1e-30){ E += cxS*(-Cx*two13*pow(ra,4.0/3.0)); d_ra += cxS*(4.0/3.0)*(-Cx*two13)*cbrt(ra); }
      if (rb>1e-30){ E += cxS*(-Cx*two13*pow(rb,4.0/3.0)); d_rb += cxS*(4.0/3.0)*(-Cx*two13)*cbrt(rb); }
    }
    double ea,dna,dsa; b88_x_spin(ra, saa, ea, dna, dsa);
    double eb,dnb,dsb; b88_x_spin(rb, sbb, eb, dnb, dsb);
    E += cxB*(ea+eb); d_ra += cxB*dna; d_rb += cxB*dnb; d_saa += cxB*dsa; d_sbb += cxB*dsb;
  }
  // ----- correlation -----
  if (func == FUNC_SVWN) {
    double ec,vca,vcb;
    vwn_spin_param(ra,rb, 0.0310907,3.72744,12.9352,-0.10498,
                          0.01554535,7.06042,18.0578,-0.32500,
                          -1.0/(6.0*pi*pi),1.13107,13.0045,-0.00475840, ec,vca,vcb);
    E += (ra+rb)*ec; d_ra += vca; d_rb += vcb;
  } else if (func == FUNC_B3LYP) {
    double ec,vca,vcb;            // 0.19 VWN1-RPA (formula I: para+ferro, no spin-stiffness)
    vwn1_spin(ra,rb, 0.0310907,13.0720,42.7198,-0.409286,
                     0.01554535,20.1231,101.578,-0.743294, ec,vca,vcb);
    E += 0.19*(ra+rb)*ec; d_ra += 0.19*vca; d_rb += 0.19*vcb;
    double L,dla,dlb,lsaa,lsab,lsbb;  lyp_c(ra,rb,saa,sab,sbb, L,dla,dlb,lsaa,lsab,lsbb);
    E += 0.81*L; d_ra += 0.81*dla; d_rb += 0.81*dlb; d_saa += 0.81*lsaa; d_sbb += 0.81*lsbb; d_sab += 0.81*lsab;
  } else if (func == FUNC_BLYP || func == FUNC_BHHLYP) {
    double L,dla,dlb,lsaa,lsab,lsbb;  lyp_c(ra,rb,saa,sab,sbb, L,dla,dlb,lsaa,lsab,lsbb);
    E += L; d_ra += dla; d_rb += dlb; d_saa += lsaa; d_sbb += lsbb; d_sab += lsab;
  }
  excd[p]   = E;
  vra_w[p]  = w[p]*d_ra;  vrb_w[p]  = w[p]*d_rb;
  vsaa_w[p] = w[p]*d_saa; vsbb_w[p] = w[p]*d_sbb; vsab_w[p] = w[p]*d_sab;
}

// ---- one UKS spin-channel contraction vector ----
// Z[u] = 0.5*vrho_w*phi[u] + (2*vs_same_w*drho_same + vs_cross_w*drho_cross).gphi[u]
// UKS same/cross variant, coalesced element-wise (see k_makeZ_gga).
__global__ void k_makeZ_gga_ch(int npts, int nbf, int gga,
        const double* vrho_w, const double* vs_same_w, const double* vs_cross_w,
        const double* drho_same, const double* drho_cross,
        const double* phi, const double* gphi, long gstr, double* Z) {
  const long nbL=nbf, npL=npts, tot=nbL*npL;
  const long idx = (long)blockIdx.x*blockDim.x + threadIdx.x;
  if (idx >= tot) return;
  const int p = (int)(idx / nbL);
  double z = 0.5*vrho_w[p]*phi[idx];
  if (gga) {
    const double a = 2.0*vs_same_w[p], b = vs_cross_w[p];
    const double cx = a*drho_same[p]       + b*drho_cross[p];
    const double cy = a*drho_same[p+npL]   + b*drho_cross[p+npL];
    const double cz = a*drho_same[p+2*npL] + b*drho_cross[p+2*npL];
    z += cx*gphi[idx] + cy*gphi[idx+gstr] + cz*gphi[idx+2*gstr];
  }
  Z[idx] = z;
}

// UKS energy + electron-count reduction: Exc = sum w*excd ; totele = sum w*(ra+rb)
__global__ void k_reduce_uks(int npts, const double* excd, const double* rho_a,
                             const double* rho_b, const double* w,
                             double* eexc, double* totele) {
  __shared__ double se[256], st[256];
  int t = threadIdx.x;
  double e=0.0, n=0.0;
  for (int p = t; p < npts; p += blockDim.x) {
    double wp = w[p];
    e += wp * excd[p];
    n += wp * (rho_a[p] + rho_b[p]);
  }
  se[t]=e; st[t]=n; __syncthreads();
  for (int s=blockDim.x/2; s>0; s>>=1){ if(t<s){se[t]+=se[t+s]; st[t]+=st[t+s];} __syncthreads(); }
  if (t==0){ atomicAdd(eexc, se[0]); atomicAdd(totele, st[0]); }  // += so grid blocks accumulate
}

// Z[u,p] = 0.5*vrho*phi[u,p]  +  (c.aoG1)[u,p]
//   c = 2*d1ds(ga)*drhoa + d1ds(gc)*drhob ; closed shell drhob=drhoa, so
//   c = (2*vsa + vsc)*drhoa  (weights already folded into vsa,vsc).
// One thread per (u,p) ELEMENT (idx = u + nbf*p) so the phi/gphi/Z accesses are
// COALESCED (adjacent threads -> adjacent memory). The per-grid-point coeffs
// (vrho_w[p], drho[p]) are broadcast-read across the nbf threads sharing a p.
// (Was: one thread per grid point looping u -> stride-nbf uncoalesced, the DFT
// bottleneck -- ~96 ms/iter at (H2O)16.)
__global__ void k_makeZ_gga(int npts, int nbf, int gga,
        const double* vrho_w, const double* vsa_w, const double* vsc_w,
        const double* drho, const double* phi, const double* gphi, long gstr,
        double* Z) {
  const long nbL=nbf, npL=npts, tot=nbL*npL;
  const long idx = (long)blockIdx.x*blockDim.x + threadIdx.x;
  if (idx >= tot) return;
  const int p = (int)(idx / nbL);
  double z = 0.5*vrho_w[p]*phi[idx];
  if (gga) {
    const double cc = 2.0*vsa_w[p] + vsc_w[p];
    z += cc*(drho[p]*gphi[idx] + drho[p+npL]*gphi[idx+gstr] + drho[p+2*npL]*gphi[idx+2*gstr]);
  }
  Z[idx] = z;
}

__global__ void k_reduce(int npts, const double* eps, const double* rho_tot,
                         const double* w, double* eexc, double* totele) {
  __shared__ double se[256], st[256];
  int t = threadIdx.x;
  double e=0.0, n=0.0;
  for (int p = t; p < npts; p += blockDim.x) {
    double wp = w[p], rt = rho_tot[p];
    e += wp * rt * eps[p];
    n += wp * rt;
  }
  se[t]=e; st[t]=n; __syncthreads();
  for (int s=blockDim.x/2; s>0; s>>=1){ if(t<s){se[t]+=se[t+s]; st[t]+=st[t+s];} __syncthreads(); }
  if (t==0){ atomicAdd(eexc, se[0]); atomicAdd(totele, st[0]); }  // += so grid blocks accumulate
}

// in-place symmetrization V <- V + V^T. The Vxc assembly is V = phi Z^T + Z phi^T
// = A + A^T with A = phi Z^T, so accumulating ONLY A per grid block (one GEMM
// instead of two) and symmetrizing once after the loop halves the V-assembly
// GEMM cost -- exactly, no approximation.
__global__ void k_sym_inplace(int nbf, double* V) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  if (i >= nbf || j >= nbf || j > i) return;    // pairs j <= i, each touched once
  double s = V[i + (long)nbf*j] + V[j + (long)nbf*i];
  V[i + (long)nbf*j] = s;
  V[j + (long)nbf*i] = s;
}

__global__ void k_unpack_fold(int nbf, const double* dpack, const double* bfnrm, double* Dfull) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  if (i >= nbf || j >= nbf) return;
  int a = i, b = j;
  if (a < b) { int t=a; a=b; b=t; }
  long idx = (long)a*(a+1)/2 + b;
  Dfull[i + (long)nbf*j] = bfnrm[i]*bfnrm[j]*dpack[idx];
}

__global__ void k_pack_fold(int nbf, const double* Vfull, const double* bfnrm, double* fxc) {
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i >= nbf) return;
  long base = (long)i*(i+1)/2;
  for (int jr = 0; jr <= i; jr++) {
    fxc[base + jr] = Vfull[jr + (long)nbf*i] * bfnrm[i]*bfnrm[jr];
  }
}

// ---- SPARSE: collocate ONLY a block's active shells into COMPACT rows --------
// Same math as k_collocate_g, but the output row is the shell's compact row
// (bsh_row) not its global AO, and strides are nact/bp (not nbf/npts).
__global__ void k_collocate_compact(int bp, int nact, int gga,
    const double* __restrict__ gx, const double* __restrict__ gy, const double* __restrict__ gz,
    const int* bsh_gid, const int* bsh_row, int nsh_b,
    const int* sh_am, const int* sh_g0, const int* sh_nc,
    const double* sh_cx, const double* sh_cy, const double* sh_cz, const double* sh_md2,
    const double* ex, const double* cc, const double* pmd2,
    const int* cartx, const int* carty, const int* cartz, int maxcart,
    double* __restrict__ phi, double* __restrict__ gphi)
{
  int p = blockIdx.x*blockDim.x + threadIdx.x;
  if (p >= bp) return;
  const long bpL=bp, nactL=nact;
  double px=gx[p], py=gy[p], pz=gz[p];
  for (int s=0; s<nsh_b; s++) {
    int ish=bsh_gid[s], row=bsh_row[s];
    int am=sh_am[ish], g0=sh_g0[ish], nc=sh_nc[ish];
    double dx=px-sh_cx[ish], dy=py-sh_cy[ish], dz=pz-sh_cz[ish];
    double r2=dx*dx+dy*dy+dz*dz;
    int ncart=(am+1)*(am+2)/2;
    if (r2 > sh_md2[ish]) {
      for (int ic=0; ic<ncart; ic++){ long o=(row+ic)+nactL*p; phi[o]=0.0;
        if(gga){gphi[o]=0.0; gphi[o+nactL*bpL]=0.0; gphi[o+2*nactL*bpL]=0.0;} }
      continue;
    }
    double vexp1=0.0, vexp2=0.0;
    for (int k=0;k<nc;k++){ int ig=g0-1+k; if(r2>pmd2[ig]) continue;
      double v=exp(-ex[ig]*r2)*cc[ig]; vexp1+=v; vexp2+=2.0*ex[ig]*v; }
    if (am==0) {
      long o=row+nactL*p; phi[o]=vexp1;
      if(gga){gphi[o]=-vexp2*dx; gphi[o+nactL*bpL]=-vexp2*dy; gphi[o+2*nactL*bpL]=-vexp2*dz;}
    } else if (am==1) {
      long o0=row+nactL*p, o1=o0+1, o2=o0+2;
      phi[o0]=vexp1*dx; phi[o1]=vexp1*dy; phi[o2]=vexp1*dz;
      if(gga){
        gphi[o0]=vexp1-vexp2*dx*dx; gphi[o1]=-vexp2*dx*dy; gphi[o2]=-vexp2*dx*dz;
        gphi[o0+nactL*bpL]=-vexp2*dx*dy; gphi[o1+nactL*bpL]=vexp1-vexp2*dy*dy; gphi[o2+nactL*bpL]=-vexp2*dy*dz;
        gphi[o0+2*nactL*bpL]=-vexp2*dx*dz; gphi[o1+2*nactL*bpL]=-vexp2*dy*dz; gphi[o2+2*nactL*bpL]=vexp1-vexp2*dz*dz;
      }
    } else {
      double drx[12],dry[12],drz[12];
      drx[0]=dry[0]=drz[0]=0.0; drx[1]=dry[1]=drz[1]=1.0; drx[2]=dx; dry[2]=dy; drz[2]=dz;
      for(int i=2;i<=am+1;i++){ drx[i+1]=drx[i]*dx; dry[i+1]=dry[i]*dy; drz[i+1]=drz[i]*dz; }
      for(int ic=0; ic<ncart; ic++){
        int ix=cartx[ic+maxcart*am], iy=carty[ic+maxcart*am], iz=cartz[ic+maxcart*am];
        double X=drx[ix+1], Y=dry[iy+1], Z=drz[iz+1];
        long o=(row+ic)+nactL*p;
        phi[o]=vexp1*X*Y*Z;
        if(gga){
          double xm=ix*drx[ix], ym=iy*dry[iy], zm=iz*drz[iz];
          double xp=drx[ix+2]*vexp2, yp=dry[iy+2]*vexp2, zp=drz[iz+2]*vexp2;
          gphi[o]           =(-xp+vexp1*xm)*Y*Z;
          gphi[o+nactL*bpL] =(-yp+vexp1*ym)*X*Z;
          gphi[o+2*nactL*bpL]=(-zp+vexp1*zm)*X*Y;
        }
      }
    }
  }
}

// Dc[i,j] = Dfold[ao[i], ao[j]]  (gather the active nact x nact sub-block).
__global__ void k_gather_D(int nact, int nb, const int* ao, const double* Dfold, double* Dc){
  int i=blockIdx.x*blockDim.x+threadIdx.x, j=blockIdx.y*blockDim.y+threadIdx.y;
  if(i>=nact||j>=nact) return;
  Dc[i+(long)nact*j] = Dfold[ao[i]+(long)nb*ao[j]];
}

// Vfull[ao[i], ao[j]] += Vc[i,j]. Grid blocks run sequentially and each (i,j)
// maps to a unique global pair within a block, so no atomics are needed.
__global__ void k_scatter_add_V(int nact, int nb, const int* ao, const double* Vc, double* Vfull){
  int i=blockIdx.x*blockDim.x+threadIdx.x, j=blockIdx.y*blockDim.y+threadIdx.y;
  if(i>=nact||j>=nact) return;
  Vfull[ao[i]+(long)nb*ao[j]] += Vc[i+(long)nact*j];
}

// =============================================================================
static long long rd_i8(FILE* f){ long long v; size_t r=fread(&v,8,1,f); (void)r; return v; }
static void rd_dvec(FILE* f, std::vector<double>& v, long n){ v.resize(n); size_t r=fread(v.data(),8,n,f); (void)r; }

static int parse_func(const char* fn) {
  if (!fn) return FUNC_SLATER;
  if (strcmp(fn,"SVWN")==0 || strcmp(fn,"SVWN5")==0) return FUNC_SVWN;
  if (strcmp(fn,"BLYP")==0) return FUNC_BLYP;
  if (strcmp(fn,"B3LYP")==0 || strcmp(fn,"B3LYPV1R")==0) return FUNC_B3LYP;
  if (strcmp(fn,"BHHLYP")==0 || strcmp(fn,"BHANDHLYP")==0) return FUNC_BHHLYP;
  return FUNC_SLATER;
}

static bool init_ctx(Ctx& g, int nfocks) {   // g SHADOWS the global: per-device ctx
  const char* dir = getenv("OQP_OWNXC_DIR");
  if (!dir) { fprintf(stderr,"[ownxc] OQP_OWNXC_DIR not set\n"); return false; }
  if (cudaSetDevice(g.dev)!=cudaSuccess) { fprintf(stderr,"[ownxc] cudaSetDevice(%d) failed\n",g.dev); return false; }
  g.func = parse_func(getenv("OQP_OWNXC_FUNC"));
  // host copies kept for the sparse block-screening step (built after nb/np known)
  std::vector<double> H_gx,H_gy,H_gz,H_gw, H_cx,H_cy,H_cz,H_md2, H_ex,H_cc;
  std::vector<int>    H_am,H_ao,H_nc,H_g0;
  std::string D = dir;
  {
    FILE* f = fopen((D+"/grid.bin").c_str(),"rb");
    if(!f){ fprintf(stderr,"[ownxc] no grid.bin\n"); return false; }
    long npts = (long)rd_i8(f);
    std::vector<double> xyz; rd_dvec(f, xyz, 3L*npts);
    std::vector<double> w;   rd_dvec(f, w, npts);
    fclose(f);
    // multi-GPU: this ctx integrates only its slice of the grid (exact split --
    // the XC energy/potential are plain sums over points).
    const long pa = (npts*(long)g.rank)/g.nranks, pb = (npts*(long)(g.rank+1))/g.nranks;
    const long nps = pb - pa;
    g.npts = (int)nps;
    H_gx.resize(nps); H_gy.resize(nps); H_gz.resize(nps); H_gw.resize(nps);
    for (long p=0;p<nps;p++){ H_gx[p]=xyz[3*(pa+p)+0]; H_gy[p]=xyz[3*(pa+p)+1]; H_gz[p]=xyz[3*(pa+p)+2]; H_gw[p]=w[pa+p]; }
    CK(cudaMalloc(&g.d_gx, nps*8)); CK(cudaMemcpy(g.d_gx,H_gx.data(),nps*8,cudaMemcpyHostToDevice));
    CK(cudaMalloc(&g.d_gy, nps*8)); CK(cudaMemcpy(g.d_gy,H_gy.data(),nps*8,cudaMemcpyHostToDevice));
    CK(cudaMalloc(&g.d_gz, nps*8)); CK(cudaMemcpy(g.d_gz,H_gz.data(),nps*8,cudaMemcpyHostToDevice));
    CK(cudaMalloc(&g.d_gw, nps*8)); CK(cudaMemcpy(g.d_gw,H_gw.data(),nps*8,cudaMemcpyHostToDevice));
  }
  {
    FILE* f = fopen((D+"/basis.bin").c_str(),"rb");
    if(!f){ fprintf(stderr,"[ownxc] no basis.bin\n"); return false; }
    g.nshell=(int)rd_i8(f); g.nbf=(int)rd_i8(f); g.nprim=(int)rd_i8(f);
    H_am.resize(g.nshell); H_ao.resize(g.nshell); H_nc.resize(g.nshell); H_g0.resize(g.nshell);
    H_cx.resize(g.nshell); H_cy.resize(g.nshell); H_cz.resize(g.nshell); H_md2.resize(g.nshell);
    auto& am=H_am; auto& nc=H_nc; auto& ao=H_ao; auto& g0=H_g0;
    auto& cx=H_cx; auto& cy=H_cy; auto& cz=H_cz; auto& md2=H_md2;
    for (int s=0;s<g.nshell;s++){
      am[s]=(int)rd_i8(f); rd_i8(f); g0[s]=(int)rd_i8(f); nc[s]=(int)rd_i8(f);
      ao[s]=(int)rd_i8(f); rd_i8(f);
      double c3[3]; size_t r=fread(c3,8,3,f);(void)r; cx[s]=c3[0]; cy[s]=c3[1]; cz[s]=c3[2];
      size_t r2=fread(&md2[s],8,1,f);(void)r2;
    }
    H_ex.clear(); H_cc.clear();
    std::vector<double> pmd2,bfn;
    rd_dvec(f,H_ex,g.nprim); rd_dvec(f,H_cc,g.nprim); rd_dvec(f,pmd2,g.nprim); rd_dvec(f,bfn,g.nbf);
    fclose(f);
    auto upI=[&](int**dp,std::vector<int>&v){ CK(cudaMalloc(dp,v.size()*4)); CK(cudaMemcpy(*dp,v.data(),v.size()*4,cudaMemcpyHostToDevice)); return true; };
    auto upD=[&](double**dp,std::vector<double>&v){ CK(cudaMalloc(dp,v.size()*8)); CK(cudaMemcpy(*dp,v.data(),v.size()*8,cudaMemcpyHostToDevice)); return true; };
    if(!upI(&g.d_sh_am,am)||!upI(&g.d_sh_g0,g0)||!upI(&g.d_sh_nc,nc)||!upI(&g.d_sh_ao,ao)) return false;
    if(!upD(&g.d_sh_cx,cx)||!upD(&g.d_sh_cy,cy)||!upD(&g.d_sh_cz,cz)||!upD(&g.d_sh_md2,md2)) return false;
    if(!upD(&g.d_ex,H_ex)||!upD(&g.d_cc,H_cc)||!upD(&g.d_pmd2,pmd2)||!upD(&g.d_bfnrm,bfn)) return false;
  }
  {
    FILE* f = fopen((D+"/cart.bin").c_str(),"rb");
    if(!f){ fprintf(stderr,"[ownxc] no cart.bin\n"); return false; }
    g.maxang=(int)rd_i8(f); g.maxcart=(int)rd_i8(f);
    long n = (long)g.maxcart*(g.maxang+1);
    std::vector<long long> tx(n),ty(n),tz(n);
    size_t r1=fread(tx.data(),8,n,f); size_t r2=fread(ty.data(),8,n,f); size_t r3=fread(tz.data(),8,n,f);
    (void)r1;(void)r2;(void)r3; fclose(f);
    std::vector<int> ix(n),iy(n),iz(n);
    for(long k=0;k<n;k++){ ix[k]=(int)tx[k]; iy[k]=(int)ty[k]; iz[k]=(int)tz[k]; }
    CK(cudaMalloc(&g.d_cartx,n*4)); CK(cudaMemcpy(g.d_cartx,ix.data(),n*4,cudaMemcpyHostToDevice));
    CK(cudaMalloc(&g.d_carty,n*4)); CK(cudaMemcpy(g.d_carty,iy.data(),n*4,cudaMemcpyHostToDevice));
    CK(cudaMalloc(&g.d_cartz,n*4)); CK(cudaMemcpy(g.d_cartz,iz.data(),n*4,cudaMemcpyHostToDevice));
  }
  long nb=g.nbf, np=g.npts;
  bool gga = func_is_gga(g.func);
  // SELECT the XC path. Small systems: dense with the FULL grid cached once
  // (collocate-once, big efficient GEMMs -- fastest when it fits). Large systems
  // that would otherwise TILE (re-collocate the full basis every iteration): the
  // active-AO SPARSE path wins (only a block's few reachable shells are touched,
  // O(N^3)->O(N^2); measured 2.6x at (H2O)32 and growing with size). AUTO picks
  // sparse exactly when the dense cache would not fit; env forces either way.
  { size_t freeb=0, totb=0; cudaMemGetInfo(&freeb,&totb);
    const double full = (double)nb*np*8.0 * (gga?6.0:4.0);
    const bool fits = full < 0.60*(double)freeb;
    if (getenv("OQP_OWNXC_SPARSE"))     g.sparse = true;
    else if (getenv("OQP_OWNXC_DENSE")) g.sparse = false;
    else                                g.sparse = !fits;   // auto (RKS + UKS/ROHF)
    (void)nfocks;
  }
  if (g.sparse) {
    long sb = getenv("OQP_OWNXC_SBLK") ? atol(getenv("OQP_OWNXC_SBLK")) : 4096;
    if (sb < 256) sb = 256;
    g.sblk = (int)std::min<long>(sb, np);
    g.blk = g.sblk; g.bwork = g.sblk; g.phi_cached = false;
  } else {
    // Cache the FULL grid (collocate once) if phi/gphi/X/Z fit comfortably; else
    // tile. BLK default 65536; OQP_OWNXC_TILE forces tiling.
    size_t freeb=0, totb=0; cudaMemGetInfo(&freeb,&totb);
    const double full = (double)nb*np*8.0 * (gga?6.0:4.0);
    long bkenv = getenv("OQP_OWNXC_BLK") ? atol(getenv("OQP_OWNXC_BLK")) : 65536;
    if (bkenv < 1024) bkenv = 1024;
    if (!getenv("OQP_OWNXC_TILE") && full < 0.60*(double)freeb) { g.blk=(int)np; g.phi_cached=true; }
    else { g.blk=(int)std::min<long>(bkenv, np); g.phi_cached=false; }
    g.bwork = (int)std::min<long>(bkenv, np);
  }
  const long bw = g.bwork;      // work buffers (rho/vrho/... and X/Z tiles) are bw wide
  // buffers used by BOTH paths
  CK(cudaMalloc(&g.d_Dfold,nb*nb*8));
  CK(cudaMalloc(&g.d_rho,  bw*8));
  if (gga) CK(cudaMalloc(&g.d_drho, bw*3*8));
  CK(cudaMalloc(&g.d_vrho, bw*8));
  CK(cudaMalloc(&g.d_vsa,  bw*8));
  CK(cudaMalloc(&g.d_vsc,  bw*8));
  CK(cudaMalloc(&g.d_eps,  bw*8));
  CK(cudaMalloc(&g.d_Vfull,nb*nb*8));
  CK(cudaMalloc(&g.d_eexc, 8));
  CK(cudaMalloc(&g.d_totele,8));
  cublasCreate(&g.cub);
  if (!g.sparse) {
    const long bk = g.blk;
    CK(cudaMalloc(&g.d_phi,  nb*bk*8));
    if (gga) CK(cudaMalloc(&g.d_gphi, nb*bk*3*8));
    CK(cudaMalloc(&g.d_X,    nb*bw*8));
    CK(cudaMalloc(&g.d_Z,    nb*bw*8));
    if (g.phi_cached) {   // collocate the whole grid once (geometry-only)
      int TB=128, nblk=(int)((np+TB-1)/TB);
      k_collocate_g<<<nblk,TB>>>((int)np,g.nshell,nb, gga?1:0, g.d_gx,g.d_gy,g.d_gz,
         g.d_sh_am,g.d_sh_g0,g.d_sh_nc,g.d_sh_ao, g.d_sh_cx,g.d_sh_cy,g.d_sh_cz,g.d_sh_md2,
         g.d_ex,g.d_cc,g.d_pmd2, g.d_cartx,g.d_carty,g.d_cartz,g.maxcart, g.d_phi, g.d_gphi);
      CKv(cudaDeviceSynchronize());
    }
  } else {
    // ---- (a) real per-shell cutoff r^2 (the prep writes md2=1e30 = no screen).
    // A shell reaches r where its slowest primitive |c| r^l exp(-a r^2) > tol. --
    { double tol = getenv("OQP_OWNXC_SCRTOL") ? atof(getenv("OQP_OWNXC_SCRTOL")) : 1e-11;
      for (int s=0;s<g.nshell;s++){
        double amin=1e300, cabs=0.0;
        for (int k=0;k<H_nc[s];k++){ int ig=H_g0[s]-1+k; if(H_ex[ig]<amin){amin=H_ex[ig]; cabs=fabs(H_cc[ig]);} }
        if (amin>1e299 || cabs<=0.0) { H_md2[s]=1e30; continue; }
        double L=log(cabs/tol), l=(double)H_am[s];
        double u=L/amin;
        for (int it=0; it<4; it++) u=(L + 0.5*l*log(u>1e-6?u:1e-6))/amin;
        H_md2[s]=(u>0?u:0.0)*1.20;
      }
      CK(cudaMemcpy(g.d_sh_md2, H_md2.data(), (size_t)g.nshell*8, cudaMemcpyHostToDevice)); }
    // ---- (b) reorder the grid into spatially-local blocks (coarse-cell sort) so
    // a contiguous sblk-block touches few shells. XC = sum over points -> exact. --
    { const long npL=np; double mn[3]={1e300,1e300,1e300};
      for (long p=0;p<npL;p++){ mn[0]=std::min(mn[0],H_gx[p]); mn[1]=std::min(mn[1],H_gy[p]); mn[2]=std::min(mn[2],H_gz[p]); }
      const double cs=3.0; std::vector<long long> key(npL);
      for (long p=0;p<npL;p++){ long long ix=(long long)((H_gx[p]-mn[0])/cs), iy=(long long)((H_gy[p]-mn[1])/cs), iz=(long long)((H_gz[p]-mn[2])/cs);
        key[p]=(ix*100000LL+iy)*100000LL+iz; }
      std::vector<long> perm(npL); for(long p=0;p<npL;p++) perm[p]=p;
      std::stable_sort(perm.begin(),perm.end(),[&](long a,long b){return key[a]<key[b];});
      std::vector<double> tx(npL),ty(npL),tz(npL),tw(npL);
      for(long p=0;p<npL;p++){ long s=perm[p]; tx[p]=H_gx[s]; ty[p]=H_gy[s]; tz[p]=H_gz[s]; tw[p]=H_gw[s]; }
      H_gx.swap(tx); H_gy.swap(ty); H_gz.swap(tz); H_gw.swap(tw);
      CK(cudaMemcpy(g.d_gx,H_gx.data(),npL*8,cudaMemcpyHostToDevice));
      CK(cudaMemcpy(g.d_gy,H_gy.data(),npL*8,cudaMemcpyHostToDevice));
      CK(cudaMemcpy(g.d_gz,H_gz.data(),npL*8,cudaMemcpyHostToDevice));
      CK(cudaMemcpy(g.d_gw,H_gw.data(),npL*8,cudaMemcpyHostToDevice)); }
    // ---- (c) build per-block active-shell lists (host screening) --------------
    const int nsh = g.nshell, sblk = g.sblk;
    g.nblock = (int)((np + sblk - 1) / sblk);
    g.h_bsh_off.assign(g.nblock+1, 0);
    g.h_ao_off.assign(g.nblock+1, 0);
    std::vector<int> bsh_gid, bsh_row, ao_idx;
    for (int b=0; b<g.nblock; b++) {
      long p0=(long)b*sblk, p1=std::min<long>(p0+sblk, np);
      int row=0;
      for (int s=0; s<nsh; s++) {
        double cxs=H_cx[s], cys=H_cy[s], czs=H_cz[s], m2=H_md2[s];
        double minr2=1e300;
        for (long p=p0; p<p1; p++) {
          double dx=H_gx[p]-cxs, dy=H_gy[p]-cys, dz=H_gz[p]-czs;
          double r2=dx*dx+dy*dy+dz*dz; if (r2<minr2) minr2=r2;
          if (minr2<=m2) break;                       // already active, stop scanning
        }
        if (minr2<=m2) {
          int ncart=(H_am[s]+1)*(H_am[s]+2)/2;
          bsh_gid.push_back(s); bsh_row.push_back(row);
          for (int ic=0; ic<ncart; ic++) ao_idx.push_back(H_ao[s]-1+ic);
          row += ncart;
        }
      }
      g.h_bsh_off[b+1]=(int)bsh_gid.size();
      g.h_ao_off[b+1]=(int)ao_idx.size();
      if (row>g.nact_max) g.nact_max=row;
    }
    auto upI=[&](int**dp,std::vector<int>&v){ if(v.empty())v.push_back(0);
      CK(cudaMalloc(dp,v.size()*4)); CK(cudaMemcpy(*dp,v.data(),v.size()*4,cudaMemcpyHostToDevice)); return true; };
    if(!upI(&g.d_bsh_gid,bsh_gid)||!upI(&g.d_bsh_row,bsh_row)||!upI(&g.d_ao_idx,ao_idx)) return false;
    const long nm=g.nact_max, sb=g.sblk;
    CK(cudaMalloc(&g.d_phic, nm*sb*8));
    if (gga) CK(cudaMalloc(&g.d_gphic, nm*sb*3*8));
    CK(cudaMalloc(&g.d_Xc,  nm*sb*8));
    CK(cudaMalloc(&g.d_Zc,  nm*sb*8));
    CK(cudaMalloc(&g.d_Dc,  nm*nm*8));
    CK(cudaMalloc(&g.d_Vc,  nm*nm*8));
    double meanact=0; for(int b=0;b<g.nblock;b++) meanact+=(g.h_ao_off[b+1]-g.h_ao_off[b]);
    meanact/=std::max(1,g.nblock);
    fprintf(stderr,"[ownxc-sparse] nblock=%d sblk=%d nact_max=%d mean_nact=%.0f/%d (%.1f%%)\n",
            g.nblock,g.sblk,g.nact_max,meanact,g.nbf,100.0*meanact/g.nbf);
  }
  g.ready=true;
  const char* names[5]={"SLATER","SVWN","BLYP","B3LYP","BHHLYP"};
  fprintf(stderr,"[ownxc] init ok: nbf=%d nshell=%d nprim=%d npts=%d func=%s (gga=%d) blk=%d bwork=%d phi_cached=%d\n",
          g.nbf,g.nshell,g.nprim,g.npts, names[g.func], gga?1:0, g.blk, g.bwork, g.phi_cached?1:0);
  return true;
}

// beta-channel buffers, only when a UKS/ROHF (nfocks==2) call actually happens.
static bool ensure_uks(Ctx& g) {   // g shadows the global (per-device ctx)
  if (g.uks_alloc) return true;
  const long nb=g.nbf, bw=g.bwork;
  const bool gga = func_is_gga(g.func);
  CK(cudaMalloc(&g.d_Dfold_b, nb*nb*8));      // both paths: gather source / dense D
  CK(cudaMalloc(&g.d_rho_b,   bw*8));
  if (gga) CK(cudaMalloc(&g.d_drho_b, bw*3*8));
  CK(cudaMalloc(&g.d_vrb,     bw*8));
  CK(cudaMalloc(&g.d_vsbb,    bw*8));
  CK(cudaMalloc(&g.d_Vfull_b, nb*nb*8));      // both paths: scatter target / dense V
  CK(cudaMalloc(&g.d_excd,    bw*8));
  if (!g.sparse) {
    CK(cudaMalloc(&g.d_Xb, nb*bw*8));
    CK(cudaMalloc(&g.d_Zb, nb*bw*8));
  } else {
    const long nm=g.nact_max, sb=g.sblk;
    CK(cudaMalloc(&g.d_Xcb, nm*sb*8));
    CK(cudaMalloc(&g.d_Zcb, nm*sb*8));
    CK(cudaMalloc(&g.d_Dcb, nm*nm*8));
    CK(cudaMalloc(&g.d_Vcb, nm*nm*8));
  }
  g.uks_alloc = true;
  return true;
}

static void vxc_one(Ctx& g, const double* d, double* fxc, const int* nbf,
                    const int* nfocks, double* eexc, double* totele, int* info) {
  *info = 1;
  if (*nfocks != 1 && *nfocks != 2) return;   // RHF (1) or UKS/ROHF (2)
  if (cudaSetDevice(g.dev)!=cudaSuccess) return;
  if (!g.ready) { if(!init_ctx(g,*nfocks)){ return; } }
  if (*nbf != g.nbf) { fprintf(stderr,"[ownxc] nbf mismatch %d vs %d\n",*nbf,g.nbf); return; }
  int nb=g.nbf, np=g.npts;
  bool gga = func_is_gga(g.func);
  long ntri = (long)nb*(nb+1)/2;

  // ===========================================================================
  // open-shell UKS (nfocks==2): d = [d_alpha_tri | d_beta_tri], fxc=[Va|Vb].
  // GGA cross-coupling: c_a = 2*vsaa*drho_a + vsab*drho_b (and a<->b for c_b).
  // RKS (nfocks==1) path below is byte-identical.
  // ===========================================================================
  if (*nfocks == 2) {
    if (!ensure_uks(g)) return;
    double* d2=nullptr;
    if (cudaMalloc(&d2, 2*ntri*8)!=cudaSuccess) return;
    cudaMemcpy(d2, d, 2*ntri*8, cudaMemcpyHostToDevice);
    { dim3 TB(16,16), nb2((nb+15)/16,(nb+15)/16);
      k_unpack_fold<<<nb2,TB>>>(nb, d2,        g.d_bfnrm, g.d_Dfold);     // alpha
      k_unpack_fold<<<nb2,TB>>>(nb, d2+ntri,   g.d_bfnrm, g.d_Dfold_b);   // beta
    }
    cudaFree(d2);
    cudaMemset(g.d_Vfull,0,(size_t)nb*nb*8); cudaMemset(g.d_Vfull_b,0,(size_t)nb*nb*8);
    cudaMemset(g.d_eexc,0,8); cudaMemset(g.d_totele,0,8);
    if (g.sparse) {
      // ---- active-AO O(N^2) UKS: two density channels share the compact phi ----
      for (int b=0; b<g.nblock; b++) {
        const long p0=(long)b*g.sblk; const int bp=(int)std::min<long>(g.sblk, np-p0);
        const int aoff=g.h_ao_off[b], nact=g.h_ao_off[b+1]-aoff;
        const int soff=g.h_bsh_off[b], nsh_b=g.h_bsh_off[b+1]-soff;
        if (nact==0) continue;
        const long gstr=(long)nact*bp; const int* aoi=g.d_ao_idx+aoff;
        { int TB=128, nblk=(bp+TB-1)/TB;
          k_collocate_compact<<<nblk,TB>>>(bp,nact,gga?1:0, g.d_gx+p0,g.d_gy+p0,g.d_gz+p0,
             g.d_bsh_gid+soff, g.d_bsh_row+soff, nsh_b,
             g.d_sh_am,g.d_sh_g0,g.d_sh_nc, g.d_sh_cx,g.d_sh_cy,g.d_sh_cz,g.d_sh_md2,
             g.d_ex,g.d_cc,g.d_pmd2, g.d_cartx,g.d_carty,g.d_cartz,g.maxcart, g.d_phic, g.d_gphic); }
        { dim3 TB(16,16), gb((nact+15)/16,(nact+15)/16);
          k_gather_D<<<gb,TB>>>(nact, nb, aoi, g.d_Dfold,   g.d_Dc);
          k_gather_D<<<gb,TB>>>(nact, nb, aoi, g.d_Dfold_b, g.d_Dcb); }
        { const double a=1.0, z=0.0;
          cublasDgemm(g.cub,CUBLAS_OP_N,CUBLAS_OP_N, nact,bp,nact, &a, g.d_Dc, nact, g.d_phic,nact, &z, g.d_Xc, nact);
          cublasDgemm(g.cub,CUBLAS_OP_N,CUBLAS_OP_N, nact,bp,nact, &a, g.d_Dcb,nact, g.d_phic,nact, &z, g.d_Xcb,nact); }
        { int TB=256; unsigned nw=(unsigned)(((long)bp*32+TB-1)/TB);
          k_rho<<<nw,TB>>>(bp, nact, g.d_phic, g.d_Xc,  g.d_rho);
          k_rho<<<nw,TB>>>(bp, nact, g.d_phic, g.d_Xcb, g.d_rho_b);
          if (gga) {
            k_drho<<<nw,TB>>>(bp, nact, g.d_gphic, gstr, g.d_Xc,  g.d_drho);
            k_drho<<<nw,TB>>>(bp, nact, g.d_gphic, gstr, g.d_Xcb, g.d_drho_b);
            const double two=2.0; cublasDscal(g.cub,3*bp,&two,g.d_drho,1); cublasDscal(g.cub,3*bp,&two,g.d_drho_b,1);
          } }
        { int TB=128, nblk=(bp+TB-1)/TB;
          k_func_gga_uks<<<nblk,TB>>>(bp, g.func, g.d_rho, g.d_rho_b, g.d_drho, g.d_drho_b, g.d_gw+p0,
                                      g.d_excd, g.d_vrho, g.d_vrb, g.d_vsa, g.d_vsbb, g.d_vsc); }
        { int TB=256; unsigned nz=(unsigned)(((long)nact*bp+TB-1)/TB);
          k_makeZ_gga_ch<<<nz,TB>>>(bp, nact, gga?1:0, g.d_vrho, g.d_vsa, g.d_vsc, g.d_drho, g.d_drho_b, g.d_phic, g.d_gphic, gstr, g.d_Zc);
          k_makeZ_gga_ch<<<nz,TB>>>(bp, nact, gga?1:0, g.d_vrb, g.d_vsbb, g.d_vsc, g.d_drho_b, g.d_drho, g.d_phic, g.d_gphic, gstr, g.d_Zcb); }
        { const double a=1.0, z=0.0;
          cublasDgemm(g.cub,CUBLAS_OP_N,CUBLAS_OP_T, nact,nact,bp, &a, g.d_phic,nact, g.d_Zc, nact, &z, g.d_Vc,  nact);
          cublasDgemm(g.cub,CUBLAS_OP_N,CUBLAS_OP_T, nact,nact,bp, &a, g.d_phic,nact, g.d_Zcb,nact, &z, g.d_Vcb, nact); }
        { dim3 TB(16,16), gb((nact+15)/16,(nact+15)/16);
          k_scatter_add_V<<<gb,TB>>>(nact, nb, aoi, g.d_Vc,  g.d_Vfull);
          k_scatter_add_V<<<gb,TB>>>(nact, nb, aoi, g.d_Vcb, g.d_Vfull_b); }
        k_reduce_uks<<<1,256>>>(bp, g.d_excd, g.d_rho, g.d_rho_b, g.d_gw+p0, g.d_eexc, g.d_totele);
      }
    } else
    for (long p0=0; p0<np; p0+=g.bwork) {
      const int bp = (int)std::min<long>(g.bwork, np-p0);
      // cached: phi/gphi hold the FULL grid, this iteration reads the bp-wide slice
      // at column p0 (contiguous, ld nb); gphi components are nb*np apart.
      // tiled: collocate this block into the [nb x bwork] buffers (components nb*bp apart).
      const double* phiB  = g.d_phi  + (g.phi_cached ? (size_t)p0*nb : 0);
      const double* gphiB = g.d_gphi ? g.d_gphi + (g.phi_cached ? (size_t)p0*nb : 0) : nullptr;
      const long    gstr  = (long)nb * (g.phi_cached ? (long)np : (long)bp);
      if (!g.phi_cached) { int TB=128, nblk=(bp+TB-1)/TB;
        k_collocate_g<<<nblk,TB>>>(bp,g.nshell,nb, gga?1:0, g.d_gx+p0,g.d_gy+p0,g.d_gz+p0,
           g.d_sh_am,g.d_sh_g0,g.d_sh_nc,g.d_sh_ao, g.d_sh_cx,g.d_sh_cy,g.d_sh_cz,g.d_sh_md2,
           g.d_ex,g.d_cc,g.d_pmd2, g.d_cartx,g.d_carty,g.d_cartz,g.maxcart, g.d_phi, g.d_gphi); }
      { const double a=1.0, b=0.0;
        cublasDgemm(g.cub, CUBLAS_OP_N, CUBLAS_OP_N, nb, bp, nb, &a, g.d_Dfold,   nb, phiB, nb, &b, g.d_X,  nb);
        cublasDgemm(g.cub, CUBLAS_OP_N, CUBLAS_OP_N, nb, bp, nb, &a, g.d_Dfold_b, nb, phiB, nb, &b, g.d_Xb, nb); }
      { int TB=256; unsigned nw=(unsigned)(((long)bp*32+TB-1)/TB);
        k_rho<<<nw,TB>>>(bp, nb, phiB, g.d_X,  g.d_rho);
        k_rho<<<nw,TB>>>(bp, nb, phiB, g.d_Xb, g.d_rho_b);
        if (gga) {
          k_drho<<<nw,TB>>>(bp, nb, gphiB, gstr, g.d_X,  g.d_drho);
          k_drho<<<nw,TB>>>(bp, nb, gphiB, gstr, g.d_Xb, g.d_drho_b);
          const double two=2.0;
          cublasDscal(g.cub, 3*bp, &two, g.d_drho,   1);
          cublasDscal(g.cub, 3*bp, &two, g.d_drho_b, 1);
        } }
      { int TB=128, nblk=(bp+TB-1)/TB;
        k_func_gga_uks<<<nblk,TB>>>(bp, g.func, g.d_rho, g.d_rho_b, g.d_drho, g.d_drho_b, g.d_gw+p0,
                                    g.d_excd, g.d_vrho, g.d_vrb, g.d_vsa, g.d_vsbb, g.d_vsc); }
      { int TB=256; unsigned nz=(unsigned)(((long)nb*bp+TB-1)/TB);
        k_makeZ_gga_ch<<<nz,TB>>>(bp, nb, gga?1:0, g.d_vrho, g.d_vsa, g.d_vsc, g.d_drho, g.d_drho_b, phiB, gphiB, gstr, g.d_Z);
        k_makeZ_gga_ch<<<nz,TB>>>(bp, nb, gga?1:0, g.d_vrb, g.d_vsbb, g.d_vsc, g.d_drho_b, g.d_drho, phiB, gphiB, gstr, g.d_Zb); }
      { const double a=1.0, b1=1.0;   // accumulate A = phi Z^T only; V = A + A^T after the loop
        cublasDgemm(g.cub, CUBLAS_OP_N, CUBLAS_OP_T, nb, nb, bp, &a, phiB, nb, g.d_Z,  nb, &b1, g.d_Vfull,   nb);
        cublasDgemm(g.cub, CUBLAS_OP_N, CUBLAS_OP_T, nb, nb, bp, &a, phiB, nb, g.d_Zb, nb, &b1, g.d_Vfull_b, nb); }
      k_reduce_uks<<<1,256>>>(bp, g.d_excd, g.d_rho, g.d_rho_b, g.d_gw+p0, g.d_eexc, g.d_totele);
    }
    { dim3 TB(16,16), nb2((nb+15)/16,(nb+15)/16);
      k_sym_inplace<<<nb2,TB>>>(nb, g.d_Vfull);       // V   <- A   + A^T
      k_sym_inplace<<<nb2,TB>>>(nb, g.d_Vfull_b); }   // V_b <- A_b + A_b^T
    double* d_fxc=nullptr;
    if (cudaMalloc(&d_fxc, 2*ntri*8)!=cudaSuccess) return;
    { int TB=128, nblk=(nb+TB-1)/TB;
      k_pack_fold<<<nblk,TB>>>(nb, g.d_Vfull,   g.d_bfnrm, d_fxc);
      k_pack_fold<<<nblk,TB>>>(nb, g.d_Vfull_b, g.d_bfnrm, d_fxc+ntri);
    }
    CKv(cudaDeviceSynchronize());
    cudaMemcpy(fxc,    d_fxc,     2*ntri*8, cudaMemcpyDeviceToHost);
    cudaMemcpy(eexc,   g.d_eexc,  8, cudaMemcpyDeviceToHost);
    cudaMemcpy(totele, g.d_totele,8, cudaMemcpyDeviceToHost);
    cudaFree(d_fxc);
    *info = 0;
    return;
  }

  double* d_dpack=nullptr;
  if (cudaMalloc(&d_dpack, ntri*8)!=cudaSuccess) return;
  cudaMemcpy(d_dpack, d, ntri*8, cudaMemcpyHostToDevice);
  {
    dim3 TB(16,16), nb2((nb+15)/16,(nb+15)/16);
    k_unpack_fold<<<nb2,TB>>>(nb, d_dpack, g.d_bfnrm, g.d_Dfold);
  }
  cudaFree(d_dpack);
  // ---- grid-blocked accumulation: Vfull and energy summed over bwork-wide grid
  // slices. Cached: phi/gphi hold the full grid and each slice is read in place
  // (X/Z/rho stay [nbf x bwork] -- no full-grid work buffers). Tiled: collocate
  // the slice into the [nbf x bwork] phi/gphi. ----
  cudaMemset(g.d_Vfull, 0, (size_t)nb*nb*8);
  cudaMemset(g.d_eexc, 0, 8); cudaMemset(g.d_totele, 0, 8);
  if (g.sparse) {
    // ---- active-AO O(N^2) path: per block collocate/contract only its shells --
    for (int b=0; b<g.nblock; b++) {
      const long p0=(long)b*g.sblk; const int bp=(int)std::min<long>(g.sblk, np-p0);
      const int aoff=g.h_ao_off[b], nact=g.h_ao_off[b+1]-aoff;
      const int soff=g.h_bsh_off[b], nsh_b=g.h_bsh_off[b+1]-soff;
      if (nact==0) continue;
      const long gstr=(long)nact*bp;
      { int TB=128, nblk=(bp+TB-1)/TB;
        k_collocate_compact<<<nblk,TB>>>(bp,nact,gga?1:0, g.d_gx+p0,g.d_gy+p0,g.d_gz+p0,
           g.d_bsh_gid+soff, g.d_bsh_row+soff, nsh_b,
           g.d_sh_am,g.d_sh_g0,g.d_sh_nc, g.d_sh_cx,g.d_sh_cy,g.d_sh_cz,g.d_sh_md2,
           g.d_ex,g.d_cc,g.d_pmd2, g.d_cartx,g.d_carty,g.d_cartz,g.maxcart, g.d_phic, g.d_gphic); }
      { dim3 TB(16,16), gb((nact+15)/16,(nact+15)/16);
        k_gather_D<<<gb,TB>>>(nact, nb, g.d_ao_idx+aoff, g.d_Dfold, g.d_Dc); }
      { const double a=1.0, z=0.0;   // Xc = Dc @ phic   (nact x bp)
        cublasDgemm(g.cub,CUBLAS_OP_N,CUBLAS_OP_N, nact,bp,nact, &a, g.d_Dc,nact, g.d_phic,nact, &z, g.d_Xc,nact); }
      { int TB=256; unsigned nw=(unsigned)(((long)bp*32+TB-1)/TB);
        k_rho<<<nw,TB>>>(bp, nact, g.d_phic, g.d_Xc, g.d_rho);
        if (gga) k_drho<<<nw,TB>>>(bp, nact, g.d_gphic, gstr, g.d_Xc, g.d_drho); }
      { int TB=128, nblk=(bp+TB-1)/TB;
        k_func_gga<<<nblk,TB>>>(bp, g.func, g.d_rho, g.d_drho, g.d_gw+p0, g.d_eps, g.d_vrho, g.d_vsa, g.d_vsc); }
      { int TB=256; unsigned nz=(unsigned)(((long)nact*bp+TB-1)/TB);
        k_makeZ_gga<<<nz,TB>>>(bp, nact, gga?1:0, g.d_vrho, g.d_vsa, g.d_vsc, g.d_drho, g.d_phic, g.d_gphic, gstr, g.d_Zc); }
      { const double a=1.0, z=0.0;   // Vc = phic Zc^T  (nact x nact, this block only)
        cublasDgemm(g.cub,CUBLAS_OP_N,CUBLAS_OP_T, nact,nact,bp, &a, g.d_phic,nact, g.d_Zc,nact, &z, g.d_Vc,nact); }
      { dim3 TB(16,16), gb((nact+15)/16,(nact+15)/16);
        k_scatter_add_V<<<gb,TB>>>(nact, nb, g.d_ao_idx+aoff, g.d_Vc, g.d_Vfull); }
      k_reduce<<<1,256>>>(bp, g.d_eps, g.d_rho, g.d_gw+p0, g.d_eexc, g.d_totele);
    }
  } else
  for (long p0=0; p0<np; p0+=g.bwork) {
    const int bp = (int)std::min<long>(g.bwork, np-p0);
    const double* phiB  = g.d_phi  + (g.phi_cached ? (size_t)p0*nb : 0);
    const double* gphiB = g.d_gphi ? g.d_gphi + (g.phi_cached ? (size_t)p0*nb : 0) : nullptr;
    const long    gstr  = (long)nb * (g.phi_cached ? (long)np : (long)bp);
    if (!g.phi_cached) { int TB=128, nblk=(bp+TB-1)/TB;
      k_collocate_g<<<nblk,TB>>>(bp,g.nshell,nb, gga?1:0, g.d_gx+p0,g.d_gy+p0,g.d_gz+p0,
         g.d_sh_am,g.d_sh_g0,g.d_sh_nc,g.d_sh_ao, g.d_sh_cx,g.d_sh_cy,g.d_sh_cz,g.d_sh_md2,
         g.d_ex,g.d_cc,g.d_pmd2, g.d_cartx,g.d_carty,g.d_cartz,g.maxcart, g.d_phi, g.d_gphi); }
    { const double a=1.0, b=0.0;
      cublasDgemm(g.cub, CUBLAS_OP_N, CUBLAS_OP_N, nb, bp, nb, &a, g.d_Dfold, nb, phiB, nb, &b, g.d_X, nb); }
    { int TB=256; unsigned nw=(unsigned)(((long)bp*32+TB-1)/TB);
      k_rho<<<nw,TB>>>(bp, nb, phiB, g.d_X, g.d_rho);
      if (gga) k_drho<<<nw,TB>>>(bp, nb, gphiB, gstr, g.d_X, g.d_drho); }
    { int TB=128, nblk=(bp+TB-1)/TB;
      k_func_gga<<<nblk,TB>>>(bp, g.func, g.d_rho, g.d_drho, g.d_gw+p0, g.d_eps, g.d_vrho, g.d_vsa, g.d_vsc); }
    { int TB=256; unsigned nz=(unsigned)(((long)nb*bp+TB-1)/TB);
      k_makeZ_gga<<<nz,TB>>>(bp, nb, gga?1:0, g.d_vrho, g.d_vsa, g.d_vsc, g.d_drho, phiB, gphiB, gstr, g.d_Z); }
    { const double a=1.0, b1=1.0;   // accumulate A = phi Z^T only; V = A + A^T after the loop
      cublasDgemm(g.cub, CUBLAS_OP_N, CUBLAS_OP_T, nb, nb, bp, &a, phiB, nb, g.d_Z, nb, &b1, g.d_Vfull, nb); }
    k_reduce<<<1,256>>>(bp, g.d_eps, g.d_rho, g.d_gw+p0, g.d_eexc, g.d_totele);
  }
  { dim3 TB(16,16), nb2((nb+15)/16,(nb+15)/16);
    k_sym_inplace<<<nb2,TB>>>(nb, g.d_Vfull); }   // V <- A + A^T
  double* d_fxc=nullptr;
  if (cudaMalloc(&d_fxc, ntri*8)!=cudaSuccess) return;
  {
    int TB=128, nblk=(nb+TB-1)/TB;
    k_pack_fold<<<nblk,TB>>>(nb, g.d_Vfull, g.d_bfnrm, d_fxc);
  }
  CKv(cudaDeviceSynchronize());
  cudaMemcpy(fxc, d_fxc, ntri*8, cudaMemcpyDeviceToHost);
  cudaMemcpy(eexc, g.d_eexc, 8, cudaMemcpyDeviceToHost);
  cudaMemcpy(totele, g.d_totele, 8, cudaMemcpyDeviceToHost);
  cudaFree(d_fxc);
  *info = 0;
}

// public entry: single-GPU -> one pass; OQP_MULTI_GPU with >=2 devices -> one
// pass per device on ITS OWN HOST THREAD (concurrent), each integrating its
// half of the grid, partials summed here on the host. Exact: XC is a plain sum
// over grid points, so splitting the grid splits nothing else.
extern "C" void routec_vxc(const double* d, double* fxc, const int* nbf,
                           const int* nfocks, double* eexc, double* totele, int* info) {
  *info = 1;
  static int mg = -1;
  if (mg < 0) {
    mg = 0;
    if (getenv("OQP_MULTI_GPU") && atoi(getenv("OQP_MULTI_GPU")) != 0) {
      int nd = 0;
      if (cudaGetDeviceCount(&nd) == cudaSuccess && nd >= 2) mg = 1;
      else fprintf(stderr, "[ownxc] OQP_MULTI_GPU set but <2 devices; single-GPU XC\n");
    }
    g.dev = 0; g.rank = 0; g.nranks = mg ? 2 : 1;
    if (mg) { g_xc1.dev = 1; g_xc1.rank = 1; g_xc1.nranks = 2;
              fprintf(stderr, "[ownxc] multi-GPU XC: grid halves on dev0/dev1\n"); }
  }
  if (!mg) { vxc_one(g, d, fxc, nbf, nfocks, eexc, totele, info); return; }
  const long ntri = (long)(*nbf)*((*nbf)+1)/2 * ((*nfocks == 2) ? 2 : 1);
  std::vector<double> fxc1(ntri, 0.0);
  double e1 = 0.0, t1 = 0.0; int info1 = 1;
  std::thread th([&]{ vxc_one(g_xc1, d, fxc1.data(), nbf, nfocks, &e1, &t1, &info1); });
  vxc_one(g, d, fxc, nbf, nfocks, eexc, totele, info);
  th.join();
  cudaSetDevice(0);                       // the SCF continues on device 0
  if (*info != 0 || info1 != 0) { if (*info == 0) *info = info1; return; }
  for (long t = 0; t < ntri; ++t) fxc[t] += fxc1[t];
  *eexc += e1; *totele += t1;
}
