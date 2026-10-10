// UKS-capable validator: feed a density through routec_vxc at nfocks=1 or 2.
// For nf=2 we build a SPIN-SYMMETRIC guess alpha=beta=D/2, whose XC energy and
// alpha-channel Vxc must equal the closed-shell (RKS) reference in ref.bin.
// Compare a dense run vs a sparse run (OQP_OWNXC_DENSE / OQP_OWNXC_SPARSE) to
// validate the active-AO UKS path.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <dlfcn.h>

typedef void (*vxc_fn)(const double*, double*, const int*, const int*, double*, double*, int*);
static long long rd_i8(FILE* f){ long long v; size_t r=fread(&v,8,1,f); (void)r; return v; }

int main(int argc, char** argv) {
  if (argc < 3) { printf("usage: %s <so> <dumpdir> [nf=1]\n", argv[0]); return 1; }
  std::string so = argv[1], D = argv[2];
  int nf = (argc>3) ? atoi(argv[3]) : 1;

  FILE* fb = fopen((D+"/basis.bin").c_str(),"rb");
  rd_i8(fb); long nbf = (long)rd_i8(fb); fclose(fb);
  long ntri = nbf*(nbf+1)/2;

  std::vector<double> dpack(ntri);
  { FILE* f=fopen((D+"/dens.bin").c_str(),"rb"); long n=(long)rd_i8(f);
    size_t r=fread(dpack.data(),8,n,f); (void)r; fclose(f); }

  std::vector<double> fxc_ref(ntri); double eexc_ref=0, totele_ref=0;
  { FILE* f=fopen((D+"/ref.bin").c_str(),"rb"); long n=(long)rd_i8(f);
    size_t r=fread(fxc_ref.data(),8,n,f); r=fread(&eexc_ref,8,1,f); r=fread(&totele_ref,8,1,f); (void)r; fclose(f); }

  void* h = dlopen(so.c_str(), RTLD_NOW);
  if (!h) { printf("dlopen fail: %s\n", dlerror()); return 1; }
  vxc_fn fn = (vxc_fn)dlsym(h, "routec_vxc");
  if (!fn) { printf("dlsym fail\n"); return 1; }

  int nbf_c=(int)nbf, info=99;
  double eexc=0, totele=0;
  std::vector<double> din, fxc;
  if (nf==1) { din=dpack; fxc.assign(ntri,0.0); }
  else {       // alpha = beta = D/2  -> [a_tri | b_tri]
    din.assign(2*ntri,0.0); fxc.assign(2*ntri,0.0);
    for (long t=0;t<ntri;t++){ din[t]=0.5*dpack[t]; din[ntri+t]=0.5*dpack[t]; }
  }
  fn(din.data(), fxc.data(), &nbf_c, &nf, &eexc, &totele, &info);

  printf("nf=%d info=%d\n", nf, info);
  printf("Exc:    ours = %.12f   ref = %.12f   diff = %.3e\n", eexc, eexc_ref, fabs(eexc-eexc_ref));
  printf("totele: ours = %.10f   ref = %.10f   diff = %.3e\n", totele, totele_ref, fabs(totele-totele_ref));

  // alpha channel Vxc vs ref (== RKS Vxc for the spin-symmetric density)
  double maxabs=0,sum2=0,refn=0,abdiff=0;
  for (long t=0;t<ntri;t++){
    double dd=fabs(fxc[t]-fxc_ref[t]); if(dd>maxabs)maxabs=dd; sum2+=dd*dd; refn+=fxc_ref[t]*fxc_ref[t];
    if (nf==2) abdiff += fabs(fxc[t]-fxc[ntri+t]);   // alpha vs beta (should be 0)
  }
  printf("Vxc(alpha): max|ours-ref|=%.3e  rel||diff||=%.3e\n", maxabs, sqrt(sum2)/sqrt(refn));
  if (nf==2) printf("Vxc: sum|alpha-beta| = %.3e (spin-symmetric -> must be ~0)\n", abdiff);
  return 0;
}
