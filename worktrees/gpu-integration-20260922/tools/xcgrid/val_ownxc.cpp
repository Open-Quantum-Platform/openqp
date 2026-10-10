// Standalone validator: dlopen libroutec_ownxc_vxc.so, feed the lane-dumped
// packed density (dens.bin) through routec_vxc, compare Exc + Vxc vs the native
// reference (ref.bin). All counts are int64 (lane is -fdefault-integer-8).
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>
#include <dlfcn.h>

typedef void (*vxc_fn)(const double*, double*, const int*, const int*, double*, double*, int*);

static long long rd_i8(FILE* f){ long long v; fread(&v,8,1,f); return v; }

int main(int argc, char** argv) {
  if (argc < 3) { printf("usage: %s <so> <dumpdir>\n", argv[0]); return 1; }
  std::string so = argv[1], D = argv[2];

  // read nbf from basis.bin
  FILE* fb = fopen((D+"/basis.bin").c_str(),"rb");
  rd_i8(fb); long nbf = (long)rd_i8(fb); fclose(fb);
  long ntri = nbf*(nbf+1)/2;

  // read dens.bin
  std::vector<double> dpack(ntri);
  { FILE* f=fopen((D+"/dens.bin").c_str(),"rb"); long n=(long)rd_i8(f);
    fread(dpack.data(),8,n,f); fclose(f); }

  // read ref.bin (native Vxc, Exc, totele)
  std::vector<double> fxc_ref(ntri);
  double eexc_ref, totele_ref;
  { FILE* f=fopen((D+"/ref.bin").c_str(),"rb"); long n=(long)rd_i8(f);
    fread(fxc_ref.data(),8,n,f); fread(&eexc_ref,8,1,f); fread(&totele_ref,8,1,f); fclose(f); }

  void* h = dlopen(so.c_str(), RTLD_NOW);
  if (!h) { printf("dlopen fail: %s\n", dlerror()); return 1; }
  vxc_fn fn = (vxc_fn)dlsym(h, "routec_vxc");
  if (!fn) { printf("dlsym fail\n"); return 1; }

  std::vector<double> fxc(ntri, 0.0);
  int nbf_c = (int)nbf, nf = 1, info = 99;
  double eexc=0, totele=0;
  fn(dpack.data(), fxc.data(), &nbf_c, &nf, &eexc, &totele, &info);

  printf("info = %d  (0 = ok)\n", info);
  printf("Exc:    ours = %.12f   native = %.12f   abs diff = %.3e\n",
         eexc, eexc_ref, fabs(eexc-eexc_ref));
  printf("totele: ours = %.10f   native = %.10f   abs diff = %.3e\n",
         totele, totele_ref, fabs(totele-totele_ref));

  // Vxc matrix comparison
  double maxabs=0, sum2=0, refnorm=0;
  long imax=0;
  for (long t=0;t<ntri;t++){
    double dd = fabs(fxc[t]-fxc_ref[t]);
    if (dd>maxabs){maxabs=dd;imax=t;}
    sum2 += dd*dd; refnorm += fxc_ref[t]*fxc_ref[t];
  }
  printf("Vxc: max|ours-native| = %.3e   ||diff||_2 = %.3e   (||native||_2 = %.3e)\n",
         maxabs, sqrt(sum2), sqrt(refnorm));
  printf("Vxc: rel ||diff||/||native|| = %.3e\n", sqrt(sum2)/sqrt(refnorm));
  printf("at max element t=%ld: ours=%.10e native=%.10e\n", imax, fxc[imax], fxc_ref[imax]);
  return 0;
}
