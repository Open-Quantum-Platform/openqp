// M3 piece: real spherical harmonics Y_lm(theta,phi) on device, matching the
// host convention (no Condon-Shortley phase; Ivanic-Ruedenberg/Choi rotation
// convention). N_l^m = sqrt((2l+1)/(4pi) (l-|m|)!/(l+|m|)!).
#pragma once
#include "sph_eri/port.h"
#include <cmath>

namespace sph_eri {

// associated Legendre P_l^m(x), m>=0, WITHOUT the (-1)^m Condon-Shortley phase
template <typename T>
__host__ __device__ inline T assoc_legendre(int l, int m, T x) {
  T s = (T)sqrt(fmax(0.0, 1.0 - (double)x * (double)x));
  T pmm = (T)1;
  for (int k = 1; k <= m; ++k) pmm *= (T)(2 * k - 1) * s;   // (2m-1)!! s^m, no (-1)^m
  if (l == m) return pmm;
  T pmmp1 = x * (T)(2 * m + 1) * pmm;
  if (l == m + 1) return pmmp1;
  T pll = (T)0;
  for (int ll = m + 2; ll <= l; ++ll) {
    pll = (x * (T)(2 * ll - 1) * pmmp1 - (T)(ll + m - 1) * pmm) / (T)(ll - m);
    pmm = pmmp1; pmmp1 = pll;
  }
  return pll;
}

template <typename T>
__host__ __device__ inline T real_sph_harm(int l, int m, T theta, T phi) {
  T x = (T)cos((double)theta);
  int am = m < 0 ? -m : m;
  T plm = assoc_legendre<T>(l, am, x);
  double N = sqrt((2.0 * l + 1.0) / (4.0 * M_PI)
                  * tgamma((double)(l - am + 1)) / tgamma((double)(l + am + 1)));
  if (m == 0) return (T)N * plm;
  T ph = (T)(m > 0 ? cos((double)(m * phi)) : sin((double)(am * phi)));
  return (T)(sqrt(2.0) * N) * plm * ph;
}

} // namespace sph_eri
