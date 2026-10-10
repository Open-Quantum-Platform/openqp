// M2: McMurchie-Davidson Hermite/Gaussian-product coefficients E^{la,lb}_t on
// device. Stable upward recursion (no cancellation). Fills out[t], t=0..la+lb.
//   E^{00}_0 = exp(-mu Q^2),  mu=ab/p, p=a+b
//   j=0:  E^{i,0}_t = 1/(2p) E^{i-1,0}_{t-1} - (mu Q/a) E^{i-1,0}_t + (t+1) E^{i-1,0}_{t+1}
//   else: E^{i,j}_t = 1/(2p) E^{i,j-1}_{t-1} + (mu Q/b) E^{i,j-1}_t + (t+1) E^{i,j-1}_{t+1}
#pragma once
#include "sph_eri/port.h"

namespace sph_eri {

// LA,LB <= 6 -> t up to 12; flat DP buffers sized for that.
template <typename T>
__host__ __device__ inline void hermite_E_axis(int la, int lb, T Q, T a, T b, T* out) {
  T p = a + b, mu = a * b / p;
  const int TM = 14;            // need index t+1 up to 13 (la+lb <= 12)
  T cur[TM], nxt[TM];
  for (int t = 0; t < TM; ++t) cur[t] = (T)0;
  cur[0] = (T)exp(-(double)(mu * Q * Q));      // E^{0,0}_0
  // build i = 1..la at j=0
  for (int i = 1; i <= la; ++i) {
    int tmax = i;                              // t <= i + 0
    for (int t = 0; t <= tmax; ++t) {
      T v = (T)0;
      if (t-1 >= 0)      v += (T)(1.0/(2.0*(double)p)) * cur[t-1];
      v += -(mu*Q/a) * cur[t];
      v += (T)(t+1) * cur[t+1];
      nxt[t] = v;
    }
    for (int t = tmax+1; t < TM; ++t) nxt[t] = (T)0;
    for (int t = 0; t < TM; ++t) cur[t] = nxt[t];
  }
  // build j = 1..lb at i=la
  for (int j = 1; j <= lb; ++j) {
    int tmax = la + j;
    for (int t = 0; t <= tmax; ++t) {
      T v = (T)0;
      if (t-1 >= 0)      v += (T)(1.0/(2.0*(double)p)) * cur[t-1];
      v += (mu*Q/b) * cur[t];
      v += (T)(t+1) * cur[t+1];
      nxt[t] = v;
    }
    for (int t = tmax+1; t < TM; ++t) nxt[t] = (T)0;
    for (int t = 0; t < TM; ++t) cur[t] = nxt[t];
  }
  for (int t = 0; t <= la+lb; ++t) out[t] = cur[t];
}

} // namespace sph_eri
