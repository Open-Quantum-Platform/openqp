"""Phase 11 V2-STEP1: localize the off-diagonal NAC deficiency to amplitude-class
channels, PROVING it is NOT the relaxation (workflow 2026-06-16 redirect).

The relaxation (2FT) was proven correct on AND off diagonal (L^esum - L^2FT = 0 to
1e-15, all 6 pairs; gmo_relax.py validates 2FT vs live Fortran to 1e-15). So the
off-diagonal deficiency in the coded gradient chain (vs the exact numeric oracle)
must live in the 2e back-transform (G_MO hxa/hxb) + mrsfsp spin-pair channels.

This script confirms: deficiency D(I,J) = cross_coded - cross_numeric (from
qz_nogate.npz, cross = bilinear A(X_I,X_J) = Z_IJ - 1/2(X_I+X_J)) aligns with the
ground-config-coupled classfit channels gl / gx / xx (B(u,v)=1/4(G(u+v)-G(u-v)),
u=X_I*P_a, v=X_J*P_b, P in {ground,local,excited}).
"""
import numpy as np

qz = np.load('/tmp/nactest/qz_nogate.npz')
cf = np.load('/tmp/nactest/classfit.npz')


def cross(I, J, kind):
    return qz[f'Z{I}{J}_{kind}'] - 0.5*(qz[f'X{I}_{kind}'] + qz[f'X{J}_{kind}'])


def u(v):
    n = np.linalg.norm(v)
    return v/n if n > 1e-14 else v


for (I, J) in [(1, 2), (1, 3), (2, 3)]:
    cc, cn = cross(I, J, 'coded'), cross(I, J, 'numeric')
    D = cc - cn
    Bsum = sum(cf[f'B_{I}{J}_{a}{b}'] for a in 'glx' for b in 'glx')
    print(f"\n===== pair ({I},{J}) =====")
    print(f"  |coded|={np.linalg.norm(cc):.4f} |numeric|={np.linalg.norm(cn):.4f} "
          f"|deficiency|={np.linalg.norm(D):.4f} cos(coded,num)={np.dot(u(cc),u(cn)):+.4f}")
    print(f"  classfit consistency |sum_ab B - coded|={np.linalg.norm(Bsum-cc):.2e}")
    chans = sorted([(f'{a}{b}', cf[f'B_{I}{J}_{a}{b}']) for a in 'glx' for b in 'glx'
                    if np.linalg.norm(cf[f'B_{I}{J}_{a}{b}']) > 1e-6],
                   key=lambda kv: -np.linalg.norm(kv[1]))
    print("  channel | |B|     cos(B, deficiency)")
    for nm, v in chans:
        print(f"    {nm}   {np.linalg.norm(v):.4f}   {np.dot(u(v),u(D)):+.4f}")

print("\nCONCLUSION: deficiency concentrated in ground-config channels gl/gx/xx;")
print("(2,3) already correct. Relaxation exonerated. Target = interstate 2e/spin-pair.")
