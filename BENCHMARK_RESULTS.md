# Numerical NAC factor-2 / sign fix — benchmark results

Validation of the fix in `pyoqp/oqp/library/single_point.py` (`NACME.nacme`,
`NAC.numerical_nac`). Two independent levels: an exact convention proof that
needs no backend, and a live MRSF run through the compiled Fortran backend.

Definitions (real adiabatic states): derivative coupling
`d_ij = <Psi_i|grad Psi_j>` is **antisymmetric**; interstate coupling
`h_ij = (E_j - E_i) d_ij` is **symmetric**. The overlap-based estimator is
`d = (S - S^T)/(2*dt)` (Hammes-Schiffer-Tully). Per displacement,
`current` denotes the OLD `numerical_nac` output and `corrected` the FIXED one:

```
current   = [(S+ - S+^T) - (S- - S-^T)] / (2*Delta)   = 2*d   (OLD, factor-2 bug)
corrected = [(S+ - S+^T) - (S- - S-^T)] / (4*Delta)   =   d   (FIXED, HST 1/2)
```

## 1. Exact convention proof — `fd_nac_logic_check.py` (no backend)

Synthetic overlaps built from a known antisymmetric `d` plus O(Delta^2)
curvature noise; exercises the exact assembly logic of both code paths.

```
[PASS] nacme dc == true d                 rel=1.33e-17
[PASS] numerical_nac dc == true d         rel=1.33e-17
[PASS] d antisymmetric (both paths)
[PASS] h symmetric (both paths)
[PASS] nacme h == (E_j-E_i) d             rel=1.15e-17
[PASS] numerical_nac h == (E_j-E_i) d     rel=1.15e-17
[PASS] cross-path dc agree                rel=0.00e+00
[PASS] cross-path h agree (sign+mag)      rel=0.00e+00
[PASS] OLD nacme dc == 2 x true d         ratio=2.0000
[PASS] OLD numerical_nac dc == 2 x true d ratio=2.0000
[PASS] OLD numerical_nac h == -2 h_std (factor 2 + sign flip)
[PASS] OLD nacme h == +2 h_std (factor 2, std sign)

pair 1->2:  current(old)/hand = 2.0000   corrected(new)/hand = 1.0000
```

Reproduces the handoff identity exactly: OLD = 2x, FIXED = 1x; both paths now
agree in magnitude **and** sign.

## 2. Live MRSF backend — `fd_nac_benchmark.py`

H2O (C2v, no spatial degeneracy), MRSF/BHHLYP/6-31G, ROHF triplet reference,
O atom displaced along z, displacement scan Delta = 0.02 ... 0.0025 Angstrom.
Run single-threaded (`OMP_NUM_THREADS=1`) for bit-reproducible overlaps.

**The factor-2 (every element, every Delta):**

| Delta (A) | current/corrected (min) | (max) |
|---|---|---|
| 0.0200 | 2.000000 | 2.000000 |
| 0.0100 | 2.000000 | 2.000000 |
| 0.0050 | 2.000000 | 2.000000 |
| 0.0025 | 2.000000 | 2.000000 |

**`corrected` vs the true HST `d` (= antisymmetric part of the FD), clean pairs:**

| Delta (A) | pair | gap (Ha) | corrected | corr/h_anti | cur/h_anti | sym/anti |
|---|---|---|---|---|---|---|
| 0.0025 | S1->S4 | 0.3937 | -0.052552 | **1.0000** | **2.0000** | 0.556 |
| 0.0025 | S1->S6 | 0.8104 |  0.003511 | **1.0000** | **2.0000** | 0.390 |
| 0.0025 | S4->S6 | 0.4167 | -0.023359 | **1.0000** | **2.0000** | 0.048 |

(`corr/h_anti` and `cur/h_anti` are 1.0000 / 2.0000 at every Delta, not just the
smallest.)

**Structural antisymmetry of the output (defining property of `d`):**

```
Delta=0.0200..0.0025:  ||cur+cur^T||/||cur|| = 0.00e+00   ||corr+corr^T||/||corr|| = 0.00e+00
```

### Note on the `sym/anti` column

The *raw* MRSF state overlap is not antisymmetric off-diagonal: e.g. for S1->S4,
`S+[1,4] = -0.000768` but `S+[4,1] = +0.000231`. It carries a symmetric component
(basis-following / phase) of comparable size (`sym/anti` up to 0.56). The
derivative coupling is **defined** as the antisymmetrized overlap, and both code
paths antisymmetrize (`S - S^T`), correctly discarding this symmetric part — as
the HST convention requires. This is why a naive `(S+ - S-)/(2*Delta)` total
difference is **not** the right reference; the antisymmetric part is, and against
it the fixed output matches exactly (`corr/h_anti = 1.0000`).

## 3. Reproduce

```bash
# build the backend once
cmake -B build -G Ninja -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=OFF
ninja -C build install

export OPENQP_ROOT=$PWD
export LD_LIBRARY_PATH=$OPENQP_ROOT/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$OPENQP_ROOT/pyoqp:$PYTHONPATH

python3 fd_nac_logic_check.py              # exact convention proof (no backend)
OMP_NUM_THREADS=1 python3 fd_nac_benchmark.py   # live H2O MRSF backend
```

## 4. Scope / not validated here

- The fix changes the **exported** DCME/NACME (and the `back_door` API) consumed
  by the external PyRAI2MD NAMD driver, halving the exported time-derivative
  coupling. Confirm/de-compensate PyRAI2MD's denominator before this affects
  production dynamics.
- An absolute cross-check of `h_ij` against an independent analytic NAC is not
  possible in OpenQP today (analytic MRSF NAC is not implemented); it is the
  recommended next validation once that backend exists.
