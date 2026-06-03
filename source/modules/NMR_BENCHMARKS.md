# RHF/DFT NMR Shielding — Benchmark & Regression Report

**Scope:** benchmark and regression of the *existing, validated* RHF/closed-shell-DFT
common-gauge-origin (CGO) NMR shielding feature (`scf_prop=nmr`). This is **not**
Gate 2b-B and introduces no new theory; it exercises and documents the shipped
ground-state code (`source/modules/nmr_shielding.F90`, `int1.F90`) against the
committed independent common-gauge reference.

## 1. Build & environment tested

| Item | Value |
|------|-------|
| Branch | `claude/cool-ritchie-vaMB6` |
| Commit | `421d8fb` (`docs(nmr): add Gate 2b-A formal derivation handoff`) |
| Host | Linux 6.18.5 x86_64 (managed remote container) |
| Build type | `Release` |
| Fortran / C / C++ | GNU 13.3.0 (`gfortran`/`gcc`/`g++`) |
| BLAS/LAPACK | OpenBLAS 0.3.26 (pthread), 32-bit int (`LINALG_LIB_INT64=OFF`) |
| 2e integral backend | **native Rys** (`USE_LIBINT=OFF`) — see §7 reproducibility note |
| CMake / Ninja | 3.28.3 / 1.x |
| Python | 3.11.15, NumPy 2.4.6, SciPy 1.17.x |

Configure / build / install (exact commands):

```bash
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_INSTALL_PREFIX=$(pwd) \
  -DENABLE_PYTHON=ON -DENABLE_OPENMP=ON -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=ON \
  -DLINALG_LIB_INT64=OFF -DUSE_LIBINT=OFF
ninja -C build && ninja -C build install
export OPENQP_ROOT=$(pwd) PYTHONPATH=$(pwd)/pyoqp
```

## 2. Regression tests — all pass

NMR tests are Python `unittest` modules (not CTest-registered; `ctest -N` lists
none). Run:

```bash
export OPENQP_ROOT=$(pwd) PYTHONPATH=$(pwd)/pyoqp
python3 -m unittest tests.test_nmr_shielding tests.test_nmr_coupled -v
```

Result: **8 / 8 passed** in ~2.0 s wall.

| Test | Checks | Status |
|------|--------|--------|
| `test_nmr_shielding.NMRShieldingTests.test_nmr_shielding` | PSO `max\|diag\|`=0, `max\|A+Aᵀ\|`=0; O/H dia+uncoupled-para vs the independent reference | **ok** |
| `test_nmr_coupled.…test_gate0_density_antisymmetric` | `max\|P^B+P^Bᵀ\|` < 1e-9 (all functionals) | **ok** |
| `…test_gate1_coulomb_response_vanishes` | `\|\|J(P^B)\|\|` < 1e-9 | **ok** |
| `…test_gate2_exchange_response_nonzero` | `\|\|K(P^B)\|\|` > 1e-3 | **ok** |
| `…test_gate3_pure_dft_coupled_equals_uncoupled` | PBE coupled ≡ uncoupled (Δ < 1e-6) | **ok** |
| `…test_gate4_hf_matches_oracle_absolutely` | HF coupled/uncoupled para vs oracle (≤ 5e-2 ppm) | **ok** |
| `…test_gate4_coupling_delta_matches_oracle` | Δ=coupled−uncoupled vs oracle (≤ 0.12 ppm) | **ok** |
| `…test_gate6_delta_scales_with_exact_exchange` | `\|Δ(O)\|`: HF > BHHLYP > PBE0 > PBE=0 | **ok** |

## 3. Benchmark matrix

System: **H₂O**, basis **STO-3G**, closed-shell RHF reference, gauge origin =
center of mass `(0, 0, −0.19886422)` Bohr (same geometry as the oracle). Methods:
HF and three functionals spanning the exact-exchange fraction `c_x`.

```bash
# per case (HF has no [dftgrid]/functional block):
python3 -m oqp.pyoqp <case>.inp      # writes <case>.log with the CGO NMR table
```

| Case | `c_x` | Functional keyword | Grid |
|------|------|--------------------|------|
| HF | 1.00 | — | — |
| BHHLYP | 0.50 | `bhhlyp` | becke |
| PBE0 | 0.25 | `pbe0` | becke |
| PBE | 0.00 | `pbe` | becke |

## 4. Isotropic shielding results (ppm) vs independent common-gauge reference

Oracle: `tests/fixtures/nmr/cgo_reference.json` (independent NMR reference). OQP values
from the `.log` NMR table; `absErr` = |OQP − oracle|. O = oxygen, H = either
(equivalent) hydrogen.

**Oxygen:**

| Case | σ_dia (OQP/oracle) | σ_para^unc (OQP/oracle) | σ_para^cpl (OQP/oracle) | σ_tot^cpl (OQP/oracle) |
|------|--------------------|-------------------------|-------------------------|------------------------|
| HF | 411.4176 / 411.4176 | −113.6305 / −113.6305 | −230.6295 / −230.6294 | 180.7881 / 180.7882 |
| BHHLYP | 411.4664 / 411.4683 | −156.7630 / −156.7420 | −242.1108 / −242.0640 | 169.3556 / 169.4043 |
| PBE0 | 411.4438 / 411.4463 | −194.2959 / −194.2395 | −248.9695 / −248.8778 | 162.4743 / 162.5686 |
| PBE | 411.3452 / 411.3488 | −258.4350 / −258.3205 | −258.4350 / −258.3205 | 152.9101 / 153.0282 |

**Hydrogen:**

| Case | σ_dia (OQP/oracle) | σ_para^unc (OQP/oracle) | σ_para^cpl (OQP/oracle) | σ_tot^cpl (OQP/oracle) |
|------|--------------------|-------------------------|-------------------------|------------------------|
| HF | 28.0618 / 28.0618 | 1.7849 / 1.7849 | 3.5059 / 3.5059 | 31.5677 / 31.5677 |
| BHHLYP | 28.0175 / 28.0157 | 2.3889 / 2.3886 | 3.5808 / 3.5802 | 31.5983 / 31.5959 |
| PBE0 | 28.0581 / 28.0558 | 2.8754 / 2.8747 | 3.6073 / 3.6061 | 31.6654 / 31.6619 |
| PBE | 28.1952 / 28.1920 | 3.6522 / 3.6507 | 3.6522 / 3.6507 | 31.8474 / 31.8427 |

### Max absolute error vs oracle (over all reported components, O & H)

| Case | max\|abs err\| (ppm) | max\|Δ-coupling err\| (ppm) | Note |
|------|---------------------|------------------------------|------|
| HF | **0.00008** | 0.00008 | grid-free → matches oracle essentially exactly |
| BHHLYP | 0.04874 | 0.02582 | abs offset = cross-code SCF/grid; coupling reproduced |
| PBE0 | 0.09431 | 0.03526 | abs offset = cross-code SCF/grid; coupling reproduced |
| PBE | 0.11806 | **0.00000** | coupled ≡ uncoupled exactly; Δ-coupling = 0 |

`Δ-coupling = σ_para^coupled − σ_para^uncoupled`; comparing OQP's Δ to the oracle's
Δ removes the cross-code DFT-SCF/grid difference in the *absolute* numbers and
isolates the exact-exchange coupling machinery. HF (no grid) agrees with the oracle
to ~1e-4 ppm; the DFT absolute offsets (≤ 0.12 ppm) are SCF/grid cross-code
differences (the diamagnetic part alone differs only ≤ 0.004 ppm), while the
coupling Δ matches to ≤ 0.035 ppm and the per-atom errors stay within the test
tolerances (`TOL_HF_ABS=5e-2`, `TOL_DELTA=0.12`).

### Phase-0 gate diagnostics (from the logs)

| Case | gate0 `max\|P^B+P^Bᵀ\|` | gate1 `\|\|J(P^B)\|\|` | gate2 `\|\|K(P^B)\|\|` | PSO `max\|diag\|`, `max\|A+Aᵀ\|` |
|------|-------------------------|------------------------|------------------------|----------------------------------|
| HF | 5.3e-17 | 2.0e-16 | 0.5417 | 0, 0 |
| BHHLYP | 1.5e-17 | 1.7e-17 | 0.5611 | 0, 0 |
| PBE0 | 1.2e-16 | 4.2e-16 | 0.5737 | 0, 0 |
| PBE | 1.1e-16 | 7.2e-17 | 0.5904 | 0, 0 |

`P^B` antisymmetric and `J(P^B)=0` to machine precision; exchange response nonzero;
PSO exactly antisymmetric — consistent with `NMR_MAGNETIC_SYMMETRY_NOTE.md`.

## 5. Runtime

Single-point energy + NMR property, H₂O/STO-3G, wall time (whole driver
subprocess, cold per run; container CPU, OpenMP on). Tiny system → dominated by
process/SCF startup, not the NMR step.

| Case | Wall (s) |
|------|----------|
| HF | 0.38 |
| BHHLYP | 0.42 |
| PBE0 | 0.42 |
| PBE | 0.40 |
| **Test suite (8 tests, both modules)** | **~2.0** |

## 6. Limitations

- **CGO only** — common gauge origin at the center of mass; results are
  gauge-origin dependent (no GIAO/London orbitals). STO-3G numbers are for
  cross-code *consistency*, not basis-converged physical shieldings.
- **Closed-shell RHF / pure-or-hybrid DFT** references only; UHF/ROHF and MRSF are
  out of scope here.
- DFT absolute shieldings carry a cross-code SCF/grid offset (≤ ~0.12 ppm at
  STO-3G) vs the independent reference; the coupling Δ is the SCF-robust cross-code metric.
- Semi-local `f_xc[P^B]=0` and `J(P^B)=0` are exact within the conventional
  (non-current-density) functional class; current-DFT terms are not included.

## 7. Reproducibility notes

- **`USE_LIBINT=OFF` in this run.** The default production build fetches a custom
  libint (`oqp-libint`) from `qchemlab.knu.ac.kr`, which is **blocked by this
  container's network policy** (HTTP 403 `host_not_allowed`). The build therefore
  used OpenQP's **native Rys** 2e backend. This is benign for these benchmarks: the
  1e NMR operators (`L_O`, PSO, diamagnetic) are native Rys/Gauss–Hermite
  regardless of backend, and the only 2e coupling is the exact-exchange image of
  `P^B`; HF reproduced the independent libint-based reference to ~1e-4 ppm, confirming the
  Rys 2e path is equivalent here. On an unrestricted network, build with
  `-DUSE_LIBINT=ON` for the production integral path.
- Toolchain (`gfortran`, OpenBLAS) and Python deps (`numpy`, `scipy`, `cffi`,
  `basis_set_exchange`, `libdlfind`) were installed into the container; `geometric`
  is **not** required for single-point NMR (it is imported lazily only on the
  geometry-optimization path).
- The oracle JSON is committed; do not recompute ad hoc (regeneration needs an independent NMR reference implementation).
- Driver invocation: `python3 -m oqp.pyoqp <input>.inp` with `OPENQP_ROOT` and
  `PYTHONPATH=$OPENQP_ROOT/pyoqp` set; the NMR table is written to `<input>.log`.

## 8. Cross-references

- `source/modules/NMR_SHIELDING_STATUS.md` — feature status & validation.
- `source/modules/NMR_MAGNETIC_SYMMETRY_NOTE.md` — `J=0`, `f_xc=0`, exchange-only.
- `tests/test_nmr_shielding.py`, `tests/test_nmr_coupled.py` — the regression tests
  exercised above; `tests/fixtures/nmr/cgo_reference.json` — the oracle.
- `examples/NMR/H2O_RHF-NMR.inp` — example input.
