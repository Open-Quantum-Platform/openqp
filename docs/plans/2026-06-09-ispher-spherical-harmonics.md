# Pure spherical-harmonic (ISPHER) basis support

Status: ENERGY PATH + ANALYTIC GRADIENTS (HF/DFT ground state, TDA, SF, and
focused MRSF checks) IMPLEMENTED. EKT smoke passes; ECP and analytic Hessian
have explicit spherical guards until their spherical kernels are implemented.
Gate: `[input] ispher` (runtime, default `true`) sets
`constants::HARMONIC_ACTIVE` before basis construction. With `ispher=false`,
`num_ao() == NUM_CART_BF` everywhere and every code path is Cartesian-only.

Branch: `feat/molecular-symmetry`. ~20 commits 4e3629d4..HEAD implement this.
Last validated state before §7 work: gate `.false.`, build green, tree clean.

## 0. RESUME HERE — build, gate, and validation recipe

The mechanics that took effort to work out; reuse them to continue.

Build (macOS, gfortran-15, ILP64 OpenBLAS):
```
wt=/Users/cheolhochoi/clone/openqp-symmetry
ob=/Users/cheolhochoi/clone/perf_int/openblas64        # prebuilt ILP64 OpenBLAS
cmake -S "$wt" -B "$wt/build" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran-15 \
  -DENABLE_OPENMP=ON -DENABLE_PYTHON=ON -DENABLE_MPI=OFF \
  -DUSE_LIBINT=OFF -DLINALG_LIB=OpenBLAS -DLINALG_LIB_INT64=ON \
  -DCMAKE_PREFIX_PATH="$ob" -DCMAKE_LIBRARY_PATH="$ob/lib" \
  -DBLAS_LIBRARIES="$ob/lib/libopenblas64.dylib" \
  -DLAPACK_LIBRARIES="$ob/lib/libopenblas64.dylib"
cmake --build "$wt/build" -j4
```
Note: the OpenTrustRegion ExternalProject only inherited a BLAS hint via
`-DBLA_VENDOR` or the NetLib path. This branch now forwards explicit
`BLAS_LIBRARIES`/`LAPACK_LIBRARIES` to OTR so the OpenBLAS configuration above
works.

Toggle the gate at runtime: set `[input] ispher=true/false` (default true).
Python calls `oqp_set_harmonic_active` before basis construction, so no rebuild
is needed to compare spherical vs Cartesian behavior.

Run a calculation (makeshift OPENQP_ROOT):
```
mkdir -p /tmp/oqp_root/include /tmp/oqp_root/lib /tmp/oqp_root/share
ln -sf "$wt/include/oqp.h"               /tmp/oqp_root/include/oqp.h
ln -sf "$wt/build/source/liboqp.dylib"   /tmp/oqp_root/lib/liboqp.dylib
ln -sfn "$wt/basis_sets"                 /tmp/oqp_root/share/basis_sets
export OPENQP_ROOT=/tmp/oqp_root PYTHONPATH="$wt/pyoqp"
export DYLD_LIBRARY_PATH="$wt/build/source:$ob/lib:/Users/cheolhochoi/clone/perf_int/OpenBLAS/exports:/opt/homebrew/opt/gcc/lib/gcc/current"
python3 -m oqp.pyoqp input.inp        # writes input.log / input.json
```
Reference: `pyscf` 2.13 is installed; use `cart=False` for spherical. Example
inputs live in `/tmp/oqp_validate/` (cc-pVDZ water/CH2O for HF/DFT/CIS/MRSF).

Validation snapshot (gate on, vs pyscf cart=False):
- RHF/cc-pVDZ water energy: -76.0266525950 (1e-10); gradient ~1e-8.
- PBE/cc-pVDZ energy 5e-6, gradient grid-level + own FD.
- RO-BHHLYP MRSF reference 8e-7; MRSF energy runs/converges.
- TDA-HF (CIS) state-1 gradient ~1e-5, dE/dx=0 restored; exc 9.0034 eV.
- planar CH2O MRSF gradient: dE/dx = 0.000000 all atoms (symmetry check).
- nbf reduced (cc-pVDZ water 24 vs 25 Cartesian).

Gotcha: MRSF finite-difference is unreliable because the triplet-ROHF
reference SCF has multiple solutions and jumps between them across displaced
geometries (seen as ROHF energy -117.9 vs -114.3). This is a reference-SCF
stability issue, NOT a gradient bug. Use closed-shell TDA/CIS vs pyscf, or
the planar-symmetry (dE/dx=0) test, to validate the response gradient.

## 1. What this delivers

Per-shell automatic selection of pure spherical (5d/7f/9g) vs Cartesian
(6d/10f/15g) AOs, matching each basis set's published convention, with the
AO dimension actually reduced (e.g. cc-pVDZ water: 24 spherical vs 25
Cartesian). Energies, the SCF guess, MOLDEN output, and log/JSON reporting
all work; gradients do not yet.

## 2. Architecture

- **Auto-selection**: `set_basis.py` reads BSE per-shell `function_type`
  (`gto_spherical`/`gto_cartesian`/`gto`) into a per-shell `harmonic` flag,
  threaded Python -> C struct `electron_shell` -> Fortran `basis_set%harmonic`
  (incl. MPI bcast and dump/load).
- **Dimension**: `constants::num_ao(l, harmonic)` returns `2l+1` for a pure
  shell when `HARMONIC_ACTIVE`, else `NUM_CART_BF(l)`. Drives `basis%naos`,
  `basis%nbf`, `basis%ao_offset`.
- **c2s transforms** (`integrals/cart2sph.F90`): matrices `B(l)` for d/f/g
  generated from and matching `pyoqp/oqp/library/symmetry.py`
  (`_solid_harmonic_coefficients`), CCA order m=-l..+l, `B S B^T = I`
  (self-test `c2s_selftest`, 1e-15). Entry points:
  - `cart2sph_eri`   — 4-index ERI shell-quartet block (full, unit-normalized in).
  - `cart2sph_mat`   — 2-index 1e block; folds `shells_pnrm2` (pure-power in);
                       `iandj` flag for same-shell lower-triangle packing.
  - `cart2sph_vec`   — 1-index AO value/derivative vector on the grid.

## 3. Normalization contract (important)

- 2e blocks reach the c2s hook AFTER backend normalization (genr22 bakes the
  sqrt(3)-type factors; libint/rys call `normalize_ints`), so they are
  unit-normalized Cartesian -> `cart2sph_eri` uses plain `B`.
- 1e blocks and grid AO values are pure-power Cartesian (`shells_pnrm2` is
  applied later, by `bas_norm_matrix` / the grid density build), so
  `cart2sph_mat`/`cart2sph_vec` FOLD `shells_pnrm2` along each transformed
  index, producing unit-normalized spherical components.
- `set_bfnorms` therefore writes `bfnrm = 1` for spherical components so the
  downstream normalization does not double-apply. Cartesian, spherical, and
  mixed shell pairs all end up unit-normalized exactly once.

## 4. Hook sites (implemented)

- `int1.F90` 1e: every `update_(triang|rectangular)_matrix` site (S, T,
  nuclear/Hcore, multipole, Lz). Same-shell sites pass `iandj`.
- `mod_shell_tools.F90`: `shell_t%harmonic`; shell-pair `%inao/%jnao =
  NUM_CART_BF(ang)` so the integral engine stays Cartesian while `nao`
  (placement) is spherical.
- `int2.F90` 2e: one hook after backend normalization covering rotspd, rys,
  and libint (libint result staged into `eri_data%ints` first).
- `basis_tools.F90` grid: `compAOv` (LDA) and `compAOvg` (GGA value+grad).

## 5. Validation (HARMONIC_ACTIVE temporarily on)

| Case | OpenQP | Reference (pyscf, cart=False) | Δ |
|------|--------|-------------------------------|---|
| RHF/cc-pVDZ water | -76.0266525950 | -76.0266525950 | 1e-10 |
| RHF/6-31G water (no d) | -75.98393862 | -75.98393863 | 1e-8 |
| PBE/cc-pVDZ water | -76.33357068 | -76.33357610 | 5e-6 (grid) |
| RO-BHHLYP ref (MRSF, CH2O/cc-pVDZ) | -114.3212163055 | -114.32121550 | 8e-7 (grid) |
| MRSF-TDDFT/BHHLYP/cc-pVDZ CH2O | runs, converges; states -2.951/1.024/5.029/5.952/7.309 eV | (integrals validated above) | — |

hcore and huckel guesses, MOLDEN ([5D7F9G] + m-order), and log/JSON
reporting (`nbf`, `nbf_cartesian`, `spherical_harmonics`) all verified.

## 6. Remaining work

### 6a. Gradients (the main gap)

The analytic gradient contracts per-shell-pair Cartesian derivative integrals
with a density slice `dens(ao_offset(ii):, ao_offset(jj):)`. With the gate on,
`dens` is spherical-dimensioned but the derivative routines iterate Cartesian
`inao/jnao`, so the slice is mis-shaped.

Recommended design — **Cartesian-effective density** (dual addressing):

1. Add a Cartesian addressing alongside the spherical one: `basis%ao_offset_cart`
   and `nbf_cart` (always Cartesian sizes). Cheap to compute next to `naos`.
2. Before `grd1`/`grd2`, build `D_cart(nbf_cart, nbf_cart)` from the spherical
   density `D_sph` block-by-block: `D_cart_block = M D_sph_block M^T` where
   `M(c,s) = B(c,s) * shells_pnrm2(c,l)` (the same pnrm-folded c2s used for 1e
   integrals; this is the transpose direction). Do the same for the
   energy-weighted density `W` used by the Pulay/overlap-derivative term.
3. Run the existing Cartesian gradient routines using `ao_offset_cart` +
   `D_cart`/`W_cart`. No change to the derivative kernels themselves.
4. XC gradient: hook `compAOvgg` (nDer=2) exactly like `compAOvg`
   (`cart2sph_vec` on the value, gradient, and the 6 second-derivative
   vectors), and confirm `dft_gridint_grad` consumes spherical AO arrays.

Validate by central finite difference of the (already-validated) energy:
`dE/dx ≈ [E(x+h) - E(x-h)] / 2h` per Cartesian coordinate, agreement to ~1e-6.
Gradients are all-or-nothing for validation (1e + 2e + XC together), so
implement and validate them as one unit.

### 6b. Smaller items

- `compAOvgg` (also needed for analytic Hessians).
- External-charge / Ewald 1e variants in `int1.F90` (PBC only) — same
  `cart2sph_mat` hook with `iandj`.
- MOLDEN spherical sign/normalization: order + [5D7F9G] markers are emitted
  and dimensions are correct; verify component signs against a MOLDEN viewer.

## 6c. Complete integral-consumer inventory

Sweep of every file that indexes per-shell AO blocks (`%ao_offset`/`%naos`)
or evaluates AOs (`aoval`), classified by status.

DONE / auto-correct (validated):
- `integrals/int1.F90`, `integrals/int2.F90` — 1e/2e integrals (hooked).
- `basis_tools.F90` (compAOv/compAOvg), `dftlib/dft_gridint.F90` — grid.
- `tdhf_lib.F90`, `tdhf_mrsf_lib.F90` — TDDFT/MRSF response: index by
  `naos`, so they auto-follow; MRSF energy validated.
- `guess.F90`, `modules/guess_minao.F90` — guesses (`naos`; huckel validated).
- `modules/electric_moments.F90`, `modules/get_basis_overlap.F90`,
  `modules/get_states_overlap.F90` — use hooked integrals + full-matrix
  density, auto-follow.
- `modules/population_analysis.F90` — fixed (naos, not NUM_CART_BF).

GROUND-STATE GRADIENT — DONE AND VALIDATED:
- `integrals/grd1.F90` — 1e gradient (prepare_grad_density, Cartesian-
  effective density + Cartesian offsets).
- `modules/hf_gradient.F90` — 2e gradient (d2a_cart/d2b_cart via build_cart;
  get_density branches on HARMONIC_ACTIVE; the 2-particle density factorizes
  so df1 is unchanged). RHF + UHF.
- `basis_tools.F90::compAOvgg` — grid 2nd derivatives (c2s on all 10 AO
  vectors); `dftlib/dft_gridint_grad.F90` auto-follows.
- Validation: water RHF/cc-pVDZ analytic gradient = pyscf to ~1e-8;
  PBE/cc-pVDZ = pyscf to ~1e-4 (grid) and OpenQP own finite difference to
  the FD noise floor.

EXCITED-STATE GRADIENTS — PARTIAL (TDA/SF validated; MRSF focused checks pass):
- Root cause: each of tdhf_gradient / tdhf_sf_gradient / tdhf_mrsf_gradient
  carries its OWN grd2 compute type whose get_density builds the response
  2-particle density on shell offsets, contracted with Cartesian derivative
  ERIs. None were hooked -> the spherical response gradient was wrong (broke
  dE/dx=0 for a planar molecule).
- Fix: build_cart on each type's densities (relaxed p, ground d, transition
  X+Y/X-Y for TDA; alpha/beta d,p + 7 spin-pair-coupling densities for MRSF;
  alpha/beta d,p + transition v for SF) and branch get_density. fock_deriv
  probes likewise. dft_gridint_tdxc_grad auto-follows (verified by the planar
  symmetry test).
- Validation: water TDA-HF (CIS)/cc-pVDZ state-1 gradient matches pyscf
  TDA nuc_grad to ~1e-5 with dE/dx=0 restored; planar CH2O MRSF/BHHLYP
  gradient has dE/dx=0.000000 for every atom.
- MRSF Cartesian baseline (2026-06-09): H2O/MRSF-HF/6-31G `ispher=false`,
  state 3, central FD of O-z with `h=1e-3 Ang` gives
  `-0.159089164 Ha/Bohr`; analytic reports `-0.159088970 Ha/Bohr`
  (`1.9e-7` difference). H2O/MRSF-BHHLYP/6-31G* `ispher=false` gives
  `7.5e-5` difference, consistent with the DFT/grid floor.
- MRSF spherical checks (2026-06-09): H2O/MRSF-HF/cc-pVDZ `ispher=true`,
  state 3, central FD of O-z with `h=1e-3 Ang` gives
  `-0.123737560 Ha/Bohr`; analytic reports `-0.123737500 Ha/Bohr`
  (`6.0e-8` difference). H2O/MRSF-BHHLYP/cc-pVDZ `ispher=true` gives
  `1.0e-4` difference on the default SG1 grid. Increasing to an unpruned
  `rad_npts=150`, `ang_npts=590`, `grid_ao_pruned=false` grid improves the
  BHHLYP difference to `5.1e-6`; a heavier `250/974` grid gives `4.6e-6`.
- MRSF f-shell checks (2026-06-09): H2O/MRSF-HF/cc-pVTZ `ispher=true`
  reports spherical AO type `5d/7f/9g` and gives O-z FD
  `-0.089715672 Ha/Bohr` vs analytic `-0.089715340 Ha/Bohr`
  (`3.3e-7` difference). H2O/MRSF-BHHLYP/cc-pVTZ with the unpruned `150/590`
  grid gives O-z FD `-0.115002670 Ha/Bohr` vs analytic
  `-0.115003790 Ha/Bohr` (`1.1e-6` difference).
- MRSF g-shell check (2026-06-09): H2O/MRSF-HF/cc-pVQZ `ispher=true`
  reports spherical AO type `5d/7f/9g` and gives O-z FD
  `-0.074139315 Ha/Bohr` vs analytic `-0.074138590 Ha/Bohr`
  (`7.2e-7` difference).
- Stress diagnostic: CH2O/MRSF-BHHLYP/cc-pVDZ `ispher=true`, target state 3,
  central FD of the state-3 energy for O-z displacement `h=1e-3 Ang` gives
  `-0.004349043 Ha/Bohr`, while the analytic gradient reports
  `+0.014762226 Ha/Bohr`. The same displaced-input harness is not a clean
  spherical-only regression because its `ispher=false` comparison also fails
  badly (`+5.600337987 Ha/Bohr` FD vs `+0.65990379 Ha/Bohr` analytic). Keep this
  CH2O case as a root/reference-state tracking diagnostic, not as the gate for
  spherical MRSF support.

GUARDED/PENDING — Hessian (2nd derivatives, der2 + compAOvgg):
- `integrals/grd1.F90` hess_* subroutines, `modules/hf_hessian.F90`,
  `modules/tdhf_hessian.F90`, `modules/tdhf_sf_hessian.F90`,
  `modules/hess1_selftest.F90`.
- `hf_hessian` and the Python analytic-Hessian driver now reject spherical AO
  dimensions with a clear message; use `[input] ispher=false` or numerical
  Hessian until spherical Hessian kernels are implemented.

GUARDED/VALIDATED — other:
- `ecp.F90` — `map_canonical` AO labeling assumes Cartesian; ECP basis sets
  with spherical shells need a spherical label map. Spherical ECP now aborts
  clearly with `set [input] ispher=false`; the Cartesian override was smoke-
  tested on HBr/aug-cc-pVDZ-PP.
- `modules/tdhf_mrsf_ekt.F90` — EKT works in MO/full-matrix space and was
  smoke-tested with H2O/MRSF-EKT-IP/cc-pVDZ `ispher=true` (24 spherical AOs).

NOT IN THIS BRANCH (handle when merged): SOC integrals and PCM/DDX solvent
have no source here (separate branches / external project).

BENIGN: `dftlib/dft_fuzzycell.F90` `num_cart_bf` is a grid-buffer upper
bound (Cartesian >= spherical); `scf_addons.F90` skeleton symmetrization is
the molecular-symmetry feature (off by default, indexes by `naos`).

## 7. Flipping the gate on by default

IN PROGRESS: `HARMONIC_ACTIVE` has been promoted to runtime `[input] ispher`,
default `true`, with `ispher=false` as the Cartesian override. The auto-
selection makes cc-pVDZ/def2/6-311G** spherical by default while 6-31G* stays
Cartesian because its BSE shells are Cartesian.
