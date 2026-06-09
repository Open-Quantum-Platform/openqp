# Pure spherical-harmonic (ISPHER) basis support

Status: ENERGY PATH IMPLEMENTED AND VALIDATED. Gradients/Hessian pending.
Gate: `constants::HARMONIC_ACTIVE` (compile-time, default `.false.`). With the
gate off, `num_ao() == NUM_CART_BF` everywhere and every code path is
byte-identical to the Cartesian-only behavior.

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

## 7. Flipping the gate on by default

Only after 6a lands and finite-difference-validates. At that point
`HARMONIC_ACTIVE` can default `.true.` (or, better, be promoted to a runtime
input keyword `ispher` with the BSE default per basis and an explicit
override), and the auto-selection makes cc-pVDZ/def2/6-311G** correct by
default while 6-31G* stays Cartesian.
