# NMR Nuclear Magnetic Shielding — Implementation Status

**Branch:** `feat/dft-nmr`
**Scope (v1):** Common Gauge Origin (CGO), closed-shell RHF / pure DFT, *uncoupled*
paramagnetic term (exact CPKS for pure functionals; HF-uncoupled approximation
for Hartree–Fock). Implemented entirely in the native OpenQP Fortran core — no
PySCF dependency at run time (PySCF was used only as an external validation
oracle during development).

Requested via the input file:

```ini
[properties]
scf_prop=nmr
```

Output: per-atom isotropic shielding `sigma = sigma_dia + sigma_para` (ppm),
gauge origin defaulting to the molecular center of mass.

## Components

| Piece | File | Status |
|-------|------|--------|
| Angular-momentum integrals `(r-O)x∇` | `mod_1e_primitives.F90::comp_amom_int1_prim`, `int1.F90::angular_momentum_integrals` | Validated |
| Diamagnetic integrals `(r-O)·(r-R_N)/r_N³` | `mod_1e_primitives.F90::comp_nmr_dia_int1_prim`, `int1.F90::nmr_dia_shielding` | Validated |
| PSO integrals `(r-R_N)x∇/r_N³` | `mod_1e_primitives.F90::comp_pso_int1_prim`, `int1.F90::pso_integrals` | **Preliminary** |
| Shielding assembly + output | `modules/nmr_shielding.F90` | Validated (dia) / preliminary (para) |
| Registration | `include/oqp.h`, `runfunc.py`, `input_checker.py` | Done |

## What has been validated

All numbers below are **H2O / STO-3G, RHF**, gauge origin = center of mass
`(0, 0, -0.19886422)` Bohr, compared against PySCF (`pyscf-properties`,
`pyscf.prop.nmr.rhf`, common-gauge, unit factor `alpha²·1e6`).

1. **Angular-momentum integrals** reproduce PySCF `int1e_cg_irxp` to ~8
   significant figures (Frobenius norms per component).
2. **Diamagnetic shielding** reproduces PySCF `dia()` to ~6 significant figures
   for *all* atoms:

   | Atom | OpenQP σ_dia (ppm) | PySCF σ_dia (ppm) |
   |------|--------------------|-------------------|
   | O    | 411.417615         | 411.417599        |
   | H    | 28.061826          | 28.061820         |

3. **Paramagnetic assembly** (MO transform of the orbital-Zeeman and PSO
   operators, occupied–virtual sum-over-states, `2·alpha²` prefactor) reproduces
   the PySCF *uncoupled* reference for the **oxygen** atom:

   | Atom | OpenQP σ_para | reference σ_para |
   |------|---------------|------------------|
   | O    | -113.630502   | -113.630         |
   | H    | 1.236891      | 1.784938 (✗)     |

## What remains unresolved

The **PSO integral is incorrect for nuclei that are not centered on a basis
function** (e.g. hydrogen in H2O). The computed matrix acquires a small
*non-antisymmetric* component (a nonzero diagonal of order 0.04 a.u.), whereas
the PSO operator `L_N/r_N³` is anti-Hermitian and its real representation must be
exactly antisymmetric. Consequences:

- Atom-centered nuclei (oxygen) are correct.
- Off-center nuclei (hydrogen) give the wrong paramagnetic shielding
  (1.24 vs reference 1.78 ppm → total 29.30 vs 29.85 ppm).

### Diagnosis (to avoid re-deriving)

- The `(r-R_N)/r_N³` **field factor is confirmed correct** — the OpenQP/PySCF
  ratio for the field integral is exactly 2.0 (an understood `fac` normalization),
  i.e. the field itself matches to full precision.
- The `alpha²` **prefactors are confirmed correct** (diamagnetic matches the
  reference, which shares the same scale).
- More Rys roots do **not** change the result (the quadrature is converged).
- The principled `½(ket-derivative − bra-derivative)` antisymmetrization does
  **not** remove the contamination.
- Conclusion: the defect is in the **combination of the field factor with the
  ket-derivative on the Rys (`DQGaussRys`) kernel** for off-center nuclei — the
  discrete curl identity that enforces antisymmetry is not satisfied.

### Recommended future work

Validate or replace the PSO field+derivative **recurrence** in
`comp_pso_int1_prim` against a reference integral formulation (or compute the PSO
integral on the atom-centered DFT grid). **Do not** spend further effort tuning
prefactors or the field factor — both are already verified correct.

## Reproducing the reference

The PySCF reference values were generated with `pyscf.prop.nmr.rhf` using
`gauge_orig` set to the same CGO origin, `_solve_mo1_uncoupled` for the
uncoupled response, and the `dia()` / `para()` helpers, with the unit conversion
`unit = pyscf.data.nist.ALPHA**2 * 1e6`.
