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
| PSO integrals `(r-R_N)x∇/r_N³` | `mod_1e_primitives.F90::comp_pso_int1_prim`, `int1.F90::pso_integrals` | Validated |
| Shielding assembly + output | `modules/nmr_shielding.F90` | Validated |
| Registration | `include/oqp.h`, `runfunc.py`, `input_checker.py` | Done |

## Validation

All numbers below are **H2O / STO-3G, RHF**, gauge origin = center of mass
`(0, 0, -0.19886422)` Bohr, compared against PySCF (`pyscf-properties`,
`pyscf.prop.nmr.rhf`, common-gauge, unit factor `alpha²·1e6`).

1. **Angular-momentum integrals** reproduce PySCF `int1e_cg_irxp` to ~8
   significant figures.
2. **Diamagnetic shielding** reproduces PySCF `dia()` to ~6 significant figures.
3. **Paramagnetic shielding** (uncoupled) reproduces the PySCF uncoupled
   reference for **all** atoms.

| Atom | σ_dia (OQP / PySCF) | σ_para (OQP / PySCF) | σ_total (OQP / ref) |
|------|---------------------|----------------------|---------------------|
| O    | 411.4176 / 411.4176 | -113.6305 / -113.6305 | 297.787 / 297.787   |
| H    | 28.0618 / 28.0618   | 1.78494 / 1.78494     | 29.847 / 29.847     |

## PSO antisymmetry (resolved)

The PSO operator `L_N/r_N³` is anti-Hermitian, so its real matrix representation
must be antisymmetric with a zero diagonal. The raw field+ket-derivative Rys
product acquired a small spurious *symmetric* component for nuclei not centered
on a basis function (a nonzero diagonal of order 0.04 a.u.), which previously
gave an incorrect hydrogen paramagnetic shielding (1.24 vs 1.78 ppm).

**Fix:** `pso_integrals` now assembles the *full* (both-triangle) matrix `M` and
returns the exact antisymmetric part `A = (M − Mᵀ)/2`. Because the exact
operator is antisymmetric, any symmetric component of `M` is pure error and is
removed exactly by this projection (`Mᵀ_symmetric = M_symmetric` cancels, while
`Mᵀ_antisymmetric = −M_antisymmetric` is preserved). After the fix the reported
diagnostics are `max|diag| = 0` and `max|A + Aᵀ| = 0`, and the hydrogen
paramagnetic shielding matches the reference (1.785 ppm).

The `(r-c)/r³` field factor and the `alpha²` prefactors were independently
confirmed correct and were not changed.

## Diagnostic checks

`nmr_shielding` reports a `PSO diagnostics  max|diag|, max|A+A^T|` line in the
log (both ~0). The pytest test `tests/test_nmr_shielding.py` asserts:

1. PSO matrix diagonal is ~0 (`max|diag|` from the log).
2. PSO + transpose(PSO) is ~0 (`max|A+A^T|` from the log).
3. H2O/STO-3G CGO@COM oxygen and hydrogen paramagnetic (and total) shieldings
   match the PySCF uncoupled reference within tolerance.

## Reproducing the reference

The PySCF reference values were generated with `pyscf.prop.nmr.rhf` using
`gauge_orig` set to the same CGO origin, `_solve_mo1_uncoupled` for the
uncoupled response, and the `dia()` / `para()` helpers, with the unit conversion
`unit = pyscf.data.nist.ALPHA**2 * 1e6`.

## Out of scope (future work)

GIAO (gauge-including atomic orbitals), hybrid functionals and the coupled HF
exchange response, UHF/ROHF open-shell references, and spin–spin coupling
constants.
