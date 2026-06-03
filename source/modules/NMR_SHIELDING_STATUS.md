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

## Phase 0 — coupled HF/hybrid magnetic response (ground state)

Extends the paramagnetic term from uncoupled to the coupled CPHF/CPKS response.
The first-order magnetic density `P^B` is imaginary/antisymmetric, so its Coulomb
response and the semi-local XC-kernel response vanish; only the exact-exchange
response survives, scaled by the functional's exact-exchange fraction `c_x`. The
coupled response is solved by fixed-point iteration of
`(eps_a-eps_i)R + c_x*K[P^B(R)] = b` (b = orbital-Zeeman occ-vir block; `K` the
antisymmetric exchange image built via the `int2` A-B path). For `c_x = 0` the
loop is skipped and the result equals the uncoupled response.

Validation vs the PySCF common-gauge oracle (H2O/STO-3G, CGO at COM):

| Functional | c_x | O para uncoupled | O para coupled (OQP / PySCF) |
|------------|-----|------------------|------------------------------|
| HF         | 1.00 | -113.63 | **-230.63 / -230.63** (exact) |
| BHHLYP     | 0.50 | -156.76 | -242.11 / -242.06 |
| PBE0       | 0.25 | -194.30 | -248.97 / -248.88 |
| PBE        | 0.00 | -258.44 | **-258.44 / -258.32** (coupled == uncoupled) |

HF (grid-free) matches the oracle exactly. For the DFT functionals the absolute
numbers differ from PySCF by ~0.1 ppm due to cross-code DFT-SCF/grid differences
(the diamagnetic part differs only ~0.003 ppm); the **coupling contribution**
`Delta = coupled - uncoupled` matches the oracle to ~0.03 ppm and scales
monotonically with `c_x` — the SCF-robust validation of the coupling itself.

**Gates (automated in `tests/test_nmr_coupled.py`):** gate 0 `max|P^B+P^B^T|~1e-16`;
gate 1 `||J(P^B)||~1e-16` (Coulomb vanishes); gate 2 `||K(P^B)||~0.54` (exchange
nonzero); gate 3 PBE coupled == uncoupled (exact); gate 4 HF matches oracle and the
coupling `Delta` matches for all functionals; gate 6 `|Delta|` scales with `c_x`.
Oracle: `tests/fixtures/nmr/generate_pyscf_cgo_reference.py` -> `pyscf_cgo_reference.json`.

## GIAO development checkpoint (branch `feat/giao-nmr`)

CGO remains the only validated production NMR path.  The Python interface now
recognizes `properties.nmr_gauge=cgo|giao`; `cgo` dispatches to this validated
CGO implementation, while `giao` is deliberately gated with a `NotImplementedError`
until the complete native London-orbital integral and response path passes the
benchmark matrix.

Current native GIAO implementation status:

| Piece | File / symbol | Status |
|-------|---------------|--------|
| Gauge selector and runtime gate | `oqpdata.py`, `input_checker.py`, `runfunc.py` | Implemented; GIAO still gated |
| Benchmark/oracle scaffold | `scripts/nmr_giao_benchmark_matrix.py`, `tests/fixtures/nmr/benchmark_results/` | Implemented; PySCF CGO/GIAO oracle + OpenQP CGO rows populated |
| GIAO overlap magnetic derivative `S10` | `mod_1e_primitives.F90::comp_giao_overlap_deriv_prim`, `int1.F90::giao_overlap_derivative` | Native one-electron building block implemented |
| GIAO first-order core Hamiltonian `h10` | `mod_1e_primitives.F90::comp_giao_h10_core_prim`, `int1.F90::giao_h10_core`; PySCF oracle calls `make_h10(..., gauge_orig=None)` | Native one-electron building block implemented; validation in progress |
| GIAO two-electron magnetic derivative contractions | `nmr_giao_debug.F90::giao_h10_twoe_matrix` (+ `nmr_giao_h10_twoe_debug` emitter); PySCF oracle calls `pyscf.prop.nmr.rhf.get_jk` | Native RHF debug contraction implemented and PySCF-validated; refactored into the reusable `giao_h10_twoe_matrix` builder |
| GIAO CPHF/CPKS RHS/response assembly (paramagnetic) | `nmr_giao_shielding.F90::nmr_giao_shielding_debug` | **Implemented and PySCF-validated** (uncoupled + coupled). h1=h10(1e+2e), s1=S10, MO transform, GIAO CPHF/CPKS first-order solve (exchange-only coupled response scaled by `c_x`), PSO contraction. See note below. |
| GIAO diamagnetic term (second-order GIAO integrals) | pending (`int1e_giao_a11part`, `int1e_a01gp` analogues) | Not implemented |
| OpenQP native GIAO shielding output (total) | pending | Not implemented; no CGO fallback allowed (needs the diamagnetic term) |

### GIAO paramagnetic checkpoint — VALIDATED

`nmr_giao_shielding.F90::nmr_giao_shielding_debug` assembles the native GIAO
paramagnetic nuclear shielding and emits machine-parseable `GIAO_SHIELDING_DEBUG_*`
records (still a debug emitter; `nmr_gauge=giao` stays gated until the diamagnetic
term completes the total).  It mirrors the PySCF GIAO algorithm exactly: the
first-order magnetic Hamiltonian `h1 = h10(one-electron) + h10(two-electron)`, the
first-order overlap `s1 = S10`, transformed to the MO basis; the GIAO CPHF/CPKS
first-order equation solved both uncoupled (`mo1[vir]=-hs/(e_a-e_i)`,
`mo1[occ]=-0.5 s1`) and coupled (fixed point with the exact-exchange image of the
imaginary/antisymmetric first-order density, scaled by `c_x`, via the validated
`int2` A-B path and `mntoia`); and the paramagnetic tensor formed by contracting
the resulting first-order density with the PSO operator (note OpenQP
`pso_integrals` is the negative of PySCF `int1e_prinvxp`).

Validation (H2O/STO-3G, isotropic paramagnetic shielding, ppm) vs the committed
PySCF GIAO oracle (`tests/test_nmr_giao_para_live.py`,
`tests/fixtures/nmr/pyscf_giao_reference.json`):

| Functional | c_x | O para uncoupled (OQP/ref) | O para coupled (OQP/ref) |
|------------|-----|----------------------------|--------------------------|
| HF     | 1.00 | 27.906 / 27.906 | **-6.3026 / -6.3026** (exact) |
| BHHLYP | 0.50 | 37.953 / 37.950 | 25.490 / 25.487 |
| PBE0   | 0.25 | 51.534 / 51.514 | 48.993 / 48.967 |
| PBE    | 0.00 | 73.848 / 73.811 | 73.848 / 73.811 (coupled == uncoupled) |

HF (grid-free) matches the oracle to ~1e-4 ppm; the DFT absolute offsets (≤0.04 ppm)
are cross-code DFT-SCF/grid differences and the coupling `Delta = coupled −
uncoupled` matches to ≤0.03 ppm — the same SCF-robust pattern established for CGO.
Pure PBE coupled ≡ uncoupled exactly (`c_x = 0`).  The GIAO paramagnetic numbers are
small/positive (very different from CGO's), so the **diamagnetic GIAO term carries
the bulk of the total** and is the remaining piece before the total shielding and
`nmr_gauge=giao` ungating.

### Remaining diamagnetic GIAO term — precise spec for the next coding gate

The GIAO diamagnetic shielding (per nucleus N, contracted with the ground-state
density `dm0`, no response) follows `pyscf.prop.nmr.rhf.dia(gauge_orig=None)`:

```
e11_{ab} = sum_{mu,nu} dm0_{mu,nu} <mu| A11part_{ab} |nu>      (a,b = 1..3, 9 comp)
e11      = e11 - I*trace(e11)                                  (trace correction)
e11_{ab} += sum_{mu,nu} dm0_{mu,nu} <mu| A01gp_{ab} |nu>
sigma_dia(N) = e11 * alpha^2 * 1e6     [ppm]
```

Two new second-order GIAO one-electron integral blocks are required (libcint
operator strings, both 9-component, `rinv` origin at nucleus N):

- `int1e_giao_a11part` = `(-.5 | nabla-rinv | r)` — GIAO diamagnetic.
  **Verified NOT equal to** the existing CGO `int1e_cg_a11part` at common origin 0
  (`||giao - cg(O=0)|| = 0.80` on H2O/STO-3G), i.e. the London-orbital prefix adds
  genuine structure; it cannot be obtained by calling `nmr_dia_shielding` with the
  gauge origin at 0.  A native kernel is needed.
- `int1e_a01gp` = `(g | nabla-rinv cross p |)` — GIAO gauge-correction: the
  London-derivative-weighted PSO operator `nabla(1/r_N) x nabla` (per-component
  mixed symmetry in mu,nu).

Magnitudes (H2O/STO-3G, isotropic, ppm) to target during validation, from the
PySCF oracle decomposition:

| atom | a11part iso | a01gp iso | GIAO dia total | (CGO dia for reference) |
|------|-------------|-----------|----------------|--------------------------|
| O | 386.94 | -8.95 | 377.99 | 411.42 |
| H | 26.76 |  4.22 |  30.98 |  28.06 |

Reusable building blocks: the GIAO position-raising trick already used in
`giao_h10_twoe_matrix` (raise the bra Cartesian power by one, combine with the
shell-center to form `r`-weighted integrals), the `rinv`/PSO primitive path
(`comp_pso_int1_prim`), and the CGO diamagnetic primitive
(`comp_nmr_dia_int1_prim`).  Validation oracle:
`tests/fixtures/nmr/pyscf_giao_reference.json` (`sigma_dia*`,
`sigma_total_*` tensors) and `generate_pyscf_giao_reference.py`.

Once both integrals are implemented and `sigma_dia + sigma_para` matches the GIAO
oracle total to the CGO-grade tolerance (HF ~1e-4 ppm; DFT Delta ≤ ~0.1 ppm),
wire the total into `nmr_giao_shielding.F90`, route `nmr_gauge=giao` through it in
`runfunc.py`, and remove the `NotImplementedError` gate.  Until then, the gate
**stays closed**: no CGO fallback, no partial-dia totals.

The implemented `S10` block returns the real coefficient of the imaginary first
magnetic-field derivative of the AO overlap matrix,
`S10_a(mu,nu) = 0.5 * [(R_mu - R_nu) x <mu|r|nu>]_a`, and is intentionally not
called by `runfunc.py` or the production CGO shielding routine.  The native
`h10` one-electron block is now implemented alongside it as real-valued storage
for the imaginary first-order GIAO one-electron operator.  Its current scope is
the PySCF/libcint one-electron convention
`h10_onee = -0.5*int1e_giao_irjxp - int1e_ignuc(asym) - int1e_igkin`.
These blocks become useful for production only after the validated two-electron
derivative image is connected to the GIAO response terms and the resulting
shieldings are benchmarked against the PySCF GIAO oracle.

`tests/test_nmr_giao_h10_twoe_live.py` is the live RHF two-electron checkpoint:
it compares the native `vj`, `vk`, and `twoe_h10 = vj - 0.5*vk` debug matrices
for H2 and H2O/STO-3G against PySCF `pyscf.prop.nmr.rhf.get_jk(mol, dm0)` and
keeps the production `nmr_gauge=giao` route gated.

Conservative `h10` implementation note: in this codebase `h10` means the native
first-order core-Hamiltonian magnetic derivative for the GIAO/London-orbital
response path.  It is a one-electron AO magnetic-perturbation block and must not
be confused with production CGO angular-momentum, PSO, or diamagnetic shielding
integrals.  The existing AO core Hamiltonian is formed in `int1.F90::omp_hst` as
nuclear attraction plus kinetic energy (`h = h + t` after separate `nuc_ints` and
`kin_ovl_ints` assembly), so the native `h10` block includes the GIAO irjxp
one-electron term plus both the kinetic magnetic-derivative contribution and the
nuclear-attraction magnetic-derivative contribution.  Any overlap/metric term needed by the response equation must be
handled consistently with the existing `S10` overlap derivative and the eventual
GIAO response convention.

The `h10` source locations should stay beside the existing one-electron integral
code, not in the production CGO shielding path.  Primitive kernels naturally
belong in `mod_1e_primitives.F90`, near `comp_kin_ovl_int1_prim`,
`comp_coulomb_int1_prim`, and `comp_giao_overlap_deriv_prim`.  Assembled AO
matrices naturally belong in `int1.F90`, parallel to `giao_overlap_derivative`
and the current one-electron assembly routines.  The current magnetic-derivative
convention is real-valued storage for the real coefficient of an imaginary
operator/derivative: `S10` omits the common factor `i`, angular momentum stores
real antisymmetric `A` with physical `L = -i A`, and PSO stores a real
antisymmetric block with physical operator `-i A`.  Future `h10` work must match
that phase/sign convention or explicitly document a tested conversion.

Component ordering is Cartesian `(x, y, z)` in the second dimension or third
array dimension for magnetic perturbation blocks: examples include packed
`ints(nbf*(nbf+1)/2,3)` for `S10`/angular momentum and full
`pso_full(nbf,nbf,3)` for PSO.  The `h10` matrix tests compare AO and, if needed,
MO-transformed component blocks against the PySCF oracle; they must not generate
native shielding output.  `h10` alone still does not enable native
GIAO shielding: `nmr_gauge=giao` remains gated until `h10`, two-electron magnetic
derivatives, the GIAO CPHF/CPKS RHS/response assembly, and PySCF-oracle benchmark
validation are all complete.  Do not use placeholder zero `h10` matrices, CGO
values as GIAO `h10`, fake shielding rows, benchmark rows, or a CGO fallback to
make native GIAO appear available.

## ROHF implementation / audit checkpoint

ROHF ground-state SCF support and ROHF NMR response support are separate tasks.
The current code recognizes `scf.type=rohf` as SCF type 3 and contains a native
ROHF SCF path in `source/scf.F90`, including alpha/beta electron counts,
ROHF occupations, shared spatial orbitals, UROHF two-electron Fock construction,
and a Guest-Saunders effective ROHF Fock transformation.  Existing SF/MRSF input
paths also require an ROHF reference.  This is ground-state/reference support; it
is not a claim of validated ROHF NMR shielding.

Current conservative ROHF/NMR status:

| Piece | Status |
|-------|--------|
| ROHF input/interface recognition | Present (`OQPData._scftypes`, input checker, examples) |
| ROHF ground-state SCF driver | Present but should remain under audit/regression coverage before being used as an NMR benchmark foundation |
| ROHF-CGO NMR response | Not implemented / not validated |
| ROHF-GIAO NMR response | Not implemented / not validated |
| ROHF-GIAO shielding | Not validated; `nmr_gauge=giao` remains gated with `NotImplementedError` |
| RHF/UHF fallback for requested ROHF | Not allowed |

ROHF is required before open-shell NMR benchmarks that need a restricted
open-shell reference, but ROHF SCF alone is insufficient for open-shell NMR.  A
future ROHF NMR implementation must explicitly add the response equations for
closed-shell/doubly occupied and singly occupied open-shell subspaces, alpha/beta
occupation handling, total and spin density requirements, effective ROHF Fock
contributions in the response, DIIS/convergence integration where applicable, and
a benchmark against an external oracle.  Until then, requesting NMR shielding with
`scf.type=rohf` raises a clear `NotImplementedError`; no requested ROHF case may
be silently routed through RHF, UHF, or the CGO fallback for GIAO.

ABI/build caution: ROHF and NMR fixes must preserve the integer-sensitive
BLAS/LAPACK and compiler-default integer behavior in the current CMake context;
do not hide type mismatches behind global compiler flags.

## Out of scope (future work)

Validated OpenQP GIAO shieldings, UHF/ROHF open-shell references, MRSF-TDDFT NMR
(see `MRSF_NMR_DESIGN.md`), and spin-spin coupling constants.
