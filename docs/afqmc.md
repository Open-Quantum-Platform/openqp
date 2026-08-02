# AFQMC

Phaseless auxiliary-field quantum Monte Carlo (ph-AFQMC) for molecular
Hamiltonians, built into `liboqp`.

The whole method runs inside the library. `[input] method = afqmc` converges the
SCF reference, builds the Cholesky-factorized molecular Hamiltonian in memory
from the tagarray, and propagates the walkers — there is no separate executable,
no `HAMILTONIAN`/`FCIDUMP` hand-off, and no return to Python inside the loop.
Nothing extra needs to be installed: `pip install .` ships it.

The kernels are vendored from
[Open-Quantum-Platform/AFQMC](https://github.com/Open-Quantum-Platform/AFQMC)
under `source/afqmc/`, with `afqmc_host.F90` supplying the handful of host
symbols they expect. `afqmc_oqp.F90` is the OpenQP entry point.

## Quick start

```sh
openqp examples/AFQMC/CH2_ROHF_STO3G_AFQMC.inp     # single-determinant trial
openqp examples/AFQMC/CH2_MRSF_STO3G_AFQMC.inp     # MRSF-TDDFT root as trial
```

## Input

### Sectioned (`.inp`)

```text
[input]
method=afqmc
runtype=energy
basis=sto-3g
functional=

[scf]
type=rohf
multiplicity=3

[tdhf]              # only when the trial is an MRSF root
type=mrsf
nstate=3

[afqmc]
trial=mrsf_state
state=1
nwalkers=640
nsteps=1000
```

### Concise (`.oqp`)

An `afqmc(...)` call is what makes the job an AFQMC job. The route then only says
what is solved first and supplies the reference:

```text
rohf/sto-3g geom="ch2.xyz" mult=3 energy afqmc(trial=mean_field)
mrsf(nstate=3)/bhhlyp/sto-3g geom="ch2.xyz" energy afqmc(trial=mrsf_state,state=1)
afqmc/sto-3g geom="ch2.xyz" mult=3 energy afqmc(walkers=256)
```

The bare `afqmc/basis` route is the mean-field-trial shorthand: it takes a basis
but no functional (the reference is HF), defaults to `S0`, and supports `energy`
only. Short spellings accepted inside `afqmc(...)`: `walkers`, `steps`, `dt`,
`chol` for `nwalkers`, `nsteps`, `timestep`, `chol_tol`.

## `[afqmc]` keywords

| keyword | default | meaning |
|---|---|---|
| `trial` | `mean_field` | trial wavefunction; see below |
| `trial_file` | — | determinant expansion for `trial = sf` / `mrsf` |
| `state` | `1` | 1-based MRSF root used by `trial = mrsf_state` |
| `trial_threshold` | `1.0e-8` | drop trial determinants with \|c\| below this |
| `chol_tol` | `1.0e-10` | Cholesky truncation; also bounds the reconstruction error |
| `nwalkers` | `640` | walker population |
| `nsteps` | `1000` | projection steps |
| `timestep` | `0.005` | imaginary-time step Δτ (Hartree⁻¹) |
| `seed` | `1` | RNG seed for the auxiliary fields |
| `stabilize_every` | `5` | steps between walker re-orthogonalizations (QR) |
| `population_control_every` | `5` | steps between branching / weight renormalization |
| `estimate_every` | `25` | steps between mixed-estimator evaluations |
| `accumulate_after` | `100` | start block-averaging after this many steps |
| `force_bias_bound` | `1.0` | clip on \|x̄\|, the force-bias magnitude |
| `oo_orbitals` | `False` | apply an imported state-specific orbital rotation |
| `oo_file` | — | that rotation (OO-DFT / ROKS / MOM / SGM), staged as `OOORBDAT` |
| `nlow` | `0` | number of lower roots to deflate out of the trial |
| `low_file` | — | their coefficient vectors, staged as `LOWSTATE` |
| `low_max` | `0.0` | stop if coherent lower-state leakage exceeds this; ≤ 0 reports only |
| `low_cap` | `10.0` | per-walker ratio cap in the coherent leakage mean; ≤ 0 disables |
| `low_start` | `0` | first step at which `low_max` is enforced |

Thread count comes from the usual `[input] omp_threads` / `OMP_NUM_THREADS`;
there is no AFQMC-specific thread keyword.

The file formats for `trial_file`, `oo_file` and `low_file` are documented in the
AFQMC repository's `MRSFCIS_FORMAT.md`. Files are copied into the working
directory under the fixed names the readers expect (`SFDAT`, `MRSFDAT`,
`OOORBDAT`, `LOWSTATE`).

## Trial wavefunctions

| `trial` | source |
|---|---|
| `mean_field` | the single converged SCF determinant |
| `sf` | an SF-CIS determinant expansion read from `trial_file` |
| `mrsf` | an MRSF-CIS determinant expansion read from `trial_file` |
| `mrsf_state` | an MRSF-TDDFT root that OpenQP computes in the same job |

### `mrsf_state`

Requires `[tdhf] type = mrsf` (or `umrsf`) and the high-spin `multiplicity = 3`
reference MRSF already uses. `[afqmc] state` picks the root.

MRSF works from an M_S=+1 reference with `nocca` alpha and `noccb` beta occupied
orbitals; the two open-shell orbitals are `noccb+1` and `noccb+2 = nocca`. A
response amplitude X(i,j) moves one alpha electron out of occupied orbital `i`
into beta orbital `j > noccb`, so it is exactly one M_S=0 determinant:

```text
alpha occupied : 1..nocca without i        (nocca-1 electrons)
beta  occupied : 1..noccb together with j  (noccb+1 electrons)
```

The stored eigenvector is expanded into those determinant amplitudes with
`tdhf_mrsf_lib`'s own `mrsfxvec`, which is what supplies the spin-complete
open-shell partner — the ±1/√2 pair that distinguishes MRSF from a plain
spin-flip vector. Reusing OpenQP's routine is what makes this a genuine MRSF
trial rather than a relabelled SF one.

Because a spin-flip trial lives in a different M_S sector than the reference, the
AFQMC electron counts come from the trial, not from the SCF.

**Check the `<S^2>` that AFQMC prints for the imported trial.** It is the
quickest confirmation that the requested root is in the spin sector you wanted;
it also validates the determinant phase convention, since a wrong relative sign
in the open-shell pair turns a singlet into a triplet.

## Hamiltonian construction

`h_pq` is the AO core Hamiltonian transformed with the converged MO
coefficients. The two-electron part is factorized as
`(pq|rs) = Σ_γ L^γ_pq L^γ_rs` by a **pivoted incomplete Cholesky in the AO
basis**, after which each vector is rotated to MOs. That is exact rather than an
approximation: the factorization is preserved term by term under
`L^γ → Cᵀ L^γ C`, and the AO route never forms the `nmo^4` MO tensor.

Only the residual diagonal and the accepted vectors are carried, with the inner
update done by DGEMV, so the cost is `O(nbf² · nchol²)`. Since the ERI
supermatrix is positive semidefinite, the largest residual diagonal at exit
bounds every residual element — the "residual diagonal bound" printed in the log
is therefore a rigorous bound on the reconstruction error, obtained for free.

The AO ERIs come from `int2e`, which materializes the full `nbf^4` tensor in
core. That is what bounds the system size.

## Output

```text
 Orbitals (nmo)          : 24
 Electrons (alpha/beta)  : 5/3
 Cholesky vectors        : 227 (threshold  1.00E-08)
 Residual diagonal bound :   4.9E-09

 AFQMC propagation
   step        E(mixed)      E(hybrid)   total weight|  E(mixed) avg  E(hybrid) avg
      1  -38.428996645781  ...
```

`E(mixed)` is the mixed estimator ⟨Ψ_T|H|Φ⟩/⟨Ψ_T|Φ⟩ and is the physical answer;
`E(hybrid)` is the hybrid estimator used to propagate the weights. The two
right-hand columns are running averages over the steps after
`accumulate_after`. The reported molecular energy is the block-averaged mixed
estimator.

A useful sanity check: with `trial = mean_field`, the step-1 mixed estimator must
equal the SCF total energy, because at step 0 every walker is the reference
determinant. On CH2/STO-3G ROHF that is `-38.428996645781` against a printed SCF
energy of `-38.4289966458` — agreement there means `h1`, the Cholesky vectors
and the core energy are all correct.

## Threading and performance

The walker loop, the re-orthogonalization and the local-energy estimator are all
OpenMP-parallel over walkers. The auxiliary fields are drawn in one block before
the threaded region, and the estimator accumulates in walker order rather than
through a reduction clause, so **results do not depend on `OMP_NUM_THREADS`**.

Set the BLAS to one thread (`MKL_NUM_THREADS=1` / `OPENBLAS_NUM_THREADS=1`) and
give the threads to AFQMC: the per-walker matrices are small, so nested BLAS
threading only oversubscribes.

Scaling on CH2/cc-pVDZ (nbf=24, nchol=227), 640 walkers × 100 steps,
propagation only:

| threads | wall |
|--------:|-----:|
| 1  | 42.0 s |
| 8  |  5.4 s |
| 32 |  3.3 s |

The dominant per-walker cost is streaming the Cholesky tensor
(`nbf² × nchol` doubles) in the force bias and in the two-body propagator; both
contract it in a single DGEMM pass. After that the largest remaining knobs are
`nwalkers`, `nsteps`, and — for multideterminant trials — the determinant count,
which `trial_threshold` controls.

## Limitations

- Energies only. There are no AFQMC gradients; the `.oqp` route rejects
  non-`energy` drivers.
- The AO ERI tensor is held in core (`O(nbf^4)`), which caps the basis size.
- `nlow` deflation is rigorous only within the supplied determinant trial space.
  It does not modify the Hamiltonian and is not a full multistate projector.
- The phaseless bias of a truncated MRSF trial is not measured by any of the
  shipped examples. Production excited-state work should converge time step,
  walker population, projection length and trial space.
