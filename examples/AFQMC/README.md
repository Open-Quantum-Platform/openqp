# AFQMC examples

Phaseless AFQMC runs inside liboqp: `[input] method = afqmc` converges the SCF
reference, builds the Cholesky-factorized molecular Hamiltonian in memory, and
propagates the walkers. There is no separate executable and no
HAMILTONIAN/FCIDUMP file hand-off.

The examples are deliberately tiny end-to-end smoke tests on triplet
CH2/STO-3G: six projection steps with 16 walkers.

```sh
openqp CH2_ROHF_STO3G_AFQMC.inp     # single-determinant (mean-field) trial
openqp CH2_MRSF_STO3G_AFQMC.inp     # MRSF-TDDFT root 1 as the trial
openqp CH2_MRSF_STO3G_AFQMC.oqp     # the same run in the concise syntax
```

## Input forms

Sectioned `.inp` is explicit: `[input] method = afqmc`, plus `[tdhf] type =
mrsf` when the trial is an MRSF root, plus the `[afqmc]` section.

In the concise `.oqp` syntax an `afqmc(...)` call is what makes the job an
AFQMC job; the route then only says what is solved first and supplies the
reference:

```text
rohf/sto-3g geom="ch2.xyz" mult=3 energy afqmc(trial=mean_field)
mrsf(nstate=3)/bhhlyp/sto-3g geom="ch2.xyz" energy afqmc(trial=mrsf_state,state=1)
afqmc/sto-3g geom="ch2.xyz" mult=3 energy afqmc(walkers=256)
```

The bare `afqmc/basis` route is the mean-field-trial shorthand; it takes no
functional (the reference is HF) and supports `energy` only. `walkers`,
`steps`, `dt` and `chol` are accepted as short spellings of `nwalkers`,
`nsteps`, `timestep` and `chol_tol`.

## Trial wavefunctions

`[afqmc] trial` selects the trial:

| value        | trial                                                            |
|--------------|------------------------------------------------------------------|
| `mean_field` | the single SCF determinant                                       |
| `sf`         | an SF-CIS determinant expansion read from `trial_file`           |
| `mrsf`       | an MRSF-CIS determinant expansion read from `trial_file`         |
| `mrsf_state` | an MRSF-TDDFT root that OpenQP computes in the same job          |

`mrsf_state` needs `[tdhf] type = mrsf` (or `umrsf`) and a high-spin `mult = 3`
reference; `[afqmc] state` picks the root. The response vector is expanded into
its spin-complete determinant basis internally, so the trial carries the MRSF
open-shell partner rather than a plain spin-flip vector. Determinants below
`[afqmc] trial_threshold` are dropped, which is the main cost knob for the
multideterminant overlap, force bias and local energy.

Watch the `<S**2>` that AFQMC prints for the imported trial: it is the quickest
check that the requested root is in the spin sector you wanted.

## Optional state-specific inputs

`oo_orbitals`/`oo_file` import an OO-DFT, ROKS, MOM or SGM orbital rotation;
`nlow`/`low_file` deflate lower roots and monitor excited-state collapse. Both
are documented in the AFQMC repository's `MRSFCIS_FORMAT.md`.
