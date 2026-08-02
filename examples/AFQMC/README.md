# AFQMC examples

Phaseless AFQMC runs inside liboqp: `[input] method = afqmc` converges the SCF
reference, builds the Cholesky-factorized molecular Hamiltonian in memory, and
propagates the walkers. There is no separate executable and no
HAMILTONIAN/FCIDUMP file hand-off.

Both examples are deliberately tiny end-to-end smoke tests on triplet
CH2/STO-3G: six projection steps with 16 walkers.

```sh
openqp CH2_ROHF_STO3G_AFQMC.inp     # single-determinant (mean-field) trial
openqp CH2_MRSF_STO3G_AFQMC.inp     # MRSF-TDDFT root 1 as the trial
```

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
