# `[cc]`

!!! note "Available in OpenQP 1.3.0"
    Energy-only CCSD and CCSD(T) are included in OpenQP 1.3.0.

The `[cc]` section controls the frozen core and the coupled-cluster solver for
`[input] method=ccsd` and `[input] method=ccsd(t)`. The defaults are usable as
they stand; the section is only needed to freeze core orbitals, to change the
convergence behaviour, or to override how the ladder integrals are stored.

Coupled cluster is a post-SCF energy workflow. It requires an HF reference, so
`[input] functional` must be empty and `[input] runtype` must be `energy`. In
concise input, select the reference with
`ccsd_t(reference=rhf|rohf|uhf)/BASIS`; the value lowers to `[scf] type` and is
not a keyword in this section.

## Minimal Examples

CCSD(T) in `.oqp`:

```text
ccsd_t/6-31g
geom="h2o.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("h2o_ccsd_t")
job.molecule(geometry="water")
job.theory.ccsd_t(basis="6-31g")
mol = job.run()
```

Legacy `.inp`:

```ini
[input]
method=ccsd(t)
runtype=energy
functional=
```

To freeze core orbitals, add an exact `.oqp` section call:

```text
ccsd_t/6-31g cc(nfzc=1)
geom="h2o.xyz"
```

The legacy section is:

```ini
[cc]
nfzc=1
```

## Keywords

### `nfzc`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` |
| Used by | CCSD and CCSD(T) |

Number of lowest orbitals excluded from the correlation treatment. Must be
non-negative and smaller than the number of occupied orbitals; the module
aborts otherwise.

For an ROHF reference the core is removed before semicanonicalisation, so the
correlated space is the span of the reference orbitals that were kept. See the
[coupled-cluster workflow](../workflows/coupled-cluster.md#frozen-core).

### `conv`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1e-7` |
| Used by | CCSD and CCSD(T) |

Convergence threshold. The iteration stops when both the amplitude RMS change
and the change in the correlation energy fall below this value.

### `maxit`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `50` |
| Used by | CCSD and CCSD(T) |

Maximum number of CCSD iterations. A run that does not converge within this
many iterations aborts rather than reporting an unconverged energy.

### `ndiis`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `8` |
| Used by | CCSD and CCSD(T) |

DIIS subspace size for the amplitude iteration. Set `0` to disable DIIS; the
same tolerance then typically needs two to three times as many iterations.

### `cholesky`

| Field | Value |
| --- | --- |
| Type | `auto`, `true`, or `false` |
| Default | `auto` |
| Used by | CCSD and CCSD(T), closed-shell reference only |

Whether to hold the ladder integrals `(ab|cd)` as an explicit `v^4` array or to
rebuild them from Cholesky vectors as the iteration needs them.

Factorising trades arithmetic for memory. The vectors are far smaller than
`v^4`, but every block of ladder integrals has to be reassembled from them, and
that assembly costs `nchol / no^2` times the ladder contraction it feeds. The
number of vectors tracks the basis set while `no^2` tracks the correlated
electrons, so a small molecule in a large basis is the worst case: H<sub>2</sub>O
in cc-pVQZ has 1680 vectors and four correlated occupied orbitals, where
factorising made the CCSD iterations about five times slower for the same
energy.

`auto` therefore takes the vectors only when the explicit route will not fit,
sizing it against the memory actually available (see
[`cholesky_direct`](#cholesky_direct) for how that is probed). On a 500 GB node
the crossover is near 460 basis functions; on a 16 GB laptop it is near 200.
Set `true` to factorise regardless, or `false` to force the explicit array and
let the memory guard refuse the job if it does not fit.

### `cholesky_tol`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1e-10` |
| Used by | CCSD and CCSD(T), when Cholesky vectors are in use |

Truncation threshold for the decomposition. The pivoted Cholesky stops once the
largest remaining diagonal falls below this value. Must be positive and finite;
the module aborts otherwise, because a non-positive or infinite threshold either
never terminates or silently produces zero vectors.

Looser values trade accuracy for both memory and time. The default is tight
enough that the correlation energy matches the explicit ladder to better than
`1e-9` Hartree. If the vector cap is reached before the threshold is met, the
log warns that the result is less accurate than requested.

### `cholesky_direct`

| Field | Value |
| --- | --- |
| Type | `auto`, `true`, or `false` |
| Default | `auto` |
| Used by | CCSD and CCSD(T), when Cholesky vectors are in use |

Whether to build the vectors straight from recomputed AO integrals instead of
from the packed AO integral store.

The direct route never allocates the packed store, which is the larger of the
two for a big basis, but it sweeps the shell-pair list once per pivot block and
is measurably slower wherever both fit. It is therefore chosen on memory, never
on speed: `auto` takes it only when the packed store would not fit.

Available memory is the tightest of physical RAM, the kernel's `MemAvailable`,
and the cgroup limit — the last being what actually binds under SLURM or in a
container. `OQP_MEMORY_LIMIT_GB` overrides the probe when the automatic answer
is wrong.

## Notes

- `method=ccsd` and `method=ccsd(t)` accept `runtype=energy` only.
- Both reject non-empty `[input] functional` values.
- The reference is selected through `[scf] type`; RHF, UHF, and ROHF are
  supported. Open-shell references run through a spin-orbital solver whose
  storage grows as `(2*nmo)^4`, so keep those systems small.
- The integrals are held in memory on both paths. The module prints the storage
  it needs and refuses when that exceeds what the machine can give, measured at
  run time rather than against a fixed ceiling — see
  [`cholesky_direct`](#cholesky_direct) for what is probed and how to override
  it.

See the [coupled-cluster workflow](../workflows/coupled-cluster.md) for
complete input and Python examples.
