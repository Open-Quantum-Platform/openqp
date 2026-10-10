# `[input]`

The `[input]` section defines the molecule, electronic-structure method, run
type, basis convention, and process-level threading. Most workflows require this
section.

## Minimal Example

```ini
[input]
system=
   O   0.000000000   0.000000000  -0.041061554
   H  -0.533194329   0.533194329  -0.614469223
   H   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
method=hf
basis=6-31g*
```

## Keywords

### `system`

| Field | Value |
| --- | --- |
| Type | string or multiline coordinate block |
| Default | empty |
| Values | XYZ file path, or inline atom coordinates |
| Used by | all molecular workflows |

Defines the molecular geometry. When `system` is a single non-empty line, OpenQP
interprets it as a file path. When `system=` is followed by indented atom lines,
OpenQP reads those lines as inline coordinates.

Inline coordinates use one atom per line:

```ini
[input]
system=
   O   0.000000000   0.000000000  -0.041061554
   H  -0.533194329   0.533194329  -0.614469223
   H   0.533194329  -0.533194329  -0.614469223
```

Each atom line must contain at least `symbol x y z`. Numeric atomic labels are
also accepted by existing examples. When `basis=library`, each atom line must
also include a tag column.

### `system2`

| Field | Value |
| --- | --- |
| Type | string or multiline coordinate block |
| Default | empty |
| Values | second XYZ file path, or second inline geometry |
| Used by | NACME and previous-geometry workflows |

Provides a second geometry. NACME uses it as the displaced or previous-step
geometry when `[guess] file2` is not supplied.

### `charge`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` |
| Used by | electron count for all workflows |

Sets the total molecular charge. Keep `charge`, `[scf] multiplicity`, and the
chosen reference type physically consistent. When DFT-D4 is enabled, this same
charge is passed to DFT-D4's charge model; charged systems are therefore no
longer evaluated as neutral dispersion models.

### `method`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `hf` |
| Values | `hf`, `tdhf`, `mp2`, `ccsd`, `ccsd(t)`, `fci`, `casci`, `casscf`, `sa-casscf`, `caspt2`, `ms-caspt2`, `xms-caspt2`, `mrmp2`, `mcqdpt2`, `xmcqdpt2` |
| Used by | workflow dispatch |

Selects the electronic-structure driver. Use `method=hf` for HF and DFT
reference calculations. Use `method=mp2` for standalone ground-state MP2
correlation. Use [`method=ccsd` or `method=ccsd(t)`](cc.md) for energy-only
coupled cluster on an HF reference, controlled by [`[cc]`](cc.md). Use
`method=tdhf` for TDHF, TDDFT, SF-TDDFT, MRSF-TDDFT, SOC,
NACME, and MRSF-EKT workflows. The native wavefunction methods and their
required sections are summarized under [Wavefunction methods](wavefunction.md).
XTB, DFTB, and AFQMC are not distributed or supported as OpenQP 1.3.0 methods.

DFT calculations still use `method=hf`; the functional is selected separately
with `functional`. MP2, coupled-cluster, and native wavefunction calculations
require `functional` to be empty.

### `functional`

| Field | Value |
| --- | --- |
| Type | string |
| Default | empty |
| Used by | DFT and TDDFT-style calculations |

Selects a density functional. An empty value means Hartree-Fock. Examples use
values such as `bhhlyp`, `pbe0`, and DTCAM-series functionals where supported.

Some property implementations have functional restrictions. For example, NMR
shielding does not support range-separated or meta-GGA functionals in the
current checker.

Standalone MP2 uses an HF reference. With `method=mp2`, leave `functional`
empty; non-empty values are rejected before the calculation starts.

### `basis`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `6-31g*` |
| Values | basis name, semicolon-separated per-atom names, or `library` |
| Used by | AO basis setup |

Sets the orbital basis. A single basis name applies to every atom:

```ini
basis=6-31g*
```

Per-atom basis names can be supplied in atom order with semicolons:

```ini
basis=aug-cc-pVDZ-PP;aug-cc-pVDZ
```

For tagged basis assignment, set `basis=library`, add a tag after each atom, and
define the tag mapping with `library`.

### `library`

| Field | Value |
| --- | --- |
| Type | multiline string |
| Default | empty |
| Used by | tagged basis assignment when `basis=library` |

Maps atom tags to basis names:

```ini
[input]
system=
 C    1.6062782722   1.5141391221  -1.8538091464  c1
 H    0.7846511041   1.8564598303  -1.2260835006  h1
basis=library
library=
 c1 6-31g
 h1 6-31g*
```

When `basis=library`, every atom line must include a tag and every tag must have
a mapping.

### `runtype`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `energy` |
| Values | `energy`, `grad`, `hess`, `nac`, `nacme`, `bp`, `optimize`, `meci`, `mecp`, `tci`, `mep`, `ts`, `irc`, `neb`, `prop`, `data`, `ekt`, `soc` |
| Used by | top-level workflow dispatch |

Selects the calculation workflow.

Common values:

| Value | Meaning |
| --- | --- |
| `energy` | Single-point energy and requested properties. |
| `grad` | Energy plus gradient. |
| `hess` | Hessian and frequency workflow. |
| `nacme` | Nonadiabatic coupling matrix element workflow. |
| `soc` | Spin-orbit coupling workflow. |
| `ekt` | MRSF-EKT ionization/electron-affinity workflow. |
| `optimize` | Geometry optimization. |
| `meci`, `mecp` | Crossing-point searches. `[optimize] meci_search=baeka` makes MECI a two-or-more-state adaptive search. |
| `tci` | Legacy three-state adaptive-penalty runtype, retained unchanged for existing inputs. New BaekA calculations use `runtype=meci` with `meci_search=baeka`. |
| `ts`, `irc`, `neb`, `mep` | Reaction-path workflows. |
| `prop`, `data` | Multi-state property/data workflows for downstream drivers. |

`runtype=namd` is available in OpenQP 1.3.0. It selects surface-hopping
molecular dynamics configured by the [`[md]`](md.md) section and requires an
MRSF-TDDFT setup.

With an RHF reference, `method=mp2` accepts `energy`, `grad`, `optimize`, `ts`,
`mep`, and `irc`; the latter five use the analytic RHF-MP2 gradient. UHF and
ROHF MP2 remain energy-only. MP2 Hessians and other derivative workflows are
not implemented.

### `ispher`

| Field | Value |
| --- | --- |
| Type | string mode |
| Default | `auto` |
| Values | `auto`, `true`, `false` |
| Used by | AO basis shell convention |

Controls pure spherical harmonic versus Cartesian AO shells.

| Value | Meaning |
| --- | --- |
| `auto` | Follow basis-set metadata where possible. |
| `true` | Force pure spherical shells, such as 5d and 7f. |
| `false` | Force Cartesian shells, such as 6d and 10f. |

Use this explicitly when reproducing a calculation from another program or when
a workflow is validated for one convention.

### `d4`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | DFT-D4 dispersion correction |

Enables the DFT-D4 dispersion correction where supported. The input checker
requires a DFT functional when `d4=true`. By default, DFT-D4 loads the published
rational-damping parameters for that functional.

To supply an explicit rational-damping parameter set, provide all six values in
a `[d4]` section:

```ini
[input]
functional=pbe
d4=true

[d4]
s6=1.0
s8=0.95948085
s9=1.0
a1=0.38574991
a2=4.80688534
alp=16.0
```

| `[d4]` keyword | Meaning |
| --- | --- |
| `s6`, `s8`, `s9` | two-body and three-body dispersion scale factors |
| `a1`, `a2` | rational-damping length parameters |
| `alp` | three-body damping exponent |

The six values form one parameterization: specifying only a subset is an input
error. All values must be finite. In concise `.oqp` input the equivalent form is
`d4(s6=1.0,s8=0.95948085,s9=1.0,a1=0.38574991,a2=4.80688534,alp=16.0)`.

### `soc_2e`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Values | `0`, `1` |
| Used by | `runtype=soc` |

Controls whether mean-field two-electron SOC terms are included.

| Value | Meaning |
| --- | --- |
| `0` | One-electron SOC terms only. |
| `1` | One-electron plus mean-field two-electron SOC terms. |

The option lives in `[input]` because it gates the whole SOC workflow rather
than a response-solver detail.

Spin-orbit coupling is a relativistic interaction that mixes spin-free states of
different spin character. OpenQP's documented SOC workflow is an MRSF-TDDFT
workflow selected with `runtype=soc`; `soc_2e=1` adds the mean-field
two-electron SOC contribution used in practical molecular SOC calculations.
Scalar relativistic DKH correction is a separate spin-free Hamiltonian option
controlled by `[scf] scal_rel`. See
[References](../references.md#spin-orbit-coupling) for OpenQP's relativistic
MRSF-TDDFT SOC method and mean-field SOC operator background.

Python style:

```python
from oqp.openqp import OpenQP

job = OpenQP("soc_keywords")
job.molecule(geometry="water", charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31G(2df,p)", nstate=12)
job.workflow.soc(soc_2e=1, scal_rel=2)
```

For Python SOC workflows, `job.workflow.soc(...)` sets `[scf] scal_rel=2` by
default. Override it with `scal_rel=0`, `1`, or `2` when needed.

### `qmmm_flag`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | QM/MM dispatch |

Enables hybrid QM/MM calculations. When `true`, the QM region is embedded in a
classical (OpenMM) MM environment through the ESPF operator, configured by the
[`[qmmm]`](qmmm.md) section. It applies to single-point QM/MM energies,
ground-state QM/MM molecular dynamics, and nonadiabatic
[SOC-NAMD-QMMM](../workflows/soc-namd-qmmm.md) dynamics (`runtype=namd` with
`[md] soc=true`). Without `qmmm_flag=true` the `[qmmm]` section is ignored.

`qmmm_flag` and the `[qmmm]` schema are available in OpenQP 1.3.0.

### `omp_threads`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` |
| Used by | OpenMP runtime setup |

Sets OpenMP threads per process or MPI rank. `0` means leave the environment or
compiled default unchanged.

Precedence is:

1. command-line `--omp`
2. `[input] omp_threads`
3. `OMP_NUM_THREADS`
4. OpenQP built-in default

Example:

```ini
[input]
omp_threads=16
```

### `perf`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` (recommended production) |
| Values | `0`, `1`, `2`, `3` |
| Used by | the performance preset (see [Performance](../performance.md)) |

Opt-in performance preset that bundles the performance input keys into one
accuracy↔speed dial: `0` strict reference, `1` recommended production (exact),
`2` faster (tiny degradation), `3` aggressive (small degradation allowed).
Explicit performance input keys override the preset. `perf=-1` disables the preset.
