# `[input]`

| Keyword | Default | Notes |
| --- | --- | --- |
| `charge` | `0` | Total molecular charge. |
| `basis` | `6-31g*` | Basis-set name. |
| `library` | empty | Optional custom basis mapping. |
| `functional` | empty | Empty means Hartree-Fock. |
| `method` | `hf` | Use `hf` for HF/DFT, `tdhf` for response workflows. |
| `runtype` | `energy` | Selects the calculation workflow. |
| `system` | empty | Inline geometry or XYZ path. |
| `system2` | empty | Second geometry for workflows such as NACME. |
| `ispher` | `auto` | AO convention: `auto`, `true`, or `false`. |
| `d4` | `False` | DFT-D4 correction where supported. |
| `soc_2e` | `1` | SOC terms for `runtype=soc`. |
| `omp_threads` | `0` | Threads per MPI rank; `0` leaves environment/defaults alone. |

## Run Types

Common values are `energy`, `grad`, `hess`, `nac`, `nacme`, `soc`, `ekt`,
`optimize`, `meci`, `mecp`, `tci`, `mep`, `ts`, `irc`, `neb`, `prop`, and
`data`.

`md` may be recognized by validation code but should not be documented as the
main production workflow in this repository.

## AO Convention

`ispher=auto` follows basis-set metadata where possible. Use `ispher=true` to
force pure spherical harmonic shells and `ispher=false` to force Cartesian
shells.

## OpenMP Threads

`omp_threads` is lower priority than the command-line `--omp` flag and higher
priority than the existing `OMP_NUM_THREADS` environment variable.
