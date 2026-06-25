# Input Files

OpenQP inputs are INI-like text files. Options are grouped by section:

```ini
[input]
runtype=energy
method=hf
basis=6-31g*

[scf]
type=rhf
```

Lines beginning with `#` are comments. Keyword names are case-insensitive in
normal use, but this manual uses lower-case names to match the Python schema.

## Geometry

Inline coordinates are written under `[input] system` with indented atom lines:

```ini
[input]
system=
   O   0.000000000   0.000000000  -0.041061554
   H  -0.533194329   0.533194329  -0.614469223
   H   0.533194329  -0.533194329  -0.614469223
```

An external XYZ file can be used instead:

```ini
[input]
system=h2o.xyz
```

Some workflows, such as NACME, also use `[input] system2` for the displaced or
previous geometry.

## Core Sections

| Section | Purpose |
| --- | --- |
| `[input]` | Charge, basis, method, run type, geometry, AO convention, threading. |
| `[guess]` | Initial orbitals and restart data. |
| `[scf]` | RHF/ROHF/UHF reference and SCF convergence controls. |
| `[dftgrid]` | DFT functional/grid controls. |
| `[tdhf]` | TDHF, TDDFT, SF-TDDFT, MRSF-TDDFT, and UMRSF settings. |
| `[properties]` | Gradients, NAC, NMR, export, and property requests. |
| `[optimize]` | Geometry-optimization target and convergence controls. |
| `[oqp]` | Native optimizer controls. |
| `[geometric]` | geomeTRIC backend controls. |
| `[pcm]` | Reference-SCF PCM/ddX energy settings. |
| `[symmetry]` | Point-group metadata and optional symmetry reductions. |
| `[hess]` | Hessian and frequency controls. |
| `[nac]` | NAC/NACME controls. |
| `[ekt]` | MRSF-EKT IP/EA channel selection. |
| `[neb]` | NEB product and image controls. |

## Run Types

Common `[input] runtype` values:

| Run type | Meaning |
| --- | --- |
| `energy` | Single-point energy and requested properties. |
| `grad` | Energy plus gradient for the requested state. |
| `hess` | Hessian/frequency workflow. |
| `nacme` | Time/geometric derivative coupling between MRSF states. |
| `soc` | MRSF-TDDFT spin-orbit coupling workflow. |
| `ekt` | MRSF-EKT ionization-potential/electron-affinity workflow. |
| `optimize` | Geometry optimization. |
| `meci`, `mecp`, `tci` | Crossing-point searches. |
| `ts`, `irc`, `neb`, `mep` | Reaction-path workflows. |
| `prop`, `data` | Multi-state property/gradient workflows for downstream drivers. |

`md` is recognized by validation code but is not the current production
workflow in this repository.
