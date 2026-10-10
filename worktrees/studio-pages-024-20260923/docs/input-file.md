# Legacy `.inp` Input

The sectioned `.inp` format remains supported for existing input decks and
controls that do not have a concise spelling. New calculations should begin
with the recommended [`.oqp` format](oqp-input.md) or the
[Python API](python-scripting.md); workflow pages show those forms before their
legacy `.inp` equivalents.

A compact `.oqp` file describes the same calculation in one readable line:

```text
mrsf/bhhlyp/6-31g* opt
geom="h2o.xyz"
```

This means an MRSF-TDDFT optimization of `S0`; OpenQP selects the required
working reference automatically. Start with the [`.oqp` Quick
Start](oqp-input.md#quick-start) for state labels, SOC counts, and more
examples.

Legacy `.inp` files use the sectioned format documented on this page and
remain supported unchanged. Keep sectioned syntax in `.inp` and compact syntax
in `.oqp`; changing formats is optional. The correction assistant is a
secondary aid that produces an inspectable `.resolved.oqp` file, while the
resolved canonical file remains the authoritative calculation record.
OpenQP renders the short positional geometry above as the explicit canonical
spelling `geom="h2o.xyz"`.

The OpenQP repository ships a same-stem `.oqp` companion for every legacy
example `.inp`. This makes the complete legacy example inventory available
in both formats without removing the established input system.

The two formats intentionally expose different optimization detail. Concise
`.oqp` geometry drivers always use the native OpenQP engine and have no `lib`
selector. Traditional `.inp` files retain `[optimize] lib=oqp`,
`lib=geometric`, or `lib=scipy` for compatibility; geomeTRIC is an optional
legacy dependency used chiefly for advanced constraints beyond native frozen
distances.

Legacy `.inp` files are INI-like text files. Options are grouped by section:

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
| [`[d4]`](keywords/input.md#d4) | Complete explicit DFT-D4 rational-damping parameter set. |
| `[guess]` | Initial orbitals and restart data. |
| `[scf]` | RHF/ROHF/UHF reference and SCF convergence controls. |
| `[mp2]` | Standalone MP2 spin-scaling controls. |
| `[cc]` | Coupled-cluster frozen core and solver controls. |
| `[dftgrid]` | DFT functional/grid controls. |
| `[tdhf]` | TDHF, TDDFT, SF-TDDFT, MRSF-TDDFT, and UMRSF settings. |
| `[md]` | Nonadiabatic-dynamics controls used by `runtype=namd`. |
| `[qmmm]` | OpenMM QM/MM system and molecular-dynamics controls. |
| `[properties]` | Gradients, NAC, NMR, export, and property requests. |
| `[optimize]` | Geometry target and convergence controls; backend selection is retained for traditional `.inp` and Python compatibility. |
| `[oqp]` | Native optimizer, TS-Hessian, IRC, MEP, and NEB controls. |
| `[geometric]` | Optional legacy geomeTRIC controls for traditional `.inp` workflows. |
| `[pcm]` | Reference-SCF PCM/ddX energy settings. |
| `[symmetry]` | Point-group metadata and optional symmetry reductions. |
| `[hess]` | Hessian and frequency controls. |
| `[nac]` | NAC/NACME controls. |
| `[ekt]` | MRSF-EKT IP/EA channel selection. |
| `[neb]` | NEB product/image controls plus optional legacy geomeTRIC compatibility keys. |
| `[json]` | Advanced JSON/restart metadata. |
| `[tests]` | Internal regression-test expectations. |

## Run Types

Common `[input] runtype` values:

| Run type | Meaning |
| --- | --- |
| `energy` | Single-point energy and requested properties. |
| `grad` | Energy plus gradient for the requested state. |
| `hess` | Hessian/frequency workflow. |
| `nac`, `bp` | Numerical nonadiabatic-coupling vector and branching-plane workflows. |
| `nacme` | Time/geometric derivative coupling between MRSF states. |
| `soc` | MRSF-TDDFT or MRSF-TDHF spin-orbit coupling workflow. |
| `ekt` | MRSF-EKT ionization-potential/electron-affinity workflow. |
| `md` | Ground-state QM/MM molecular dynamics. The command-line runner dispatches this OpenMM path specially. |
| `namd` | Nonadiabatic molecular dynamics using `[md]` controls. |
| `optimize` | Geometry optimization. |
| `meci`, `mecp` | Crossing-point searches. `meci_search=baeka` selects the two-or-more-state adaptive MECI algorithm. |
| `tci` | Existing three-state adaptive-penalty workflow, retained for backward compatibility. It is distinct from the new general `meci_search=baeka` option. |
| `ts`, `irc`, `neb`, `mep` | Reaction-path workflows. |
| `prop`, `data` | Multi-state property/gradient workflows for downstream drivers. |

For QM/MM MD, use `[qmmm] n_steps=N`. The older `[qmmm] nsteps=N` keyword
remains available for legacy bookkeeping. Canonical `qmmm(...)` may accompany
`energy`, `md`, or `namd`; `md` requires it, while `namd` may also run gas
phase. Canonical QM/MM gradients and optimizations are rejected until their
active backends provide the required assembled gradient.

Bare `runtype=md` without QM/MM remains invalid. With `qmmm_flag=true`, both
the command-line path and programmatic `Runner` dispatch ground-state MD to the
OpenMM `QMMM_MD` driver. This applies after a concise `.oqp` request has been
lowered as well as to a traditional sectioned `.inp`. Run it without MPI, for
example `openqp file.oqp --nompi` or `openqp file.inp --nompi`.

Standalone MP2 is selected with `[input] method=mp2` and requires an empty
`[input] functional`. RHF references support `runtype=energy`, `grad`,
`optimize`, `ts`, `mep`, and `irc`; the derivative runtypes use the analytic
RHF-MP2 gradient. UHF and ROHF references remain energy-only. Spin-scaled MP2
variants are controlled by the optional `[mp2]` section and use the same
analytic derivative when the reference is RHF.
