# `.oqp` Input

The `.oqp` format describes one calculation as a short readable sequence:
the method route, the task, and the geometry. Write everything except the
geometry on one line and put the geometry last; splitting the same items over
several lines is also accepted and means exactly the same thing. Write only
what differs from the defaults: a bare `opt` with no SCF, guess, or optimizer
keywords is a complete, production-ready request. It has no leading `#` and no
`[input]` or `[scf]` blocks.

!!! note "Available in OpenQP 1.3.0"

    `.oqp` is a supported input style in OpenQP 1.3.0. Sectioned `.inp` files
    remain supported for compatibility and advanced exact-section controls.

## Quick Start

Put `h2o.xyz` beside a file named, for example, `h2o-opt.oqp`, and write:

```text
mrsf/bhhlyp/6-31g* opt
geom="h2o.xyz"
```

For a keyword with no arguments, empty parentheses are optional: `opt` and
`opt()` mean the same thing. The examples use the shorter form when possible.
An `.xyz` or `.pdb` geometry may also immediately follow the route in a
one-line shorthand, `mrsf/bhhlyp/6-31g* h2o.xyz opt`; the examples put an
explicit `geom` last.

Use quotes for a path containing spaces, for example
`geom="structures/water molecule.xyz"`.

Read the input from top to bottom:

| Text | Meaning |
| --- | --- |
| `mrsf` | Use MRSF-TDDFT. |
| `bhhlyp` | Use the BHHLYP functional. |
| `6-31g*` | Use the 6-31G* basis. |
| `h2o.xyz` | Read the molecular geometry from `h2o.xyz`. |
| `opt` | Optimize the geometry. With MRSF, an omitted state means `S0`. |

Run it as usual:

```bash
openqp h2o-opt.oqp
```

Useful command-line variations are:

```bash
openqp h2o-opt.oqp --nompi
openqp h2o-opt.oqp --silent
openqp h2o-opt.oqp --omp 8
```

`--nompi` is required for ground-state OpenMM QM/MM molecular dynamics.

To optimize a particular excited state with tighter controls, add details only
where they are needed:

```text
mrsf(nstate=3)/bhhlyp/6-31g* opt(S1,maxit=100) scf(conv=1e-08)
geom="h2o.xyz"
```

This explicitly requests three states, the physical `S1` surface, at most 100
optimization steps, and an SCF convergence threshold of `1e-8`.

### Common MRSF calculations

```text
mrsf(nstate=3)/bhhlyp/6-31g*
geom="h2o.xyz"
```

```text
mrsf(nstate=3)/bhhlyp/6-31g* grad(S1)
geom="h2o.xyz"
```

```text
mrsf(nstate=3)/bhhlyp/6-31g* opt(T0)
geom="h2o.xyz"
```

```text
mrsf/bhhlyp/6-31g* meci(S1,S2)
geom="guess.xyz"
```

### Four rules to remember

1. **Write one main task (primary driver) per file.** Choose one of `energy`,
   `grad`, `opt`, `meci`, `soc`, and the other main tasks. For example,
   `grad(S1) opt(S1)` is an error. Modifiers such as `pcm(water)` and exact
   section controls such as `scf(conv=1e-8)` do not count as another task.
2. **MRSF state labels start at zero in each spin group (manifold).** `S0`,
   `T0`, and `Q0` are the lowest singlet, triplet, and quintet states,
   respectively. `S1`, `T1`, and `Q1` are the next states. These labels
   describe ordering *within each spin group*, not a combined energy ordering
   between singlets and triplets.
3. **Do not specify the internal reference for MRSF.** Write `opt(S0)` or
   `opt(T0)`, not `mult=...`, `spin=...`, or a ROHF multiplicity. OpenQP selects
   the required high-spin working reference automatically.
4. **Existing `.inp` files still work unchanged.** Use the `.oqp` syntax in
   a new `.oqp` file. Keep the traditional sectioned syntax in `.inp`; do not
   mix the two formats in one file.

The repository's `examples` tree provides a same-stem `.oqp` companion for
every `.inp`, covering the full legacy example inventory rather than only a
small hand-written subset.

### SOC state counts

One count requests the same number of singlets and triplets:

```text
mrsf(nstate=3)/bhhlyp/6-31g* soc
geom="h2o.xyz"
```

This calculates `S0`--`S2` and `T0`--`T2`. To use different counts, put both
counts in `soc(...)`:

```text
mrsf/bhhlyp/6-31g* soc(ns=3,nt=5)
geom="h2o.xyz"
```

This calculates `S0`--`S2` and `T0`--`T4`. Do not combine route `nstate` with
`soc(ns=...,nt=...)`; choose one form or the other.

## General Form

The complete pattern is:

```text
ROUTE GLOBAL... PRIMARY_DRIVER MODIFIER... SECTION_CALL...
```

The route comes first. Whitespace outside quotes and parentheses separates the
remaining parts, so spaces and line breaks are interchangeable. The canonical
layout — what OpenQP writes when it renders an input — puts the route, the
global options, the driver, and the modifier/section calls on one line,
wrapping to a new line only when that line would exceed about 100 characters
(long multiconfigurational inputs therefore span a few lines), and moves
`geom`/`geom2` to the end. It also omits every option that merely restates a
runtime default. OpenQP converts any layout to its established configuration
and applies the same schema validation used for traditional input.

## Route and Model Reference

Use one of the following route forms. `N` is a positive number of response
roots.

| Calculation model | Canonical route |
| --- | --- |
| Hartree--Fock with automatic restricted/unrestricted reference | `hf/BASIS` |
| Explicit HF reference | `rhf/BASIS`, `uhf/BASIS`, or `rohf/BASIS` |
| Kohn--Sham DFT with automatic restricted/unrestricted reference | `dft/FUNCTIONAL/BASIS` |
| Explicit KS reference | `rks/FUNCTIONAL/BASIS`, `uks/FUNCTIONAL/BASIS`, or `roks/FUNCTIONAL/BASIS` |
| MP2 | <code>mp2(reference=rhf&#124;rohf&#124;uhf,variant=...,same_spin_scale=...,opposite_spin_scale=...)/BASIS</code> |
| CCSD | <code>ccsd(reference=rhf&#124;rohf&#124;uhf)/BASIS</code> |
| CCSD(T) | <code>ccsd_t(reference=rhf&#124;rohf&#124;uhf)/BASIS</code>; `ccsd-t` and `ccsdt` are accepted input aliases |
| TDHF | `tdhf(nstate=N)/BASIS` |
| TDDFT | `tddft(nstate=N)/FUNCTIONAL/BASIS` |
| TDA-TDDFT | `tda(nstate=N)/FUNCTIONAL/BASIS` |
| CIS / TDA-TDHF | `tda-tdhf(nstate=N)/BASIS`; `cis(...)` is an accepted input alias |
| SF-TDDFT | `sf(nstate=N)/FUNCTIONAL/BASIS` |
| MRSF-TDDFT | `mrsf(nstate=N)/FUNCTIONAL/BASIS` |
| UMRSF-TDDFT | `umrsf(nstate=N)/FUNCTIONAL/BASIS` |
| SF-TDHF | `sf-tdhf(nstate=N)/BASIS` |
| MRSF-TDHF | `mrsf-tdhf(nstate=N)/BASIS` |
| UMRSF-TDHF | `umrsf-tdhf(nstate=N)/BASIS` |

The parenthesized route options are optional; omit the parentheses when none
are needed.

The table shows the canonical spellings. The parser also accepts these
compatibility aliases: `mrsf-tddft`, `mrsftddft`, `umrsf-tddft`,
`umrsftddft`, `sf-tddft`, `sftddft`, `td-dft`, and `ks-dft`. Canonical
rendering normalizes them to the routes in the table.

The explicit-reference routes never silently change the requested reference.
`rhf` and `rks` use the default `mult=1`; writing it is unnecessary.
`rohf` and `roks` require an open-shell `mult` of at least 2. `uhf` and `uks`
may also be used deliberately with the omitted singlet default. The TDHF-family
routes contain no functional component, whereas the corresponding
TDDFT-family routes do.

All response-route parentheses accept only `nstate`. Select a physical spin
manifold with the primary-driver state, not with a route option. Thus
`mrsf(nstate=3)/... energy(T0)` requests the three triplet states `T0`--`T2`.
Numeric route `mult=...`, `multiplicity=...`, and `spin=...` are not public
spellings. Put advanced response controls in their exact section call, for
example:

```text
mrsf(nstate=5)/bhhlyp/6-31g* energy(S0) tdhf(nvdav=30)
geom="h2o.xyz"
```

MP2 route parentheses likewise accept only `reference`, `variant`,
`same_spin_scale`, and `opposite_spin_scale`; other MP2 settings use
`mp2(...)` as an exact section call. `reference=uhf`, for example, preserves a
deliberate UMP2 reference even when `mult=1`. Use physical state labels in the
driver for state-specific work.
`mrsf(...) energy` defaults to the singlet target manifold. `energy(T0)`
selects the triplet manifold, and `energy(Q0)` selects the quintet manifold for
all-electron MRSF. The high-spin working reference is supplied internally.

## Top-Level Globals

These common values follow the route without a section name:

| Keyword | Meaning |
| --- | --- |
| `geom="FILE"` | Required primary geometry. Relative paths are resolved from the `.oqp` file directory. |
| `geom2="FILE"` | Previous or second geometry for workflows that require one. |
| `charge=N` | Molecular charge; default `0`. |
| `mult=N` | Ordinary SCF-reference multiplicity for non-mixed-reference models. It is forbidden for MRSF, UMRSF, and SF routes. |
| `library=VALUE` | Basis-library mapping used by the established input schema. |
| `ispher=VALUE` | Spherical/Cartesian AO convention. |
| `perf=N` | OpenQP performance preset. |
| `omp_threads=N` | OpenMP threads per MPI rank. `threads=N` is an alias. |

Put geometry last in a readable `.oqp` file. For an external coordinate file,
use `geom="h2o.xyz"`. To keep coordinates inside the `.oqp` file, use a
triple-quoted block; each atom occupies its own line:

```text
hf/sto-3g
geom="""
O   0.000000   0.000000   0.000000
H   0.000000   0.000000   0.960000
H   0.000000   0.750000  -0.240000
"""
```

The same coordinates may stay on one source line by escaping atom separators:

```text
hf/sto-3g
energy
geom="O 0 0 0\nH 0 0 0.96\nH 0 0.75 -0.24"
```

The escaped one-line value and the triple-quoted multiline block are parsed
identically. Canonical rendering chooses the multiline block for readability.

`geometry`/`system` are accepted aliases for `geom`, and
`geometry2`/`system2` for `geom2`. `basis=...` and `functional=...` are accepted
only when they agree with the route; the slash route is the canonical spelling.
Top-level `multiplicity=N` remains an accepted input alias for canonical
`mult=N`.

For MRSF-family calculations, select the target with a physical driver label
such as `S0`, `S1`, `T0`, or `Q0`; never expose the automatic high-spin
reference through a top-level `mult`.

## Exactly One Primary Driver

Every `.oqp` calculation has exactly one primary driver. If it is omitted,
OpenQP inserts `energy`. Square brackets below mean an optional argument;
`STATE` is a physical label such as `S0`, `S1`, `T0`, or `Q0`.

In the signatures, `OPT` means the common native controls `maxit`, `rmsd_grad`,
`rmsd_step`, `max_grad`, `max_step`, `energy_shift`, `energy_gap`, and
`init_scf` used by the established optimizers. `BAEKA_OPT` is the narrower
BaekA convergence set `maxit`, `rmsd_grad`, `energy_shift`, and `init_scf`;
BaekA does not terminate on `rmsd_step`, `max_step`, or `max_grad`. `BAEKA`
means the public adaptive-penalty controls `sigma`,
`alpha`, `delta_beta`, `beta_schedule`, and `gap`. The MECI and MECP
`algorithm` options map to `meci_search` and `mecp_search`, and `gap_sigma`
scales the `auglag` gap term. The MECI `algorithm` option
lowers to the traditional
[`meci_search`](keywords/optimize.md#meci_search) selector. The BaekA controls
lower respectively to `pen_sigma`, `pen_alpha`, `pen_delta`, `pen_jump`, and
`energy_gap`; the historical multiplicative `pen_incre` key is not a BaekA
control. `TCI` means the established legacy three-state controls `pen_sigma`,
`pen_alpha`, `pen_incre`, and `gap_weight`. `ENGINE` means `coordsys`, `trust`,
`trust_max`, `auto_recovery`, `recovery_maxit`, and `recovery_trust`; these
apply to the native minimum, crossing-point, and transition-state optimizers.
MEP and IRC own their path step and gradient threshold, while NEB owns its FIRE
band controls; those drivers do not use `ENGINE`. `NAMD`
means the current `[md]` controls `nstep`,
`dt`, `active`, `substep`, `decoherence`, `edc_c`, `thrshe`, `tdc`, `trivial`,
`trivial_thresh`, `init_temp`, `velocity`, `seed`, `rng_stream`,
`first_hop_step`, `nacme_check`, `ba_gap_max`, `nacme_gate`,
`nacme_gate_invariant_tol`, `nacme_gate_abs_tol`, `nacme_gate_rel_tol`,
`nacme_gate_consecutive`, `nve_gate`, `nve_gate_abs_tol`,
`nve_gate_step_tol`, `nve_gate_transition_tol`, `nve_gate_consecutive`,
`trajectory_interval`, `restart_interval`,
`trajectory_file`, `restart_file`, `restart`, `soc`,
`soc_basis`, `soc_du_dt_corr`, `soc_tdc_grad_corr`, `grad_wthr`, `init_state`,
`econs`, `dt_adaptive`, `dt_min`, and `dx_max`.
`NEB` means the native options `product`, `images`/`nimage`, `spring`, `climb`,
`fmax`, `frms`, `climb_fmax`, `dt`/`neb_dt`, `maxmove`, `align`, `opt_ends`,
`end_fmax`, and `output`.

| Driver signature | Lowered option family |
| --- | --- |
| <code>energy([S0&#124;T0&#124;Q0])</code> | Single-point energy. On a response route, the optional zero-state label selects a spin manifold; no other state or options are accepted. MRSF defaults to singlet. |
| `grad([STATE],td_prop=...,export=...,title=...)` | Gradient target plus concise `[properties]` controls. Default target is `S0`, except that SF routes require `root=N`. |
| `opt([STATE],OPT...,ENGINE...,freeze="distance(i,j)")` | Native minimum optimization. Default target is `S0`, except that SF routes require `root=N`. `freeze` holds one or more semicolon-separated atom-pair distances at their initial values. |
| <code>meci(STATE1,STATE2,[algorithm=auglag&#124;penalty&#124;ubp&#124;hybrid],[gap_sigma=...],OPT...,ENGINE...)</code><br><code>meci(STATE1,STATE2[,STATE3...],algorithm=baeka,BAEKA_OPT...,BAEKA...,ENGINE...)</code> | Native intersection search for states of the same multiplicity. A two-state call defaults to `auglag`, the augmented Lagrangian, whose gap term is scaled by `gap_sigma` (default `10.0`); a call with three or more states selects `baeka`, the only N-state algorithm. Write `algorithm=baeka` explicitly only to select BaekA for two states; with three or more states it is implied and the canonical rendering omits it. State order is normalized. Use public `gap` rather than the internal `energy_gap` spelling in a BaekA call, and do not supply both. See [BaekA Multistate MECI](workflows/baeka-multistate-meci.md). |
| `tci(STATE1,STATE2,STATE3,OPT...,TCI...,ENGINE...)` | Backward-compatible three-state adaptive-penalty driver. It preserves the established `pen_sigma`/`pen_alpha`/multiplicative `pen_incre` behavior and is not an alias for BaekA. New work should select the intended MECI algorithm explicitly. |
| <code>mecp(STATE1,STATE2,[algorithm=auto&#124;sqp&#124;auglag&#124;penalty&#124;quad],[gap_sigma=...],OPT...,ENGINE...)</code> | Native crossing search for two states of different multiplicity. Defaults to `algorithm=auto`, which selects `sqp` with the native optimizer and `auglag` elsewhere. `sqp` solves the KKT equations of the constrained problem for the step and the multiplier together and reads no penalty parameter; `auglag` scales its gap term by `gap_sigma` (default `10.0`). `quad` is the legacy fixed-weight quadratic penalty; it settles at a residual gap of order `1/gap_weight` and is kept only to reproduce older runs. |
| <code>mep([STATE],maxit=...&#124;points=...,step=...,gtol=...)</code> | Native minimum-energy path with a path limit, step size, and gradient stopping threshold. |
| <code>ts([STATE],OPT...,ENGINE...,follow=N,hessian=model&#124;numerical&#124;analytical)</code> | Native P-RFO transition-state search. `follow` chooses the initial mode and `hessian` chooses the initial Hessian policy. |
| <code>irc([STATE],maxit=...,direction=forward&#124;backward,step=...,hessian=numerical&#124;analytical,gtol=...)</code> | Native IRC with an explicit branch, step size, Hessian type, and gradient stopping threshold. |
| `neb([STATE],maxit=...,NEB...)` | Native NEB; `product="FILE"` is required and the final band is written as a multi-frame XYZ file. |
| <code>hess([STATE],type=numerical&#124;analytical,dx=...,nproc=...,read=...,restart=...,temperature=...,clean=...)</code> | Hessian/frequency calculation. |
| `nac(STATE1,STATE2,type=numerical,dx=...,nproc=...,restart=...,clean=...,align=...)` | Numerical nonadiabatic-coupling vector between two states in the same spin manifold. |
| `bp(STATE1,STATE2,type=numerical,dx=...,nproc=...,restart=...,clean=...,align=...)` | Numerical branching-plane calculation between two states in the same spin manifold. |
| `nacme(STATE1,STATE2,dt=...,align=...)` | Coupling matrix element; requires `geom2` or `guess(file2=...)`. |
| `soc(soc_2e=...,ns=...,nt=...)` | Spin-orbit coupling; accepts no single target state. `ns` and `nt` must be supplied together. |
| `md([S0]) qmmm(...)` | Ground-state QM/MM molecular dynamics. `qmmm(...)` is mandatory and owns the OpenMM controls. |
| `namd([STATE],NAMD...)` | MRSF nonadiabatic molecular dynamics using the `[md]` controls listed above; defaults to `S1`. |
| <code>ekt([STATE],ip=true&#124;false,ea=true&#124;false)</code> | MRSF extended Koopmans IP/EA options; the parent state defaults to `S0`. |
| <code>thermo([STATE],type=numerical&#124;analytical,dx=...,nproc=...,read=...,restart=...,temperature=...,clean=...)</code> | Alias that lowers to the supported Hessian path; the state defaults to `S0`. |
| `prop([STATE],scf_prop=...,nmr_gauge=...,td_prop=...,export=...,title=...)` | MRSF-TDDFT/MRSF-TDHF property driver; defaults to `S0`. |
| `data([STATE],scf_prop=...,nmr_gauge=...,td_prop=...,export=...,title=...)` | Multi-state data/gradient export. A response route defaults to `S1`; an SF route defaults to `root=1`. |

Aliases such as `sp`, `gradient`, `optimize`, `optimization`, and `hessian` are
accepted, but the spellings in the table are the canonical generated forms.
The older `tci(STATE1,STATE2,STATE3,...)` command remains available with its
established behavior. It is deliberately separate from
`meci(STATE1,STATE2,STATE3,algorithm=baeka,...)`, so existing inputs do not
silently change algorithms.
Method and workflow availability, together with option values, remain subject
to the existing OpenQP validator.

SF state character is not known before diagonalization. Therefore every
state-specific SF request must write `root=N`; do not omit the state and do not
write `S0`, `S1`, or `T0` on an SF route.

The NAC family currently requires an MRSF route and two distinct states in the
same spin manifold. `nac` and `bp` support numerical finite differences only;
analytical NAC is rejected. `nacme` accepts only its time-step/alignment
controls and additionally requires a previous geometry through `geom2` or
`guess(file2=...)`. Branching-plane analysis is not available for
The `prop` driver is currently limited to all-electron MRSF-TDDFT/MRSF-TDHF;
use `prop(S1)` or omit the state for its `S0` default.

For MRSF SOC, the route count is applied equally to singlets and triplets:

```text
mrsf(nstate=3)/bhhlyp/6-31g* soc
geom="h2o.xyz"
```

This requests `S0`--`S2` and `T0`--`T2`. Use explicit counts when the two
manifolds need different sizes:

```text
mrsf/bhhlyp/6-31g* soc(ns=3,nt=5)
geom="h2o.xyz"
```

`ns` and `nt` are inseparable, and `soc(ns=...,nt=...)` cannot be combined with
route `nstate`. For a crossing between different multiplicities, use physical
labels directly, for example `mecp(S0,T0,maxit=100)`.

Two primary drivers never imply a sequence. This is an error:

```text
mrsf(nstate=5)/bhhlyp/6-31g*
grad(S1)
opt(S1)
geom="h2o.xyz"
```

Use separate `.oqp` files for separate calculation steps.

### Native Geometry Drivers

Concise `.oqp` geometry and reaction-path drivers always use the native OpenQP
engine. There is no backend selector in this format, so do not write `lib`,
`optimizer`, `step_size`, `step_tol`, or `mep_maxit`. Traditional sectioned
`.inp` files and the Python workflow API retain their existing optional
geomeTRIC and SciPy backends for compatibility; see
[Legacy `.inp`](input-file.md) when that escape hatch is required.

MEP uses `points` as the path limit, `step` as the native path step, and `gtol`
as its gradient stopping threshold. Native TS accepts `hessian=model`,
`numerical`, or `analytical`; `model` is the fast default, while the other
choices calculate a real Cartesian initial Hessian. `follow=N` selects a
non-negative initial P-RFO mode index. Native IRC accepts only numerical or
analytical Hessians and lowers its branch, step, and `gtol` controls to the
native IRC engine.
Native minima also accept `freeze="distance(1,2);distance(2,3)"`; atom indices
are one-based and the current constraint surface is limited to frozen
distances.

Native NEB keeps backend details out of the command:

```text
dft/pbe0/def2-svp neb(S0,product="product.xyz",images=7,spring=0.08,frms=0.001,output="path.xyz")
geom="reactant.xyz"
```

`climb`, `align`, and `opt_ends` are booleans. `fmax` and `frms` are the maximum
and RMS force thresholds, and `dt` is the concise spelling of the native FIRE
time step. Unless `output` is supplied, OpenQP writes `<project>_neb.xyz` in the
log directory. The file is a multi-frame XYZ trajectory containing every final
image and its energy in Hartree.

## Modifiers

Modifiers may accompany the one primary driver:

| Modifier | Meaning |
| --- | --- |
| `d4` | Enable DFT-D4. `d4()` uses functional defaults; `d4(s6=...,s8=...,s9=...,a1=...,a2=...,alp=...)` supplies a complete explicit rational-damping set. |
| `pcm([SOLVENT], pcm...)` | Enable PCM and optionally name a solvent, for example `pcm(water)`. The current production path is an RHF/ROHF, ddX, reference-SCF single-point energy. |
| `nmr([gauge=cgo|giao])` | Request NMR shielding. Bare `nmr` defaults to GIAO. |
| `ir` | Record that IR intensities are requested; valid only with `hess(...)` or `thermo()`. |
| `raman` | Record that Raman activities are requested; valid only with `hess(...)` or `thermo()`. |
| `qmmm(qmmm...)` | Supply QM/MM options and enable `qmmm_flag` automatically. It may accompany `energy`, `md`, or `namd`; `md` requires it. |

The compatibility spellings `d4=true` and `qmmm=true` remain accepted, but
`d4` and `qmmm(...)` are preferred in canonical files.

DFT-D4 receives the molecular `charge` from the route or `[input]` section.
Explicit damping requires all six named values; partial and non-finite sets are
rejected before the calculation starts.

In canonical `.oqp` input, `nmr` explicitly lowers to the GIAO gauge; use
`nmr(gauge=cgo)` to request CGO. This canonical default does not alter the
unchanged defaults of a traditional sectioned `.inp` file.

Use `qmmm(n_steps=N,...)` for the QM/MM MD step count. `qmmm.n_steps` is the
first-class schema key used by the OpenMM MD engine. Concise `.oqp` also accepts
`qmmm(nsteps=N,...)` as a compatibility alias and lowers it to `n_steps`; do not
write both spellings. In a traditional sectioned `.inp`, `[qmmm] nsteps` remains
the separate legacy static-driver key.

`energy qmmm(...)` selects the active QM/MM single-point path. `md` is the
ground-state QM/MM molecular-dynamics driver and therefore requires
`qmmm(...)`. MRSF `namd(...)` may run gas phase or with `qmmm(...)`. Attaching
`qmmm(...)` to `grad`, `opt`, or another driver is rejected because those
generic backends do not yet provide a verified QM/MM gradient assembly.
QM/MM/OpenMM controls belong in `qmmm(...)`; nonadiabatic-dynamics controls
belong in `namd(...)` and lower to the legacy `[md]` section.

## Advanced Section Calls

Non-owned legacy keywords keep their section and option names. A section becomes
a keyword-only call:

```text
scf(conv=1e-8,maxit=100)
cc(nfzc=1,conv=1e-7,maxit=60,ndiis=8)
properties(scf_prop=mulliken)
tdhf(nvdav=30,zvconv=1e-7)
```

Advanced exact calls include non-driver sections such as `input`, `guess`,
`scf`, `mp2`, `cc`, `dftgrid`, `tdhf`, `properties`, `pcm`, `symmetry`, `json`,
and `tests`. `oqp(...)` is not a normal advanced call: it is accepted only as a
read-time compatibility form for older concise inputs and is folded into the
primary geometry driver. When a schema section is represented by a primary
driver, put its public controls in that driver instead. The established schema
remains authoritative for keyword spelling, type, allowed values, and
cross-section constraints.

The backend selectors in `[optimize]` and the entire `[geometric]` section are
available only in traditional sectioned `.inp` files and the compatible Python
API. They are deliberately rejected as exact calls in concise `.oqp` input.

This is not a promise that bookkeeping keys can be repeated verbatim. Method,
state, spin, and workflow selectors are owned by the route and primary driver;
the reserved-key table below is authoritative. In particular, use
`soc(ns=...,nt=...)`, not `tdhf(nstate_s=...,nstate_t=...)`, in `.oqp`, and do
not combine those two spellings. Likewise, put `soc_2e` in `soc(...)` rather
than repeating it through `input(soc_2e=...)`.

When a section name is also a primary driver, put its options in that driver.
For example, write `opt(S0,maxit=100,coordsys=tric,trust=0.2)`, not
`opt(S0) optimize(maxit=100)`. Older `oqp(...)` calls are accepted only as a
compatibility form and are rewritten into the primary driver when OpenQP
renders canonical input; do not write them in new `.oqp` files.

## Physical States and Reserved Internal Keys

State labels always describe the physical target:

| Model and label | Physical meaning | Internal lowering |
| --- | --- | --- |
| HF/DFT `S0` | Ground-state SCF surface | State index 0 |
| Conventional TDDFT/TDHF `S0` | Ground-state DFT/HF surface for a state-specific derivative driver | State index 0 |
| Conventional TDDFT/TDHF `Sn`, `n >= 1` | Singlet excited state `n` | Response root `n` |
| Conventional TDDFT/TDHF `Tn`, `n >= 0` | Triplet response state `n` | Response root `n + 1`, target multiplicity 3 |
| SF family `root=N` | Response root whose spin character is assigned after diagonalization | Response root `N` |
| All-electron MRSF `Sn`, `Tn`, or `Qn`, `n >= 0` | State `n` within the selected singlet, triplet, or quintet manifold | Response root `n + 1`, target multiplicity 1, 3, or 5 |

Every available MRSF spin manifold is therefore zero-based. In all-electron
MRSF, `S0`, `T0`, and `Q0` all lower to internal response root 1, while `S1`,
`T1`, and `Q1` lower to root 2.
`opt(S0)` therefore lowers to internal `istate=1`. On non-SF routes, omitted
states for `grad`, `opt`, `hess`, and similar state-specific drivers normally
default to `S0`. SF routes are the exception and require an explicit `root=N`.
OpenQP also supplies the triplet ROHF reference for MRSF and the triplet UHF
reference for UMRSF. These are implementation references, not requested
triplet surfaces.
MRSF NAMD uses exactly the same zero-based physical labels. For example,
`namd(T0,soc=true)` lowers to legacy `[md] init_state=T0`; the first triplet is
`T0`, not `T1`. A non-SOC `namd(T0)` selects internal active root 1, while an
omitted NAMD state defaults to `S1`. Do not write the internal `active` or
`init_state` selector alongside a physical driver state.

SF-family state character cannot be assigned safely before diagonalization, so
SF-TDDFT and SF-TDHF state-specific drivers require `root=N`
instead of an `S` or `T` label.
Conventional TD calculations accept `S` and `T` labels but reject `Q` labels;

The following bookkeeping keys are reserved and cannot contradict the route or
physical state syntax:

| Legacy/internal keys | Canonical source of truth |
| --- | --- |
| `input.method`, `input.runtype`, `input.system`, `input.system2`, `input.basis`, `input.functional`, `input.charge` | Route, globals, and primary driver |
| `scf.type`, `scf.multiplicity` | Route-selected reference for every model; MRSF/SF/UMRSF choose their high-spin reference automatically |
| `tdhf.type`, `tdhf.multiplicity`, `tdhf.nstate`, `tdhf.target` | Response route `nstate` and physical state labels |
| `properties.grad` | `grad(STATE)`, `prop(STATE)`, or `data(STATE)` |
| `optimize.states`, `optimize.istate/jstate/kstate`, `optimize.imult/jmult` | `opt`, `meci`, `mecp`, the legacy `tci` driver, and reaction-path state labels |
| `hess.state`, `nac.states`, `md.active/init_state` | Corresponding state-aware driver |
| `qmmm.istate` | No canonical equivalent; obsolete disconnected-path selector |
| `tdhf.nstate_s/nstate_t`, `input.soc_2e` | `soc(ns=...,nt=...,soc_2e=...)`; do not combine the two sources |

An attempt such as `mrsf/... opt(S0) scf(multiplicity=1)` is rejected instead
of allowing the last spelling to win.

`qmmm.istate` remains in the legacy schema for compatibility, but it is
an obsolete selector from the disconnected legacy `libopenmm` path. It is
reserved in `.oqp`. Ground-state `md` uses `S0`; use physical driver labels
where a state-aware canonical workflow supports them.

## Canonical Examples

1. HF energy with a tighter SCF threshold:

    ```text
    hf/6-31g* scf(conv=1e-08)
    geom="h2o.xyz"
    ```

2. DFT-D4/PCM single point:

    ```text
    dft/pbe0/def2-svp d4 pcm(water)
    geom="h2o.xyz"
    ```

3. MRSF gradient on the physical first singlet excited state:

    ```text
    mrsf(nstate=5)/bhhlyp/6-31g* grad(S1)
    geom="h2o.xyz"
    ```

4. MRSF optimization on the physical singlet ground state:

    ```text
    mrsf(nstate=5)/bhhlyp/6-31g* opt
    geom="h2o.xyz"
    ```

5. BaekA multistate MRSF conical-intersection search:

    ```text
    mrsf(nstate=5)/bhhlyp/6-31g* meci(S0,S1,S2)
    geom="guess.xyz"
    ```

6. Ground-state QM/MM molecular dynamics:

    ```text
    dft/pbe0/def2-svp md
    qmmm(pdb_file="system.pdb",forcefield_files="amber14-all.xml",qm_atoms="0-2",n_steps=100)
    geom="qm.xyz"
    ```

7. QM/MM single-point energy with the QM atom selection in the PDB geometry:

    ```text
    dft/pbe0/def2-svp qmmm()
    geom="ala.pdb 9 10 17 18 19"
    ```

8. Ground-state NEB calculation:

    ```text
    dft/pbe0/def2-svp neb(S0,product="product.xyz",images=7)
    geom="reactant.xyz"
    ```

9. NMR shielding with the canonical GIAO default:

    ```text
    dft/pbe0/def2-svp nmr
    geom="h2o.xyz"
    ```

10. Analytic Hessian with IR and Raman intent:

    ```text
    dft/pbe0/def2-svp hess(S0,type=analytical) ir raman
    geom="h2o.xyz"
    ```

11. Forward native IRC:

    ```text
    dft/pbe0/def2-svp irc(S0,hessian=analytical)
    geom="ts.xyz"
    ```

Run any canonical file normally:

```bash
openqp calculation.oqp
```

## Errors and Corrections

OpenQP reports canonical-looking mistakes as errors; it does not reinterpret
them as prose. Common corrections are:

| Error | Correction |
| --- | --- |
| `# mrsf/...` | Remove the leading `#`. It is not a route marker in `.oqp`. |
| `[input]` or `[scf]` in an `.oqp` file | Keep sectioned input in a `.inp` file, or rewrite it in compact `.oqp` syntax. |
| `grad(S1) opt(S1)` | Put the two primary calculations in separate `.oqp` files. |
| `mrsf/... mult=3` | Remove `mult`; choose the physical target in the driver, such as `opt(T0)`. |
| `sf/... grad(S1)` | Use an unlabeled response root, for example `grad(root=1)`. |
| `opt(S0,lib=oqp)` | Remove `lib`; concise geometry drivers select the native engine automatically. |
| `geometric(...)` or `lib=geometric` | Use a traditional sectioned `.inp` file and install the optional geomeTRIC extra. |
| `optimizer`, `step_size`, `step_tol`, or `mep_maxit` | These are legacy SciPy controls for traditional `.inp` or Python workflows, not concise `.oqp`. |
| NEB `k`, `maxg`, `avgg`, or `optep` | Use native `spring`, `fmax`, `frms`, or `opt_ends`; semantics and units differ, so do not copy values mechanically. |
| route `nstate` together with `soc(ns=...,nt=...)` | Choose equal counts with route `nstate`, or unequal counts with `ns` and `nt`. |
| `qmmm(istate=...)` | The key belongs to an obsolete disconnected path; use the physical-state driver contract instead. |

### Optional correction assistant

Canonical syntax is the authoritative runnable input. As a secondary aid,
OpenQP can translate a short Korean or English request into canonical syntax.
It never executes the prose directly: `request.oqp` produces
`request.resolved.oqp`, then that generated file is reparsed and validated
before execution. The original is not overwritten, and ambiguous requests are
rejected. Canonical-looking text with a syntax error is reported as a syntax
error rather than reinterpreted as prose.

For reproducible production work, inspect and retain the resolved canonical
file or write the canonical form directly.

## MRSF Log Presentation

MRSF logs present the physical request first and implementation bookkeeping
second. For example, an MRSF ground-state optimization begins with a summary
equivalent to:

```text
Method:                       MRSF-TDDFT
Calculation:                  Geometry optimization
Physical target state(s):     S0
Target spin:                  singlet
Reference:                    triplet ROHF (internal working reference)
State labels:                 physical labels shown; engine root numbers are internal
```

The SCF detail block likewise says `reference type (internal)` and `reference
multiplicity (internal)` instead of presenting the ROHF triplet as the target
state. Gradient, optimization, intersection, and dynamics progress headings
use `S0`, `S1`, `T0`, or `Q0`; an engine root number is shown only as internal
bookkeeping when needed. SOC summaries show both requested ranges, such as
`S0-S2 and T0-T2`.

The same presentation is used for existing `.inp` calculations; their input
syntax is not changed. In a traditional MRSF input, `[scf] multiplicity=3`
still describes the internal high-spin reference, while `[tdhf]
multiplicity=1|3|5` selects the physical response manifold.
