# `[md]`

The `[md]` section controls nonadiabatic molecular dynamics (`runtype=namd`):
Tully fewest-switches surface hopping (FSSH) on MRSF states, optionally
with spin-orbit-coupled intersystem crossing (SOC-NAMD) and ESPF QM/MM
embedding. It is used together with [`[input] runtype=namd`](input.md#runtype)
and an all-electron MRSF-TDDFT theory block (`method=tdhf`,
`[tdhf] type=mrsf`). For embedded dynamics, use the [`[qmmm]`](qmmm.md)
section. See the
[SOC-NAMD-QMMM workflow](../workflows/soc-namd-qmmm.md) for complete decks and
theory.

!!! note "Available in OpenQP 1.3.0"
    NAMD, validation gates, packed trajectory/restart records, and the SOC
    trajectory/restart extensions are included in OpenQP 1.3.0.

## Background

Surface-hopping dynamics propagates classical nuclei on one active
Born-Oppenheimer (or spin-adiabatic) potential energy surface while the
electronic amplitudes evolve; stochastic hops between surfaces reproduce
nonadiabatic transitions. OpenQP implements FSSH on MRSF-TDDFT states with
energy-based decoherence (EDC), time-derivative couplings, and trivial-crossing
following. Enabling `soc=true` extends the dynamics to the spin-adiabatic
manifold so intersystem crossing (ISC) between singlet and triplet MRSF states
is described (SOC-NAMD). Enabling [`[input] qmmm_flag=true`](input.md#qmmm_flag)
embeds the MRSF-TDDFT QM region in an OpenMM MM environment via the ESPF
operator.

## Minimal NAMD Example

Gas-phase FSSH on MRSF-TDDFT states:

`.oqp`:

```text
mrsf(nstate=5)/bhhlyp/6-31g* namd(S1)
geom="molecule.xyz"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("molecule_namd")
job.molecule("molecule.xyz")
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=5)
job.workflow.namd(init_state="S1")
mol = job.run()
```

Legacy `.inp`:

```ini
[input]
runtype    = namd
method     = tdhf
functional = bhhlyp
basis      = 6-31g*
system     = molecule.xyz

[scf]
type         = rohf
multiplicity = 3

[tdhf]
type   = mrsf
nstate = 5

[md]
active    = 2
nstep     = 200
dt        = 0.5
init_temp = 300.0
```

## Core Dynamics Keywords

### `nstep`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `100` |
| Used by | nuclear propagation |

Number of nuclear (velocity-Verlet) steps.

### `dt`

| Field | Value |
| --- | --- |
| Type | float (fs) |
| Default | `0.5` |
| Used by | nuclear propagation |

Nuclear timestep in femtoseconds.

### `active`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | initial active surface |

Initial active state (1-based). For plain FSSH this indexes the MRSF states
(`1 <= active <= [tdhf] nstate`). For SOC-NAMD it indexes the spin-adiabatic
manifold (`1 <= active <= ns + 3*nt`; see [`soc`](#soc)). For SOC runs,
[`init_state`](#init_state) can override `active` by MCH character.

### `substep`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `200` |
| Used by | electronic propagation |

Number of electronic sub-steps integrated per nuclear step.

### `decoherence`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `edc` |
| Values | `edc`, `off` |
| Used by | electronic propagation |

Decoherence correction. `edc` applies the energy-based decoherence correction
(EDC) of Granucci & Persico (see
[References](../references.md#nonadiabatic-dynamics)), the SHARC default; `off`
disables it. Energy-based decoherence is recommended for surface hopping.

### `edc_c`

| Field | Value |
| --- | --- |
| Type | float (Ha) |
| Default | `0.1` |
| Used by | EDC decoherence |

The EDC constant `C` (in Hartree) in the energy-based decoherence rate. Only
used when `decoherence=edc`.

### `thrshe`

| Field | Value |
| --- | --- |
| Type | float (Ha) |
| Default | `0.1` |
| Used by | hop gating |

Energy-gap gate for hops: a hop is blocked when the state gap exceeds `thrshe`.
The `0.1` Ha default applies to both same-spin and SOC dynamics and blocks
large-gap transitions outside the intended local crossing region.

### `tdc`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `fd` |
| Values | `fd`, `baeck_an` (`npi` pending) |
| Used by | time-derivative couplings |

Time-derivative coupling scheme. `fd` uses the finite-difference
(Hammes-Schiffer / Tully) overlap form. The norm-preserving interpolation
(`npi`) option is pending.

`baeck_an` propagates the electronic amplitudes with the lagged
time-dependent Baeck-An (TD-BA) coupling instead of the overlap TDC. From the
adiabatic gap `dE_ij(t)`, the magnitude at the centre of three consecutive
energy points is

\[
\left|\tau^{\mathrm{BA}}_{ij}(t_n)\right|
=\frac{1}{2}\sqrt{
\frac{\mathrm d^2 \Delta E_{ij}(t_n)/\mathrm dt^2}
     {\Delta E_{ij}(t_n)}}.
\]

A pair is set to zero when the radicand is nonpositive, the centre gap is zero,
or the centre gap exceeds [`ba_gap_max`](#ba_gap_max).

**The coupling is one nuclear step lagged.** The centred curvature at `t_n`
becomes available only after the energies at `t_(n+1)` have been evaluated, so
causal dynamics applies it during the electronic propagation at `t_(n+1)`. It
is not an instantaneous Baeck-An coupling.

Energies alone do not fix the wavefunction gauge, so the magnitude takes the
pairwise sign of the phase-tracked, centred overlap TDC; a pair whose overlap
sign is exactly indeterminate is set to zero. The first interval, and any
interval following a discontinuous history, fall back to the overlap NPI TDC as
a warm-up value. Dense trajectory records distinguish the two sources:

| `tdc_source` | Meaning |
| --- | --- |
| `1` | Overlap NPI warm-up |
| `3` | Lagged Baeck-An magnitude with overlap-transported sign |

Scope: same-spin FSSH only. It supplies no `3N` coupling vector and therefore no
direction-specific momentum adjustment, and it makes no claim to supply a Berry
phase or a signed electronic-state gauge. The one-step lag and the `ba_gap_max`
pair selection are part of the named approximation and must be preserved in
comparisons and restart signatures. Pair it with `rescale=isotropic` to test the
approximate electronic coupling without adding an analytic NAC-vector
calculation at a hop.

### `trivial`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | trivial-crossing handling |

Enable trivial- (weakly avoided) crossing detection and diabatic following, so
the active surface tracks state character through sharp crossings instead of
hopping. This is an opt-in heuristic rather than part of standard FSSH; leave
it off unless the chosen protocol has been validated with it.

### `trivial_thresh`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.5` |
| Used by | trivial-crossing handling |

State-overlap threshold that flags a trivial crossing. Only used when
`trivial=True`.

## Initial Conditions

### `init_temp`

| Field | Value |
| --- | --- |
| Type | float (K) |
| Default | `300.0` |
| Used by | initial velocities |

Temperature for Maxwell-Boltzmann initial velocities (used when
`velocity=maxwell`).

### `velocity`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `maxwell` |
| Values | `maxwell`, `zero`, *(file path)* |
| Used by | initial velocities |

Initial velocity source: `maxwell` samples a Maxwell-Boltzmann distribution at
`init_temp`, `zero` starts from rest, or a file path reads velocities from a
file. `maxwell` is a classical distribution, not a vibrational Wigner sample.
One `.oqp` NAMD request uses one initial geometry; a Wigner ensemble must supply
one independently sampled geometry/velocity pair per trajectory.

### `seed`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` (resolved once to the local date as `YYYYMMDD`) |
| Used by | trajectory counter-RNG |

Seed for the resident Fortran counter-RNG that draws Maxwell initial velocities
and hopping random numbers. The zero sentinel is resolved once when a run starts,
and the resulting integer is frozen in the generated restart manifest. Set an
explicit nonzero campaign seed for reproducible ensembles. A hopping
draw is a pure function of `(seed, rng_stream, physical MD step)`, so worker
scheduling and unrelated calls cannot shift the hopping sequence.

### `rng_stream`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | initial velocities and hopping counter-RNG |

Non-negative trajectory stream identifier. Use a distinct value for every
trajectory in one ensemble while keeping `seed` fixed as the campaign seed.
It separates both Maxwell initial velocities and hopping draws. The same
`(seed, rng_stream, step)` triple always gives the same full-precision uniform
value, which permits exact two-code replay. Do not reuse one stream for two
nominally independent trajectories.

### `first_hop_step`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | active-state transitions and hopping RNG |

First physical nuclear step at which an active-state transition is permitted.
Step 0 is the initial electronic structure; step 1 has both endpoint structures
and therefore defines the first overlap/TLF2 interval. Electronic coefficients
and hop probabilities are propagated at every such interval, beginning at step
1. With the default `first_hop_step=1`, the first FSSH decision is also made at
step 1. If `first_hop_step=2` is selected explicitly, step 1 still propagates
the electronic coefficients but cannot change the active state, rescale the
velocity, or consume a hopping random number.

For a strict OpenQP/KNU comparison, use the same full-precision random tape and
the same `first_hop_step`; note that a code which labels the initial structure
as step 1 may call OpenQP's first interval “step 2.” Rounded values copied from
ordinary text output can change a hop when the probability lies close to the
random threshold.

### `nacme_check`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `baeck_an` for same-spin NAMD; contextually `off` for SOC-NAMD |
| Values | `off`, `baeck_an` |
| Used by | independent NACME validation |

Enable an energy-only time-dependent Baeck–An (TD-BA) diagnostic alongside the
overlap/TLF coupling used by NAMD. With `baeck_an`, OpenQP uses three consecutive
energy points to evaluate the nonuniform central curvature of every energy gap.
It compares the resulting TD-BA coupling magnitude with the overlap-derived TDC
interpolated to the same central time and logs both matrices plus RMS and maximum
magnitude errors.

TD-BA does not use wavefunctions and therefore cannot check MO/root phase or the
signed NACME gauge. Treat it as an independent check of coupling magnitude and
peak location, not as an oracle or a replacement for TLF. It is based on a
two-state near-crossing approximation and may overestimate couplings, especially
outside its intended region. The current implementation supports same-spin
NAMD. SOC-NAMD records its full complex spin-adiabatic overlap and anti-Hermitian
TDC instead, so its inherited default is disabled contextually. An explicit
non-`off` request through the Python workflow API is rejected rather than
silently ignored. See the
[Baeck-An references](../references.md#nonadiabatic-dynamics).

### `ba_gap_max`

| Field | Value |
| --- | --- |
| Type | float (Ha) |
| Default | `0.0734986443513` (2 eV) |
| Used by | TD-Baeck–An NACME validation |

Maximum central energy gap included in the TD-BA diagnostic. Pairs above this
gap, or pairs without a positive TD-BA curvature radicand, are not evaluated by
the reference-comparison gate.

With `nacme_check=baeck_an` the setting is diagnostic only and does not modify
the overlap/TLF coupling, electronic propagation, or hopping probabilities. With
[`tdc=baeck_an`](#tdc) the same selection governs the coupling that is actually
propagated, so a pair above this gap contributes no electronic coupling at all.

### `nacme_gate`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `off` |
| Values | `off`, `warn`, `error` |
| Used by | MD NACME validation policy |

Policy applied to the common resident-Fortran NACME gate. The gate always
reports matrix invariants and reference-comparison metrics when a check is
enabled. `off` records diagnostics only, `warn` logs failed gates without
stopping dynamics, and `error` stops immediately for a finite-value or matrix
invariant failure and stops after `nacme_gate_consecutive` consecutive reference
failures.

The exact invariants are a zero diagonal and antisymmetry of both the MD TDC and
the supplied reference. TD-BA is compared by magnitude. A future phase-aligned
analytic NAC reference can use the same gate in signed mode after contracting
the analytic vector with the nuclear velocity, `d_IJ . v`, at the matching time.

### `nacme_gate_invariant_tol`

| Field | Value |
| --- | --- |
| Type | float (au^-1) |
| Default | `1.0e-10` |
| Used by | diagonal and antisymmetry checks |

Absolute tolerance for exact NACME matrix invariants.

### `nacme_gate_abs_tol`

| Field | Value |
| --- | --- |
| Type | float (au^-1) |
| Default | `1.0e-4` |
| Used by | reference comparison |

Absolute component of the pair acceptance threshold.

### `nacme_gate_rel_tol`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `1.0` |
| Used by | reference comparison |

Relative component of the pair acceptance threshold. Pair `IJ` passes when
`error <= nacme_gate_abs_tol + nacme_gate_rel_tol * abs(reference_IJ)`.
The deliberately broad default reflects that TD-BA is an approximation; choose
thresholds from a validated system before using `nacme_gate=error` for a
production campaign.

### `nacme_gate_consecutive`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `3` |
| Used by | `nacme_gate=error` |

Number of consecutive time points with at least one failed reference pair
required before aborting. A passing point resets the count. Exact invariant or
non-finite failures are not delayed.

### `nve_gate`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `warn` for NAMD |
| Values | `off`, `warn`, `error` |
| Used by | same-spin NVE/FSSH energy validation |

Validate the nominally microcanonical gas-phase or QM/MM trajectory, including
same-spin and SOC drivers.
The driver records total-energy drift from step zero, the change from the
previous step, the energy discontinuity at a successful hop or trivial state
change, and drift per femtosecond. `warn` prints the NVE table without stopping;
`error` records the failing point and then aborts for a failed transition-energy
check or after `nve_gate_consecutive` consecutive drift/step failures. The
restart checkpoint is not advanced past the rejected point.

This is a quantum-classical FSSH energy validation, not a claim that the
electronic subsystem alone is microcanonical. Decoherence, frustrated hops,
finite time steps, and electronic-structure convergence can contribute to
drift. NAMD drivers always own their velocity-Verlet propagation and do not
inherit `[qmmm] ensemble`; therefore the default remains `warn` even if that
unrelated ground-state QM/MM setting says `NVT` or `NPT`. The current
surface-hopping QM/MM driver uses OpenMM for forces but
performs its own velocity-Verlet + SHAKE/RATTLE propagation; the ground-state
[`[qmmm] ensemble`](qmmm.md#ensemble) NVT/NPT integrators do not control NAMD.

### `nve_gate_abs_tol`

| Field | Value |
| --- | --- |
| Type | float (Ha) |
| Default | `5.0e-3` |
| Used by | total drift from the initial energy |

### `nve_gate_step_tol`

| Field | Value |
| --- | --- |
| Type | float (Ha) |
| Default | `1.0e-3` |
| Used by | change in total energy between adjacent integration steps |

The comparison is made at every physical MD step, not between saved trajectory
records — `trajectory_interval` does not widen it. Calibrate the tolerance
against one integration step, or an `error` gate will be far more permissive
than intended once the automatic ~10 fs write cadence is resolved.

### `nve_gate_transition_tol`

| Field | Value |
| --- | --- |
| Type | float (Ha) |
| Default | `1.0e-6` |
| Used by | same-geometry energy discontinuity across a state change |

This local quantity is evaluated immediately before and after hop velocity
rescaling (or trivial state following), so it is a stricter check than ordinary
integrator drift.

### `nve_gate_consecutive`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `3` |
| Used by | `nve_gate=error` drift/step policy |

### `trajectory_interval`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` (automatic, approximately 10 fs) |
| Used by | dense NAMD trajectory |

Positive values write every Nth MD step to the dense binary trajectory. Zero
chooses `round(10 fs / dt)`, with a minimum of one step. The final point and a
point that triggers either strict NACME or NVE validation are written even when
they are not on the regular interval.

### `trajectory_file`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `<project>.namd.trj` |
| Used by | dense NAMD trajectory |

Appendable, packed fixed-record binary trajectory for machine analysis. Numeric
values are stored directly rather than expanded as repeated decimal text, so it
is substantially more compact and faster to scan than a text trajectory holding
the same matrices. Every record contains coordinates, velocities, energies,
populations, complex electronic coefficients, hop decision and full-precision
random value, state overlap, overlap TDC, the active reference TDC/mask, gate
metrics, and root/phase tracking order, phase, matched overlap, and margin. It
also stores the NVE drift, step change, transition jump, drift rate, verdict,
and failure streak. SOC records additionally retain the complex spin-adiabatic
overlap and anti-Hermitian TDC as real/imaginary components plus the active
representation (`adiabatic` or `mch`). It is not an NPZ/ZIP archive: those formats cannot be
appended safely without rewriting the complete trajectory and cannot be
memory-mapped record by record. Compress a completed trajectory only as an
archival/post-processing step.

Read it without loading the complete trajectory into memory:

```python
from oqp.library.namd import read_namd_trajectory

header, trajectory = read_namd_trajectory("job.namd.trj")
print(trajectory["active"])
print(trajectory["overlap_tdc_au"][:, 0, 1])
```

The human-readable NACME validation table is written to the main log; no
separate tabulated NACME file is produced. On restart, OpenQP validates the
trajectory schema and calculation identity,
removes only records newer than the last atomic checkpoint, and appends the
continued trajectory. Committed records are not rewritten.

### `restart_interval`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` (automatic, approximately 10 fs) |
| Used by | atomic NAMD checkpoint |

Positive values write the restart checkpoint every Nth step. Zero chooses
`round(10 fs / dt)`, with a minimum of one step. The final MD step is always
saved. Increasing the effective interval reduces checkpoint I/O but also
increases the maximum amount of accepted dynamics that must be recomputed after
an error. The generated restart input is written once because its contents do
not change between checkpoints.

### `restart_file`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `<project>.namd.restart.npz` |
| Used by | exact NAMD continuation |

Compressed, non-pickle numerical checkpoint containing coordinates, velocities,
acceleration, electronic coefficients, the previous electronic-structure tag
bundle required by MO/root/phase tracking, counter-RNG identity, NACME gate
streak, and TD-BA three-point history. SOC checkpoints additionally contain the
previous SOC eigensystem and singlet/triplet response vectors needed for the
next spin-adiabatic overlap. Their exact shapes, dtypes, and finite values are
validated before restoration. The file is written to a temporary file, flushed,
and atomically replaced.

### `restart`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | trajectory restart |

Restart the trajectory from a saved state.

Every canonical `.oqp` NAMD run also writes a directly runnable file named
`<project>.namd.restart.oqp` beside the main log. It preserves the original
request, resolves input-owned paths against the original input directory, and
adds `restart=true` plus explicit checkpoint and trajectory paths and freezes a
date-derived seed. Run this manifest to continue toward the original final `nstep`;
`nstep` is not interpreted as an additional number of steps. The manifest and
its NPZ numerical checkpoint must remain together. Deriving the manifest name
from the project/log stem prevents simultaneous trajectories in one output
directory from overwriting each other.

Restart, packed trajectory/checkpoint output, and NVE gating support all four
same-spin/SOC and gas-phase/QM/MM driver combinations. The independent TD-BA
NACME comparison remains same-spin only; SOC stores its complex overlap/TDC but
does not reinterpret TD-BA as a spin-adiabatic reference.

## SOC-NAMD (Intersystem Crossing)

### `soc`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | SOC-NAMD dispatch |

Enable SOC-NAMD: surface hopping on the **spin-adiabatic** manifold so
intersystem crossing between MRSF singlet and triplet states is described. When
`soc=true`, the manifold has `ns + 3*nt` states (`ns` singlets and `nt` triplets,
each triplet contributing three `Ms` sublevels, with `ns = nt = [tdhf] nstate`).
Combined with [`[input] qmmm_flag=true`](input.md#qmmm_flag), this selects an
SOC-QM/MM driver; [`soc_basis`](#soc_basis) chooses between the spin-adiabatic
and MCH-basis variants (see the
[dispatch table](../workflows/soc-namd-qmmm.md#how-the-driver-is-selected)).

### `soc_basis`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `adiabatic` |
| Values | `adiabatic`, `mch` |
| Used by | SOC-NAMD propagation and force basis |

Selects the SOC-NAMD representation.

| Value | Meaning |
| --- | --- |
| `adiabatic` | Propagate on spin-adiabatic SOC eigenstates and use the weighted-MCH diagonal gradient controlled by [`grad_wthr`](#grad_wthr). |
| `mch` | Propagate in the spin-pure MCH basis with exact active-root MCH gradients. With QM/MM, this selects `NAMD_SOC_MCH_QMMM`. |

The `mch` basis is the recommended production mode from the current validation
work because it avoids the approximate weighted-gradient force used by the
spin-adiabatic path.

### `soc_du_dt_corr`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | spin-adiabatic SOC-NAMD force correction |

For `soc_basis=adiabatic`, add a finite-difference `dU/dt` force correction to
the weighted-MCH diagonal gradient. This is a diagnostic/validation option for
the spin-adiabatic force path and is ignored by the MCH-basis driver.

### `soc_tdc_grad_corr`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | spin-adiabatic SOC-NAMD force correction |

For `soc_basis=adiabatic`, add an approximate MCH time-derivative-coupling
projected gradient correction. It can be combined with
[`soc_du_dt_corr`](#soc_du_dt_corr) for force-basis testing, and is ignored by
the MCH-basis driver.

### `grad_wthr`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.001` |
| Used by | SOC-NAMD active-surface force |

Weight threshold for the spin-adiabatic weighted-MCH diagonal gradient. Only
the spin-pure (MCH) components whose weight in the active spin-adiabatic state
exceeds `grad_wthr` contribute to the active-surface force (the three `Ms`
sublevels of a triplet share a summed weight). A small value keeps the force
continuous across regions of strong spin mixing.

### `init_state`

| Field | Value |
| --- | --- |
| Type | string |
| Default | *(empty)* |
| Used by | SOC-NAMD initial surface |

Start SOC-NAMD on the spin-adiabatic state whose dominant character matches this
MCH label (`S0`, `S1`, `T0`, `T1`, ...). MRSF labels are zero-based within
each spin manifold in concise `.oqp`, so `T0` is the first triplet and `T1` is
the second. Historical sectioned `.inp` and Python configurations used `T1`
for the first triplet; that spelling remains a compatibility alias there, and
`T0` is also accepted unambiguously. When empty, the initial surface is taken
from the [`active`](#active) index. `init_state` overrides `active` for SOC runs.

### `econs`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | SOC-NAMD energy conservation |

Rescale velocities each step to conserve the total energy. This is a temporary
stabilizer (band-aid) for residual drift of the spin-adiabatic weighted-MCH
diagonal gradient; leave it off unless a trajectory shows systematic energy
drift.

## Adaptive Timestep

### `dt_adaptive`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | nuclear propagation |

Shrink the timestep automatically when atoms move fast or the surface is stiff,
down to `dt_min`.

### `dt_min`

| Field | Value |
| --- | --- |
| Type | float (fs) |
| Default | `0.05` |
| Used by | adaptive timestep |

Minimum timestep for the adaptive scheme. Only used when `dt_adaptive=True`.

### `dx_max`

| Field | Value |
| --- | --- |
| Type | float (bohr) |
| Default | `0.02` |
| Used by | adaptive timestep |

Maximum per-step atomic displacement used as the adaptive-timestep criterion.
Only used when `dt_adaptive=True`.

## Notes

- **NAMD requires MRSF-TDDFT.** The input checker accepts `runtype=namd` only
  with [`[input] method=tdhf`](input.md#method) and
  [`[tdhf] type=mrsf`](tdhf.md#type).
- **Active-state range.** For plain FSSH, `1 <= active <= [tdhf] nstate`. For
  SOC-NAMD (`soc=true`), the spin-adiabatic manifold has `ns + 3*nt` states, so
  `1 <= active <= ns + 3*nt` (with `ns = nt = [tdhf] nstate`).
- **SOC-only keywords.** [`soc_basis`](#soc_basis), [`soc_du_dt_corr`](#soc_du_dt_corr),
  [`soc_tdc_grad_corr`](#soc_tdc_grad_corr), [`grad_wthr`](#grad_wthr),
  [`init_state`](#init_state), and [`econs`](#econs) apply only when
  `soc=true`. `init_state` overrides `active`; `econs` is a temporary
  stabilizer.
- **QM/MM dynamics.** Combine `[md]` with [`[qmmm]`](qmmm.md) and
  `[input] qmmm_flag=true` for embedded dynamics; see the
  [SOC-NAMD-QMMM workflow](../workflows/soc-namd-qmmm.md).

## Python API

In the compact `OpenQP` Python API,
[`job.workflow.namd(...)`](../python-scripting.md#qmmm-and-nonadiabatic-dynamics)
selects the surface-hopping run (`runtype=namd`) and sets `[md]` keywords; it
requires an MRSF-TDDFT theory. Pass `soc=True` (with an optional `soc_basis`) for
SOC-NAMD, and combine with [`job.qmmm(...)`](qmmm.md#python-api) for QM/MM
dynamics.

```python
from oqp.openqp import OpenQP

job = OpenQP("gas_socnamd", silent=1)
job.molecule(geometry="water", charge=0)
job.theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=3)

# SOC-NAMD (intersystem crossing); drop soc=... for internal-conversion FSSH
job.workflow.namd(soc=True, soc_basis="mch", nstep=200, dt=0.5,
                  init_state="S1", thrshe=0.1, init_temp=300.0)

mol = job.run()
```
