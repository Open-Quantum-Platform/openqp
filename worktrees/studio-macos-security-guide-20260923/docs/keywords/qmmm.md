# `[qmmm]`

The `[qmmm]` section configures hybrid quantum-mechanics/molecular-mechanics
(QM/MM) calculations. A quantum region described by an OpenQP method is
embedded in a classical force-field environment handled by
[OpenMM](https://openmm.org). The established `.inp` surface includes
single-point QM/MM energy and ground-state QM/MM molecular dynamics, while
NAMD uses the separate embedded surface-hopping driver.

In `.oqp`, `qmmm(...)` is accepted with `energy`, ground-state `md`,
and embedded `namd`; writing it also enables `qmmm_flag`. QM/MM gradients and
geometry optimizations are rejected because those generic backends do not yet
provide the assembled QM/MM gradient.

```text
dft/pbe0/def2-svp qmmm()
geom="ala.pdb 9 10 17 18 19"
```

!!! note "Available in OpenQP 1.3.0"
    The QM/MM implementation and its covalent-boundary `frontier_scheme`
    controls are included in OpenQP 1.3.0.

## Background

The QM subsystem is polarized by the MM electrostatic potential through the
**electrostatic potential-fitted (ESPF) operator**: the MM point charges enter
the QM core Hamiltonian, and the reaction of the QM density on the MM atoms is
represented by ESPF-fitted atomic charges. This gives an analytic,
energy-conserving QM/MM gradient. See
[References](../references.md#qmmm-espf-embedding) for the ESPF operator and its
periodic (particle-mesh Ewald) extension.

Two ways of defining the QM region are supported, matching the driver paths:

- **Single-point QM/MM energy** reads the geometry and QM selection from
  `[input] system = file.pdb <indices>` (see [`[input] system`](input.md#system)).
  Dangling covalent bonds crossing the QM/MM boundary are capped automatically
  (see [Link atoms](#link-atoms)).
- **Ground-state QM/MM molecular dynamics and NAMD-QMMM** read the PDB, force
  field, and QM selection from the `[qmmm]` keys `pdb_file`,
  `forcefield_files`, and `qm_atoms` below. NAMD additionally uses
  `runtype=namd` and the `[md]` section.

## Minimal QM/MM Example

Single-point QM/MM energy in `.oqp`:

```text
dft/bhhlyp/6-31g* qmmm()
geom="ala.pdb 9 10 17 18 19"
```

Python:

```python
from oqp.openqp import OpenQP

job = OpenQP("ala_qmmm", silent=1)
job.molecule("ala.pdb 9 10 17 18 19", basis="6-31g*")
job.theory("hf", functional="bhhlyp")
job.qmmm(embedding="electrostatic")
mol = job.run()
```

Legacy `.inp` (QM selection inline in `[input] system`):

```ini
[input]
qmmm_flag  = true
runtype    = energy
method     = hf
functional = bhhlyp
basis      = 6-31g*
system     = ala.pdb 9 10 17 18 19

[scf]
type = rhf
```

Embedded nonadiabatic molecular dynamics in `.oqp`:

```text
mrsf(nstate=5)/bhhlyp/6-31g namd
qmmm(forcefield_files="amber14-all.xml,amber14/tip3p.xml",qm_atoms="0-2",cutoff=PME)
geom="water_box.pdb 0-2"
```

The legacy `.inp` spelling is:

```ini
[input]
qmmm_flag  = true
runtype    = namd
method     = tdhf
functional = bhhlyp
basis      = 6-31g

[scf]
type = rohf
multiplicity = 3

[tdhf]
type   = mrsf
nstate = 5

[qmmm]
pdb_file         = water_box.pdb
forcefield_files = amber14-all.xml,amber14/tip3p.xml
qm_atoms         = 0-2
cutoff           = PME
embedding        = electrostatic
```

## Python API

In the compact `OpenQP` Python API, `job.qmmm(...)` enables QM/MM: it sets
[`[input] qmmm_flag=true`](input.md#qmmm_flag) and the `[qmmm]` section in one
call. `forcefield` is an alias for [`forcefield_files`](#forcefield_files); a
list is joined into the comma-separated string OpenQP expects, and `qm_atoms`
accepts a string (`"0-2"`) or a list of indices. Any other `[qmmm]` keyword can
be passed through as a keyword argument.

```python
from oqp.openqp import OpenQP

# Single-point QM/MM energy (QM selection inline in job.molecule)
job = OpenQP("ala_qmmm", silent=1)
job.molecule("ala.pdb 9 10 17 18 19", basis="6-31g*", charge=0)
job.theory("hf", functional="bhhlyp")
job.qmmm(embedding="electrostatic")
mol = job.run()
```

```python
# QM/MM molecular dynamics: PDB, force field, and QM atoms in job.qmmm(...)
job = OpenQP("water_box_qmmm", silent=1)
job.molecule("water_box.pdb 0 1 2", basis="6-31g")
job.theory.mrsf(functional="bhhlyp", nstate=5)
job.qmmm(
    pdb_file="water_box.pdb",
    forcefield=["amber14-all.xml", "amber14/tip3p.xml"],
    qm_atoms="0-2",
    cutoff="PME",
    embedding="electrostatic",
)
job.workflow.namd(nstep=200, dt=0.5)   # add soc=True for SOC-NAMD-QMMM
mol = job.run()
```

See the [SOC-NAMD-QMMM workflow](../workflows/soc-namd-qmmm.md) and
[Run from Python](../python-scripting.md#qmmm-and-nonadiabatic-dynamics) for the
full nonadiabatic QM/MM setup.

## Keywords

### `pdb_file`

| Field | Value |
| --- | --- |
| Type | string (path) |
| Default | *(empty)* |
| Used by | QM/MM molecular dynamics and SOC-NAMD-QMMM |

Path to the PDB file that defines the full QM+MM system (coordinates and
topology) for ground-state `runtype=md` and nonadiabatic `runtype=namd`.
Single-point QM/MM energy instead takes its PDB path and QM indices together
from `[input] system = file.pdb <indices>`.

### `forcefield_files`

| Field | Value |
| --- | --- |
| Type | string (comma- or space-separated list) |
| Default | *(empty)* |
| Used by | QM/MM molecular dynamics and SOC-NAMD-QMMM |

OpenMM force-field XML files applied to the MM region, e.g.
`amber14-all.xml,amber14/tip3p.xml` for a protein/water system or `tip3p.xml`
for a pure water box. Multiple files are combined by OpenMM in order.

### `forcefield`

| Field | Value |
| --- | --- |
| Type | string list |
| Default | `amber14-all.xml,amber14/tip3p.xml` |
| Used by | ground-state QM/MM setup |

Default force-field list used when a driver builds the MM system without an
explicit `forcefield_files` value. New QM/MM-MD decks should set
`forcefield_files` explicitly.

### `qm_atoms`

| Field | Value |
| --- | --- |
| Type | string (index list) |
| Default | *(empty)* |
| Used by | QM/MM molecular dynamics and SOC-NAMD-QMMM |

Zero-based indices of the atoms placed in the QM region, as individual indices
and/or ranges, e.g. `0 1 2` or `0-2` or `0-8 12 15`. Give the indices in
**ascending order**. This key is required by ground-state QM/MM MD and NAMD;
single-point energy writes the equivalent selection after its PDB path in
`[input] system`. Whole-molecule QM selections (e.g. a solute in a solvent box)
are the common case, and the only case supported by the nonadiabatic
(`runtype=namd`) path. In single-point energy and ground-state QM/MM MD, a
selection that cuts a covalent bond is capped with a hydrogen [link
atom](#link-atoms) and the MM frontier charge is treated per
[`frontier_scheme`](#frontier_scheme); see the [SOC-NAMD-QMMM scope
note](../workflows/soc-namd-qmmm.md#scope-and-limitations).

### `cutoff`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `NoCutoff` |
| Values | `NoCutoff`, `PME`, `Ewald`, `CutoffNonPeriodic`, `CutoffPeriodic` |
| Used by | QM/MM molecular dynamics and SOC-NAMD-QMMM |

OpenMM nonbonded method for the MM region. `NoCutoff` is used for isolated
(non-periodic) clusters. `PME` (particle-mesh Ewald) or `Ewald` select a
periodic box and enable ESPF-PME electrostatics for a solvated/periodic system;
these turn on the periodic branch of the driver.

### `embedding`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `electrostatic` |
| Values | `electrostatic`, `mechanical` |
| Used by | QM/MM electrostatic coupling |

Selects how the MM environment couples to the QM subsystem.

| Value | Meaning |
| --- | --- |
| `electrostatic` | Full ESPF electrostatic embedding. The MM charges polarize the QM density through the ESPF operator, and the QM density reacts on the MM atoms via ESPF-fitted charges. This is the production value and gives the analytic, energy-conserving QM/MM gradient. |
| `mechanical` | No electrostatic coupling into the QM Hamiltonian; the QM/MM interaction is mechanical (bonded/van der Waals) only. |

Legacy spellings such as `espf` and `split` appear in older decks; new inputs
should use `electrostatic`.

### `frontier_scheme`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `none` |
| Values | `none`, `rcd`, `rc`, `z1` |
| Used by | ESPF electrostatics at a covalent QM/MM boundary (ground-state QM/MM MD) |

When the QM/MM partition cuts a covalent bond, the MM host atom (`M1`, the MM end
of the severed bond) sits ~1.5 Å from the QM density. `frontier_scheme` selects
how that frontier charge is treated in the ESPF embedding. It is a **no-op for
whole-molecule QM regions** (no cut bond). Covalent QM/MM boundaries are handled
by the single-point and ground-state QM/MM MD paths; the nonadiabatic
(`runtype=namd`) path builds its QM molecule from `qm_atoms` only and does not
support a covalent cut.

| Value | Meaning |
| --- | --- |
| `none` | Full-field embedding: the QM density sees the complete MM charge set. This is the **default and the validated ESPF baseline**. ESPF couples the MM potential to the QM *atomic-charge operators* (`h += Σ_A φ_A Q̂_A`), not to the raw density via `1/\|r−R_M\|` integrals, which structurally suppresses the electron spill-out that motivates frontier redistribution in density-based embedding — so the ESPF method uses the full MM charges even at a covalent boundary. |
| `rcd` | Redistributed charge and dipole: delete `M1`'s charge and place virtual point charges at the `M1–M2` bond midpoints (`2·q₁/N`) plus `−q₁/N` on each MM neighbour `M2`, conserving the **total charge and the dipole about `M1`**. The virtual charges sit at bond midpoints (linear in the real atom positions), so the analytic gradient stays exact. |
| `rc` | Redistributed charge: midpoint charges only (conserves the total charge, not the dipole). |
| `z1` | Charge deletion: remove `M1`'s charge (conserves neither; provided for comparison). |

`rcd`/`rc`/`z1` are **optional refinements**, not the ESPF default. See
[References](../references.md#qmmm-espf-embedding) for the ESPF charge-operator
formulation and the redistribution schemes.

### `rigidwater`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | MM constraint setup |

Apply rigid-water (SHAKE/RATTLE) constraints to MM water molecules. QM atoms are
never constrained. Rigid water removes the stiff O-H stretch from the MM region
and allows a normal MD timestep (~0.5-1 fs). The nonadiabatic velocity-Verlet
loop always constrains MM rigid water; QM atoms move under the QM forces.

### `nonbondedmethod`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `NoCutoff` |
| Used by | ground-state QM/MM setup |

OpenMM nonbonded method for the ground-state QM/MM path. The nonadiabatic and
newer MD paths use `cutoff` instead.

### `constraints`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `None` |
| Used by | ground-state QM/MM MD |

OpenMM bond/angle constraint setting for the MM region in the ground-state MD
path (for example, constrain X-H bonds).

### `temperature`

| Field | Value |
| --- | --- |
| Type | float (K) |
| Default | `300.0` |
| Used by | ground-state QM/MM MD |

Target/initial temperature for the ground-state QM/MM MD path. The nonadiabatic
path sets the initial temperature from [`[md] init_temp`](md.md#init_temp).

### `timestep`

| Field | Value |
| --- | --- |
| Type | number (fs) |
| Default | `1` |
| Used by | ground-state QM/MM MD |

MD timestep for the ground-state QM/MM MD path. The nonadiabatic path uses
[`[md] dt`](md.md#dt).

### `nsteps`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | legacy static QM/MM bookkeeping |

Retained for compatibility with the older static-driver configuration. The
current ground-state OpenMM MD engine uses [`n_steps`](#n_steps); the
nonadiabatic path uses [`[md] nstep`](md.md#nstep). In concise `.oqp`, however,
`qmmm(nsteps=N)` is accepted as an alias and lowered to the active `n_steps`
key; this does not change the meaning of `nsteps` in a sectioned `.inp`.

### `istate`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` |
| Used by | ground-state QM/MM |

Obsolete numeric selector from the disconnected legacy `libopenmm` path. It is
retained in the sectioned schema for compatibility but is reserved in `.oqp`;
do not write `qmmm(istate=...)`. Ground-state canonical `md` uses `S0`, and
state-aware drivers own their physical state labels.

## Ground-state OpenMM MD Keywords

These keys are consumed by the command-line ground-state QM/MM-MD driver. The
command must run without MPI.

### `n_steps`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1000` |
| Used by | number of ground-state QM/MM-MD integration steps |

This is the preferred MD spelling. Concise `.oqp` accepts `nsteps` as an alias
for this key, while a traditional sectioned `.inp` keeps `nsteps` as separate
legacy bookkeeping with its own default of one.

### `ensemble`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `nve` |
| Values | `nve`, `nvt`, `npt` |

### `friction`

| Field | Value |
| --- | --- |
| Type | float, ps^-1 |
| Default | `1.0` |
| Used by | Langevin integration for NVT/NPT |

### `pressure`

| Field | Value |
| --- | --- |
| Type | float, bar |
| Default | `1.0` |
| Used by | NPT barostat target |

### `barostat_interval`

| Field | Value |
| --- | --- |
| Type | integer steps |
| Default | `25` |
| Used by | NPT barostat update interval |

NPT requires a periodic nonbonded method rather than `NoCutoff`.

### `trajectory_format`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `pdb` |
| Values | `pdb`, `dcd` |

### `trajectory_file`

| Field | Value |
| --- | --- |
| Type | path |
| Schema default | empty |
| Used by | coordinate trajectory output |

When omitted from the active MD configuration, the driver uses
`qmmm_trajectory.<trajectory_format>`.

### `log_file`

| Field | Value |
| --- | --- |
| Type | path |
| Schema default | empty |
| Used by | OpenMM state-data output |

The active MD driver uses `qmmm_trajectory.dat` when no path is supplied.

### `report_interval`

| Field | Value |
| --- | --- |
| Type | integer steps |
| Default | `1` |
| Used by | trajectory, log, and energy reporting frequency |

### `energy_file`

| Field | Value |
| --- | --- |
| Type | path |
| Schema default | empty |
| Used by | saved QM, MM, kinetic, and total-energy arrays |

The active MD driver uses `total_energy.npz` when no path is supplied.

### `qm_atoms_xyz`

| Field | Value |
| --- | --- |
| Type | path |
| Default | empty |
| Used by | optional replacement coordinates for the QM atoms |

### `qm_list`

| Field | Value |
| --- | --- |
| Type | integer list |
| Default | empty |
| Used by | mapping QM atoms to rows of `qm_atoms_xyz` |

When `qm_atoms_xyz` is present and `qm_list` is omitted, rows are used in QM
atom order. If supplied, `qm_list` must have the same length as the QM atom
selection and contain valid XYZ indices.

## Covalent QM/MM boundaries

Most QM/MM calculations put a whole molecule (or a set of whole molecules) in the
QM region — a chromophore in a solvent box, a substrate in a binding pocket. When
the QM/MM partition instead **cuts through a covalent bond** (e.g. carving a side
chain out of a protein backbone), two problems appear at the frontier, and OpenQP
handles each automatically.

### 1. The dangling bond → hydrogen link atom

Cutting a bond leaves the QM frontier atom with an unsatisfied valence. OpenQP
caps it with a **hydrogen link atom** on the broken QM–MM bond at the scaled
(IMOMM) position `R_L = R_QM + g·(R_MM − R_QM)`, with `g = (r_H + r_QM)/(r_QM +
r_MM)` from covalent radii. Because the cap position is a fixed linear function of
its two real host atoms, the force computed on the link atom is redistributed
onto them by the chain rule — no extra degrees of freedom, and the analytic
gradient stays exact. Details in [Link atoms](#link-atoms).

### 2. The frontier MM charge → ESPF is robust; `frontier_scheme` refines it

The MM host atom (`M1`, the MM end of the severed bond) sits only ~1.5 Å from the
QM density. In a **density-based** electrostatic embedding, embedding `M1`'s full
point charge that close over-polarizes the QM density (electron spill-out), which
is why such methods delete or shift the frontier charge.

OpenQP's embedding is **ESPF**, which is structurally more robust here: the MM
potential couples to the QM *atomic-charge operators*, `h += Σ_A φ_A Q̂_A`
(Huix-Rotllant & Ferré; see [Background](#background)), **not** to the raw density
via `1/|r−R_M|` integrals. Because the coupling runs through bounded atomic-charge
operators evaluated at the QM nuclei, the spill-out that plagues density-based
embedding is suppressed — so the full MM charges can be used "with no scaling"
even at a covalent boundary. That is the default,
[`frontier_scheme = none`](#frontier_scheme).

For a still-more-faithful boundary field, the optional redistributions delete the
raw `M1` charge and spread it over **virtual point charges** at the `M1–M2` bond
midpoints (`M2` = the MM neighbours of `M1`):

- **`rcd`** — `2·q₁/N` at each of the `N` midpoints plus `−q₁/N` on each `M2`,
  conserving both the total charge **and** the dipole about `M1`.
- **`rc`** — `q₁/N` at each midpoint (total charge only).
- **`z1`** — delete `M1` (conserves neither; for comparison).

Each virtual charge sits at a bond midpoint — a linear function of the real
`M1`/`M2` positions — so the electrostatic force on it redistributes onto its
hosts by the chain rule, exactly like the link-atom cap, keeping the analytic
gradient exact.

### Gradient accuracy

The full analytic QM/MM gradient (link-atom projection + ESPF density response +
frontier redistribution) matches a finite-difference of the self-consistent QM/MM
energy to ~1e-3 au — the ESPF grid-derivative floor, a property of the embedding
independent of the boundary treatment. The redistribution does not degrade this;
at the deleted host `M1` the redistributed schemes are in fact *more* accurate
than full-field, because they remove the raw close-range point charge.

### Which path supports it

Covalent-boundary QM/MM is available in the **single-point** and **ground-state
QM/MM MD** (`runtype=md`) paths. The nonadiabatic (`runtype=namd`,
[SOC-NAMD-QMMM](../workflows/soc-namd-qmmm.md)) path builds its QM molecule from
`qm_atoms` only and **raises on a covalent cut** — use the ground-state MD path.
A runnable deck is `examples/QMMM/ala-dipeptide_BHHLYP-QMMM-MD-RCD.inp` (alanine
dipeptide, QM = the C-terminal amide cutting the `ALA C–CA` bond,
`frontier_scheme=rcd`); see [QM/MM examples](../examples/index.md#qmmm-examples).

## Link atoms

When the QM region defined through `[input] system = file.pdb <indices>` cuts a
covalent bond, OpenQP caps each dangling bond automatically with a **hydrogen
link atom** — there is no keyword to enable it; the boundary bonds are detected
from the PDB topology. The link hydrogen is placed along the broken QM-MM bond at
a scaled distance set by the covalent-radius factor

```
g = (r_H + r_QM) / (r_QM + r_MM)
```

where `r_H`, `r_QM`, and `r_MM` are the covalent radii of hydrogen and of the QM
and MM boundary atoms. Only hydrogen capping is currently supported. The
link-atom energy gradient is redistributed onto its two real host atoms by the
chain rule of this scaled position, so no extra degrees of freedom are added.

Across a covalent boundary the MM frontier-host charge is treated per
[`frontier_scheme`](#frontier_scheme) (default `none` = full-field ESPF), and
the ESPF grid switching width is selected automatically — see
[ESPF grid switching](#espf-grid-switching) below.

!!! note "Covalent boundaries are not supported in nonadiabatic dynamics"
    Automatic link-atom capping applies to the single-point and ground-state
    QM/MM MD paths. The nonadiabatic (`runtype=namd`,
    [SOC-NAMD-QMMM](../workflows/soc-namd-qmmm.md)) path builds the QM molecule
    from `qm_atoms` only, so it supports **whole-molecule** QM regions and raises
    on a covalent cut; use the ground-state QM/MM MD path for covalent-boundary
    QM/MM.

## ESPF grid switching

The ESPF embedding gradient carries a small finite-difference residual (the
grid-derivative floor, ~1e-3 a.u.). At a **covalent** boundary that residual
does not spread over the QM region: it concentrates as an over-shoot on the MM
host atom, and it is governed by the smooth-switching width `ESPF_SWSCALE`.

The shipped whole-molecule value, `1.8`, over-smooths at a covalent boundary,
and the tighter `1.5` is about 20% worse without one — neither is a safe
global default. OpenQP therefore selects it **per system**, from whether the
QM/MM partition actually cut a bond:

| QM/MM boundary | `ESPF_SWSCALE` |
| --- | --- |
| covalent (a bond is cut; link atoms present) | **1.5**, selected automatically |
| whole-molecule (no bond cut) | `1.8` |

The choice is made by the QM/MM driver when it builds the partition, because
that is where it is known whether a bond was cut: the native ESPF gradient is
handed a grid and a density and cannot tell. A covalent-boundary run reports
the selection:

```
[QM/MM] 1 covalent boundary bond(s); frontier-charge embedding = ...
[QM/MM] covalent boundary detected; ESPF_SWSCALE=1.5 (whole-molecule default
        is 1.8). Set ESPF_SWSCALE to override.
```

Setting `ESPF_SWSCALE` yourself always wins and the automatic choice stands
down entirely, so a parameter sweep or a reproduction is never silently
overridden:

```bash
export ESPF_SWSCALE=1.8      # force the whole-molecule width at a boundary
```

### ESPF environment variables

These are runtime **environment variables**, not `[qmmm]` keywords, read by the
native ESPF module:

| Variable | Default | Meaning |
| --- | --- | --- |
| `ESPF_SMOOTH` | on | enable smooth grid switching at all |
| `ESPF_SWSCALE` | selected (1.5 / 1.8) | switching width scale, see above |
| `ESPF_SWDELTA` | 0.7 | switching offset |
| `ESPF_WSCALE` | 1.0 | grid weight scale; keep at 1.0 |
| `ESPF_WDERIV` | on | include the grid weight-derivative gradient term |

!!! note "Applies to every electrostatically embedded QM/MM path"
    Grid switching is a property of baseline ESPF, so it affects every path
    that runs [`embedding = electrostatic`](#embedding), including
    [SOC-NAMD-QMMM](../workflows/soc-namd-qmmm.md), and it is independent of
    [`frontier_scheme`](#frontier_scheme). Under `embedding = mechanical`
    there is no electrostatic coupling into the QM Hamiltonian and therefore
    no ESPF grid at all, so none of these variables — including the automatic
    `ESPF_SWSCALE` selection — has any effect.

## Notes

- Set `[input] qmmm_flag=true` to activate any QM/MM path; without it the
  `[qmmm]` section is ignored.
- For nonadiabatic QM/MM dynamics, combine this section with the
  [`[md]`](md.md) section and see the
  [SOC-NAMD-QMMM workflow](../workflows/soc-namd-qmmm.md).
- Use `cutoff = PME` (or `Ewald`) with a solvated periodic water box for
  production QM/MM-MD; `NoCutoff` is for isolated clusters.
