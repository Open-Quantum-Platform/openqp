# QM/MM examples

ESPF electrostatic QM/MM: a quantum region (HF/DFT/MRSF-TDDFT) embedded in a
classical (OpenMM) MM environment. These examples require the optional **OpenMM**
backend (`pip install openmm`) and read auxiliary topology/force-field files
(`*.pdb`, `*.xml`) from this directory.

The NAMD-QMMM examples (`runtype=namd`) resolve their auxiliary files relative
to the input file, so they **are part of `openqp --run_tests all`**. When OpenMM
is not installed they are reported **SKIPPED** (like the ddX/PCM examples on a
build without ddX), so the suite stays green either way. Run one directly:

```bash
openqp examples/QMMM/H2CO-water_BHHLYP-MRSF-NAMD-QMMM.inp
```

## Defining the QM region

Name the QM atoms **once**, with `qmmm(qm_atoms=...)`: 0-based indices in PDB
order, as ranges or lists (`"0-3"`, `"8,9,16,17,18"`, `"0-10,21,22"`). With
`geom="file.pdb"` the QM molecule is built from that list, including one
hydrogen link atom per covalent bond the selection cuts.

```
qmmm(forcefield_files="amber14-all.xml",qm_atoms="8,9,16,17,18",cutoff=NoCutoff)
geom="ala.pdb"
```

- `charge=` is the charge of the QM region, not of the whole system.
- The older spelling that repeats the atoms 1-based after the PDB name
  (`geom="ala.pdb 9 10 17 18 19"`, or `[input] system = ala.pdb 9 10 17 18 19`)
  is still accepted, alone or together with `qm_atoms`. Given together, the two
  lists must name the same atoms; a mismatch is an error rather than a QM
  molecule that differs from the embedded region
  (`ala-dipeptide_RHF-QMMM-OPT-constraints.oqp` keeps both).
- A QM region of whole molecules may instead take its atoms from an xyz file
  (`geom="solute.xyz"`), as the formaldehyde NAMD decks do. The coordinates
  are then taken from the PDB, so the xyz only supplies the elements, in PDB
  order. It cannot add link atoms.

## The four dynamics combinations

`md(...)` owns the nuclear propagation, `namd(...)` the electronic/hopping side
and `qmmm(...)` the embedding, so the same controls read the same way whether
the dynamics is ground state or nonadiabatic, gas phase or embedded. One
runnable two-step `.oqp` deck per combination:

| | gas phase | QM/MM |
| --- | --- | --- |
| **ground state** (`md`) | [`examples/MD/Mg_RHF-GROUND-STATE-MD.oqp`](../MD/Mg_RHF-GROUND-STATE-MD.oqp) | [`ala-dipeptide_RKS-QMMM-GROUND-STATE-MD.oqp`](ala-dipeptide_RKS-QMMM-GROUND-STATE-MD.oqp) |
| **surface hopping** (`namd` + `md`) | [`examples/NAMD/C2H4_BHHLYP-MRSF-NAMD-SCHEME.oqp`](../NAMD/C2H4_BHHLYP-MRSF-NAMD-SCHEME.oqp) | [`H2CO-water_BHHLYP-MRSF-NAMD-QMMM-SCHEME.oqp`](H2CO-water_BHHLYP-MRSF-NAMD-QMMM-SCHEME.oqp) |

Adding `qmmm(...)` is the only difference between a column, and adding
`namd(...)` the only difference between a row:

```
rks/bhhlyp/6-31g md(nstep=2,dt=0.5,velocity=zero,ensemble=nve,trajectory_file="ala-gsmd.pdb")
qmmm(forcefield_files="amber14-all.xml",qm_atoms="8,9,16,17,18",cutoff=NoCutoff)
geom="ala.pdb 9 10 17 18 19"
```

```
mrsf(nstate=2)/bhhlyp/6-31g* namd(S1,scheme=Overlap,nacme_check=baeck_an,nacme_policy=warn)
md(nstep=2,dt=0.5,velocity=zero,seed=3,trajectory_interval=1,trajectory_file="H2CO-water-scheme.namd.trj")
qmmm(pdb_file="formaldehyde_water.pdb",forcefield_files="formaldehyde.xml tip3p.xml",qm_atoms="0-3",cutoff=NoCutoff)
geom="../geometries/CH2O-2bc62dda4b8a.xyz"
```

Check any of them with the regression runner (`Passed: 1` means the deck ran):

```bash
openqp --run_tests examples/QMMM/ala-dipeptide_RKS-QMMM-GROUND-STATE-MD.oqp
```

Both QM/MM decks need OpenMM; the two gas-phase ones do not. `cutoff=NoCutoff`
is the schema default and is spelled out here only because these systems are
non-periodic clusters — see `ala-box_BHHLYP-MRSF-NAMD-QMMM-PME.inp` for the
periodic counterpart.

## NAMD-QMMM (surface-hopping dynamics)

Minimal nonadiabatic-dynamics demonstrations on formaldehyde (QM) solvated by
5 TIP3P waters (MM), `NoCutoff` (non-periodic cluster):

| Input | What it shows |
| --- | --- |
| `H2CO-water_BHHLYP-MRSF-NAMD-QMMM.inp` | Two-step NVE MRSF-TDDFT FSSH (internal conversion, `[md] soc=false`) with ESPF QM/MM, an independent conservative water-droplet boundary, and a solute-COM restraint. |
| `H2CO-water_BHHLYP-MRSF-NAMD-QMMM.oqp` | Semantic-input version of the same two-step calculation that writes a restart checkpoint. |
| `H2CO-water_BHHLYP-MRSF-NAMD-QMMM.restart.oqp` | Paired continuation that loads the step-2 checkpoint and advances through step 3. `openqp --run_tests all` schedules it after the producer and reuses the same isolated run directory. |
| `H2CO-water_BHHLYP-MRSF-NAMD-QMMM-NVT.inp` | One-step NVT smoke run with the independent Langevin thermostat and separately recorded energy exchange. |
| `H2CO-water_BHHLYP-SOC-NAMD-QMMM.inp` | SOC-NAMD (intersystem crossing, `[md] soc=true`) on the spin-adiabatic manifold with ESPF QM/MM. |
| `ala-dipeptide_BHHLYP-MRSF-NAMD-QMMM-linkatom.inp` | Two-step MRSF-TDDFT FSSH across a **covalent QM/MM boundary** (hydrogen link atom): alanine dipeptide, QM = the C-terminal amide, `NoCutoff`. |
| `ala-dipeptide_RHF-QMMM-OPT-linkatom.inp` | **QM/MM geometry optimisation** across the same covalent boundary (RHF/6-31G): minimises the embedded QM/MM energy over the QM atoms with the MM fixed (`[optimize] qmmm_radius=0`), 12 steps, writes the full-system PDB. |
| `ala-dipeptide_RHF-QMMM-OPT-constraints.inp` | **QM/MM optimisation with held bonds**: as above, with the MM atoms within 3 Å of the QM region free to move (`[optimize] qmmm_radius=3.0`) and their X–H bond lengths held by `[qmmm] constraints=HBonds` (frozen distances in the native optimizer), 6 steps; writes the geometry to the file named by `[optimize] qmmm_output`. |
| `ala-box_BHHLYP-MRSF-NAMD-QMMM-PME.inp` | The same boundary in a **periodic TIP3P box** (`cutoff=PME`, Ewald QM/MM electrostatics with the self-consistent QM-image term) exercising `ewald_tol`, `lj_switch`, `h_lj` and `mm_charge_width`. |
| `ala-box_BHHLYP-MRSF-NAMD-QMMM-PME-FLEXIBLE-ODP.oqp` | The periodic box with water HOH 76 moved into the QM region, kept there by the FlexiBLE boundary (`flexible(...)`), sampled as one ODP umbrella window (`odp(...)`) in Langevin NVT: the setup for QM/MM free energies of solvent-involving reactions. |

To exercise checkpoint loading, select the semantic examples in the regression
runner. It gives the producer and continuation the same project directory and
runs the continuation only after the checkpoint-producing job finishes:

```bash
openqp --run_tests examples/QMMM --input-format oqp
```

Auxiliary files: `formaldehyde_water.pdb` (QM+MM coordinates/topology),
`formaldehyde.xml` (minimal QM-residue force field — only the Lennard-Jones
parameters matter; QM electrostatics come from ESPF), `tip3p.xml` (water),
`ala.pdb` (alanine dipeptide in vacuum) and `ala_box.pdb` (the dipeptide with
106 TIP3P waters in a 16 Å cubic box, AMBER-14 + `amber14/tip3p.xml`). A PDB
named in `[input] system = file.pdb <indices>` is looked up relative to the
working directory first and then next to the input file, like the `[qmmm]`
auxiliary files, so every NAMD deck runs from any directory.

### Solvent molecules in the QM region — `[flexible]`

QM waters diffuse away from the solute and MM waters diffuse in. The FlexiBLE
boundary potential of Shen and Glover (J. Chem. Phys. 155, 224112 (2021))
holds the separation: each solvent molecule is placed by its anchor atom
(`anchor=O`, or `com`) at its distance from a centre, and

    V = -kT ln[ h(QM, MM) / sum over partitions h ],
    h = exp(-sum_{i in QM, j in MM} g(r_i - r_j)),  g(R) = (aR)^3/(1 + aR) for R > 0

so the actual partition costs nothing when every QM molecule is closer than
every MM molecule. A whole residue named in `resnames` is a QM molecule when
all its atoms are in `[qmmm] qm_atoms` and an MM molecule otherwise; a residue
split between the regions is an error, and so is a massless anchor (a
virtual site cannot carry the boundary force). `[flexible]` needs QM/MM and is
refused in SOC-NAMD, whose propagators do not apply the conservative
restraints. For a solute in water:

```text
flexible(enabled=true, alpha=5.0, temperature=300, center=solute, atoms="0-11", anchor=O)
```

`alpha` is per angstrom (the plugin's tests use 5.0, i.e. 50 nm^-1); a small
alpha lets very many partitions contribute and the run stops with an error
rather than enumerating them. `temperature` is the kT of the potential and is
set explicitly. `center=solute` uses the centre of mass of `atoms` (0-based;
by default the QM atoms that are not solvent; in a periodic cell made whole
along the topology bonds), `center=point` the fixed
`point` in angstrom, and `center=molecules` the centre of mass of all
solvent molecules (non-periodic only). Periodic runs (`cutoff=PME`) use the
minimum-image distance and stop when a QM solvent molecule is not closer than
half the box (MM molecules farther out contribute nothing). The partition sum follows the authors' OpenMM reference kernel
(`threshold`, `scale`, `max_iterations`, `keep_below_threshold`), except
that a first sum of exactly one is re-evaluated at the lowered threshold
instead of being accepted (the kernel's shortcut drops every partition below
the initial threshold and can return zero boundary energy and force), and
`tests/test_flexible_boundary.py` checks energies and forces against values
computed with that plugin (`tests/data/flexible_boundary_plugin_reference.json`;
the generator is kept in openqp-devkit).

NAMD writes a trajectory log (`<project>.log`), not a regression `.json`, so
these serve as runnable demonstrations rather than numeric regression tests. See
the [SOC-NAMD-QMMM workflow](https://open-quantum-platform.github.io/openqp-docs/workflows/soc-namd-qmmm/)
and the `[md]` / `[qmmm]` keyword pages in the manual for the full input
contract and the compact `job.qmmm(...)` / `job.workflow.namd(...)` Python API.

## Ground-state QM/MM energy at a fixed geometry

`ala.oqp` (the alanine dipeptide) and `2E4E_RHF-DFT-QMMM_energy.oqp` (a 129-atom
peptide) take **one** `md(nstep=1)` step, which evaluates the embedded QM/MM
energy and forces at the input geometry; `run.oqp` is a two-step ground-state
QM/MM MD of a water dimer in a periodic box.

### Restart, and snapshots for surface hopping

Initial conditions for QM/MM surface hopping come from an equilibrated
ground-state trajectory. The usual protocol and the piece of input that
connects each stage:

| stage | level | ensemble | connects to the next stage by |
| --- | --- | --- | --- |
| 1 | classical (pure MM, run in OpenMM) | NPT, ~ns | a snapshot written with `oqp.utils.md_snapshot.write_snapshot`, read by `md(snapshot=...)` |
| 2 | QM/MM ground state | NVT, a few ps | `md(restart=true)` to continue; `md(snapshot_interval=N)` to save a point every N steps |
| 3 | QM/MM surface hopping | NVE | `md(snapshot="...")`, one trajectory per saved point |

A **snapshot** is one small `.npz` file holding the positions, the velocities
and the cell of every atom, QM and MM (`pyoqp/oqp/utils/md_snapshot.py`).

- **Checkpoint and restart.** Ground-state QM/MM MD rewrites
  `md(restart_file=...)` (default `qmmm_md.restart.npz`) every
  `restart_interval` steps and at the end of the run. Adding `restart=true` to
  the same deck continues it: same step counter, trajectory and energy log
  appended, rows and frames written after the checkpoint dropped. `nstep` is
  the **total** length, so raise it. The ensemble may change between stages
  (NVT, then NVE). An NVE run stopped and continued reproduces the
  uninterrupted run.
- **Numbered snapshots.** `md(snapshot_interval=N)` also writes
  `<name>.snapshot.<step>.npz` every N steps
  (`H2CO-water_RKS-QMMM-MD-SNAPSHOTS.oqp`). Space them 50-100 fs apart.
- **Starting from a snapshot.** `md(snapshot="file.npz")` starts a new
  trajectory from that point, for ground-state QM/MM MD and for QM/MM surface
  hopping (`H2CO-water_BHHLYP-MRSF-NAMD-QMMM-SNAPSHOT.oqp`). No velocities are
  drawn and none are rescaled to a temperature; without it QM/MM NAMD always
  starts from fresh Maxwell-Boltzmann velocities. It cannot be combined with
  `velocity=`.
- **Rigid water.** QM/MM surface hopping keeps MM water rigid, so equilibrate
  with `qmmm(rigidwater=true)`; a snapshot with flexible water is refused.
- **From a classical run.** In the OpenMM script of stage 1:

  ```python
  from oqp.utils.md_snapshot import write_snapshot
  state = simulation.context.getState(getPositions=True, getVelocities=True)
  box = state.getPeriodicBoxVectors()
  write_snapshot("classical.npz",
      positions_nm=state.getPositions(asNumpy=True).value_in_unit(nanometer),
      velocities_nm_ps=state.getVelocities(asNumpy=True).value_in_unit(nanometer/picosecond),
      masses_dalton=[system.getParticleMass(i).value_in_unit(dalton)
                     for i in range(system.getNumParticles())],
      box_nm=[box[i][i].value_in_unit(nanometer) for i in range(3)])
  ```

  The atom order must be that of the PDB given to `qmmm(pdb_file=...)`; a
  snapshot of a different system is refused. The cell of the snapshot replaces
  the `CRYST1` cell of the PDB.

### Constant pressure

There is no `ensemble=npt`, for any driver. A QM/MM barostat would need the
full QM/MM energy at every trial cell, and a QM/MM run covers a few
picoseconds, far too short for a density to equilibrate. Equilibrate the cell
in the classical stage and hand it over with `md(snapshot=...)`, as above; the
snapshot's cell replaces the `CRYST1` cell of the PDB.

There is deliberately no `runtype=energy` QM/MM deck: the plain single-point
route never builds the ESPF embedding records (`OQP::ESPF_CORR` / `OQP::POTQM`)
— only the force-based driver behind `md`, `optimize` and `namd` does — so it
aborts inside the SCF. One `md` step is the supported spelling.

`run.inp` and `ala-dipeptide_BHHLYP-QMMM-MD-RCD.inp` omit `[input] system`
on purpose: they are the fixtures for the legacy converter's
synthesize-`geom`-from-`[qmmm] pdb_file` path, and they run, since the QM
molecule is built from `qm_atoms` (see below).

## Covalent QM/MM boundary — `[qmmm] frontier_scheme`

When the QM/MM partition cuts a covalent bond, the dangling QM bond is capped
with a hydrogen link atom and the MM host atom (`M1`) sits ~1.5 Å from the QM
density. `[qmmm] frontier_scheme` selects how that frontier charge is treated in
the ESPF electrostatics. Covalent QM/MM boundaries are handled by both the
ground-state QM/MM MD path (`QMMM_MD`) and the nonadiabatic `runtype=namd`
paths (FSSH, SOC-NAMD). For NAMD the QM molecule must contain the link
hydrogens, so build it from the PDB (`geom="file.pdb"` with
`qmmm(qm_atoms=...)`), which appends one
H per cut bond in the order the driver detects them. The link atoms carry no
dynamical degrees of freedom: their positions follow the two host atoms, their
forces are chain-ruled onto the hosts, and the surface-hopping velocity
rescaling acts on the real QM atoms only.

| value | meaning |
| --- | --- |
| `rcd` (default) | Delete `M1`'s charge and redistribute it to virtual point charges at the `M1–M2` bond midpoints, conserving the **total charge and the dipole about `M1`** (Lin & Truhlar, *J. Phys. Chem. A* 2005, 109, 3991). Gradient-consistent (the midpoints are linear in the real atom positions). |
| `rc` | As `rcd` but conserving only the total charge. |
| `z1` | Delete `M1`'s charge (conserves neither; for comparison). |
| `none` | Keep the full `M1` charge (for comparison only; a warning is printed). |

`none` is not a usable model at a covalent cut. The link hydrogen sits on the
cut bond 0.4–0.5 Å from `M1`, so the full `M1` charge is a point charge inside
the link-atom density and the ESPF energy is unbounded below: optimization and
MD collapse the QM frontier atom onto `M1` (issue #551, bovine rhodopsin
Lys296 CB|CA). On the ACE CH3|C cut of the alanine dipeptide (`M1` charge
+0.60), RHF/6-31G with `none` falls from −41.5 to −249 Eh in ten optimization
steps, while `rcd` converges at −40.872 Eh. Select another scheme with
`[qmmm] frontier_scheme = rc` or `job.qmmm(..., frontier_scheme="rc")`. The
setting is a no-op for whole-molecule QM regions (no cut bond).

A runnable covalent-boundary deck is
`ala-dipeptide_BHHLYP-QMMM-MD-RCD.oqp` — the alanine dipeptide (ACE-ALA-NH2) with
AMBER-14, QM = the C-terminal amide so the QM/MM partition cuts the `ALA C–CA`
backbone bond, run as ground-state QM/MM MD (`md(...)`) with
`frontier_scheme=rcd`:

```bash
openqp --run_tests examples/QMMM/ala-dipeptide_BHHLYP-QMMM-MD-RCD.oqp
```

Like the other ground-state QM/MM decks it is skipped by `openqp --run_tests all`.
The nonadiabatic decks `ala-dipeptide_BHHLYP-MRSF-NAMD-QMMM-linkatom.inp`
(vacuum) and `ala-box_BHHLYP-MRSF-NAMD-QMMM-PME.inp` (periodic box) run the
same covalent boundary with MRSF-TDDFT surface hopping and are part of the
suite. The same alanine boundary is exercised automatically — link-atom detection +
frontier-charge conservation on the real AMBER-14 charges — in
`tests/test_qmmm_frontier_openmm.py` (OpenMM-gated), and the pure redistribution
math (including a finite-difference gradient check) in
`tests/test_qmmm_frontier.py`.

## ODP umbrella windows in ground-state QM/MM MD

`H2CO-water_HF-QMMM-MD-ODP.oqp` runs one ODP umbrella window in ground-state
QM/MM MD (`md(...) qmmm(...)`), with the `droplet` wall and the `solute_com`
restraint, as in the original ODP paper (Baek and Choi, JCTC 14, 2312 (2018)).
Each step of `total_energy.npz` holds the unbiased potential, every bias term
and the ODP coordinate; `oqp.library.odp.odp_wham` (or `examples/ODP/odp_wham.py`)
combines the tables of all windows. A restart (`md(restart=true)`) must keep
the bias, the restraints, the QM/MM Hamiltonian and the temperature of the
table it continues; periodic runs need an orthorhombic cell.

