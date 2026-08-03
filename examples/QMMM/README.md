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

## NAMD-QMMM (surface-hopping dynamics)

Minimal nonadiabatic-dynamics demonstrations on formaldehyde (QM) solvated by
5 TIP3P waters (MM), `NoCutoff` (non-periodic cluster):

| Input | What it shows |
| --- | --- |
| `H2CO-water_BHHLYP-MRSF-NAMD-QMMM.inp` | MRSF-TDDFT FSSH (internal conversion, `[md] soc=false`) with ESPF QM/MM. |
| `H2CO-water_BHHLYP-MRSF-NAMD-QMMM.oqp` | Two-step semantic-input version that writes a restart checkpoint. |
| `H2CO-water_BHHLYP-MRSF-NAMD-QMMM.restart.oqp` | Paired continuation that loads the step-2 checkpoint and advances through step 3. `openqp --run_tests all` schedules it after the producer and reuses the same isolated run directory. |
| `H2CO-water_BHHLYP-SOC-NAMD-QMMM.inp` | SOC-NAMD (intersystem crossing, `[md] soc=true`) on the spin-adiabatic manifold with ESPF QM/MM. |

To exercise checkpoint loading, select the semantic examples in the regression
runner. It gives the producer and continuation the same project directory and
runs the continuation only after the checkpoint-producing job finishes:

```bash
openqp --run_tests examples/QMMM --input-format oqp
```

Auxiliary files: `formaldehyde_water.pdb` (QM+MM coordinates/topology),
`formaldehyde.xml` (minimal QM-residue force field — only the Lennard-Jones
parameters matter; QM electrostatics come from ESPF), and `tip3p.xml` (water).

NAMD writes a trajectory log (`<project>.log`), not a regression `.json`, so
these serve as runnable demonstrations rather than numeric regression tests. See
the [SOC-NAMD-QMMM workflow](https://open-quantum-platform.github.io/openqp-docs/workflows/soc-namd-qmmm/)
and the `[md]` / `[qmmm]` keyword pages in the manual for the full input
contract and the compact `job.qmmm(...)` / `job.workflow.namd(...)` Python API.

## Single-point / ground-state QM/MM

`ala.inp` and `2E4E_RHF-DFT-QMMM_energy.inp` are QM/MM single-point energy decks
(QM selection via `[input] system = file.pdb <indices>`); `run.inp` is a
ground-state OpenMM-integrator QM/MM MD deck.

## Covalent QM/MM boundary — `[qmmm] frontier_scheme`

When the QM/MM partition cuts a covalent bond, the dangling QM bond is capped
with a hydrogen link atom and the MM host atom (`M1`) sits ~1.5 Å from the QM
density. `[qmmm] frontier_scheme` selects how that frontier charge is treated in
the ESPF electrostatics. Covalent QM/MM boundaries are handled by the
ground-state QM/MM MD path (`QMMM_MD`); the nonadiabatic `runtype=namd` path does
not yet append link atoms to its QM molecule and raises on a covalent cut.

| value | meaning |
| --- | --- |
| `none` (default) | Full-field: the QM density sees the complete MM charge set. This is the **validated ESPF baseline** — ESPF couples the MM potential to QM *atomic-charge operators* (`h += Σ_A φ_A Q̂_A`, Huix-Rotllant & Ferré, *JCTC* 2021, 17, 538, eq 6), which already suppresses the electron spill-out that motivates redistribution in density-based embedding, so the ESPF papers use full MM charges even at a covalent protein boundary. |
| `rcd` | Delete `M1`'s charge and redistribute it to virtual point charges at the `M1–M2` bond midpoints, conserving the **total charge and the dipole about `M1`**. Gradient-consistent (the midpoints are linear in the real atom positions). |
| `rc` | As `rcd` but conserving only the total charge. |
| `z1` | Delete `M1`'s charge (conserves neither; for comparison). |

`rcd`/`rc`/`z1` are **optional refinements**, not the ESPF default. Enable via the
input (`[qmmm] frontier_scheme = rcd`) or the Python API
(`job.qmmm(..., frontier_scheme="rcd")`). It is a no-op for whole-molecule QM
regions (no cut bond).

A runnable covalent-boundary deck is
`ala-dipeptide_BHHLYP-QMMM-MD-RCD.inp` — the alanine dipeptide (ACE-ALA-NH2) with
AMBER-14, QM = the C-terminal amide so the QM/MM partition cuts the `ALA C–CA`
backbone bond, run as ground-state QM/MM MD (`runtype=md`) with
`frontier_scheme=rcd`:

```bash
cd examples/QMMM && openqp ala-dipeptide_BHHLYP-QMMM-MD-RCD.inp
```

Like the other ground-state QM/MM decks it is skipped by `openqp --run_tests all`
(covalent-boundary QM/MM is the ground-state MD path, not `runtype=namd`). The
same alanine boundary is exercised automatically — link-atom detection +
frontier-charge conservation on the real AMBER-14 charges — in
`tests/test_qmmm_frontier_openmm.py` (OpenMM-gated), and the pure redistribution
math (including a finite-difference gradient check) in
`tests/test_qmmm_frontier.py`.
