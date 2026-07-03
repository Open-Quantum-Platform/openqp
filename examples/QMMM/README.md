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

Minimal (`[md] nstep=1`) nonadiabatic-dynamics demonstrations on formaldehyde
(QM) solvated by 5 TIP3P waters (MM), `NoCutoff` (non-periodic cluster):

| Input | What it shows |
| --- | --- |
| `H2CO-water_BHHLYP-MRSF-NAMD-QMMM.inp` | MRSF-TDDFT FSSH (internal conversion, `[md] soc=false`) with ESPF QM/MM. |
| `H2CO-water_BHHLYP-SOC-NAMD-QMMM.inp` | SOC-NAMD (intersystem crossing, `[md] soc=true`) on the spin-adiabatic manifold with ESPF QM/MM. |

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
