# OpenQP

Open Quantum Platform ([OpenQP](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01117))
is a source-available quantum chemistry program centered on
[mixed-reference spin-flip TDDFT (MRSF-TDDFT)](https://doi.org/10.1021/acs.jpclett.3c02296).
It provides HF and DFT references, correlated wavefunction methods, linear-response
excited states, nuclear derivatives, reaction-path calculations, nonadiabatic
dynamics, and QM/MM calculations in one program.

MRSF-TDDFT retains the linear-response structure of TDDFT while removing the spin
contamination that limits conventional spin-flip TDDFT. This permits consistent
descriptions of multiconfigurational ground and excited states, including
diradicals, bond breaking, conical intersections, photochemical processes, and
spin-orbit coupling.

Use the **[OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/)**
for methods, input keywords, and build options, and the
**[OpenQP tutorials](https://open-quantum-platform.github.io/openqp-tutorials/)**
for runnable examples.

### Functionality

#### Electronic-Structure Methods

| Method | References / variants | Notes |
| --- | --- | --- |
| **MRSF-TDDFT** | [Mixed-reference spin-flip](https://doi.org/10.1021/acs.jpclett.3c02296) with [DTCAM-series functionals](https://doi.org/10.1021/acs.jctc.4c00640) | Spin-adapted multiconfigurational ground and excited states within linear response |
| UMRSF-TDDFT | MRSF excitation energies from a UHF reference | Energy-only |
| MRSF-EKT | [IP/EA via Extended Koopmans' Theorem](https://doi.org/10.1021/acs.jpclett.1c02494) | Dyson orbitals and pole strengths (`runtype=ekt`) |
| Hartree–Fock | RHF, ROHF, UHF | Closed- and open-shell SCF foundations |
| DFT | RKS / UKS / ROKS via [LibXC](https://gitlab.com/libxc/libxc) | Hundreds of LCAO functionals; range-separated (CAM/LRC) support |
| MP2 | RHF, UHF, and ROHF energies; MP2, SCS-MP2, SOS-MP2, OS/SS-MP2, SCS-MI-MP2, and custom spin scaling | Analytic nuclear gradients for RHF references; UHF/ROHF nuclear derivatives are not yet available |
| Coupled cluster | CCSD and CCSD(T) with RHF, UHF, and ROHF references | Energy-only; frozen-core, in-core, Cholesky-factorized, and integral-direct options are available as appropriate to the reference |
| TDHF / TDDFT | RPA, TDA | Conventional linear-response excited states |
| SF-TDDFT | Spin-flip TDA | Spin-flip excited states from a high-spin reference |
| CI and CASSCF | FCI, CASCI, CASSCF, SA-CASSCF | Native determinant CI and orbital optimization; CASSCF/SA-CASSCF have central-difference nuclear gradients |
| Multireference PT2 | CASPT2 (SS/MS/XMS), NEVPT2 (uncontracted/strongly contracted), QDPT2 (MRMP2/MCQDPT2/XMCQDPT2) | Native energy and central-difference nuclear-gradient calculations; see [`examples/WF_methods`](examples/WF_methods) for variants, controls, and present limits |

**Tutorials:** [MRSF-TDDFT](https://open-quantum-platform.github.io/openqp-tutorials/mrsf-tddft/) · [UMRSF-TDDFT](https://open-quantum-platform.github.io/openqp-tutorials/umrsf-tddft/) · [Spin-flip TDDFT](https://open-quantum-platform.github.io/openqp-tutorials/sf-tddft/) · [TDDFT/TDHF](https://open-quantum-platform.github.io/openqp-tutorials/tddft-and-tdhf/) · [Hartree–Fock & DFT](https://open-quantum-platform.github.io/openqp-tutorials/hf-and-dft/) · [MP2 & spin-scaled MP2](https://open-quantum-platform.github.io/openqp-tutorials/mp2/)

#### Properties & Spectroscopy

| Capability | Scope | Notes |
| --- | --- | --- |
| Analytic gradients | HF, DFT, RHF-based MP2 variants, TDDFT, SF/MRSF-TDDFT | State-specific gradients for optimization and dynamics |
| Numerical gradients | CASSCF, SA-CASSCF, CASPT2, NEVPT2, QDPT2 | Central differences of converged multireference energies for `grad`, `optimize`, `ts`, `mep`, and `irc` |
| Hessians | Native **analytic** HF/DFT Hessians + numerical Hessians | Covers UHF/ROHF references, ECPs, and CAM/LRC functionals |
| Vibrational analysis | Frequencies, normal modes, thermochemistry, **IR and Raman intensities** | Native dipole / CPHF-polarizability kernels |
| **NMR shieldings** | CGO and GIAO (London-orbital) gauges | HF and DFT, closed- and open-shell |
| **Nonadiabatic couplings** | NAC / NACME between MRSF-TDDFT states | [TLF technology](https://doi.org/10.1021/acs.jpclett.1c00932) for dynamics workflows |
| **Spin-orbit coupling** | SOC between MRSF-TDDFT states | One- and two-electron contributions ([Relativistic MRSF-TDDFT](https://doi.org/10.1021/acs.jctc.2c01036)) |
| **X-ray absorption** | XAS / core-excitation workflows (incl. ΔCHP-MRSF) | Core-level excited states |
| **Implicit solvation** | PCM via the ddX backend (ddCOSMO / ddPCM / ddLPB) | Energy-only continuum solvent on RHF/ROHF references |
| Population & moments | Mulliken, Löwdin, RESP charges; electric multipole moments | `runtype=prop` |
| Dispersion | [DFT-D4](https://dftd4.readthedocs.io/en/latest/) correction | — |

**Tutorials:** [Hessians, frequencies & IR/Raman](https://open-quantum-platform.github.io/openqp-tutorials/vibrational-analysis/) · [NMR shielding](https://open-quantum-platform.github.io/openqp-tutorials/nmr-shielding/) · [Spin–orbit coupling](https://open-quantum-platform.github.io/openqp-tutorials/spin-orbit-coupling/) · [Population, moments & charges](https://open-quantum-platform.github.io/openqp-tutorials/properties-and-population/) · [PCM/ddX solvation](https://open-quantum-platform.github.io/openqp-tutorials/pcm-solvation/)

#### Geometry & Reaction Paths

| Workflow | `runtype` | Default execution |
| --- | --- | --- |
| Energy / gradient / Hessian | `energy`, `grad`, `hess` | native |
| Minimization & transition states | `optimize`, `ts` | native OQP |
| Conical intersections | `meci`, `mecp`, `tci` | native OQP |
| Reaction paths | `irc`, `mep`, `neb` | native OQP |
| Nonadiabatic data | `nac`, `nacme` | native |

The built-in optimizer uses redundant-internal / DLC / TRIC coordinates with
restricted-step RFO/P-RFO and needs no external optimizer package.  It covers
all primary geometry and reaction-path workflows above, including aligned,
endpoint-optimized climbing-image NEB and optional numerical/analytical initial
Hessians for transition-state searches. SciPy and
[geomeTRIC](https://github.com/leeping/geomeTRIC) remain optional compatibility
backends for advanced constrained optimization.

**Tutorials:** [Geometry optimization & TS](https://open-quantum-platform.github.io/openqp-tutorials/geometry-optimization/) · [Conical intersections](https://open-quantum-platform.github.io/openqp-tutorials/conical-intersections/)

#### Dynamics & QM/MM

Nonadiabatic molecular dynamics (`runtype=namd`) with Tully fewest-switches
surface hopping on MRSF-TDDFT states, spin-orbit-coupled intersystem crossing,
and ESPF QM/MM embedding. These compose into **[SOC-NAMD-QMMM](https://open-quantum-platform.github.io/openqp-tutorials/soc-namd-qmmm/)**: excited-state
surface-hopping dynamics of an MRSF-TDDFT chromophore, with singlet/triplet
intersystem crossing, embedded in an explicit (OpenMM) MM solvent.

- Native fewest-switches surface hopping (`runtype=namd`) for gas-phase MRSF-TDDFT internal conversion.
- SOC-NAMD for intersystem crossing: SHARC-like spin-adiabatic propagation and an MCH-basis mode with exact active-root MCH gradients (`[md] soc_basis=mch`).
- ESPF electrostatic QM/MM via OpenMM (PME periodic electrostatics, smooth ESPF grid forces, rigid-water constraints); QM/MM composes with both FSSH and SOC-NAMD to give SOC-NAMD-QMMM.
- Overlap-based MRSF state tracking, finite-difference NAC/TDC propagation, and energy-based decoherence (EDC).

**Tutorials:** [SOC-NAMD-QMMM](https://open-quantum-platform.github.io/openqp-tutorials/soc-namd-qmmm/) · [ESPF QM/MM embedding](https://open-quantum-platform.github.io/openqp-tutorials/qmmm-embedding/)

See the **SOC-NAMD-QMMM** guide and the `[md]` / `[qmmm]` keyword pages in the
[manual](https://open-quantum-platform.github.io/openqp-docs/workflows/soc-namd-qmmm/)
for complete input decks and the compact `job.qmmm(...)` / `job.workflow.namd(...)`
Python API.

#### SCF, Initial Guesses & Performance

| Area | What OpenQP provides |
| --- | --- |
| Initial guesses | Native `hcore`, `huckel`, `modhuckel`, `minao`, `sap`; `json` restart and `auto`; optional PySCF (`sad`/`sap`/`pyscf`) guesses |
| SCF convergence | DIIS family (C/E/A/V-DIIS), SOSCF, and OpenQP's **own native TRAH** (Trust-Region Augmented Hessian) solver, with the external [OpenTrustRegion](https://github.com/eriksen-lab/opentrustregion) library as an optional alternative |
| Symmetry | Point-group detection; MO/state/mode labels; petite-list reductions accelerating integrals, XC, gradients, and response |
| DFT grids | Lebedev plus SG-0/SG-1/SG-2/SG-3 pruned grids with per-element DE2 radial quadrature; OpenMP-parallel XC kernels |
| Excited-state robustness | Davidson auto-restart; MINRES/AUTO Z-vector fallbacks |
| Parallelism & deployment | OpenMP and MPI; BLAS/LAPACK optimization; pip install and Docker images |
| BLAS/LAPACK threading | Parallel regions declare their own width, so a global BLAS setter cannot serialize them; dense eigensolves size their thread count to the matrix, since LAPACK's eigensolvers do not thread the way GEMM does — see [docs/blas_threading.md](docs/blas_threading.md) |

**Tutorials:** [SCF convergence & guesses](https://open-quantum-platform.github.io/openqp-tutorials/scf-convergence/) · [Effective core potentials](https://open-quantum-platform.github.io/openqp-tutorials/effective-core-potentials/)

#### Ecosystem & Integrations

| Integration | Purpose |
| --- | --- |
| [LibXC](https://gitlab.com/libxc/libxc) | Wide library of exchange-correlation functionals |
| [basis_set_exchange](https://github.com/MolSSI-BSE/basis_set_exchange) | Standard basis sets |
| [libecpint](https://github.com/robashaw/libecpint) | Effective Core Potentials |
| [DFT-D4](https://dftd4.readthedocs.io/en/latest/) | Dispersion correction |
| [PyRAI2MD](https://github.com/mlcclab/PyRAI2MD-hiam) | AI-driven ab initio molecular dynamics |
| [Molden](https://www.theochem.ru.nl/molden/) format | Standards-oriented geometry, basis, SCF/Dyson orbitals, and optional frequency sections for common graphics tools |
| [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) | Browser-based inspection of log, JSON, Molden, cube, and XYZ outputs |
| Optional [DFTB+](https://dftbplus.org/) backend | Ground-state energy, gradient, and geometry optimization |
| Optional [MOKIT](https://github.com/1234zou/MOKIT) | Broader external wavefunction conversion workflows |

### Upcoming Features
- Full analytic spin-adiabatic SOC gradients, requiring MCH derivative-coupling vectors and SOC-gradient matrix elements.
- **Scalar-relativistic (X2C) framework** extending the relativistic MRSF-TDDFT treatment

### Install

```bash
pip install openqp
```

For a source checkout:

```bash
git clone https://github.com/Open-Quantum-Platform/openqp.git
cd openqp
pip install .
```

The package install keeps the Python wrapper, native library, headers, and data files together for normal `openqp` command-line use. A ready-to-use [Docker image](https://github.com/Open-Quantum-Platform/openqp/wiki/OpenQP_Docker_Image) is also available. The image is a distribution of OpenQP, so the same research/commercial licensing terms apply to its use. Build options (MPI, LibXC/ERI backends, BLAS/LAPACK selection) are documented in the [Build options](https://open-quantum-platform.github.io/openqp-docs/build-options/) guide.

### First Run

```bash
openqp examples/HF/H2O_RHF-HF_ENERGY.oqp          # OpenMP / sequential run
```

Every legacy example has a concise `.oqp` companion. To run the regression
suite exclusively with concise inputs, use `--input-format`:

```bash
openqp --run_tests all --input-format oqp         # standard suite through .oqp
```

Omitting the selector uses `auto`: it prefers `.oqp`, retains any legacy-only
inputs, and keeps a small representative legacy compatibility set. The
historical `all` scope still excludes unusually slow or non-self-contained
examples; selecting an explicit directory applies the requested format to every
input in that directory. Each calculation receives an isolated output folder,
so paired optimization artifacts cannot overwrite one another.

Control OpenMP threads per process or MPI rank with `--omp 16` or `[input] omp_threads=16`.

### Documentation

- [OpenQP Manual](https://open-quantum-platform.github.io/openqp-docs/) (reference docs; source: [openqp-docs](https://github.com/Open-Quantum-Platform/openqp-docs))
- [OpenQP Tutorials](https://open-quantum-platform.github.io/openqp-tutorials/) (guided, runnable walkthroughs)
- [Build options](https://open-quantum-platform.github.io/openqp-docs/build-options/)
- [API guide](https://open-quantum-platform.github.io/openqp-docs/api/)
- [Example inputs](examples)
- [BLAS/LAPACK threading](docs/blas_threading.md) — how OpenQP's OpenMP regions and a threaded BLAS are kept out of each other's way, and how to check which BLAS you actually linked

### Graphic Web Tools

- [OpenQP Web](https://app.openqp.org/) — prepare inputs and preview structures locally in the browser.
- [OpenQP Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/) — browser-based input builder.
- [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) — inspect OpenQP log, JSON, Molden, cube, and XYZ outputs in the browser; files are processed locally and never uploaded.

Current full JSON output includes a portable Molden-ordered AO basis and SCF orbital block. Hessian JSON adds the same orbital data beside frequencies and normal modes, and MRSF-EKT JSON adds state-specific IP/EA Dyson orbitals transformed to the AO basis. Frequency Molden output combines `[Atoms]`, `[GTO]`, `[MO]`, `[FREQ]`, standard one-value-per-mode `[INT]`, optional `[RAMAN]`, `[FR-COORD]`, and `[FR-NORM-COORD]` sections in one file. EKT runs with `save_molden=True` also write a `dyson` Molden file whose labeled orbital records correspond to individual IP/EA roots. The small `examples/HESS/H2O_RHF-DFT_VIEWER_EXPORT` example exercises the combined JSON/Molden export path.

### Citing OpenQP
If you use OpenQP in your research, please cite the OpenQP platform paper:

- **Mironov V, Komarov K, Li J, Gerasimov I, Mazaheri M, Park W, Lashkaripour A, Oh M, Nakata H, Ishimura K, Huix-Rotllant M, Lee S, and Choi CH.** "OpenQP: A Quantum Chemical Platform Featuring MRSF-TDDFT with an Emphasis on Open-source Ecosystem" [Journal of Chemical Theory and Computation, 2024](https://doi.org/10.1021/acs.jctc.4c01117)

Original MRSF-TDDFT theory and analytic-gradient papers:

- **Lee S, Filatov M, Lee S, and Choi CH.** "Eliminating Spin-Contamination of Spin-Flip Time-Dependent Density Functional Theory Within Linear Response Formalism by the Use of Zeroth-Order Mixed-Reference (MR) Reduced Density Matrix." [The Journal of Chemical Physics, vol. 149, no. 10, 2018.](https://doi.org/10.1063/1.5044202)
- **Lee S, Kim EE, Nakata H, Lee S, and Choi CH.** "Efficient Implementations of Analytic Energy Gradient for Mixed-Reference Spin-Flip Time-Dependent Density Functional Theory (MRSF-TDDFT)." [The Journal of Chemical Physics, vol. 150, no. 18, 2019.](https://doi.org/10.1063/1.5086895)

Recent MRSF-TDDFT accounts and overview papers:

- **Park W, Komarov K, Lee S, and Choi CH.** "Mixed-Reference Spin-Flip Time-Dependent Density Functional Theory: Multireference Advantages with the Practicality of Linear Response Theory." [The Journal of Physical Chemistry Letters. 2023 Sep 28;14(39):8896-908.](https://doi.org/10.1021/acs.jpclett.3c02296)
- **Lee S, Park W, and Choi CH.** "Expanding Horizons in Quantum Chemical Studies: The Versatile Power of MRSF-TDDFT." [Accounts of Chemical Research, 2025.](https://doi.org/10.1021/acs.accounts.4c00640)
- **Park W, Lee S, Komarov K, Mironov V, Nakata H, Zeng T, Huix-Rotllant M, and Choi CH.** "MRSF-TDDFT: A New Tool in Quantum Chemistry for Better Understanding Molecules and Materials." [Bulletin of the Korean Chemical Society, 2025.](https://doi.org/10.1002/bkcs.70011)

### Contributors

**Principal Investigator**

- **Cheol Ho Choi** (PI), Kyungpook National University, South Korea, [cheolho.choi@gmail.com](mailto:cheolho.choi@gmail.com), [https://www.openqp.org](https://www.openqp.org)

**Development team**

- **Seunghoon Lee**, Seoul National University, South Korea, [seunghoonlee89@gmail.com](mailto:seunghoonlee89@gmail.com)
- **Vladimir Mironov**, [vladimir.a.mironov@gmail.com](mailto:vladimir.a.mironov@gmail.com)
- **Konstantin Komarov**, [constlike@gmail.com](mailto:constlike@gmail.com)
- **Jingbai Li**, Hoffmann Institute of Advanced Materials, China, [lijingbai2009@gmail.com](mailto:lijingbai2009@gmail.com)
- **Igor Gerasimov**, [i.s.ger@yandex.ru](mailto:i.s.ger@yandex.ru)
- **Hiroya Nakata**, Fukui Institute for Fundamental Chemistry, Japan, [nakata.hiro07@gmail.com](mailto:nakata.hiro07@gmail.com)
- **Mohsen Mazaherifar**, Kyungpook National University, South Korea, [moh.mazaheri@gmail.com](mailto:moh.mazaheri@gmail.com)
- **Vladimir Makhnev** ([VladimirMakhnev](https://github.com/VladimirMakhnev))
- **[Alireza Lashkaripour](https://github.com/Alireza-Lashkaripour)**

### Legal Notice

OpenQP v1.3.0 and later versions distributed with the current license are
source-available under a dual-licensing model:
qualifying academic and nonprofit users may use OpenQP under the free OpenQP
Research License 1.0, while every use by or for a commercial entity requires a
separate written commercial license from Open Quantum Inc. OpenQP v1.2.1 and
earlier release tags remain under the GPLv3 terms that accompanied them, but those
terms do not extend to current or future versions distributed under the new
license. See the complete [license](LICENSE), [licensing guide](LICENSING.md),
and [community and sustainable development statement](SUSTAINABILITY.md).
