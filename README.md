## Open Quantum Platform: OpenQP

Open Quantum Platform ([OpenQP](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01117)) is a quantum chemical platform built around [Mixed-Reference Spin-Flip (MRSF)-TDDFT](https://doi.org/10.1021/acs.jpclett.3c02296) with an emphasis on an open-source ecosystem. It combines conventional HF/DFT and TDHF/TDDFT with MRSF-TDDFT to treat multiconfigurational ground and excited states — diradicals, bond breaking, conical intersections, nonadiabatic dynamics, and spin-orbit coupling — using autonomous, interoperable modules driven through the **PyOQP** Python wrapper.

MRSF-TDDFT is the central scientific feature of OpenQP: it retains the practical linear-response structure of TDDFT while removing the spin contamination that limits conventional spin-flip TDDFT, making it useful for multiconfigurational ground-state surfaces as well as excited-state and photochemical workflows.

### Functionality

#### Electronic-Structure Methods

| Method | References / variants | Notes |
| --- | --- | --- |
| Hartree–Fock | RHF, ROHF, UHF | Closed- and open-shell SCF foundations |
| DFT | RKS / UKS / ROKS via [LibXC](https://gitlab.com/libxc/libxc) | Hundreds of LCAO functionals; range-separated (CAM/LRC) support |
| TDHF / TDDFT | RPA, TDA | Conventional linear-response excited states |
| SF-TDDFT | Spin-flip TDA | Spin-flip excited states from a high-spin reference |
| **MRSF-TDDFT** | [Mixed-Reference Spin-Flip](https://doi.org/10.1021/acs.jpclett.3c02296) + [DTCAM-series functionals](https://doi.org/10.1021/acs.jctc.4c00640) | Main production method; multireference accuracy with LR practicality |
| UMRSF-TDDFT | MRSF excitation energies from a UHF reference | Energy-only |
| MRSF-EKT | [IP/EA via Extended Koopmans' Theorem](https://doi.org/10.1021/acs.jpclett.1c02494) | Dyson orbitals and pole strengths (`runtype=ekt`) |

#### Properties & Spectroscopy

| Capability | Scope | Notes |
| --- | --- | --- |
| Analytic gradients | HF, DFT, TDDFT, SF/MRSF-TDDFT | State-specific gradients for optimization and dynamics |
| Hessians | Native **analytic** HF/DFT Hessians + numerical Hessians | Covers UHF/ROHF references, ECPs, and CAM/LRC functionals |
| Vibrational analysis | Frequencies, normal modes, thermochemistry, **IR and Raman intensities** | Native dipole / CPHF-polarizability kernels |
| **NMR shieldings** | CGO and GIAO (London-orbital) gauges | HF and DFT, closed- and open-shell |
| **Nonadiabatic couplings** | NAC / NACME between MRSF-TDDFT states | [TLF technology](https://doi.org/10.1021/acs.jpclett.1c00932) for dynamics workflows |
| **Spin-orbit coupling** | SOC between MRSF-TDDFT states | One- and two-electron contributions ([Relativistic MRSF-TDDFT](https://doi.org/10.1021/acs.jctc.2c01036)) |
| **X-ray absorption** | XAS / core-excitation workflows (incl. ΔCHP-MRSF) | Core-level excited states |
| **Implicit solvation** | PCM via the ddX backend (ddCOSMO / ddPCM / ddLPB) | Energy-only continuum solvent on RHF/ROHF references |
| Population & moments | Mulliken, Löwdin, RESP charges; electric multipole moments | `runtype=prop` |
| Dispersion | [DFT-D4](https://dftd4.readthedocs.io/en/latest/) correction | — |

#### Geometry & Reaction Paths

| Workflow | `runtype` | Backends |
| --- | --- | --- |
| Energy / gradient / Hessian | `energy`, `grad`, `hess` | native |
| Minimization & transition states | `optimize`, `ts` | `oqp` (native), [geomeTRIC](https://github.com/leeping/geomeTRIC), SciPy |
| Conical intersections | `meci`, `mecp`, `tci` | `oqp`, geomeTRIC, SciPy |
| Reaction paths | `irc`, `mep`, `neb` | `oqp`, geomeTRIC, SciPy |
| Nonadiabatic data | `nac`, `nacme` | native |

The built-in native optimizer (`lib=oqp`) uses redundant-internal / DLC / TRIC coordinates with a restricted-step RFO step and needs no external optimizer package.

#### SCF, Initial Guesses & Performance

| Area | What OpenQP provides |
| --- | --- |
| Initial guesses | Native `hcore`, `huckel`, `modhuckel`, `minao`, `sap`; `json` restart and `auto`; optional PySCF (`sad`/`sap`/`pyscf`) guesses |
| SCF convergence | DIIS family (C/E/A/V-DIIS), SOSCF, and OpenQP's **own native TRAH** (Trust-Region Augmented Hessian) solver, with the external [OpenTrustRegion](https://github.com/eriksen-lab/opentrustregion) library as an optional alternative |
| Symmetry | Point-group detection; MO/state/mode labels; petite-list reductions accelerating integrals, XC, gradients, and response |
| DFT grids | Lebedev plus SG-0/SG-1/SG-2/SG-3 pruned grids with per-element DE2 radial quadrature; OpenMP-parallel XC kernels |
| Excited-state robustness | Davidson auto-restart; MINRES/AUTO Z-vector fallbacks |
| Parallelism & deployment | OpenMP and MPI; BLAS/LAPACK optimization; pip install and Docker images |

#### Ecosystem & Integrations

| Integration | Purpose |
| --- | --- |
| [LibXC](https://gitlab.com/libxc/libxc) | Wide library of exchange-correlation functionals |
| [basis_set_exchange](https://github.com/MolSSI-BSE/basis_set_exchange) | Standard basis sets |
| [libecpint](https://github.com/robashaw/libecpint) | Effective Core Potentials |
| [DFT-D4](https://dftd4.readthedocs.io/en/latest/) | Dispersion correction |
| [PyRAI2MD](https://github.com/mlcclab/PyRAI2MD-hiam) | AI-driven ab initio molecular dynamics |
| [Molden](https://www.theochem.ru.nl/molden/) format | Visualization compatible with common graphics tools |
| [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) | Browser-based inspection of log, JSON, Molden, cube, and XYZ outputs |
| Optional [DFTB+](https://dftbplus.org/) backend | Ground-state energy, gradient, and geometry optimization |
| Optional [MOKIT](https://github.com/1234zou/MOKIT) | Broader external wavefunction conversion workflows |

### Upcoming Features
- **Efficient electrostatic embedding QM/MM** by [ESPF QM/MM](https://doi.org/10.1063/5.0133646)
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

The package install keeps the Python wrapper, native library, headers, and data files together for normal `openqp` command-line use. A ready-to-use [Docker image](https://github.com/Open-Quantum-Platform/openqp/wiki/OpenQP_Docker_Image) is also available. Build options (MPI, LibXC/ERI backends, BLAS/LAPACK selection) are documented in the [Build options](https://open-quantum-platform.github.io/openqp-docs/build-options/) guide.

### First Run

```bash
openqp examples/HF/H2O_RHF-HF_ENERGY.inp          # OpenMP / sequential run
mpirun -np <n> openqp any_example_file.inp        # MPI run
openqp --run_tests all                            # run the packaged example tests
```

Control OpenMP threads per process or MPI rank with `--omp 16` or `[input] omp_threads=16`.

### Documentation

- [OpenQP Manual](https://open-quantum-platform.github.io/openqp-docs/)
- [Build options](https://open-quantum-platform.github.io/openqp-docs/build-options/)
- [API guide](https://open-quantum-platform.github.io/openqp-docs/api/)
- [Example inputs](examples)

### Graphic Web Tools

- [OpenQP Web](https://app.openqp.org/) — prepare inputs and preview structures locally in the browser.
- [OpenQP Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/) — browser-based input builder.
- [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) — inspect OpenQP log, JSON, Molden, cube, and XYZ outputs in the browser; files are processed locally and never uploaded.

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

### Legal Notice

See the separate LICENSE file.
