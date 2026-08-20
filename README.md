# OpenQP

Open Quantum Platform ([OpenQP](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01117);
abbreviated **OQP** and pronounced “oh-coop”) is a source-available quantum
chemistry program centered on
[mixed-reference spin-flip TDDFT (MRSF-TDDFT)](https://doi.org/10.1021/acs.jpclett.3c02296)
that combines HF and DFT references, correlated wavefunction methods,
linear-response excited states, nuclear derivatives, reaction paths,
nonadiabatic dynamics, and QM/MM calculations in one program.

MRSF-TDDFT retains the linear-response structure of TDDFT while removing the spin
contamination that limits conventional spin-flip TDDFT. This permits consistent
descriptions of multiconfigurational ground and excited states, including
diradicals, bond breaking, conical intersections, photochemical processes, and
spin-orbit coupling.

Use the **[OpenQP tutorials](https://open-quantum-platform.github.io/openqp-tutorials/)**
for guided calculations and the
**[OpenQP manual](https://open-quantum-platform.github.io/openqp-docs/)**
for method, keyword, and build references.

### Functionality

#### Electronic-Structure Methods

| Method | References / variants | Available calculations | Learn |
| --- | --- | --- | --- |
| **MRSF-TDDFT** | ROHF mixed reference with [DTCAM-series functionals](https://doi.org/10.1021/acs.jctc.4c00640) | Ground and excited states, analytic gradients, properties, NAC/NACME, SOC, optimization, and dynamics | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/mrsf-tddft/) |
| **UMRSF-TDDFT** | Unrestricted mixed reference based on UHF orbitals | Excitation energies for unrestricted mixed-reference calculations | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/umrsf-tddft/) |
| **MRSF-EKT** | [Extended Koopmans' theorem](https://doi.org/10.1021/acs.jpclett.1c02494) applied to MRSF states | Ionization and electron-attachment energies, Dyson orbitals, and pole strengths | [Guide](https://open-quantum-platform.github.io/openqp-docs/workflows/ekt/) |
| **SF-TDDFT** | Spin-flip TDA from a high-spin reference | Excitation energies and analytic state-specific nuclear gradients | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/sf-tddft/) |
| TDHF / TDDFT | RPA and TDA response from HF and DFT references | Excitation energies, analytic state-specific nuclear gradients, and excited-state properties | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/tddft-and-tdhf/) |
| Hartree–Fock | RHF, ROHF, and UHF references | Energies, analytic nuclear gradients and Hessians, and molecular properties | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/hf-and-dft/) |
| Density functional theory | RKS, ROKS, and UKS through [LibXC](https://gitlab.com/libxc/libxc) | Energies, analytic nuclear gradients and Hessians, properties, and CAM/LRC functionals | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/hf-and-dft/) |
| MP2 | RHF, ROHF, and UHF with conventional and spin-scaled variants | Correlation energies for all references and analytic nuclear gradients for RHF-based variants | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/mp2/) |
| Coupled cluster | RHF-, ROHF-, and UHF-based CCSD and CCSD(T) | Correlation energies with reference-appropriate frozen-core, Cholesky, in-core, or direct algorithms | [Guide](https://open-quantum-platform.github.io/openqp-docs/workflows/coupled-cluster/) |
| Configuration interaction | FCI and CASCI in determinant active spaces | State energies and CI vectors at fixed molecular orbitals | [Examples](examples/WF_methods) |
| CASSCF | State-specific CASSCF and state-averaged SA-CASSCF | Energies and analytic nuclear gradients: state specific, the weighted SA-CASSCF objective, and an individual averaged root through the coupled orbital/CI Z-vector; the legacy `method=casscf` plus `[state_average]` spelling stays on central differences | [Guide](https://open-quantum-platform.github.io/openqp-docs/workflows/casscf-gradient/) |
| Multireference PT2 | SS/MS/XMS-CASPT2, `nevpt2` / `sc-nevpt2`, and MRMP2/MCQDPT2/XMCQDPT2 | Correlated energies and central-difference nuclear gradients for supported geometry workflows | [Examples](examples/WF_methods) |

#### Capabilities

| Capability | Scope | Available calculations | Learn |
| --- | --- | --- | --- |
| Analytic nuclear gradients | MRSF/SF/TDDFT, HF/DFT, RHF-based MP2, and state-specific CASSCF and SA-CASSCF | State-specific derivatives for `grad`, optimization, and reaction paths; `runtype=namd` dynamics is MRSF-TDDFT only, and the CASSCF/SA-CASSCF and PT2 routes accept grad/optimize/ts/mep/irc; SA-CASSCF additionally differentiates the weighted objective, which is not a state and is selected by `[casscf] gradient_state` rather than by `[properties] grad` | [Overview](https://open-quantum-platform.github.io/openqp-docs/capabilities/) |
| Numerical nuclear gradients | The legacy `method=casscf` plus `[state_average]` spelling, and the CASPT2/NEVPT2/QDPT2 families | Central differences of converged energies for `grad`, `optimize`, `ts`, `mep`, and `irc` | [Examples](examples/WF_methods) |
| Hessians and vibrations | Analytic HF/DFT Hessians and numerical Hessians from supported gradients | Frequencies, normal modes, thermochemistry, IR intensities, and Raman intensities | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/vibrational-analysis/) |
| NMR shieldings | HF and DFT with CGO or GIAO (London-orbital) gauges | Closed- and open-shell nuclear magnetic shielding tensors | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/nmr-shielding/) |
| Nonadiabatic couplings | NAC and NACME between MRSF-TDDFT states | State-pair couplings for optimization and nonadiabatic dynamics | [Guide](https://open-quantum-platform.github.io/openqp-docs/workflows/nacme/) |
| Spin-orbit coupling | One- and two-electron SOC between MRSF-TDDFT states | Spin-mixed energies and coupling matrices for intersystem crossing | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/spin-orbit-coupling/) |
| X-ray absorption | Core-excitation workflows including ΔCHP-MRSF | Core-level excitation energies and oscillator strengths | [Examples](examples/XAS) |
| Implicit solvation | ddCOSMO, ddPCM, and ddLPB through the ddX backend | Continuum-solvent energies for RHF and ROHF references | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/pcm-solvation/) |
| Populations and moments | Mulliken, Löwdin, and RESP analyses | Atomic charges and electric multipole moments from `runtype=prop` | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/properties-and-population/) |
| Dispersion correction | [DFT-D4](https://dftd4.readthedocs.io/en/latest/) with supported DFT references | Dispersion-corrected energies, gradients, and Hessians | [Example](examples/DFT/H2O_RHF-PBE_D4_EXPLICIT.oqp) |
| Geometry optimization and transition states | Native redundant-internal, DLC, and TRIC coordinates with RFO/P-RFO steps | Minima and transition states through `runtype=optimize` and `runtype=ts` | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/geometry-optimization/) |
| Conical intersections | MECI, MECP, and three-state intersection searches | Crossing structures through `runtype=meci`, `mecp`, and `tci` | [Tutorial](https://open-quantum-platform.github.io/openqp-tutorials/conical-intersections/) |
| Reaction paths | IRC, MEP, and climbing-image NEB calculations | Intrinsic reaction coordinates, minimum-energy paths, and endpoint-connected paths | [Guide](https://open-quantum-platform.github.io/openqp-docs/workflows/optimization/) |
| Nonadiabatic dynamics | MRSF-TDDFT fewest-switches surface hopping with state tracking and decoherence | Internal-conversion trajectories through `runtype=namd` | [Examples](examples/ODP) |
| QM/MM and SOC-NAMD-QMMM | ESPF electrostatic embedding through OpenMM with MRSF-TDDFT | Embedded energies and forces, QM/MM dynamics, and spin-orbit-coupled surface hopping | [QM/MM tutorial](https://open-quantum-platform.github.io/openqp-tutorials/qmmm-embedding/) · [Dynamics tutorial](https://open-quantum-platform.github.io/openqp-tutorials/soc-namd-qmmm/) |

#### Ecosystem & Integrations

| Integration | Purpose |
| --- | --- |
| [LibXC](https://gitlab.com/libxc/libxc) | Wide library of exchange-correlation functionals |
| [basis_set_exchange](https://github.com/MolSSI-BSE/basis_set_exchange) | Standard basis sets |
| [libecpint](https://github.com/robashaw/libecpint) | Effective Core Potentials |
| [DFT-D4](https://dftd4.readthedocs.io/en/latest/) | Dispersion correction |
| [PyRAI2MD](https://github.com/mlcclab/PyRAI2MD-hiam) | Machine-learning-assisted nonadiabatic molecular dynamics |
| [Molden](https://www.theochem.ru.nl/molden/) format | Standards-oriented geometry, basis, SCF/Dyson orbitals, and optional frequency sections for common graphics tools |
| [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) | Browser-based inspection of log, JSON, Molden, cube, and XYZ outputs |
| Optional [MOKIT](https://github.com/1234zou/MOKIT) | Broader external wavefunction conversion workflows |

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

Windows, macOS and Linux are all supported; per-platform requirements and
source-build options are documented in the
[installation guide](https://open-quantum-platform.github.io/openqp-docs/installation/)
and the [build guide](https://open-quantum-platform.github.io/openqp-docs/build-options/).
A ready-to-use [Docker image](https://github.com/Open-Quantum-Platform/openqp/wiki/OpenQP_Docker_Image)
is also available.

#### Standalone builds (no Python required)

Self-contained archives for Windows, macOS and Linux — unzip and run, with no
Python, compiler or BLAS on the machine — are produced and published by
[OQP Studio](https://github.com/Open-Quantum-Platform/oqp-studio), which
rebuilds them weekly from OpenQP's latest release tag and ships them with the
desktop application:

```bash
./openqp examples/HF/H2O_RHF-HF_ENERGY.inp
```

`pip install openqp` remains the usual route when a Python environment is
available; the archives exist for machines without one, for clusters where pip
is not an option, and for embedding in other applications.

### Quick Start

```bash
openqp examples/HF/H2O_RHF-HF_ENERGY.oqp
```

```bash
openqp --run_tests all --input-format oqp
```

See the [quick-start guide](https://open-quantum-platform.github.io/openqp-docs/quickstart/)
for calculation setup and the
[test-runner reference](https://open-quantum-platform.github.io/openqp-docs/keywords/tests/)
for selectors, input formats, and parallel execution.

### Documentation

- [OpenQP Manual](https://open-quantum-platform.github.io/openqp-docs/) (reference docs; source: [openqp-docs](https://github.com/Open-Quantum-Platform/openqp-docs))
- [OpenQP Tutorials](https://open-quantum-platform.github.io/openqp-tutorials/) (guided, runnable walkthroughs)
- [SCF](https://open-quantum-platform.github.io/openqp-docs/keywords/scf/), [initial guesses](https://open-quantum-platform.github.io/openqp-docs/keywords/guess/), and [performance](https://open-quantum-platform.github.io/openqp-docs/performance/)
- [Build options](https://open-quantum-platform.github.io/openqp-docs/build-options/)
- [API guide](https://open-quantum-platform.github.io/openqp-docs/api/)
- [Example inputs](examples)

### Graphic Web Tools

- [OpenQP Web](https://app.openqp.org/) — prepare inputs and preview structures locally in the browser.
- [OpenQP Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/) — browser-based input builder.
- [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) — inspect OpenQP log, JSON, Molden, cube, and XYZ outputs in the browser; files are processed locally and never uploaded.
- [Orbital and vibrational output guide](https://open-quantum-platform.github.io/openqp-docs/workflows/orbital-vibrational-output/) — JSON and Molden exports for visualization and analysis.

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
