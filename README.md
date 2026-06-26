# OpenQP

Open Quantum Platform (OpenQP) is an open-source quantum chemistry package built
around Mixed-Reference Spin-Flip TDDFT (MRSF-TDDFT) and practical workflows for
multiconfigurational ground and excited states. It combines ordinary HF/DFT and
TDHF/TDDFT capabilities with tools for diradicals, bond breaking, conical
intersections, nonadiabatic couplings, and spin-orbit coupling.

MRSF-TDDFT is the central scientific feature of OpenQP. It keeps the practical
linear-response structure of TDDFT while reducing the spin contamination that
limits conventional spin-flip TDDFT, making it useful for multiconfigurational
ground-state surfaces as well as excited-state and photochemical workflows.

## Highlights

- MRSF-TDDFT as the main production method for multiconfigurational ground and
  excited states, designed to reduce spin contamination while retaining
  TDDFT-like practicality
- MRSF-TDDFT energies, analytic gradients, geometry optimization, MECI/MECP/TCI,
  NACME, spin-orbit coupling, and MRSF-EKT workflows
- DTCAM-series functionals, UMRSF-TDDFT, and newer MRSF extensions for broader
  multistate coverage
- Practical photochemistry targets: diradicals, bond breaking, conical
  intersections, nonadiabatic dynamics, X-ray absorption, SOC, and inverted
  singlet-triplet materials
- HF/DFT, TDHF/TDDFT, and SF-TDDFT foundations with RHF, ROHF, and UHF
  references
- Native optimizer (`lib=oqp`), geomeTRIC/SciPy backends, HF/DFT Hessians,
  NMR/IR/Raman, PCM/ddX energies, OpenMP/MPI, pip installs, and Docker support

## Install

```bash
pip install openqp
```

For a source checkout:

```bash
git clone https://github.com/Open-Quantum-Platform/openqp.git
cd openqp
pip install .
```

The package install keeps the Python wrapper, native library, headers, and data
files together for normal `openqp` command-line use.

## First Run

```bash
openqp examples/HF/H2O_RHF-HF_ENERGY.inp
```

Run the packaged example tests:

```bash
openqp --run_tests all
```

Control OpenMP threads per process or MPI rank with:

```bash
openqp input.inp --omp 16
```

or:

```ini
[input]
omp_threads=16
```

## Documentation

The user manual lives in the separate documentation repository:

- [OpenQP Manual](https://open-quantum-platform.github.io/openqp-docs/)
- [Manual source](https://github.com/Open-Quantum-Platform/openqp-docs)
- [Build options](https://open-quantum-platform.github.io/openqp-docs/build-options/)
- [Example inputs](examples)
- [API guide](https://open-quantum-platform.github.io/openqp-docs/api/)
- [Generated API reference](https://open-quantum-platform.github.io/openqp)

## Graphic Web Tools

- [OpenQP Web](https://app.openqp.org/) prepares inputs and previews structures
  locally in the browser.
- [OpenQP Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/)
  provides a browser input builder.
- [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) inspects
  OpenQP logs, JSON, Molden, cube, and XYZ data in the browser.

## Key Papers

If you use OpenQP in research, cite the OpenQP platform paper:

- Mironov V, Komarov K, Li J, Gerasimov I, Mazaheri M, Park W, Lashkaripour A,
  Oh M, Nakata H, Ishimura K, Huix-Rotllant M, Lee S, and Choi CH. "OpenQP: A
  Quantum Chemical Platform Featuring MRSF-TDDFT with an Emphasis on Open-source
  Ecosystem." Journal of Chemical Theory and Computation, 2024.
  https://doi.org/10.1021/acs.jctc.4c01117

Original MRSF-TDDFT theory and implementation papers:

- Lee S, Filatov M, Lee S, and Choi CH. "Eliminating Spin-Contamination of
  Spin-Flip Time-Dependent Density Functional Theory Within Linear Response
  Formalism by the Use of Zeroth-Order Mixed-Reference Reduced Density Matrix."
  Journal of Chemical Physics, 2018. https://doi.org/10.1063/1.5044202
- Lee S, Kim EE, Nakata H, Lee S, and Choi CH. "Efficient Implementations of
  Analytic Energy Gradient for Mixed-Reference Spin-Flip Time-Dependent Density
  Functional Theory (MRSF-TDDFT)." Journal of Chemical Physics, 2019.
  https://doi.org/10.1063/1.5086895

Recent MRSF-TDDFT accounts and overview papers:

- Park W, Komarov K, Lee S, and Choi CH. "Mixed-Reference Spin-Flip
  Time-Dependent Density Functional Theory: Multireference Advantages with the
  Practicality of Linear Response Theory." Journal of Physical Chemistry
  Letters, 2023. https://doi.org/10.1021/acs.jpclett.3c02296
- Lee S, Park W, and Choi CH. "Expanding Horizons in Quantum Chemical Studies:
  The Versatile Power of MRSF-TDDFT." Accounts of Chemical Research, 2025.
  https://doi.org/10.1021/acs.accounts.4c00640
- Park W, Lee S, Komarov K, Mironov V, Nakata H, Zeng T, Huix-Rotllant M, and
  Choi CH. "MRSF-TDDFT: A New Tool in Quantum Chemistry for Better
  Understanding Molecules and Materials." Bulletin of the Korean Chemical
  Society, 2025. https://doi.org/10.1002/bkcs.70011

## License

See [`LICENSE`](LICENSE).
