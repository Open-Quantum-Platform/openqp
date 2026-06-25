# OpenQP

Open Quantum Platform (OpenQP) is an open-source quantum chemistry package built
around HF/DFT, TDHF/TDDFT, SF-TDDFT, MRSF-TDDFT, and related excited-state
workflows.

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
files together. Normal `openqp` command-line use does not require
`OPENQP_ROOT`.

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

The repo-owned manual starts in [`docs/index.md`](docs/index.md).

Important entry points:

- [Installation](docs/installation.md)
- [Quickstart](docs/quickstart.md)
- [Input files](docs/input-file.md)
- [Capabilities](docs/capabilities.md)
- [Workflow examples](docs/examples/index.md)
- [Keyword reference](docs/keywords/index.md)
- [Troubleshooting](docs/troubleshooting.md)

To preview the manual locally:

```bash
pip install -r docs/requirements.txt
mkdocs serve
```

## Highlights

- RHF, ROHF, and UHF references for HF and DFT
- TDHF/TDDFT, SF-TDDFT, MRSF-TDDFT, and energy-only UMRSF-TDDFT
- MRSF-EKT ionization potentials and electron affinities
- Analytic gradients and native HF/DFT Hessians for supported workflows
- MRSF-TDDFT NACME and spin-orbit couplings
- Energy-only reference-SCF PCM/ddX
- Native geometry optimizer (`lib=oqp`) plus geomeTRIC and SciPy backends
- NMR shielding, vibrational frequencies, IR, and Raman intensities
- OpenMP/MPI parallelism, ILP64 BLAS/LAPACK, pip installs, and Docker support

## Web Tools

- [OpenQP Web](https://app.openqp.org/) prepares inputs and previews structures
  locally in the browser.
- [OpenQP Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/)
  provides a browser input builder.
- [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) inspects
  OpenQP logs, JSON, Molden, cube, and XYZ data in the browser.

## Citation

If you use OpenQP in research, cite:

Mironov V, Komarov K, Li J, Gerasimov I, Mazaheri M, Park W, Lashkaripour A,
Oh M, Nakata H, Ishimura K, Huix-Rotllant M, Lee S, and Choi CH. "OpenQP: A
Quantum Chemical Platform Featuring MRSF-TDDFT with an Emphasis on Open-source
Ecosystem." Journal of Chemical Theory and Computation, 2024.
https://doi.org/10.1021/acs.jctc.4c01117

## License

See [`LICENSE`](LICENSE).
