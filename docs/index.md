# OpenQP Manual

Open Quantum Platform (OpenQP) is a quantum chemistry package centered on
HF/DFT, TDHF/TDDFT, SF-TDDFT, MRSF-TDDFT, and related excited-state workflows.
The manual is organized around what users usually need first: install OpenQP,
prepare an input, run a calculation, then look up keywords when a workflow needs
more control.

## Start Here

- [Installation](installation.md): package install, source build, BLAS/LAPACK,
  OpenMP, MPI, and runtime-file layout.
- [Quickstart](quickstart.md): the shortest path from a molecule to an OpenQP
  output file.
- [Input Files](input-file.md): section layout, geometry input, run types, and
  output files.
- [Examples](examples/index.md): runnable inputs already stored in the repository.

## Common Workflows

- Ground-state [HF and DFT](workflows/hf-dft.md)
- [MRSF-TDDFT](workflows/mrsf-tddft.md) energies and gradients
- Geometry [optimization](workflows/optimization.md)
- [Spin-orbit coupling](workflows/soc.md)
- [NACME](workflows/nacme.md)
- Energy-only [PCM/ddX](workflows/pcm.md)
- [MRSF-EKT](workflows/ekt.md) ionization potentials and electron affinities
- [NMR, IR, and Raman](workflows/nmr-ir-raman.md)

## Keyword Reference

The [keyword reference](keywords/index.md) is the code-aligned lookup layer for
high-drift input sections such as `[input]`, `[scf]`, `[optimize]`, `[oqp]`,
`[pcm]`, and `[symmetry]`. It should be updated together with
`pyoqp/oqp/molecule/oqpdata.py` and `pyoqp/oqp/utils/input_checker.py`.

## Web Tools

- [OpenQP Web](https://app.openqp.org/) prepares inputs in the browser and
  previews structures locally.
- [OpenQP Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/)
  provides a browser-based input builder.
- [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) inspects
  OpenQP logs, JSON, Molden, cube, and XYZ data in the browser.
