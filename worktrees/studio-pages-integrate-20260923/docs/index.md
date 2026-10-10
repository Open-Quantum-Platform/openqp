# OpenQP Manual

**Manual version:** 1.3.0; **OpenQP version:** 1.3.0

Open Quantum Platform (OpenQP) is a quantum chemistry package centered on
HF/DFT, MP2, TDHF/TDDFT, SF-TDDFT, MRSF-TDDFT, and related workflows for
multiconfigurational ground and excited states.
The manual is organized around what users usually need first: install OpenQP,
prepare an input, run a calculation, then look up keywords when a workflow needs
more control. The lower-level API reference is intentionally placed last, after
workflows and keywords.

## Start Here

- [Installation](installation.md): package install, source build, BLAS/LAPACK,
  OpenMP, MPI, and runtime-file layout.
- [Build Options](build-options.md): complete CMake options, package-build
  overrides, BLAS/LAPACK choices, and external dependency cache behavior.
- [Quickstart](quickstart.md): the shortest path from a molecule to an OpenQP
  output file.
- [OQP Studio](studio/index.md): install and use the desktop interface for
  molecular building, execution, analysis, and ray-traced rendering.
- [`.oqp` input](oqp-input.md): the recommended input style used first
  throughout the workflow chapters.
- [Run OpenQP from Python](python-scripting.md): the second, programmatic style
  for scripts, notebooks, and workflow managers.
- [Legacy `.inp` input](input-file.md): the supported sectioned format for
  existing decks and legacy-only controls.
- [Examples](examples/index.md): runnable inputs stored in the OpenQP code
  repository.
- [References](references.md): platform, MP2, MRSF-TDDFT, SOC, scalar
  relativistic, and PCM/ddX papers
  cited by the manual.

## Common Workflows

- Ground-state [HF and DFT](workflows/hf-dft.md)
- [MP2](workflows/mp2.md) and spin-scaled MP2 energies, plus analytic RHF gradients
- [SC-NEVPT2](workflows/sc-nevpt2-gradient.md) analytic gradients and geometry calculations
- [TDDFT and TDHF](workflows/tddft.md) energies and gradients
- [SF-TDDFT](workflows/sf-tddft.md) spin-flip energies and gradients
- [MRSF-TDDFT](workflows/mrsf-tddft.md) energies and gradients
- [QMRSF-DK](workflows/qmrsf-dk.md) quintet-reference CAS(4,4) spectra
- [Hessian and frequencies](workflows/hessian.md)
- Geometry [optimization](workflows/optimization.md)
- [Spin-orbit coupling](workflows/soc.md) and scalar relativistic DKH correction
- [NACME](workflows/nacme.md)
- Energy-only [PCM/ddX](workflows/pcm.md)
- [MRSF-EKT](workflows/ekt.md) ionization potentials and electron affinities
- [NMR, IR, and Raman](workflows/nmr-ir-raman.md)

## Keyword Reference

The [keyword reference](keywords/index.md) is the code-aligned lookup layer for
high-drift input sections such as `[input]`, `[scf]`, `[optimize]`, `[oqp]`,
`[mp2]`, `[pcm]`, and `[symmetry]`. It should be checked against
[`pyoqp/oqp/molecule/oqpdata.py`](https://github.com/Open-Quantum-Platform/openqp/blob/main/pyoqp/oqp/molecule/oqpdata.py)
and
[`pyoqp/oqp/utils/input_checker.py`](https://github.com/Open-Quantum-Platform/openqp/blob/main/pyoqp/oqp/utils/input_checker.py)
when OpenQP changes.

## API Documentation

Workflow pages use the same order throughout the manual: recommended `.oqp`,
Python, then legacy `.inp`. The [Run OpenQP from Python](python-scripting.md)
chapter expands the script-based form. The
[API chapter](api/index.md) appears last in the manual and is the lower-level
reference for `oqp.pyoqp.Runner`, in-memory inputs, result extraction through
`runner.mol`, and the input-checking API used by front ends and automated
workflows.

## Web Tools

- [OpenQP Web](https://app.openqp.org/) prepares inputs in the browser and
  previews structures locally.
- [OpenQP Input Generator](https://open-quantum-platform.github.io/OpenQP_Input_Generator/)
  provides a browser-based input builder.
- [OpenqpView](https://open-quantum-platform.github.io/OpenqpView/) inspects
  OpenQP logs, JSON, Molden, cube, and XYZ data in the browser.
