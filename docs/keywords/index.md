# Keyword Reference

This section documents the high-drift input sections that should stay aligned
with:

- `pyoqp/oqp/molecule/oqpdata.py`
- `pyoqp/oqp/utils/input_checker.py`

The keyword pages are intentionally compact. Workflow pages show complete input
decks; keyword pages show defaults, allowed values, and cross-section behavior.

## Sections

| Section | Purpose |
| --- | --- |
| `[input]` | Global calculation setup, geometry, run type, AO convention, threading. |
| `[guess]` | Initial orbitals, restart files, and MO swaps. |
| `[scf]` | Reference type, convergence controls, pFON, TRAH, fallback manager. |
| `[dftgrid]` | DFT quadrature and hybrid/range-separated functional controls. |
| `[tdhf]` | TDHF/TDDFT/SF/MRSF/UMRSF response settings. |
| `[properties]` | Gradients, NMR, export, and property requests. |
| `[hess]` | Hessian, frequency, and thermochemistry controls. |
| `[nac]` | NAC and NACME state-pair controls. |
| `[ekt]` | MRSF-EKT IP/EA channel selection. |
| `[optimize]` | Target state and optimization convergence controls. |
| `[oqp]` | Native optimizer controls. |
| `[geometric]` | geomeTRIC optimizer backend controls. |
| `[neb]` | NEB product endpoint and image count. |
| `[pcm]` | Reference-SCF PCM/ddX settings. |
| `[symmetry]` | Point-group metadata and optional symmetry reductions. |
