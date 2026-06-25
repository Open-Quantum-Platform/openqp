# Examples

The repository examples are the preferred source of runnable input decks. Use
them as templates before writing a new input from scratch.

| Folder | Contents |
| --- | --- |
| `examples/HF` | RHF, ROHF, and UHF Hartree-Fock energies and gradients. |
| `examples/DFT` | DFT energies, gradients, and optimization examples. |
| `examples/SCF` | SCF convergence controls such as DIIS variants, MOM, pFON, and SOSCF. |
| `examples/TDHF` | TDHF energy and gradient examples. |
| `examples/TDDFT` | TDDFT energy and gradient examples. |
| `examples/SF-TDDFT` | Spin-flip TDDFT examples. |
| `examples/MRSF-TDDFT` | MRSF-TDDFT energies, gradients, and optimization data. |
| `examples/OPT` | Native and geomeTRIC optimization, MECI, MECP, TCI, TS, IRC, NEB, and MEP. |
| `examples/HESS` | Analytic and numerical Hessian workflows. |
| `examples/PCM` | ddX reference-SCF PCM energy cases. |
| `examples/SOC` | MRSF-TDDFT SOC cases. |
| `examples/NMR` | NMR shielding examples. |
| `examples/TRAH` | TRAH SCF examples. |
| `examples/ISPHER` | Spherical-harmonic AO convention examples. |
| `examples/ECP` | Effective-core-potential examples. |
| `examples/UMRSF-TDDFT` | UMRSF-TDDFT energy examples. |
| `examples/XAS` | X-ray absorption examples. |

Run a single example:

```bash
openqp examples/HF/H2O_RHF-HF_ENERGY.inp
```

Run the packaged example tests:

```bash
openqp --run_tests all
```

When adding a new manual page, link to a repository example whenever possible
instead of pasting a long input deck into prose.
