# PCM validation matrix

This matrix defines the minimum evidence required before enabling energy-only PCM runtime support in OpenQP. The first scientific target is **MRSF-TDDFT with PCM-solvated ROHF reference**: self-consistent PCM on the high-spin ROHF reference, followed by the existing gas-phase MRSF response machinery on those solvated reference orbitals.

## Scope guardrails

In scope for the first runtime milestone:

- `runtype=energy` only.
- RHF closed-shell reference SCF and ROHF high-spin reference SCF.
- MRSF-TDDFT only through a PCM-solvated ROHF reference.
- PySCF reference calculations for independent energy/sign convention checks.
- ddX as the preferred optional native backend candidate; PCMSolver remains a fallback candidate.

Out of scope until separate implementation and validation:

- analytic PCM gradients and any PCM geometry optimization.
- state-specific excited-state PCM.
- nonequilibrium or linear-response solvent in the MRSF kernel.
- solvent response terms in MRSF sigma vectors, Z-vectors, or transition properties.

## Evidence levels

| Level | Evidence | Required before claiming |
| --- | --- | --- |
| 1 | Input/schema guardrails | Reserved `[pcm]` options parse, invalid backend/model/mode/scope combinations fail clearly. |
| 2 | Backend smoke | Optional ddX/PCMSolver configuration, link, lifecycle, and adapter smoke tests pass without enabling runtime SCF coupling. |
| 3 | SCF seam smoke | OpenQP can build unweighted MEP data and an AO reaction-field matrix from backend cavity charges. |
| 4 | Energy validation | RHF/ROHF `runtype=energy` with PCM matches PySCF/reference backend energies within documented tolerances. |
| 5 | MRSF reference validation | MRSF response runs on PCM-solvated ROHF orbitals, with gas-phase vs solvated reference energies and response roots reported separately. |

## Initial calculation set

| Case | Method target | Backend/reference | Purpose |
| --- | --- | --- | --- |
| H2O | RHF/6-31g* water PCM, `runtype=energy` | PySCF + ddX candidate | Small closed-shell sign/unit check. |
| CH2O | BHHLYP/6-31g* water PCM, `runtype=energy` | PySCF + ddX candidate | DFT closed-shell energy and dipole-sensitive case. |
| small triplet/ROHF molecule | ROHF/6-31g* water PCM, `runtype=energy` | PySCF/reference package if available | High-spin reference-SCF validation before MRSF. |
| tiny MRSF case | MRSF-TDDFT/BHHLYP/6-31g* on PCM-solvated ROHF reference | OpenQP runtime vs gas-phase control | First MRSF scope check; not state-specific PCM. |

## ddX reaction-field caveat

The current adapter smoke exposes cavity coordinates and a projected cavity quantity suitable for the future `external_charge_potential` seam. The **ddX q_cav sign/scale is provisional**: finite-difference smoke data suggests testing `chg = -0.5*q_cav`, but this must be cross-checked against PySCF/ddX/reference-package data before runtime PCM is enabled.

## Incremental Fock guard

Reviewed reference-PCM reaction potentials may only enter the first prototype through the non-incremental SCF Fock path. If OpenQP old-buffer state is present, the guard must fail fast with the diagnostic `reference PCM incremental Fock is not validated` and preserve whether `dens_old`, `f_old`, or both triggered the incremental Fock shortcut. This keeps the candidate PCM energy bookkeeping from silently mixing with an unvalidated incremental-Fock reuse path.

## Reporting rules

For validation reports, include:

- gas-phase and solvated total energies in Hartree;
- polarization energy if exposed by the backend;
- backend, model, dielectric, cavity/radii settings;
- SCF iteration count and convergence status;
- whether TRAH or another fallback occurred;
- exact code commit, backend version, and artifact path.

Do not report solution-phase excitation energies as validated until the PCM-solvated ROHF reference path has live OpenQP evidence and the MRSF response-root mapping is checked separately from the reference SCF energy shift.
