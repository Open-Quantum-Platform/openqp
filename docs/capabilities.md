# Capabilities

This page summarizes what the current user manual documents. It is not a
promise that every combination of method, property, and backend is available.
When a workflow has limits, use the linked workflow page or keyword page for the
specific input contract.

## Electronic Structure

| Area | Status |
| --- | --- |
| HF and DFT | RHF, ROHF, and UHF references. |
| TDHF/TDDFT | Energy and gradient workflows for supported references. |
| SF-TDDFT and MRSF-TDDFT | Energy, gradients, NACME, SOC, optimization workflows. |
| UMRSF-TDDFT | Energy-only UHF-reference workflow. |
| MRSF-EKT | IP/EA analysis with Dyson-like orbital data. |

## Properties

| Property | Status |
| --- | --- |
| Analytic gradients | Available for the supported HF/DFT and response workflows. |
| HF/DFT Hessians | Native analytic path for supported HF/DFT references. |
| Numerical Hessians | Available through the Hessian workflow. |
| NACME | MRSF-TDDFT state-coupling workflow. |
| SOC | MRSF-TDDFT one-electron and mean-field two-electron SOC. |
| PCM/ddX | Energy-only reference-SCF path for RHF/ROHF. |
| NMR | Nuclear magnetic shielding via `[properties] scf_prop=nmr`. |
| IR/Raman | Frequency-analysis intensities from supported Hessian workflows. |

## Geometry and Paths

The native optimizer is selected with `[optimize] lib=oqp` and supports
`optimize`, `ts`, `meci`, `mecp`, `tci`, `neb`, `irc`, and `mep`. geomeTRIC and
SciPy remain optional backends for the workflows wired to them.

## Upcoming or Limited Areas

- Production electrostatic embedding QM/MM is an active development direction.
- PCM gradients, PCM optimizations, and state-specific excited-state PCM are not
  part of the first ddX energy path.
- UMRSF-TDDFT gradients and Hessians are not part of the documented production
  surface yet.
