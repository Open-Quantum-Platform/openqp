# Capabilities

This page summarizes what the current user manual documents. It is not a
promise that every combination of method, property, and backend is available.
When a workflow has limits, use the linked workflow page or keyword page for the
specific input contract.

## Electronic Structure

| Area | Status |
| --- | --- |
| HF and DFT | RHF, ROHF, and UHF references. |
| MP2 | RHF, UHF, and ROHF energies with spin-scaled variants; analytic RHF gradients and RHF gradient-driven geometry calculations. |
| Coupled cluster | Energy-only CCSD and CCSD(T) on RHF, UHF, and ROHF references, with a frozen core. In-core integrals; the open-shell path is a spin-orbital solver for small systems. See [Coupled Cluster](workflows/coupled-cluster.md). |
| Wavefunction methods | Native determinant-space FCI, CASCI, CASSCF, state-averaged CASSCF, CASPT2, NEVPT2, and QDPT-family energies on an RHF reference. State-specific CASSCF, dedicated SA-CASSCF, single-state CASPT2/MRMP2, MCQDPT2, XMS-CASPT2/XMCQDPT2, and single-state SC-NEVPT2 have analytic nuclear gradients; the SA-CASSCF compatibility spelling, multi-set MS-CASPT2, and partially contracted NEVPT2 use central differences. See [Wavefunction methods](keywords/wavefunction.md). |
| CASSCF nuclear gradient | Analytic gradients for state-specific `casscf`, the weighted SA-CASSCF objective, and an individual averaged root through a coupled orbital and CI Z-vector. Individual-state gradients support gradient-driven optimization; the weighted objective currently supports direct gradients only. CASCI and FCI remain energy-only. See [CASSCF Nuclear Gradient](workflows/casscf-gradient.md). |
| CASPT2 nuclear gradient | Analytic gradients for single-state CASPT2/MRMP2, MCQDPT2, and XMS-CASPT2/XMCQDPT2 (`[pt2] gradient=auto` prefers them and falls back to central differences when a precondition fails). See [CASPT2 Nuclear Gradient](workflows/caspt2-gradient.md). |
| SC-NEVPT2 nuclear gradient | Analytic gradient and gradient-driven geometry calculations for single-state, strongly contracted Dyall-Hamiltonian NEVPT2 on a state-specific CASSCF reference. See [SC-NEVPT2 Nuclear Gradient](workflows/sc-nevpt2-gradient.md). |
| TDHF/TDDFT | Energy and gradient workflows for supported references. |
| SF-TDDFT and MRSF-TDDFT | Multiconfigurational ground- and excited-state energies, gradients, NACME, SOC, and optimization workflows. |
| UMRSF-TDDFT | Energy-only UHF-reference workflow. |
| MRSF-EKT | IP/EA analysis with Dyson-like orbital data. |

## Properties

| Property | Status |
| --- | --- |
| Analytic gradients | Available for the supported HF/DFT, RHF MP2 and spin-scaled MP2, state-specific CASSCF, dedicated SA-CASSCF weighted and individual-root derivatives, the analytic PT2 variants, in-scope SC-NEVPT2, and response methods. |
| Numerical multireference gradients | Central differences for `method=casscf` with state averaging enabled, multi-set MS-CASPT2, and NEVPT2 outside the analytic SC-NEVPT2 scope, with `grad`, `optimize`, `ts`, `mep`, and `irc`. |
| HF/DFT Hessians | Native analytic path for supported HF/DFT references. |
| Numerical Hessians | Available through the Hessian workflow. |
| NACME | MRSF-TDDFT state-coupling workflow. |
| SOC | MRSF-TDDFT one-electron and mean-field two-electron SOC. |
| MRSF excited-state analysis | NTOs, attachment/detachment densities, state-to-state transition densities, cube export, QCSchema export, FCIDUMP export, and external-code comparisons through `oqp.interop`. |
| Nonadiabatic MD (NAMD) | Tully surface hopping, SOC-NAMD, and ESPF QM/MM embedding. See [SOC-NAMD-QMMM](workflows/soc-namd-qmmm.md). |
| QM/MM | ESPF electrostatic embedding for single-point energies, ground-state MD, and nonadiabatic dynamics. See [`[qmmm]`](keywords/qmmm.md). |
| Scalar relativistic correction | Spin-free DKH correction through `[scf] scal_rel=1` or `2`. |
| PCM/ddX | Energy-only reference-SCF path for RHF/ROHF. |
| NMR | Nuclear magnetic shielding via `[properties] scf_prop=nmr`. |
| IR/Raman | Frequency-analysis intensities from supported Hessian workflows. |

## Geometry and Paths

The native optimizer supports `optimize`, `ts`, `meci`, `mecp`, `neb`, `irc`,
and `mep`; legacy `tci` inputs remain compatible. MECI includes the BaekA
adaptive algorithm for two or more same-manifold states. Concise `.oqp`
geometry drivers select the native engine automatically and do not expose a
backend selector. Native TS supports model, numerical, or analytical initial
Hessians and mode following. Native NEB supports endpoint alignment and
relaxation, climbing images, maximum and RMS force convergence, and final
multi-frame XYZ path output.

Method availability is narrower than the optimizer's complete run-type list.
The native multireference wavefunction methods support `optimize`, `ts`,
`mep`, and `irc`: state-specific CASSCF, an individual root of dedicated
SA-CASSCF, the analytic PT2 variants (single-state CASPT2/MRMP2, MCQDPT2,
XMS-CASPT2/XMCQDPT2), and in-scope SC-NEVPT2 use analytic gradients, while the
SA-CASSCF compatibility spelling, multi-set MS-CASPT2, and partially contracted
NEVPT2 use numerical gradients. The weighted SA objective is available for
direct gradients only because it has no state-energy entry for the optimizer.
These methods do not support `meci`, `mecp`, or `neb`; FCI and
CASCI remain energy-only.

geomeTRIC and SciPy remain optional compatibility backends for traditional
sectioned `.inp` files and the Python API. Native `.oqp` supports frozen-distance
minimum searches; legacy geomeTRIC remains an escape hatch for advanced
constraint types that the concise grammar does not yet support.

## Upcoming or Limited Areas

- XTB, DFTB, and AFQMC backends are not distributed or supported as OpenQP
  1.3.0 capabilities. Their separate implementation repositories and earlier
  documentation proposals are outside this release's public support surface.
- Electrostatic embedding QM/MM is an active development direction. Nonadiabatic
  QM/MM dynamics currently supports whole-molecule QM regions only; covalent
  QM/MM boundaries (link atoms) in dynamics are not yet available.
- PCM gradients, PCM optimizations, and state-specific excited-state PCM are not
  part of the first ddX energy path.
- Open-shell MP2 gradients, MP2 Hessians, RI/Laplace/local MP2 kernels, and
  periodic MP2 are not part of the documented standalone MP2 path.
- Coupled-cluster gradients and derivative workflows are not implemented, and
  there is no density-fitted or disk-based coupled-cluster mode; the CC
  iteration itself is not integral-direct and its integrals are held in memory.
  The closed-shell Cholesky path can still *build* its vectors directly from
  recomputed AO integrals, skipping the packed AO store — see
  [`cholesky_direct`](keywords/cc.md#cholesky_direct).
- UMRSF-TDDFT gradients and Hessians are not part of the documented production
  surface yet.
