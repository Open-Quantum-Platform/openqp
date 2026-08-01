# OpenQP WaveFunction examples

This directory contains native OpenQP wavefunction examples for the FCI/CAS/CASPT2 project branch.

Current native support snapshot:

| Method | Native status | Main limitations |
| --- | --- | --- |
| FCI | determinant CI energy | energy-only; spin-adapted CSF solver is not implemented |
| CASCI | fixed-orbital active-space CI energy | energy-only; orbital localization is not implemented |
| CASSCF | closed-shell RHF macroiterations | energy-only; open-shell and analytic nuclear gradients are not implemented |
| SA-CASSCF | closed-shell RHF state-averaged macroiterations | energy-only; root tracking is guarded by overlap diagnostics |
| CASPT2 | uncontracted determinant-space state-specific PT2 | internally contracted and noncanonical PT2 are not implemented |
| MS-CASPT2 | uncontracted determinant-space multistate PT2 | internally contracted and noncanonical PT2 are not implemented |
| XMS-CASPT2 | XMS state rotation plus uncontracted determinant-space multistate PT2 | internally contracted and noncanonical PT2 are not implemented |

Retained runnable examples:

- `H2_RHF-CASCI_ENERGY.inp` — H2/STO-3G CASCI using the native OpenQP determinant-CI engine.
- `H2_RHF-CASCI_CASPT2_REFERENCE_PREFLIGHT.inp` — H2/STO-3G CASCI that writes a diagnostic-only CASPT2 reference preflight JSON sidecar from native OpenQP CASCI tensors. This is not a native PT2 energy calculation.
- `H2_RHF-CASCI_DAVIDSON.inp` — H2/STO-3G CAS(2e,2o) CASCI with `nroot=1` and `solver=davidson`; this is example/interface coverage only, not a new method acceptance or performance claim.
- `H2_RHF-CASCI_NROOT2.inp` — H2/STO-3G CAS(2e,2o) CASCI with `nroot=2`, reporting the singlet ground root and Ms=0 triplet component.
- `H2O_RHF-CASCI_CAS2E2O_FROZEN_CORE.inp` — H2O/STO-3G frozen-core CASCI-style active-space calculation.
- `H2O_RHF-CASCI_CAS2E2O_MO_INDICES.inp` — H2O/STO-3G CAS(2e,2o) CASCI using explicit disjoint 1-based core and active MO labels.
- `H2_RHF-CASCI_JSON_MO_INDICES.inp` — CASCI with OpenQP JSON MO coefficients and explicit 1-based active-orbital labels.
- `H2_RHF-CASSCF_ZERO_MACRO_DIAGNOSTIC.inp` — CASSCF zero-orbital-update scaffold example that writes a diagnostic-only native CASSCF microiteration JSON sidecar and compares the diagnostic electronic energy to a PySCF CASCI reference. This is not an accepted orbital-optimized CASSCF calculation.
- `LiH_RHF-CASSCF_ZERO_MACRO_DIAGNOSTIC_NONZERO_GRADIENT.inp` — CASSCF zero-orbital-update scaffold example for LiH/STO-3G frozen-core CAS(2e,2o), retaining the PySCF CASCI benchmark gate while recording a nonzero native orbital-gradient diagnostic trial. This is not an accepted orbital-optimized CASSCF calculation.
- `LiH_RHF-CASSCF_ONE_MACROITERATION.inp` — Native OpenQP state-specific CASSCF smoke example for LiH/STO-3G frozen-core CAS(2e,2o). It accepts one orbital macroiteration, commits the live MO update, recomputes CASSCF CI/RDM tensors, and writes a strict macroiteration evidence JSON sidecar.
- `LiH_RHF-SA-CASSCF_ONE_MACROITERATION.inp` — Native OpenQP equal-weight two-state SA-CASSCF smoke example for LiH/STO-3G frozen-core CAS(2e,2o). It runs one weighted state-average macroiteration and writes a strict macroiteration evidence JSON sidecar. The paired JSON/log files store a fully converged PySCF SA-CASSCF reference for the same molecule/basis/active space.
- `LiH_RHF-CASSCF_POWELL_CONVERGED.inp` — Native OpenQP state-specific CASSCF example for LiH/STO-3G inactive-core CAS(2e,2o), converged with `optimizer=powell` using OpenQP-native FCI/RDM objective evaluations and benchmarked against PySCF.
- `LiH_RHF-SA-CASSCF_POWELL_CONVERGED.inp` — Native OpenQP equal-weight two-state SA-CASSCF example for LiH/STO-3G inactive-core CAS(2e,2o), converged with `optimizer=powell` and benchmarked against the matching PySCF inactive-core state-average reference.
- `H4_RHF-CASCI_CASPT2_NATIVE.inp` — Native OpenQP state-specific uncontracted determinant-space CASPT2 for linear H4/STO-3G with a CASCI reference, IPEA shift, and amplitude save/print controls.
- `H4_RHF-CASCI_MS-CASPT2_NATIVE.inp` — Native OpenQP uncontracted determinant-space MS-CASPT2 for selected CASCI roots of linear H4/STO-3G.
- `H4_RHF-CASCI_XMS-CASPT2_NATIVE.inp` — Native OpenQP XMS state rotation plus uncontracted determinant-space multistate PT2 for selected CASCI roots of linear H4/STO-3G.

Removed by scientific-scope correction:

- NEVPT2 examples were removed because the previous implementation used PySCF MRPT rather than native OpenQP perturbation theory.
- Molden-import CASCI/CASSCF examples were removed where the import path required PySCF's Molden parser.

External programs may still be used in tests as references, but no user-facing method/example in this directory should require a non-native backend or be presented as an OpenQP implementation unless the corresponding native OpenQP algorithm has been implemented and validated.
