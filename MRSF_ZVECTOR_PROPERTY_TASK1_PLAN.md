# Task 1: Reusable Z-vector / Property-Response Interface Plan for MRSF-TDDFT

## Goal
Create a response infrastructure where new MRSF-TDDFT properties (NMR shielding, polarizability, hyperpolarizability, etc.) can be added by supplying property-specific operators/RHS while reusing the existing Z-vector solver and contraction pathways.

## Current code map (what exists now)

### MRSF-TDDFT energy path
- `source/modules/tdhf_mrsf_energy.F90` computes MRSF/UMRSF energies and excitation vectors, and stores intermediate arrays in tagarray records such as `OQP_td_bvec_mo`, `OQP_td_t`, and `OQP_td_energies` for downstream stages.
- This module is effectively the source of many excited-state quantities required by gradient/response code.

### MRSF-TDDFT gradient path
- `source/modules/tdhf_mrsf_gradient.F90` orchestrates one- and two-electron gradient terms.
- It consumes tagarray data (`OQP_DM_A/B`, `OQP_td_abxc`, `OQP_td_p`, `OQP_td_mrsf_density`) and dispatches MRSF-specific two-electron work through `mrsf_2e_grad`.
- The design is property-specific and tightly coupled to geometry derivatives.

### MRSF Z-vector solver
- `source/modules/tdhf_mrsf_z_vector.F90` already contains robust iterative linear-solver machinery (GMRES with restart/preconditioning and reusable workspace).
- This is the main reusable asset for future analytic property development.

### Existing property modules (ground-state style)
- `source/modules/electric_moments.F90` computes multipole moments from density and one-electron integrals.
- `source/modules/resp.F90` computes ESP/RESP charges from SCF density and electrostatic potentials.
- These modules show how to expose property drivers via C bindings and tagarray I/O, but they do not use excited-state Z-vector response.

## Proposed architecture

### 1) Introduce a response “contract” module
Add a new module (proposed file):
- `source/modules/tdhf_mrsf_response_api.F90`

Define types/interfaces to separate common solver flow from property-specific details:
- `type :: mrsf_response_context_t`
  - Holds references to `infos`, `basis`, MO blocks, state metadata, and cached intermediates.
- `type :: mrsf_property_rhs_t`
  - Stores property RHS vectors for alpha/beta/spin blocks and optional metadata (component index, frequency).
- Deferred procedures / callback interfaces:
  - `build_property_rhs(context, rhs)`
  - `apply_property_operator(context, zvec, contrib)`
  - `accumulate_property_value(context, rhs, zvec, value)`

Design rule: the Z-vector linear operator stays property-agnostic; only RHS and final contractions are property-specific.

### 2) Refactor the existing MRSF Z-vector module behind a stable driver
In `source/modules/tdhf_mrsf_z_vector.F90`, introduce a high-level entry point:
- `solve_mrsf_zvector(context, rhs, zvec, opts)`

Keep existing GMRES/preconditioner internals, but make the external interface independent from gradient naming assumptions. This avoids duplicating linear-solver code in every new property module.

### 3) Add a shared “response common” utility layer
Add a utility module (proposed file):
- `source/modules/tdhf_mrsf_response_common.F90`

Responsibilities:
- Loading/validating required tagarray records once.
- Building common intermediates (MO partition dimensions, reference density pieces, optional DFT objects).
- Providing deterministic memory ownership/cleanup helpers.

This removes repeated boilerplate currently spread across energy/gradient code paths.

### 4) Add the first property plugin as a thin vertical slice
Add a first response-style property module (proposed file):
- `source/modules/tdhf_mrsf_property_dipole.F90`

Scope:
- Static excited-state dipole correction through Z-vector using one-electron dipole integrals.
- Reuse `electric_moments`-like integral machinery (`multipole_integrals`) where possible.
- Output per-state/per-component contributions in the log and tagarray.

Reason to start with dipole plugin:
- Lower integral complexity than NMR.
- Strongly validates the response API and data flow before magnetic response.

### 5) Preserve current gradient behavior while migrating incrementally
- Keep `tdhf_mrsf_gradient.F90` behavior unchanged initially.
- Optionally (later) migrate parts of gradient Z-vector usage to call the same response driver, reducing duplicate implementations and ensuring one solver path.

## Recommended implementation sequence (small PR slices)

### PR-A (infrastructure only)
1. Add `tdhf_mrsf_response_api.F90` with core types/interfaces.
2. Add `tdhf_mrsf_response_common.F90` with context assembly and tag checks.
3. Add non-invasive wrapper entry point in `tdhf_mrsf_z_vector.F90`.
4. No scientific behavior change expected.

### PR-B (first property)
1. Add `tdhf_mrsf_property_dipole.F90` implementing RHS and contraction callbacks.
2. Add C-bind wrapper and runtime dispatch hook.
3. Add at least one small MRSF example input/output check for regression.

### PR-C (extend to NMR)
1. Implement magnetic perturbation integral path and gauge strategy.
2. Reuse response API + solver from PR-A.
3. Add targeted reference tests.

## Data/API conventions to lock early
- **State indexing**: define a single convention for ground/excited state IDs between tagarray, logs, and wrappers.
- **Spin-block ordering**: codify alpha/beta and occupied/virtual flattening in one place.
- **Tag naming**: reserve new tags with property namespace, e.g. `OQP_td_prop_<name>`.
- **Numerical controls**: centralize GMRES tolerance/restart/maxiter settings and expose in input schema.

## Performance opportunities aligned with this refactor
- Reuse GMRES workspace already present in `tdhf_mrsf_z_vector.F90` across property components/states.
- Batch property components (e.g., X/Y/Z) when applying common kernels.
- Cache transformed one-electron operators in MO basis for all states/components.
- Ensure DFT response objects are initialized once per job where feasible.

## Validation plan for Task 1 outputs
For the infrastructure PR (PR-A), verify:
1. Existing MRSF energy + gradient examples still reproduce baseline values.
2. New response API wrappers compile without changing existing execution paths.
3. Logging prints dimensions and solver options for traceability.

## Definition of done for Task 1
Task 1 is complete when:
- A documented response architecture exists (this plan).
- A concrete PR-A file-level change list is agreed.
- Interfaces are defined to let future properties be implemented without touching GMRES internals.

## Proposed file-level change list for PR-A
- **Modify** `source/modules/tdhf_mrsf_z_vector.F90`
  - Add stable high-level `solve_mrsf_zvector` driver and option struct.
- **Add** `source/modules/tdhf_mrsf_response_api.F90`
  - Property callback interfaces and core response types.
- **Add** `source/modules/tdhf_mrsf_response_common.F90`
  - Context initialization, tag validation, and shared helpers.
- **(Optional) Modify** `source/modules/tdhf_mrsf_gradient.F90`
  - Only minimal wiring if needed; otherwise defer to later migration PR.

## Coordination update (Apr 24, 2026): PR #124 status and integration guidance

A separate PR is already in progress:
- **PR #124**: `feat: add excited-state electric moments`.
- Scope observed from PR metadata: new C API symbol (`electric_moments_excited`), Python call-site wiring in `single_point.py`, and a large addition to `electric_moments.F90` implementing excited-state multipole analysis from relaxed density (`OQP_TD_P`).

### What this means for our sequence
- Treat PR #124 as the functional equivalent of the originally proposed PR-B dipole vertical slice.
- Keep this Task-1 infrastructure work focused on PR-A (response API + common context + stable Z-vector driver), and **do not** open a second competing dipole implementation.

### PR #124 quick review checklist (to merge safely)
1. **State correctness check**
   - Confirm `target_state` is used to select the correct per-state relaxed density contribution, not only printed.
2. **Data contract check**
   - Ensure tag naming is consistent (`OQP_td_p` vs `OQP_TD_P`) and documented in one place.
3. **Code-duplication check**
   - Prefer extracting shared multipole contraction/output helpers instead of duplicating the whole ground-state routine.
4. **User-facing behavior check**
   - Gate excited-state moment printing behind an explicit input flag or clear workflow stage to avoid surprising output during gradients.
5. **Requested follow-up from review**
   - A maintainer comment requested adding **Mulliken charge** support; define whether this belongs in the same PR or an immediate follow-up PR.

### Updated next actionable order
1. Finalize/review PR #124 (electric moments for excited states).
2. Implement PR-A infrastructure (response API/common/driver).
3. Rebase/refactor PR #124 onto new infrastructure only if the churn is small; otherwise merge first and migrate in a dedicated cleanup PR.
4. Start NMR shielding on top of the stabilized response driver.
