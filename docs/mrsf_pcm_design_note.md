# MRSF-PCM design note

Status: design gate only. Do not implement MRSF-PCM until HF/DFT PCM diagnostics
remain stable and an identical-protocol ground-state PCM reference strategy is in
place.

## Why MRSF-PCM is a separate design problem

The current PCM path is an equilibrium ground-state SCF coupling. MRSF-TDDFT adds
state-specific and response questions that cannot be answered by simply calling
`add_pcm_reaction_field` in the existing SCF loop.

Key choices must be explicit before coding:

1. Which density sources the reaction field?
   - triplet reference density,
   - state-averaged MRSF density,
   - state-specific relaxed density,
   - transition density or difference density contribution.

2. Which solvation regime is modeled?
   - equilibrium state-specific solvation,
   - nonequilibrium vertical excitation solvation,
   - frozen ground-state reaction field plus fast electronic response.

3. Where does the reaction field enter?
   - reference SCF/orbital optimization only,
   - MRSF Hamiltonian diagonal terms,
   - off-diagonal state couplings,
   - response equations / transition-density terms.

4. What energy is reported?
   - ddX scalar for the source density,
   - state-resolved correction,
   - difference between state-specific PCM energies,
   - nonequilibrium correction split into slow/fast components.

## Proposed staged implementation

### Stage 0: guardrail and status ledger

Add an explicit runtime log notice for any MRSF calculation with `[pcm] enabled`
that says MRSF-PCM is design-gated/deferred unless a specific stage below is
implemented. Add a regression test so manuscript/report text cannot accidentally
claim MRSF-PCM is active.

### Stage 1: frozen-reference reaction field only

Use the converged reference SCF density to build the same ddX source as the HF/DFT
path. Insert the reaction field only into the reference SCF/Fock path. Do not
modify the MRSF state Hamiltonian or report state-specific solvation energies.

Validation:

- PCM-off MRSF energies unchanged.
- PCM-on log states `frozen-reference PCM only`.
- Ground-state/reference density diagnostics match the HF/DFT PCM invariants:
  `phi_source_vs_exact`, `q_cav_sum`, and fixed-density `-0.5*q_cav` derivative.

This is useful only as an orbital-environment experiment, not a physical
MRSF-PCM model.

### Stage 2: state-diagonal density correction prototype

For each requested MRSF state, build a state-specific one-particle density (or
well-defined approximation) and evaluate a ddX scalar correction with the same
source convention. Keep it post-SCF/post-diagonalization first; do not feed back
into orbital or CI response.

Validation:

- exact density/ESP diagnostic per state,
- invariant source-vs-exact cavity potential check,
- symmetry and state-index stability checks,
- no claims about gradients or nonequilibrium response.

### Stage 3: nonequilibrium/state-response formulation

Only after Stage 2 is numerically stable, design the slow/fast solvent split and
transition-density/state-coupling terms. This stage requires theory review before
implementation.

## Initial tests to add before coding physics

- Parser/config test that `[pcm] enabled` with an MRSF runtype produces an
  explicit deferred/staged status line.
- PCM-off MRSF regression remains unchanged.
- Stage 1, if enabled, must print `frozen-reference` and must not print
  state-specific PCM energies.
- Any future state-specific mode must print the density/source used for each
  state and must emit `phi_source_vs_exact` diagnostics.

## Non-goals for the next increment

- gradients,
- TDDFT/MRSF response solvent terms,
- nonequilibrium PCM,
- literature benchmark claims,
- benchmark-table population without an identical-protocol reference.
