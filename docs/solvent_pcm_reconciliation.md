# PCM dual-path reconciliation task

Recorded after integrating the local ddX/Fortran PCM energy path onto the
private solvent line. This note exists so the duplication below is reconciled
**before** the PCM backend is treated as production-ready.

## Current state

- Private branch: `feat/solvent-backend-spike` (in repo `openqp-private`).
- Merge commit: `cdddfba` — a `--no-ff` merge preserving both histories:
  - private tip parent: `1f5b173`
  - local PCM parent: `9d3026e`
- The push to `private` was **fast-forward**. `origin` and `upstream` were **not touched**.
- Build/test status at the merge: **compile + link clean; full Python suite
  passing; the ddX-gated PCM-on test (`tests/test_pcm_energy_path.py`) is
  skipped as expected** because this environment has no ddX (`OQP_ENABLE_DDX`).

## The duplication — two PCM paths now coexist

Both are present in `source/scf_addons.F90 :: calc_jk_xc` and in
`scf_energy_t`. They are **independently gated**, so they do not double-fire and
tests pass — this is **not an immediate correctness failure** — but it is a
design redundancy that must not remain permanently.

| | Private path | New ddX/Fortran path |
| --- | --- | --- |
| Energy field | `epcm` | `e_pcm` |
| Fock hook | `add_reference_pcm_reaction_field` | `add_pcm_reaction_field` (`source/solvent_pcm.F90`) |
| Gate | `present(pcm_reaction_potential)` (optional arg) | `infos%control%pcm_enabled` |
| Reaction field source | supplied from `pyoqp/oqp/library/solvent.py` | computed in-Fortran via the Fortran/C ddX adapter (`solvent_ddx_adapter.c`) |

Because the gates are disjoint, neither path fires unless its own trigger is
set, so the merged tree builds and the combined test suite (private's
`solvent.py` helpers plus the new energy-path tests) passes.

## The decision to make (do NOT implement yet)

Pick the canonical architecture, then collapse the other path into it:

1. **Python-supplied reaction potential** — Python (`solvent.py`) owns the ddX
   solve and hands a packed AO reaction potential to Fortran
   (`add_reference_pcm_reaction_field`). Fortran stays thin.
2. **Fortran-side ddX solve** — Fortran owns the full loop through the ddX C
   adapter (`add_pcm_reaction_field`); Python only sets `pcm_enabled`/`epsilon`.
3. **Layered** — Python configures (solvent, radii, model, epsilon); Fortran
   executes the per-SCF-iteration solve. Likely the cleanest, but defines a
   single owner for each responsibility and one energy field.

Whichever is chosen: unify to **one** energy field and **one** hook, and fold
the redundant tests/docs together. Until then, treat PCM as **not
production-ready**, consistent with the provisional-convention inventory in
`docs/solvent_pcm_validation_gate.md`.

## Scope guardrails (unchanged)

- No unification performed here — this note only records the task.
- No PR, no GitHub issue, no push to `origin`/`upstream`.
- Gradients, TDDFT, MRSF, and excited-state solvation remain out of scope until
  the single canonical ground-state energy path is validated.

## Related

- `docs/solvent_pcm_validation_gate.md` — provisional conventions + ddX validation gate.
- `docs/solvent_ddx_scf_integration_seam.md` — seam mapping and ddX API findings.
