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

Reviewed reference-PCM reaction potentials may only enter the first prototype through the non-incremental SCF Fock path. If OpenQP old-buffer state is present, the guard must fail fast with the diagnostic `reference PCM incremental Fock is not validated` and preserve whether `dens_old`, `f_old`, or both triggered the incremental Fock shortcut. Caller-facing and native diagnostics should carry exact old-buffer provenance fields (`dens_old_present=true/false`, `f_old_present=true/false`) plus ordered trigger labels (`incremental_trigger_fields=dens_old,f_old`, `incremental_trigger_fields=dens_old`, or `incremental_trigger_fields=f_old`) before any reviewed `pcm_reaction_potential_in` reaches `calc_jk_xc`. This keeps the candidate PCM energy bookkeeping from silently mixing with an unvalidated incremental-Fock reuse path and lets later SCF-call-site wiring identify which stale buffer caused the block without rerunning the audit helper.

## SCF call-site bridge contract

`reference_scf_pcm_calc_fock_call_site_bridge()` is the dependency-light staging helper for the future SCF caller. A molecule with no reviewed runtime payload remains an explicit disabled/no-payload request even when old SCF buffers are present. A reviewed payload may forward only `pcm_reaction_potential_in` through the non-incremental path, with compact `nbf` and packed AO length metadata preserved so the native caller can enforce `size(pcm_reaction_potential_in) == nbf * (nbf + 1) / 2` before any reaction field reaches the Fock builder.

The native `calc_fock(..., pcm_reaction_potential_in=...)` guard must independently validate the same shape before `calc_jk_xc`: `size(pcm_reaction_potential_in)` must match `nbf_tri`, and malformed prototype callers must fail with `reference PCM reaction potential length must match packed AO dimension` rather than forwarding a wrong-length AO reaction field into the Fock build.

## Runtime payload shape contract

`reference_scf_pcm_runtime_payload()` must produce the packed AO shape metadata (`nbf`, `packed_ao_length`, `expected_packed_ao_length`, and `packed_ao_shape_formula`) together with the reviewed `OQP::pcm_reaction_potential`. The consumer `reference_scf_pcm_reaction_potential_from_payload()` must not trust stale metadata: the consumer must recompute and validate the triangular packed AO length before exposing `pcm_reaction_potential_in` to any `calc_fock` bridge. This keeps every no-runtime handoff boundary aligned on the same `nbf * (nbf + 1) / 2` contract before native SCF wiring.

Those shape fields are required, not optional defaults. A payload with missing `nbf`, missing `packed_ao_length`, missing `expected_packed_ao_length`, or missing `packed_ao_shape_formula` must fail fast before any reviewed potential can be exposed. The consumer must also reject malformed numeric provenance such as non-integer/bool `nbf` metadata and boolean numeric entries in `OQP::pcm_reaction_potential` or `OQP::pcm_epcm`, so Python truthiness cannot masquerade as a valid packed AO payload.

The same producer/consumer boundary must also preserve backend-validation provenance while runtime PCM remains disabled. The payload field `backend_validation_status` is required before exposing `pcm_reaction_potential_in`, and its only accepted scaffold value is `pending PySCF/ddX/reference cross-check`; hand-written or stale payloads that omit this status must fail before the reviewed reaction potential can reach any SCF Fock handoff.

The reviewed consumer boundary must also reject non-mapping payload objects before field-level checks. `None`, list, tuple, or other non-mapping values must fail with the domain-specific diagnostic `PCM runtime payload must be a mapping` before exposing `pcm_reaction_potential_in`, so malformed caller state cannot fall through to misleading missing-field or shape errors.

## Molecule JSON round-trip contract

`Molecule.set_pcm_runtime_payload()` stores reviewed PCM runtime payload fields outside native Fortran tag arrays and uses an explicit allowlist for JSON round-trip persistence. The Molecule-level setter must reject non-mapping restored/prototype payloads with `PCM runtime payload must be a mapping` before applying the allowlist, so malformed JSON/restored data cannot be silently iterated or stripped before the reviewed payload consumer sees it. The allowlist must preserve `OQP::pcm_reaction_potential`, `OQP::pcm_epcm`, `nbf`, `packed_ao_length`, `expected_packed_ao_length`, `packed_ao_shape_formula`, `pcm_runtime_payload_version`, first-scope labels, and `backend_validation_status`; these fields must not be stripped during save/load because the reviewed payload consumer requires them before exposing `pcm_reaction_potential_in`.

## Reporting rules

For validation reports, include:

- gas-phase and solvated total energies in Hartree;
- polarization energy if exposed by the backend;
- backend, model, dielectric, cavity/radii settings;
- SCF iteration count and convergence status;
- whether TRAH or another fallback occurred;
- exact code commit, backend version, and artifact path.

Do not report solution-phase excitation energies as validated until the PCM-solvated ROHF reference path has live OpenQP evidence and the MRSF response-root mapping is checked separately from the reference SCF energy shift.
