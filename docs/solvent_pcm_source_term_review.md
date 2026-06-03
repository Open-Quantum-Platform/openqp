# Review — QM source-term fix (commits `0f66593…7e93ca8`)

Reviewer pass over hermes's 5 commits that land the QM-SCF PCM source-term fix,
checked against `docs/solvent_pcm_source_term_handoff.md` and the
self-consistency validation specified there. **Review only — no behavior was
changed by this note.**

## Verdict

The **core source-term fix is real and faithful to the handoff** and resolves
the nuclear-only bug for small molecules. The **benchmark matrix overreaches**
and should not be treated as scientific validation: its references are OpenQP's
own output (a regression lock, honestly labeled), several entries are
physically wrong yet marked `verified`, and the scope went well past the
requested "H₂O first".

## What matches the spec (accept)

- **Approach B, exactly as specified.** A per-atom real-solid-harmonic source
  (monopole = net atomic charge via a Mulliken row partition, plus electronic
  dipole + quadrupole, `l ≤ 2`, `nmultipoles = 9`) is fed to **one** ddX call
  that drives **both** `phi` and `psi` (`ddx_multipole_electrostatics` +
  `ddx_multipole_psi`). Production stays Fortran-side; the new C entry points
  (`oqp_ddx_pcm_solve_multipole_source[_with_phi]`) only wrap ddX.
- **The self-consistency gate is implemented and enforced.**
  `phi_source_vs_exact_{rms,max}` compares ddX's multipole-`phi` against
  OpenQP's exact `phi_cav`, and the diagnostics test asserts `rms < 0.05`. For
  H₂O it is `0.0076` (~8 %, consistent with `l = 2` truncation; normalization
  sound).
- **The central bug is fixed.** `q_cav_sum` went from `≈ +10` (total nuclear
  charge) to `≈ 1e-4`, and `source_charge_sum ≈ 1e-15` (net molecular charge),
  for every molecule. `e_pcm` is now `~1e-2` Eh, not tens of Hartree.
  H₂O: `e_pcm = -0.0094 Eh` (~ -5.9 kcal/mol).
- **FMM was kept ON** (the forbidden FMM-off "fix" was not committed);
  convergence comes from routing the QM solve through the consistent
  `ddx_ddrun` multipole path (the one Tier-1 already used).
- **Tier-1 preserved** (nuclear-only `psi` retained for the point-charge
  regression); **diagnostic label updated** to
  `psi_source=total_qm_atom_multipoles_l2`; **FD-derived Fock scale `-0.5`**
  confirmed (`fd_fock_scale_maxerr ≈ 1.6e-7`).
- **Two genuine bug fixes:** `int1.F90` now zeroes `pot` before accumulation
  (the exact `phi_elec` was previously garbage-accumulated), and the `phi_cav`
  sign was corrected to `phin + phi_elec` with rationale.
- **No out-of-scope code.** `docs/mrsf_pcm_design_note.md` is explicitly
  "design gate only, do not implement"; no MRSF/TDDFT/gradient code was added.
  Python (pyscf scripts) is reference/diagnostic only.

## Concerns (do not treat the matrix as validated)

1. **References are circular, not independent.** All 30 `reference_value`s have
   `reference_kind = openqp_ddx_closed_shell_diagnostic_baseline` — OpenQP's own
   ddX output recorded at `1e-7` tolerance. Honestly labeled ("not an
   independent literature value"), but it is a **regression lock**, not a check
   against literature or a trusted external code. A protocol-matched independent
   reference is still missing (the pyscf path is explicitly "protocol-near, not
   identical").

2. **Physically wrong values marked `verified`.** `e_pcm` is **positive** for
   CH₃CN (`+0.0338`), CH₃OH (`+0.0145`), CO₂ (`+0.0007`). Electrostatic
   solvation of polar neutrals (CH₃CN, CH₃OH) must be negative — the `l = 2`
   atom-centered Mulliken-multipole source is too crude for these and gives the
   wrong sign, yet they are blessed `verified`. This violates the acceptance
   rule "no benchmark reference unless physically meaningful."

3. **One entry fails its own gate but is still `verified`.** `ch3oh_rhf_hf` has
   `phi_source_vs_exact_rms = 0.0527 > 0.05`, so its diagnostics test would
   **fail** with a real ddX build (others sit just under threshold). The matrix
   was populated past where the self-consistency gate holds.

4. **Scope creep.** The handoff said **H₂O only; add NH₃ only after H₂O is
   clean.** This is 10 molecules × {HF, **BHHLYP**, **PBE**} = 30 entries,
   including DFT — beyond the RHF-energy-first scope.

5. **Energy convention still unreconciled (expected).** H₂O: `e_pcm = -0.0094`
   vs `½ Tr[D·V] = -0.035` (ratio ~3.7). The `e_pcm` vs `½ Tr[D·V]` bookkeeping
   remains open, and the baselines anchor the still-unvalidated `e_pcm`.

6. **Minor.** `enable_forces` flipped `0 → 1` (a ddX `ddx_ddrun` flag; forces
   are computed and discarded — slight tension with "no gradients", not an
   OpenQP gradient). Commit messages have no bodies.

## Diagnostic evidence (from `tests/data/pcm_literature_benchmarks.json`)

| molecule (RHF/HF/6-31g*) | `e_pcm` (Eh) | `q_cav_sum` | `phi_src_vs_exact_rms` | note |
| --- | --- | --- | --- | --- |
| H₂O  | -0.0094 | -1.8e-4 | 0.0076 | healthy |
| NH₃  | -0.0124 |  1.3e-4 | 0.0164 | healthy |
| HF   | -0.0042 | -6.8e-5 | 0.0058 | healthy |
| CH₄  | -0.0030 |  1.5e-4 | 0.0167 | healthy |
| CO   | -0.0011 | -2.7e-5 | 0.0060 | small (nonpolar) |
| CO₂  | **+0.0007** | -7.9e-5 | 0.0052 | wrong sign |
| HCN  | -0.0131 |  1.1e-5 | 0.0096 | ok |
| H₂CO | -0.0110 | -1.7e-5 | 0.0110 | ok |
| CH₃OH| **+0.0145** |  2.5e-4 | **0.0527** | wrong sign + fails rms gate |
| CH₃CN| **+0.0338** |  2.2e-5 | 0.0127 | wrong sign (large) |

`q_cav_sum` and `source_charge_sum` are `≈ 0` for all molecules — the
nuclear-only signature is gone. The positive `e_pcm` for the more polar
molecules is the `l = 2` truncation limit.

## Recommendation

- **Accept the source-term fix itself** as a genuine, self-consistency-validated
  improvement that resolves the nuclear-only bug (H₂O / NH₃ / HF / CH₄ are
  solid).
- **Do not present the 30-entry matrix as validated.** Trim to the in-scope
  small molecules, stop marking the positive / threshold-failing entries
  `verified`, keep the honest "not independent" labeling, and defer DFT and the
  larger polar set.
- **Pursue an actual protocol-matched independent reference** (pyddx or another
  ddPCM code with identical geometry / radii / discretization) for the real
  Tier-2 validation.
- The `l = 2` truncation is the next physics limiter for polar molecules
  (needs higher `l` or a density-projected `psi`).

## Update — trim applied

The benchmark matrix was relabeled (no physics changed):

- The 9 entries with **positive `e_pcm`** (CO₂, CH₃OH, CH₃CN × HF/BHHLYP/PBE)
  were demoted from `status=verified` to `diagnostic_unphysical_pending`, each
  with a `physical_review` field. Both the diagnostics and reference tests now
  **skip** these (the diagnostics test gained a `status != verified` skip), so
  no unphysical value is locked in as a passing gate. `ch3oh_rhf_hf` (which also
  exceeded the `phi_source_vs_exact_rms < 0.05` self-consistency gate) is in
  this set.
- The 21 physically-sound, gate-passing negative-`e_pcm` entries (H₂O, NH₃, HF,
  CH₄, CO, HCN, H₂CO × 3 methods) remain `verified` **regression baselines**.
- A top-level `_review` block in the JSON records the regression-baseline
  (circular, not independent) caveat and the scope note (DFT + extended molecule
  set exceed the H₂O-first instruction).

Validated with a real ddX build: H₂O diagnostics + reference **pass**, the
demoted entries **skip** with an explanatory message, and the Tier-1 adapter
regression still **passes**. A protocol-matched independent reference is still
the outstanding requirement for true Tier-2 validation; a stricter trim to the
H₂O/NH₃-only scope is available if desired.

### DFT deferral

The closed-shell **DFT** entries (BHHLYP/PBE × the 7 sound molecules, 14 rows)
were demoted `verified → deferred_out_of_scope_dft` with a `scope_review` field.
They are physically sound (negative `e_pcm`, pass the self-consistency gate) and
kept as regression baselines, but DFT-PCM is beyond the handoff's
RHF-energy-first scope and the values remain self-baselines. After this, the
`verified` gates are: the Tier-1 point-charge adapter regression **plus 7
RHF/HF closed-shell rows** (H₂O, NH₃, HF, CH₄, CO, HCN, H₂CO). Confirmed with a
real ddX build: a DFT row now skips, the RHF H₂O row still passes, Tier-1 still
passes.
