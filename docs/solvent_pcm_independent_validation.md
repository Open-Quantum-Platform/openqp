# Independent ddPCM validation — OpenQP/ddX vs PySCF ddPCM

True Tier-2 validation against an **independent** ddPCM implementation
(PySCF's native `ddPCM`, a different codebase from OpenQP/ddX), at a matched
protocol. Reproduce with `scripts/pcm_independent_ddpcm_validation.py`
(reference/diagnostic Python only; the production PCM path stays Fortran-side).

## Protocol (matched)

`model = ddPCM`, `eps = 78.3553`, `lmax = 8`, Lebedev order 29 (302 cavity
points), `eta = 0.1`, Bondi-by-Z radii identical to
`source/solvent_ddx_adapter.c::vdw_radius_angstrom`. RHF/6-31G*.

**Caveats (documented, not hidden):** PySCF ddcosmo/ddPCM requires spherical
orbitals (5d) while OpenQP 6-31G* is Cartesian (6d) — negligible for the
electrostatic shift; PySCF-ddPCM and ddX-ddPCM are independent implementations
of the same model; OpenQP currently builds the ddX source from an `l ≤ 2`
atom-centered Mulliken multipole expansion while PySCF uses the full QM density,
so a large gap signals the OpenQP source truncation.

## Observable

The unambiguous quantity both codes define identically: the total SCF
solvation shift `ΔE = E_solvated − E_vacuum`.

## Result

| molecule | OpenQP/ddX `ΔE` | PySCF ddPCM `ΔE` (independent) | gap | ratio | verdict |
| --- | --- | --- | --- | --- | --- |
| **NH₃** | −7.42 kcal/mol | −7.83 kcal/mol | 0.41 | 0.95 | within ~5 % — consistent |
| **H₂O** | −5.55 kcal/mol | −9.63 kcal/mol | **4.08** | **0.58** | **discrepant — NOT validated** |

(OpenQP: H₂O `−76.010747 → −76.019593`; NH₃ `−56.183840 → −56.195661`, from
real ddX-enabled runs on this branch.)

## Conclusion

**OpenQP's ddPCM is not yet physically validated.**

- **NH₃ agrees within ~5 %** with an independent ddPCM at matched protocol. This
  also confirms the protocol alignment and the OpenQP↔ddX coupling are sound for
  that case — the comparison is meaningful, not a gross mismatch.
- **H₂O is ~42 % too small** (a 4 kcal/mol gap). Because NH₃ agrees, this is an
  OpenQP-side problem specific to water's charge distribution, not a protocol
  artifact. The most likely cause is the **`l ≤ 2` atom-centered Mulliken
  multipole source** under-representing the water reaction field (and/or the
  energy/scale convention). This is exactly the truncation flagged in
  `docs/solvent_pcm_source_term_review.md`.

## Status of the benchmark matrix

The tier-2 `verified` rows are **regression baselines** (they lock OpenQP's own
ddX output), **not** physically validated values — H₂O disagrees with an
independent ddPCM by 4 kcal/mol. None of the rows should be presented as a
validated solvation energy.

## Next step (not done here)

Replace the `l ≤ 2` Mulliken source with a higher-`l` and/or full
density-projected `psi` so the H₂O reaction field matches an independent ddPCM,
then re-run this comparison. A pass criterion of ratio ∈ [0.95, 1.05] vs PySCF
ddPCM (as NH₃ already meets) would constitute genuine Tier-2 validation.

## Root-cause refinement (investigation update)

Further investigation changes the recommended fix. Three results:

1. **The reference is solid.** PySCF `ddPCM` ("under testing") was cross-checked
   against PySCF's mature `ddCOSMO` at the same protocol:

   | molecule | OpenQP/ddX | PySCF ddPCM | PySCF ddCOSMO | OpenQP/ddCOSMO |
   | --- | --- | --- | --- | --- |
   | H₂O | −5.55 | −9.63 | −10.41 | 0.53 |
   | NH₃ | −7.42 | −7.83 | −8.33 | 0.89 |

   ddPCM and ddCOSMO agree to ~7 %, so the ~2× H₂O shortfall is real, not a
   reference artifact.

2. **It is NOT a source-moment (dipole) error.** OpenQP's `l ≤ 1` atom-centered
   source (net Mulliken charge + atomic dipole, row partition) reproduces the
   **exact** H₂O RHF dipole to ratio **1.000** (net charges alone give 1.088×).
   So a "higher-`l` density-projected `psi`" does **not** address the gap — the
   molecular multipoles are already represented correctly.

3. **Root cause: the point-multipole source is fundamentally insufficient for a
   molecular cavity.** ddX's multipole interface places the entire solute as
   point multipoles **at the atom centers**, then evaluates `phi` at the cavity
   points from those. An atom-centered multipole expansion only converges to the
   true potential *outside* the sphere enclosing the density; the ddPCM cavity
   points sit close to (and around) neighbouring atoms, inside that region, so
   the point-multipole `phi` is wrong there regardless of `l`. PySCF instead
   evaluates the **exact** density potential at every cavity point (3-center
   integrals) and projects the density onto a volume grid for `psi`. That the
   error is larger for the compact H₂O cavity than for NH₃ is consistent with
   this near-field explanation.

### Revised recommendation

Abandon the multipole-source route; do **not** pursue higher-`l` multipoles
(they will not converge at the near cavity). The correct fix is the
**exact-density coupling** PySCF/ddCOSMO use, which OpenQP is already half-way
to: `electrostatic_potential_unweighted` already gives the **exact** `phi` at
cavity points (the primitive the retired `oqp_ddx_pcm_solve` path used). What
remains is (a) a correct density-derived `psi` (the ddCOSMO volume-grid
projection, not a multipole expansion) and (b) re-solving the ddX convergence
that drove the earlier switch to the multipole path. This is research-grade
Fortran work, not a small source-term tweak; it was therefore **not** attempted
here in preference to shipping an unvalidated change. The PySCF ddPCM/ddCOSMO
gate above is the acceptance test for it (ratio ∈ [0.95, 1.05]).

