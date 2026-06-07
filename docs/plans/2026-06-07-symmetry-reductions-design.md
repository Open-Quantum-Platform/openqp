# Symmetry reductions design (integral + response space)

Status: DESIGN ONLY — no implementation authorized by this document alone.
Prerequisite: the metadata/labeling layer (PR #184) is merged and validated.
Scope: abelian (D2h-family) subgroups only, exactly what the detection layer
provides. Reductions stay behind `use_integral_symmetry` /
`use_response_symmetry`, default off; `enabled=false` and C1 molecules must
take byte-identical code paths to today.

## 0. Expected gains

For a molecule with abelian group order |G| (2, 4, or 8):

- Two-electron Fock build: up to ~|G|x fewer shell quartets via the
  petite-list/skeleton-Fock method (Dupuis & King, IJQC 11, 613 (1977)).
  Benzene (D2h, |G|=8) is the canonical win; water (C2v) gives ~4x.
- XC grid: integrate only the grids of symmetry-unique atoms and scale by
  orbit size — direct multiplier on the dftlib work in PR #182.
- Response: per-irrep Davidson subspaces — smaller subspaces, no
  cross-irrep contamination, and state targeting by symmetry.

## 1. Frame decision (must be settled first)

The labeling layer conjugates operations into the input frame so it never
moves the molecule. Reductions want the opposite: in the standard
orientation every abelian operation is a signed shell permutation (cheap,
exact); in a general frame p/d/f/g blocks are dense.

DECIDED (2026-06-07, CHC): when `use_integral_symmetry=on`, reorient the
input geometry to the detection's standard orientation at load time
(GAMESS-style). Gradients/modes are rotated back on output. When the
flag is off (default), geometry is never touched — labeling continues to
use the input-frame machinery. This keeps "off == today" exactly.

## 2. Prerequisite metadata (Gate A — Python side, small) [DONE: build_reduction_maps + tests, stored under symmetry_metadata['reduction_maps']]

From the existing detection payload, build and pass to Fortran (tagarray):

- `OQP::sym_nops` and, per operation: shell permutation map
  `OQP::sym_op_shell_map(nshell, nops)` (shell s on atom a -> same-rank
  shell on atom perm(a); the validation in `_ao_operator_maps` already
  guarantees identical shell signatures on equivalent atoms) and per-AO
  sign vector `OQP::sym_op_ao_sign(nbf, nops)` (from `_component_signs`,
  valid because standard-orientation abelian ops are sign-diagonal).
- Orbit representatives: `OQP::sym_shell_rep(nshell)` (lowest-index shell
  in each orbit) and orbit sizes.

Validation gate: round-trip test — apply the Fortran-side maps to the
overlap matrix and confirm invariance (reuse the one-electron block
diagnostics in Python on the same data).

## 3. Phase I: petite-list skeleton Fock (the core) [VALIDATED 2026-06-07]

Validation results (tests/smoke_petite_list_validation.py, RHF/ROHF/UHF
6-31G*): all 7 cases pass the <= 1e-10 Ha gate vs C1 (worst 4.1e-12,
benzene |G|=8). Engagement proven by a corrupted-map probe (wrong Fock ->
SCF divergence). UHF triplet water exercises the fail-safe: the TRAH
stability stage keeps a broken-symmetry solution, petite disables itself
(status disabled_symmetry_broken_scf) and the energy matches C1 exactly.
Screening uses the unweighted cutoff so the integral set is identical to
the C1 path.

Anchor points in `source/integrals/int2.F90` (module `int2_compute`):

- `int2_twoei` — OpenMP quartet driver: outer loop over `ij_pair` from
  `int2_build_shell_pair_map` (triangular shell pairs), inner k,l loops,
  Schwarz test per quartet, `shellquartet` evaluation, then
  `int2_consumer%update`.
- Consumers `int2_rhf_data_t` / `int2_urohf_data_t` accumulate the Fock
  matrix from quartet batches.

Plan:

1. Petite list: skip quartet (ij|kl) unless it is the canonical
   representative of its orbit under G (lexicographic comparison of the
   image quartets using the shell permutation maps; sign bookkeeping not
   needed for the skip decision). Multiply the quartet contribution by
   the orbit size q4 in the consumer update. This is a pure filter +
   weight: the Rys/libint kernels and screening are untouched.
2. Skeleton Fock symmetrization after the quartet loop:
   `F = (1/|G|) sum_op T_op F_skel T_op^T` with T_op the signed
   permutation — O(|G| * nbf^2), negligible. Density entering `shlden`
   (shell-density screening) is already totally symmetric, so
   density-weighted Schwarz screening commutes with the petite list.
3. Keep the existing OpenMP schedule; orbit-representative pairs are
   distributed like any other `ij_pair`. Load imbalance from q4 weights
   is bounded by |G| and amortizes over many quartets; measure before
   optimizing.

Validation gates (each vs the C1 reference, same build):

- RHF/ROHF/UHF energies: water (c2v), ethylene (d2h), trans-difluoro-
  ethylene (c2h), H2O2 (c2), benzene (d2h subgroup of d6h) — agreement
  to <= 1e-10 Hartree, plus iteration-count parity.
- Quartet-count telemetry: measured reduction factor reported in the log
  (claim-honest: report achieved factor, not |G|).
- A deliberately desymmetrized geometry (1e-3 distortion) must fall back
  to C1 automatically via detection tolerance.

## 4. Phase II: one-electron integrals and gradients [GRADIENTS VALIDATED 2026-06-07]

Implemented: petite filter + q4 weight in grd2_driver_gen (skeleton 2e
gradient), final-gradient projection onto the totally symmetric component
in pyoqp (exact for 1-dim irreps; all abelian irreps are 1-dim, so TD
state gradients qualify). Validated at machine precision (<= 6.1e-14
Ha/bohr vs C1): water/ethylene RHF, water BHHLYP, water MRSF S1.

Two robustness fixes from validation: (1) reorientation now ITERATES
detection until the current geometry's detection returns the identity
frame -- a single rotate-then-redetect can land in a different but
equivalent degenerate frame (e.g. the three C2 axes of d2h ethylene),
leaving maps and geometry inconsistent (this produced a variationally
impossible SCF that the stability stage then masked); (2) staging now
verifies T S T^T = S against the real overlap for every operation and
falls back to C1 on any deviation (defense-in-depth). One-electron
gradient parts are not reduced (cheap); they are covered by the final
projection.

- One-electron matrices: cheap already; symmetrize for consistency
  (Hcore from unique atom pairs is optional, low value).
- Gradient quartet loop (`grd2*.F90`) admits the same petite list with
  contributions scattered to orbit atoms via the permutation maps; the
  symmetric-group projection of forces also cleans numerical noise.
  Gate: analytic vs numerical gradient on c2v/d2h cases with the flag on.

## 5. Phase III: XC grid reduction (ties into PR #182) [VALIDATED 2026-06-07]

Implemented: xc_options_t%symAtomWeight (set only by the SCF dftexcor
path -- response/gradient consumers never see it), slice skip + orbit
weighting in both run_xc slice loops, skeleton-XC symmetrization in
calc_dft_xc, per-atom orbit weights staged as OQP::sym_atom_weight.
Validated (BHHLYP/6-31G*): water c2v 4.3e-14, benzene d2h 7.0e-12 Ha vs
orientation-matched C1. Note: XC quadrature error is orientation
dependent (~1e-6 for default grids), so C1 references must be computed
in the standard orientation -- the validation script does this.

`dftlib` builds atom-centered grids (`dft_molgrid.F90`). With symmetry:
process only symmetry-unique atoms' slices and multiply their quadrature
weights by the atom-orbit size. The density and potential are totally
symmetric, so this is exact for the XC energy and the symmetrized Fock.
Gate: XC energy match to <= 1e-10 on the Phase-I molecule set; interacts
with the slice infrastructure touched by perf/xc-numerical-kernel, so
rebase on whichever lands first.

## 6. Phase IV: response-space blocking

- The MO irrep labels already exist. Excitation pairs (i,a) partition by
  irrep product; the A/B matrices are block-diagonal across irreps.
- Davidson (`tdhf_lib.F90` / `tdhf_mrsf_lib.F90`): restrict trial vectors
  to one irrep block per requested state symmetry; orthogonalization and
  preconditioning stay unchanged within a block.
- MRSF caveat: the special spin-pairing rows (`mrsfxvec` lr1/lr2/g/d
  components) couple specific SOMO pairs; their irrep assignment must be
  derived from the SOMO labels before blocking — treat as its own gate
  with 'mixed'-orbital bailout (fall back to unblocked solver whenever
  any involved MO is 'mixed').
- Gate: state energies vs unblocked solver to <= 1e-8 on water MRSF
  (states 1A1/1B1/1A2 from the smoke test), plus a c1 fallback test.

## 7. Failure-mode policy

Any inconsistency detected at runtime (orbit map mismatch, 'mixed'
occupied orbital, detection failure, distorted geometry) silently
disables the reduction for that run and records the reason in
`symmetry_metadata` — never abort, never produce unsymmetrized/skeleton
results without the final symmetrization.

## 8. Suggested order of implementation

Gate A (maps to Fortran) -> Phase I (RHF first, then UHF/ROHF) ->
Phase III (XC, cheap and high value for DFT) -> Phase II (gradients) ->
Phase IV (response). Each phase lands behind the off-by-default flags
with its validation gate green before the next starts.
