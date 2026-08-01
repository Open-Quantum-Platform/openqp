# MRSF-TDDFT Analytic NAC — the nac-lagrangian rebuild — HANDOFF

**Branch:** `nac-lagrangian` (from Alireza's `nac` tip b71b864; local
worktree `sessions/20260731_nac_audit/repo_nac` — NB the plain `repo`
clone was taken over by a parallel session, do not touch it; chc3 mirror
`~/nac_audit/repo`, synced by git bundles; publish target is
`karmachoi/openqp-private` branch `nac-lagrangian` — NOT pushed to
origin/Alireza). Checkpoint `ae0bf33` was 65 commits on top of the `nac`
tip; the 7.54–7.56 work described below is the next private-branch
checkpoint. **Full log:** `MRSF_NAC_DERIVATION.md` (theory Secs. 0–6,
campaign 7.1–7.56).


## CAMPAIGN STATE AT 2026-08-01 — 7.54–7.56 RESPONSE CLOSED; SEAM OPEN

READ 7.52-7.56 FIRST. The former 0.0102 “fold-sector residual” is now
identified exactly: v19's DM-only `hf_energy` probe included JK response
but omitted the KS XC-kernel response because J/K reads `DM_A/B` while
`calc_dft_xc -> dftexcor` rebuilds its density from `VEC_MO_A/B`.

Source-level result: `mrsfcbc -> int2_driver -> mrsfmntoia` and the SPC/
`trans` machinery have no ground-state Fock/density dependency. The ONLY
one is `mrsfesum(wrk1,fa,fb,amo)`. Hence there is no new sigma channel:

```
P^{s,x} = C(U^x O_s + O_s U^{x,T})C^T
F^{s,x} = G^s_JK[P^x] + sum_t f_xc^{s,t} P^{t,x}
Delta sigma_xc = mrsfesum_X(Delta F^alpha_xc, Delta F^beta_xc)
```

v20 H2O/c5 gate: actual-F injection closes `trueF` to 2.23e-7; the
missing XC contribution has norm 1.016837e-2 and closes the v19 residual
to 2.65e-7. Existing analytic `get_response_packed` agrees with the
independent MO+DM central difference at 2.06e-10 (alpha) / 2.05e-10
(beta), and its sigma injection closes to 2.65e-7. The theory-level
F-channel derivation is DONE and the analytic JK+XC response has been
propagated into the v2 reference assembly.

The first H2O v21 verdict against the old `dx=1e-3` numerical referee
looked closed (5.37e-4, 5.35e-5, 7.80e-4), but that referee was not
converged. New `dx=5e-4` and `2.5e-4` freezes agree to 3.08e-7,
1.32e-9, and 9.23e-6 for pairs (1,2), (1,3), and (2,3). Against the
converged `dx=2.5e-4` referee, v2 errors are 5.813e-4, 1.958e-6, and
3.419e-3. A same-process direct-U assembly remains tight at 3.192e-5,
1.992e-5, and 2.645e-4. Therefore the remaining structural error is NOT
the amplitude/F-response formula: it is introduced when direct U is
replaced by the z-vector interchange seam, dominated by ordered pair
(3,2), whose J4 mismatch is 6.551e-3. Forcing NAC-only MINRES to a
1e-10 residual leaves the result unchanged and falsifies loose CG
convergence as the cause.

CURRENT WORKING-TREE ENTRY POINTS:
- `v20_fold_audit.py` — private-worker decomposition + actual Fock,
  DM-only JK, MO+DM, and analytic JK+XC gates.
- `mrsf_nac_response` in `tdhf_mrsf_energy.F90` and `include/oqp.h` —
  packed `nac_dm1_a/b -> nac_v1_a/b` JK+XC response using
  `scf_addons:get_response_packed`, now consumed by v2.
- `mrsf_nac_wpair_impl` — resident Fortran reference harvest. It removes
  the Python O(nbf^2) orbital-generator loop and reproduces the previous
  Python vector to 1.27e-9 on H2O. It still uses central orbital
  generators internally and is NOT the final closed-form adjoint.
- `v21_production_gate.py` — sign-resolved production/referee gate;
  latest H2O resident-Fortran artifact is
  `H2O_energy_tlf0_v25_fortran_wpair.npz` on chc3.

NEXT SESSION START SENTENCE:
"7.54–7.56을 읽고 수렴된 H2O dx=2.5e-4 기준으로 (3,2) seam/interchange
오차를 Fortran에서 유도·수정한 뒤, wpair 중앙차분 수확기를 닫힌형식
adjoint로 교체하고 ethylene/Acrolein을 게이트해."

THE COMPLETE CERTIFIED ACCOUNTING (all referees frozen):
  d = ampdir(dX) + gamma:(Sk+U)          [machine, both molecules]
  U-channel = stagedC                    [closed 4e-4]
  F-channel = mrsfesum_X(JK[P^x]+f_xc[P^x]) [closed 2.65e-7 on v20 c5]
  production v2 reference: full response propagated; direct-U closes,
  but seam replacement still leaves 3.419e-3 on H2O (2,3).

THREE REAL BUGS FIXED/IDENTIFIED THIS CAMPAIGN: tlf=0 diagonal minors
(2242bff; duplicate fix exists elsewhere -- diff before upstreaming);
the Sk-record diagonal (7.49; rebuild from overlap_ao); and the
`mrsf_nac_wpair` C/Fortran symbol collision (internal scaffold renamed).

## What is PROVEN (do not re-litigate; regression tests enforce)

| Result | Evidence |
|---|---|
| Conventions: `dcv[i,j]=d_ij` antisym, `nacv=(E_j−E_i)d` sym | absolute-orientation gate vs a code-independent exact biorthogonal oracle; suite green |
| numpy 2-D tagarrays are TRANSPOSED Fortran matrices | the root of a whole trap family; storage-boundary fix 60c5412 |
| Davidson state phases random per run | 3-run experiment; gauge-product rule (=+1) in the tests |
| compute_states_overlap exactly replicated | 1.9e-16 (`nac_formula_kernel.replica_S`) |
| γ^formula closed form (one pass) | == generator sweep to 1e-13 (`gamma_closed`) |
| **Master decomposition** `d_num = antisym[∂S/∂X̃·dX̃ + γ^formula:T]`, `T=dM/dx` | H2O 0.1%; **C1 ethylene 0.02–1.2%, cos +0.9996..+1.0000000, signed, per-component** |
| tlf=0 exact-overlap path fixed (GAMESS-inherited diagonal-minor bug) | bounds-check + valgrind; tlf0 vs tlf2 numerical NACs agree to 1e-10 |
| Raw-frame matvec PT numerically impossible | Rayleigh 2nd-order: ‖δ‖²×(spectral radius ~300 Ha); quantifies why transport/interchange is mandatory |
| z-vector interchange OPERATIONAL | polarized-L RHS → sfrolhs → gZ−gS seam == direct L∘U^x to ~1% |
| Amplitude response closes scan-free on H2O | zB + same-space L + inter-sym term: 95/64/93%, all signs +1 as derived |
| Same-space rotations need NO Fock rebuild | D invariant under same-space rotations (derived + used) |

## Benchmark table (v3g ANALYTIC-STRUCTURE assembly, 7.27-7.28, signed)

| system | pair | cos(pred,num) | maxdiff |
|---|---|---|---|
| H2O (C2v) | (1,2) | phase-gauge exact | \|pred\|==\|num\| all digits |
| H2O | (1,3) | **+1.00000000** | **4.2e-9 (machine)** |
| H2O | (2,3) | -0.99999999 (phase) | 0.12% |
| ethylene (C1) | (1,2) | +0.99999979 | 9.6e-5 |
| ethylene | (1,3) | **+1.00000000** | 1.2e-5 |
| ethylene | (2,3) gap=10.2mHa | +0.99999997 | 4.6e-4 (0.012%) |

(previous v4 all-FD gate: 0.02-1.2%; the analytic structure is BETTER
because the resolvent removes the stacked-FD near-degeneracy noise)

## Deliverables

**PRODUCTION REFERENCE (rewired, v2; not final):**
- `pyoqp/oqp/library/nac_analytic.py` — `analytic_nac(mol)`: closed-form
  γ + Sk export (rebuilt from overlap_ao per 7.49!), full-A-free MINRES
  ỹ (certified 2.6e-9), slot-injection T1 engines (certified 1e-5),
  resident-Fortran MT harvest, full JK+XC response, direct-injection seam,
  V-mask S^x eliminations, γ:Sk, antisymmetrize. Returns (nacv, dcv).
  Honest current H2O/dx=2.5e-4 max errors: 5.813e-4, 1.958e-6,
  3.419e-3; the seam/interchange defect keeps this a reference path.
- `pyoqp/oqp/library/nac_kernel.py` — copy of `nac_formula_kernel.py`.
- `pyoqp/oqp/library/single_point.py` — `analytical_nac()` REPLACED
  (old pseudo-state polarization removed) with a wrapper calling
  `nac_analytic.analytic_nac`; reports `analytic-v2-reference`.
- `source/modules/tdhf_mrsf_gradient.F90` — NAC_DUMP_DS (dbg_dsket/
  dbg_dsfull, bfnrm applied), NAC_DUMP_PIJ (dbg_pij_a/b), and the
  resident `mrsf_nac_wpair` reference harvest. Its internal implementation
  name is `mrsf_nac_wpair_impl` to avoid collision with the C symbol.
- `source/modules/tdhf_mrsf_energy.F90` + `include/oqp.h` —
  `mrsf_nac_response`: apply existing analytic JK+XC response to packed
  first-order alpha/beta densities and export packed response Focks.
- `source/modules/tdhf_mrsf_z_vector.F90` — NAC interchange solves use
  MINRES with residual tolerance at most 1e-10. This improves solver
  safety but did not change the H2O seam error.

**Kernel + gates (this directory):**
- `nac_formula_kernel.py` — replica, staged Fortran driver, one-pass
  γ^formula. Imported by the tests.
- `ROUTE_A_SPEC.md` — the derivative-sigma amp-channel blueprint,
  rewiring plan, validation ladder, landmine checklist.
- Gate harnesses: `assembly_gate.py` (master referee), `assembly_v3*.py`
  (v3–v3h incl. `v3g` = the analytic-structure benchmark), `v4a_adjoint`,
  `v5*`–`v18*` (the forensic chain; `v7j` = unified one-process gate,
  `v17_fix` = the Sk-fix closure test, `v19_fchan.py` = the original
  F-channel referee, `v20_fold_audit.py` = the decisive JK-vs-XC audit
  and analytic-response gate), `skel_gate.py` (subprocess worker), plus the
  ladder/γ/orientation gates from the first half.
- `ETH_energy.inp` — the C1 validation input.
- Frozen data (session dir + chc3 `~/nac_audit/probe/`):
  `H2O_energy_tlf0_{v7h,v7i,v7o}.npz`, `ETH_energy_{v7o,ctx,dnum}.npz`,
  `H2O_energy_tlf0_dnum.npz`, `H2O_energy_tlf0_c5_v20.npz` and
  `v20_h2o_analytic.out`; worker dirs `H2O_energy_tlf0_skel/`,
  `ETH_energy_skel/` (displaced p{idx}.inp reusable; regenerate ref.npz
  per process for phase consistency). Plus the first-half references
  (`h2o_nac_reference.npz`, `h2o_nac_ref_tlf0.npz`, `ladderA8_data.npz`,
  `gamma_formula_h2o.npz`).
- Repo tests: `tests/test_nac_convention_numeric.py`,
  `tests/test_nac_gamma_exact.py`, `tests/test_nac_formula_kernel.py`
  (suite: 11 passed).

## Remaining work (in order; every step has a frozen referee)

1. **DONE — derive/gate/propagate the fold-sector term (7.54–7.55).** It is
   `f_xc[P^x]` passed through the existing `mrsfesum` fold, not a new
   MRSF/SPC channel. v20 analytic closure: 2.65e-7.
2. **Fix the seam/interchange contraction in Fortran.** Direct-U is
   already at 2.65e-4 against the converged referee; seam assembly is
   3.419e-3. Start with ordered (3,2) and the three ROHF rotation blocks.
3. **Replace `mrsf_nac_wpair` central harvest with the closed-form
   bilinear adjoint in Fortran.** Preserve its C/tagarray interface and
   gate against the resident reference result (1.27e-9 current match).
4. Re-judge/refreeze H2O and ethylene, then sum-rule/translational checks;
   Acrolein NAMD shakedown (ndtlf=2,
   /bighome/jin/Projects/MRSF_SOC_NAMD/Acrolein/).
5. Before upstreaming: diff our 2242bff tlf=0 fix against the duplicate
   fixes in the 07-17/18 sessions (see the related-session map).

Historical context for (1)-(2): C2's empirical program closed 7.22–7.23
(vv/dd and target-ordering ruled out; near-degeneracy sharp; production
prescription = ONE z-vector with the combined RHS Ltot = L + gap·γ_a,
never term-by-term splitting — the split gates A2–A11 certify
ingredients only).

## Landmines (all cost a wrong result once; all now understood)

- Transpose EVERY 2-D tagarray read/write (`.T` or explicit F-order).
- `_run_oqp_external` prefers PATH `openqp` over the launching venv.
- The matvec needs `int2e_cutoff=1e-20` for column-exact linearity.
- Symmetric-generator staging: `(Ce^{tS})^T = e^{+tS}W` (opposite sign
  from the antisymmetric case).
- The literal ov_exact index layouts (overwrite semantics) differ from
  the clean set-theoretic minors at socc edges — use `_ia_maps`.
- C2v H2O makes cosines degenerate (1-D irreps) — judge by ratios there;
  C1 ethylene is the real test.
- Davidson phases are random per run — gauge-resolve, enforce the
  pair-sign product rule.
- **Sk record (7.49, a REAL bug):** `OQP::overlap_mo_non_orthogonal`
  from a SECOND same-geometry call has a wrong DIAGONAL. Always rebuild
  `Sk = C0^T · (OQP::overlap_ao_non_orthogonal) · C0`. One line; 50x.
- Fortran-CREATED tagarray records reach numpy transposed; PYTHON-created
  records keep their dims verbatim (certified via nac_trden_mo echo).
  Know which side created the record before transposing.
- FD sweep data (Ux, w_ref, dX) is valid ONLY in-process: displaced-
  geometry orbital gauge is run-nondeterministic. Never mix npz sweeps
  from different processes at the vector level.
- In-process 1-iter SCF emulation is broken (trueG blowup) — displaced
  F-channel measurements need SEPARATE-PROCESS workers with phases
  generated by the parent process (v19/v20 pattern).
- **DM-only DFT response is incomplete (7.54):** `fock_jk` reads
  `DM_A/B`, but `calc_dft_xc -> dftexcor` reads `VEC_MO_A/B`. Perturbing
  DM and calling one-iteration `hf_energy` gives JK response only and
  silently omits `f_xc[P^x]`. Use `get_response_packed` (analytic) or a
  consistent MO+DM central difference. Never label the DM-only result
  the full KS response.
- `sfrorhs` adds 2FT terms absent from the matvec — NAC_ZERO_2FT env.
- `int2e_cutoff=1e-20` for matvec linearity; get_jacobi is umrsf-only;
  the matvec depends ONLY on VEC_MO+FOCK records.
- Do NOT polarize the gradient chain state-specifically (4-term
  polarization needs +1/2·g(0); ruled out quantitatively in 7.30).

## Related-session map (2026-08-01, CORRECTED after a wider sweep)

An extended sweep (keywords: nacme/tlf/analytic NAC/Alireza) DID find
substantive related chats missed by the first pass:

1. "Fix same ov_exact OOB bug in native OpenQP" (2026-07-17) and
   "Fix states_overlap heap abort on gfortran-11/Linux" (2026-07-18):
   THE SAME tlf=0 diagonal-minor bug this campaign re-discovered and
   fixed (2242bff) was found and fixed there two weeks earlier (DFTB
   worktrees / native OpenQP). ACTION ITEM: before upstreaming or
   merging, diff our 2242bff against those fixes (possible duplicate/
   conflicting patches in openqp-dftb worktrees or upstream PRs).
2. "NAMD-QMMM" (2026-06-17, /Volumes/External_Storage/claude/NAMD-QMMM):
   origin of the canonical NAC conventions used by this campaign's
   audit brief (d_ij antisym, h_ij = gap*d_ij, (S-S^T)/2dt), plus real
   NAMD trajectories: Acrolein ndtlf=2 runs (with/without alignment,
   /bighome/jin/Projects/MRSF_SOC_NAMD/Acrolein/) -- a natural
   REAL-WORLD VALIDATION TARGET for the analytic NAC once rewired.
3. "Audit OpenQP NAMD/SOC/QMMM energy conservation" (2026-07-10):
   NACME/state-overlap machinery audit (Python layer, md.soc flags).
4. DFTB chats: NACME plumbing analysis (dc=(S-S^T)/dt, single_point
   L1669-1711) -- the production numerical pipeline map.

Clarification kept from the first pass: the "Phase 11/12" results the
derivation cites are ALIREZA'S `nac` BRANCH COMMITS -- in-repo, not a
chat. Key clarification: the "Phase 11/12" results the
derivation cites (fcac55a, 90943f8, the engine FD-validations) are
ALIREZA'S `nac` BRANCH COMMITS -- already in this repo's history.

Also adjacent (pointers): the two bounds-fix spin-off sessions
(upstreaming, independent); "Hessian and MRSF-TDDFT status"
(2026-07-01: CPHF infrastructure exists -- cphf_* selftests,
hf_hessian); QMRSF chats are unrelated (QMRSF_NACT is an orbital
count).
