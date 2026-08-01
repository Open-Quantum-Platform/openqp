# MRSF-TDDFT Analytic NAC — the nac-lagrangian rebuild — HANDOFF

**Branch:** `nac-lagrangian` (from Alireza's `nac` tip b71b864; local
worktree `sessions/20260731_nac_audit/repo_nac` — NB the plain `repo`
clone was taken over by a parallel session, do not touch it; chc3 mirror
`~/nac_audit/repo`, synced by git bundles; **PUSHED to
`karmachoi/openqp-dev-private` branch `nac-lagrangian`** — NOT pushed to
origin/Alireza). 64 commits on top of the `nac` tip, every claim
numeric-gated.
**Full log:** `MRSF_NAC_DERIVATION.md` (theory Secs. 0–6, campaign 7.1–7.53).


## CAMPAIGN STATE AT 2026-08-01 SESSION CLOSE (read 7.44-7.53 first)

THE SINGLE REMAINING ITEM for theory-level closure: the fold-sector
(~1%-of-channel; ~2e-2 on d for the worst pair) Fock-response term --
7.52. Derive it ON PAPER from the sigma source's SOCC/spc usage
(tdhf_mrsf_energy.F90 fa/fb socc combinations, mrsfmntoia/mrsfcbc spc
channels, the `trans(xvec_dim,2)` U-matrix pairing of JCP 158, 194105,
and the lr1/lr2 fold application), then gate with v19_fchan.py
(slot-resolved; targets: the 0.0102 residual on slot 4 = LR1 (i=5,a=5)
and the i=6 socc2 rows).

7.53 (the last act of this session) DERIVED AND FALSIFIED the one
shortcut candidate: the antisym-dD hypothesis (dD = C[U_a,D°]C^T from a
Loewdin-guess premise). Offline test on the frozen v7h matrices: (1,3)
improved 8.7e-3->7.0e-3 and assembly-(2,3) 0.205->0.120, BUT (1,2)
worsened 1.9e-2->1.3e-1 and J1-(2,3)->0.54; the premise is independently
contradicted by the G-A gate bound (engines == worker-w_skel at
1e-5..3e-4, so no O(h) orthonormalization shift exists in the workers).
CONCLUSION: the symmetric-dD channel is REAL under the frozen referee
protocol; there is no shortcut around reading the sigma source. That
reading is a focused fresh-context task — start the next session with:
"Read MRSF_NAC_DERIVATION.md 7.52-7.53, then read the socc/spc sigma
assembly in tdhf_mrsf_energy.F90 + tdhf_mrsf_lib.F90 (mrsfmntoia,
mrsfcbc, spc scaling, trans pairing, the fold) and derive the
fold-sector response term; gate with v19_fchan.py."

THE COMPLETE CERTIFIED ACCOUNTING (all referees frozen):
  d = ampdir(dX) + gamma:(Sk+U)          [machine, both molecules]
  U-channel = stagedC [closed 4e-4] + F-channel [dD-model+G, 89%]
            + fold-sector residual [0.010, THE item]
  production: T1 engines + direct-injection seam + V-mask + gamma:Sk
  (v1 rewired & verified); MINRES ytil; G[P] build; all certified.

TWO REAL BUGS FIXED THIS CAMPAIGN: tlf=0 diagonal minors (2242bff, NB
duplicate fix exists in other sessions -- diff before upstreaming);
the Sk-record diagonal (7.49: NEVER use overlap_mo_non_orthogonal's
second-call output for Sk -- rebuild from overlap_ao. One-line fix,
50x vector-channel improvement).

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

**PRODUCTION (rewired, v1):**
- `pyoqp/oqp/library/nac_analytic.py` — `analytic_nac(mol)`: closed-form
  γ + Sk export (rebuilt from overlap_ao per 7.49!), full-A-free MINRES
  ỹ (certified 2.6e-9), slot-injection T1 engines (certified 1e-5),
  polarized-combined seam T2 (seam(e_pq) = −U_full), V-mask S^x
  eliminations, γ:Sk, antisymmetrize. Returns (nacv, dcv).
- `pyoqp/oqp/library/nac_kernel.py` — copy of `nac_formula_kernel.py`.
- `pyoqp/oqp/library/single_point.py` — `analytical_nac()` REPLACED
  (old pseudo-state polarization removed) with a wrapper calling
  `nac_analytic.analytic_nac`. Known v1 accuracy (H2O, vs numerical):
  (1,3) 8.7e-3, (1,2) 1.9e-2, (2,3) 0.21 — limited ONLY by the
  fold-sector term; propagate the derived term here after closure.
- `source/modules/tdhf_mrsf_gradient.F90` — NAC_DUMP_DS (dbg_dsket/
  dbg_dsfull, bfnrm applied), NAC_DUMP_PIJ (dbg_pij_a/b), and the
  `mrsf_nac_wpair` scaffold (aborts by design until the fold term is
  derived; term checklist in comments).

**Kernel + gates (this directory):**
- `nac_formula_kernel.py` — replica, staged Fortran driver, one-pass
  γ^formula. Imported by the tests.
- `ROUTE_A_SPEC.md` — the derivative-sigma amp-channel blueprint,
  rewiring plan, validation ladder, landmine checklist.
- Gate harnesses: `assembly_gate.py` (master referee), `assembly_v3*.py`
  (v3–v3h incl. `v3g` = the analytic-structure benchmark), `v4a_adjoint`,
  `v5*`–`v18*` (the forensic chain; `v7j` = unified one-process gate,
  `v17_fix` = the Sk-fix closure test, `v19_fchan.py` = THE live referee
  for the fold term), `skel_gate.py` (subprocess worker), plus the
  ladder/γ/orientation gates from the first half.
- `ETH_energy.inp` — the C1 validation input.
- Frozen data (session dir + chc3 `~/nac_audit/probe/`):
  `H2O_energy_tlf0_{v7h,v7i,v7o}.npz`, `ETH_energy_{v7o,ctx,dnum}.npz`,
  `H2O_energy_tlf0_dnum.npz`; worker dirs `H2O_energy_tlf0_skel/`,
  `ETH_energy_skel/` (displaced p{idx}.inp reusable; regenerate ref.npz
  per process for phase consistency). Plus the first-half references
  (`h2o_nac_reference.npz`, `h2o_nac_ref_tlf0.npz`, `ladderA8_data.npz`,
  `gamma_formula_h2o.npz`).
- Repo tests: `tests/test_nac_convention_numeric.py`,
  `tests/test_nac_gamma_exact.py`, `tests/test_nac_formula_kernel.py`
  (suite: 11 passed).

## Remaining work (in order; every step has a frozen referee)

1. **THE fold-sector term (theory closure — the only derivation left).**
   Read the sigma source (tdhf_mrsf_energy.F90 socc/fa/fb assembly,
   tdhf_mrsf_lib.F90 mrsfmntoia/mrsfcbc/spc, the trans pairing, the
   lr1/lr2 fold), derive the Fock-response term the dD-model misses
   (7.52: F-channel 89% closed, 0.0102 vector residual on LR1 + socc2
   rows), gate slot-resolved with `v19_fchan.py`. 7.53's antisym-dD
   shortcut is FALSIFIED — do not retry it.
2. **Propagate into `nac_analytic.py`** (production), re-judge H2O AND
   ethylene against d_num, refreeze references, extend the suite.
   (The v1 rewiring is already in place with documented accuracy; this
   step upgrades it to closure. The original branch's mistake — rewiring
   BEFORE closure — is avoided because v1's accuracy is honestly
   documented in-code and here.)
3. **`mrsf_nac_wpair`**: fill the Fortran scaffold with the derived
   term(s) for a resident implementation.
4. Sum-rule + translational checks; Acrolein NAMD shakedown (ndtlf=2,
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
- In-process 1-iter SCF emulation is broken (trueG blowup) — G[P] and
  F-channel measurements need SEPARATE-PROCESS workers (v19 pattern:
  push DM_A/B ← D+εP, oqp.hf_energy at control.maxit=1 in a fresh
  process, G=(F_eps−F0)/eps; certified linear/self-adjoint).
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
