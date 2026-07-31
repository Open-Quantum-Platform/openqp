# MRSF-TDDFT Analytic NAC — the nac-lagrangian rebuild — HANDOFF

**Branch:** `nac-lagrangian` (from Alireza's `nac` tip b71b864; local clone
`sessions/20260731_nac_audit/repo`, chc3 mirror `~/nac_audit/repo`, synced
by git bundles; pushed nowhere). 25+ commits, every claim numeric-gated.
**Full log:** `MRSF_NAC_DERIVATION.md` (theory Secs. 0–6, campaign 7.1–7.21).

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

## Benchmark table (v4 assembly-gate, signed, no alignment)

| system | pair | cos(pred,num) | \|d_pred\| / \|d_num\| |
|---|---|---|---|
| H2O (C2v) | (1,2) | +0.99999 | 0.12661 / 0.12636 |
| H2O | (1,3) | +1.00000 | 0.03685 / 0.03679 |
| H2O | (2,3) | +0.99999 | 0.52900 / 0.52861 |
| ethylene (C1) | (1,2) | +0.99998 | 0.21751 / 0.21489 |
| ethylene | (1,3) | **+1.00000000** | 0.95505 / 0.95519 |
| ethylene | (2,3) | +0.99961 | 3.71479 / 3.68050 |

## Deliverables in this directory

- `nac_formula_kernel.py` — the production kernel library (replica,
  staged Fortran driver, one-pass γ^formula). Imported by the tests.
- Gate harnesses: `assembly_gate.py` (the master referee), `gamma_gate.py`,
  `formula_gamma.py`, `orientation_gate.py`, `ladderA*_gate.py` (A2–A10),
  `gamma_closed.py`, `diagA.py`, `conv_check.py`, `freeze_ref.py`.
- `ETH_energy.inp` — the C1 validation input.
- Frozen references (chc3 `~/nac_audit/`): `h2o_nac_reference.npz`,
  `h2o_nac_ref_tlf0.npz`, `ladderA8_data.npz` (H2O La/Ls/Ux/residuals),
  `gamma_formula_h2o.npz`.
- Repo tests: `tests/test_nac_convention_numeric.py`,
  `tests/test_nac_gamma_exact.py`, `tests/test_nac_formula_kernel.py`.

## Remaining work (assembly only; every step has a frozen referee)

1. **C2 finish** — ethylene response-closure bookkeeping. Item (a) is
   RESOLVED-NEGATIVE (the unrestricted sweep, `ladderA10.out` +
   `ladderA10_eth_full.npz`, changes nothing: vv/dd are not the missing
   content). Open, sharpened: (b) near-degenerate conditioning — (2,3)'s
   zB is 2.2x the residual regardless of sign; (c) a MOLECULE-LEVEL zB
   sign in the seam chain — ethylene wants s1=-1 on ALL pairs (93.4% for
   (1,2)) while H2O wants +1 on all; inspect the polarized-RHS push/pack
   orientation vs the seam gradient sign and the set_mrsf_nac_cphf
   target ordering. All sweep data saved in the npz — no re-sweep needed.
2. **Full analytic d assembly** — combine: skeleton (`mrsf_nac_amp` +
   `mrsf_nac_esum`, both FD-validated) + zB (operational) + same-space L
   (from `mrsf_matvec_apply`, restricted generators; Fortran-resident
   version = assemble from the G_MO/gchan/fa/fb exports) + inter-sym term
   + γ^closed:Sk (push γ as `OQP::nac_gamma_tlf` → `mrsf_nac_overlap`'s
   dSket contraction). Gate every sub-swap against `assembly_gate`.
3. **Rewire `analytical_nac()`** once (2) closes on BOTH H2O and ethylene.
   Do NOT rewire before — shipping a partially-closed response was the
   original branch's central mistake.
4. Sum-rule + translational checks; refreeze references; extend the suite.

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
