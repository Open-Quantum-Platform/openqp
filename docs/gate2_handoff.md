# Gate 2 handoff — native Rys 1e nuclear-attraction 2nd derivative

Status as of this commit: **build GREEN, primary integral oracle RED-but-real
(executes, does not silently skip), overlap basis-bridge SOLVED.** This document
is the single source of truth for finishing Gate 2. It is written so the work can
be continued directly.

## Branch / remote facts (verified)

- Work branch: `feat/hf-dft-hessian`. `HEAD == origin/feat/hf-dft-hessian`.
- Gate-1 (`source/integrals/rys_deriv.F90`, validated `RYS_DERIV_SELFTEST PASS`,
  max abs 4.7e-11) exists ONLY on this branch. `claude/practical-hypatia-4tcBa`
  lacks it and reflects the superseded `der2_coul_xyz` direction — do NOT build
  Gate 2 there.
- Fix-forward only. No force-push, no history rewrite.

## Design contract (frozen, do not relitigate)

- 2e / F^X / CPHF response stay on the **libint** path. Untouched.
- 1e nuclear-attraction energy, gradient, **and 2nd derivative** stay **native
  Gauss-Rys**.
- The production basis-basis 2nd derivative uses **angular-momentum (AM) shift
  identities** and therefore does **NOT** differentiate Rys roots/weights.
  `rys_deriv.F90` is a validated auxiliary layer, not on this path.
- Per shell pair: `nroots_der2 = floor((Li+Lj+2)/2)+1` (set explicitly; never
  inherit `cp%nroots`). Abort cleanly for `nroots_der2 > 5` (any L>=4 combination
  reaching `rys_general`).
- Per-nucleus single-center operator `-Z_C/|r-C|`. Compute `p_AA(C)`, `p_AB(C)`,
  `p_BB(C)` directly; then translational invariance:
  - `p_AC = -(p_AA + p_AB)`
  - `p_BC = -(p_BB + p_AB^T)`
  - `p_CC = -(p_AC^T + p_BC^T) = p_AA + p_AB + p_AB^T + p_BB`
- Frozen/forbidden to touch: `fock_deriv.F90`, `cphf.F90`, the `*_selftest`
  response files, the validated overlap/kinetic 2nd-deriv path, the `rys_rtN`
  value routines, libint files, and guarded `hf_hessian.F90`.

## Where the code is

- `source/integrals/mod_1e_primitives.F90`
  - `comp_coulomb_der2_blocks(cp, c, znuc, pAA, pAB)` — shared AM-shift kernel,
    returns **uncontracted** per-AO-pair `(3,3,inao,jnao)` blocks
    `pAA=d2/dA_adA_b`, `pAB=d2/dA_a dB_b`. `nroots_der2` per Eq. above. PUBLIC.
  - `comp_coulomb_der2_braC(...)` — thin density contraction wrapper over the
    kernel; returns `p_XX` and `p_XC=-(p_XX+p_AB)`. Used by `hess_en`.
- `source/integrals/grd1.F90::hess_en` — assembles the 9 atom blocks + optional
  `hess_cc` diagnostic. Structure validated against the TI relations above.
- `source/modules/hess1_selftest.F90` — emits the oracle dump
  `/tmp/hess_nuc_blocks.txt` when `nbf<=64`. Dump layout (current):
  1. `nbf natom`
  2. `nshell`, then per shell: `atom0based L nprim` + `nprim` lines of `exp cc`
     (**cc are OpenQP's primitive-normalized stored coefficients**)
  3. `nbf` lines: `lx ly lz cx cy cz` (AO cart powers + center, bohr)
  4. `natom` lines: `Znuc x y z` (bohr)
  5. `nbf x nbf` overlap `OQP_SM` (OpenQP normalized)
  6. per nucleus C: `PAA[a,b]` then `PAB[a,b]` (each 3x3 of nbf x nbf),
     **chargeless (znuc=1), bfnrm-normalized**.
- `tests/test_hess_nuc_oracle.py` — primary oracle (currently RED; see below).

## THE BASIS BRIDGE (this was the blocker — now SOLVED)

The oracle must build PySCF from the **OpenQP-exported basis**, not PySCF's
library, and compare in **OpenQP's AO normalization**. Two corrections are both
required; with both, the overlap fingerprint passes to **7.5e-8** (NOT 1e-10 —
see precision note):

1. **Raw contraction coefficients.** OpenQP's stored `cc` are already
   primitive-normalized. PySCF re-normalizes, so feed it
   `cc_raw = cc_oqp / gto_norm(L, exp)` (`from pyscf.gto.mole import gto_norm`).
   This makes the radial functions identical (fixes the s-s block, which the
   √3-only d-bridge can never fix).
2. **Cartesian normalization bridge `D`.** PySCF cart self-overlaps are not 1
   (d_xx = 2.513, etc.); OpenQP AOs are fully unit-normalized. Use the
   self-overlap bridge `D_ii = 1/sqrt(S_py_ii)` (empirically exact and
   convention-proof), then compare
   `S_oqp  vs  D (P^T S_py P) D`.

   NOTE: the originally-proposed hand-coded `bf` factors (axial d = 1,
   off-axis = √3, etc.) do **not** reproduce PySCF's cart normalization and give
   a 1.5 residual. Use the self-overlap diagonal instead. (Both are "OpenQP
   normalization as target"; the self-overlap form just reads the factor off
   PySCF directly so it is exact regardless of PySCF's internal cart convention.)

3. **AO permutation `P`** by (center, lmn, radial-rank-within-(atom,L)). PySCF
   cart order per L is `lx desc, then ly desc`; OpenQP d order is
   `xx,yy,zz,xy,xz,yz`. Matching on the full (center,lmn,radial-rank) key is
   robust and was verified unique for HOF/6-31G*.

A self-contained, working reference implementation of steps 1-3 (parses the dump,
builds PySCF, proves the fingerprint) is `/tmp/stage1_raw2.py` in the dev
container — fold it into `tests/test_hess_nuc_oracle.py`.

### Precision note (why 1e-10 is currently unreachable)

The exported exponents/coeffs are **float32-origin** (e.g. `5484.671875` is an
exact float32; the `cc` are not float64-exact). So a PySCF mol built from the
dump differs from OpenQP's internal double math at the ~1e-7 level. The overlap
fingerprint bottoms out at **7.5e-8**, and that floor is OpenQP's stored-basis
precision, not the bridge. Options for hermes (pick one, document it):
- (a) accept a `1e-6` overlap-fingerprint gate, justified by the float32 basis
  provenance, and keep derivative thresholds tight relative to that; or
- (b) export the basis in full float64 from OpenQP's source-of-truth (basis
  library before single-precision storage) and keep the strict `1e-10` gate.
Do NOT silently loosen the derivative tolerance to mask a real error — only the
overlap-fingerprint floor is provenance-limited.

## PySCF derivative conventions (verified by finite difference, sign +, order [a,b])

- `int1e_ipiprinv[a,b] = + d2/dA_a dA_b <mu|1/|r-C||nu>`  (bra on A) → matches `pAA`
- `int1e_iprinvip[a,b] = + d2/dA_a dB_b <mu|1/|r-C||nu>`  (bra A, ket B) → matches `pAB`
- Evaluate inside `with mol.with_rinv_at_nucleus(C):`. Apply `-Z_C` explicitly
  (the bare intors are chargeless; the dump is also chargeless).
- The same `D (P^T · P) D` bridge applies unchanged to `int1e_rinv`,
  `int1e_ipiprinv`, `int1e_iprinvip` (normalization is displacement-independent).

## Remaining work (ordered)

1. Fold `/tmp/stage1_raw2.py` logic into `tests/test_hess_nuc_oracle.py`:
   build PySCF from exported basis (raw cc), prove overlap fingerprint FIRST
   (gate; skip→fail), then compare per nucleus:
   base `int1e_rinv`, then `pAA` vs `ipiprinv`, `pAB` vs `iprinvip`, with explicit
   `a!=b` component checks and the `-Z_C` factor. `p_BB` is covered by the full
   nbf×nbf sweep (it is `pAA` of the swapped pair).
2. Make skip impossible to misreport: the runner must emit
   `build_rc / oracle_rc / skipped / failed / errored`; `skipped>0` on the
   primary oracle ⇒ Gate RED.
3. Only after the integral oracle passes element-wise, run the contracted
   nucattr FD (`hess_en` vs central FD of `grad_en_pulay+grad_en_hellman_feynman`)
   on the **low-symmetry HOF** geometry; require O(h^2), abs < 1e-6. H2O is
   diagnostic only, never primary.
4. Re-validate the `hess_en` 9-block scatter end-to-end on HOF.
5. `pyoqp/oqp/utils/input_checker.py`: reject `[hess] type=analytical` for basis
   with L>=4 (mirror the Fortran abort).
6. Commit + push to `feat/hf-dft-hessian` only when: clean rebuild from source,
   freshly synced `lib/liboqp.so`, oracle executed (not skipped), element-wise
   PySCF pass, contracted FD O(h^2), no stale artifact.

## Build / run recipe (verified)

```
cd build && ninja oqp           # build_rc must be 0
cp build/source/liboqp.so ../lib/liboqp.so   # MUST sync; runtime loads lib/, not build/
# run oracle:
OPENQP_ROOT=$PWD OMP_NUM_THREADS=1 PYTHONPATH=$PWD/pyoqp python3 <driver>
# driver runs an energy calc on HOF then calls oqp.hess1_selftest(runner.mol),
# which writes /tmp/hess_nuc_blocks.txt
```

Gotchas confirmed in this environment:
- Never leave a file named `/tmp/inspect.py` around — it shadows stdlib
  `inspect` and breaks numpy import.
- The runtime loads `lib/liboqp.{so}` (cffi `dlopen`), NOT `build/source`. A
  stale `lib/` copy is the classic "passed but didn't really run the new code"
  trap — always re-sync after `ninja`.
- pytest is not installed; drive via a plain python script + unittest or direct
  calls.
