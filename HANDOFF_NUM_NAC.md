# Handoff — numerical NAC factor-2 / sign bugfix

**Branch:** `bugfix-num-NAC` (off `origin/main` @ `1890fef`), worktree
`/Users/cheolhochoi/Documents/claude/openqp-numnac`.
**Status:** investigation complete, **no code changed yet**. Proceed only after the
convention in `FD_NAC_CONVENTION_AND_FIX.md` is agreed. This session continues NAC work that
was scoped out of the `gradient-fix` (NMR) session.

## TL;DR
OpenQP's overlap-based numerical NAC produces a **derivative coupling 2× too large**, and the
two code paths use **opposite energy-gap signs**. Both are confirmed analytically and by a LiH
MRSF run (S₁/S₂ pair: `current/hand = 2.0000`, `corrected/hand = 1.0000`). The fix is one line
plus a sign unification, but it also changes the value handed to the **external** NAMD driver
(PyRAI2MD), so the export side is gated on that contract.

## The two bugs (line refs verified on this branch)

`pyoqp/oqp/library/single_point.py`:
- **Factor 2** — `NACME.nacme:979`: `dc_matrix = (state_overlap - state_overlap.T) / self.dt`.
  The antisymmetrized numerator `S−Sᵀ ≈ ±2Δ·d_ij`, so the HST convention needs `/(2·dt)`.
  Missing ½ ⇒ `dc = 2·d_ij`. `NAC.numerical_nac` inherits it (its own `/2` at `:1182` is the
  legitimate central-difference average, not the culprit).
- **Opposite gap signs** — `nacme` (`:980–982`, `e_i=reshape(-1,1)`, `e_j=reshape(1,-1)`) →
  `gap[i,j]=E_j−E_i`; `numerical_nac` (`:1183–1185`, reshapes swapped) → `gap[i,j]=E_i−E_j`.
  So `nacme` gives `+2·h_std`, `numerical_nac` gives `−2·h_std`.
- `nac.dt` default `= 1.0` (`pyoqp/oqp/molecule/oqpdata.py:220`); `numerical_nac` sets `dt=dx`
  per point (`single_point.py:1239`).
- Fortran `get_states_overlap.F90:884–893` (`get_dcv`) returns **raw `F−B`** (no division) and
  its docstring says the denominator "MUST be provided ... can be 2·a for second-order" — i.e.
  the 2 is intended to be supplied by the caller. The Python layer is where it's wrong.

## Correct conventions
`d_ij = ⟨Ψ_i|∇Ψ_j⟩` **antisymmetric**; `h_ij = ⟨Ψ_i|∇H|Ψ_j⟩ = (E_j−E_i)·d_ij` **symmetric**.
(Confirmed numerically: antisymmetrized `dcv` had `‖M+Mᵀ‖/‖M‖ = 0.00`.) The private analytic
branch had both backwards (symmetric `d`, antisymmetric `h`).

## Proposed fix (see FD_NAC_CONVENTION_AND_FIX.md §3)
1. `NACME.nacme:979`: `/ self.dt` → `/ (2.0 * self.dt)` (corrects both oracle and the MD path).
2. Unify gap to `E_j − E_i` everywhere (make `numerical_nac:1183–1185` match `nacme`).
3. Store `d` as canonical, derive `h = (E_j−E_i)·d`.

## NAMD audit result (gating) — see §4a
- **No NAMD propagator inside OpenQP** (`md`/`soc` ∈ `NOT_AVAILABLE_RUNTYPES`). No internal
  consumer compensates for the 2×.
- **Production dynamics is external = PyRAI2MD** (`README.md:15`,
  `github.com/mlcclab/PyRAI2MD-hiam`), fed via file export (`file_utils.py:611`) and the
  `back_door` API (`pyoqp.py:163`).
- ⇒ Internal oracle fix is **safe now**; changing the **exported** time-derivative coupling
  requires checking/coordinating PyRAI2MD's expected denominator first.

## Validation (FD_NAC_CONVENTION_AND_FIX.md §5) — use `fd_nac_verify.py`
LiH MRSF/BHHLYP/6-31G*, single H–z displacement. Targets: S₁/S₂ ratio → 1.0 after fix; **add a
non-degenerate pair**; assert `d` antisym & `h` sym; `nacv == (E_j−E_i)·dcv`; cross-path sign
agreement. Δ-convergence study = supporting evidence only, run afterward.
Build/run recipe: `cd build && ninja`; then
`OPENQP_ROOT=<repo> PYTHONPATH=$OPENQP_ROOT/pyoqp $OPENQP_ROOT/.venv/bin/python fd_nac_verify.py`
(may need to re-point paths to this worktree and re-create/symlink `.venv`).

## Artifacts in this worktree
- `FD_NAC_CONVENTION_AND_FIX.md` — the convention + fix proposal (review this first).
- `FD_NAC_ORACLE_DESIGN.md` — FD-oracle design note (overlap/phase/benchmark/acceptance).
- `fd_nac_verify.py` — the LiH verification script (paths point at the gradient-fix clone;
  update to this worktree before running).

## Open decisions before coding
1. Fix at source (`nacme`) vs localize to `numerical_nac` only (leaves MD at 2×).
2. Confirm PyRAI2MD's denominator convention (external).
3. Approve storing `d` canonical / deriving `h`.

## Not in scope here
The **analytic** MRSF NAC rebuild lives on `karmachoi/openqp-private:feat/mrsf-analytic-nac-private`
(scaffold, not trustworthy — separate effort). This branch is the numerical-oracle bugfix only.
