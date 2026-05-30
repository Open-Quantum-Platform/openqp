# Finite-Difference NAC Oracle — Design Note (for review)

**Status:** proposal, not implemented. No code to be written until this is approved.
**Goal:** a trusted, reproducible finite-difference (FD) nonadiabatic-coupling reference
that becomes the oracle for *any* future analytic NAC implementation (incl. the
`feat/mrsf-analytic-nac-private` rebuild).

---

## 0. Key finding that reframes Phase A

A numerical NAC path **already exists in mainline** (`feat/dft-nmr`, and presumably
`main`/`dev`). We are not building the overlap method from scratch — we are
**validating and hardening** existing code into an oracle. Relevant pieces in
`pyoqp/oqp/library/single_point.py`:

| Component | Location | Role |
|---|---|---|
| `NAC.numerical_nac()` | `single_point.py:1113` | 6N central-difference driver over Cartesian coords; converts `d → h` via energy gap (`:1182–1192`) |
| `nacme_wrapper()` | `single_point.py:1201` | per-displacement worker: SCF → `BasisOverlap` → excitation → `NACME` |
| `NACME.nacme()` | `single_point.py:962` | state overlap → `dc = (S − Sᵀ)/dt` → `h = dc·ΔE` (`:979–986`) |
| `NACME.align_x()` | `single_point.py:934` | amplitude (X-vector) phase/order alignment |
| `BasisOverlap.align_mo()` / `find_vec_order()` | `single_point.py:867`/`903` | MO-level phase/sign + reorder, block-wise (occ / SOMO / virt) |
| `get_structures_ao_overlap` (Fortran) | `source/modules/get_basis_overlap.F90` | cross-geometry AO overlap `⟨χμ(Rₐ)|χν(R_b)⟩`; MO overlap via `Cₐᵀ S C_b` |
| `get_states_overlap` (Fortran) | `source/modules/get_states_overlap.F90` | many-electron state overlap `⟨Ψᵢ(R)|Ψⱼ(R′)⟩` |

**Implication:** Phase A = *characterize + validate + pin conventions* of the above,
plus a clean benchmark driver. We only write new physics code in Phase C.

> ⚠️ **Q1 RESOLVED by trace (see below).** The two FD layers (`nacme` divides by `dt`,
> `numerical_nac` then does `(fwd−bwd)/2`) net to dividing the antisymmetrized
> overlap difference by `2Δ`. Since `S±−S±ᵀ ≈ ±2Δ·d_ij`, the implemented
> **`dcv = 2·d_ij + O(Δ²)`** — a **factor-of-2 overcount**, traceable to `nacme`
> using `(S−Sᵀ)/dt` instead of the HST `(S−Sᵀ)/(2·dt)`. Additionally
> **`nacv = −2·h_standard`** (overall sign from `gap = E_i−E_j` in `numerical_nac`).
> **The oracle magnitude/sign must be corrected before use.** Fix is one line in
> `NACME.nacme` (`/dt → /(2·dt)`) but touches the dynamics path too — decision pending,
> confirm with a single numerical run in Phase A.
>
> Overlap orientation is fixed by `get_states_overlap.F90:915` (`get_dcv` comment):
> `s_st(i,j) = ⟨Ψ_i(old)|Ψ_j(new)⟩`, so the differentiated state is index `j`.

---

## 1. Overlap convention

- **Many-electron state overlap** between displaced geometries is the primitive:
  `S_IJ = ⟨Ψ_I(R) | Ψ_J(R′)⟩`, computed by `get_states_overlap` from the MRSF state
  vectors and the **non-orthogonal cross-geometry AO overlap**
  `S^AO_{μν} = ⟨χ_μ(R) | χ_ν(R′)⟩` (`get_structures_ao_overlap`), transformed to MO
  basis as `S^MO = C(R)ᵀ S^AO C(R′)`.
- **Derivative coupling (DC)** uses the antisymmetric finite difference already in
  `nacme()`:  `d_IJ = (S_IJ − S_JI) / (2 Δ)` — antisymmetrization cancels the
  `O(Δ²)` symmetric (overlap-normalization) part and isolates the coupling.
  *(Code currently divides by `dt`; confirm `dt == 2Δ` vs `Δ` in the wrapper.)*
- **NAC vector (energy-scaled):** `h_IJ = d_IJ · (E_J − E_I)`, matching the analytic
  target `h_IJ = ⟨Ψ_I|∂H/∂R|Ψ_J⟩`. Diagonal `d_II ≡ 0`.
- **Symmetry conventions (to be enforced & asserted) — CORRECTED after the Q1 trace:**
  - `d_IJ = −d_JI` (**derivative coupling: antisymmetric**).
  - `h_IJ = +h_JI` (**NAC vector `h = d·ΔE = ⟨i|∇H|j⟩`: SYMMETRIC**, because `∇H`
    is a Hermitian operator; `h_ij = (E_j−E_i)d_ij` is invariant under i↔j).
  - The private branch is wrong **both** ways: `tdhf_mrsf_nac_backend.F90` made `d`
    symmetric and `h` antisymmetric — exactly the labels swapped. Correct is
    antisymmetric `d`, symmetric `h`. Mainline `numerical_nac` already produces
    symmetric `h` (`dcv`·`gap` = antisym·antisym = symmetric). The analytic rebuild
    must match this.
- **Units:** coordinates in Bohr; energies in Hartree; `d` in 1/Bohr; `h` in Ha/Bohr.

## 2. Phase-tracking strategy

Wavefunction phase/order is arbitrary per geometry and per SCF; without tracking the
FD DC is meaningless. The existing code aligns at **two levels** — keep both:

1. **MO alignment** (`align_mo` + `find_vec_order`): greedy max-|overlap| matching of
   current↔previous MOs with sign fix, done **block-wise** (occupied, the two SOMOs
   `nocc-2:nocc`, virtual) so MRSF's open-shell reference is respected. Optional
   reorder (`nac.align = reorder|sign|no`).
2. **State/amplitude alignment** (`align_x`): greedy max-|overlap| matching + sign fix
   of the MRSF X-amplitude vectors across geometries.

**Protocol additions for the oracle (proposed):**
- Reference phase = the **origin geometry** for every displacement (single anchor),
  rather than chaining, to avoid drift over the 6N points.
- After alignment, **assert** `|S_II| > 0.9` for every tracked state at every
  displacement; if any falls below, the point is rejected (too large `Δ`, or a
  near-crossing) and logged — never silently averaged in.
- Record the applied sign/permutation per point in the benchmark log for
  reproducibility.

## 3. Benchmark molecules

Deliberately boring, far from conical intersections, cheap, reproducible:

| System | State pair | Geometry | Basis | Rationale |
|---|---|---|---|---|
| **H₂** | S₁/S₂ | r ≈ 0.74 Å (near-eq), single stretch coord | 6-31G* | Smallest possible; DC essentially 1-D; trivial to reason about |
| **LiH** | S₁/S₂ | r ≈ 1.6 Å, away from the well-known LiH avoided crossing (~3–4 Å) | 6-31G* | Classic NAC testbed; clean, non-degenerate gap at chosen R |

- MRSF-TDDFT, same functional/settings as the intended analytic target.
- Gaps chosen so `E_J − E_I` is comfortably non-zero (no `1/ΔE` blow-up).
- Geometries, basis, functional, and all `nac`/`tdhf` config frozen in the repo
  (§Phase B) so values are bit-reproducible.

## 4. Displacement protocol

- Central difference, per existing `numerical_nac`: 6N points `R ± Δ eᵢ`.
- **Step study:** run `Δ ∈ {1e-2, 5e-3, 2e-3, 1e-3, 5e-4}` Bohr (mirrors
  `grad_check.py`'s `h = 2e-3` default) to locate the FD plateau.
- Use a **consistent SCF guess** across displacements (the wrapper already passes a
  `guess_file`) so the FD signal is the geometry derivative, not guess noise.
- Tight thresholds: SCF and MRSF response convergence ≥ the planned analytic tol
  (e.g. `1e-8`/`1e-7`) so FD noise floor sits below the acceptance tolerance.
- Single-thread / fixed `OMP_NUM_THREADS` for the reference run to keep values
  deterministic.

## 5. Acceptance criteria (oracle is "trusted" only if all pass)

1. **Phase stability:** with anchored alignment, the DC **direction** (unit vector
   over the 3N components) is invariant to (a) repeated runs and (b) deliberate
   sign/order scrambling of the per-geometry wavefunctions. Cosine similarity to the
   canonical direction `> 1 − 1e-6`.
2. **FD convergence / plateau:** `||d_IJ||` vs `Δ` shows a flat central region; Richardson
   check between `Δ` and `Δ/2` agrees to the target tolerance (e.g. `< 1%`). Establishes
   the recommended production `Δ`.
3. **Finite limit, no `1/Δ` divergence:** as `Δ → 0` the DC **converges to a finite
   value** (confirms we are far from a CI). Explicitly demonstrate the *absence* of
   `1/Δ` growth — `||d_IJ||` must not scale inversely with `Δ`.
4. **Symmetry:** `d` antisymmetric — `||d + dᵀ|| / ||d|| < 1e-6`, diagonal `< 1e-12`.
   `h` **symmetric** — `||h − hᵀ|| / ||h|| < 1e-6` (NOT antisymmetric; see §1).
5. **Energy-scaling consistency:** `h_IJ` reconstructed as `d_IJ·ΔE` matches the
   directly-formed `h` within tolerance (internal cross-check of the `d↔h` path).
6. **Reproducibility:** frozen geometries + config reproduce stored reference
   `d`/`h` to `< 1e-6` on a clean checkout (regression test, Phase B).

Only after 1–6 hold do the stored values become the oracle for Phase C, where the
analytic NAC must match the FD reference in **direction (cosine > 0.999)** and
**magnitude (within a few %)** at the benchmark geometries.

---

## Phase plan (for context — not part of this note's deliverable)

- **Phase A** — validate/harden the existing numerical NAC into a reproducible
  benchmark driver on H₂ + LiH; confirm conventions (§1–2), pass §5.
- **Phase B** — freeze geometries + reference `d`/`h` values; add automated
  regression tests = the permanent oracle.
- **Phase C** — *only then* rebuild analytic machinery: correct interstate 1-PDM
  (oo/vv blocks + Z-vector occ-virt), interstate Z-vector/Lagrangian, antisymmetry
  by formalism (not sign-stamping); validate against the Phase B oracle.

## Open questions for the reviewer

1. **`dt`/`Δ` scaling — RESOLVED (§0):** implemented `dcv = 2·d_ij`, `nacv = −2·h_std`.
   Decision needed: absorb the missing ½ by fixing `NACME.nacme` (`/dt → /(2·dt)`,
   cleanest but affects the MD time-derivative-coupling consumer) **or** divide by 2
   inside `numerical_nac` only (localized to the oracle). Confirm numerically before
   committing. Also pick the global sign convention (`numerical_nac` and `nacme`
   currently disagree on the gap sign).
2. **Independent cross-check:** do you want a *second* independent FD route
   (e.g. direct `(S_IJ−S_JI)/2Δ` between origin and a single displaced geom, bypassing
   the 6N driver) to confirm the production driver — or is the within-method
   convergence study (§5.2) sufficient?
3. **Reference data home:** store frozen geometries + values under `tests/` (pytest
   regression) and/or a `benchmarks/nac/` dir? 
4. **Molecule/state choices:** OK with H₂ and LiH at the proposed non-CI geometries
   and S₁/S₂, or do you want a specific functional/basis to match a literature value?
