# NAC Convention & Fix Proposal — OpenQP finite-difference / overlap NAC

**Status:** proposal for agreement. **No code changes are to be made until this is approved.**
**Scope:** the numerical/overlap NAC path (`NACME.nacme`, `NAC.numerical_nac`) and its
consumers (FD oracle, `nacme` time-derivative coupling, NAMD propagator). The analytic
MRSF NAC backend is out of scope here.

**Evidence base:** analytic trace (prior note) + LiH MRSF run. The S₁/S₂ degenerate pair
gives `current/hand = 2.0000` and `corrected/hand = 1.0000` exactly — the factor of two
and its removal are proven, not hypothesized. (`fd_nac_verify.py`, logs in `/tmp/fd_nac_lih/`.)

---

## 1. Current implementation (as written today)

All line numbers are `pyoqp/oqp/library/single_point.py` unless noted.

### 1a. `NACME.nacme()` — single state-overlap pair (also the MD path)

```python
# state_overlap = td_states_overlap, with (Fortran get_dcv comment, get_states_overlap.F90:915):
#   state_overlap[i,j] = <Psi_i(old) | Psi_j(new)>      # column j = "new"/variable geometry
dc_matrix = (state_overlap - state_overlap.T) / self.dt          # :979
e_i = energies[1:].reshape((-1, 1))                              # :980  row index  -> i
e_j = energies[1:].reshape((1, -1))                              # :981  col index  -> j
gap = e_j - e_i                                                  # :982  gap[i,j] = E_j - E_i
nac_matrix = dc_matrix * gap                                     # :983
```
- `self.dt = config['nac']['dt]`, **default = 1.0** (`oqpdata.py:228`).
- In dynamics, `dt` = the MD time step; `state_overlap` is between consecutive MD frames.

### 1b. `NAC.numerical_nac()` — 6N central-difference Cartesian driver

Per displacement, `nacme_wrapper` sets `config['nac']['dt'] = dx` (`:1239`) and runs `nacme`,
producing `dc_matrix = (S±−S±ᵀ)/dx` for each ±Δ point. Then:

```python
forward  = dcm[0:ncoord]            # +Δ points : (S⁺−S⁺ᵀ)/Δ
backward = dcm[ncoord:]             # −Δ points : (S⁻−S⁻ᵀ)/Δ
dcm = (forward - backward) / 2                                   # :1182
e_i = energies[1:].reshape((1, -1))                             # :1183  col index -> i
e_j = energies[1:].reshape((-1, 1))                             # :1184  row index -> j
gap = e_j - e_i                                                 # :1185  gap[i,j] = E_i - E_j
nacm = dcm * gap                                               # :1188
dcv  = dcm.T  -> (nstate,nstate,natom,3)                        # :1191  stored derivative coupling
nacv = nacm.T -> (nstate,nstate,natom,3)                        # :1192  stored NAC vector
```

### 1c. Where the factor of 2 enters

Net assembled derivative coupling (substituting 1b into 1a, `dt = Δ`):
```
dcv[i,j] = { (S⁺−S⁺ᵀ) − (S⁻−S⁻ᵀ) }_{ij} / (2Δ)
```
Taylor (`S±(i,j) = δ_ij ± Δ·d_ij + O(Δ²)`, `d_ij ≡ <Ψ_i|∇Ψ_j>`, `d_ji = −d_ij`):
```
S±−S±ᵀ = ±Δ·(d_ij − d_ji) = ±2Δ·d_ij      ⇒   dcv[i,j] = (2Δ d_ij −(−2Δ d_ij))/(2Δ) = 2·d_ij
```
**The 2× is structural:** `nacme` divides the antisymmetrized overlap by `dt`, where the
Hammes-Schiffer–Tully (HST) convention requires `2·dt`. The missing ½ is the *entire*
discrepancy. (`numerical_nac`'s own `/2` is the legitimate central-difference average and is
not the source.)

### 1d. Where the gap sign enters — and the internal inconsistency

| Path | gap formula | `nac = dc·gap` |
|---|---|---|
| `nacme` (`:982`) | `gap[i,j] = E_j − E_i` | `+ (E_j−E_i)·dc` |
| `numerical_nac` (`:1185`) | `gap[i,j] = E_i − E_j` | `− (E_j−E_i)·dc` |

The two paths use **opposite** energy-gap signs. Relative to the standard
`h_ij = (E_j−E_i)·d_ij`, and with `dc = 2·d_ij`:
- `numerical_nac`: `nacv = 2·d_ij·(E_i−E_j) = −2·h_standard`  (factor 2 **and** sign flip).
- `nacme`:        `nac_matrix = 2·d_ij·(E_j−E_i) = +2·h_standard` (factor 2, standard sign).

---

## 2. Standard convention (target definitions)

For real, orthonormal adiabatic states `Ψ_i(R)`:

```
Derivative coupling:   d_ij = <Ψ_i | ∇_R Ψ_j>                       (vector, per Bohr)
NAC / interstate force: h_ij = <Ψ_i | ∇_R H | Ψ_j> = (E_j − E_i)·d_ij   (Ha/Bohr)
```

Symmetry (both provable, both confirmed numerically in the LiH run):

| object | symmetry | reason |
|---|---|---|
| `d_ij` | **antisymmetric** `d_ij = −d_ji`, `d_ii = 0` | `∇⟨Ψ_i|Ψ_j⟩ = 0` ⇒ `⟨∇Ψ_i|Ψ_j⟩ = −⟨Ψ_i|∇Ψ_j⟩` |
| `h_ij` | **symmetric** `h_ij = +h_ji` | `∇H` Hermitian ⇒ `⟨i|∇H|j⟩=⟨j|∇H|i⟩`; equivalently `(E_j−E_i)d_ij` is invariant under i↔j |

Note this means the private analytic branch is wrong **both** ways (it made `d` symmetric and
`h` antisymmetric — labels swapped); the correct pairing is antisymmetric `d`, symmetric `h`.

Index/orientation convention to standardize on (matches the existing Fortran overlap
orientation, `get_states_overlap.F90:915`): **row = bra state `i`, column = differentiated
ket state `j`**, i.e. stored `d[i,j] = <Ψ_i|∇Ψ_j>`. The overall sign of `d` (and hence `h`)
is gauge-dependent on state phases; phase tracking (`align_x`/`align_mo`) fixes the gauge per
run, and the convention below fixes the *relative* sign of `h` vs `d`.

---

## 3. Proposed fix

### 3.1 Restore the HST ½ at the source (`NACME.nacme:979`)
```
dc_matrix = (state_overlap - state_overlap.T) / self.dt
        →   dc_matrix = (state_overlap - state_overlap.T) / (2.0 * self.dt)
```
This single change yields the **true** `d_ij` from *both* paths:
- `nacme` (MD): `dc = (S−Sᵀ)/(2·dt)` = HST time-derivative coupling (no longer 2×).
- `numerical_nac`: per-point `dc± = (S±−S±ᵀ)/(2Δ) = ±d_ij`, then `(forward−backward)/2 = d_ij`.
  Its existing `/2` is kept (it is the central-difference average, now correct).

### 3.2 Unify the energy-gap convention to `E_j − E_i` everywhere
Make `numerical_nac` (`:1183–1185`) match `nacme` (`:980–982`):
```
e_i = energies[1:].reshape((-1, 1))   # row -> i
e_j = energies[1:].reshape(( 1,-1))   # col -> j
gap = e_j - e_i                        # gap[i,j] = E_j − E_i
```
Result after 3.1 + 3.2: `nacv[i,j] = d_ij·(E_j−E_i) = h_ij^standard` (symmetric, standard
sign) on **both** paths; `dcv[i,j] = d_ij` (antisymmetric). Keep `np.fill_diagonal(gap,1)`
as a divide-by-zero guard (diagonal `d`/`h` are zero regardless).

### 3.3 Canonical internal quantity: store `d`, derive `h`
Recommendation: treat **`d_ij` (antisymmetric) as the canonical stored object**
(`OQP::dc_matrix` / `dcv`) and compute `h_ij = (E_j−E_i)·d_ij` on demand
(`OQP::nac_matrix` / `nacv`). Rationale: `d` is the gauge-fundamental quantity; `h` is a
trivial energy-weighting; storing both independently is exactly what allowed the gap-sign
divergence in §1d. Document the symmetries (`d` antisym, `h` sym) at the storage site and add
cheap assertions (below).

> **Single decision point for the reviewer:** adopt §3.1 (fix at source in `nacme`, corrects
> MD too) **or** the localized alternative — divide by 2 only inside `numerical_nac`
> (`(forward−backward)/2 → /4`) and leave `nacme` at 2×. The localized option avoids touching
> the MD path but *perpetuates* the 2× in `nacme`/dynamics and keeps the two paths physically
> inconsistent. **Recommended: §3.1**, contingent on the NAMD-propagator audit in §4.

---

## 4. Impact analysis

| Consumer | Today | After §3.1+§3.2 | Action required |
|---|---|---|---|
| **FD NAC oracle** (this work) | `dcv = 2·d`, `nacv = −2·h` | `dcv = d`, `nacv = +h` | None beyond the fix; re-baseline reference values |
| **`nacme` time-derivative coupling** (MD per-frame `σ_ij`) | `σ = (S−Sᵀ)/dt` = **2× HST** | `σ = (S−Sᵀ)/(2dt)` = HST | Halves `σ`; see propagator row |
| **NAMD propagator** (electronic TDSE / hopping) | **external (PyRAI2MD)** consumes the 2× `σ` | would receive ½× the previous `σ` | **External to this repo — see §4a.** Coordinate with PyRAI2MD before changing the *exported* `nacme`/`dcme`. |
| **Any saved reference / regression baselines** keyed to current `nacv`/`dcv` | magnitudes 2×, `numerical_nac` sign flipped | corrected | Regenerate baselines; flag in changelog |
| **Branching-plane code** (`compute_bp`, uses `nacv[i,j]`) | uses `−2·h` | uses `+h` | Verify BP geometry/scaling unaffected (overall scale of `h` enters `pitch`/`tilt`) |
| **Analytic MRSF NAC rebuild** (future) | — | must match corrected oracle | Target: `d` antisym, `h` sym, `h=(E_j−E_i)d` |

### 4a. NAMD-propagator audit result (gating item — RESOLVED for the OpenQP side)

**OpenQP contains no internal NAMD propagator.** `md` and `soc` are in
`NOT_AVAILABLE_RUNTYPES` ("recognized but not implemented", `input_checker.py:17,55`); there
is no Verlet/velocity integrator, no electronic-amplitude propagation, no FSSH/Tully hopping
anywhere in `pyoqp` or `source/`. (`population_analysis.F90` is Mulliken/Löwdin charges,
unrelated.) The only internal consumers of `OQP::dc_matrix`/`OQP::nac_matrix` are the producer
(`nacme`) and `numerical_nac` (reads the per-point `dcme` files) — **neither compensates**;
`numerical_nac` simply inherits and propagates the 2× (exactly what the LiH run measured).

**Production dynamics is external — PyRAI2MD** (`README.md:15`,
`github.com/mlcclab/PyRAI2MD-hiam`; `pyoqp/README.md:11` "interface for nonadiabatic molecular
dynamics"). OpenQP is a NAC/gradient **provider**; it hands `σ` to PyRAI2MD via (i) file export
(`np.savetxt` of `nacme`/`dcme`, `file_utils.py:611–620`) and (ii) the `back_door` API
(`pyoqp.py:163`) that injects previous-step data and drives `nacme`.

**Documented design intent confirms the bug location.** The Fortran `get_dcv`
(`get_states_overlap.F90:884–893`) returns the *raw* `F − B` (antisymmetrized overlap, **no
division**) and its docstring states: *"This routine returns d_IJ^a · a = F − B. The
denominator MUST be provided in subsequent calculations because it can be 2·a in the case of
second-order numerical differentiation."* So the Fortran is convention-neutral **by design**,
and the intended denominator for the antisymmetric `F−B` form is **2·a** (= HST `2·dt`). The
Python `NACME.nacme:979` applies only `/dt` — i.e. it contradicts its own documented design
intent. Reference: MRSF-TLF NAC paper, JCTC 2018 (`acs.jctc.8b01049`).

**Conclusion:** No compensation exists anywhere inside OpenQP, so restoring `/(2·dt)` in
`nacme` is unambiguously correct for the FD oracle and *aligns with the documented intent*.
The only open risk is the **external** PyRAI2MD contract: if it was developed against OpenQP's
present 2× export (or compensates internally), changing the *exported* `σ` could shift
dynamics. That code is not in this repo and cannot be verified here.

**Recommended sequencing given the audit:**
- **Now (safe):** apply the ½ on the *internal* oracle path so `numerical_nac`/`nacv`/`dcv` are
  physically correct (`§3.1` effect, scoped to internal NAC). This unblocks the analytic-NAC
  oracle with no external exposure.
- **Before touching the export:** verify PyRAI2MD's expected denominator (read its OpenQP
  interface, or ask mlcclab/lijingbai2009). Then either (a) source-fix `nacme` and
  de-compensate PyRAI2MD in lockstep, or (b) keep the exported time-derivative `σ` legacy and
  document the convention at the boundary.

---

## 5. Validation plan (run after the convention is agreed, before/with the code change)

Reuse `fd_nac_verify.py` (LiH, MRSF/BHHLYP/6-31G*, single H–z displacement).

1. **LiH S₁/S₂ (degenerate pair):** confirm `corrected/hand = 1.0000` (already observed);
   after the code change, `numerical_nac` `dcv` must equal the hand `d` directly (ratio 1.0),
   not 2.0.
2. **One non-degenerate pair:** pick a geometry/state pair with a clean, non-zero gap
   (e.g. S₀→S₁ at a stretched but non-crossing R, or LiH at a geometry lifting the Π
   degeneracy) so `h = (E_j−E_i)d` is non-trivial. Confirm the factor-2 removal there too,
   not only on the degenerate pair. *(Addresses the one contaminated pair seen at single Δ.)*
3. **Symmetry assertions:**
   - `‖d + dᵀ‖ / ‖d‖ < 1e-6` (antisymmetry of `d`) — already 0.00 for the antisymmetrized form.
   - `‖h − hᵀ‖ / ‖h‖ < 1e-6` (symmetry of `h`).
   - `|diag(d)|, |diag(h)| < 1e-12`.
4. **Factor-2 removal / consistency:** corrected `numerical_nac dcv` equals an independent
   hand-built `(S⁺−S⁻)/(2Δ)`; and `nacv == (E_j−E_i)·dcv` elementwise.
5. **Cross-path sign agreement:** `nacme` and `numerical_nac` now report the same sign of `h`
   for the same pair/geometry.

**Δ-convergence study:** to be run *afterward* as supporting evidence (plateau / Richardson),
**not** as a decision input — the factor-2 is already settled by the degenerate-pair identity.

---

## Sign-off checklist (what "agreed" means)
- [ ] Standard definitions in §2 accepted (`d` antisym, `h = (E_j−E_i)d` sym).
- [ ] Fix location chosen: **§3.1 source fix** (recommended) vs localized `numerical_nac`-only.
- [ ] Gap convention unified to `E_j − E_i` (§3.2).
- [ ] Canonical stored object chosen (recommended: store `d`, derive `h`, §3.3).
- [ ] **NAMD propagator audited** for existing 2× compensation (§4) — gating item.
- [ ] Validation plan §5 accepted.

Only after these boxes are checked should code changes begin.
