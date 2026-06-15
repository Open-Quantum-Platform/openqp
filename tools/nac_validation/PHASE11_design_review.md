# NAC Phase 11 — interstate design spec + adversarial review (2026-06-15)

> KEY OUTCOME: the designed 'Route B' (interstate densities via the gradient chain `sfrorhs`)
> is mathematically identical to the polarization route and would REPRODUCE the documented
> ground×ijlr deficiency (Z12 37%, Z13 16.6%). Phase 11 must use the Davidson MATVEC derivative
> (mrsfmntoia/mrsfesum). Step-1 grd2 bilinear is a valid contraction ENGINE; correctness lives
> in the relaxed density fed to it (must come from the matvec, not sfrorhs).

================================================================================
## DESIGN SPEC
================================================================================

I have what I need; the spec is the deliverable, not a task-status change. Here is the implementation spec.

---

# NAC Phase 11 — Implementation Spec (Steps 2–4: interstate machinery)

## Orienting reality check (read first — it changes the plan)

Two distinct numerator routes exist in this tree, and the reader reports describe the **explicit-interstate** route while the production driver already runs the **polarization** route. The spec below specs the explicit route (as requested) but flags where it is redundant with code that already works, so you don't rebuild validated machinery.

- **Route A — pseudo-state polarization (already wired, already gives the full G_IJ numerator).** `analytical_nac` (`pyoqp/oqp/library/single_point.py:1747-1796`) computes `h_ci = 1/2[G(X_+) - G(X_-)]` with `X_± = (X_I±X_J)/√2`, where `G` is the *unmodified* `tdhf_mrsf_z_vector` + `tdhf_mrsf_gradient` path. Because `(X_I+X_J)⊗(X_I+X_J) - (X_I-X_J)⊗(X_I-X_J) = 2(X_I⊗X_J + X_J⊗X_I)`, every quadratic-in-X object (tij/tab, hx, rhs, xk, td_p, wao, the 7 channels' bilinear products inside grd2) is automatically replaced by its symmetrized interstate bilinear. This is the **full** numerator G_IJ = X_I^T(d_x A)X_J + relaxation, with the SCF/ground term cancelling in the difference. No interstate density code is required for the numerator at all.
- **Route B — explicit interstate densities (what Steps 2–4 below build).** Build `spcI`/`spcJ`, `tij^{IJ}/tab^{IJ}`, `Z^{IJ}`, `p2^{IJ}`, `W^{IJ}` directly and call the bilinear `grd2_mrsf_nac_compute_data_t` (Step 1, done) + a 1e/W mirror in one shot. This avoids the two-gradient subtraction (cheaper, no `1/2[G_+ - G_-]` cancellation noise) and is the clean analytic form for the paper.

**Decision the implementer must make explicitly (RISK-0):** Route A already produces the validated numerator. Route B's *only* deliverables of new value are (i) a single-pass G_IJ that avoids polarization round-off, and (ii) the explicit `p2^{IJ}`/`W^{IJ}` objects. Build Route B as an *independent verifier* of Route A, gated behind a flag, and treat **Route A's `h_ci` as the production numerator**. The I=J self-test in Step 1 (`mrsf_nac_amp_selftest`, `tdhf_mrsf_gradient.F90:928-1035`, asserts `max|de_nac(I=J) − de_prod| ≈ 0`) is the diagonal anchor for Route B's 2e block; Steps 2–4 extend that anchor to 1e/W and to the off-diagonal.

---

## Step 2 — interstate transition + relaxed densities

### 2a. Build `spcI` and `spcJ` (the 7 channels for states I and J)

Per-state amplitude slice is column `bvec_mo(:,state)` of `OQP_td_bvec_mo`, shape `(xvec_dim, nstate)`, `xvec_dim = nocca*nvirb` (`tdhf_mrsf_z_vector.F90:1051,979`). The channel builder is `mrsfcbc` (`tdhf_mrsf_lib.F90:490`), fed a *square* MO-amplitude matrix from `iatogen` (`tdhf_lib.F90:381-399`). Production wiring to replicate: `tdhf_mrsf_z_vector.F90:1675-1684`.

Recipe (two private buffers, do **not** reuse `td_mrsf_den` — the CPHF override zeros it at `:1818`):

```fortran
real(kind=dp), allocatable :: spcI(:,:,:), spcJ(:,:,:)   ! (7, nbf, nbf)
real(kind=dp), allocatable :: wsq(:,:)                    ! (nbf, nbf) scratch
allocate(spcI(7,nbf,nbf), spcJ(7,nbf,nbf), wsq(nbf,nbf), source=0.0_dp)

call iatogen(bvec_mo(:,I), wsq, nocca, noccb)
call mrsfcbc(infos, mo_a, mo_a, wsq, spcI)      ! va=vb=mo_a (shared MRSF ref)

call iatogen(bvec_mo(:,J), wsq, nocca, noccb)
call mrsfcbc(infos, mo_a, mo_a, wsq, spcJ)
```

Channel order is fixed by `tdhf_mrsf_lib.F90:513-519`: `spc(1..7) = {bo2v,bo1v,bco1,bco2,o21v,co12,ball}`. Each channel is **linear** in its state's amplitude; the bilinear cross-product `spcI×spcJ` is the correct order-2 object, and `grd2_mrsf_nac_compute_data_t_get_density` (`tdhf_mrsf_gradient.F90:719-918`) already symmetrizes it as `1/2[f_I g_J + f_J g_I]`.

**RISK-2a-1 (channel-7 ball, from report open-question):** production overwrites `fmrst1(1,7,:,:)` with `td_abxc` from `sfdmat` (`:1682`), noting `ball == td_abxc` to 1e-15 *for the target state only*. For arbitrary I, use `mrsfcbc`'s own `spcI(7,:,:)` directly (the algebra is state-independent — `ball` is assembled purely from `bvec` rows/cols at `tdhf_mrsf_lib.F90:727-794`). Add a one-time assertion in the I=J self-test extension: build `spc` via `mrsfcbc` and via the `sfdmat` path and assert `max|spc(7)_cbc − td_abxc| < 1e-12`. If it holds at I=J it holds for all I by the same algebra. **This is low risk but must be asserted once.**

**RISK-2a-2:** `mrsfcbc` reads `mrst` from `infos%tddft%mult` internally (singlet/triplet sign pattern at `:780-794`). Confirm `infos%tddft%mult ∈ {1,3}` before the bilinear path; mrst==5 (quintet) has no interstate tden (`mrsf_interstate_tden` aborts, `tdhf_mrsf_lib.F90:1737-1738`) — **Phase-11 off-diagonal is singlet/triplet only.**

### 2b. Interstate UNRELAXED difference density → `p2`'s unrelaxed part

Call `mrsf_interstate_tden(infos, bvec_mo, I, J, tij, tab)` (`tdhf_mrsf_lib.F90:1716-1770`). It is already general:
- `tij(noca,noca) = −1/2(Xi Xj^T + Xj Xi^T)` (occ-occ, alpha block)
- `tab(nvirb,nvirb) = +1/2(Xi^T Xj + Xj^T Xi)` (vir-vir, beta block)

Both factors run through `mrsfxvec` (`tdhf_mrsf_lib.F90:1643-1693`) so the SOMO O→O folding (`ijlr1/ijlr2`, 1/√2) is consistent.

**Mapping into the `(nbf,nbf,2)` p2 layout the bilinear `get_density` expects.** The production unrelaxed→AO→pack path is `sfropcal` (with `z=0`) → `orthogonal_transform` → `symmetrize` → `pack` → `×0.5` (`tdhf_mrsf_z_vector.F90:1493-1530`, `tdhf_sf_lib.F90:704-763`). The relaxed density `td_p` is `(nbf_tri, 2)` packed, **alpha in slot 1, beta in slot 2** (`:1036`). When read by the gradient driver it is unpacked into `p(nbf,nbf,2)` (`tdhf_mrsf_gradient.F90:164-165`), and the bilinear type's `init` (`:700-701`) folds α/β → tot/spin **identically** to the production type (`:478-479`):

```fortran
p2(:,:,1) = p2(:,:,1) +   p2(:,:,2)      ! total
p2(:,:,2) = p2(:,:,1) - 2*p2(:,:,2)      ! spin
```

So the **unrelaxed part of `p2^{IJ}` enters through the same `td_p` tagarray** — you do not hand a separate p2 to `grd2`; you let the z-vector pipeline write `td_p` as `0.5*(T^{IJ}_unrelaxed + relaxation(Z^{IJ}))` and the gradient reads it. The α/β→tot/spin fold is the same matrix algebra regardless of I=J vs I≠J, so **no new fold logic is needed**; you only need `td_p` to *contain* the interstate object.

**I=J reduction (exact):** at I=J, `mrsf_interstate_tden` takes `xj=xi` (`:1744-1745`) → `tij=−Xi Xi^T`, `tab=+Xi^T Xi`, **bit-identical** to the production unrelaxed difference density (the driver itself calls it with `target_state,target_state` at `:1735`). Confirmed exact by construction.

---

## Step 3 — interstate relaxed density via the z-vector

### 3a. Interstate RHS — keep the gradient (`sfrorhs`) branch, NOT the CPHF override

The production RHS is `sfrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, Tij, Tab, Fa, Fb, ...)` (`tdhf_sf_lib.F90:345-421`), assembled at `tdhf_mrsf_z_vector.F90:1758`. It is **state-index-free** — correct for I≠J the moment its inputs (`tij/tab`, `hxa/hxb`, `ab1_mo_*`) are the interstate versions. The interstate RHS formula is exactly `sfrorhs` with:

| input | production (I=J) source | interstate (I≠J) source |
|---|---|---|
| `tij/tab` | `mrsf_interstate_tden(...,target,target,...)` `:1735` | `mrsf_interstate_tden(...,I,J,...)` |
| `pa` (for `ab1_mo_*`) | `unpack(ta/tb)` `:1609-1610`, ta/tb from `sfdmat(target)` `:1080,1082` | `unpack(tij/tab^{IJ})` — i.e. feed the *interstate* T into the apb response `int2_tdgrd_data_t(d2=pa)` `:1639` |
| `hxa/hxb` | `mrsfcbc(target)→fmrst2→mrsfsp` `:1675-1727` | **see RISK-3a-1 below** |

**The two `sfrorhs`/`sfropcal` consumers of tij/tab are state-free; only the *producers* change.** Replace the `:1735` call's args with `(I,J)`, and feed `pa` (the apb-response density at `:1609-1610` / `:1639`) from the interstate `tij/tab` instead of `sfdmat`'s `ta/tb`.

**RISK-3a-1 (the central open question — hxa/hxb interstate factorization).** `hxa/hxb` (`:1707-1727`) are built from `fmrst2` (the ERI response of state-`target` channels via `mrsfcbc`→`int2_data_st`, `:1676,1687-1694`) contracted with `iatogen(mrsfxvec(target))` (`:1709-1716`), then `mrsfsp` (`:1727`). This is **bilinear in the single target X** (one X in `fmrst2`'s channel build, one in `wrk3`). For I≠J the bilinear must pair **one factor from I and one from J**, symmetrized. Concretely the interstate `hx^{IJ}` should be `1/2[ h(spcI, X_J) + h(spcJ, X_I) ]`:
- build `fmrst2` from `spcI` (i.e. `int2_data_st` fed with `spcI` instead of `fmrst1` from target) and contract with `wrk3=iatogen(mrsfxvec(J))`,
- plus the transpose term `fmrst2` from `spcJ` contracted with `iatogen(mrsfxvec(I))`,
- each `mrsfsp` call (`:1727`) similarly fed the cross channel set.

This is **not yet implemented or proven** in the tree. **The cleanest validation is Route A:** the polarization path generates exactly this symmetrized `hx^{IJ}` automatically (since `hx(X_+) − hx(X_−)` is bilinear). So: implement the explicit `hx^{IJ}` per the formula above, then verify the resulting `G_IJ` against Route A's `h_ci`. If they disagree, the `hx` symmetrization is the prime suspect. **Do not assume `hx` factorizes as "only tij/tab need the J index" — it does not; `hx` is an independent bilinear that needs its own I/J symmetrization.**

**I=J reduction (exact):** with I=J, `mrsf_interstate_tden` → production tij/tab; `spcI=spcJ=spc(target)` → `hx^{IJ}=h(spc,X)` = production hx; `pa` from interstate T = production T. So `sfrorhs` returns the production RHS bit-for-bit. Verified structurally because the diagonal calls are literally the production calls.

**Do NOT take the CPHF override branch** (`:1769-1823`). That branch (`mrsf_nac_cphf_mode`) zeros `rhs` and rebuilds it from the *antisymmetrized* `gamma^IJ = wrk1(p,q)−wrk1(q,p)` via `get_mrsf_transition_density` (`:1782`). That is the **orbital-overlap (Pulay/U^x) RHS for `mrsf_nac_overlap`**, a *different physical contribution* (the `d^ov` term in `analytical_nac:1799`), NOT the amplitude relaxed-difference density. The amplitude path needs the **symmetric** interstate T + nonzero hx; the override needs the **antisymmetric** gamma + zeroed T/hx. **These are mutually exclusive modes** — gate them with separate flags (see Code Plan, RISK-CP-1).

### 3b. Solve and assemble `p2^{IJ} = T^{IJ}_unrelaxed + relaxation(Z^{IJ})`

The LHS operator `apply_z_operator` (`:747-803`), preconditioner `apply_z_precond`/`sfromcal` (`:811-823,1120`), and all solvers (`:1127-1136`) are **state-independent** (built from orbital energies + fa/fb only) and reused **verbatim**. `(A+B) Z^{IJ} = rhs^{IJ}` with the interstate RHS from 3a.

Assembly via `sfropcal(wrk1, wrk2, tij, tab, xk, nocca, noccb)` (`:1496`, def `tdhf_sf_lib.F90:704-763`): copies `T^{IJ}` into the occ-occ/vir-vir blocks, then adds `Z^{IJ}` with weight 0.5 over the three ROHF rotation blocks. Then AO-transform (`:1505-1508`), `symmetrize_matrix` (`:1525-1526`), `pack` into `td_p` (`:1527-1528`), global `×0.5` (`:1530`). All reused unchanged — `sfropcal` is state-free.

**RISK-3b-1 (td_p has only 2 spin slots, report open-question).** `OQP_td_p` shape `(nbf_tri,2)` (`:1036`). The interstate `p2^{IJ}` from `mrsf_interstate_tden` is **already symmetric** (the 1/2 transpose-sum at `tdhf_mrsf_lib.F90:1750-1768`), so the `symmetrize_matrix` (`:1525-1526`) and `0.5` packing (`:1530`) preserve normalization exactly as in the diagonal case. No extra slot needed. Confirm by the I=J self-test: the assembled `td_p` must equal the production `td_p` (already true since the diagonal call is the production call).

### The differencing doctrine (-1 prefactor / zero-densities)

This is **only** for the CPHF/overlap override path, NOT the amplitude path. In `mrsf_nac_cphf_mode` the code zeros `td_abxc, td_mrsf_den, tij, tab, hxa, hxb` (`:1817-1822`) so the gradient contracts *only* the orbital-response Z, and the ground/SCF part is removed by the Python-side differencing (the `1/2[G_+−G_−]` cancellation, or the explicit subtraction of the SCF gradient). **For the amplitude (Route B) path, do the OPPOSITE: keep tij/tab/hx nonzero, and let the α/β→tot/spin fold + the bilinear grd2 handle the contraction.** The SCF/ground cancellation in Route B happens because `d2` (reference density) is state-independent and the bilinear `get_density` only ever multiplies it against the *interstate* `p2`/channels (the pure `d2×d2` term is absent from `get_density` — check: `grd2_mrsf_nac_compute_data_t_get_density:787` has `d2×d2` only inside the `(d2+p2)×d2` Coulomb structure, identical to production, so the ground part is handled by the *production* gradient being subtracted, i.e. Route A, or by `G_IJ` being inherently a difference). **RISK-3b-2: in single-pass Route B you must still subtract the pure-SCF gradient** (the `d2`-only piece), because `grd2_mrsf_nac_compute_data_t` includes `d2×d2`-containing terms via `(d2+p2)`. The polarization identity (Route A) removes it automatically; single-pass Route B needs an explicit ground-gradient subtraction or a verification that the SCF part is X-independent and cancels in the J-axis antisymmetry. **Flag this — it is the subtlest correctness point.**

---

## Step 4 — the 1e/W (mrsfesum-derived) term

`sf_1e_grad` (`tdhf_sf_gradient.F90:179-267`) is the production 1e/W routine, called at `tdhf_mrsf_gradient.F90:143`. It reads `OQP_td_p` (`p`) and `OQP_WAO` (`w`) tagarrays and does:

1. **W·S^[x] term:** `eijden(dens)` (SCF Lagrangian, `grd1.F90:114-138`) → `dens = dens + 2*w` (`:242`) → `grad_ee_overlap(basis, dens, grad)` (`:245`).
2. **1e term:** `dens = dmat_a + p(:,1) (+ dmat_b + p(:,2))` (`:248-251`) → `grad_en_hellman_feynman` + `grad_ee_kinetic` + `grad_en_pulay` + `grad_1e_ecp` (`:254-263`).

**The interstate 1e/W term needs NO new routine** — it is `sf_1e_grad` reading interstate `td_p` (= `p2^{IJ}` from Step 3b) and interstate `wao` (= `W^{IJ}`). The density object for the 1e contractions is exactly `dmat_a + p^{IJ}(:,1) (+β)`; the W object for `grad_ee_overlap` is `eijden + 2·W^{IJ}`. **Both are written by the z-vector step before the gradient runs.** So Step 4 reduces to: make `build_mrsf_relaxed_density_and_w` (`:1493-1587`) write the interstate `td_p` and `wao`, then call the *stock* `sf_1e_grad` and `mrsf_2e_grad`.

`W^{IJ}` is built by `mrsfrowcal(wmo, mo_energy_a, fa, fb, xk, hxa, hxb, ppija, ppijb, ...)` (`:1570`, def `tdhf_mrsf_lib.F90:2242-2383`), AO-transformed, packed, then **quartered** (`×0.5×0.5` at `:1585,1587` — gradient 1/2 × ROHF 1/2). This is state-free in its formula: feed it the interstate `xk=Z^{IJ}`, interstate `hxa/hxb=hx^{IJ}`, and `ppija/ppijb` from the interstate relaxed-density apb response (`:1546-1563`, also state-free given interstate `pa`). Reused unchanged.

**I=J reduction (exact):** at I=J, `Z^{IJ}→Z`, `hx^{IJ}→hx`, `p2^{IJ}→td_p`, `W^{IJ}→WAO`. Then `sf_1e_grad`'s `eijden+2·WAO` and `dmat+td_p` are the production inputs → production 1e/W gradient. This is the exact 1e/W analogue of the passing 2e self-test (`mrsf_nac_amp_selftest:1010-1031`).

**RISK-4-1 (SOMO corrections, report open-question).** `mrsfesum`'s SOMO special-cases (`lr1=nocca-1, lr2=nocca`, `xlr/√2` terms, `tdhf_mrsf_lib.F90:1474-1544`) live in the **Davidson matvec**, not the gradient. The question is whether they are already absorbed into `WAO`/`td_p` before the gradient stage. They are: `W` is built by `mrsfrowcal` (which carries its own SOMO handling via `fa/fb/xk/hx`), and the SOMO folding in the *amplitude* enters through `mrsfxvec` inside `mrsf_interstate_tden` and inside the `hx` build. So the 1e/W gradient term assembles purely from `td_p^{IJ} + W^{IJ}` **provided `hx^{IJ}` (RISK-3a-1) carries the correct SOMO pairing**. The SOMO risk is therefore *entirely localized in the hx interstate build* (Step 3a), not in Step 4. **Verify via the I=J self-test extended to 1e/W (see Code Plan).**

**RISK-4-2 (mult=1 vs 3 SOMO signs).** `mrsfesum`/`mrsfrowcal`/`mrsfcbc` all branch on `mrst` with different SOMO sign patterns (singlet `tdhf_mrsf_lib.F90:780-787` vs triplet `:788-794`; `sgnk=±1` in `get_density:534,766`). Since the interstate build reuses these same routines with the same `mrst`, the signs propagate correctly — **but the bilinear cross-term must use the same `mrst` for both spcI and spcJ** (both states share the MRSF reference and multiplicity). Assert `infos%tddft%mult` is the single mult for the pair.

---

## Code plan

### Files and routines

**`tdhf_mrsf_z_vector.F90`:**

1. Add module flag `mrsf_nac_p2_mode` + states `mrsf_nac_p2_istate/jstate` (mirror `:27-57`), with a C setter `set_mrsf_nac_amplitude` bind(C). This selects the **amplitude (Route B)** path, **distinct from** `mrsf_nac_cphf_mode` (RISK-CP-1: the two NAC modes are mutually exclusive — assert not both set).
2. In `build_mrsf_zvector_rhs` (`:1594`): when `mrsf_nac_p2_mode`, replace the `mrsf_interstate_tden(...,target,target,...)` call (`:1735`) with `(...,p2_istate,p2_jstate,...)`; build `spcI/spcJ` per Step 2a into `fmrst1`-shaped buffers and feed the apb response `pa` (`:1609-1610`) and the `hx^{IJ}` build (`:1707-1727`) with the **symmetrized interstate** versions per RISK-3a-1. Take the `sfrorhs` branch (`:1758`), **skip** the CPHF override (`:1769-1823`).
3. `build_mrsf_relaxed_density_and_w` (`:1493`): unchanged — it consumes interstate `tij/tab/xk/hx` transparently and writes interstate `td_p`/`wao`.

**`tdhf_mrsf_gradient.F90`:**

4. Add a new C-bound entry `mrsf_nac_amplitude` bind(C, name="mrsf_nac_amplitude") that runs the **single-pass Route B**: sets `mrsf_nac_p2_mode`, calls `tdhf_mrsf_z_vector`, then runs `mrsf_2e_grad` with `grd2_mrsf_nac_compute_data_t(spcI,spcJ,...)` (Step 1 type, `:36-48`) instead of `grd2_mrsf_compute_data_t`, plus the stock `sf_1e_grad`. Subtract the SCF/ground gradient per RISK-3b-2.
5. Extend `mrsf_nac_amp_selftest` (`:928-1035`) to also run `sf_1e_grad` and compare the **full numerator** (1e+2e+W) at I=J against the production excited-state gradient — this is the **G_22 == production gradient** check.

**`pyoqp/oqp/library/single_point.py`:** `analytical_nac` (`:1695`) already implements Route A. Add an optional `oqp.mrsf_nac_amplitude(mol)` single-pass path gated by config, cross-checked against `h_ci`.

### Validation ladder (each sub-step has a hard gate)

1. **Step 2 (spcI build):** at I=J, assert `mrsfcbc(bvec_mo(:,target))` channels == production `td_mrsf_den` to 1e-13, and channel-7 == `td_abxc` to 1e-12 (RISK-2a-1).
2. **Step 3 (z-vector RHS/p2):** at I=J, assert interstate `rhs == sfrorhs` production rhs and assembled `td_p == production td_p` to 1e-12 (structural — same calls).
3. **Step 4 (1e/W):** extend `mrsf_nac_amp_selftest` so the **full** G_22 (1e+2e+W) at I=J equals the stock `tdhf_mrsf_gradient` output to ~1e-12 (the headline diagonal anchor).
4. **Off-diagonal (I≠J):** assert Route B `G_IJ` == Route A `h_ci = 1/2[G_+−G_−]` to ~1e-9 per atom/component, for at least one HF and one BHHLYP molecule. **This cross-check is the off-diagonal correctness gate** — Route A is the trusted reference.
5. **Benchmark npz:** run the existing analytical-vs-numerical NAC comparison (the `numerical_nac` path at `:1810`, `TLF(1)`, `dx=1e-4`) and store `(nstate,nstate,natom,3)` arrays; target the documented 0.00–0.08% agreement.

### Consolidated risk register

- **RISK-0:** Route A already gives the production numerator; build Route B as a single-pass verifier, not a replacement. Don't rebuild validated machinery.
- **RISK-3a-1 (highest):** `hx^{IJ}` is an independent bilinear needing its own I/J symmetrization `1/2[h(spcI,X_J)+h(spcJ,X_I)]`; the reports' guess that "only tij/tab need the J index" is **wrong**. Validate against Route A.
- **RISK-3b-2:** single-pass Route B must explicitly subtract the X-independent SCF/ground gradient (the `(d2+p2)·d2` Coulomb/`dq1` pieces in `get_density` contain `d2·d2`); Route A removes it automatically via polarization. Prove cancellation or subtract.
- **RISK-CP-1:** `mrsf_nac_cphf_mode` (overlap/U^x, antisymmetric gamma, zeroed T/hx) and `mrsf_nac_p2_mode` (amplitude, symmetric T, nonzero hx) are mutually exclusive — separate flags, assert not-both.
- **RISK-2a-1:** channel-7 `ball` from `mrsfcbc` assumed == `sfdmat` `td_abxc` for arbitrary I — assert once at I=J.
- **RISK-4-1/4-2:** SOMO corrections + mult=1/3 signs are correct iff `hx^{IJ}` (Step 3a) carries the right SOMO pairing; localized entirely in the hx build.
- **RISK-quintet:** off-diagonal is singlet/triplet (mult 1,3) only; `mrsf_interstate_tden` aborts for mult=5 (`tdhf_mrsf_lib.F90:1737-1738`).

### Key file:line anchors
- Bilinear 2e type (done): `tdhf_mrsf_gradient.F90:36-48, 689-918`; I=J self-test `:928-1035`.
- Interstate tden: `tdhf_mrsf_lib.F90:1716-1770` (I=J reduces exact, `:1744-1745`).
- 7-channel build: `mrsfcbc` `tdhf_mrsf_lib.F90:490-811`; `iatogen` `tdhf_lib.F90:381-399`; driver `tdhf_mrsf_z_vector.F90:1675-1684`.
- RHS: `sfrorhs` `tdhf_sf_lib.F90:345-421`; driver call `tdhf_mrsf_z_vector.F90:1758`; CPHF override (do NOT use) `:1769-1823`.
- Relaxed density: `sfropcal` `tdhf_sf_lib.F90:704-763`; driver `tdhf_mrsf_z_vector.F90:1493-1530`. W: `mrsfrowcal` `tdhf_mrsf_lib.F90:2242-2383`; driver `:1565-1587`.
- 1e/W gradient: `sf_1e_grad` `tdhf_sf_gradient.F90:179-267`; `eijden` `grd1.F90:114-138`.
- Polarization numerator (Route A, working): `single_point.py:1695-1808`.
- C setters pattern: `tdhf_mrsf_z_vector.F90:39-57`.

================================================================================
## ADVERSARIAL REVIEW (verdict)
================================================================================

I now have everything I need to write the adversarial review. The key findings are grounded in the actual source. Let me verify one final critical point: the spec's claim that `eijden` provides a state-independent SCF Lagrangian (which would mean the W·S^x term in single-pass Route B contains an SCF piece that must be retained, not subtracted) — and confirm `mrsfsp` is fed `wrk3` from the single target state (the bilinear-in-one-X structure for hx).

The `eijden` reads `OQP_FOCK_A/B` and `OQP_DM_A/B` — confirmed **state-independent** SCF Lagrangian. This is important for RISK-3b-2: the 1e/W path's `eijden` term is pure-SCF and is the same in production and NAC, so in single-pass Route B it does NOT cancel by itself.

I have enough to deliver the review.

---

# Adversarial review — NAC Phase 11 Steps 2–4 spec

I read the spec against the actual tree at `~/openqp-nac` @ `9f6c2c4` (the branch the spec targets). Findings are grounded in source. Ranked by severity.

## RISK-A (BLOCKER) — The whole premise contradicts a documented, proven-dead route. Route B ≡ Route A *exactly*, so it cannot rescue the one pair Phase 11 exists to fix.

**Why.** Your own project memory (`nac-project.md`, "RETIRED ROUTES — DO NOT REVISIT") records: *"Hand-symmetrized bilinear z-vector mode ≡ polarization identity exactly (bilinear of a quadratic form is unique; parallelogram law proved the chain quadratic to 1e-4) — would reproduce the deficiency."* Route B is *precisely* that hand-symmetrized bilinear mode. The spec even concedes this in RISK-3a-1: it validates `G_IJ` against Route A's `h_ci` and says "if they disagree, the hx symmetrization is the prime suspect." But a correct Route B must **agree** with Route A — they are two evaluations of the same unique bilinear form of the same quadratic chain `G(X)`. The spec's validation ladder gate #4 ("Route B == Route A to 1e-9") is therefore a *tautology check*, not a correctness check: passing it proves only that you re-derived Route A, and Route A is already known to be **structurally deficient for the ground×ijlr (S0↔open-shell-singlet) pair** — the one pair that motivated Phase 11 (the key NAMD pair, per memory). For the pairs where Route A works (non-ground, ground×generic), Route B adds nothing but round-off avoidance.

**The spec has the value proposition inverted.** It says (RISK-0) "Route A already gives the production numerator; build Route B as a verifier." But the actual open problem is the ground×ijlr pair where Route A's *quadratic chain itself* is deficient — and Route B inherits that exact chain (same `mrsfcbc`, same `int2_data_st`, same `mrsfsp`, same `sfrorhs`, same `apply_z_operator`). Memory is explicit that Phase 11's fix must come from the **Davidson matvec derivative** (`X_I·mrsfesum(X_J)` and `mrsfmntoia` channel pairing), which re-derives the *A the states actually diagonalize*, so "the class-mixture deficiency cannot arise." Steps 2–4 as written rebuild the gradient-side chain, not the matvec-side `A`.

**Decisive test.** Before implementing anything: take the existing `(1, ijlr-state)` pair where Route A is known deficient, and hand-evaluate the proposed Route B `G_IJ` for *that* pair (you can do it numerically via the polarization identity, which the spec admits equals Route B). If `d_amp = G_IJ/ΔΩ` still mismatches `d_num` for ground×ijlr — which the parallelogram-law proof says it will — Route B is dead on arrival and Steps 2–4 must be re-scoped onto the matvec derivative. **Run this test first; it is one displaced-amplitude gradient pair and gates the entire 3-step effort.**

**Resolution.** Either (a) re-scope Steps 2–4 to differentiate the *Davidson matvec* (`mrsfesum`, `mrsfmntoia` pairing) per the memory's Phase-11 spec, not the gradient chain; or (b) explicitly accept Route B as "single-pass verifier for the pairs Route A already solves" and state in the spec that it does NOT address ground×ijlr — but then Phase 11's headline deliverable is unmet.

## RISK-B (HIGH) — `hx^{IJ}` symmetrization (RISK-3a-1) is under-specified and the spec's own formula omits a channel.

**Why.** I confirmed `hxa/hxb` at `tdhf_mrsf_z_vector.F90:1707-1727` are built from *two* X-dependent objects: (i) `fmrst2` = `int2_data_st%f3`, whose ERI response is built from `fmrst1` = `mrsfcbc(target)` channels (one X), contracted at `:1713-1720` with `wrk3 = iatogen(mrsfxvec(target))` (a second X); **and** (ii) `mrsfsp(hxa,hxb, mo_a, mo_a, wrk3, fmrst2,...)` at `:1727`, which contracts `fmrst2` channels *again* with `xv=wrk3` (a third appearance of the target X — see `tdhf_mrsf_lib.F90:1820` `scr(:,j)*xv(lr1,j)` etc.). So `hx` is bilinear in X across **two independent factors**: the channel set (in `fmrst2`, ultimately from `spc`) and the `wrk3`/`xv` amplitude. The spec's formula `½[h(spcI,X_J)+h(spcJ,X_I)]` is the right *idea*, but:
- It must symmetrize **both** the dgemm contraction (`:1713-1720`, `wrk3`) **and** the `mrsfsp` call (`:1727`, `xv=wrk3`) — the spec mentions `mrsfsp` "similarly fed the cross channel set" but `mrsfsp` takes `fmrst2` (channels) *and* `xv` (amplitude) as **separate** arguments; you must pass `fmrst2(spcI)` with `xv=mrsfxvec(J)` and add the transpose. The spec does not pin this down and it is the single most error-prone edit.
- `fmrst2` is the **ERI response of the channels**, not the channels themselves. Feeding `int2_data_st` with `spcI` (per spec Step 2a/3a) gives `fmrst2(spcI)`; but note `:1697` flips sign for triplet and `:1700-1705` rescales channels 5/6/1-4 by `spc_*` ratios. These rescalings must be applied identically in the interstate build or the I≠J `hx` is silently wrong while I=J stays exact (the rescale is X-order-independent, so the self-test won't catch a missing rescale if you accidentally drop it only in one of the two cross terms).

**Resolution.** Spec the `hx^{IJ}` build as a literal two-term sum, each term a *full* replay of `:1687-1727` (int2 run on `spcX` → sign/rescale → dgemm with `mrsfxvec(Y)` → `mrsfsp(fmrst2(spcX), xv=mrsfxvec(Y))`), then `½(term_IJ + term_JI)`. Add an I≠J unit test that perturbs `spc_ovov` away from `hfscale` to exercise the `:1702` rescale path (otherwise the self-test runs with all ratios = 1 and never tests it).

## RISK-C (HIGH) — Single-pass Route B does NOT avoid the SCF-reference cancellation; the spec's RISK-3b-2 is correct but its "maybe it cancels by antisymmetry" hedge is wrong.

**Why.** I confirmed `grd2_mrsf_nac_compute_data_t_get_density` keeps the pure reference term: `df1 = (d2+p2)·d2 + d2·p2` at `:787-788` contains `d2·d2` (pure SCF Coulomb), and `dq1` at `:792-799` contains `d2·d2` exchange. These are **X-independent** and present in every call. Likewise `sf_1e_grad`'s `eijden` (grd1.F90, reads `OQP_FOCK_A/B`+`OQP_DM_A/B` — state-independent) produces the **full SCF Lagrangian**, and `:248 dens = dmat_a + p(:,1)` always includes the bare SCF density `dmat_a`. So single-pass Route B's raw output = `G_SCF + G_amp^{IJ}`, **not** `G_amp^{IJ}`. You must subtract the pure SCF gradient.

The spec floats "or a verification that the SCF part is X-independent and cancels in the J-axis antisymmetry." **It does not cancel.** `d_IJ` is antisymmetric in I↔J, but a *single-pass* evaluation of `G_IJ` at fixed (I,J) produces `G_SCF` as an additive constant offset — antisymmetry would only remove it if you also computed `G_JI` and subtracted, which is a *second* pass (and `G_IJ` is already symmetric in the bilinear, so `G_JI=G_IJ` and the SCF part doesn't cancel at all). Route A removes `G_SCF` because `½[G(X₊)−G(X₋)]` differences two gradients with identical `G_SCF`. Single-pass Route B has nothing to difference against.

**Resolution.** Mandate an explicit ground-state SCF gradient subtraction in the new `mrsf_nac_amplitude` entry (Code Plan item 4 says "Subtract the SCF/ground gradient per RISK-3b-2" — make this non-optional and specify *which* gradient: the `grd2`+`sf_1e_grad` evaluated with `p2=0, td_p=dmat only, W=0, spc=0`, i.e. the stock ground ROHF gradient). Promote RISK-3b-2 from "flag this" to a hard precondition. Delete the "cancels by antisymmetry" alternative — it is false.

## RISK-D (MEDIUM) — `W^{IJ}` is not obviously well-defined off-diagonal; the spec asserts state-freedom without checking `mrsfrowcal`'s internal target dependence.

**Why.** The spec claims `mrsfrowcal` is "state-free in its formula." The energy-weighted/relaxed W (Handy–Schaefer Lagrangian) is conventionally a *single-state* object; for I≠J the interstate W is `½[W(I,J)+W(J,I)]` and must be built from interstate `xk=Z^{IJ}`, interstate `hx^{IJ}`, **and** the interstate `ppija/ppijb` (apb response of the interstate relaxed density). The spec asserts all three feed in cleanly, but `mrsfrowcal` (`tdhf_mrsf_lib.F90:2242-2383`) also takes `mo_energy_a` and may contain SOMO/`xlr` special-casing analogous to `mrsfesum` (`:1474-1544`). The spec's RISK-4-1 *assumes* all SOMO handling is "entirely localized in the hx build" and that W "carries its own SOMO handling via fa/fb/xk/hx." That is an assertion, not a verification — I did not see the body of `mrsfrowcal` checked in the spec, and if it has its own `lr1/lr2` folding driven by a *single* amplitude (like `mrsfesum` does at `:1474`), then feeding it a single `xk=Z^{IJ}` may not be the correct symmetric interstate W.

**Resolution.** Read `mrsfrowcal:2242-2383` and confirm every X-dependence enters only through `xk/hxa/hxb/ppij*` (which you control) and not through a hidden `bvec(:,target)` or SOMO fold keyed to one state. If it has internal SOMO logic, the interstate W needs the explicit `½[W(spcI,X_J)+W(spcJ,X_I)]` symmetrization, same as `hx`. The I=J self-test will NOT catch this (it passes by construction at I=J) — you need an independent I≠J check of W against Route A's implied W (extractable from the `½[G₊−G₋]` W-channel).

## RISK-E (MEDIUM) — Double-counting risk between the 2e bilinear `p2` and the z-vector relaxed density is real and the spec's "let td_p carry it" plan is the right mitigation but must be asserted.

**Why.** The 2e bilinear `get_density` consumes `p2` (via `OQP_td_p`) in `df1`/`dq1` (`:787-799`) AND the channel cross-products `spcI×spcJ` (`:801-909`). The z-vector writes `td_p = 0.5·(T^{IJ}_unrel + relax(Z^{IJ}))` via `sfropcal` (`:1496`/`tdhf_sf_lib.F90:704-763`). The danger: `T^{IJ}_unrel` (the unrelaxed difference density from `mrsf_interstate_tden`) appears **both** as the occ-occ/vir-vir blocks of `td_p` (sfropcal `:724-733`) **and** is the seed that built the `spcI/spcJ` channels — but these are *different* objects (T is the difference density; spc are the 7 transition channels), so there is no direct double count. The genuine risk is that `sfrorhs` *also* adds `2·Fa·Tij + 2·Fb·Tab` into the RHS (`:375-390`), so `Tij/Tab` influence Z^{IJ}, which then re-enters `td_p` via `sfropcal`'s Z blocks. This is correct relaxation coupling, not double-counting — **provided** you do not ALSO add a separate `p2^{IJ}` to `grd2` (the spec correctly says "you do not hand a separate p2 to grd2"). 

**Resolution.** The spec's plan (let the z-vector pipeline write `td_p`, grd2 reads it) is correct and avoids double-counting. But assert it: the I=J self-test (gate #3) must compare the *assembled* `td_p` against production `td_p` to 1e-12 (the spec lists this) AND verify the gradient is run with exactly ONE `td_p` source. Add an assertion that `grd2_mrsf_nac_compute_data_t` is never passed a non-zero `p2` argument *in addition to* reading `td_p` — the current self-test (`:1011`) passes `p2 = pB` built from `td_p`, which is correct, but the single-pass entry (Code Plan item 4) must use the identical wiring, not a second p2.

## RISK-F (LOW, but assert-once) — Channel-7 `ball` and the mult-pairing assertions are fine but the self-test doesn't currently cover them.

**Why.** RISK-2a-1 (mrsfcbc `ball` == `td_abxc`) is already *exonerated in the tree* — the in-code note at `tdhf_mrsf_z_vector.F90:1678-1682` and memory both confirm `ball ≡ td_abxc` to 1e-15, and commit `604501a` is literally "channel-7 exonerated." So the spec's "must assert once" is satisfiable but the current `mrsf_nac_amp_selftest` (`:928-1035`) does **not** build `spc` via `mrsfcbc` independently — it reads `td_mrsf_density` from the tagarray (`:989`), which already had `ball` overwritten by `td_abxc` at `:1682`. So the self-test as written cannot detect a `mrsfcbc`-vs-`sfdmat` channel-7 divergence for arbitrary I. RISK-4-2 (mult sign) is structurally safe since `sgnk` (`:766`) and the `:1697` triplet flip use the single shared `mrst`.

**Resolution.** Add the one-time `max|spc(7)_cbc − td_abxc| < 1e-12` assertion *inside* a fresh `mrsfcbc` call in the extended self-test, not from the pre-overwritten tagarray. Low effort, closes the open question for arbitrary I.

---

## Answers to the five posed questions

1. **Does interstate RHS reduce to `sfrorhs` at I=J?** Yes for `tij/tab` (bit-exact: `mrsf_interstate_tden` `:1744-1745` sets `xj=xi`, and the driver already calls it with `target,target` at `:1735` — it IS the production call). **But the spec's bilinear `hx^{IJ}` is NOT proven to reduce, and is the gap (RISK-B):** production `hx` is bilinear across two factors (channels in `fmrst2` + `wrk3`/`xv` in both the dgemm `:1713` and `mrsfsp` `:1727`); the symmetrization must cover both, including the `:1697`/`:1700-1705` sign/rescale. At I=J it collapses correctly by construction, but a partial symmetrization passes I=J while being wrong at I≠J. No hidden factor in `sfrorhs` itself (it's state-index-free), but the inputs are not all proven.

2. **Is `p2^{IJ} = T_unrel + relax(Z)` symmetric under I↔J and does it reduce to `td_p`?** The *unrelaxed* part is symmetric (`mrsf_interstate_tden`'s explicit `½(XiXj^T+XjXi^T)` at `:1750-1768`). The *relaxed* part inherits I↔J symmetry only if `Z^{IJ}` does, i.e. only if `rhs^{IJ}` is symmetric — which requires `hx^{IJ}` symmetric (RISK-B). Reduces to `td_p` at I=J exactly (sfropcal is state-free, the calls are the production calls). So symmetry is *contingent on RISK-B being done correctly*, not automatic.

3. **1e/W reduction and is `W^{IJ}` well-defined off-diagonal?** The `½[T^{IJ}+(I↔J)]` density reduces to production `td_p` at I=J exactly (per the bilinear → quadratic collapse). `W^{IJ}` reduces to `WAO` at I=J. **But "well-defined for I≠J" is asserted, not verified (RISK-D):** the interchange/Handy–Schaefer W is conventionally single-state; the interstate symmetric W requires `mrsfrowcal` to have *no* internal single-target SOMO fold. That needs checking in `mrsfrowcal:2242-2383` before claiming Step 4 is "no new routine."

4. **Double-counting between the done 2e bilinear and the new z-vector/relaxed pieces?** No double-count **if** the spec's plan is followed literally: `td_p` is the *only* relaxed-density channel into `grd2` (not a separate `p2`), and `Tij/Tab` entering both `sfrorhs` (RHS, `:375-390`) and `sfropcal` (td_p blocks) is correct relaxation coupling, not duplication. The risk is an implementation slip where the single-pass entry passes a second `p2` alongside `td_p` (RISK-E). Assert one-`td_p`-source.

5. **Is the zero-densities differencing doctrine needed?** It is needed for the **CPHF/overlap override** path (which zeros `td_abxc/td_mrsf_den/tij/tab/hxa/hxb` at `:1817-1822`) — correctly flagged mutually exclusive (RISK-CP-1, real: confirmed both modes touch `rhs` and the densities). For the **amplitude single-pass Route B**, the differencing doctrine does NOT save you: the SCF reference (`d2·d2` in `df1`/`dq1` `:787-799`, and `eijden`+`dmat_a` in `sf_1e_grad`) is X-independent but **additively present** and does **not** cancel by I↔J antisymmetry in a single pass (RISK-C). You must explicitly subtract the ground SCF gradient. So: doctrine not directly applicable, but its *purpose* (remove the SCF reference) still must be achieved by an explicit subtraction.

---

## Overall verdict

**Not ready as-is — RISK-A must be resolved first.** Before any coding, run the single decisive test: evaluate the proposed Route B `G_IJ/ΔΩ` for the ground×ijlr (S0↔open-shell-singlet) pair and compare to `d_num`. Because Route B is the hand-symmetrized bilinear of the same quadratic chain that the project already proved deficient for exactly this pair (parallelogram law, 1e-4), it will almost certainly fail there — which means Steps 2–4 are differentiating the *wrong operator* (the gradient chain) and must be re-scoped onto the **Davidson matvec derivative** (`mrsfesum`/`mrsfmntoia` pairing) that Phase 11 was defined to use. If that test unexpectedly passes, then proceed but treat RISK-B (hx bilinear, two-factor + rescale) and RISK-C (mandatory SCF subtraction) as hard gates, since those are the two places a partial implementation passes the I=J self-test while being silently wrong off-diagonal.

Key file:line anchors I verified: interstate tden `tdhf_mrsf_lib.F90:1716-1770`; hx two-factor bilinear `tdhf_mrsf_z_vector.F90:1697-1727` + `mrsfsp` SOMO fold `tdhf_mrsf_lib.F90:1800-1890`; `sfrorhs` adds `2·F·T` `tdhf_sf_lib.F90:375-419`; bilinear `get_density` keeps `d2·d2` SCF reference `tdhf_mrsf_gradient.F90:787-799`; `sf_1e_grad` state-independent `eijden`+`dmat` `tdhf_sf_gradient.F90:239-254`; self-test reads pre-overwritten `td_mrsf_density` `tdhf_mrsf_gradient.F90:989`; CPHF override mutual-exclusion `tdhf_mrsf_z_vector.F90:1769-1823`; Route A polarization `single_point.py:1747-1808`. Contradicting memory: `~/.claude/projects/-bighome-alireza/memory/nac-project.md` "RETIRED ROUTES" (hand-symmetrized bilinear ≡ polarization, reproduces the deficiency) and the Phase-11 mandate (matvec derivative, not gradient chain).