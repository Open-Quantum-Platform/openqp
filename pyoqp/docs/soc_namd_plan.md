# SOC-NAMD Energy Conservation: Diagnosis and Improvement Plan

**Branch:** `namd-qmmm`  
**Date:** 2026-06-09  
**Status:** Option 2 and Option 3 planned for implementation

---

## Background: Completed fixes

| Commit | Fix |
|--------|-----|
| `41fc343b` | Include ΔE_ESPF in velocity rescaling at ISC hops |
| `58b37cb9` | Recommend `thrshe=0.1` Ha for SOC-NAMD (docs/comment) |
| `f786cc9d` | Honour `ESPF_ROHF=1` in SOC gradient loop |

These fixes resolved issues in the ESPF energy bookkeeping at hops and prevent
unphysical large-gap hops at the FC geometry when `thrshe=0.1` is set.

---

## Root cause of residual SOC energy drift

### IC-NAMD (no SOC): exact

For internal-conversion NAMD the active surface is a single adiabatic state _a_.
The force equals the exact gradient of that state:

```
F = -dE_a/dR
```

Because F = -∇V exactly, total energy is conserved to integration accuracy.
Expected drift for H₂CO + 5 TIP3P over 500 fs: **< 0.1 kJ/mol**.

### SOC-NAMD (spin-adiabatic representation): approximate

The SOC states are spin-adiabats obtained by diagonalising the SOC matrix **H**:

```
H = U† · diag(E_i^MCH) · U
```

where `E_i^MCH` are the MCH (spin-pure: S0, T1_ms, S1, T2_ms) energies and **U**
is the transformation.  The active spin-adiabat _a_ is a mixture of MCH states,
so its exact gradient is:

```
dE_a/dR = Σ_i |U_ia|² dE_i^MCH/dR   +   Σ_{i≠j} 2Re[U*_ia U_ja](E_j-E_i) dU_ia/dR
```

**Current implementation (SHARC diagonal approximation):**

```
F ≈ -Σ_i |U_ia|² dE_i^MCH/dR        # only the first term
```

The missing second term is:

```
ΔF_corr = -Σ_{i≠j} 2Re[U*_ia U_ja](E_j - E_i) dU_ia/dR
```

This term is nonzero whenever:
- The active state is a strong mix of MCH states (|U_ia| ~ |U_ja|), AND
- The SOC matrix changes along the trajectory (dU/dt ≠ 0).

This is the **irreducible source of energy non-conservation** in the SHARC/spin-adiabatic
representation.  It is not a bug; it is a structural approximation of the method.

Expected drift with diagonal approximation for H₂CO + 5 TIP3P over 500 fs:
**~1–20 kJ/mol** depending on how close the trajectory comes to a spin-adiabat crossing.

---

## Planned Solutions

### Option 2: Finite-difference dU/dR correction (cheap)

**Idea:** Approximate `dU_ia/dR` via a time-finite-difference using the stored
`prev_u` (U at step _t-dt_):

```python
dU_dt = (u - self.prev_u) / dt_au          # shape (nstate, nstate), a.u.
```

Then add the correction force for each atom _I_:

```
ΔF_I,corr = -Σ_{i≠j} 2Re[U*_ia U_ja] (E_j - E_i) (dU_ia/dR_I)
```

Since `dU/dt ≈ (dU/dR)·(dR/dt)` and `dR/dt = v`, one can recover a
direction-projected correction; or, at each time step, use the scalar:

```python
correction_factor = 2 * Re[ (u.conj() @ (u - prev_u) / dt) * delta_E_matrix ]
```

**Location:** `_soc_gradient_qmmm` or at the end of `NAMD_SOC_QMMM.run()`, after
computing the weighted gradient `g_qm`.

**Pros:**
- No extra QM calls
- Uses already-stored `self.prev_u`
- Straightforward ~20-line addition

**Cons:**
- Time-derivative dU/dt is a proxy for the spatial gradient dU/dR; valid when
  the nuclear velocity direction is roughly constant (short dt)
- Does not fix the gradient during the very first step (prev_u not yet available)

**Implementation outline:**

```python
# In _soc_gradient_qmmm (or run()), after computing g_weighted:
if hasattr(self, 'prev_u') and self.prev_u is not None and dt_au > 0:
    dU_dt = (u - self.prev_u) / dt_au          # (nstate, nstate) complex
    nst = eval_ha.shape[0]
    a = self.active - 1                         # 0-indexed active state
    for i in range(nst):
        for j in range(nst):
            if i == j:
                continue
            dE = eval_ha[j] - eval_ha[i]       # Hartree
            coeff = 2.0 * np.real(u[i, a].conj() * u[j, a] * dU_dt[i, a])
            # coeff has units of 1/time; multiply by |v| to get force direction
            # approximation: project along nuclear velocity direction
            v_norm = self.vel / (np.linalg.norm(self.vel) + 1e-30)
            g_weighted -= coeff * dE * v_norm   # shape (natom, 3), Hartree/Bohr
```

---

### Option 3: MCH-basis FSSH (rigorous, matches GAMESS)

**Idea:** Run the FSSH entirely in the MCH (spin-pure) basis.  The active state
is one of the MCH states (S0, T1_ms, S1, T2_ms), and the force is the exact
single-state MCH gradient — no mixing, no missing dU/dR term.

**How GAMESS does it:**
- States: S0, T1 (×3), S1, T2 (×3)
- Active state: one MCH state label (mult, state)
- Force: exact gradient of the active MCH state
- Electronic propagation: TDSE in MCH basis using H_SOC as the off-diagonal coupling
- Hop decision: compare densities `d_{jj}` in MCH basis; hop from MCH state _a_ → _b_
- Velocity rescaling: along NAC direction (IC) or isotropic (ISC); gap = E_b^MCH − E_a^MCH

**Implementation changes in `namd.py`:**

1. **`self.active`** changes meaning from spin-adiabat index (1–8) to an MCH label
   `(mult, state)` — e.g., `(1, 1)` = S0, `(3, 1)` = T1, `(1, 2)` = S1.

2. **Gradient:** Replace the weighted-MCH loop in `_soc_gradient_qmmm` with a
   single MCH gradient call for `(self.active_mult, self.active_state)`.

3. **Electronic propagation:**  Propagate the coefficient vector in the MCH basis
   using the full SOC Hamiltonian `H_SOC` as the coupling:
   ```
   iħ dc_j/dt = Σ_k H_SOC[j,k] c_k
   ```
   No need for local diabatization or `prev_u`.

4. **Hop decision:** Compute transition probability from density matrix:
   ```
   g_{aj} = -2 Re[c_a* c_j] (H_SOC[a,j] / ħ) * dt
   ```
   Accept hop `a → j` if cumulative random number threshold is crossed.

5. **Velocity rescaling:** Use MCH energy gap `ΔE = E_j^MCH − E_a^MCH`.

**Pros:**
- Exact force (no diagonal approximation)
- Direct comparison to GAMESS FSSH output
- Conceptually simpler: each state has one gradient

**Cons:**
- Requires computing `_mch_gradient(mult, state)` — currently `_soc_gradient_qmmm`
  computes all MCH gradients; need to expose single-state version
- Electronic propagation logic changes substantially (no LD propagator)
- `self.active` semantics change — careful not to break IC branch

**Phasing:** Implement as a new class `NAMD_SOC_MCH_QMMM` (alongside the existing
`NAMD_SOC_QMMM`) controlled by an input keyword, e.g., `soc_basis = mch` (default
`adiabatic` for current behaviour).

---

## Testing plan

After implementing Option 2:
1. Run `form_soc_verify.inp` (H₂CO + 5 TIP3P, 30 steps at dt=0.5 fs)
2. Confirm total-energy drift per step < 0.5 kJ/mol (vs ~5 kJ/mol without correction)
3. Compare to Option 3 on the same trajectory for cross-check

After implementing Option 3:
1. Run same test; force should be identical to GAMESS single-state force at each step
2. Cross-check hop sequence and rescaling with a GAMESS reference trajectory

---

## Files to modify

| File | Change |
|------|--------|
| `pyoqp/oqp/library/namd.py` | Option 2: add dU/dt correction in `_soc_gradient_qmmm` |
| `pyoqp/oqp/library/namd.py` | Option 3: add `NAMD_SOC_MCH_QMMM` class |
| `pyoqp/oqp/molecule/oqpdata.py` | Add `soc_basis` keyword (`adiabatic`/`mch`) |
| `pyoqp/docs/soc_namd_plan.md` | This file (update as work progresses) |
| `overleaf_namd/main.tex` | Add planned improvements section |

---

## Summary of energy drift expectations (H₂CO + 5 TIP3P, 500 fs)

| Method | Expected drift |
|--------|---------------|
| IC-NAMD (no SOC) | < 0.1 kJ/mol |
| SOC-NAMD, diagonal approx (current) | 1–20 kJ/mol |
| SOC-NAMD + Option 2 (dU/dt correction) | ~0.1–1 kJ/mol |
| SOC-NAMD + Option 3 (MCH basis, exact) | < 0.1 kJ/mol |
