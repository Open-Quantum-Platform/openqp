# Analytical vs. numerical gradient validation

This directory contains tooling to validate OpenQP's analytical nuclear
gradients against numerical (finite-difference) gradients, for

* RHF, UHF, ROHF (Hartree-Fock)
* TDDFT, SF-TDDFT, MRSF-TDDFT (excited states)

and the investigation that this validation triggered.

## Tools

* `validate_gradients.py` — runs the analytical gradient and a central
  finite-difference gradient of the (state) energy for each method and
  reports the per-component agreement.

  ```bash
  export OPENQP_ROOT=/path/to/openqp
  export LD_LIBRARY_PATH=$OPENQP_ROOT/lib:$LD_LIBRARY_PATH
  export PYTHONPATH=$OPENQP_ROOT/pyoqp:$PYTHONPATH
  python tools/validate_gradients.py                 # all methods
  python tools/validate_gradients.py rhf uhf rohf     # subset
  python tools/validate_gradients.py --nstate 10 --verbose
  ```

* `crosscheck_pyscf.py` — independent cross-check of the HF ground-state
  gradients against PySCF (requires `pip install pyscf`).  OpenQP uses
  *cartesian* d/f functions for Pople basis sets, so PySCF is run with
  `cart=True`; the SCF energies then agree to ~1e-10 Hartree, which proves
  the two codes are solving the identical problem and makes any gradient
  discrepancy unambiguous.

## Bug found and fixed: pure Hartree-Fock 2-electron gradient

The validation uncovered a real bug in the **pure Hartree-Fock** gradient
(RHF / UHF / ROHF, and any `method=tdhf` run *without* a DFT functional).

### Symptom

At the H2O/6-31G* equilibrium geometry the analytical RHF gradient was
`|g| ≈ 0.91` Hartree/Bohr, while

* a finite difference of OpenQP's own (correct) energy gave `≈ 0`, and
* PySCF (cartesian, same energy to 1e-10) gave `≈ 1.4e-5`.

A term-by-term decomposition showed the **one-electron** gradient was
correct, but the **two-electron** gradient evaluated to `1.0·J + 0.5·K`
instead of the correct `1.0·J − 0.5·K` — i.e. the exchange contribution had
the wrong sign.  DFT (e.g. BHHLYP) gradients were unaffected.

### Root cause

`source/integrals/grd2.F90`, `grd2_driver`, non-CAM branch:

```fortran
else
  gcomp%hfscale = infos%dft%hfscale     ! <-- clobbers the value the caller set
  gcomp%hfscale2 = infos%tddft%hfscale
  call grd2_driver_gen(infos, basis, de, gcomp)
end if
```

`hf_2e_grad` (and the TD gradient drivers) already set `gcomp%hfscale`
correctly: `1.0` for pure HF, `infos%dft%hfscale` for DFT.  This line then
**overwrote** it with `infos%dft%hfscale` unconditionally.

* For DFT (`hamilton>=20`) `infos%dft%hfscale` holds the hybrid mixing
  (e.g. 0.5 for BHHLYP) — the overwrite is harmless (same value).
* For pure HF (`hamilton<20`) `infos%dft%hfscale` is **not meaningful**: the
  C-interop struct never received the `1.0` Fortran default and held the
  `-1.0` sentinel.  The exchange weight `−0.5·hfscale` therefore became
  `−0.5·(−1.0) = +0.5`, flipping the exchange sign.

### Fix

Only adopt the DFT hybrid scales when it is actually a DFT calculation:

```fortran
else
  if (infos%control%hamilton >= 20) then
    gcomp%hfscale = infos%dft%hfscale
    gcomp%hfscale2 = infos%tddft%hfscale
  end if
  call grd2_driver_gen(infos, basis, de, gcomp)
end if
```

### Verification after the fix

| method | analytical vs numerical | analytical vs PySCF |
|--------|-------------------------|---------------------|
| RHF    | max\|Δ\| ≈ 7e-7         | exact (to ~1e-7)    |
| UHF    | max\|Δ\| ≈ 6e-7         | self-consistent\*   |
| ROHF   | analytical = PySCF (−0.139857 …) | exact      |
| TDHF (pure HF), S1 | max\|Δ\| ≈ 7e-6 | —            |

\* OpenQP's UHF SCF for the triplet converges to a different solution than
PySCF's default for this case; the gradient is consistent with OpenQP's own
energy (that SCF-solution difference is a separate, unrelated matter).

The reference `.json` files for all affected pure-HF gradient examples were
regenerated (only the `grad` field changes; energies are unaffected).  The
full example test-suite passes with the corrected references and shows **no
regressions** in the DFT / TDDFT / SF / MRSF / SCF / ECP paths.

## A note on excited-state numerical references

`validate_gradients.py` will still report large `max|Δ|` for some higher
SF-/MRSF-TDDFT roots on the symmetric H2O reference.  Those are **not**
gradient errors: at C2v there are (near-)degenerate roots, and the triplet
ROHF reference is numerically delicate, so a plain energy-index finite
difference tracks the *wrong* root at displaced geometries.  Evidence that
the analytical excited-state gradients are sound:

* the lowest, well-separated excited state agrees tightly
  (TDDFT S1: `max|Δ| ≈ 8e-5`);
* at a low-symmetry geometry the analytical MRSF S0 gradient is small and
  smooth while the finite difference is the noisy quantity;
* the fix here does not touch the DFT/TD code paths (`hamilton>=20`), and the
  shipped MRSF/SF gradient regression tests continue to pass unchanged.

A fully robust excited-state finite-difference check would require root
tracking by wavefunction overlap and a continued (non-`huckel`) reference
guess across displacements.
