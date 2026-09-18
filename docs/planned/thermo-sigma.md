# Rotational symmetry number is missing from the thermochemistry

## The defect

`pyoqp/oqp/library/frequency.py:198-213` builds the rotational partition function
with **no symmetry number**. `grep` for `symmetry_number` / `sigma_rot` across
`pyoqp/` and `source/` returns nothing.

Every gas-phase entropy and Gibbs energy for a symmetric molecule is therefore
wrong by `R ln σ`:

| molecule | σ | error in G(298.15) |
|---|---|---|
| H2O | 2 | 0.411 kcal/mol |
| CH4, C6H6 | 12 | 1.47 kcal/mol |

Every reaction and activation free energy is off by the net difference.

## Two coupled bugs found alongside it

This is **not** a one-line patch.

1. **Linear molecules print `st_rot = -inf`.** `rt` comes from `np.linalg.eigh`
   in ascending order (`frequency.py:100`), so `rt[0]` is the near-zero moment
   and `qr = 0` in both the `len(atoms)==2` branch and the `else` branch.
2. **The `linear` flag is never passed.** `single_point.py:1128-1136` does not
   supply it, so `u_rot` is 1.5RT for linear species and the `+1` entropy branch
   at `:211` is dead code.

## Fix

1. Divide `qr` by σ, counted from the `det > 0` operations in
   `enumerate_full_group` for nonlinear cases, and from the `dooh` (2) / `coov`
   (1) label for linear ones. σ has to be computed at the call site —
   `thermal_analysis` has no coordinates.
2. Select the nonzero principal moments for linear molecules.
3. Thread the `linear` flag through.
4. Print σ next to the detected point group, so a tolerance-induced wrong σ is
   visible rather than silent.

## Verification

Textbook σ for H2O (2), NH3 (3), CH4 (12), C6H6 (12), staggered C2H6 (6), CO2
(2), CO (1), H2 (2); total S and G cross-checked against ORCA or Gaussian for the
same geometry and frequencies.

**Existing thermochemistry regression numbers will all shift and must be
rebaselined.**

## Risk

A geometry converged to only 1e-5 may detect as C1 and silently take σ = 1.
