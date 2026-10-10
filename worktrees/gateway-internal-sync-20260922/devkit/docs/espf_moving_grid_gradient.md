# ESPF gradient: the fitting grid moves with the atoms

The ESPF charge operators are fitted on an atom-centred Lebedev grid: every
grid point `g_k` is placed at `R_J + r_J n_k` on its parent atom `J` and
therefore translates rigidly with that atom. The embedding energy at a fixed
MM potential,

    E = sum_i a_i sum_k Z_{i,k}(R, g) u_k(g),     Z = (T diag(s) T^T)^-1 T diag(s)

depends on the grid positions through the kernel `T_{i,k} = 1/|R_i - g_k|`,
the smooth switching weights `s_k(g)`, and the electronic ESP at the points,
`u_k = Tr[P V(g_k)]`.

Until 2026-09 the analytic gradient (`grad_esp_qmmm`, `grad_esp_qmmm_excited`,
`espf_grad_weight` in `source/modules/qmmm.F90`) differentiated all of these
with respect to the atomic positions **while holding the grid fixed in
space**. The energy that is integrated in MD, however, is computed on a grid
that follows the atoms, so the force was not the derivative of the energy:

| check (alanine dipeptide, QM = C-terminal amide + link H, HF/6-31G) | before | after |
| --- | --- | --- |
| fixed-potential QM gradient vs finite difference (Ha/bohr) | 6e-5 | 2.4e-7 |
| sum of the fixed-potential QM gradient (must vanish) | 1e-4 | 1e-12 |
| total QM/MM force vs finite difference, worst atom (kJ/mol/nm) | 1.8 | 0.24 (0.08 with the 12-decimal geometry string) |
| net force on the whole system (kJ/mol/nm) | 7.2 | 1.5e-5 |
| whole-molecule control (H2O in 3 point charges, Ha/bohr) | 1.6e-4 | 3.5e-9 |

The same defect was present for whole-molecule QM regions; it was merely
larger at a covalent boundary because the link hydrogen sits ~0.4 A from the
MM host atom and its grid overlaps a large potential.

## The added terms

Every term depends on `g_k` only through differences `g_k - R_i` (kernel and
switching) or through the point at which the ESP is evaluated (integral
term). By translational invariance the derivative with respect to `g_k` is
minus the sum over atoms of the corresponding atom derivatives, and since
`g_k` rides on its parent atom `J = parent(k)` that counter-term is charged
to `J`:

* integral term: `grad_elpot` returns the basis-centre derivatives of
  `w_k Tr[P V(g_k)]`; their negative sum is added to `parent(k)`;
* kernel term `dT/dR` in `espf_grad_weight`: the `(g_k - R_i)/r^3`
  contribution to atom `i` is subtracted from `parent(k)`;
* switching term `ds_k/dR_J`: likewise, the contribution to atom `J` is
  subtracted from `parent(k)`.

`form_espf_grid` now returns the parent index of every kept grid point
(optional argument, maintained through the smooth-weight compaction), and the
gradient routines pass it down. The overlap (charge-conservation) term and
the classical `q dphi/dR` coupling force are unchanged.

This is the "grid derivative" of Huix-Rotllant & Ferre, JCTC 2021, 17, 538,
and the `T^x` term of Bonfrate, Ferre & Huix-Rotllant, JCTC 2024, 20, 4338,
eq 13, which the link-atom gradients of that paper rely on.

MRSF-TDDFT (BHHLYP) excited-state forces through the NAMD driver agree with finite
differences to ~1e-4 Ha/bohr on the QM atoms with an exactly vanishing net force; the same residual
is present for the gas-phase MRSF gradient (2.6e-5 Ha/bohr at fixed zero potential), so it is the
MRSF gradient's own accuracy, not an ESPF or link-atom term.

Regression test: `tests/test_espf_moving_grid_gradient.py` (needs the compiled
runtime) checks translational invariance and finite-difference agreement of the
fixed-potential gradient.
