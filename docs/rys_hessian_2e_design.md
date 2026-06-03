# Analytic two-electron Rys second-derivative (Hessian) kernel — design

Status: **2e ∂² kernel IMPLEMENTED & VALIDATED** (H2O/6-31G*, rel err 1.4e-7
vs finite difference of the analytic 2e gradient; symmetric to 5e-16). Driver
`grd2_hess_driver` + per-quartet `grd2_rys_hess_compute` are live; self-test
`grd2_hess_selftest` / `tests/test_hess_2e_selftest.py`. Remaining: wire into
`hf_hessian`, then finish 1e `hess_en` and the XC ∂² skeleton.

## Goal

Add the two-electron ERI second-derivative skeleton to the analytic HF/DFT
Hessian:

    H2e[Aa,Bb] = ½ Σ_μνλσ P_μν P_λσ · ∂²/∂A_a∂B_b [ (μν|λσ) − ¼ c_x (μλ|νσ) ]

This is the dominant missing term (isolation experiment: response + nuclear +
1e skeleton still leaves ~560% rel. error on pure-HF H₂O/6-31G*, because this
term is absent).

## What already exists (reuse)

- `grd2_driver` / `grd2_driver_gen` (grd2.F90): the 2e **gradient** driver.
  Loops shell quartets (i,j,k,l), builds the 2-body density block via
  `gcomp%get_density(basis, gdat%id, dab, dabmax)` (the exact J−K combination,
  already scaled by `hfscale`), calls `grd2_rys_compute`, scatters per-quartet
  `gdat%fd(3,4)` into `de(:,gdat%at)` by atom index.
- `grd2_rys_compute` (grd2_rys.F90): per-quartet engine. Builds Rys roots,
  the 2D integrals `gnm` (`compute_xyz_p0q0`), transfers to 4-center 1D
  integrals `gijkl` (`compute_xyz_ijkl`), forms first derivatives
  `fi,fj,fk,fl` (`compute_der_xyz_ijkl`), contracts with density into
  `fd(3,4)` (`compute_der_ijkl`), then `apply_translation_invariance`.
- `mod_1e_primitives::der2_coul_xyz` — validated 1D Rys **second**-derivative
  recursion (`test_der2_recursion`: der2 == der(der) to 1e-10).
- `gdat_init` already pads angular ranges by `nder` (`mxcart = maxang+1+nder`,
  `nroots = (ΣL+2+nder)/2`), so `nder=2` allocates enough room.

## Key simplification (makes the kernel tractable)

An ERI factorizes into independent x, y, z 1D integrals; each 1D factor in
direction d couples all four centers. A derivative ∂/∂(center c, dir a) acts
only on the **a-direction** 1D factor, replacing it with its c-derivative.
Therefore the per-quartet second-derivative block decomposes as:

- **a1 ≠ a2** (different Cartesian directions): the two derivatives hit two
  *different* 1D factors →
      D_{a1}^{c1} D_{a2}^{c2} (gx·gy·gz)
        = f{c1}(a1) · f{c2}(a2) · g0(third)
  i.e. a **product of the existing first-derivative arrays** fi/fj/fk/fl.
- **a1 == a2 == d** (same direction): both derivatives hit the *same* 1D
  factor → need the 1D **second** derivative w.r.t. the center pair (c1,c2):
      f{c1 c2}(d)  =  der_{c2}( f{c1} )(d)
  built by re-applying the per-center first-derivative operator (the same
  `der2 == der(der)` identity, generalized to mixed centers c1≠c2).

So the only genuinely new integral arrays are the 10 center-pair second
derivatives {ii,ij,ik,il,jj,jk,jl,kk,kl,ll} per direction, each obtained by
applying the existing single-center derivative operator a second time:

    der_i:  A(...,i)   -> A(...,i+1)·(2α_i) − A(...,i−1)·(i−1)
    der_j:  A(...,j,i) -> A(...,j+1,i)·(2α_j) − A(...,j−1,i)·(j−1)   (etc.)

    fii = der_i(fi),  fij = der_j(fi),  fik = der_k(fi),  fil = der_l(fi),
    fjj = der_j(fj),  fjk = der_k(fj),  fjl = der_l(fj),
    fkk = der_k(fk),  fkl = der_l(fk),  fll = der_l(fl)

(operators commute; fij via der_j(fi) == der_i(fj), a free internal check.)

## Per-quartet contraction (compute_der2_ijkl)

For each AO quartet (i,j,k,l) with direction offsets nx,ny,nz and density df:

    for centers c1≤c2 in {1..4}, dirs a1,a2 in {x,y,z}:
      if a1==a2==d:
        fd2[a1,c1,a2,c2] += df · Σ_roots f{c1c2}(d)[n_d] · g0(o1)[n_o1] · g0(o2)[n_o2]
      else:
        fd2[a1,c1,a2,c2] += df · Σ_roots f{c1}(a1)[n_a1] · f{c2}(a2)[n_a2] · g0(a3)[n_a3]

`fd2(3,4,3,4)` is symmetric under (a1,c1)↔(a2,c2).

## Driver (grd2_hess_driver, mirrors grd2_driver_gen)

- Same quartet loop + screening + `gcomp%get_density` (identical density block).
- `gdat%init(mxam, nder=2, …)`; differentiate **all four centers**
  (skip(:) = .false.) for a correctness-first version; translation invariance
  (Σ_c fd2[:,c,:,:] = 0) becomes a *check*, not a shortcut.
- Scatter: for c1,c2:
      hess(3(at(c1)−1)+a1, 3(at(c2)−1)+a2) += fd2(a1,c1,a2,c2)
- ½ prefactor + permutational multiplicity handled exactly as the gradient
  (`iandj`/`kandl`/`same` halving) but applied to the 12×12 block.

## Validation strategy (gate before wiring into hf_hessian)

Fortran self-test (mirror `hess1_selftest`): with the converged density fixed,
compare the analytic 2e Hessian from `grd2_hess_driver` against a central
finite difference of `grd2_driver` (the validated 2e gradient), displacing
geometry in-process via `init_shell_centers`. Identity:
`d/dR [ Σ M·dX/dR ] = Σ M·d²X/dR²` for fixed M. Tolerance ~1e-6.

Then add `grd2_hess_driver` into `hf_hessian` and validate the **total**
analytic Hessian against the numerical reference (pure HF first, then DFT once
the XC ∂² skeleton is added).

## Status after wiring into hf_hessian (H2O/6-31G*, pure HF)

Total analytic-vs-numerical HF Hessian error dropped **627% → 27%** once the
2e skeleton was added. Contribution norms (‖·‖_F) reveal heavy physical
cancellation:

    response            0.205
    + nn                9.30      (nuclear repulsion ∂²)
    + 1e skeleton      10.24      (overlap·W + kinetic·P + en·P)
    + 2e skeleton       1.40      (ERI ∂², this kernel)  vs numerical 1.77

All four integral second-derivative pieces are individually validated against
finite difference:
  - overlap ∂²   1e-10   (hess1_selftest)
  - kinetic ∂²   4e-9
  - nucattr ∂²   1.9e-8  (hess_en — NOT WIP after all; matches FD of grad_en)
  - 2e ERI ∂²    1.4e-7  (grd2_hess_selftest)

So the residual 27% is **not** in the integral kernels — it is an assembly
issue surfacing after the ~10→1.4 cancellation. Prime suspects, in order:
  1. CPHF **response** scale/sign — norm 0.205 looks small for H2O; check the
     closed-shell occupancy factor (2/4) and sign of matmul(bvecᵀ,uvec).
  2. Energy-weighted (Lagrangian) density **W** from `eijden`, and the sign
     convention of `hess_ee_overlap(W)` (Pulay term −Σ W d²S).
A few-percent error in either large term explains the 0.37 residual.

## Remaining
- Resolve the assembly factor above (focused: instrument each term vs a
  reference decomposition, e.g. PySCF skeleton/response split).
- XC ∂² skeleton (DFT/bhhlyp only).
