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

## RESOLVED: full analytic HF Hessian matches PySCF (0.56%, = OQP/PySCF baseline)

The response term is now correct. Learning from the PySCF bridge
(`external.analytic_hessian_from_pyscf` delegates to `pyscf/hessian/rhf.py`),
the response is `hess_elec`'s three-term form:

    H^resp_xy = 4 Tr[F^x dm1^y] - 4 Tr[S^x (eps.dm1^y)] - 2 Tr[s1oo^x mo_e1^y]
    dm1^y_pq   = sum_k dC^y_pk C_qk
    mo_e1^y_kl = (h^y + G[P]^y + G[dP^y])^MO_kl - 1/2(eps_k+eps_l) s1oo^y_kl   (full occ-occ)

Two bugs were fixed:
1. **Off-diagonal occ-occ term.** The third term must contract the FULL occ×occ
   `s1oo^x_kl mo_e1^y_kl`; the earlier diagonal-only `dε_k S^x_kk` dropped the
   off-diagonal Lagrangian-derivative contribution (the missing structural term).
2. **d-function normalization.** `der_overlap/kinetic/nucattr_matrix` return the
   UNNORMALIZED basis; they must be scaled by `bfnrm_mu bfnrm_nu` to match the
   normalized MO coefficients/density before use in the CPHF RHS and the
   response. Invisible for s/p-only bases (STO-3G), ~10% error with d in 6-31G*.

Result on H2O/6-31G* (RHF): analytic vs OQP numerical = 6e-5; analytic vs PySCF
analytic = 0.56% (identical to OQP-numerical vs PySCF, i.e. the residual is the
OQP/PySCF integral-convention baseline, not a Hessian error). Frequencies match
to <0.1 cm^-1.

Implemented in `hf_hessian` using `fock_deriv_contract` (=1/2 Tr[M G[P]^x]) and
`fock_jk` (G[dP^y]) for all 2e traces.

Validated (analytic vs PySCF, the OQP/PySCF integral baseline):
  - H2O/6-31G* RHF : 0.56%   (vs OQP numerical 0.006%)
  - NH3/6-31G* RHF : 0.31%

## DFT (RKS) Hessian -- working (bhhlyp H2O/6-31G* 9.3% vs OpenQP numerical)

The RKS Hessian reuses the RHF response (the CPHF solver is already XC-aware via
`tddft_fxc`, so dP^y is the correct DFT relaxed density) and adds the XC
contribution in two pieces, mirroring PySCF `hessian/rks` but realised entirely
through the OpenQP **moving-grid** XC machinery so it stays consistent with the
OpenQP numerical Hessian (which differentiates the same moving-grid XC gradient):

  - **dHse (skeleton + term-1)** = central FD of the analytic XC gradient
    `derexc_blk` along the relaxed path `R+lambda, P+lambda*dP^y`. This is the
    genuine total derivative `d/dR[g_XC(R,P(R))]` of OpenQP's XC gradient, so the
    grid-weight derivatives are treated exactly as in the SCF/numerical gradient.
    (It equals PySCF's `vxc_deriv2`/`veff_diag` skeleton + the XC part of
    `4 Tr[h1ao dm1]` combined.)
  - **dHt3 (energy-weighted term-3)** = `+2 Tr[s1oo^x (vxc^y+fxc[dP^y])_oo]`,
    from the FD of the XC Fock matrix `dftexcor` along the same relaxed orbital
    path `C_occ +/- h dC^y`. This is the XC part of the `mo_e1` (Lagrangian)
    contribution that the Pulay/energy-weighted term carries.
  - a warm-up `dft_initialize`/`dftclean` before the loop flushes stale grid
    state left by the CPHF solver (else the first column is garbage, ~5000%).

What PySCF does (hessian/rks.py + hessian/rhf.py), and how it maps here:
  - `partial_hess_elec`: `_get_vxc_diag` + `_get_vxc_deriv2` (XC skeleton) on a
    grid that moves with the nuclei -> matches our dHse skeleton part (norms
    agree: PySCF 0.579 vs ours 0.581).
  - `make_h1`/`_get_vxc_deriv1`: XC part of `h1ao` (term-1) built on a **static**
    grid (no weight derivatives -- those live in the skeleton). PySCF norm 0.284.
    Our moving-grid `vxc^x` is ~2x larger (0.593) because it carries the weight
    derivatives; that is the CORRECT split for OpenQP because our skeleton is the
    *moving*-grid `derexc` FD -- the two pieces sum to the same total, only the
    decomposition differs. (A static-grid `dftexcor` derivative is not available:
    OpenQP's grid carries equilibrium AO-screening tables that are invalid for a
    displaced basis, so a fixed-grid/moving-basis evaluation returns garbage.)
  - `hess_elec`: `-2 Tr[s1oo^x mo_e1^y]` with `mo_e1` carrying the XC Fock
    derivative -> our dHt3. The OpenQP moving-grid sign convention makes this
    `+2` (validated against the OpenQP numerical Hessian, the ground truth).

Validation (bhhlyp/6-31G*, H2O):
  - analytic vs OpenQP numerical Hessian: 9.3% max / 11.3% Frobenius.
  - analytic vs PySCF RKS: 9.5% (6d-vs-5d basis + grid-convention baseline).
  - frequencies: stretches 4073/4147 vs num 4032/4161 (<2%); bend 1561 vs 1724
    (~9% soft -- the residual is concentrated in the bend, i.e. the angular block
    of the moving-grid FD).

Remaining: the ~9% residual is the moving-grid finite-difference floor (the
analytic XC pieces are FD of grid quantities at a single step; the numerical
Hessian re-grids per geometry). Closing it needs genuine analytic XC second-
derivative grid integrals (PySCF's `_get_vxc_deriv2`/`_get_vxc_deriv1` kernels:
contractions of v2rho2/v2rhosigma/... and v3rho3/... against AO derivatives),
which OpenQP does not yet expose. Also pending: UHF/ROHF response.

## CPHF response term — detailed findings (historical, now resolved)

The skeleton (1e + 2e + nn) is now **proven exact**: `hess_skel_selftest`
finite-differences the full frozen-density gradient (physical energy-weighted
density W, not a placeholder) and matches the analytic skeleton to 1.9e-7. So
the entire residual error of the total HF Hessian is localised to the **CPHF
orbital-relaxation response** term in `hf_hessian` (`matmul(bvec,uvec)`).

Established facts:
- `bvec`/`uvec` are correct: `cphf_dpdx_selftest` validates that U^y
  reproduces the relaxed density derivative dP^y (O(h²), ~1e-6).
- `cphf_solve` is `(A+B)U=B` with `A_diag=(e_a−e_i)`; calibrated by the
  polarizability path `α=−4 Σ μ U` (factor 4 for this solver).
- `fock_deriv_contract(P,M)` returns **½·Tr[M G[P]^x]** (the validated `bvec`
  uses `2*gx`); any full trace built from it needs a factor 2.
- `eijden` (OQP energy-weighted density): `W_OQP = −2 Σ_i ε_i C_iμ C_iν`, and
  the gradient overlap term is `+Tr[W_OQP S^x]`.

The correct closed-shell form is the density-derivative interchange
    H^resp_xy = Tr[dP^y F^x] + Tr[dW^y S^x],   F^x = h^x + G[P]^x
(derived from differentiating the converged gradient
 G^x = Tr[P h^x] + ½Tr[P G[P]^x] + Vnn^x + Tr[W_OQP S^x] through P and W).
This was implemented (build dP^y, dW^y from the validated dC^y; all 2e traces
via fock_deriv_contract×2 and fock_jk) and reproduces the individual pieces,
but a least-squares fit of the **known** target response (Hn − exact skeleton)
to all six computed sub-pieces still leaves a **21% structural residual** — i.e.
a 7th term, orthogonal to every U/dP-based piece, is missing. The pieces are
each ≈ a multiple of `sym(bvec^T uvec)` (e.g. `Tr[dP^y F^x] ≈ −4·R`), and the
target ≈ `+2.5·R` with corr 0.81, so neither the bare `matmul(bvec,uvec)`
(currently shipped, ~27%) nor the density-derivative form (~37%) is complete.

Next step: identify the missing structural term (candidate: a pure
overlap-derivative reorthonormalization contribution Σ_ij ε_i S^x_ik S^y_kj not
captured by dW^y, or the second-order dP^x·dP^y coupling through the orbital
Hessian). The numerical target response is fully available
(`Hn − hess_skel_selftest`), so this is a bounded reverse-engineering task.
