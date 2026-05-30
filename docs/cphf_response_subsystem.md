# CPHF/CPKS Response Subsystem — Validated, Frozen

This records the validated native orbital-response subsystem as a frozen
checkpoint. It is an input to the analytic Hessian but is **independent of** the
second-derivative integral work (see `second_derivative_integral_design.md`).
Do not modify these modules as part of derivative-2 development.

Checkpoint commit: `e14a75e` (tag intended: `cphf-response-validated`).

## Scope

Closed-shell (RHF/RKS) coupled-perturbed orbital response for nuclear
perturbations, the chain

```
F^x  ->  B^x  ->  U^x  ->  dC/dx  ->  dP/dx
```

with no libint dependency (uses the native Rys 2e engine and the DFT XC kernel).

## Modules (frozen)

| Module | Public entry | Role |
|---|---|---|
| `source/integrals/grd1.F90` | `der_overlap_matrix`, `der_kinetic_matrix`, `der_nucattr_matrix` | AO first-derivative matrices dS/dx, dT/dx, dV/dx (nbf,nbf,3,natom) |
| `source/modules/fock_deriv.F90` | `fock_deriv_contract` | 2e derivative-Fock contraction g_x = sum_uv M_uv F^x_uv[P] (F^x builder) |
| `source/modules/cphf.F90` | `cphf_solve`, `cphf_polarizability_selftest` | CPHF/CPKS PCG solver A U = B (A = orbital Hessian via the native response Fock) |
| `source/modules/fock_deriv_selftest.F90` | `fockx_selftest` | F^x M=P trace-identity validation |
| `source/modules/cphf_nuclear_selftest.F90` | `cphf_f0x_selftest` | F^x M!=P gateway validation |
| `source/modules/cphf_dpdx_selftest.F90` | `cphf_dpdx_selftest` | end-to-end dP/dx validation |

## Validation evidence

All validated against exact/non-iterative references with O(h^2) scaling
(truncation-limited, not reference-limited):

- **F^x, M=P trace identity** (`fockx_selftest`): contraction(P,P) = dE_2e/dx vs
  FD of tr(P G[P]) at frozen P. max|an-fd| = 1.4e-8; O(h^2)
  (5.1e-8, 1.4e-8, 4.0e-9 at h=2e-4,1e-4,5e-5). H2O/6-31G*.
- **F^x, M!=P gateway** (`cphf_f0x_selftest`): d/dx tr(M G[P]) = 2
  fock_deriv_contract(M,P) for a fixed AO probe M != P, vs FD of the
  basis-consistent scalar tr(M G[P]). relative error 3.0e-9; O(h^2). H2O/6-31G*.
  (This caught and fixed a Coulomb-symmetrization bug: the energy-gradient
  density product 4 c D D generalizes to the symmetrized 2 c (M P + P M) for
  M != P.)
- **dS/dx, dT/dx, dV/dx** (`hess1_selftest`): contracting each with the
  normalized density reproduces the production gradients
  (grad_ee_overlap, grad_ee_kinetic, grad_en_pulay + grad_en_hellman_feynman) to
  ~1e-14 - ~1e-16. H2O/6-31G*.
- **cphf_solve** (`cphf_polarizability_selftest`): static dipole polarizability,
  symmetric and positive, magnitudes matching a PySCF field finite-difference
  reference (signed agreement at the FD reference's noise floor).
- **End-to-end dP/dx** (`cphf_dpdx_selftest`): analytic d(P_alpha)/dx vs central
  FD of DM_A/2 (alpha density), H2O/STO-3G:

  | h (bohr) | abs_err | rel_err | ratio |
  |---|---|---|---|
  | 8.0e-3 | 9.88e-6 | 2.61e-5 | - |
  | 4.0e-3 | 2.46e-6 | 6.49e-6 | 4.02 |
  | 2.0e-3 | 6.14e-7 | 1.62e-6 | 4.00 |
  | 1.0e-3 | 2.13e-7 | 5.62e-7 | (noise floor) |

  Textbook O(h^2) (error /4 per halving) into the SCF/FD noise floor (~2e-7).

## Key conventions (verified from production code; do not re-derive by guessing)

- **OQP::DM_A is the TOTAL closed-shell density**: DM_A = 2 P_alpha,
  P_alpha = sum_occ C C^T (idempotent, P_alpha S P_alpha = P_alpha).
- **fock_jk expects the TOTAL density**: production Fock F = Hcore + fock_jk(D),
  D = DM_A (`scf_addons` calc_jk_xc).
- **Fortran mo_a is indexed (AO, MO)**; C^T S C = I (proven in hess1_selftest).
- **cphf_solve operates in bare/alpha amplitude scale** (trial density from the
  bare occ-vir amplitude; validated by the polarizability test). The RHS must be
  in the same scale.
- **Occupied-occupied overlap response** (from C^T S C = I, U + U^T + S^x = 0):
  canonical gauge U_ij = -1/2 S^x_ij; the relaxed density occ-occ block is
  dP_alpha^{oo} = -sum_ij S^x_ij C_i C_j^T (the two mirror terms double -1/2 to
  -1); the response density fed to fock_jk (total scale) is
  d0 = -2 sum_ij S^x_ij C_i C_j^T.

## Boundary rule

Second-derivative integral development (the Hessian's d2h, d2S, d2(uv|ls) terms)
is a separate workstream. It must not modify the modules above. It consumes the
response subsystem only at final Hessian assembly, via `cphf_solve` and the
existing public derivative-matrix builders.
