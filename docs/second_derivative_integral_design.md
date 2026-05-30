# Second-Derivative Integral Strategy — Design Document

Status: **design only**. No second-derivative integral code has been written
against this document yet. This is the plan to be reviewed before implementation.

This document covers ONLY the missing second-derivative integral contributions
to the analytic HF/DFT Hessian. The orbital-response machinery
(F^x -> B^x -> U^x -> dP/dx) is a separate, **validated and frozen** subsystem;
see Section 1 and `cphf_response_subsystem.md`. Derivative-2 development must not
modify the response-subsystem modules.

---

## 1. Boundary with the validated response subsystem (do not cross)

The following modules are validated (see `cphf_response_subsystem.md`) and are
**off-limits** to second-derivative integral development except through their
existing public interfaces:

- `source/modules/fock_deriv.F90` (`fock_deriv_contract`, the F^x builder)
- `source/modules/cphf.F90` (`cphf_solve`, the CPHF/CPKS solver)
- `source/integrals/grd1.F90` first-derivative *matrix* builders
  (`der_overlap_matrix`, `der_kinetic_matrix`, `der_nucattr_matrix`)
- the validation harnesses `fock_deriv_selftest`, `cphf_nuclear_selftest`,
  `cphf_dpdx_selftest`, `cphf` polarizability self-test.

The second-derivative integrals are an **independent** input to the final
Hessian assembly. They do not feed the response chain and must not be tangled
with it. Final Hessian assembly (a later, separate step) combines:

```
H = H_nuc(2)            (done: grd1::hess_nn)
  + sum P d2h           (1e core-Hamiltonian 2nd deriv)  <-- this document
  - sum W d2S           (overlap 2nd deriv, energy-weighted) <-- this document
  + 1/2 sum PP d2(uv|ls) (2e 2nd deriv)                  <-- this document
  + (XC 2nd geometric deriv, DFT only)                   <-- later
  + response term from U^x (validated subsystem)          (done)
```

---

## 2. The mathematical problem (why fixed-root recursion failed)

OpenQP integrals over 1/r (1e nuclear attraction) and ERIs use **Rys
quadrature**. The Rys argument is center-dependent:

- 1e nuclear attraction: `X = pp%aa * sum((pp%r - c)**2)`
  (`mod_1e_primitives.F90:383`), `pp%r` = Gaussian-product center of the bra-ket
  pair, which depends on the basis centers.
- 2e ERI: `X = rho * sum((P - Q)**2)` (`grd2_rys.F90:309`), with `P`, `Q` the
  bra and ket Gaussian-product centers, depending on all four basis centers.

The Rys roots `u_t(X)` and weights `w_t(X)` are functions of `X`. A first
derivative w.r.t. a center `A` is

```
d/dA [I] = (d/dA of polynomial xyz tables, at FIXED roots)
         + (dX/dA) * (d/dX of roots/weights)
```

OpenQP's **first**-derivative code (`compute_der_xyz_ijkl`,
`der_coul_xyz`, `der_kinovl_xyz`) captures both contributions because it rebuilds
the full integral implicitly and the single angular-momentum shift with the
existing roots is sufficient at first order. The **fixed-root recursion-twice**
approach for the second derivative is mathematically invalid: it applies the
polynomial shift twice while holding the roots frozen, dropping the
`(dX/dA) d(roots)/dX` contribution to the second derivative. This was confirmed
empirically (1e nuclear attraction `der2_coul`: on-axis components wrong; bumping
the root count moved the answer *away* from the FD value) and is recorded in
`analytic_hessian_design.md`.

Consequence: native second-derivative nuclear-attraction and ERI integrals
require either (a) a Rys engine extended with root/weight derivatives, or (b) a
different integral scheme (Obara-Saika / McMurchie-Davidson) for the 2nd
derivative. Overlap and kinetic second derivatives have **no Rys roots** and are
already done and validated (`hess_ee_overlap`, `hess_ee_kinetic`).

---

## 3. Strategy options

Three routes, assessed against effort, risk, and the project's existing
infrastructure.

### Option A — libint deriv-2 integrals (pragmatic, chosen baseline)

OpenQP already links a custom libint (`oqp-libint`), and the wrapper
(`int_libint.F90`) already exposes a `deriv_order` interface with a
`libint2_build_eri2` (second-derivative ERI) code path gated by
`#if INCLUDE_ERI >= 2`. The libint engine computes the root-derivative
contributions correctly internally.

- **1e (Task A):** libint 1-body deriv-2: `int1e_ipipnuc` (d2/dA2),
  `int1e_ipnucip` (d2/dA dB), and per-nucleus `int1e_ipiprinv` / `int1e_iprinvip`
  (charge center). Feed into the validated `grd1::hess_en` block bookkeeping
  (the nine-block translational-invariance scatter is already correct; only the
  integrals it was fed were wrong).
- **2e (Task B):** `libint2_build_eri2` for `d2(uv|ls)`; contract
  `1/2 sum P_uv P_ls d2(uv|ls)` (Coulomb + exchange, scaled by HFscale for
  hybrids), paralleling `grd2_driver`.
- **Build:** regenerate `oqp-libint` with deriv order 2
  (`INCLUDE_ERI=2`, `INCLUDE_ONEBODY=2`); build with `USE_LIBINT=ON`. Rys and
  libint 2e engines coexist (runtime `libint2_active`), so energy/gradient are
  unchanged.
- **Effort:** medium. **Risk:** low (battle-tested integrals; the assembly
  bookkeeping is already validated). **Environment:** requires a libint-capable
  build machine (the libint tarball host is unreachable from the current
  sandbox).

### Option B — native Rys root-derivative extension (self-contained)

Extend `rys.F90` / `grd2_rys.F90` to provide `du_t/dX`, `dw_t/dX` and chain them
through the existing xyz-table recursion to form second derivatives.

- **Effort:** high (root-derivative recurrences, careful validation per angular
  momentum). **Risk:** medium-high (this is the part that already broke once with
  the fixed-root shortcut). **Environment:** buildable anywhere (no libint).
- **Value:** fully self-contained native Hessian; no libint runtime dependency.

### Option C — native Obara-Saika / McMurchie-Davidson second derivatives

Introduce a parallel auxiliary-integral scheme (Hermite/OS) for the 2nd
derivative, avoiding Rys-root derivatives.

- **Effort:** high (a second integral scheme alongside Rys). **Risk:** medium.
- **Value:** native; but duplicates integral capability rather than extending the
  existing Rys path.

**Recommendation:** **Option A** as the baseline to reach a working, validated
native Hessian fastest, with **Option B** as a later "no-libint" hardening goal.
The two share the entire assembly/validation layer, so Option A's work is not
wasted if Option B is pursued later.

---

## 4. Validation strategy (mirrors the response subsystem's rigor)

Every second-derivative integral term is validated bottom-up against an
exact/non-iterative reference, with O(h^2) scaling demonstrated, before it is
allowed into the final Hessian assembly. Two independent reference classes:

1. **Element-wise vs PySCF/libcint** (test-time oracle, not a runtime
   dependency): `int1e_ipipnuc`, `int1e_ipnucip`, `int1e_ipiprinv`,
   `int2e_ipip1`/`int2e_ipvip1`, contracted the same way as OpenQP. PySCF is
   already wired in the bridge.
2. **Finite difference of the validated first-derivative routines** (the
   `hess1_selftest` pattern): FD the production gradient/derivative-matrix
   routines (`grad_ee_*`, `grad_en_*`, `grd2_driver`) and compare to the analytic
   second derivative, contracting with a fixed matrix so basis normalization is
   handled consistently. The separated-FD trick (basis vs charge positions
   independently) localizes failures by sub-block.

Acceptance per term: max-abs and relative error reported; O(h^2) ratio shown
across >= 3 step sizes; truncation floor vs reference floor separated explicitly.
This is the same discipline that validated the response subsystem; do not relax
it for the integral terms.

---

## 5. Task breakdown and dependency graph

```
Task A  1e nuclear-attraction d2  (libint ipipnuc/ipnucip/ipiprinv)
            -> validate vs PySCF + FD of grad_en_*  -> feed grd1::hess_en
Task B  2e d2(uv|ls)              (libint eri2)
            -> validate vs PySCF + FD of grd2_driver -> new hess_2e driver
Task D  XC 2nd geometric deriv    (GGA/meta-GGA only; method-specific)

Final assembly (separate step, AFTER A+B validated):
   H = hess_nn + (sum P d2h - sum W d2S) + 1/2 sum PP d2(uv|ls)
       + [DFT: XC 2nd deriv] + [response term from U^x: validated subsystem]
   -> validate full H vs central FD of the analytic HF/DFT gradient on
      H2/H2O/formaldehyde; optional PySCF mf.Hessian() cross-check
   -> only then unguard source/modules/hf_hessian.F90
```

Tasks A, B, D are mutually independent and independent of the response
subsystem. The final assembly depends on all of them plus the (done) response
term. RHF Hessian needs A + B only (no XC). GGA adds D (LDA/GGA grid); meta-GGA
adds the tau/laplacian grid terms.

---

## 6. Method separation (RHF vs GGA vs meta-GGA)

| Component                  | RHF | GGA/RKS | meta-GGA |
|----------------------------|-----|---------|----------|
| nuclear-repulsion d2       | done (shared) | done | done |
| 1e overlap/kinetic d2      | done (shared) | done | done |
| 1e nuclear-attraction d2 (A) | needed (shared) | needed | needed |
| 2e d2 (B)                  | needed (shared) | needed (HFx-scaled) | needed |
| CPHF/CPKS response (U^x)   | done (shared) | done (fxc in A-matrix) | tau path to re-check |
| XC 2nd geometric deriv (D) | n/a | needed (hard) | needed (harder: tau, lapl) |
| functional-class gating    | n/a | low | low |

Everything except Task D is shared and method-independent. RHF is the natural
first complete native Hessian milestone (A + B + assembly, no XC).

---

## 7. Open questions to resolve before implementation

1. Confirm the regenerated `oqp-libint` exposes `INCLUDE_ONEBODY >= 2` (1-body
   deriv-2) in addition to `INCLUDE_ERI >= 2`, and the exact component ordering
   libint returns for the 9-component (3x3) blocks, so the `hess_en` scatter maps
   correctly.
2. Decide the 2e screening/permutational-symmetry handling for the deriv-2
   contraction (reuse `grd2_driver`'s Schwarz machinery vs libint's own).
3. Decide whether to keep `USE_LIBINT=OFF` builds able to compute the analytic
   Hessian at all (they cannot, by this design) -> they must fall back, with an
   explicit error, to the PySCF bridge or numerical Hessian (no silent fallback).
4. For Option B (later): enumerate the root-derivative formulas and the per-L
   validation plan.
