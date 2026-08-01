# MRSF-TDDFT analytic NAC Lagrangian — final handoff

## Status block (2026-08-01)

- Checkout: `sessions/20260731_nac_audit/repo_nac`
- Branch: `nac-lagrangian`
- Private publication target: remote `private`,
  `git@github.com:karmachoi/openqp-private.git`, branch `nac-lagrangian`
- Never push this private NAC work to `origin` (Alireza) unless the owner
  explicitly authorizes an upstream submission.
- H2O final-binary v94 one-call gate: maximum pair component error
  `2.30391559e-6`, **PASS** at `1e-5`.
- H2O final-binary v94 Lee diagonal gate: value error `6.62971067e-7` and
  native/legacy Z closure `4.20265769e-7`, both **PASS** at `5e-6`.
- C1 ethylene final-binary v94 one-call gate: maximum pair component error
  `7.55836774e-6`, **PASS** at `1e-5`.
- Acrolein final-binary v94 one-call gate: pair errors
  `(1,2)=2.13672772e-6`, `(1,3)=2.30452788e-6`, and
  `(2,3)=2.21028518e-6`; **PASS** at `2e-5` against the numerical `tlf=2`
  HST reference. The analytic metric remains the exact `ndtlf=0` identity
  derivative; there is no analytic `tlf=2` branch.
- Final Z-vector-only/log guard rebuild: `H2O_fort6_guard_v96.npz` reproduces
  the v94 H2O errors, logs six ordered-pair Z-vector solves and zero forward
  CPHF solves, and creates no `fort.6`.
- Publication status is defined by the live branch-HEAD equality below. The
  authoritative revision is the branch `HEAD`, not a prospective hash
  embedded in that same commit.

### Publication verification record

After the closeout commit is made, verify all three commands report the same
object:

```
git rev-parse HEAD
git rev-parse private/nac-lagrangian
ssh chc3 'git -C ~/nac_audit/repo rev-parse HEAD'
```

Fetch the private branch before the second check and bundle-sync the chc3
mirror before the third. Equality of those live branch heads is the
publication record. Do not insert a future commit's own hash into the file
that the commit is supposed to contain.

The old 7.52--7.53 conclusion that a fold-sector Fock-response term still had
to be derived is superseded. Reading the `socc/fa/fb`, `mrsfmntoia`,
`mrsfcbc`, `spc`, `trans`, and fold paths produced the closed resident
`mrsf_nac_wpair_impl`; the old antisymmetric-`dD` shortcut remains falsified
and must not be revived. The last large diagonal residual was then localized
to, and fixed in, the XC moving-grid derivative of the relaxed probe `P=T+Z`;
it was not caused by ROKS versus MRSF or by an omitted coordinate CPHF
calculation.

The intended production architecture is now:

1. Python validates scope and calls one C entry point.
2. Fortran owns the complete ordered-state-pair Lagrangian, exact metric,
   amplitude source, one-RHS adjoint Z-vector, coordinate contractions,
   antisymmetrization, and gap scaling.
3. Python only reshapes the final `dcv` and `nacv` records.

This is deliberately a **Z-vector implementation**. The routine named
`cphf_solve_rohf` is reused only as the native ROHF orbital-Hessian linear
solver. Production never constructs or solves the `3N` coordinate-response
block.

## Exact next-session entry sentence

"`HANDOFF_NAC_LAGRANGIAN.md`의 v94 H2O·ethylene·Acrolein 최종 통과를
읽고, 완료된 정적 singlet raw-electronic analytic NAC 범위를 유지한 채
`openqp-private` push와 chc3 branch-HEAD 일치만 확인해."

## Scientific anchor: Lee gradient Lagrangian and its NAC continuation

Primary reference: S. Lee, E. E. Kim, H. Nakata, S. Lee, and C. H. Choi,
*J. Chem. Phys.* **150**, 184111 (2019),
doi:[10.1063/1.5086895](https://doi.org/10.1063/1.5086895).

Use the paper's index sets `C` (doubly occupied), `O` (the two singly
occupied orbitals), and `V` (virtual):

- Eq. (3.3) is the MRSF excited-state Lagrangian. The C--V multiplier has
  the explicit leading factor two; C--O and O--V do not.
- Eq. (3.6) is `Jbar Zbar = -Rbar`.
- Eqs. (3.8)--(3.10) define the SF/ROHF orbital Hessian, the MRSF source,
  and `H+`. The left-hand operator is the SF/ROHF orbital Hessian; MRSF
  changes the source.
- Eq. (3.16) gives the relaxed density `P=T+Z`.
- Eq. (3.21) gives the full gradient structure
  `Omega^R = h^R:P - S^R:W + ERI^R:Gamma`.
- Published Eq. (3.8f) has a `-1/c_H` typo. Eq. (3.10), SI Eq. (S26),
  direct differentiation, and the OpenQP implementation require `-c_H/2`.

The paper proves the **state-diagonal energy-gradient** Lagrangian. It does
not by itself prove the interstate polarization of every quadratic MRSF
term. The off-diagonal continuation implemented here is therefore supported
by the exact-overlap derivation, the duplicated-slot Lee limit, adjoint
identities, and independent numerical NAC references; it must not be
described as a formula quoted directly from the paper.

For a ket column `J`, the ordered HST derivative is

```
Dord_IJ = Gmet_IJ^T X_J^R + gamma_IJ : (Sk^R + U^R),
d_IJ    = 1/2 (Dord_IJ - Dord_JI),
h_IJ    = (E_J-E_I) d_IJ.
```

On the physical folded response space,

```
(Omega_J-A) y_IJ = Q_J Gmet_IJ,
Q_J = 1-X_J X_J^T,
Gmet_IJ^T X_J^R = y_IJ^T A^R X_J.
```

The resident implementation splits the last bilinear derivative into

```
T1   = mrsf_nac_amp(y_IJ,X_J) + mrsf_nac_esum(y_IJ,X_J),
M_IJ = MT_frozen + MT_response + gamma_IJ,

Dord_IJ = T1
        + zeta_IJ^T B_HF/JK/Pulay^R
        + zeta_IJ^T B_XC^R
        + M_IJ : V^R
        + gamma_IJ : Sk^R.
```

`V^R` is the overlap-fixed symmetric/reorthonormalization response. The
independent ROHF rotations are eliminated by the adjoint equation, not by
coordinate-wise response:

```
H_native^T zeta_IJ = E^T (M_IJ-M_IJ^T).
```

### Z-vector, not coordinate CPHF

The production solve has exactly one RHS for each ordered state pair and
contracts `zeta^T B^R` analytically for all nuclear coordinates. It does not
solve `H U^R=B^R` independently for `3N` coordinates. The implementation
reuses `cphf_solve_rohf` because that routine already provides the correct
native symmetric-indefinite ROHF Hessian action and MINRES convergence
check. Its historical name does not change the algorithm into forward CPHF.
`source/modules/hf_hessian.F90` and forward CPHF records are diagnostic
oracles only and must stay out of the production path.

### Sign and normalization that must not change

In the native OpenQP independent-rotation coordinates,

```
H_native = 2 Jbar,
ell_pair = E^T(M_II-M_II^T) = +2 Rbar,
H_native zeta = ell_pair  =>  zeta = -Zbar.
```

The legacy diagonal gradient uses a different bookkeeping seam:

```
sfrorhs = -2 Rbar,
sfrolhs = Jbar,
xk      = 2 Zbar,
sfropcal inserts xk/2.
```

Therefore the diagonal convention closure is

```
zeta_native + xk_legacy/2 = 0.
```

Never feed native `zeta=-Zbar` through `sfropcal`, and never feed legacy
`xk=2Zbar` directly to a native adjoint contraction. The earlier handoff
statement `native z=Zbar` was wrong and is superseded by the duplicated-slot
source and Z-closure gates.

The native tangent and dual maps are also not ordinary array inverses:

```
unpack(kappa): SD -> beta; DV -> alpha and beta; SV -> alpha
pack(g):       SD <- beta; DV <- alpha + beta; SV <- alpha
pack(unpack(kappa)) = diag(1,2,1) kappa.
```

The DV factor is already present in both sides of the stationarity equation.
Do not add another `1/2` or `2`.

## Production: one-call resident Fortran path

The public entry point is
`include/oqp.h:mrsf_nac_lagrangian(struct oqp_handle_t *inf)`. The source is
`source/modules/mrsf_nac_driver.F90`, module `mrsf_nac_driver_mod`. The module
source glob in `source/modules/CMakeLists.txt` picks up the new `.F90` file.

The driver uses the **actual resident** `infos%tddft%nstate`, not a Python
configuration guess. It checks finite resident amplitudes and energies,
converged SCF and Davidson states, and a numerically resolved nonzero pair
gap. It currently hard-gates:

- ROHF/ROKS (`scftype==3`);
- restricted MRSF, not UMRSF;
- singlet target states;
- the two-SOMO reference (`N_alpha-N_beta=2`);
- `scf.conv` and TD convergence thresholds no looser than `1e-8`
  (`1e-10` is recommended near crossings).

The full resident traversal is:

1. Save `td_bvec_mo` and `int2e_cutoff`; set the NAC ERI cutoff to `1e-20`.
2. For each ket column `J`, build one exact streamed metric column with
   `mrsf_nac_metric_column`.
3. For each `I != J`, form `ytil=X_I/(E_J-E_I)` and zero the redundant
   folded coordinate.
4. Run `mrsf_nac_wpair_impl`, inject `ytil` only into the selected amplitude
   slot, run `mrsf_nac_amp`, `mrsf_nac_esum`, restore amplitudes, and run
   `mrsf_nac_response`.
5. Stream the current `gamma_IJ`, form the overlap/native source, solve the
   one-RHS `mrsf_nac_rohf_zvector`, then evaluate HF/JK/Pulay and XC
   adjoints.
6. Accumulate the ordered pair, then finalize all pairs by HST
   antisymmetrization and energy-gap scaling.
7. Restore amplitudes and the original ERI cutoff. The wrapper keeps the
   requested log unit open so a normal run does not create `fort.6`.

The streamed work records are `OQP::nac_ytil`, `OQP::nac_xstate`,
`OQP::nac_gamma_pair`, and `OQP::nac_rohf_z`; only the final
`OQP::nac_dcv`/`OQP::nac_nacv` tensors cross back to Python in production.

### TagArray pointer lifetime and dimension hardening

A TagArray `remove_records` or `reserve_data` operation may relocate the
entire backing store and invalidate pointers to records that were not
themselves changed. The final driver is hardened around that ownership rule:

- `mrsf_nac_lagrangian` copies both the resident amplitudes and state energies
  into owned `bvec_saved`/`energies_saved` arrays before it creates streamed
  records. Pair gaps are always read from the owned energy copy.
- `mrsf_nac_amp` owns copies of its alpha/beta densities, MO coefficients,
  and amplitude matrix before reserving `OQP::nac_amp` or entering nested
  gradient routines.
- `mrsf_nac_esum` likewise owns its amplitude, alpha/beta MO, and density
  inputs before repeatedly reserving pair and diagnostic outputs.
- `mrsf_nac_rohf_pair_overlap` reacquires the unchanged MO record after its
  final output reservation. `mrsf_nac_pair_finalize` reacquires both the
  ordered-pair tensor and state energies after creating `dcv` and `nacv`.

The one-call driver performs dimension arithmetic in 64-bit integers before
converting to default Fortran integers. It guards `3*natom`, `nstate^2`,
`ncoord*nstate^2`, `noca*(nbf-nocb)`, `nbf^2`, and the resident
state-amplitude/metric products against overflow and TagArray dimension
limits. Invalid orbital occupations and inconsistent resident tensor shapes
abort before allocation or reshape.

### Source map

| Source | Production responsibility |
|---|---|
| `source/modules/mrsf_nac_driver.F90` | One-call ordered-pair loop, scope/state guards, streaming, restoration, finalization |
| `source/modules/mrsf_nac_metric_data.F90` | Exact `tlf=0` metric oracle and `O(nstate*nbf^2)` production column streamer |
| `source/modules/tdhf_mrsf_gradient.F90` | `mrsf_nac_amp`, `mrsf_nac_esum`, closed analytic `mrsf_nac_wpair_impl`, standard-gradient moving-grid correction |
| `source/modules/tdhf_mrsf_energy.F90` | Resident MRSF response/source assembly |
| `source/modules/mrsf_nac_interchange.F90` | Pair overlap/source, one-RHS Z-vector, HF/JK/Pulay and XC adjoints, accumulation/finalization |
| `source/modules/cphf.F90` | Reused native ROHF Hessian action, tangent/dual maps, MINRES; not a production `3N` solve |
| `source/dftlib/dft_gridint_tdxc_grad.F90` | Linear-probe XC derivative and probe-only partition/owner moving-grid response |
| `pyoqp/oqp/library/nac_analytic.py` | Thin scope wrapper, exactly one `oqp.mrsf_nac_lagrangian(mol)`, optional observational debug save, output reshape |
| `include/oqp.h` | C declarations for the one-call driver and diagnostic component entries |

`nac_analytic.py` must not regain state-pair scientific loops. Future OpenQP
numerical kernels are Fortran-first; Python is for API orchestration,
validation, and artifact inspection.

## Exact metric streaming

`mrsf_nac_metric_data` is the full diagnostic oracle and may materialize all
pair metric derivatives at `O(nstate^2*nbf^2)` storage. Production instead
calls

```
mrsf_nac_metric_column(infos, jstate, gamma_column)
```

once for a fixed normalized ket column and consumes each `I` immediately.
Its workspace is `O(nstate*nbf^2)`.

For exact `ndtlf=0`, the implementation differentiates the seven determinant
blocks used by the state-overlap code. It uses the same physical MRSF fold as
`mrsfxvec`, zeros the special 2-by-2 SOMO coefficients in the generic blocks,
and constructs the `sij`, `sab`, and `sia` minor families. The normalized
column derivative includes

```
delta_IK/n_J - R_IJ R_KJ/n_J^3.
```

Direct identity-minor cofactors are required; inverse/Jacobi formulas are
invalid because some relevant minors are singular. The orbital source stores
the half-antisymmetric derivative in opposite slots,
`0.5*(g_pq-g_qp)`. Because the overlap normalization is one-sided by ket
column, `gamma_IJ` must **not** be assumed equal to `-gamma_JI`; HST
antisymmetrization happens only after both ordered kernels have been built.

Primary checks are `tests/test_mrsf_nac_metric_data.py`,
`tests/test_nac_gamma_resident_gate.py`, and
`tests/test_mrsf_nac_fortran_driver.py`. The closed-metric artifacts are
`H2O_closed_metric_v81.npz` and `H2O_closed_metric_v81_debug.npz` on chc3.

## XC moving-grid correction and the large-error diagnosis

The large residual was not evidence that ROKS is incompatible with MRSF and
was not repaired by adding coordinate CPHF. A finite atom-centred quadrature
for a relaxed linear probe is

```
Q_XC[P;R] = sum_g w_g(R) q_P(r_g(R);R),     P=T+Z.
```

Its derivative contains all of

```
sum_g w_g dq_P/dR
+ sum_g q_P dw_g/dR
+ owner-point motion of q_P.
```

The pair `mrsf_nac_esum(T)` and `mrsf_nac_xc_adjoint(Z)` already carried
their moving-grid pieces. The ordinary diagonal MRSF gradient did not carry
the corresponding relaxed-probe `P=T+Z` weight/owner derivative. That omission
was why the full duplicated-slot Lee value did not close even after the
source and Z-vector were correct.

The standard gradient now retains its established ground-plus-probe XC call
and adds a second correction-only resident sweep with

```
include_ground_state=.false.
include_weight_derivative=.true.
weight_derivative_only=.true.
```

The second sweep uses relaxed `P=T+Z`, omits the ordinary AO/basis-gradient
part, and adds only the partition-weight and owner-point correction. A naive
single combined call is mathematically mismatched: it adds ground-state owner
motion while the partition probe contains only `P`. The deliberately failed
v89 experiment produced a state-independent error of about `0.533282` and
must never be used as a shortcut.

The implemented moving-grid branch is currently restricted to a single
linear probe (`nMtx=1`, `doFxc=false`). It aborts for the unsupported `xa/xb`
third-derivative branch rather than silently returning an incomplete answer.

### Decisive duplicated-slot gates

- `H2O_diagonal_rhs_all_v71.npz`: `ell_pair + rhs_legacy` maximum absolute
  errors `3.22e-8`, `1.60e-7`, `3.93e-7`; relative L2
  `5.15e-8`, `4.25e-6`, `4.87e-6`; fitted sign/scale is `-1` for all roots.
  This proves the aggregate `y=X` Eq. (3.9) source and the native
  `zeta=-Zbar` convention.
- v84: Z closure passed, but the full value still failed at
  `2.894e-5`, `8.432e-6`, `1.224e-5`.
- v85: direct `-xk/2` re-contraction left the same value error, excluding the
  Z solver and exposing amplitude/esum/Z-HF/Z-XC/overlap components.
- v86: `1e-12` SCF/TD/Z thresholds did not change the residual.
- v87, no esum grid: residuals fell to
  `7.0430e-6`, `1.43475e-6`, `2.33657e-6`.
- v88, pure HF: `3.9063e-6`, `6.37e-8`, `3.58e-7`, localizing the remaining
  DFT term.
- v89: the intentionally wrong combined ground-plus-probe moving-grid call;
  error about `0.533282`.
- `H2O_diagonal_value_v90_gridfixed.npz`: correct probe-only moving-grid
  sweep. Full-value errors are `6.6297e-7`, `3.6219e-8`, `6.6150e-8`;
  native/legacy Z closure is `2.89e-8`, `2.51e-7`, `4.20e-7`; pair and
  production translational sums are at about `1e-14`.

The `T` moving-grid weight accounts for 95--99% of the squared pre-v90
residual; the full correction has cosine `0.999966`, `0.999979`, `0.999993`
with that residual. v90 establishes the complete Lee diagonal continuation
on this H2O case. It does not, by itself, certify the final one-call
off-diagonal driver on other molecules; the v94 H2O, ethylene, and Acrolein
gates in the next section now supply that separate evidence.

## Final-binary v94 numerical evidence

The current one-call Fortran binary has now been judged against independent,
separately reconverged tight references rather than only against its
predecessor implementation:

| Final artifact | Quantity | Result | Gate |
|---|---|---:|---:|
| `H2O_final_fortran_v94.npz` | maximum off-diagonal pair component error | `2.30391559e-6` | **PASS**, `<=1e-5` |
| `H2O_diagonal_value_v94_final.npz` | maximum Lee diagonal value error | `6.62971067e-7` | **PASS**, `<=5e-6` |
| `H2O_diagonal_value_v94_final.npz` | maximum native/legacy Z closure | `4.20265769e-7` | **PASS**, `<=5e-6` |
| `ETH_final_fortran_v94.npz` | maximum off-diagonal pair component error | `7.55836774e-6` | **PASS**, `<=1e-5` |
| `Acrolein_final_fortran_v94.npz` | pair `(1,2)` component error | `2.13672772e-6` | **PASS**, `<=2e-5` |
| `Acrolein_final_fortran_v94.npz` | pair `(1,3)` component error | `2.30452788e-6` | **PASS**, `<=2e-5` |
| `Acrolein_final_fortran_v94.npz` | pair `(2,3)` component error | `2.21028518e-6` | **PASS**, `<=2e-5` |

The H2O, C1 ethylene, and Acrolein results close the final-driver state-pair
ownership, streaming, TagArray lifetime, and output-reshape question across
the symmetric H2O case, a C1 small molecule, and Acrolein judged against a
numerical `tlf=2` HST reference. In every case the analytic production metric
is the exact `ndtlf=0` identity derivative described above. The v94 diagonal
result reproduces the v90 moving-grid closure with the same final binary and
confirms that the one-RHS native Z-vector remains Lee-normalized after the
pointer-lifetime hardening.

For historical comparison, the pre-one-call Fortran-kernel/Python-pair
orchestrated path had passed the following `2.5e-4 Angstrom`,
`scf.conv=tdhf.conv=1e-10` references:

| Historical artifact | Pair errors (max component) |
|---|---|
| `H2O_production_optimized_v65.npz` | `4.59e-8`, `2.15e-7`, `1.13e-6` |
| `ETH_production_movinggrid_v60.npz` | `2.71e-6`, `1.25e-7`, `7.07e-6` |
| `Acrolein_production_movinggrid_v62.npz`, `tlf=2` | `2.14e-6`, `2.31e-6`, `2.22e-6` |

Those older freezes remain useful scientific referees, but all three final
completion claims now rest on v94.

## v3--v19 forensic gate map

These scripts document how the old error was localized. They are referees,
not the current production API.

| Campaign | Files | Purpose / surviving conclusion |
|---|---|---|
| v3 | `assembly_v3.py`, `assembly_v3c.py`, `v3d`--`v3h` | Master exact-overlap decomposition, PT amplitude response, analytic `Sk`, exact displaced `w_ref`; H2O and C1 structure established |
| v4 | `v4a_adjoint.py` | Machine-level identity between direct response and amplitude adjoint |
| v5 | `v5_prod_h2o.py`, `v5b_h2o.py`, `v5c_probe.py` | First production-form source, polarization and record-scaling diagnostics; state-specific gradient-chain polarization rejected |
| v6 | `v6_h2o.py`, `v6b_h2o.py` | Split skeleton, W, Z/P, and seam channels |
| v7a--v7d | corresponding scripts | Symmetric channel, same-space response, calibrated seam, and full-record restoration |
| v7e--v7i | corresponding scripts | Unit-L seam calibration, V-mask convention, full `MT` harvest, direct-injection and two-step FD checks |
| v7j | `v7j_h2o.py` | Unified one-process referee; never cross process-local gauges |
| v7k--v7o | `v7k`, `v7l`, `v7m`, `v7o` | Fock/Jacobi record dependence, DM response, G[P] linearity/self-adjointness, consolidated freezes |
| v8 | `v8_wprobe.py`, `v8b_wprobe.py`, `v8c_wprobe.py` | Slot-resolved `w`/fold response and fixed-reference gauge |
| v9 | `v9_gvec.py` | Measured Fock-response `G[dD]` channel |
| v10 | `v10_probe.py` | Continuation of fixed-reference slot diagnostics |
| v11 | `v11_branch.py` | One-sided branch/kink check |
| v12 | `v12_tscan.py` | Staged-response step-size scan |
| v13 | `v13_final.py` | Coefficient reconstruction and mixed-order error scaling |
| v14 | `v14_swap.py` | Same-state versus restored-state bisection |
| v15 | `v15_emap.py` | Entrywise reconstruction-error scaling map |
| v16 | `v16_convcheck.py` | Located the overlap/`Sk` convention mismatch |
| v17 | `v17_fix.py` | Fixed the second-call `Sk` diagonal by rebuilding `C0^T S_cross C0` |
| v18 | `v18_h2o.py` | Re-ran the unified gate with the `Sk` fix; closed the U channel to FD level |
| v19 | `v19_fchan.py` | Separate-process F-channel referee; isolated the fold-sector residual that was subsequently closed in analytic Fortran `mrsf_nac_wpair_impl` |

`v20_fold_audit.py` and `v21_production_gate.py` are later bridge/audit
scripts. The current decisive API tests live under `tests/`; do not extend the
historical Python campaign into a new production implementation.

## Frozen artifacts and worker directories

Canonical remote artifact root is chc3 `~/nac_audit/probe/`. The repo mirror
is `~/nac_audit/repo`.

### Current decisive artifacts

- Final-binary v94 gates:
  `H2O_final_fortran_v94.npz`,
  `H2O_diagonal_value_v94_final.npz`, and
  `ETH_final_fortran_v94.npz` are passed freezes with the exact errors listed
  above. `~/nac_audit/probe/Acrolein_final_fortran_v94.npz` is also passed;
  its three pair errors are `2.13672772e-6`, `2.30452788e-6`, and
  `2.21028518e-6` at the `2e-5` gate.
- Prior tight off-diagonal production referees:
  `H2O_production_optimized_v65.{npz,out}`,
  `ETH_production_movinggrid_v60.{npz,out}`,
  `Acrolein_production_movinggrid_v62.{npz,out}`.
- Intermediate resident/one-call development freezes:
  `H2O_production_fortranfirst_v74.{npz,out}`,
  `H2O_production_fortranfirst_v76.{npz,out}`,
  `H2O_resident_finalize_v80.{npz,log}`,
  `H2O_closed_metric_v81{,_debug}.npz`.
  These explain the route to v94 but no longer define the final H2O verdict.
- Diagonal source/value campaign:
  `H2O_diagonal_rhs_all_v71.npz`,
  `H2O_diagonal_value_v84.npz`, `v85`, `v86`, `v87_nogrid`, `v88_hf`,
  deliberately failed `v89_gridfixed`, and decisive
  `H2O_diagonal_value_v90_gridfixed.npz`. The final-binary successor is
  `H2O_diagonal_value_v94_final.npz`.
- Tight numerical references:
  - `H2O_energy_tlf0_tight_dx25e5_dnum.npz`, SHA256
    `6c889520914e61c2ad2e5687875d3496471f87513424391645dfed9b51a8e450`,
    worker `~/nac_audit/probe/H2O_energy_tlf0_tight_num_nacv/`;
  - `ETH_energy_tight_dx25e5_dnum.npz`, SHA256
    `5bd3e1fcc1b796176be68b22a4cd030b0076d8c7cd6942757114fb77d372a086`,
    worker
    `~/nac_audit/repo/tools/nac_lagrangian/ETH_energy_tight_dx25e5_num_nacv/`;
  - `Acrolein_S2_tlf2_dx25e5_dnum.npz`, SHA256
    `d5a8eac3ef9c83231f5180efd8188a00a780b3d76629f476c1837c3f62532d55`,
    worker
    `~/nac_audit/repo/tools/nac_lagrangian/Acrolein_S2_tlf2_dx25e5_num_nacv/`.
  The old `H2O_energy_tlf0_tight_v36_z.npz` has a stale embedded reference;
  compare its `dcv` only against the separate tight file.

### Historical v3--v19 freezes

- Session parent and/or chc3:
  `H2O_energy_tlf0_{v7h,v7i,v7o}.npz`,
  `ETH_energy_{v7o,ctx,dnum}.npz`, and
  `H2O_energy_tlf0_dnum.npz`.
- The early campaign also cited `h2o_nac_reference.npz`,
  `h2o_nac_ref_tlf0.npz`, `ladderA8_data.npz`, and
  `gamma_formula_h2o.npz`. Only `h2o_nac_reference.npz` is present in the
  current session parent; verify the physical file before depending on the
  other three historical names.
- Verified chc3 worker directories:
  `~/nac_audit/probe/H2O_energy_tlf0_skel/` and
  `~/nac_audit/probe/ETH_energy_skel/`.
- Verified numerical displacement directories include
  `H2O_energy_tlf0_dx25e5_num_nacv/`,
  `H2O_energy_tlf0_dx5e4_num_nacv/`, and tight H2O variants.
- Worker amplitudes, state phases, and displaced-frame vectors are
  process-local. Recreate `ref.npz` in the same worker/gauge cycle; never
  splice vector arrays from independent processes just because filenames
  match.

## Completed stages

1. Fixed the exact `tlf=0` diagonal-minor path and established the signed HST
   `dcv`/`nacv` conventions against independent overlap derivatives.
2. Derived and gated the master exact-overlap decomposition, normalized
   metric derivative, amplitude adjoint, and one-sided state gauge.
3. Read the MRSF sigma implementation and replaced the old fold/F-channel
   hypothesis with closed analytic `mrsf_nac_wpair_impl`; the O(`nbf^2`)
   orbital-generator FD harvest is gone from production.
4. Established the native tangent/dual metric and the one-RHS ROHF/ROKS
   Z-vector convention; fixed DFT initialize/clean lifecycle and eliminated
   coordinate CPHF from production.
5. Moved exact metric work to resident Fortran and reduced production metric
   workspace to `O(nstate*nbf^2)` by column streaming.
6. Moved the ordered-pair traversal and final assembly to the single Fortran
   `mrsf_nac_lagrangian` entry; reduced Python to one call and reshape.
7. Closed the duplicated-slot Lee source and Z sign with v71, then closed the
   full H2O diagonal value with the probe-only `P=T+Z` moving-grid fix in v90.
8. Hardened every final-driver TagArray lifetime boundary with owned
   amplitude/esum inputs, copied driver amplitudes and energies, and explicit
   pointer reacquisition in pair-overlap and finalization. Added 64-bit
   dimension-product and occupation/shape guards before allocation.
9. Rebuilt the final one-call binary and passed H2O off-diagonal, H2O full
   Lee diagonal/Z closure, C1 ethylene, and the Acrolein v94 gate against its
   numerical `tlf=2` HST reference at the criteria recorded above. The
   analytic metric is exact `ndtlf=0` in all four gates.

## Original five-stage analytic-NAC closeout

| Original stage | Status | Closure evidence |
|---|---|---|
| 1. Derive the fold/F-channel term | **COMPLETE** | Sigma-source audit and closed analytic `mrsf_nac_wpair_impl`; the 7.53 antisymmetric-`dD` shortcut remains falsified |
| 2. Propagate, rejudge, and refreeze | **COMPLETE** | Final one-call H2O, diagonal, ethylene, and Acrolein v94 freezes pass |
| 3. Fill the resident wpair implementation | **COMPLETE** | No production orbital-generator FD harvest; resident analytic source is called by the one-call driver |
| 4. Validate Acrolein against the numerical `tlf=2` HST reference | **COMPLETE** | Exact-`ndtlf=0` analytic v94 pair errors `2.13672772e-6`, `2.30452788e-6`, `2.21028518e-6`; pass criterion `2e-5` |
| 5. Audit the duplicate `tlf=0` fix before upstreaming | **COMPLETE for this private closeout** | `e5eeed07` on `private/main` and this branch's historical `2242bff` implement the same principal-minor correction; the randomized comparison differed by zero |

All five stages are complete. The original task is therefore
**COMPLETE: static singlet MRSF-TDDFT raw-electronic analytic NAC**. Do not
keep it artificially open for publication mechanics or the extensions below.

## Publication operation — separate from scientific completion

To publish or re-verify the result, review the intended dirty-tree scope,
commit and push to `openqp-private`, bundle-sync `~/nac_audit/repo`, and
execute the three branch-HEAD checks in the publication record. This branch
operation is required before saying a particular revision is published, but
it is not a missing derivation, code path, or numerical gate.

## Future extensions — not blockers to the static analytic-NAC closeout

- Consecutive-frame Acrolein SOC-NAMD validation with one causal
  state/subspace gauge applied consistently to NAC, SOC, gradients, and
  coefficients.
- An explicitly derived ETF-corrected coupling mode. The present result is
  raw electronic NAC and obeys its nonzero raw translation identity.
- Triplet, UMRSF, arbitrary/non-two-SOMO open-shell references, and
  exact-degenerate subspace response.
- The `doFxc`/third-functional-derivative moving-grid branch and multi-probe
  grid batching.
- Additional term-resolved diagnostics such as
  `grid(P)=grid(T)+grid(Z)`, broader method/basis coverage, and performance
  optimization. These can strengthen or extend the result but are not
  missing terms in the gated singlet raw-electronic formula.

## Landmine registry

1. **Z-vector versus CPHF:** production is one adjoint RHS per pair.
   `cphf_solve_rohf` is only the reused Hessian solver; a `3N` forward block
   is a diagnostic and is neither necessary nor acceptable in production.
2. **Two Hessians:** the MRSF Davidson matrix `A` is not the ROHF orbital
   Hessian `Jbar`/`H_native`.
3. **Two Z conventions:** native computational `zeta=-Zbar`; legacy gradient
   `xk=2Zbar` followed by `xk/2`. Mixing them produces a sign/factor error.
4. **Tangent/dual metric:** `pack(unpack)=diag(1,2,1)`, not identity. The DV
   factor two is not missing.
5. **TagArray layout:** Fortran-created non-square 2-D records arrive in
   NumPy with the Fortran storage interpretation; recover `(ltot,ncart)` as
   `raw.reshape(ncart,ltot).T`. Python-created and Fortran-created records do
   not have one universal transpose rule.
6. **Actual state count:** use resident `infos%tddft%nstate` and Fortran-first
   state indexing. Never reshape outputs from a requested/configured state
   count without checking what Davidson actually produced.
7. **State gauge:** Davidson phases are process-random. Solve one phase per
   state and enforce the cycle product; never choose signs independently per
   pair.
8. **Process-local vectors:** `Ux`, `w_ref`, amplitudes, and displaced-state
   frames cannot be mixed across processes. Scalar contractions may be
   compared after gauge resolution; raw arrays may not.
9. **Separate-process workers:** DM/Fock perturbation and `G[P]` referees need
   fresh worker processes. The in-process one-iteration SCF emulation has a
   known `trueG` blow-up and is not a valid referee.
10. **`Sk` diagonal bug:** the second same-geometry
    `overlap_mo_non_orthogonal` record has a bad diagonal. Rebuild
    `Sk=C0^T S_cross,AO C0` from the AO overlap record.
11. **DFT lifecycle:** every `dft_initialize` must be paired with `dftclean`.
    Otherwise functionals append to process-global state and XC can double.
12. **Moving-grid pairing:** partition-weight and owner-point terms must be
    differentiated together. The v89 combined ground-plus-probe shortcut is
    wrong; use the probe-only correction sweep. Production hard-wires the
    required weight derivative on; never add an environment kill switch that
    silently changes the scientific result.
13. **Moving-grid scope:** the present implementation is only for one linear
    probe with `doFxc=false`. Abort unsupported branches.
14. **Exact metric direction:** normalized overlap is one-sided by ket
    column. Do not impose `gamma_IJ=-gamma_JI` before ordered assembly.
15. **Singular minors:** exact `tlf=0` metric derivatives need direct
    identity-minor cofactors; inverse/Jacobi formulas fail at singular
    minors.
16. **Fold conventions:** `mrsfxvec` produces the physical `sqrt(2)` unfolded
    amplitudes, while the response-space fold uses `1/sqrt(2)`. These are
    adjoint roles, not interchangeable arrays. Always zero the redundant
    folded coordinate in `ytil`.
17. **ERI cutoff:** the response/matvec path needs `int2e_cutoff=1e-20` for
    column-exact linearity. Save and restore the user's cutoff.
18. **Log unit:** resident kernels called after the energy driver closes `IW`
    can create `fort.6`. Keep the requested log open around the one-call
    driver and gate that no normal path leaves `fort.6`.
19. **Loose stationarity:** a displacement-converged result at
    `scf.conv=tdhf.conv=1e-6` can still be wrong near a small state gap. Use
    `1e-10` for certification; tightening only the Z solve cannot repair
    unconverged reference orbitals.
20. **Exact degeneracy:** `y=X_I/(E_J-E_I)` is singular at an exact or
    numerically unresolved gap. Abort; a degenerate-subspace theory is a
    separate derivation.
21. **Raw translation is not zero:** without ETF,
    `sum_A d_IJ^A = antisym sum_A gamma_IJ:Sk_A`, not zero. Never subtract an
    atomic mean and call it ETF.
22. **Paper typo and scope:** use `-c_H/2`, not published `-1/c_H`; the Lee
    gradient paper anchors the diagonal limit but does not alone prove the
    interstate formula.
23. **Stale freezes:** v36 embeds the wrong reference, v89 is deliberately
    wrong, and historical vector freezes are gauge/process dependent. Artifact
    names are not proof; always state the exact comparison.
24. **Duplicate `tlf=0` fixes:** `e5eeed07` is already on `private/main`, while
    this branch carries historical `2242bff`. They implement the same
    principal-minor correction and the randomized comparison differed by
    zero. When integrating with `private/main`, drop/rebase the functional
    `2242bff` patch instead of stacking it; preserve only useful comments.
25. **TagArray pointer lifetime:** any reserve/remove can relocate the entire
    arena, not only the named record. Never retain resident pointers across a
    reservation; copy long-lived scientific inputs into owned arrays or
    reacquire unchanged records after the final reservation. This was a real
    final-driver hardening issue in amp, esum, pair overlap, and finalize.
26. **Dimension-product overflow:** validate occupations and compute all
    `3*natom`, `nstate^2`, orbital-matrix, amplitude, metric, and final-output
    products in a wider integer kind before default-integer conversion. A
    successful small-factor check does not prove the multiplied TagArray size
    is representable.
27. **Production environment switches:** do not put coordinate-FD self-tests or
    scientific kill switches inside `mrsf_nac_lagrangian` or its resident
    kernels. The hidden Python forward-CPHF gate, the esum moving-grid disable
    switch, and the in-routine `/tmp` FD self-test were removed at closeout.
    Observational record exports may remain only when they cannot change the
    returned coupling or mutate geometry-dependent state.

## Honest scope and completion claim

What is established now:

- the Lee-consistent diagonal source, native Z sign, and final-binary full H2O
  diagonal value including relaxed XC moving-grid response;
- final one-call off-diagonal accuracy on H2O, C1 ethylene, and Acrolein,
  where the Acrolein exact-`ndtlf=0` analytic result is judged against an
  independent tight numerical `tlf=2` HST reference in v94;
- a Fortran-resident exact metric, closed MRSF source, one-RHS ROHF/ROKS
  adjoint, analytic coordinate contractions, pointer-safe record ownership,
  overflow guards, and one-call driver architecture;
- completion of the original static singlet raw-electronic analytic-NAC task.

What remains outside that completed scope:

- triplet, UMRSF, arbitrary open-shell references, exact-degenerate
  subspaces, or `doFxc` moving-grid response;
- a trajectory-continuous NAC/SOC/gradient/coefficient gauge and a complete
  production SOC-NAMD shakedown;
- an ETF-corrected coupling mode;
- publication is an external branch-state assertion and must be read from the
  live private/chc3 branch-HEAD checks rather than an embedded self-hash.

Therefore the correct closeout statement is: **static singlet MRSF-TDDFT
raw-electronic analytic NAC is complete and final-binary gated on H2O, C1
ethylene, and Acrolein; publication status is determined operationally from
the live branch heads, while NAMD, ETF, triplet/UMRSF, `doFxc`, and
degenerate-subspace work are explicit future extensions rather than
blockers.**
