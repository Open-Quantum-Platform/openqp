# Gate 3 MRSF Hessian design

This document records the Gate 3A inventory and conservative implementation plan for MRSF-TDDFT Hessians in OpenQP. It is intentionally strict: an incomplete analytical MRSF Hessian must never be enabled silently, and the generic state-index-only numerical Hessian path is not considered safe for MRSF until Gate 3B root tracking is implemented.

## Current support matrix

| Request | Status | Rationale |
| --- | --- | --- |
| `method=hf`, `runtype=hess`, `hess.type=analytical`, `hess.state=0` | Supported scaffold via explicit external PySCF bridge | Gate 2 HF/DFT bridge is guarded and labeled; no numerical fallback is used for analytical requests. |
| `method=tdhf`, `tdhf.type=mrsf`, `runtype=hess`, `hess.type=analytical` | Unsupported and rejected | Native `tdhf_mrsf_hessian` is only a C ABI/build scaffold and aborts. Required MRSF Lagrangian/response terms are not assembled. |
| `method=tdhf`, `tdhf.type=mrsf`, `runtype=hess`, `hess.type=numerical` | Unsupported in Gate 3A and rejected | The existing generic numerical Hessian driver displaces geometries and requests the same state index. For MRSF, this can silently follow the wrong electronic character near root crossings or near-degeneracies. Gate 3B will add a root-tracked finite-difference oracle. |
| `tdhf.type=umrsf`, MRSF-style Hessian | Unsupported and rejected | UMRSF energy has a wrapper, but gradient/Z-vector/Hessian support is incomplete for a validated Hessian oracle. |
| `tdhf.type=sf`, analytical Hessian | Unsupported and rejected separately from MRSF | SF has distinct code paths and should not be routed through the private MRSF scaffold. |
| Ordinary TDDFT/TDA/RPA analytical Hessian | Unsupported and rejected | TDDFT analytical Hessian response terms are not implemented. |

## Existing Python data flow

1. `SinglePoint.energy()` in `pyoqp/oqp/library/single_point.py` computes the reference SCF energy and then dispatches TDHF excitations through `SinglePoint.excitation()` when `input.method=tdhf`.
2. `SinglePoint.__init__` maps `tdhf.type=mrsf` to `oqp.tdhf_mrsf_energy` and `tdhf.type=umrsf` to `oqp.tdhf_umrsf_energy`.
3. `SinglePoint.excitation()` calls the selected Fortran TDHF/MRSF energy function, then builds total state energies as `ref_energy[0] + OQP::td_energies` and stores them in `mol.energies`.
4. `Gradient.tddft_grad()` loops over requested `properties.grad` roots, calls `mol.data.set_tdhf_target(i)`, solves the matching Z-vector, then calls the matching gradient routine. For `tdhf.type=mrsf`, this dispatches `oqp.tdhf_mrsf_z_vector` followed by `oqp.tdhf_mrsf_gradient`.
5. `Hessian.numerical_hess()` currently uses central finite differences of gradients, spawning displaced gradient calculations through `grad_wrapper()`. This generic path rewrites each displaced calculation to `runtype=grad`, sets `properties.grad = hess.state`, and reuses the original state index without any MRSF-specific root tracking. Gate 3A disables that path for MRSF/UMRSF.
6. `Hessian.analytical_hess()` routes `tdhf.type=mrsf` and `tdhf.type=umrsf` to `analytical_mrsf_hess()`, which is a clear `NotImplementedError` guard. It reports OpenQP root mapping for MRSF (`root 1 -> physical S0`, `root 2 -> physical S1`, etc.) and `<S^2>` when `OQP::td_s2` is available.

## Existing Fortran MRSF energy flow

Primary files:

- `source/modules/tdhf_mrsf_energy.F90`
- `source/modules/tdhf_mrsf_z_vector.F90`
- `source/modules/tdhf_mrsf_gradient.F90`
- `source/modules/tdhf_mrsf_hessian.F90`

Energy flow in `tdhf_mrsf_energy.F90`:

1. C wrappers expose `tdhf_mrsf_energy` and `tdhf_umrsf_energy`.
2. The routine requires triplet reference multiplicity (`mol_mult=3`). Non-UMRSF MRSF expects ROHF (`scftype=3`); UMRSF requires UHF (`scftype=2`).
3. `mrst = infos%tddft%mult` selects MRSF spin manifold:
   - `mrst=1` or `mrst=3`: spin-adapted spin-flip/MRSF blocks with dimension `nocca*nvirb`, with special truncation of available roots near the spin-adapted reference sector.
   - `mrst=5`: beta-to-alpha spin-flip block with dimension `noccb*nvira`.
4. Davidson/Jacobi TD equations produce excitation energies and vectors.
5. Important stored records include:
   - `OQP::td_energies` (`OQP_td_energies`): excitation energies.
   - `OQP::td_bvec_mo` (`OQP_td_bvec_mo`): selected TD/MRSF eigenvectors in MO excitation space.
   - `OQP::td_t` (`OQP_td_t`): packed transition-density-like intermediate.
6. Spin/transition diagnostics are computed with routines such as `get_mrsf_transition_density`, `get_mrsf_transitions`, `get_transition_dipole`, and `get_spin_square`/`umrsfssqu` as applicable. Python-side guards can only report `<S^2>` if the value is exported to an accessible tag such as `OQP::td_s2`.
7. The Fortran energy routine records `infos%mol_energy%excited_energy = mrsf_energies(infos%tddft%target_state)`.

## Existing MRSF Z-vector and gradient flow

Z-vector flow in `tdhf_mrsf_z_vector.F90`:

1. `tdhf_mrsf_z_vector` reads `infos%tddft%target_state` and clamps it with `min(infos%tddft%target_state, infos%tddft%nstate)` in the current implementation.
2. Required records include reference Fock/MO data plus `OQP_td_bvec_mo`, `OQP_td_t`, and `OQP_td_energies`.
3. Allocated records for gradient consumption include:
   - `OQP_WAO`
   - `OQP_td_mrsf_density`
   - `OQP_td_p`
   - `OQP_td_abxc`
4. The target eigenvector `bvec_mo(:, target_state)` drives response-density and Z-vector construction.
5. The solver sets `infos%mol_energy%Z_Vector_converged` false/true based on convergence.

Gradient flow in `tdhf_mrsf_gradient.F90`:

1. `tdhf_mrsf_gradient` requires triplet multiplicity and reads `infos%tddft%target_state`.
2. The one-electron gradient is reused from spin-flip infrastructure through `sf_1e_grad(infos, basis)`.
3. It reads:
   - `OQP_DM_A`, `OQP_DM_B`: reference densities.
   - `OQP_td_abxc`: response/AB/XC intermediate.
   - `OQP_td_p`: packed relaxed/state density blocks.
   - `OQP_td_mrsf_density`: seven MRSF spin-coupling density blocks for `mrst=1` or `mrst=3`.
4. It builds AO-space density blocks `d`, `p`, and `spc` and accumulates:
   - DFT XC gradient with `utddft_xc_gradient` when `hamilton == 20`.
   - two-electron MRSF terms through `mrsf_2e_grad` for `mrst=1/3`.
   - SF-like two-electron terms through `sf_2e_grad` for `mrst=5`.
5. The final gradient is written through the normal `infos%atoms%grad`/`mol.get_grad()` path.

## Current Hessian scaffold

`source/modules/tdhf_mrsf_hessian.F90` exposes a C ABI symbol `tdhf_mrsf_hessian`, but its only runtime behavior is an aborting message:

`Analytic MRSF-TDDFT Hessian kernel scaffold reached; implementation awaits validated MRSF gradient/Z-vector finite-difference baseline.`

This is a build/link scaffold only. It must not return zeros, stale data, generic HF/DFT Hessians, or a numerical fallback.

## Required density and Lagrangian terms for analytical MRSF Hessian

A scientifically defensible analytical MRSF Hessian needs at least these categories, each compared against a numerical oracle before being enabled:

1. Ground-state/reference derivative reuse:
   - nuclear repulsion second derivatives;
   - one-electron kinetic, overlap, and nuclear-attraction second derivatives;
   - reference Coulomb/exchange/XC second-derivative contributions;
   - ROHF reference orbital-response consistency.
2. MRSF state-density and relaxed-density plumbing:
   - target-state eigenvector density for the chosen physical MRSF root;
   - relaxed density from the MRSF Z-vector;
   - `OQP_td_p`, `OQP_td_abxc`, and `OQP_td_mrsf_density` derivative contracts;
   - consistent AO/MO transformations and spin blocks for `mrst=1/3/5`.
3. Lagrangian stationarity terms:
   - reference SCF orbital-response terms;
   - TD/MRSF amplitude-response terms;
   - Z-vector derivative terms;
   - constraints for ROHF semicanonical/open-shell occupations.
4. Response and coupling terms:
   - derivative of the MRSF response matrix/operator;
   - derivative of spin-flip coupling blocks;
   - derivative of DFT XC response kernels for MRSF densities;
   - derivative of special spin-coupling (`spc`) density terms.
5. Symmetry and invariance checks:
   - Hessian symmetry;
   - rotational/translational mode behavior;
   - consistency between serial/OpenMP/Linux-MPI runs;
   - analytical-vs-numerical agreement over multiple FD steps.

## Root-selection and root-tracking hazards

MRSF root numbering is not the same as ordinary ground-state indexing:

- OpenQP root 0 is the high-spin reference in the Python guard language.
- OpenQP root 1 maps to physical S0 for MRSF Hessian-target diagnostics.
- OpenQP root 2 maps to physical S1, and so on.

Hazards:

1. The generic numerical Hessian path preserves only the requested state index. It does not verify that plus/minus displaced geometries keep the same electronic character.
2. The Z-vector code currently clamps `target_state` to `nstate` internally. A Hessian oracle must fail early instead of letting an out-of-range request silently become the highest computed state.
3. Near-degenerate MRSF roots can swap order under small displacements.
4. TD eigenvectors can change sign or rotate in a near-degenerate subspace. Sign-only checks are insufficient for multidimensional near-degenerate blocks.
5. `<S^2>`, excitation character, transition-density signatures, and vector overlaps should be recorded as diagnostics for each displaced point.

## Gate 3B numerical oracle plan

Default implementation choice: central finite difference of existing MRSF gradients.

Fallback implementation choice: central finite difference of MRSF total energies only if a gradient calculation is unavailable or fails capability checks. Energy FD is lower priority because it is more expensive and less directly comparable to the existing gradient implementation.

Required Gate 3B behavior:

1. Preserve the requested OpenQP state index and physical root label in metadata.
2. Refuse `hess.state=0` for TDHF/MRSF Hessians.
3. Refuse `hess.state > tdhf.nstate` before entering Fortran.
4. For every displaced geometry, collect:
   - total state energies for all computed roots;
   - selected-root energy;
   - target gradient when gradient FD is used;
   - available TD/MRSF eigenvector or transition-density fingerprints;
   - available `<S^2>` diagnostics;
   - convergence status for SCF, Davidson, and Z-vector.
5. Compare plus/minus displaced root character using, in priority order:
   - eigenvector overlap/fingerprint overlap when exported;
   - transition-density overlap when available;
   - energy continuity and minimum root gap as weaker diagnostics;
   - `<S^2>` continuity as an auxiliary diagnostic.
6. Detect root flips or ambiguous near-degenerate assignments and fail loudly with a diagnostic table rather than producing a Hessian.
7. Run FD step sweep at `1.0e-3`, `3.0e-4`, and `1.0e-4` Bohr.
8. Select the default step only after checking Hessian symmetry, repeatability, step convergence, and absence of root flips.

## Analytical implementation roadmap

Small commits by gate:

1. Gate 3A: design document plus guards/tests that reject unsupported MRSF Hessian modes.
2. Gate 3B: root-tracked numerical MRSF Hessian oracle from central FD of gradients; energy-FD fallback only behind explicit capability detection.
3. Gate 3C: input/API routing for supported numerical MRSF Hessian and explicit analytical rejection.
4. Gate 3D.1: reuse validated HF/DFT nuclear and one-electron second-derivative infrastructure behind comparison to numerical MRSF oracle.
5. Gate 3D.2: add MRSF target-state density and Lagrangian plumbing, tested against oracle deltas.
6. Gate 3D.3: add response/Z-vector derivative terms.
7. Gate 3D.4: add spin-flip/MRSF coupling terms and DFT XC response terms.
8. Gate 3D.5: enable full analytical MRSF Hessian only after full analytical-vs-numerical comparison is GREEN over the benchmark set.
9. Gate 3E: benchmark-quality report over stable examples, including serial/OpenMP and practical Linux-MPI representative evidence.

## Benchmark seed set

Use existing stable `examples/other/*mrsf*` water examples first because they are small and already part of the OpenQP example ecosystem:

- Minimal closed-shell-like physical root through ROHF triplet reference/MRSF singlet target: `h2o_rohf_mrsf-s_6-31g_hf.inp` or a smaller derivative of it.
- Spin-flip/MRSF state: `h2o_rohf_mrsf-t_6-31g_hf.inp` or a DFT analogue after HF is stable.
- Root-tracking stress case: a compact displaced/modified water MRSF case with adjacent roots close enough to exercise overlap diagnostics, but not so close that the oracle is unusable.

Benchmark report fields:

- state energy;
- gradient norm when available;
- Hessian symmetry max error;
- FD step convergence table;
- lowest modes and imaginary-mode count;
- root-tracking diagnostics and root-flip status;
- serial vs OpenMP reproducibility;
- Linux MPI representative result when practical.

## Gate 3A regression expectations

Gate 3A adds/keeps tests proving:

1. Analytical MRSF Hessian is rejected explicitly and includes root mapping and available `<S^2>` diagnostics.
2. Generic numerical MRSF Hessian is rejected until the Gate 3B root-tracked finite-difference oracle exists.
3. Input validation rejects unsupported MRSF numerical and analytical modes before runtime.
4. SF analytical rejection remains separate from MRSF rejection.
5. Runtime Hessian dispatch refuses MRSF numerical Hessian rather than entering generic state-index-only FD.
