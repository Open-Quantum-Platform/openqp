# Wavefunction-method examples

A single, flat directory of curated native wavefunction-method inputs — one
small, option-complete, self-documenting example per method.  Every file carries
a header comment stating what it computes and the reference energy it reproduces.

Run the whole set with:

```bash
openqp --nompi --run_tests WF_methods
```

or a single one with `openqp --nompi <file>.inp` (results land in `<file>.log`).

## Examples (one folder, one file per method)

| File | `method=` | `[pt2]` selector | What it shows |
| --- | --- | --- | --- |
| `H2_FCI.inp` | `fci` | — | full configuration interaction |
| `H2_CASCI.inp` | `casci` | — | active-space CI at fixed RHF orbitals |
| `LiH_CASSCF.inp` | `casscf` | — | state-specific CASSCF (orbital + CI optimization) |
| `LiH_CASSCF_grad.inp` | `casscf` | — | analytic nuclear gradient of the state-specific CASSCF root |
| `LiH_CASSCF_optimize.inp` | `casscf` | — | geometry optimization with the analytic state-specific CASSCF gradient |
| `H2O_CASSCF_CAS44.inp` | `casscf` | — | CASSCF, CAS(4,4) (vs PySCF −75.0085688625) |
| `LiH_SA-CASSCF.inp` | `sa-casscf` | — | state-averaged CASSCF (`[state_average]`) |
| `LiH_SA-CASSCF_grad.inp` | `sa-casscf` | — | central-difference gradient of a state on common SA-CASSCF orbitals |
| `H2O_CASSCF_CAS44_AH.inp` | `casscf` | `[casscf] converger=ah, hessian=analytic` | trust-region AH converger + exact orbital Hessian (17 vs 335 CI solves, both at 15 macroiterations) |
| `H2O_CASSCF_CAS44_TRAH.inp` | `casscf` | `[casscf] converger=trah` | matrix-free trust-region converger: the Hessian is never assembled (303 vs 1963 CI evaluations against `converger=ah` on the default FD Hessian, 13 vs 19 macroiterations, same energy to 1e-13) |
| `H4_CASPT2.inp` | `caspt2` | `h0=fock` | single-state CASPT2 (Fock H0) |
| `H4_NEVPT2.inp` | `caspt2` | `h0=dyall` | NEVPT2 (Dyall H0), uncontracted |
| `H4_SC-NEVPT2.inp` | `caspt2` | `h0=dyall, contraction=strong` | strongly-contracted NEVPT2 (vs PySCF, nEh) |
| `H4_MS-CASPT2.inp` | `ms-caspt2` | `h0=fock` | multi-set MS-CASPT2 (state-specific Fock) |
| `H4_XMS-CASPT2.inp` | `xms-caspt2` | — | extended multistate CASPT2 |
| `H4_MRMP2.inp` | `mrmp2` | — | single-state MRMP2 (Hirao; == caspt2 h0=fock) |
| `H4_MCQDPT2.inp` | `mcqdpt2` | `edshft=0.02` | GAMESS-convention multistate QDPT (Nakano, single-set) + ISA |
| `H4_XMCQDPT2.inp` | `xmcqdpt2` | — | Granovsky's extended QDPT (XZERO; invariant) |
| `H4_CASPT2_grad.inp` | `caspt2` | `gradient=analytic` | **analytic** CASPT2 nuclear gradient (Lagrangian: amplitude, CI and orbital response) |
| `H4_CASPT2_numgrad.inp` | `caspt2` | `gradient=numerical` | the central-difference route the out-of-scope variants still use |
| `H4_XMS-CASPT2_grad.inp` | `xms-caspt2` | `gradient=analytic` | analytic gradient of an XMS-CASPT2 root (state-rotation response, `H_eff` eigenvalue) |

## Method controls

- `[cas]` — `active_electrons`, `active_orbitals`, `frozen_core`, `orbital_source`.
- `[ci]` — `nroot`, `solver` (dense/davidson), `target_spin`.
- `[casscf]` — orbital optimization and convergence; `grad_step`, `grad_guess`,
  `grad_gap_warn`, and `grad_ranks_per_group` control the SA-CASSCF
  central-difference nuclear gradient. State-specific CASSCF uses the analytic
  derivative.
- `[state_average]` — `enabled`, `target_roots`, weights.
- `[pt2]` — `h0` (fock=CASPT2 / dyall=NEVPT2), `contraction` (none / strong=SC-NEVPT2),
  `frozen` (auto=standard deep cores, matches OpenMolcas), `ipea_shift`,
  `level_shift`, `imaginary_shift`, `edshft` (GAMESS ISA `d -> d + edshft/d`,
  QDPT-style intruder handling; exclusive with the level shifts); `gradient`
  selects the nuclear-gradient route (see below).

The native scope is determinant-space active-space CI and determinant-space
CASPT2/NEVPT2, plus RDM-based strongly-contracted NEVPT2 (`contraction=strong`)
and the GAMESS-convention QDPT family (`mrmp2` / `mcqdpt2` / `xmcqdpt2` =
single-set diagonal-Fock MRPT2, Granovsky extended H0 via `xmcqdpt2`).

## PT2 nuclear gradients

`[pt2] gradient` selects the route:

| value | behaviour |
| --- | --- |
| `auto` (default) | analytic where the variant has one, central differences otherwise |
| `analytic` | refuse rather than fall back |
| `numerical` | always central differences |

The analytic derivative is one PT2 evaluation plus one pass of derivative
integrals, independent of the number of nuclei; the central-difference route
costs `2 * 3 * natom` displaced PT2 energies.

Scope of the analytic route (everything else falls back under `auto`, and is
refused with a specific reason under `analytic`):

| supported | not supported |
| --- | --- |
| `caspt2`, `mrmp2` (single state) | `ms-caspt2` — the **multi-set** construction (per-state orbitals, per-state full-Fock-matrix H0, inter-state Löwdin-minor rotations); use `xms-caspt2` |
| `mcqdpt2` (single-set multistate) | `h0=dyall` (NEVPT2) |
| `xms-caspt2`, `xmcqdpt2` | `contraction=strong` (SC-NEVPT2) |
| `level_shift`, `imaginary_shift`, `edshft` — carried exactly | `ipea_shift != 0` — the active-diagonal bias is not invariant under rotations inside the active block |
| `reference=casci` (RHF orbitals) and `reference=casscf` (state-specific or state-averaged) | `[cas] orbital_source` reading orbitals from a file — imported orbitals are not a differentiable function of the geometry |
| the PT2 frozen core (`[pt2] frozen`) | |

The gradient is refused, rather than approximated, when the reference itself
cannot support it: an unconverged CASSCF, orbitals that are not semicanonical,
a degenerate effective-Hamiltonian root, or a reconstruction that does not
reproduce the reported energy. Each writes the offending number into the log.

## QDPT engines

The QDPT family's H0 is diagonal in the determinant basis, so it runs on a
matrix-free **direct engine** by default (`[pt2] engine=auto`): the
singles/doubles interacting space is streamed from the reference determinants
with a diagonal resolvent — no dense Hamiltonian, no external-space
diagonalization, no `max_det` wall.  Two backends implement the stream, both
pinned against each other and against `engine=dense` in
`tests/test_qdpt2_direct.py`:

- **NumPy streaming** (`engine=direct` forces it; `auto` uses it only when the
  native kernel is unavailable): vectorized
  class generation + lexsort merge.  Measured on the isolated streaming stage
  (CAS(8,8)-in-24 orbitals, 4900 references, 114.5M streamed terms): **~6 s**.
- **liboqp OpenMP kernel** (what `auto` prefers when built;
  `source/modules/qdpt2_kernel.F90`; threads default to min(8, cores),
  `[pt2] nproc` overrides): numerically identical.  History on the same
  benchmark: v1 hash tables 73 s (DRAM-bound, negative scaling); v2
  sort/reduce 18 s serial and **7.4 s at 8 threads** — ahead of NumPy — once
  a clamp bug was fixed (explicit thread requests were being capped by the
  core's ambient BLAS-control OMP setting of 1, which had silently
  serialized every earlier "parallel" measurement).  The serial T-way merge
  still limits scaling to ~2.4x; a key-range-parallel merge is the follow-up.

## Determinant-CI Fortran engine

`source/modules/fci_hamiltonian.F90` moves the CASCI hot loops into
Fortran/OpenMP behind Python fallbacks (equivalence-pinned to 1e-12 in
`tests/test_fortran_fci.py`): dense Hamiltonian assembly, the Davidson
matrix-free application `Y = 0.5(H+H^T)X` and diagonal, and `oqp_dsyevd` — a
dense eigensolver through liboqp's own ILP64 LAPACK, slotted between numpy
and the Jacobi fallback in `_symmetric_eigh`.  Together with a one-line fix
that vectorized the eigen-residual verification (previously an O(n^2)-Python
loop per root that alone cost ~20 minutes at 4900 determinants), H2O/6-31G
CAS(6,6) SA-2 MCQDPT2 now runs end-to-end in ~6 s.

The dense-path spin diagnostics now classify only a growing low-root prefix
(the full-spectrum sweep at 4900 determinants used to dominate everything),
and the engine's OpenMP build scales (3.5x at 8 threads) now that explicit
thread requests are honored.  Measured milestones: H2O/6-31G CAS(6,6) SA-2
MCQDPT2 end-to-end 6.2 s; CAS(8,8) SA-2 completes in ~171 s (previously a
>20-minute timeout).  Remaining CAS(8,8) time is dominated by the
semicanonicalization loop's repeated full dense eigensolves — profiling that
stage (fewer passes, or `[ci] solver=davidson` with the Fortran matvec) is
the next lever.  Measured
on this machine: H6/STO-3G SA-3 MCQDPT2 21x faster than the dense engine at
machine-precision agreement (<=1e-14 Eh); C2H4/cc-pVDZ SA-3 XMCQDPT2 (a space
the dense engine cannot even allocate, ~1e14 determinants) completes in ~0.7 s
end-to-end.  `engine=dense` keeps the reference path (equivalence-tested in
`tests/test_qdpt2_direct.py`); `[pt2] max_terms` guards the streamed-term
count.  Current limits: <= 63 orbitals per spin after the PT2 frozen core, and
work scales as N_CAS-dets x (occ^2 x virt^2) — large active spaces await the
planned Fortran-core port.
