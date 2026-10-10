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
| `LiH_SA-CASSCF_grad.inp` | `casscf` + state average | — | compatibility central-difference gradient of a state on common SA-CASSCF orbitals |
| `LiH_SA-CASSCF_ANALYTIC_grad.inp` | `sa-casscf` | — | analytic gradient of the weighted SA-CASSCF objective |
| `LiH_SA-CASSCF_ROOT1_grad.inp` | `sa-casscf` | — | analytic individual-state gradient with the coupled orbital/CI Z-vector |
| `H2O_CASSCF_CAS44_AH.inp` | `casscf` | `[casscf] converger=ah, hessian=analytic` | trust-region AH converger + exact orbital Hessian (17 vs 335 CI solves, both at 15 macroiterations) |
| `H2O_CASSCF_CAS44_TRAH.inp` | `casscf` | `[casscf] converger=trah` | matrix-free trust-region converger: the Hessian is never assembled (303 vs 1963 CI evaluations against `converger=ah` on the default FD Hessian, 13 vs 19 macroiterations, same energy to 1e-13) |
| `H4_CASPT2.inp` | `caspt2` | `h0=fock` | single-state CASPT2 (Fock H0) |
| `H4_NEVPT2.inp` | `caspt2` | `h0=dyall` | NEVPT2 (Dyall H0), uncontracted |
| `H4_SC-NEVPT2.inp` | `caspt2` | `h0=dyall, contraction=strong` | strongly-contracted NEVPT2 (vs PySCF, nEh) |
| `H4_SC-NEVPT2_grad.inp` | `caspt2` | `h0=dyall, contraction=strong, gradient=analytic` | **analytic** SC-NEVPT2 nuclear gradient (vs 5-point central differences, 5.9e-10 Eh/Bohr) |
| `LiH_SC-NEVPT2_optimize.inp` | `caspt2` | `h0=dyall, contraction=strong, gradient=analytic` | geometry optimization driven by the analytic SC-NEVPT2 gradient; also covers a PT2 frozen core and exactly degenerate (2p-pi) virtuals |
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
- `[casscf]` — orbital optimization and convergence. For the dedicated
  `method=sa-casscf` analytic derivative, `gradient_state=averaged` selects the
  weighted objective and an integer selects an averaged root; `zvector_tol`
  and `zvector_degeneracy_tol` govern the individual-state response solve. The
  legacy `method=casscf` plus `[state_average] enabled=true` spelling retains
  its central-difference controls (`grad_step`, `grad_guess`, `grad_gap_warn`,
  and `grad_ranks_per_group`). State-specific CASSCF also uses an analytic
  derivative.
- `[state_average]` — `enabled`, `target_roots`, weights.
- `[input] method` — `caspt2`, `ms-caspt2`, `xms-caspt2`, `nevpt2`,
  `sc-nevpt2`, `mrmp2`, `mcqdpt2`, `xmcqdpt2`.  NEVPT2 used to be the one
  member without its own name (`caspt2` plus `[pt2] h0=dyall`); that spelling
  is the same calculation and still works, but two different theories reaching
  the gradient dispatch under one name is what made routing them a matter of
  inspecting options instead of reading the method.
- `[pt2]` — `h0` (fock=CASPT2 / dyall=NEVPT2), `contraction` (none / strong=SC-NEVPT2);
  a method that names them (`nevpt2`, `sc-nevpt2`) rejects an option that
  contradicts it rather than letting one win,
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

Two independent analytic derivatives sit behind that one selector, and the
route picks between them from `[pt2] h0` and `contraction`. Both cost one PT2
evaluation plus one pass of derivative integrals, independent of the number of
nuclei; the central-difference alternative costs `2 * 3 * natom` displaced PT2
energies.

### CASPT2 family (`h0=fock`)

Scope of the analytic route (everything else falls back under `auto`, and is
refused with a specific reason under `analytic`):

| supported | not supported |
| --- | --- |
| `caspt2`, `mrmp2` (single state) | `ms-caspt2` — the **multi-set** construction (per-state orbitals, per-state full-Fock-matrix H0, inter-state Lowdin-minor rotations); use `xms-caspt2` |
| `mcqdpt2` (single-set multistate) | `ipea_shift != 0` — the active-diagonal bias is not invariant under rotations inside the active block |
| `xms-caspt2`, `xmcqdpt2` | `[cas] orbital_source` reading orbitals from a file — imported orbitals are not a differentiable function of the geometry |
| `level_shift`, `imaginary_shift`, `edshft` — carried exactly | |
| `reference=casci` (RHF orbitals) and `reference=casscf` (state-specific or state-averaged) | |
| the PT2 frozen core (`[pt2] frozen`) | |

`h0=dyall` with `contraction=strong` is NOT out of scope: it is the other
route, below. `h0=dyall` with `contraction=none` (uncontracted NEVPT2) has no
analytic derivative on either route and falls back.

The derivation also has preconditions that hold almost everywhere and can fail
at a particular geometry: canonical (`casci`) or stationary (`casscf`) reference
orbitals, semicanonical PT2 orbitals, a nonsingular orbital response,
non-degenerate effective-Hamiltonian roots, and a reconstruction that reproduces
the reported energy. These are preconditions of the **route**, not of the
energy, so `auto` treats a failed one exactly like an out-of-scope variant — it
falls back to central differences and writes the offending number into the log —
and `analytic` refuses. A penalty-function MECI search walks into the degenerate
case by construction; the fallback is what keeps it running.

### Strongly contracted NEVPT2 (`method=sc-nevpt2`)

Analytic when the run is strongly contracted, state specific, on a
state-specific CASSCF reference, with `runtype` in grad/optimize/ts/mep/irc. It
costs one CASSCF, one second-order pass, one coupled orbital/CI response solve
and one set of derivative integrals. Its error is FIRST order in the residual
CASSCF orbital-rotation gradient, so `[casscf] gradient_norm_tol` controls how
many digits it is worth; the run reports `|g_orb|`, the response residual and
the Lagrangian asymmetry so the size of every approximation is visible rather
than assumed.

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
