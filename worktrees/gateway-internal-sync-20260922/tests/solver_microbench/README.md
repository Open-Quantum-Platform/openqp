# Z-vector solver micro-benchmark

A self-contained harness that exercises the real `source/pcg.F90` and
`source/minres.F90` on synthetic symmetric systems. It needs **only `gfortran`**
— no BLAS/LAPACK, Libint, or the rest of OpenQP — so it runs in environments
where the full package cannot be built, and isolates the solver algorithms from
the (expensive) Fock-build operator they normally drive.

## Run

```bash
cd tests/solver_microbench
./run.sh
```

Override the compiler/flags or the baseline commit if needed:

```bash
FC=gfortran FFLAGS="-O3 -ffree-line-length-none" OLD_REF=f601f83 ./run.sh
```

## What it checks

1. **Performance (PCG).** Builds the pre-optimization `pcg.F90` (`OLD_REF`) and the
   current one, solves the same SPD system many times, and reports wall time.
   The optimization (carried `rz` numerator + scalar-reduction breakdown guards
   instead of per-iteration full-vector scans) must give the **same iteration
   count and the same residual** while running faster.

2. **Stability (PCG).** Injects a NaN into the operator output and a +Inf into
   the preconditioner; both must produce `PCG_BREAKDOWN` (errcode 3) with no
   non-finite value leaking into the returned solution, and an exact initial
   guess must converge at `init`.

3. **MINRES correctness/robustness.**
   - SPD system: MINRES must match PCG (same solution).
   - Indefinite system (`Laplacian − 2I`): MINRES must still converge.
   - `A = diag(+1, −1)`, `b = (1,1)`: the textbook case where `pᵀAp = 0` on
     step 1 — **CG breaks down, MINRES solves it in 2 steps**. This is the
     justification for offering MINRES as the symmetric-indefinite fallback.

## Representative results (gfortran 13.3, -O2)

```
PERF  OLD pcg.F90 : n=200000, 18 iters/solve, total 19.9 s
PERF  NEW pcg.F90 : n=200000, 18 iters/solve, total 12.3 s   (~1.6x, identical residual)
STAB  nan_in_Ax / inf_in_Minv : PCG_BREAKDOWN, no non-finite solution
MINRES T1 SPD        : matches PCG, max|Δx| = 4.7e-10
MINRES T2 indefinite : ||b-Ax||/||b|| = 1.9e-12
MINRES T3 diag(1,-1) : MINRES solved (2 iters), CG broke down (0 iters)
```

> The ~1.6× PERF speedup is the *cheap-matvec* upper bound: here the matvec is a
> 3-point stencil costing about as much as the vector ops, so removing redundant
> vector passes is visible. In a real z-vector solve each iteration is a full
> Fock/ERI build that dwarfs the vector work, so the wall-clock effect of this
> particular change is small there — its value is the proven equivalence and the
> preserved fail-closed stability, not raw speed.
