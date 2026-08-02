# BLAS/LAPACK threading inside OpenQP

OpenQP runs its own OpenMP regions *and* links a threaded BLAS. The two have to
be told about each other, or they either oversubscribe the machine or cancel
each other out. This note records what OpenQP does about that, why, and the one
knob you can turn.

| Knob | Default | Effect |
| --- | --- | --- |
| `OQP_EIGEN_ROWS_PER_THREAD=<n>` | `512` | rows of the matrix required before a dense eigensolve asks for one more BLAS thread; `0` (or negative) disables the policy entirely |

Nothing else here is configurable — the rest is unconditional and needs no
attention. Read on only if you are tuning, or if a number below disagrees with
what you measure on your own machine.

## Rule 1 — a parallel region declares its own width

Code that calls no BLAS but runs its own OpenMP region pins BLAS to one thread
around it, so the BLAS worker pool does not spin against the region's threads.
`int2_twoei` (the two-electron Fock build) is the important case.

That pin has a trap. `blas_thread_set` resolves at run time through
`dlsym(RTLD_DEFAULT, "openblas_set_num_threads")`, and **in the OpenMP build of
OpenBLAS that symbol is `omp_set_num_threads`**. Pinning "BLAS" therefore also
pinned OpenMP, and the region the pin existed to protect ran serial.

So every such region now **declares** its width with `num_threads`, taken from a
count read *before* the pin (after it, `omp_get_max_threads` already reads 1). A
declared width cannot be undone by a global setter, whichever library owns it.

This is invisible on a serial or pthreads OpenBLAS, where the two setters really
are different functions — which is why it survived a long time. It costs
everything on a build that links `libopenblaso64`:

Plain RHF, CH2O/cc-pVTZ, wall seconds:

| threads | 1 | 8 | 32 | 64 | 128 |
|---|---:|---:|---:|---:|---:|
| before | 8.02 | 8.00 | 8.13 | 8.33 | 8.78 |
| after | 7.98 | 1.54 | **0.98** | 1.16 | 1.44 |

Flat to within noise before — a 152-core machine was not beating one core.

## Rule 2 — dense eigensolves are sized to the problem

LAPACK's symmetric eigensolvers do not thread the way GEMM does. Inside OpenBLAS
they are largely the netlib reference implementation, so they parallelise only
through the BLAS calls they themselves make — and the tridiagonal reduction that
dominates `DSYEVD` is memory-bound level-2 work spread over many small calls.
Each of those forks and joins a full thread team.

`DSYEVD` wall time in ms, best of several repeats, on a 2 × 38-core Xeon 8368
against `libopenblaso64` (OpenBLAS 0.3.15, `Core: SkylakeX`):

| n | 1t | 2t | 4t | 8t | 16t | 32t | 64t | 128t |
|-----:|------:|------:|------:|------:|------:|------:|------:|-------:|
|  128 | **0.64** | 1.55 | 4.49 | 6.97 | 18.0 | 23.2 | 54.9 | 190.7 |
|  256 | **3.11** | 5.10 | 8.00 | 12.4 | 37.7 | 49.4 | 105.8 | 509.6 |
|  512 | **19.8** | 22.6 | 31.9 | 53.3 | 254.5 | 334.6 | 977.1 | 3630.8 |
| 1024 | 131.4 | **122.7** | 152.1 | 319.0 | 1555 | 2092 | 5547 | 30527 |
| 2048 | 935.9 | 741.8 | **583.3** | 1203 | 5497 | 7353 | 18552 | 102936 |
| 4096 | 9279 | 5887 | **4302** | 4935 | 16308 | 21139 | 48123 | 259416 |

One thread is strictly best below n = 1024, nothing above four threads wins at
any size tested, and at n = 128 the 128-thread case is **299× slower** than the
1-thread one.

**This is the LAPACK layer, not the BLAS.** On the same library, in the same
sitting, `DGEMM` scales properly — 90.6 GFlop/s at one thread to 1775 at 128
(n = 2048, 19.6×). GEMM-heavy paths are therefore left alone; only eigensolves
are gated.

What ships is one linear rule rather than that table: **one BLAS thread per
`OQP_EIGEN_ROWS_PER_THREAD` rows**, applied in `source/mathlib/eigen.F90` and
around the remaining direct `DSYEV`/`ZHEEV` call sites (`trah_converger.F90`,
`soc_mrsf.F90`). The table above is one machine and one OpenBLAS: the mechanism
generalises, the constants do not. A monotone rule with no cliff degrades
gracefully on a BLAS that *does* thread these routines well — it leaves some
speed unclaimed on very large matrices, where the eigensolve is already small
next to the O(nbf⁴) integral work around it.

Two safety properties, both worth knowing if you are reading the code:

* the policy only ever **lowers** the count, so a caller that has already pinned
  BLAS keeps its pin; and
* it changes nothing inside an OpenMP parallel region — `blas_thread_set`
  mutates process-global state, so calling it at all from several workers would
  race even when they all write the same value.

### What it buys

The methods that gain are the ones doing many small eigensolves. H2O/cc-pVTZ
CAS(6e,6o) CASSCF with the analytic Hessian, wall seconds:

| threads | 1 | 8 | 32 | 64 | 128 |
|---|---:|---:|---:|---:|---:|
| without | 37.60 | 10.42 | 14.39 | 34.73 | 125.69 |
| with | 37.76 | 9.60 | **6.57** | 7.99 | **11.77** |

The shape changes, not just the constant. Without the policy CASSCF *degrades*
past 8 threads — more hardware made it 12× slower. With it the optimum moves to
32 threads and the curve is nearly flat.

Plain SCF gains much less (see the RHF table above, 0.98 → 0.81 s at 32
threads): one diagonalisation per iteration is not where its time goes.

### Reproducibility

Changing a thread count reorders floating-point reductions inside the BLAS, so
bit-identity is not guaranteed — and measurably does not always hold. Toggling
`OQP_EIGEN_ROWS_PER_THREAD=512` versus `=0` in the same binary:

* all 16 `examples/WF_methods` cases (CASCI, CASSCF, SA-CASSCF, CASPT2, MS- and
  XMS-CASPT2, NEVPT2, SC-NEVPT2, MRMP2, MCQDPT2, XMCQDPT2) — identical energies
  to every printed digit;
* the `examples/TRAH` and `examples/SOC` cases — every non-zero quantity
  identical to every printed digit, with **the sign of some printed zeros
  differing** (`0.00000000` vs `-0.00000000`). That is a reduction-order
  artefact of exactly the kind described above, and carries no information.

If you need a hard bitwise guarantee, do not rely on this being benign in
general — fix the thread count instead, as `perf=0` does.


## Checking your own build

The behaviour above depends on *which* OpenBLAS you linked. To see:

```bash
ldd $(python -c "import oqp,os;print(os.path.dirname(oqp.__file__))")/../liboqp.so | grep -iE "openblas|gomp|omp"
```

* `libopenblaso64` **+** `libgomp` — the OpenMP build, sharing one runtime with
  liboqp. This is the configuration the notes above describe.
* `libopenblas64` alone — the **serial** build. Many distributions install this
  as the default `FindBLAS` target, and it silently caps every BLAS call at one
  thread. Point the link at the OpenMP one with `CMAKE_LIBRARY_PATH`.
* `libopenblasp64` — the pthreads build. Distinct setters, so Rule 1's pin does
  real work here; Rule 2 still applies.
* both `libgomp` **and** `libomp` — two OpenMP runtimes in one process. Fix this
  before trusting any timing.

One more environment trap worth checking on older OpenBLAS: its CPU probe may
not recognise a recent processor and can fall back to a decades-old microkernel.

```bash
OPENBLAS_VERBOSE=2 python -c "import oqp" 2>&1 | grep Core:
```

If that prints something implausible for your hardware (`Core: Prescott` on a
Xeon, say), set `OPENBLAS_CORETYPE` explicitly — it is read at load time and
cannot be changed afterwards. On the machine benchmarked here that alone was
worth 4.2–4.9× on DGEMM.
