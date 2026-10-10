# GPU branch integration — 2026-09-22

Base: `open-quantum-platform/openqp-gpu:integrate-native-pieces`,
`ad8783e469e50aa2072ffd286adefd4fe7cbcb18`.

The audit found 21 distinct GPU-related or GPU-supporting tips across personal/private OpenQP
repositories and local worktrees. `gpu-branches.json` lists every source alias,
SHA, and disposition. `source-files.json` maps imported files to immutable
source commits and records original and integrated SHA-256 values.
No source branch was deleted, archived, or rewritten.

## Implemented in this repository

- `src/metc/`: direct-integral MRSF/UMRSF contraction backend, exported as
  `libopenqp_gpu_metc`. It is independent of the existing density-fitting
  SCF/sigma engine. The later `beb51a2617ad` implementation contains reference,
  combined-Coulomb, warp-aggregation, and two-phase accumulation kernels, plus
  profiling and persistent sessions. Earlier METC development stages are
  represented by this implementation; they are not merged as full OpenQP trees.
- `include/openqp_gpu/metc.h`: public C ABI, tensor layout and ownership contract.
- `python/openqp_gpu/`: workspace planning, persistent allocation records,
  GPU configuration and XC cache planning. These planners do not allocate CUDA
  memory. `workspace_legacy.py` preserves the distinct local manager API from
  `6413215eeb41`; `gpu_workspace.py` is the later canonical planner.
- `tests/test_gpu_metc_regression.py` and `benchmarks/metc/`: recovered local
  numerical regression and kernel benchmark from `fc92929bbdbb`, adapted to
  the standalone library.

Integration fixes preserve nonzero input accumulators in the owning wrapper,
propagate CUDA copy/zero/free errors, reject invalid indices and oversized
32-bit kernel shapes, avoid int64 range-check overflow, reset diagnostic
accumulators between passes, and return an error for a host-only no-op.
The package's existing DF preparation script now runs only when explicitly
called, so importing the workspace planners no longer starts a PySCF calculation.

## Preserved integration and experimental code

- `seam/metc/`: original Fortran session/error-state bridge and host benchmark
  tools. These depend on an OpenQP host build; they are not automatically compiled
  into the standalone C/CUDA library. Tools require explicit host source/build/
  input paths. Their historical auto-build defaults are not the current package
  build policy. Do not use a tool's skipped result as validation.
- `seam/routec/`: committed Route-C host bridge plus six source-only patches for
  the external SCF J/K, XC, gradient and low-rank response connections. Apply only
  after comparing against the exact internal OpenQP target. `d9a752d4002a` adds
  README material beyond `4d8afd6e122f`; it contains no newer GPU engine.
- `experimental/xc_response/`: the latest packed elementwise response scaffold.
  This is not full XC quadrature/response and is excluded from production targets.
- `seam/mrsf/` and existing `src/sigma.cu`: retained. `gpu-mrsf-seam` was already
  integrated upstream through PR #265; importing it again would duplicate work.
- `docs/history/`: original design/validation notes, not new benchmark evidence.
  ERI kernels remain maintained independently in `libintRot`/`Spherical_ERI`.

The rotation-seam CUDA headers were byte-identical to the seven headers already
in `tools/dfbuild/sph_eri`. Their CPU host changes remain in the OpenQP/ERI
repositories. The anonymous `claude/elegant-hypatia-wZpoB` branch also supplied
a recovered Fortran/CUDA term-equivalence test, including a negative mutation
check. GPU benchmark timer parsing and its real-log fixture were recovered
from `fix/oqp-timer-parser-blockers`.

The local workspace and regression worktrees were clean at inspection. The
Route-C checkout had unrelated uncommitted Z-vector changes and untracked
files; those were not included. This is repository/backend consolidation, not
an automatic switch of the current OpenQP host to METC and not a new method.
The single-device METC interface must be called serially. Its per-thread slices
are separate accumulators, not permission for concurrent global profiler calls.
Existing multi-GPU DF functionality is unchanged.

## Build and validation

```sh
cmake -S . -B build -DOPENQP_GPU_METC_ONLY=ON -DCMAKE_CUDA_ARCHITECTURES=80
cmake --build build -j2
OPENQP_GPU_METC_LIB="$PWD/build/libopenqp_gpu_metc.so" \
OPENQP_GPU_METC_REQUIRE=1 \
python3 -m unittest discover -s tests -p test_gpu_metc_regression.py -v
```

Normal builds also enable `OPENQP_GPU_METC=ON`; use OFF to omit it. The METC-only
build does not require BLAS, generated rotation code, or a complete OpenQP build.
Use `OQP_GPU_METC_VARIANT=reference|combined_coulomb|combined_coulomb_warp|two_phase_accum`
to select the preserved UMRSF accumulation implementation.

Verified on CHC A100, CUDA 12.2.140, GCC 8.5.0, CMake 3.26.4:

- Slurm job **5450346**, `COMPLETED`, exit `0:0`, 22 seconds.
- Standalone shared library and recovered CUDA benchmark compiled successfully.
- **12 numerical tests passed**, including 16 mode/pass/variant subcases and
  resident two-slice repeated accumulation with pass reuse. No GPU skips.
- Resident diagnostic maximum absolute error was **4.4408920985006262e-16**;
  the Python comparison required `rtol=1e-9`, `atol=1e-11`.
- **79 host tests passed**: workspace/planning/runtime/import, including overflow
  and host-only contraction rejection, five source-equivalence checks and eight
  timer-parser tests. GitLab CI reruns these host checks.

The source-pinned manifest and raw logs are in `validation/`. Slurm accounting
service (`sacct`) was unavailable; the scheduler's `scontrol` terminal receipt
independently confirms success. Lease bundle
`LEASE-27E140AFE70647129589DF5A38031D04` was released after process exit;
manifest pair and submit resource each had zero active leases after completion.

These checks establish kernel equivalence to the preserved CPU contractions,
not full molecular/host integration, speedup, or validity of a new MRSF/UMRSF
model. Existing DF/gradient paths were not changed or re-benchmarked.
