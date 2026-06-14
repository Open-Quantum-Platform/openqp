# ddX-enabled PCM validation (reviewer / maintainer guide)

The Polarizable Continuum Model (PCM) energy path is **optional** and depends on
the external ddX continuum-solvation library (`-DENABLE_DDX=ON`). The default CI
job (`.github/workflows/CI.yml::gcc-build`) builds with `ENABLE_DDX=OFF`, so it
does **not** exercise the real ddX SCF path: every PCM-on benchmark in
`tests/test_pcm_literature_benchmarks.py` skips cleanly when OpenQP is built
without ddX (the C adapter reports the missing backend and the run aborts).

Building ddX from source (ILP64) plus running the SCF benchmark suite is slow, so
it is intentionally kept out of the per-push CI. This document is the exact,
reproducible recipe to validate the ddX/PCM path locally before merging or when
reviewing changes under `source/solvent_*`, `source/dftlib/`, or
`tests/*pcm*`.

## 1. Build OpenQP with ddX

```bash
cmake -B _build_ddx -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DUSE_LIBINT=OFF \
  -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=ON -DENABLE_MPI=OFF \
  -DENABLE_DDX=ON \
  -DENABLE_OPENTRAH=OFF \
  -DBUILD_TESTING=ON \
  -DCMAKE_INSTALL_PREFIX=.
cmake --build _build_ddx --target install -j
```

Notes:

* `-DENABLE_DDX=ON` builds ddX v0.8.0 from source as ILP64 (matching OpenQP's
  64-bit-integer BLAS) unless a prebuilt `-DDDX_ROOT=...` is supplied.
* `-DENABLE_OPENTRAH=OFF` avoids the OpenTrustRegion external, whose build is
  broken on some toolchains; PCM does not need it (`run_otr` falls back to the
  native solver).
* `libddx.so` is installed next to `liboqp.so` (see `source/CMakeLists.txt`), so
  no extra runtime path juggling is needed beyond pointing at `./lib`.

## 2. ddX link / adapter smoke (fast, seconds)

```bash
ctest --test-dir _build_ddx --output-on-failure -R 'ddx'
```

This runs `oqp_ddx_link_smoke` and `oqp_ddx_adapter_smoke` (real ddX through the
C adapter) — a quick check that the library links and the adapter numerics are
sane, independent of the QM SCF coupling.

## 3. ddX-enabled scientific PCM gates (minutes)

```bash
cd pyoqp && pip install -r requirements.txt && pip install . && cd ..

export OPENQP_ROOT=$PWD
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH   # finds liboqp.so + libddx.so
export OMP_NUM_THREADS=1

python -m pytest tests/test_pcm_literature_benchmarks.py \
                 tests/test_pcm_canonical_runtime_path.py -q
```

What this validates:

* Runs the real OpenQP Fortran SCF+PCM path on the closed-shell HF/ROHF and
  DFT (PBE/BHHLYP) 10-molecule tier and compares the reaction-field energy
  `e_pcm` to standard PySCF ddCOSMO references in
  `tests/data/pcm_literature_benchmarks.json` (tolerance `3e-4` Ha).
* Only rows with `status: "verified"` are enforced pass/fail gates; everything
  else skips with an explanatory message (never a false pass).

Expected result on a ddX-enabled build: all 30 verified gates pass
(`HF` MAE ≈ 0.034, `PBE` ≈ 0.027, `BHHLYP` ≈ 0.028 kcal/mol vs PySCF ddCOSMO).
On a non-ddX build the PCM-on rows skip and only the source-level canonical
checks run.
