# CI coverage and runner requirements

CI runs for merge requests, standalone branch pushes, and tags. An open MR
replaces the branch pipeline. There are deliberately no macOS jobs.

| Job | Coverage | Runner |
| --- | --- | --- |
| `linux-host` | Python planning/import/source-equivalence tests and the compiled C workspace runtime | `small-linux` |
| `windows-python` | Portable Python planning/import/source-equivalence tests | `windows`, `winserver1` |
| `linux-cuda-build` | All default CUDA targets, install step, shared-library dependencies, and a separate METC-only build | `small-linux` Docker |
| `linux-cuda-regression` | Real-device METC numerical regression using the exact libraries from the build job | `linux`, `nvidia-gpu` Docker |

Every job is required. Unit tests publish JUnit reports and fail if any test is
skipped or a selected module contains no tests. The CUDA suite additionally sets
`OPENQP_GPU_METC_REQUIRE=1`; a missing CUDA device/library is an error. CUDA
toolchain, build, install, linkage, hashes and device information are retained
for 14 days. Compilation on a CPU runner does not establish GPU correctness.

## GPU runner deployment status

At setup inspection on 2026-09-23, Linux CPU and Windows runners were online,
but no NVIDIA GPU runner was registered. The required GPU job intentionally
remains pending until one is assigned. Do not make this job optional, disable
it with a variable, or label a CPU runner as a GPU runner to turn CI green.

The GPU runner must use Linux x86_64, Docker with NVIDIA Container Toolkit,
an NVIDIA driver compatible with CUDA 12.6, and an A100 (SM80) device for this
initial validated architecture. Assign it to this project with both `linux`
and `nvidia-gpu` tags, disable untagged jobs, and expose an allocated GPU to the
container. Use one concurrent GPU job; the pipeline also serializes this
project's numerical jobs. A CHC Slurm GPU must be allocated through Slurm and
the central admission procedure, never used by an unrestricted login-shell
runner. Provisioning that runner is separate from adding these pipeline jobs.

Enable **Pipelines must succeed** before merging. A pending GPU job then blocks
merge rather than accepting a CPU-only success as full validation.

## Scope

The Windows job does not claim native Windows CUDA support. The current build
contains Linux-specific OpenMP/linker options; enabling that platform needs a
separate native build and GPU validation.

The standalone DF library currently uses 32-bit C integers at its host
BLAS/LAPACK boundary. CI explicitly selects matching OpenBLAS LP64. This does
not change OpenQP's ILP64 engine policy and does not prove the two libraries
can safely share BLAS symbols in one process. Host ABI and full molecular
integration remain separate tests in the OpenQP repository.

The GPU numerical suite currently covers METC contractions, accumulation,
layout and persistent workspace behavior. Full-library compilation covers
SCF, XC, response, gradients and DF building, but does not validate their
complete molecular results. Historical benchmark tables are not CI evidence.

To run host checks locally in a disposable Python 3.11 environment:

```sh
python3.11 -m venv .venv-ci
.venv-ci/bin/python -m pip install -r ci/requirements.txt
.venv-ci/bin/python ci/run_tests.py host
```
