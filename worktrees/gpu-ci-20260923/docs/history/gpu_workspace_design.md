# Stage-1 Unified GPU Workspace Manager — Design Note

**Module:** `pyoqp/oqp/utils/gpu_workspace.py`
**Status:** Stage-1 (planning + bookkeeping only). No CUDA, Fortran, or CMake.
**Supersedes (for direct integration):**
`perf/gpu-metc-persistent-buffers` (PMB) and `perf/tdhf-xc-response-cache` (XCC).

## Purpose

`GpuWorkspaceManager` is the single source of truth for two questions that the
two legacy planning experiments answered separately:

1. **Residency** — *where* a planned GPU buffer lives.
2. **Reuse** — *when* a previously planned allocation can be reused.

It unifies the legacy schemas without merging the off-base PMB/XCC history. The
legacy planning modules are retained verbatim as the source-of-truth shape
contracts and are mapped into the unified manager by parity adapters:

| Legacy module | Legacy type | Unified entry point |
|---|---|---|
| `gpu_metc_buffers.py` | `PersistentMetcBufferPlan` | `from_metc_plan(plan, *, f3_policy, nthreads)` |
| `gpu_metc_persistent_runtime.py` | `PersistentMetcAllocationRegistry` | `GpuWorkspaceManager` (register/lookup/release) |
| `tdhf_xc_response_cache.py` | `XcResponseCachePlan` | `from_xc_response_plan(plan, *, dtype_bytes)` |

## Residency classes

`Residency` (enum):

- `HOST_ONLY` — staged on the host; never uploaded.
- `DEVICE_RESIDENT` — lives on the device for the plan's lifetime.
- `MIRRORED` — exists on both host and device (upload/download).
- `BORROWED` — aliases memory owned by another allocation; the borrower must
  not free it.

Stage-1 keeps the *kernel* residency refit out of scope. The later refit
(`d3 → RESIDENT_INPUT`, `f3 → RESIDENT_ACCUM`) is intentionally **not** done
here; this note only establishes the residency vocabulary and bookkeeping.

## C-ABI-ready verbs

The manager exposes four verbs designed to map onto a future C ABI:

- `validate_table(target, table)` — normalize/validate a workspace table
  (rows are `WorkspaceBuffer` or `(name, bytes, role, residency)`); rejects
  type drift, bool-as-int byte counts, duplicate names, and non-`Residency`
  residency values.
- `allocate(target, reuse_key, table, *, f3_policy, nthreads, scalar_layout)` —
  validate and record a reusable allocation; re-allocating the same key with a
  divergent table raises a reuse conflict.
- `borrow(target, reuse_key)` — return a `BORROWED` alias of an owned record
  without transferring ownership.
- `release(target, reuse_key)` — release one outstanding borrow, or (once no
  borrows remain) free the owning record.

## Target namespacing

Every record is keyed by `WorkspaceKey(target, reuse_key)` where `target` is one
of `VALID_TARGETS = ("metc", "xc_response")`. Because the namespace prefixes the
subsystem reuse key, a METC `density` buffer and an XC-response `density` buffer
**cannot collide** even though they share a role name. This mirrors handoff
decision 4: the config's `target` may be multi-valued, but residency/cache
records are always stored under one concrete target string at a time.

### Reuse-key parity

- **METC** legacy reuse key `(nbf, nf, nmatrix, max_integrals, dtype_bytes)`
  maps to `WorkspaceKey("metc", (nbf, nf, nmatrix, max_integrals, dtype_bytes))`.
  The legacy allocation manifest (`ids`, `integrals`, `density`, `fock`) and its
  bytes/roles survive the mapping unchanged (with `nthreads=1`).
- **XC-response** legacy reuse key
  `(functional, basis, scf_type, response_type, nbf, ngrid, spin_channels)`
  maps to `WorkspaceKey("xc_response", …)`. The legacy scalar
  `(name, offset, length)` layout is preserved verbatim in
  `WorkspaceAllocation.scalar_layout`, so cache offsets are not lost; byte sizes
  apply the dtype width on top of the scalar lengths.

## f3 accumulator policy

`F3AccumulatorPolicy`:

- `PER_THREAD` (**default**) — one Fock/`f3` accumulator per OpenMP thread,
  mirroring the host `this%f3` residency model. `from_metc_plan` replicates the
  `output_matrix` (Fock) buffer bytes by `nthreads`.
- `SINGLE_ATOMIC` — one shared accumulator updated with atomics. Supported, but
  selected only when requested explicitly.

The manager never silently collapses the accumulator to a single global buffer:
the policy is an explicit field on every METC allocation, and `PER_THREAD` is
the default so the first refit preserves the CPU threading model.

## Determinism note (forward-looking)

When the kernel refit lands, `SINGLE_ATOMIC` implies `atomicAdd`, whose
floating-point summation order is non-deterministic. Any CPU/GPU parity gate
built on top of this manager must therefore stay tolerance-based
(`rtol≈1e-9, atol≈1e-11`), never bit-exact. This is recorded here so the policy
choice and the parity-tolerance requirement travel together.

## Explicitly out of Stage-1 scope

No CUDA kernels, no Fortran changes, no CMake changes, no METC runtime refit, no
MRT work, no DVT timer work, no ERI / HF-DFT J/K or quadrature work. Those are
later, separately-gated steps.
