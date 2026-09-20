# Contributor and AI-reviewer guide

This file is read by Codex, automated reviewers, and human contributors. Every
merge request is expected to satisfy the rules below. A reviewer should call
out, per rule, whether the merge request satisfies it or explain what is
missing.

## Code Review Rules

### 1. BLAS/LAPACK must go through the OpenQP wrapper layer

Never call a raw BLAS/LAPACK routine (`dgemm`, `dgels`, `dgetrf`, `dsyev`,
`zheev`, ...) directly. Fortran modules must `use oqp_linalg`, which renames the
`oqp_<name>_i64` wrappers back to the standard names, and then call the standard
name (`call dgels(...)` resolves to `oqp_dgels_i64`).

- Raw calls are **always allowed inside the wrapper modules**
  `source/mathlib/lapack_wrap.F90` and `source/mathlib/blas_wrap.F90` (their whole
  job is to call the raw symbol). Elsewhere, a bare call is acceptable only if the
  routine is registered (see the enforced-minimum note below).
- OpenQP is **ILP64-only** (`-fdefault-integer-8`, `BLA_SIZEOF_INTEGER=8`), so a
  bare call binds correct-width 8-byte integers on Linux. The failure a stray raw
  call causes is on **macOS Accelerate**, whose classic Fortran symbols are LP64:
  every referenced name must be interposed onto its `$NEWLAPACK$ILP64` variant via
  the alias list, so a routine that is CALLED but not in `OQP_ACCELERATE_ILP64_SYMS`
  fails `cmake/check_accelerate_aliases.cmake`.
- If a wrapper is missing, add `oqp_<name>_i64` to the wrapper module, add the
  `<name> => oqp_<name>_i64` rename in `source/mathlib/oqp_linalg.F90`, and add
  the routine to `OQP_ACCELERATE_ILP64_SYMS` in `cmake/oqp_functions.cmake`.
- The **enforced minimum** is registration in `OQP_ACCELERATE_ILP64_SYMS`, not
  literal `use oqp_linalg`: a number of long-standing sites call bare BLAS
  directly and are correct on the ILP64-only build (they predate the wrapper
  convention). New code should still prefer `use oqp_linalg` for clarity, but the
  CI gate rejects any *unregistered* BLAS/LAPACK routine regardless of call style.

**Enforced by CI:** `tools/check_blas_wrapper.py` (the `PR policy` workflow)
fails the build if any BLAS/LAPACK routine referenced from `source/` — or from the
`tests/fortran/*.F90` harnesses linked into `liboqp` — (outside the wrapper
modules) is not registered in `OQP_ACCELERATE_ILP64_SYMS`. It is the portable
complement to the macOS-only, link-time `cmake/check_accelerate_aliases.cmake`, and
CI runs the trusted base-branch copy of the script so a PR cannot weaken its own gate.

### 2. New functionality ships with a test/example under `examples/`

Any new capability — a new `runtype`, `method`, `[section]`, keyword, or opt-in
feature flag — must add at least one **small, fast** example under `examples/`
that exercises it (single-point or minimal steps; no big calculations).

- Opt-in boolean flags are already enforced by `openqp --check_feature_coverage`
  (a CI gate): a new flag must be set true by an example, or classified in
  `EXEMPT_FLAGS` / `KNOWN_UNCOVERED` in `pyoqp/oqp/utils/regression.py` with a
  reason.
- Give the example a reference where reproducible (`openqp --validate_examples`
  checks committed `.json` references); a demo that needs an optional backend
  (e.g. OpenMM) may be run-only.

**Reviewer check:** if the diff adds a new runtype/method/section/keyword but no
`examples/**` file, flag it.

### 3. User-facing additions get a pythonic Python-API addition

New user-facing features should be usable from the compact `OpenQP` Python API in
`pyoqp/oqp/openqp.py`.

- New `[section]` keywords are auto-exposed through the schema-driven proxies
  (`job.settings.<section>(...)`, `job.<section>.<key>`), so no manual work is
  usually needed for plain keywords.
- A new workflow/runtype should get a `job.workflow.<name>` entry (and a
  convenience helper like `job.qmmm(...)` where it improves ergonomics), plus a
  unit test in `tests/test_openqp_api.py`.

**Reviewer check:** if the diff adds a new workflow runtype or a new `[section]`
with no corresponding `openqp.py` handling or `tests/test_openqp_api.py` test,
flag it.

### 4. New functionality is documented in openqp-docs

User-facing keywords, sections, and workflows must be documented in the manual
repo [open-quantum-platform/openqp-docs](https://qchemlab.knu.ac.kr/open-quantum-platform/openqp-docs):
a keyword-page entry under `docs/keywords/` and/or a workflow page under
`docs/workflows/`, wired into `mkdocs.yml` nav.

- openqp-docs is a **separate repository**, so this merge request's CI cannot
  see it directly. Link the companion openqp-docs merge request in this merge
  request's description.

**Reviewer check:** if the diff adds/changes user-facing keywords or workflows,
confirm the merge request description links an openqp-docs merge request; flag
it if missing.

### 5. openqp carries product code, not development material

This repository holds the engine, its tests and examples, and the scripts the
build and CI actually consume. Method notes, derivations, validation harnesses,
performance investigations and one-off diagnostics belong in
[open-quantum-platform/openqp-devkit](https://qchemlab.knu.ac.kr/open-quantum-platform/openqp-devkit)
(private), which was split out of this repository with history preserved.

This is not housekeeping. GitHub refuses to serve a diff past 20,000 changed
lines (`406 too_large`), and Codex review reads that diff, so an oversized pull
request gets no review and **no error message** — #405, at 24,955 changed lines,
asked six times and was answered zero times. Nearly half of that PR was
development scaffolding: a 2,883-line derivation, eight validation gate scripts,
and ten tests mirroring them.

- `docs/` is an allowlist, like `tools/`: a document stays only if a test reads
  it as part of what that test checks. User documentation goes to openqp-docs
  (rule 4); design and method notes go to openqp-devkit.
- `tools/` is an allowlist. Every entry must name the path in this repository
  that consumes it.
- Markdown at the repository root is limited to the files a newcomer needs plus
  this one.

**Enforced by CI:** `tools/check_repo_layout.py` (the `PR policy` workflow)
checks the **whole tree**, not just the files a PR adds, so the split cannot
erode one merge at a time. As with rule 1, CI runs the trusted base-branch copy
of the script. If the build or CI genuinely needs a new script, add it to
`TOOLS_ALLOWED` together with its consumer.

**Reviewer check:** the gate is a dumb allowlist. Judgment calls are yours — is
a new `docs/`-style page user documentation or a design note, is a new test a
regression test or a mirror of a development-time gate? Say so.

### 6. A new two-electron gradient digest handles a spherical basis, and is tested in one

`grd2` drives every `get_density` with **Cartesian** shell extents. A digest
that indexes its densities with `basis%ao_offset` / `basis%naos`, which count
the *actual* AOs, therefore reads the wrong elements whenever the basis is
spherical — and only then, because in a Cartesian basis the two index spaces
coincide.

That failure is completely silent in a Pople basis. The analytic MRSF NAC
digest shipped without the branch and was wrong for every spherical basis: with
d functions (5 vs 6) it read neighbouring AO elements, broke molecular
symmetry, and turned the 1e-13 run-to-run noise of the threaded SCF into an
O(1) change in the answer; with f functions (7 vs 10) it wrote past the end of
the array and aborted. Every analytic-NAC example and test in the repository
used `6-31G` or `6-31G*`, so CI was green throughout.

Under `HARMONIC_ACTIVE`, build Cartesian-effective, `bfnrm`-folded copies of
the densities and address them at Cartesian offsets. See
`grd2_mrsf_build_cart` and the `usecart` branch of
`grd2_mrsf_compute_data_t_get_density` in `source/modules/tdhf_mrsf_gradient.F90`.

**Enforced by CI:** `tools/check_digest_harmonic.py` (the `PR policy` workflow)
requires the `get_density` bound to each concrete `grd2_compute_data_t`
extension to mention `HARMONIC_ACTIVE`. The check is **per type**, not per
file: the NAC digest above was added to a file whose other digest already
handled the spherical case. It follows the inheritance transitively, so a
digest that extends an intermediate subtype — `grd2_rhf_compute_data_t` under
the abstract `grd2_hf_compute_data_t`, say — is examined like any other. As
with rules 1 and 5, CI runs the trusted base-branch copy of the script.

**Reviewer check:** the gate proves only that the question was asked, never
that the branch is right. Any change touching AO-indexed code should be
exercised in a spherical basis, not just a Pople one — `cc-pVDZ` reaches the d
mismatch and `cc-pVTZ` the f one. Symmetry makes a good detector and needs no
reference value: see `tests/test_mrsf_nac_spherical_basis.py`.


### 7. Initialization and resource cleanup are part of correctness

For every changed allocation, pointer, persistent buffer, native handle, file or
worker, review its complete lifetime: creation, use, replacement, normal return,
early return and error exit. This applies to Fortran, C/C++, Python and bundled
library patches, not only to the numerical kernel that uses them.

- Initialize every value before it is read, including optional-branch outputs,
  allocation descriptors, pointer association and module/SAVE state. Declaration
  initialization of a Fortran local implies SAVE; do not use it for per-call state.
- Identify the owner of each allocation and borrowed view. After replacing or
  erasing a backing record, reacquire pointers; never retain a NumPy view past
  its native owner's lifetime. Match allocator/deallocator, alignment and ABI.
- Cleanup must release owned resources once, permit safe repeated cleanup where
  exposed, and preserve borrowed resources. Close a file only when the current
  operation opened it. Restore temporary settings and modified molecular data
  on recoverable failures as well as successful returns.
- Size persistent buffers from the current molecule, basis and state count;
  invalidate or resize them when that identity changes. Exercise repeated calls
  in one process, increasing and decreasing dimensions, and relevant early/error
  paths. One successful fresh-process calculation does not establish safety.
- Memory fixes need a regression that fails for the original defect, plus a
  suitable native memory/bounds check when feasible. Record compiler, dimensions,
  result and any remaining diagnostics. Do not silence a memory error, enable an
  allocator workaround, or loosen numerical tolerances to obtain a passing CI.
- Patches to bundled libraries must invalidate only the affected cache entries;
  validate allocation and destruction with the patched dependency actually linked.

**Reviewer check:** report PASS with the relevant tests, NOT APPLICABLE with a
reason, or identify the missing evidence. Static pattern checks and an ordinary
CI pass alone cannot prove initialization or resource lifetime correctness.
