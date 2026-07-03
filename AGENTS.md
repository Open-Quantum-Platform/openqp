# Contributor and AI-reviewer guide

This file is read by the automated PR reviewers (Codex reads `AGENTS.md`
natively; the Claude review workflow is pointed at it in
`.github/workflows/claude.yml`) and by human contributors. Every pull request is
expected to satisfy the four rules below. A reviewer should call out, per rule,
whether the PR satisfies it or explain what is missing.

## PR rules

### 1. BLAS/LAPACK must go through the OpenQP wrapper layer

Never call a raw BLAS/LAPACK routine (`dgemm`, `dgels`, `dgetrf`, `dsyev`,
`zheev`, ...) directly. Fortran modules must `use oqp_linalg`, which renames the
`oqp_<name>_i64` wrappers back to the standard names, and then call the standard
name (`call dgels(...)` resolves to `oqp_dgels_i64`).

- Raw calls are only allowed **inside the wrapper modules**
  `source/mathlib/lapack_wrap.F90` and `source/mathlib/blas_wrap.F90`.
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
fails the build if any BLAS/LAPACK routine referenced from `source/` (outside the
wrapper modules) is not registered in `OQP_ACCELERATE_ILP64_SYMS` — the portable
complement to the macOS-only, link-time `cmake/check_accelerate_aliases.cmake`.

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
repo [Open-Quantum-Platform/openqp-docs](https://github.com/Open-Quantum-Platform/openqp-docs):
a keyword-page entry under `docs/keywords/` and/or a workflow page under
`docs/workflows/`, wired into `mkdocs.yml` nav.

- openqp-docs is a **separate repository**, so this PR's CI cannot see it
  directly. Link the companion openqp-docs PR in this PR's description.

**Reviewer check:** if the diff adds/changes user-facing keywords or workflows,
confirm the PR description links an openqp-docs PR; flag it if missing.
