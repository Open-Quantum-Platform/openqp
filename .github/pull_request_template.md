<!-- Describe the change and link any related issues. -->

## Summary

## Contribution checklist

See [AGENTS.md](../AGENTS.md) for details. The automated Claude and Codex
reviewers check these too.

- [ ] **BLAS/LAPACK** — no raw `dgemm`/`dgels`/`dgetrf`/... calls; all go through
      the `oqp_linalg` wrapper (raw calls only inside
      `source/mathlib/lapack_wrap.F90` / `blas_wrap.F90`).
- [ ] **Example** — new functionality (runtype/method/`[section]`/keyword/flag)
      adds a small, fast example under `examples/` (passes
      `openqp --check_feature_coverage`).
- [ ] **Python API** — user-facing additions are usable from `oqp.openqp.OpenQP`
      and covered by a `tests/test_openqp_api.py` test.
- [ ] **Docs** — user-facing keywords/workflows documented in openqp-docs; link
      the companion PR here: <!-- openqp-docs PR # -->

<!-- If a checklist item does not apply, say why. -->
