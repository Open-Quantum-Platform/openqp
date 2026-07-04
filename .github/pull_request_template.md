<!-- Describe the change and link any related issues. -->

## Summary

## Contribution checklist

See [AGENTS.md](../AGENTS.md) for details. The automated Claude and Codex
reviewers check these too.

- [ ] **BLAS/LAPACK** — every referenced BLAS/LAPACK routine is registered in
      `OQP_ACCELERATE_ILP64_SYMS` (`cmake/oqp_functions.cmake`); new code prefers
      `use oqp_linalg`. Unregistered routines fail the `PR policy` lint.
- [ ] **Example** — new functionality (runtype/method/`[section]`/keyword/flag)
      adds a small, fast example under `examples/`. (New opt-in **boolean** flags
      are gated by `openqp --check_feature_coverage`; other additions are checked
      by the automated reviewers.)
- [ ] **Python API** — user-facing additions are usable from `oqp.openqp.OpenQP`
      and covered by a `tests/test_openqp_api.py` test.
- [ ] **Docs** — user-facing keywords/workflows documented in openqp-docs; link
      the companion PR here: <!-- openqp-docs PR # -->

<!-- If a checklist item does not apply, say why. -->
