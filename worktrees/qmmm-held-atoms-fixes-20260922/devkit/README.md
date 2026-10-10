# OpenQP private development kit

This directory contains private development material beside the OpenQP engine
in the canonical GitLab internal superset.

Nothing here is built, imported, or executed by OpenQP, and nothing in the
engine depends on it. It is kept because it is expensive to
reconstruct: derivations, validation harnesses, performance investigations and
design notes that explain *why* the engine looks the way it does.

## Why this repository exists

`openqp` accumulated documentation and one-off tooling that shares no
dependency with the product. That material made every pull request larger than
GitHub's 20,000-line diff limit, which silently disables Codex review, and it
buried the small set of scripts that CI genuinely needs.

The material was moved to a separate `openqp-devkit` repository to protect the
public GitHub review path. On 2026-09-21 its history was brought back under
`devkit/` in the private GitLab superset. `tools/check_repo_layout.py` still
keeps the engine directories clean, and public promotion must exclude this
directory completely.

## Layout

| Path | Contents |
|---|---|
| `docs/` | Method notes, design documents, release notes, `planned/` proposals |
| `tools/` | Gradient validation, cross-checks, example gradient checks |
| `tools/nac_lagrangian/` | MRSF NAC derivation (2,883 lines) and its validation gates |
| `tests/` | The gates' own tests, plus source-shape contracts read out of an `openqp` checkout |
| `GRAD_SCREENING_NOTES.md` | Gradient screening investigation |
| `MRSF_ZVECTOR_PERF_NOTES.md` | MRSF Z-vector performance investigation |

Git history is preserved: every file here carries the commits that produced it,
so `git log --follow <path>` still works.

## Running the tests

```
python -m pytest devkit/tests/
```

Some of the tests read engine source rather than running a calculation: they
assert that a guard precedes a division, that a subroutine keeps no state
between calls, that a loop nest is ordered a particular way. Those are
development gates, not engine regressions -- one of them failed on a pure hoist
that changed no arithmetic -- which is why they live here and not in `openqp`'s
own suite.

They read the engine source from the repository root automatically. A standalone
historical devkit checkout can still set `OPENQP_SOURCE_ROOT` or sit beside an
`openqp/` checkout. Without a usable engine tree they skip: there is nothing to
assert against, which is not the same as an assertion being violated.
`tests/openqp_checkout.py` performs the lookup. One test also needs a *built*
engine carrying the analytic NAC work and skips when the installed `oqp`
predates it.

## What did NOT move

These stay in `openqp` because the build or CI depends on them:
`tools/check_blas_wrapper.py` (PR policy rule 1), `tools/convert_legacy_examples.py`,
`tools/generate_int2_pure_kernels.py`, `tools/minao/`, `tools/sap/`,
`tools/scf-converger-ml/`, and the code generators `tools/gen_boys.F90`,
`tools/gen_rys.F90`, `tools/parallel_gen.fpp`, `tools/libxc/`.

Two more stayed for a reason worth recording: `tools/diagnostics/` and
`tools/validate_analytic_hessian.py` are loaded by tests
(`test_trace_namd_hop.py`, `test_namd_baeck_an.py`, `test_namd_diagnostics.py`,
`test_analytic_hessian_validator.py`). A Python import spells the path with
dots and an importlib load builds it from parts, so a search for the literal
string `tools/<name>` finds neither. Search for the stem, not the path.

## Cross-references

Engine comments should name development material as `devkit/docs/...` or
`devkit/tools/...` so the target resolves inside the private superset.

## Contributing

Changes go through a GitLab merge request in the internal OpenQP project. They
remain private unless a separate public promotion explicitly selects them.
