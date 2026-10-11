#!/usr/bin/env python3
"""Choose the pytest files a pull request runs; print ``ALL`` for the full suite.

CI used to run the whole tests/ suite (about 3100 tests, ~21 min on a hosted
runner) on every pull request, after the example suite in the same job.  Pull
requests now run the tests that the change can affect plus a fixed smoke set,
and the full suite runs after every merge to main (push trigger in CI.yml).

Selection, from ``git diff --name-only BASE...HEAD``:

* changed test files ``tests/test_*.py`` that still exist, and every test file
  that names one (tests import helpers from other test modules);
* for each changed Python module under ``pyoqp/oqp``, every test file whose
  source names the module (``oqp.a.b`` or its last component as a word);
* for each changed native source under ``source/`` or ``include/``, every test
  file that names the file stem (e.g. ``tdhf_energy``);
* for each changed document (Markdown, docs/, LICENSE, .gitignore, CITATION),
  every test file that names it, e.g. tests/test_readme_log_metadata.py for
  README.md;
* the smoke set below, always.

The full suite is chosen instead when the change touches what every test
depends on: the test harness (conftest, shared helpers, known_failures), this
script or the workflows, the build (CMakeLists, pyproject, cmake/, external/),
or when the selection would exceed ``MAX_SELECTED`` files anyway.
"""
import ast
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "tests"
MAX_SELECTED = 120

# Fast files that exercise the installed native library end to end
# (SCF, TDHF, DFT buffers, ECP, Hessian bindings, input handling, exit status).
SMOKE = [
    "tests/test_failed_run_exit_status.py",
    "tests/test_oqp_input_schema_manifest.py",
    "tests/test_tdhf_input_checker.py",
    "tests/test_bounds_regressions_scf_tdhf.py",
    "tests/test_default_guess_huckel.py",
    "tests/test_dft_buffer_lifetime.py",
    "tests/test_gradient_buffer_cleanup.py",
    "tests/test_ecp_integrals.py",
    "tests/test_analytic_hessian_bindings.py",
]

# Paths every test depends on.  Non-test modules directly in tests/ (shared
# helpers such as oqp_module_isolation.py) are handled in needs_full_suite().
FULL_SUITE_PATTERNS = [
    r"^tests/conftest\.py$",
    r"^tests/known_failures\.txt$",
    r"^\.github/",
    r"^CMakeLists\.txt$",
    r"^pyproject\.toml$",
    r"^cmake/",
    r"^external/",
    r"^pyoqp/requirements",
]


def changed_files(base):
    """Changed paths, or None when the diff cannot be computed (e.g. base not fetched)."""
    out = subprocess.run(["git", "diff", "--name-only", f"{base}...HEAD"],
                         cwd=ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        print(f"select_pr_tests: git diff {base}...HEAD failed: {out.stderr.strip()}",
              file=sys.stderr)
        return None
    return [line for line in out.stdout.splitlines() if line]


def needs_full_suite(path):
    if path.startswith("tests/") and "/" not in path[len("tests/"):]:
        name = path[len("tests/"):]
        if name.endswith(".py") and not name.startswith("test_"):
            return True
    return any(re.search(p, path) for p in FULL_SUITE_PATTERNS)


def test_sources():
    return {p: p.read_text(errors="ignore") for p in sorted(TESTS.glob("test_*.py"))}


def tests_naming(words, sources):
    hits = set()
    for word in words:
        rx = re.compile(r"(?<![\w.])" + re.escape(word) + r"(?![\w])")
        for path, text in sources.items():
            if rx.search(text):
                hits.add(path)
    return hits


def public_names(path):
    """Top-level public names a module defines (classes, functions, assignments).

    Tests often reach a module through a package re-export
    (``from oqp.interop import MRSFExcitedStates``) without naming the module,
    so the names it defines are searched for as well.
    """
    try:
        tree = ast.parse((ROOT / path).read_text(errors="ignore"))
    except (OSError, SyntaxError):
        return set()
    names = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names |= {t.id for t in node.targets if isinstance(t, ast.Name)}
    return {n for n in names if not n.startswith("_")}


# Documentation and repository metadata.  They run no code, but tests read some
# of them (README.md, docs/logging.md, .gitignore), so a change selects every
# test that names the file and is test-neutral only when none does.
DOCUMENT_PATHS = (r"\.md$", r"^docs/", r"^LICENSE", r"^\.gitignore$", r"^CITATION")


def select(paths):
    """Return the selected test files, or None for the full suite.

    Conservative: a changed path that is neither a document nor mapped to at
    least one test file selects the full suite; a document selects the tests
    that read it.
    """
    if any(needs_full_suite(p) for p in paths):
        return None
    sources = test_sources()
    selected = set()
    for p in paths:
        if any(re.search(rx, p) for rx in DOCUMENT_PATHS):
            selected |= tests_naming({p, Path(p).name}, sources)
            continue
        if re.match(r"^tests/test_[^/]+\.py$", p):
            if (ROOT / p).exists():
                selected.add(ROOT / p)
            # Test modules also serve helpers to other tests
            # (``from test_symmetry_metadata import load_molecule_module``).
            selected |= tests_naming({Path(p).stem}, sources)
            continue
        if p.startswith("pyoqp/oqp/") and p.endswith(".py"):
            mod = p[len("pyoqp/"):-len(".py")].replace("/", ".")
            if mod.endswith(".__init__"):
                mod = mod[: -len(".__init__")]
            words = {mod, mod.rsplit(".", 1)[-1]} | public_names(p)
        elif p.startswith(("source/", "include/", "tests/")):
            # native sources, and test data/drivers in tests/ subdirectories
            # (tests/fortran/x.F90 is driven by the test file that names x)
            words = {Path(p).stem}
        else:
            return None  # unclassified path (Docker, lock files, examples, ...)
        hits = tests_naming(words, sources)
        if not hits:
            return None  # nothing maps to it: run everything
        selected |= hits
    selected |= {ROOT / s for s in SMOKE if (ROOT / s).exists()}
    if len(selected) > MAX_SELECTED:
        return None
    return sorted(str(s.relative_to(ROOT)) for s in selected)


def main():
    base = os.environ.get("BASE_SHA") or (sys.argv[1] if len(sys.argv) > 1 else "")
    if not base:
        print("ALL")
        return
    paths = changed_files(base)
    # Fail safe: without a reliable diff, run everything.
    chosen = None if paths is None else select(paths)
    print("ALL" if chosen is None else "\n".join(chosen))


if __name__ == "__main__":
    main()
