#!/usr/bin/env python3
"""Choose the pytest files a pull request runs; print ``ALL`` for the full suite.

CI used to run the whole tests/ suite (about 3100 tests, ~21 min on a hosted
runner) on every pull request, after the example suite in the same job.  Pull
requests now run the tests that the change can affect plus a fixed smoke set,
and the full suite runs after every merge to main (push trigger in CI.yml).

Selection, from ``git diff --name-only BASE...HEAD``:

* changed test files ``tests/test_*.py`` that still exist;
* for each changed Python module under ``pyoqp/oqp``, every test file whose
  source names the module (``oqp.a.b`` or its last component as a word);
* for each changed native source under ``source/`` or ``include/``, every test
  file that names the file stem (e.g. ``tdhf_energy``);
* the smoke set below, always.

The full suite is chosen instead when the change touches what every test
depends on: the test harness (conftest, shared helpers, known_failures), this
script or the workflows, the build (CMakeLists, pyproject, cmake/, external/),
or when the selection would exceed ``MAX_SELECTED`` files anyway.
"""
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
    out = subprocess.run(["git", "diff", "--name-only", f"{base}...HEAD"],
                         cwd=ROOT, check=True, capture_output=True, text=True)
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


def select(paths):
    if any(needs_full_suite(p) for p in paths):
        return None
    sources = test_sources()
    selected = set()
    for p in paths:
        if re.match(r"^tests/test_[^/]+\.py$", p) and (ROOT / p).exists():
            selected.add(ROOT / p)
        elif p.startswith("pyoqp/oqp/") and p.endswith(".py"):
            mod = p[len("pyoqp/"):-len(".py")].replace("/", ".")
            if mod.endswith(".__init__"):
                mod = mod[: -len(".__init__")]
            words = {mod, mod.rsplit(".", 1)[-1]}
            selected |= tests_naming(words, sources)
        elif p.startswith(("source/", "include/")):
            selected |= tests_naming({Path(p).stem}, sources)
    selected |= {ROOT / s for s in SMOKE if (ROOT / s).exists()}
    if len(selected) > MAX_SELECTED:
        return None
    return sorted(str(s.relative_to(ROOT)) for s in selected)


def main():
    base = os.environ.get("BASE_SHA") or (sys.argv[1] if len(sys.argv) > 1 else "")
    if not base:
        print("ALL")
        return
    chosen = select(changed_files(base))
    print("ALL" if chosen is None else "\n".join(chosen))


if __name__ == "__main__":
    main()
