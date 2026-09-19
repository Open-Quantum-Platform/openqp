#!/usr/bin/env python3
"""Rule 5 gate: keep openqp to product code and what the build/CI depends on.

openqp accumulated method notes, derivations, validation harnesses and one-off
diagnostics that share no dependency with the engine. Two concrete costs:

  * Pull requests grew past GitHub's 20,000-line diff limit. The diff API then
    answers 406 too_large, and Codex review -- which reads that diff -- produces
    nothing at all, with no error on the pull request. #405 (24,955 changed
    lines) got six "@codex review" requests and never answered one of them.
  * The handful of scripts CI genuinely needs sat among a dozen that nothing
    referenced, so neither readers nor tooling could tell them apart.

That material now lives in Open-Quantum-Platform/openqp-devkit, with history.
This gate keeps it from coming back. It is deliberately dumb: an allowlist, no
heuristics. Judgment calls -- is this new docs/ page user documentation or a
design note, is this test a regression test or a mirror of a development gate --
belong to the Claude review, which comments without blocking.

Finding consumers: search for the file STEM, not the path. Two rounds of CI
failures on the split PR came from grepping "tools/<name>" and "docs/<name>",
which matches neither `from tools.diagnostics.trace_namd_hop import ...` nor
`ROOT / "docs" / "logging.md"`. The exhaustive query is for the directory name
as a quoted path component: `["\']docs["\']`, `["\']tools["\']`, plus the
slash forms.

Run: python3 tools/check_repo_layout.py [repo_root]
Exit 0 when the tree conforms, 1 otherwise. Reads the tree as text; builds
nothing and executes nothing from it.
"""

from __future__ import annotations

import sys
from pathlib import Path

DEVKIT = "Open-Quantum-Platform/openqp-devkit"

# Every entry here must have a demonstrable consumer inside this repository.
# Adding one without a consumer is how the previous accumulation started, so
# the comment naming that consumer is part of the entry, not decoration.
TOOLS_ALLOWED = {
    "check_blas_wrapper.py": ".github/workflows/pr-policy.yml (rule 1)",
    "check_repo_layout.py": ".github/workflows/pr-policy.yml (rule 5, this gate)",
    "convert_legacy_examples.py": "tests/known_failures.txt",
    "generate_int2_pure_kernels.py": "source/integrals/int2_pure_generated.F90",
    "minao": "source/minao_lut.F90",
    "sap": "source/sap_lut.F90",
    "scf-converger-ml": ".github/workflows/train-scf-selector.yml",
    # Loaded by tests, not by the build. Found the hard way: the first pass of
    # this split searched for the literal path "tools/<name>" and so missed a
    # Python import and an importlib load built from path parts, which broke CI.
    "diagnostics": "tests/test_trace_namd_hop.py, test_namd_baeck_an.py, test_namd_diagnostics.py",
    "validate_analytic_hessian.py": "tests/test_analytic_hessian_validator.py",
    # Code generators whose output is committed. The build does not invoke them,
    # but regenerating Boys/Rys tables or the libxc bindings requires them to
    # sit beside the source they produce.
    "gen_boys.F90": "regenerates the Boys table in source/",
    "gen_rys.F90": "regenerates the Rys roots in source/",
    "parallel_gen.fpp": "regenerates OpenMP scaffolding in source/",
    "libxc": "regenerates the libxc bindings in source/",
}

# Markdown at the repository root is for people arriving at the repository, plus
# AGENTS.md, which states the contribution rules this workflow enforces.
ROOT_MD_ALLOWED = {
    "README.md",
    "CONTRIBUTING.md",
    "AGENTS.md",
    "LICENSING.md",
    "SUSTAINABILITY.md",
    "THIRD_PARTY_NOTICES.md",
    "TROUBLESHOOTING.md",
}

# docs/ follows the same rule as tools/: a file stays only if something in this
# repository reads it. Two tests read their subject document directly, so those
# two documents are part of the test, not prose that happens to sit nearby.
# Everything else -- method notes, design documents, release notes -- moved.
DOCS_ALLOWED = {
    "TDDFT_HESSIAN_IMPORT.md": "tests/test_tdhf_hessian_import.py (ROOT / \"docs\" / ...)",
    "logging.md": "tests/test_log_format.py (ROOT / \"docs\" / ... .read_text())",
}


def check(root: Path) -> list[str]:
    problems: list[str] = []

    docs = root / "docs"
    if docs.is_dir():
        for path in sorted(docs.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(docs).as_posix()
            if rel in DOCS_ALLOWED:
                continue
            problems.append(
                f"docs/{rel} has no consumer in this repository.\n"
                f"    Method notes, design documents and release notes belong in\n"
                f"    {DEVKIT}; user documentation belongs in openqp-docs.\n"
                f"    If a test reads this file, add it to DOCS_ALLOWED in\n"
                f"    tools/check_repo_layout.py together with that test."
            )

    tools = root / "tools"
    if tools.is_dir():
        for entry in sorted(p.name for p in tools.iterdir()):
            if entry in TOOLS_ALLOWED or entry == "__pycache__":
                continue
            problems.append(
                f"tools/{entry} has no consumer in this repository.\n"
                f"    Development tooling belongs in {DEVKIT}.\n"
                f"    If the build or CI genuinely needs it, add it to "
                f"TOOLS_ALLOWED in tools/check_repo_layout.py together with the "
                f"path that consumes it."
            )

    for entry in sorted(p.name for p in root.iterdir() if p.is_file()):
        if entry.endswith(".md") and entry not in ROOT_MD_ALLOWED:
            problems.append(
                f"{entry} is a development note at the repository root.\n"
                f"    Investigation and performance notes belong in {DEVKIT}."
            )

    return problems


def main(argv: list[str]) -> int:
    root = Path(argv[1] if len(argv) > 1 else ".").resolve()
    if not (root / "CMakeLists.txt").is_file():
        print(f"error: {root} does not look like the openqp root", file=sys.stderr)
        return 2

    problems = check(root)
    if not problems:
        print("Repository layout gate (rule 5): PASS")
        return 0

    print("Repository layout gate (rule 5): FAIL\n")
    for problem in problems:
        print(f"  ✗ {problem}\n")
    print(
        f"openqp holds product code and what the build and CI depend on.\n"
        f"Everything else goes to {DEVKIT}, which preserves history."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
