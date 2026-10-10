"""Keep application source identity separate from packaging-tool revisions."""
from __future__ import annotations

import re
import subprocess
from pathlib import Path


def source_commits(root: Path, application_commit: str | None = None) -> tuple[str, str]:
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=root, text=True).strip()

    packaging = git("rev-parse", "HEAD")
    if git("status", "--porcelain", "--untracked-files=no"):
        raise ValueError("commit and review the candidate before packaging")
    application = application_commit or packaging
    if not re.fullmatch(r"[0-9a-f]{40}", application):
        raise ValueError("application commit must be a full Git SHA")
    subprocess.run(["git", "merge-base", "--is-ancestor", application, packaging],
                   cwd=root, check=True)
    changed = git("diff", "--name-only", application, packaging).splitlines()
    # These files are build tools/tests only. Every other tracked file, including
    # app code, assets, dependencies and engine identity, must remain identical.
    if any(not (p.startswith(("tools/release/", "tests/")) or p == ".gitlab-ci.yml")
           for p in changed):
        raise ValueError("application files changed since the requested source commit")
    return application, packaging
