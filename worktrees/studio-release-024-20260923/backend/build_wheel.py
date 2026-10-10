"""Build the pip-installable distribution (`pip install oqp-studio`).

Copies the built frontend into the package so one install ships the server,
the UI, and the results viewer:

    cd frontend && npm install && npm run build
    cd ../backend && python build_wheel.py
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


def main() -> int:
    backend = Path(__file__).resolve().parent
    dist = backend.parent / "frontend" / "dist"
    if not (dist / "index.html").is_file():
        print("frontend/dist not found — run `npm run build` in frontend/ first.")
        return 1

    web = backend / "oqp_studio" / "web"
    shutil.rmtree(web, ignore_errors=True)
    shutil.copytree(dist, web)
    print(f"copied frontend into {web.relative_to(backend)}")

    try:
        return subprocess.call([sys.executable, "-m", "build"], cwd=backend)
    finally:
        # The copy only exists to be packaged. Leaving it behind would shadow
        # frontend/dist for anyone running the server from a source checkout,
        # who would then keep seeing the UI as it was at the last wheel build.
        shutil.rmtree(web, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
