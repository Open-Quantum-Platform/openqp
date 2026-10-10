"""Build the PyInstaller-frozen backend binary for the Tauri sidecar.

Run from backend/ with the frontend already built (frontend/dist present):

    pip install -e . pyinstaller
    python build_binary.py

On macOS, writes dist/oqp-studio-backend/ with the server and its runtime
libraries. Keeping that runtime unpacked avoids PyInstaller extracting a large
archive before every launch. Other platforms keep the existing single-file
sidecar until their installers can stage the runtime directory beside it.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def main() -> int:
    backend = Path(__file__).resolve().parent
    frontend_dist = backend.parent / "frontend" / "dist"
    if not (frontend_dist / "index.html").is_file():
        print("frontend/dist not found — run `npm run build` in frontend/ first.")
        return 1

    sep = ";" if os.name == "nt" else ":"
    try:
        import rdkit  # noqa: F401

        rdkit_args = ["--collect-all", "rdkit"]
    except ImportError:
        rdkit_args = []
        print("note: rdkit not installed — sketch-to-3D disabled in this binary")
    cmd = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--noconfirm",
        "--onedir" if sys.platform == "darwin" else "--onefile",
        "--name",
        "oqp-studio-backend",
        "--add-data",
        f"{frontend_dist}{sep}frontend_dist",
        # uvicorn loads its loop/protocol classes by string name at runtime.
        "--hidden-import",
        "uvicorn.logging",
        "--collect-submodules",
        "uvicorn",
        # truststore is imported lazily and picks its backend by platform, so
        # collect it whole rather than relying on the import graph.
        "--collect-all",
        "truststore",
        *rdkit_args,
        str(backend / "oqp_studio" / "server_main.py"),
    ]
    print(" ".join(cmd))
    return subprocess.call(cmd, cwd=backend)


if __name__ == "__main__":
    raise SystemExit(main())
