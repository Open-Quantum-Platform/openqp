#!/usr/bin/env python3
"""Run a NAC command in a fresh work directory and reject a stray fort.6."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run COMMAND in a temporary working directory and fail if OpenQP "
            "creates the process-global fallback log fort.6. Use absolute paths "
            "for command arguments that name files."
        )
    )
    parser.add_argument(
        "command",
        nargs=argparse.REMAINDER,
        help="command to run, conventionally preceded by --",
    )
    args = parser.parse_args()

    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a command is required after --")

    with tempfile.TemporaryDirectory(prefix="openqp-nac-fort6-") as tmp:
        workdir = Path(tmp)
        print(f"fort.6 guard work directory: {workdir}", flush=True)
        result = subprocess.run(command, cwd=workdir, check=False)
        fallback_log = workdir / "fort.6"

        if fallback_log.exists():
            size = fallback_log.stat().st_size
            print(
                f"FAIL: command created {fallback_log.name} ({size} bytes)",
                file=sys.stderr,
            )
            return 86
        if result.returncode != 0:
            print(
                f"FAIL: guarded command exited with status {result.returncode}",
                file=sys.stderr,
            )
            return result.returncode

    print("PASS: guarded command completed without creating fort.6")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
