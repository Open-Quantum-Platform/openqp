#!/usr/bin/env python3
"""Refresh OpenQP example JSON reference files from live Runner output.

Usage examples:
  python scripts/refresh_example_json_refs.py examples/OPT/H2O_RHF-DFT_OPTIMIZE.inp
  python scripts/refresh_example_json_refs.py examples/OPT --pattern '*.inp'

The script runs each input through oqp.pyoqp.Runner, calls mol.save_data(), and
copies the produced JSON back next to the input file. It is intentionally scoped
to explicit paths; avoid refreshing the full examples tree unless that is what
you intend to review and commit.
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

from oqp.pyoqp import Runner


def discover_inputs(paths: list[Path], pattern: str) -> list[Path]:
    inputs: list[Path] = []
    for path in paths:
        if path.is_dir():
            inputs.extend(sorted(path.rglob(pattern)))
        elif path.is_file() and path.suffix == ".inp":
            inputs.append(path)
        else:
            raise FileNotFoundError(f"not an input file or directory: {path}")
    return sorted(dict.fromkeys(p.resolve() for p in inputs))


def refresh_one(inp: Path, workroot: Path, save_on_error: bool) -> tuple[str, float]:
    project = inp.stem
    work = workroot / project
    work.mkdir(parents=True, exist_ok=True)
    log = work / f"{project}.log"
    started = time.perf_counter()
    runner = Runner(project=project, input_file=str(inp), log=str(log), silent=1, usempi=False)
    try:
        runner.run(test_mod=True)
        status = "ok"
    except Exception:
        if not save_on_error:
            raise
        status = "saved_after_error"
    runner.mol.save_data()
    produced = log.with_suffix(".json")
    if not produced.exists():
        raise RuntimeError(f"Runner did not produce {produced}")
    target = inp.with_suffix(".json")
    shutil.copy2(produced, target)
    return status, time.perf_counter() - started


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path, help="Example .inp files or directories")
    parser.add_argument("--pattern", default="*.inp", help="Input glob used for directories")
    parser.add_argument("--workdir", type=Path, default=Path("refresh_example_json_tmp"))
    parser.add_argument(
        "--save-on-error",
        action="store_true",
        help="Still save the molecule JSON after a bounded optimizer raises, useful only for reviewed smoke examples.",
    )
    args = parser.parse_args(argv)

    inputs = discover_inputs(args.paths, args.pattern)
    if not inputs:
        print("No input files found", file=sys.stderr)
        return 1

    args.workdir.mkdir(parents=True, exist_ok=True)
    for inp in inputs:
        rel = inp if not inp.is_relative_to(Path.cwd()) else inp.relative_to(Path.cwd())
        print(f"RUN {rel}", flush=True)
        status, elapsed = refresh_one(inp, args.workdir, args.save_on_error)
        print(f"SAVED {inp.with_suffix('.json')} status={status} elapsed={elapsed:.2f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
