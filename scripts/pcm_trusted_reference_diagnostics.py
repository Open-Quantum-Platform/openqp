#!/usr/bin/env python3
"""Capture identical-protocol OpenQP/ddX PCM diagnostics as a JSON snapshot.

This is not a literature benchmark generator. It runs the same Fortran/C/ddX
path used by tests/test_pcm_literature_benchmarks.py and records the diagnostic
quantities needed to review or regenerate a future trusted-reference table.
Rows with reference_value=null remain pending-reference in the benchmark table.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"


def git_text(args: list[str]) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()
    except Exception:
        return "unknown"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", type=Path, default=ROOT / "tests" / "data" / "pcm_trusted_reference_diagnostics.json")
    args = parser.parse_args()

    sys.path.insert(0, str(TESTS))
    import test_pcm_literature_benchmarks as benchmod  # noqa: PLC0415

    records = []
    for bench in benchmod._load_benchmarks():
        proc, log = benchmod._run(benchmod._input_text(bench, pcm_on=True))
        record = {
            "id": bench["id"],
            "method": bench.get("method"),
            "functional": bench.get("functional"),
            "basis": bench.get("basis"),
            "scf_type": bench.get("scf_type"),
            "epsilon": bench.get("epsilon"),
            "returncode": proc.returncode,
            "ddx_unavailable": benchmod._ddx_unavailable(log),
            "ddpcm_nonconvergence": benchmod._ddpcm_nonconvergence(log),
            "pcm_total_energy_hartree": benchmod._total_energy(log),
            "pcm_solvent_energy_hartree": benchmod._pcm_energy(log),
            "diagnostics": benchmod._pcm_diag(log),
        }
        if proc.returncode != 0:
            record["log_tail"] = log[-2000:]
        records.append(record)

    snapshot = {
        "schema": "openqp-ddx-pcm-diagnostics-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo": str(ROOT),
        "branch": git_text(["branch", "--show-current"]),
        "source_tree_note": "Snapshot generated from the checked-out working tree; use the enclosing git commit for the committed version.",
        "ddx_prefix": os.environ.get("DDX_ROOT") or "/tmp/ddx-install",
        "environment": {
            "OPENQP_ROOT": os.environ.get("OPENQP_ROOT"),
            "DYLD_LIBRARY_PATH_set": bool(os.environ.get("DYLD_LIBRARY_PATH")),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
        },
        "status_note": (
            "Identical OpenQP/ddX protocol diagnostic snapshot for the closed-shell "
            "trusted-regression benchmark table. pcm_solvent_energy_hartree is a "
            "same-branch regression reference, not an independent literature value."
        ),
        "records": records,
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"json_out": str(args.json_out), "records": len(records)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
