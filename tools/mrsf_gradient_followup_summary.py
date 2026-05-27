#!/usr/bin/env python3
"""Summarize post-OV-OV-fix MRSF gradient finite-difference diagnostics.

The helper is intentionally dependency-light so cron/autopilot runs can rank
existing validation artifacts before launching any new quantum chemistry jobs.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Iterable

DEFAULT_THRESHOLD = 1.0e-3
TARGET_CASES = {
    ("ch2o", 4),
    ("ch2o", 5),
    ("ch2o", 7),
    ("h2s", 5),
    ("nh3", 3),
    ("nh3", 4),
    ("hcn", 5),
    ("hcn", 6),
    ("h2o", 4),
}


def _to_float(value: Any, default: float = 0.0) -> float:
    if value in (None, ""):
        return default
    return float(value)


def _to_int(value: Any, default: int = 0) -> int:
    if value in (None, ""):
        return default
    return int(float(value))


def physical_state_for_mrsf_root(root: int, explicit: str | None = None) -> str | None:
    """Map OpenQP MRSF response roots to physical singlets.

    Root 0 is the ROHF high-spin reference and root 1 is physical S0, so root N
    maps to S(N-1).  Compact CSV artifacts often leave physical_state blank.
    """

    if explicit:
        return explicit
    if root <= 0:
        return None
    return f"S{root - 1}"


def _classification(row: dict[str, Any], threshold: float) -> str:
    if _to_int(row.get("trah_total")) > 0 or _to_int(row.get("failed_count")) > 0:
        return "trah_or_failed"
    if _to_float(row.get("max_abs_diff_ha_per_bohr")) > threshold:
        return "fd_mismatch_no_trah"
    return "clean"


def _mechanism_hint(row: dict[str, Any], threshold: float) -> str:
    classification = _classification(row, threshold)
    if classification == "clean":
        return "clean_control"
    if classification == "trah_or_failed":
        return "not-clean: inspect SCF/TRAH before gradient algebra"

    molecule = str(row.get("molecule", "")).lower()
    root = _to_int(row.get("root"))
    component = str(row.get("worst_component", ""))
    if (molecule, root) in {("nh3", 3), ("nh3", 4)}:
        return "prioritize root-tracking/near-degeneracy state-character check"
    if component.endswith("_z") and molecule in {"hcn", "h2s", "ch2o"}:
        return "localized axial component; compare root continuity and Z-vector/operator mapping"
    return "compare state continuity, SPC toggles, and XC/Z-vector density handoff"


def normalize_row(row: dict[str, str], threshold: float = DEFAULT_THRESHOLD) -> dict[str, Any]:
    root = _to_int(row.get("root"))
    normalized: dict[str, Any] = {
        "method": row.get("method", ""),
        "molecule": str(row.get("molecule", "")).lower(),
        "root": root,
        "physical_state": physical_state_for_mrsf_root(root, row.get("physical_state") or None),
        "max_abs_diff_ha_per_bohr": _to_float(row.get("max_abs_diff_ha_per_bohr")),
        "rms_diff_ha_per_bohr": _to_float(row.get("rms_diff_ha_per_bohr")),
        "worst_component": row.get("worst_component", ""),
        "worst_analytic": _to_float(row.get("worst_analytic")),
        "worst_fd": _to_float(row.get("worst_fd")),
        "trah_total": _to_int(row.get("trah_total")),
        "failed_count": _to_int(row.get("failed_count")),
        "elapsed_s": _to_float(row.get("elapsed_s")),
    }
    normalized["target_case"] = (normalized["molecule"], root) in TARGET_CASES
    normalized["classification"] = _classification(normalized, threshold)
    normalized["mechanism_hint"] = _mechanism_hint(normalized, threshold)
    return normalized


def load_compact_csv(path: Path, threshold: float = DEFAULT_THRESHOLD) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return [normalize_row(row, threshold) for row in csv.DictReader(handle)]


def summarize_rows(rows: Iterable[dict[str, Any]], threshold: float = DEFAULT_THRESHOLD) -> dict[str, Any]:
    all_rows = list(rows)
    failures = [
        row
        for row in all_rows
        if row["classification"] != "clean" or row["max_abs_diff_ha_per_bohr"] > threshold
    ]
    failures.sort(
        key=lambda row: (
            row["classification"] != "trah_or_failed",
            -row["max_abs_diff_ha_per_bohr"],
            row["molecule"],
            row["root"],
        )
    )
    clean_controls = [
        f"{row['molecule']} root {row['root']}"
        for row in all_rows
        if row["classification"] == "clean" and not row["target_case"]
    ]
    target_failures = [row for row in failures if row["target_case"]]
    return {
        "threshold_ha_per_bohr": threshold,
        "total_cases": len(all_rows),
        "failure_count": len(failures),
        "target_failure_count": len(target_failures),
        "trah_or_failed_count": sum(1 for row in failures if row["classification"] == "trah_or_failed"),
        "failures": failures,
        "clean_controls": clean_controls,
    }


def summarize_compact_csv(path: Path | str, threshold: float = DEFAULT_THRESHOLD) -> dict[str, Any]:
    return summarize_rows(load_compact_csv(Path(path), threshold), threshold)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("compact_csv", type=Path)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--output", type=Path, help="Write JSON summary to this path")
    args = parser.parse_args(argv)

    summary = summarize_compact_csv(args.compact_csv, args.threshold)
    payload = json.dumps(summary, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
