#!/usr/bin/env python3
"""Summarize post-OV-OV-fix MRSF gradient finite-difference diagnostics.

The helper is intentionally dependency-light so cron/autopilot runs can rank
existing validation artifacts before launching any new quantum chemistry jobs.
"""

from __future__ import annotations

import argparse
import ast
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


def _to_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def _parse_s2_map(value: Any) -> dict[int, float]:
    if value in (None, ""):
        return {}
    try:
        parsed = ast.literal_eval(str(value))
    except (SyntaxError, ValueError):
        return {}
    if not isinstance(parsed, dict):
        return {}
    result: dict[int, float] = {}
    for key, item in parsed.items():
        try:
            result[int(key)] = float(item)
        except (TypeError, ValueError):
            continue
    return result


def _s2_max_delta(*maps: dict[int, float]) -> float:
    if not maps:
        return 0.0
    keys = set().union(*(mapping.keys() for mapping in maps))
    max_delta = 0.0
    for key in keys:
        values = [mapping[key] for mapping in maps if key in mapping]
        if values:
            max_delta = max(max_delta, max(values) - min(values))
    return max_delta


def _s2_evidence(*maps: dict[int, float]) -> str:
    populated = [mapping for mapping in maps if mapping]
    if not populated:
        return "unknown"
    values = [abs(value) for mapping in populated for value in mapping.values()]
    if values and max(values) == 0.0:
        return "unknown"
    return "present"


def normalize_component_row(row: dict[str, str], threshold: float = DEFAULT_THRESHOLD) -> dict[str, Any]:
    root = _to_int(row.get("root"))
    s2_grad = _parse_s2_map(row.get("s2_grad"))
    s2_plus = _parse_s2_map(row.get("s2_plus"))
    s2_minus = _parse_s2_map(row.get("s2_minus"))
    abs_diff = _to_float(row.get("abs_diff_ha_per_bohr"))
    return {
        "method": row.get("method", ""),
        "molecule": str(row.get("molecule", "")).lower(),
        "root": root,
        "physical_state": physical_state_for_mrsf_root(root, row.get("physical_state") or None),
        "component": row.get("component", ""),
        "axis": str(row.get("component", ""))[-1:] or "?",
        "analytic_ha_per_bohr": _to_float(row.get("analytic_ha_per_bohr")),
        "fd_ha_per_bohr": _to_float(row.get("fd_ha_per_bohr")),
        "diff_ha_per_bohr": _to_float(row.get("diff_ha_per_bohr")),
        "abs_diff_ha_per_bohr": abs_diff,
        "trah_count": _to_int(row.get("trah_count")),
        "failed_any": _to_bool(row.get("failed_any")),
        "s2_max_delta": _s2_max_delta(s2_grad, s2_plus, s2_minus),
        "s2_evidence": _s2_evidence(s2_grad, s2_plus, s2_minus),
        "bad_component": abs_diff > threshold,
        "target_case": (str(row.get("molecule", "")).lower(), root) in TARGET_CASES,
    }


def load_components_csv(path: Path, threshold: float = DEFAULT_THRESHOLD) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return [normalize_component_row(row, threshold) for row in csv.DictReader(handle)]


def _component_mechanism_hint(group: dict[str, Any]) -> str:
    if group["trah_or_failed"]:
        return "not-clean: inspect SCF/TRAH before gradient algebra"
    if group["possible_state_character_change"]:
        return "root_tracking_or_state_character_change"
    if group["s2_evidence"] == "unknown":
        return "state-character evidence missing; rerun/parse root-continuity metadata before algebra edits"
    if group["bad_component_count"] == 1 and group["worst_axis"] == "z":
        return "localized_z_component_z_vector_or_operator_mapping"
    return "compare_spc_toggles_xc_density_handoff_and_z_vector_mapping"


def summarize_component_rows(rows: Iterable[dict[str, Any]], threshold: float = DEFAULT_THRESHOLD) -> dict[str, Any]:
    grouped: dict[tuple[str, str, int], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault((row["method"], row["molecule"], row["root"]), []).append(row)

    summaries: list[dict[str, Any]] = []
    for (method, molecule, root), items in grouped.items():
        worst = max(items, key=lambda item: item["abs_diff_ha_per_bohr"])
        bad_items = [item for item in items if item["abs_diff_ha_per_bohr"] > threshold]
        bad_items.sort(key=lambda item: (-item["abs_diff_ha_per_bohr"], item["component"]))
        axes = sorted({item["axis"] for item in bad_items})
        bad_components = [
            {
                "component": item["component"],
                "axis": item["axis"],
                "abs_diff_ha_per_bohr": item["abs_diff_ha_per_bohr"],
                "analytic_ha_per_bohr": item["analytic_ha_per_bohr"],
                "fd_ha_per_bohr": item["fd_ha_per_bohr"],
                "s2_evidence": item["s2_evidence"],
            }
            for item in bad_items
        ]
        summary = {
            "method": method,
            "molecule": molecule,
            "root": root,
            "physical_state": worst["physical_state"],
            "max_abs_diff_ha_per_bohr": worst["abs_diff_ha_per_bohr"],
            "worst_component": worst["component"],
            "worst_axis": worst["axis"],
            "bad_component_count": len(bad_items),
            "bad_axes": axes,
            "bad_components": bad_components,
            "s2_max_delta": max(item["s2_max_delta"] for item in items),
            "s2_evidence": "present" if any(item.get("s2_evidence") == "present" for item in items) else "unknown",
            "possible_state_character_change": max(item["s2_max_delta"] for item in items) > 0.5,
            "trah_or_failed": any(item["trah_count"] > 0 or item["failed_any"] for item in items),
            "target_case": any(item.get("target_case", (molecule, root) in TARGET_CASES) for item in items),
        }
        summary["mechanism_hint"] = _component_mechanism_hint(summary)
        summaries.append(summary)

    summaries.sort(
        key=lambda item: (
            not item["possible_state_character_change"],
            -item["max_abs_diff_ha_per_bohr"],
            item["molecule"],
            item["root"],
        )
    )
    target_groups = [item for item in summaries if item["target_case"]]
    target_bad_groups = [item for item in target_groups if item["bad_component_count"] > 0 or item["trah_or_failed"]]
    return {
        "threshold_ha_per_bohr": threshold,
        "component_group_count": len(summaries),
        "target_group_count": len(target_groups),
        "target_bad_group_count": len(target_bad_groups),
        "target_bad_groups": target_bad_groups,
        "groups": summaries,
    }


def summarize_components_csv(path: Path | str, threshold: float = DEFAULT_THRESHOLD) -> dict[str, Any]:
    return summarize_component_rows(load_components_csv(Path(path), threshold), threshold)


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
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--components", action="store_true", help="Summarize per-component FD diagnostics CSV")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--output", type=Path, help="Write JSON summary to this path")
    args = parser.parse_args(argv)

    if args.components:
        summary = summarize_components_csv(args.csv_path, args.threshold)
    else:
        summary = summarize_compact_csv(args.csv_path, args.threshold)
    payload = json.dumps(summary, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
