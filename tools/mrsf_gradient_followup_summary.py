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
import re
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
    grouped: dict[tuple[str, str, int, str], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault((row["method"], row["molecule"], row["root"], str(row.get("source_csv", ""))), []).append(row)

    summaries: list[dict[str, Any]] = []
    for (method, molecule, root, source_csv), items in grouped.items():
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
            "source_csv": source_csv,
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


def summarize_component_datasets(
    datasets: Iterable[tuple[str, Iterable[dict[str, Any]]]],
    threshold: float = DEFAULT_THRESHOLD,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    dataset_count = 0
    for source_csv, dataset_rows in datasets:
        dataset_count += 1
        for row in dataset_rows:
            tagged = dict(row)
            tagged["source_csv"] = source_csv
            rows.append(tagged)
    summary = summarize_component_rows(rows, threshold)
    summary["dataset_count"] = dataset_count
    return summary


def summarize_components_csvs(paths: Iterable[Path | str], threshold: float = DEFAULT_THRESHOLD) -> dict[str, Any]:
    datasets = []
    for path in paths:
        csv_path = Path(path)
        datasets.append((str(csv_path), load_components_csv(csv_path, threshold)))
    return summarize_component_datasets(datasets, threshold)


_STATE_ENERGY_RE = re.compile(r"State #\s+(\d+)\s+Energy =\s+([-+0-9.Ee]+) eV")
_S2_RE = re.compile(r"<S\^2> =\s+([-+0-9.Ee]+)")
_SUMMARY_ROW_RE = re.compile(
    r"^\s*(\d+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$"
)
_REFERENCE_ROW_RE = re.compile(r"^\s*0\s+([-+0-9.Ee]+)\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+\(ROHF/UHF Reference state\)")
_SPCSCALE_ORDER_RE = (
    r"spcscale\s*=\s*\[\s*infos%tddft%spc_coco\s*,\s*&?\s*"
    r"infos%tddft%spc_ovov\s*,\s*&?\s*infos%tddft%spc_coov\s*\]"
)


def parse_mrsf_log_state_table(log_text: str) -> dict[int, dict[str, Any]]:
    """Extract compact MRSF root state-character evidence from an OpenQP log.

    Root 0 is the ROHF/UHF reference.  Root N maps to physical S(N-1).
    """

    states: dict[int, dict[str, Any]] = {}
    pending_root: int | None = None
    for line in log_text.splitlines():
        match = _STATE_ENERGY_RE.search(line)
        if match:
            pending_root = int(match.group(1))
            states.setdefault(pending_root, {})["root"] = pending_root
            states[pending_root]["physical_state"] = physical_state_for_mrsf_root(pending_root)
            states[pending_root]["raw_mrsf_root_value_ev"] = float(match.group(2))
            continue
        if pending_root is not None:
            s2_match = _S2_RE.search(line)
            if s2_match:
                states.setdefault(pending_root, {})["s2"] = float(s2_match.group(1))
                pending_root = None
                continue

        summary_match = _SUMMARY_ROW_RE.match(line)
        if summary_match:
            root = int(summary_match.group(1))
            state = states.setdefault(root, {"root": root})
            state.update(
                {
                    "physical_state": physical_state_for_mrsf_root(root),
                    "total_energy_ha": float(summary_match.group(2)),
                    "raw_mrsf_root_value_ev": float(summary_match.group(3)),
                    "physical_excitation_energy_ev": float(summary_match.group(4)),
                    "s2": float(summary_match.group(5)),
                    "transition_dipole_x": float(summary_match.group(6)),
                    "transition_dipole_y": float(summary_match.group(7)),
                    "transition_dipole_z": float(summary_match.group(8)),
                    "transition_dipole_abs": float(summary_match.group(9)),
                    "oscillator_strength": float(summary_match.group(10)),
                    "state_type": "response_root",
                }
            )
            continue

        reference_match = _REFERENCE_ROW_RE.match(line)
        if reference_match:
            states[0] = {
                "root": 0,
                "physical_state": None,
                "total_energy_ha": float(reference_match.group(1)),
                "state_type": "ROHF/UHF Reference state",
            }
    return states


def _read_mrsf_log_state_table(path: Path) -> dict[int, dict[str, Any]]:
    return parse_mrsf_log_state_table(path.read_text(errors="replace"))


def _target_neighbor_gap(states: dict[int, dict[str, Any]], root: int) -> float | None:
    target = states.get(root)
    if not target or "raw_mrsf_root_value_ev" not in target:
        return None
    target_energy = target["raw_mrsf_root_value_ev"]
    gaps = [
        abs(item["raw_mrsf_root_value_ev"] - target_energy)
        for item_root, item in states.items()
        if item_root > 0 and item_root != root and "raw_mrsf_root_value_ev" in item
    ]
    return min(gaps) if gaps else None


def summarize_root_continuity_dir(path: Path | str, root: int, near_degenerate_threshold_ev: float = 1.0e-3) -> dict[str, Any]:
    root_dir = Path(path)
    log_paths = sorted(root_dir.glob("grad/*.log")) + sorted(root_dir.glob("e_*/*.log"))
    entries: list[dict[str, Any]] = []
    s2_values: list[float] = []
    neighbor_gaps: list[float] = []
    trah_logs: list[str] = []
    missing_target_logs: list[str] = []
    for log_path in log_paths:
        text = log_path.read_text(errors="replace")
        if "SCF did not converge. Restarting SCF with the TRAH method." in text or "TRAH / Trust-Region Augmented Hessian Settings" in text:
            trah_logs.append(str(log_path))
        states = parse_mrsf_log_state_table(text)
        target = states.get(root)
        if not target:
            missing_target_logs.append(str(log_path))
            continue
        if "s2" in target:
            s2_values.append(target["s2"])
        gap = _target_neighbor_gap(states, root)
        if gap is not None:
            neighbor_gaps.append(gap)
        entries.append(
            {
                "log": str(log_path),
                "root": root,
                "physical_state": physical_state_for_mrsf_root(root),
                "raw_mrsf_root_value_ev": target.get("raw_mrsf_root_value_ev"),
                "physical_excitation_energy_ev": target.get("physical_excitation_energy_ev"),
                "s2": target.get("s2"),
                "oscillator_strength": target.get("oscillator_strength"),
                "transition_dipole_abs": target.get("transition_dipole_abs"),
                "nearest_neighbor_gap_ev": gap,
            }
        )
    s2_evidence = "present" if s2_values else "unknown"
    s2_max_delta = max(s2_values) - min(s2_values) if s2_values else 0.0
    min_neighbor_gap = min(neighbor_gaps) if neighbor_gaps else None
    near_degenerate = min_neighbor_gap is not None and min_neighbor_gap < near_degenerate_threshold_ev
    if trah_logs:
        evidence_hint = "not-clean: TRAH fallback present; inspect SCF before algebra edits"
    elif s2_evidence == "unknown":
        evidence_hint = "state-character evidence missing; parse/rerun root-continuity metadata before algebra edits"
    elif near_degenerate:
        evidence_hint = "root-continuity risk: target root is near-degenerate with a neighboring response root"
    elif s2_max_delta > 0.5:
        evidence_hint = "root-continuity risk: target <S^2> changes across displaced geometries"
    else:
        evidence_hint = "target <S^2> evidence is present; no large <S^2> change detected"
    return {
        "root_dir": str(root_dir),
        "root": root,
        "physical_state": physical_state_for_mrsf_root(root),
        "log_count": len(log_paths),
        "parsed_target_log_count": len(entries),
        "missing_target_logs": missing_target_logs,
        "s2_evidence": s2_evidence,
        "s2_max_delta": s2_max_delta,
        "min_neighbor_gap_ev": min_neighbor_gap,
        "near_degenerate_target": near_degenerate,
        "trah_log_count": len(trah_logs),
        "trah_logs": trah_logs,
        "evidence_hint": evidence_hint,
        "entries": entries,
    }


def _root_dir_for_component_group(group: dict[str, Any]) -> Path:
    return Path(group["source_csv"]).parent / str(group["molecule"]) / str(group["method"]) / f"root_{int(group['root'])}"


def summarize_root_continuity_targets(
    component_summary_path: Path | str,
    near_degenerate_threshold_ev: float = 1.0e-3,
) -> dict[str, Any]:
    """Summarize root-continuity evidence for target bad groups in a component summary."""

    component_summary = json.loads(Path(component_summary_path).read_text())
    target_groups = component_summary.get("target_bad_groups", [])
    seen: set[tuple[str, str, str, int]] = set()
    cases: list[dict[str, Any]] = []
    missing_cases: list[dict[str, Any]] = []
    for group in target_groups:
        root = int(group["root"])
        molecule = str(group["molecule"])
        method = str(group["method"])
        source_csv = str(group["source_csv"])
        key = (source_csv, molecule, method, root)
        if key in seen:
            continue
        seen.add(key)
        root_dir = _root_dir_for_component_group(group)
        case_meta = {
            "molecule": molecule,
            "method": method,
            "root": root,
            "physical_state": physical_state_for_mrsf_root(root),
            "source_csv": source_csv,
            "expected_root_dir": str(root_dir),
        }
        if not root_dir.exists():
            missing_cases.append(case_meta)
            continue
        case_summary = summarize_root_continuity_dir(root_dir, root, near_degenerate_threshold_ev)
        case_summary.update(case_meta)
        cases.append(case_summary)
    return {
        "source_component_summary": str(component_summary_path),
        "case_count": len(seen),
        "parsed_case_count": len(cases),
        "missing_case_count": len(missing_cases),
        "cases": cases,
        "missing_cases": missing_cases,
    }


def _case_key(item: dict[str, Any]) -> tuple[str, str, str, int]:
    return (
        str(item.get("source_csv", "")),
        str(item.get("molecule", "")).lower(),
        str(item.get("method", "")),
        int(item.get("root", 0)),
    )


def _recommended_source_check(group: dict[str, Any]) -> str:
    axes = set(group.get("bad_axes", []))
    if axes == {"z"}:
        return "Prioritize a source-level Z-vector/operator mapping diagnostic for the localized z component; no production algebra edit yet."
    return "Prioritize a source-level SPC/XC density-handoff diagnostic across the listed bad components; no production algebra edit yet."


def _source_file_presence(source_root: Path, relative_paths: Iterable[str]) -> list[dict[str, Any]]:
    files = []
    for relative_path in relative_paths:
        path = source_root / relative_path
        files.append({"path": relative_path, "exists": path.exists()})
    return files


def summarize_source_diagnostic_plan(
    source_targets: dict[str, Any],
    source_root: Path | str = Path("."),
) -> dict[str, Any]:
    """Create a conservative source-inspection plan for the top stable residual.

    The plan intentionally does not edit algebra.  It records the first source
    files and guardrails to inspect before any future production fix claim.
    """

    stable = list(source_targets.get("stable_source_candidates", []))
    stable.sort(key=lambda item: (-float(item.get("max_abs_diff_ha_per_bohr", 0.0)), str(item.get("molecule", "")), int(item.get("root", 0))))
    if not stable:
        return {
            "selected_candidate": None,
            "diagnostic_family": "none",
            "scope_guard": "no production algebra edit; no stable no-TRAH source candidate selected",
            "source_files_to_inspect": [],
            "source_checklist": [],
            "validation_required_before_fix_claim": [
                "finite-difference validation is required before claiming any MRSF gradient source fix"
            ],
        }

    candidate = stable[0]
    axes = set(candidate.get("bad_axes", []))
    if axes == {"z"}:
        diagnostic_family = "localized_z_component"
        checklist = [
            "Inspect MRSF Z-vector/operator mapping for the target root and bad Cartesian component.",
            "Compare gradient assembly terms against existing SPC toggle evidence before changing signs or densities.",
            "Preserve root-continuity classification and avoid NH3-style near-degenerate cases as source-fix evidence.",
        ]
    else:
        diagnostic_family = "multi_component_spc_xc_density_handoff"
        checklist = [
            "Inspect MRSF SPC/XC density handoff across all bad components before changing algebra.",
            "Separate SPC toggle, XC-gradient, and Z-vector hypotheses with source-level diagnostics first.",
            "Preserve root-continuity classification and avoid near-degenerate cases as source-fix evidence.",
        ]

    return {
        "selected_candidate": {
            "molecule": candidate.get("molecule"),
            "method": candidate.get("method"),
            "root": candidate.get("root"),
            "physical_state": candidate.get("physical_state"),
            "max_abs_diff_ha_per_bohr": candidate.get("max_abs_diff_ha_per_bohr"),
            "bad_axes": candidate.get("bad_axes", []),
            "bad_components": candidate.get("bad_components", []),
            "root_dir": candidate.get("root_dir"),
        },
        "diagnostic_family": diagnostic_family,
        "scope_guard": "no production algebra edit; source diagnostic plan only",
        "source_files_to_inspect": _source_file_presence(
            Path(source_root),
            [
                "source/modules/tdhf_mrsf_gradient.F90",
                "source/modules/tdhf_mrsf_z_vector.F90",
            ],
        ),
        "source_checklist": checklist,
        "validation_required_before_fix_claim": [
            "finite-difference validation on the selected stable target residual",
            "root-continuity evidence with real <S^2> and no TRAH for R/R+h/R-h logs",
            "no-fix or pre-change control artifact for the same molecule/root/component",
        ],
    }


def _read_source_text(source_root: Path, relative_path: str) -> str:
    path = source_root / relative_path
    if not path.exists():
        return ""
    return path.read_text(errors="replace")


def _xc_call_has_explicit_xa_xb_handoff(source_text: str) -> bool:
    match = re.search(r"call\s+utddft_xc_gradient\s*\((.*?)\)", source_text, flags=re.IGNORECASE | re.DOTALL)
    if not match:
        return False
    call_text = match.group(1).lower()
    return bool(re.search(r"\bxa\s*=", call_text)) and bool(re.search(r"\bxb\s*=", call_text))


def _find_line_numbers(source_text: str, pattern: str, flags: int = 0) -> list[int]:
    """Return one-based source line numbers matching a diagnostic regex."""

    compiled = re.compile(pattern, flags)
    return [idx for idx, line in enumerate(source_text.splitlines(), start=1) if compiled.search(line)]


def _find_multiline_start_lines(source_text: str, pattern: str, flags: int = 0) -> list[int]:
    """Return one-based start lines for diagnostic regexes spanning continuations."""

    compiled = re.compile(pattern, flags | re.MULTILINE | re.DOTALL)
    return [source_text.count("\n", 0, match.start()) + 1 for match in compiled.finditer(source_text)]


def _line_snippets(source_text: str, lines: Iterable[int], limit: int = 3) -> list[dict[str, Any]]:
    """Return compact one-line snippets for diagnostic source anchors."""

    source_lines = source_text.splitlines()
    snippets = []
    for line in list(lines)[:limit]:
        if 1 <= line <= len(source_lines):
            snippets.append({"line": line, "text": source_lines[line - 1].strip()})
    return snippets


def _bad_components_for_validation(plan: dict[str, Any]) -> list[dict[str, Any]]:
    candidate = plan.get("selected_candidate") or {}
    result = []
    for item in candidate.get("bad_components") or []:
        if not item.get("component"):
            continue
        result.append(
            {
                "component": str(item.get("component")),
                "axis": str(item.get("axis", "")),
                "abs_diff_ha_per_bohr": item.get("abs_diff_ha_per_bohr"),
                "analytic_ha_per_bohr": item.get("analytic_ha_per_bohr"),
                "fd_ha_per_bohr": item.get("fd_ha_per_bohr"),
            }
        )
    return result


def _next_validation_plan(source_plan: dict[str, Any]) -> dict[str, Any]:
    candidate = source_plan.get("selected_candidate") or {}
    bad_components = candidate.get("bad_components") or []
    component_names = [str(item.get("component")) for item in bad_components if item.get("component")]
    return {
        "molecule": candidate.get("molecule"),
        "method": candidate.get("method"),
        "root": candidate.get("root"),
        "physical_state": candidate.get("physical_state"),
        "diagnostic_family": source_plan.get("diagnostic_family"),
        "root_dir": candidate.get("root_dir"),
        "components_to_validate": component_names,
        "bad_components_to_validate": _bad_components_for_validation(source_plan),
        "requires_no_fix_control": True,
        "requires_finite_difference_rerun": True,
        "requires_root_continuity_evidence": True,
        "scope_guard": "diagnostic plan only; no production algebra edit until FD/no-fix/root-continuity validation exists",
    }


def summarize_source_diagnostic_evidence(
    source_plan: dict[str, Any],
    source_root: Path | str = Path("."),
) -> dict[str, Any]:
    """Record static source signals for the selected residual without editing algebra.

    This is a conservative source-diagnostic artifact: it helps choose the next
    RED diagnostic/live FD check, but it is not validation of a production fix.
    """

    root = Path(source_root)
    gradient_text = _read_source_text(root, "source/modules/tdhf_mrsf_gradient.F90")
    z_vector_text = _read_source_text(root, "source/modules/tdhf_mrsf_z_vector.F90")
    source_signals = {
        "gradient_source_present": bool(gradient_text),
        "z_vector_source_present": bool(z_vector_text),
        "ovov_gradient_sign_uses_post_pr153_plus": "df1 = df1 + sgnk*qfspcp2*db2" in gradient_text,
        "ovov_gradient_sign_uses_pre_pr153_minus": "df1 = df1 - sgnk*qfspcp2*db2" in gradient_text,
        "gradient_xc_call_has_explicit_xa_xb_handoff": _xc_call_has_explicit_xa_xb_handoff(gradient_text),
        "z_vector_channel7_overwrites_mrsfcbc_with_td_abxc": bool(
            re.search(r"fmrst1\s*\(\s*1\s*,\s*7\s*,\s*:\s*,\s*:\s*\)\s*=\s*td_abxc", z_vector_text)
        ),
        "z_vector_mrsfcbc_uses_rohf_same_mo": bool(
            re.search(r"call\s+mrsfcbc\s*\([^\n]*mo_a\s*,\s*mo_a", z_vector_text, flags=re.IGNORECASE)
        ),
        "gradient_spcscale_order_present": bool(_find_multiline_start_lines(gradient_text, _SPCSCALE_ORDER_RE, flags=re.IGNORECASE)),
        "z_vector_td_mrsf_den_consumes_seven_channels": bool(
            re.search(r"td_mrsf_den\s*\(\s*1\s*:\s*7\s*,\s*:\s*,\s*:\s*\)\s*=\s*fmrst1", z_vector_text)
        ),
    }
    source_signal_locations = {
        "ovov_gradient_sign_post_pr153_plus_lines": _find_line_numbers(
            gradient_text, r"df1\s*=\s*df1\s*\+\s*sgnk\s*\*\s*qfspcp2\s*\*\s*db2", flags=re.IGNORECASE
        ),
        "ovov_gradient_sign_pre_pr153_minus_lines": _find_line_numbers(
            gradient_text, r"df1\s*=\s*df1\s*-\s*sgnk\s*\*\s*qfspcp2\s*\*\s*db2", flags=re.IGNORECASE
        ),
        "gradient_xc_call_lines": _find_line_numbers(gradient_text, r"call\s+utddft_xc_gradient", flags=re.IGNORECASE),
        "z_vector_channel7_td_abxc_overwrite_lines": _find_line_numbers(
            z_vector_text, r"fmrst1\s*\(\s*1\s*,\s*7\s*,\s*:\s*,\s*:\s*\)\s*=\s*td_abxc", flags=re.IGNORECASE
        ),
        "z_vector_mrsfcbc_rohf_same_mo_lines": _find_line_numbers(
            z_vector_text, r"call\s+mrsfcbc\s*\([^\n]*mo_a\s*,\s*mo_a", flags=re.IGNORECASE
        ),
        "gradient_spcscale_order_lines": _find_multiline_start_lines(
            gradient_text,
            _SPCSCALE_ORDER_RE,
            flags=re.IGNORECASE,
        ),
        "z_vector_td_mrsf_den_handoff_lines": _find_line_numbers(
            z_vector_text,
            r"td_mrsf_den\s*\(\s*1\s*:\s*7\s*,\s*:\s*,\s*:\s*\)\s*=\s*fmrst1",
            flags=re.IGNORECASE,
        ),
    }
    source_signal_snippets = {
        "ovov_gradient_sign_post_pr153_plus": _line_snippets(
            gradient_text, source_signal_locations["ovov_gradient_sign_post_pr153_plus_lines"]
        ),
        "ovov_gradient_sign_pre_pr153_minus": _line_snippets(
            gradient_text, source_signal_locations["ovov_gradient_sign_pre_pr153_minus_lines"]
        ),
        "gradient_xc_call": _line_snippets(gradient_text, source_signal_locations["gradient_xc_call_lines"]),
        "z_vector_channel7_td_abxc_overwrite": _line_snippets(
            z_vector_text, source_signal_locations["z_vector_channel7_td_abxc_overwrite_lines"]
        ),
        "z_vector_mrsfcbc_rohf_same_mo": _line_snippets(
            z_vector_text, source_signal_locations["z_vector_mrsfcbc_rohf_same_mo_lines"]
        ),
        "gradient_spcscale_order": _line_snippets(
            gradient_text, source_signal_locations["gradient_spcscale_order_lines"]
        ),
        "z_vector_td_mrsf_den_handoff": _line_snippets(
            z_vector_text, source_signal_locations["z_vector_td_mrsf_den_handoff_lines"]
        ),
    }
    static_hypotheses: list[str] = []
    if source_signals["z_vector_channel7_overwrites_mrsfcbc_with_td_abxc"]:
        static_hypotheses.append(
            "channel-7 MRSF density provenance should be isolated with a RED source diagnostic and live FD control"
        )
    if not source_signals["gradient_xc_call_has_explicit_xa_xb_handoff"]:
        static_hypotheses.append(
            "MRSF XC-gradient density handoff remains a source-inspection seam; do not infer a fix without spin-resolved FD validation"
        )
    if source_signals["ovov_gradient_sign_uses_post_pr153_plus"]:
        static_hypotheses.append("OV-OV SPC gradient sign is at the merged PR #153 baseline")

    return {
        "selected_candidate": source_plan.get("selected_candidate"),
        "diagnostic_family": source_plan.get("diagnostic_family"),
        "evidence_scope": "static_source_diagnostic_only",
        "scope_guard": "no production algebra edit; static source signals only",
        "source_signals": source_signals,
        "source_signal_locations": source_signal_locations,
        "source_signal_snippets": source_signal_snippets,
        "next_validation_plan": _next_validation_plan(source_plan),
        "static_hypotheses_to_test": static_hypotheses,
        "validation_required_before_fix_claim": [
            "finite-difference validation on the selected stable target residual",
            "root-continuity evidence with real <S^2> and no TRAH for R/R+h/R-h logs",
            "no-fix or pre-change control artifact for the same molecule/root/component",
        ],
    }


def _detect_trah_in_logs(log_paths: Iterable[Path]) -> bool:
    for log_path in log_paths:
        text = log_path.read_text(errors="replace")
        if "SCF did not converge. Restarting SCF with the TRAH method." in text:
            return True
        if "TRAH / Trust-Region Augmented Hessian Settings" in text:
            return True
    return False


def summarize_source_validation_manifest(source_evidence: dict[str, Any]) -> dict[str, Any]:
    """Create a conservative manifest for the next validation controls.

    The manifest records what evidence already exists for the selected stable
    residual and what controls remain missing before any source/algebra fix can
    be claimed.  It does not launch calculations or infer a mechanism.
    """

    plan = source_evidence.get("next_validation_plan") or {}
    root_dir = Path(str(plan.get("root_dir") or ""))
    root_dir_exists = bool(plan.get("root_dir")) and root_dir.exists()
    log_paths = sorted(root_dir.glob("grad/*.log")) + sorted(root_dir.glob("e_*/*.log")) if root_dir_exists else []
    trah_detected = _detect_trah_in_logs(log_paths)
    components = [str(item) for item in plan.get("components_to_validate", [])]
    bad_components_to_validate = list(plan.get("bad_components_to_validate") or [])
    control_dir = root_dir / "validation_controls" if plan.get("root_dir") else Path("validation_controls")
    control_artifact_plan = []
    fd_artifact_exists: list[bool] = []
    no_fix_artifact_exists: list[bool] = []
    for component in components:
        fd_component_csv = control_dir / f"fd_rerun_{component}_components.csv"
        fd_summary_json = control_dir / f"fd_rerun_{component}_summary.json"
        no_fix_control_json = control_dir / f"no_fix_{component}_control.json"
        fd_component_csv_exists = fd_component_csv.exists()
        fd_summary_json_exists = fd_summary_json.exists()
        no_fix_control_json_exists = no_fix_control_json.exists()
        fd_artifact_exists.extend([fd_component_csv_exists, fd_summary_json_exists])
        no_fix_artifact_exists.append(no_fix_control_json_exists)
        control_artifact_plan.append(
            {
                "component": component,
                "fd_component_csv": str(fd_component_csv),
                "fd_component_csv_exists": fd_component_csv_exists,
                "fd_summary_json": str(fd_summary_json),
                "fd_summary_json_exists": fd_summary_json_exists,
                "no_fix_control_json": str(no_fix_control_json),
                "no_fix_control_json_exists": no_fix_control_json_exists,
            }
        )

    def _status_from_exists(exists: list[bool]) -> str:
        if not exists or not any(exists):
            return "missing"
        if all(exists):
            return "present"
        return "partial"

    fd_status = _status_from_exists(fd_artifact_exists)
    no_fix_status = _status_from_exists(no_fix_artifact_exists)
    root_continuity_status = "present" if root_dir_exists and log_paths and not trah_detected else "missing"
    required_controls = [
        {
            "control": "finite_difference_rerun_for_selected_components",
            "components": components,
            "status": fd_status,
            "reason": "static evidence is not a fresh FD validation of a source hypothesis",
        },
        {
            "control": "no_fix_or_pre_change_control_same_case",
            "status": no_fix_status,
            "reason": "needed to separate source-hypothesis effects from existing post-PR #153 residuals",
        },
        {
            "control": "root_continuity_no_trah_evidence",
            "status": root_continuity_status,
            "reason": "existing logs are acceptable only when parseable and TRAH-free",
        },
    ]
    blocking_controls = [control["control"] for control in required_controls if control["status"] != "present"]
    validation_readiness = {
        "ready_for_source_edit": not blocking_controls,
        "status": "ready_for_source_edit" if not blocking_controls else "blocked_missing_or_partial_controls",
        "blocking_controls": blocking_controls,
    }
    return {
        "selected_case": {
            "molecule": plan.get("molecule"),
            "method": plan.get("method"),
            "root": plan.get("root"),
            "physical_state": plan.get("physical_state"),
            "diagnostic_family": plan.get("diagnostic_family"),
            "root_dir": str(root_dir) if plan.get("root_dir") else None,
        },
        "components_to_validate": components,
        "bad_components_to_validate": bad_components_to_validate,
        "existing_evidence": {
            "root_dir_exists": root_dir_exists,
            "log_count": len(log_paths),
            "trah_detected": trah_detected,
            "log_paths": [str(path) for path in log_paths],
        },
        "control_artifact_plan": control_artifact_plan,
        "required_controls": required_controls,
        "validation_readiness": validation_readiness,
        "scope_guard": "diagnostic manifest only; no production algebra edit or fix claim",
    }


def summarize_source_diagnostic_targets(
    component_summary: dict[str, Any],
    root_continuity_summary: dict[str, Any],
) -> dict[str, Any]:
    """Rank target residuals that are clean enough for source-level diagnostics.

    This combines component-level FD residuals with parsed root-continuity evidence.
    Near-degenerate, TRAH, missing, or unparsed cases are deferred rather than
    labeled as algebra/source candidates.
    """

    continuity_by_key = {_case_key(case): case for case in root_continuity_summary.get("cases", [])}
    continuity_by_case = {
        (str(case.get("molecule", "")).lower(), str(case.get("method", "")), int(case.get("root", 0))): case
        for case in root_continuity_summary.get("cases", [])
    }
    missing_by_key = {_case_key(case): case for case in root_continuity_summary.get("missing_cases", [])}
    missing_by_case = {
        (str(case.get("molecule", "")).lower(), str(case.get("method", "")), int(case.get("root", 0))): case
        for case in root_continuity_summary.get("missing_cases", [])
    }
    stable: list[dict[str, Any]] = []
    deferred: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    for group in component_summary.get("target_bad_groups", []):
        key = _case_key(group)
        case_key = (str(group.get("molecule", "")).lower(), str(group.get("method", "")), int(group.get("root", 0)))
        merged = dict(group)
        continuity = continuity_by_key.get(key) or continuity_by_case.get(case_key)
        missing_case = missing_by_key.get(key) or missing_by_case.get(case_key)
        if missing_case is not None:
            merged.update(missing_case)
            merged["diagnostic_status"] = "missing_root_continuity_evidence"
            missing.append(merged)
            continue
        if continuity is None:
            merged["diagnostic_status"] = "missing_root_continuity_evidence"
            missing.append(merged)
            continue
        merged.update(
            {
                "root_dir": continuity.get("root_dir"),
                "s2_evidence": continuity.get("s2_evidence"),
                "s2_max_delta": continuity.get("s2_max_delta"),
                "min_neighbor_gap_ev": continuity.get("min_neighbor_gap_ev"),
                "near_degenerate_target": continuity.get("near_degenerate_target"),
                "trah_log_count": continuity.get("trah_log_count", 0),
                "parsed_target_log_count": continuity.get("parsed_target_log_count"),
            }
        )
        if merged["trah_log_count"]:
            merged["diagnostic_status"] = "defer_trah_or_failed"
            deferred.append(merged)
        elif merged["s2_evidence"] != "present" or merged["near_degenerate_target"]:
            merged["diagnostic_status"] = "defer_root_continuity_risk"
            deferred.append(merged)
        else:
            merged["diagnostic_status"] = "source_diagnostic_candidate"
            merged["recommended_next_check"] = _recommended_source_check(merged)
            stable.append(merged)

    stable.sort(key=lambda item: (-float(item.get("max_abs_diff_ha_per_bohr", 0.0)), item["molecule"], int(item["root"])))
    deferred.sort(key=lambda item: (-float(item.get("max_abs_diff_ha_per_bohr", 0.0)), item["molecule"], int(item["root"])))
    missing.sort(key=lambda item: (item["molecule"], int(item["root"])))
    return {
        "stable_source_candidate_count": len(stable),
        "deferred_root_continuity_risk_count": sum(
            1 for item in deferred if item.get("diagnostic_status") == "defer_root_continuity_risk"
        ),
        "deferred_trah_or_failed_count": sum(1 for item in deferred if item.get("diagnostic_status") == "defer_trah_or_failed"),
        "missing_root_continuity_count": len(missing),
        "stable_source_candidates": stable,
        "deferred_root_continuity_risks": [
            item for item in deferred if item.get("diagnostic_status") == "defer_root_continuity_risk"
        ],
        "deferred_trah_or_failed": [item for item in deferred if item.get("diagnostic_status") == "defer_trah_or_failed"],
        "missing_root_continuity": missing,
    }


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
    parser.add_argument("csv_path", type=Path, nargs="+")
    parser.add_argument("--components", action="store_true", help="Summarize per-component FD diagnostics CSV")
    parser.add_argument("--root-continuity", action="store_true", help="Summarize MRSF root-continuity evidence from a root run directory")
    parser.add_argument(
        "--root-continuity-targets",
        action="store_true",
        help="Resolve target_bad_groups from a component summary and summarize each existing root run directory",
    )
    parser.add_argument(
        "--source-diagnostic-targets",
        action="store_true",
        help="Combine component and root-continuity summaries to rank stable source-diagnostic candidates",
    )
    parser.add_argument(
        "--source-diagnostic-plan",
        action="store_true",
        help="Create a conservative source-inspection plan from ranked source-diagnostic candidates",
    )
    parser.add_argument(
        "--source-diagnostic-evidence",
        action="store_true",
        help="Record static source signals for a planned stable source diagnostic without editing algebra",
    )
    parser.add_argument(
        "--source-validation-manifest",
        action="store_true",
        help="Record existing evidence and missing controls before any source-fix claim",
    )
    parser.add_argument("--source-root", type=Path, default=Path("."), help="Repository root for source diagnostic modes")
    parser.add_argument("--root", type=int, help="Target MRSF response root for --root-continuity")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--output", type=Path, help="Write JSON summary to this path")
    args = parser.parse_args(argv)

    if args.root_continuity:
        if args.root is None:
            parser.error("--root-continuity requires --root")
        if len(args.csv_path) != 1:
            parser.error("--root-continuity accepts exactly one root directory")
        summary = summarize_root_continuity_dir(args.csv_path[0], args.root)
    elif args.root_continuity_targets:
        if len(args.csv_path) != 1:
            parser.error("--root-continuity-targets accepts exactly one component summary JSON")
        summary = summarize_root_continuity_targets(args.csv_path[0])
    elif args.source_diagnostic_targets:
        if len(args.csv_path) != 2:
            parser.error("--source-diagnostic-targets accepts component summary JSON and root-continuity summary JSON")
        component_summary = json.loads(args.csv_path[0].read_text())
        root_continuity_summary = json.loads(args.csv_path[1].read_text())
        summary = summarize_source_diagnostic_targets(component_summary, root_continuity_summary)
    elif args.source_diagnostic_plan:
        if len(args.csv_path) != 1:
            parser.error("--source-diagnostic-plan accepts exactly one source-diagnostic targets JSON")
        source_targets = json.loads(args.csv_path[0].read_text())
        summary = summarize_source_diagnostic_plan(source_targets, source_root=args.source_root)
    elif args.source_diagnostic_evidence:
        if len(args.csv_path) != 1:
            parser.error("--source-diagnostic-evidence accepts exactly one source-diagnostic plan JSON")
        source_plan = json.loads(args.csv_path[0].read_text())
        summary = summarize_source_diagnostic_evidence(source_plan, source_root=args.source_root)
    elif args.source_validation_manifest:
        if len(args.csv_path) != 1:
            parser.error("--source-validation-manifest accepts exactly one source-diagnostic evidence JSON")
        source_evidence = json.loads(args.csv_path[0].read_text())
        summary = summarize_source_validation_manifest(source_evidence)
    elif args.components:
        if len(args.csv_path) == 1:
            summary = summarize_components_csv(args.csv_path[0], args.threshold)
        else:
            summary = summarize_components_csvs(args.csv_path, args.threshold)
    else:
        if len(args.csv_path) != 1:
            parser.error("compact summary mode accepts exactly one CSV path")
        summary = summarize_compact_csv(args.csv_path[0], args.threshold)
    payload = json.dumps(summary, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
