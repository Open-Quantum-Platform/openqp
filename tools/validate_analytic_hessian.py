#!/usr/bin/env python3
"""Dependency-light helpers for analytic Hessian validation summaries.

This module intentionally does not launch OpenQP jobs yet.  Early analytic
Hessian work needs a stable reporting contract first: compare an analytic
Hessian against a finite-difference/reference Hessian, report compact max/RMS
metrics, surface the largest offending components, and never hide asymmetry by
printing only a symmetrized matrix.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

Matrix = list[list[float]]


def _as_float_matrix(name: str, value: Any) -> Matrix:
    if hasattr(value, "tolist"):
        value = value.tolist()
    if not isinstance(value, (list, tuple)) or not value:
        raise ValueError(f"{name} must be a nonempty 2-D matrix")

    matrix: Matrix = []
    width = None
    for row in value:
        if not isinstance(row, (list, tuple)):
            raise ValueError(f"{name} must be a 2-D matrix")
        if width is None:
            width = len(row)
            if width == 0:
                raise ValueError(f"{name} must be a nonempty 2-D matrix")
        elif len(row) != width:
            raise ValueError(f"{name} rows must all have the same length")
        converted_row = []
        for col_index, item in enumerate(row):
            number = float(item)
            if not math.isfinite(number):
                raise ValueError(f"{name}[{len(matrix)}][{col_index}] must be finite")
            converted_row.append(number)
        matrix.append(converted_row)
    return matrix


def _shape(matrix: Matrix) -> tuple[int, int]:
    return len(matrix), len(matrix[0])


def compare_hessians(analytic: Any, reference: Any, top_n: int = 10) -> dict[str, Any]:
    """Return compact analytic-vs-reference Hessian error metrics.

    Parameters
    ----------
    analytic, reference
        Square or rectangular Hessian-like matrices with identical shape.
    top_n
        Number of largest absolute-difference components to include.  This is
        capped at the total component count and can be zero.
    """

    analytic_matrix = _as_float_matrix("analytic", analytic)
    reference_matrix = _as_float_matrix("reference", reference)
    analytic_shape = _shape(analytic_matrix)
    reference_shape = _shape(reference_matrix)
    if analytic_shape != reference_shape:
        raise ValueError(
            "analytic and reference Hessians must have the same shape; "
            f"got {analytic_shape} and {reference_shape}"
        )

    nrow, ncol = analytic_shape
    components = []
    sum_sq = 0.0
    max_abs_diff = 0.0
    for row in range(nrow):
        for col in range(ncol):
            delta = analytic_matrix[row][col] - reference_matrix[row][col]
            abs_delta = abs(delta)
            sum_sq += delta * delta
            max_abs_diff = max(max_abs_diff, abs_delta)
            components.append((row, col, delta, abs_delta))
    rms_diff = math.sqrt(sum_sq / len(components)) if components else 0.0

    if nrow == ncol and components:
        max_asymmetry = max(
            abs(analytic_matrix[row][col] - analytic_matrix[col][row])
            for row in range(nrow)
            for col in range(ncol)
        )
    else:
        max_asymmetry = 0.0

    ranked_components = sorted(components, key=lambda item: (-item[3], item[0], item[1]))
    largest_components = []
    for row, col, delta, abs_delta in ranked_components[: max(0, min(int(top_n), len(components)))]:
        largest_components.append(
            {
                "row": int(row),
                "col": int(col),
                "computed_value": analytic_matrix[row][col],
                "reference_value": reference_matrix[row][col],
                "diff": delta,
                "abs_diff": abs_delta,
            }
        )

    return {
        "schema_version": "analytic_hessian_validation.v1",
        "report_type": "metric_comparison",
        "matrix_payload": "omitted",
        "shape": [int(nrow), int(ncol)],
        "max_abs_diff": max_abs_diff,
        "rms_diff": rms_diff,
        "max_asymmetry": max_asymmetry,
        "largest_components": largest_components,
    }


def _positive_finite_context_scalar(name: str, value: float) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return number


def _format_float(value: float) -> str:
    return f"{value:.12g}"


def _failure_summary_lines(failed_metric_details: list[dict[str, Any]]) -> list[str]:
    lines = []
    for detail in failed_metric_details:
        line = (
            f"{detail['metric']} observed {_format_float(detail['observed'])} "
            f"exceeds tolerance {_format_float(detail['tolerance'])} "
            f"by {_format_float(detail['excess'])}"
        )
        worst_component = detail.get("worst_component")
        if worst_component is not None:
            line += f" at component [{worst_component['row']},{worst_component['col']}]"
        lines.append(line)
    return lines


def build_validation_summary(
    analytic: Any,
    reference: Any,
    *,
    method: str,
    td_type: str,
    state: int,
    molecule: str,
    basis: str,
    displacement: float,
    max_tolerance: float,
    rms_tolerance: float,
    asymmetry_tolerance: float | None = None,
    top_n: int = 10,
) -> dict[str, Any]:
    """Compare Hessians and attach validation context plus pass/fail status.

    This stays deliberately dependency-light and file-format agnostic: callers can
    supply analytic Hessians from a native kernel and reference Hessians from
    central finite differences, while this helper records the scientific context
    needed to interpret the comparison without dumping full matrices.
    """

    summary = compare_hessians(analytic, reference, top_n=top_n)
    displacement_value = _positive_finite_context_scalar("displacement", displacement)
    max_tolerance_value = _positive_finite_context_scalar("max_tolerance", max_tolerance)
    rms_tolerance_value = _positive_finite_context_scalar("rms_tolerance", rms_tolerance)
    tolerances = {
        "max_abs_diff": max_tolerance_value,
        "rms_diff": rms_tolerance_value,
    }
    if asymmetry_tolerance is not None:
        tolerances["max_asymmetry"] = _positive_finite_context_scalar(
            "asymmetry_tolerance", asymmetry_tolerance
        )

    failed_metrics = []
    failed_metric_details = []
    worst_component = summary["largest_components"][0] if summary["largest_components"] else None
    for metric in ("max_abs_diff", "rms_diff", "max_asymmetry"):
        if metric not in tolerances:
            continue
        observed = float(summary[metric])
        tolerance = tolerances[metric]
        if observed > tolerance:
            failed_metrics.append(metric)
            detail = {
                "metric": metric,
                "observed": observed,
                "tolerance": tolerance,
                "excess": observed - tolerance,
            }
            if metric in {"max_abs_diff", "rms_diff"} and worst_component is not None:
                detail["worst_component"] = worst_component
            failed_metric_details.append(detail)

    return {
        "method": method,
        "td_type": td_type,
        "state": int(state),
        "molecule": molecule,
        "basis": basis,
        "displacement": displacement_value,
        "tolerances": tolerances,
        "passed": not failed_metrics,
        "failed_metrics": failed_metrics,
        "failed_metric_details": failed_metric_details,
        "failure_summary": _failure_summary_lines(failed_metric_details),
        **summary,
        "report_type": "contextual_validation",
    }


def summary_to_json(summary: dict[str, Any]) -> str:
    """Serialize a compact Hessian validation summary deterministically."""

    return json.dumps(summary, indent=2, sort_keys=True)


def _load_matrix(path: Path) -> Matrix:
    if path.suffix == ".npy":
        try:
            import numpy as np  # type: ignore
        except ModuleNotFoundError as exc:
            raise SystemExit("reading .npy Hessian matrices requires NumPy") from exc
        return _as_float_matrix(str(path), np.load(path))

    if path.suffix == ".json" or path.name.endswith(".hess.json"):
        data = json.loads(path.read_text())
        if "hessian" not in data:
            raise ValueError(f"{path} must contain a hessian field")
        return _as_float_matrix(f"{path}:hessian", data["hessian"])

    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped:
            rows.append([float(item) for item in stripped.split()])
    return _as_float_matrix(str(path), rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare analytic and reference/finite-difference Hessian matrices."
    )
    parser.add_argument("analytic", type=Path, help="Analytic Hessian matrix (.npy or text)")
    parser.add_argument("reference", type=Path, help="Reference Hessian matrix (.npy or text)")
    parser.add_argument("--top", type=int, default=10, help="Largest components to report")
    parser.add_argument("--output", type=Path, help="Optional JSON summary output path")
    parser.add_argument("--method", help="Method label for contextual validation summaries")
    parser.add_argument("--td-type", help="TDHF type label, or 'none' for ground state")
    parser.add_argument("--state", type=int, help="State/root index used for the comparison")
    parser.add_argument("--molecule", help="Molecule/system label")
    parser.add_argument("--basis", help="Basis-set label")
    parser.add_argument("--displacement", type=float, help="Finite-difference displacement used for the reference Hessian")
    parser.add_argument("--max-tolerance", type=float, help="Pass/fail threshold for max_abs_diff")
    parser.add_argument("--rms-tolerance", type=float, help="Pass/fail threshold for rms_diff")
    parser.add_argument(
        "--asymmetry-tolerance",
        type=float,
        help="Optional pass/fail threshold for analytic Hessian max_asymmetry",
    )
    args = parser.parse_args(argv)

    analytic = _load_matrix(args.analytic)
    reference = _load_matrix(args.reference)
    contextual_fields = [
        args.method,
        args.td_type,
        args.state,
        args.molecule,
        args.basis,
        args.displacement,
        args.max_tolerance,
        args.rms_tolerance,
        args.asymmetry_tolerance,
    ]
    if any(value is not None for value in contextual_fields):
        missing = [
            name
            for name, value in [
                ("--method", args.method),
                ("--td-type", args.td_type),
                ("--state", args.state),
                ("--molecule", args.molecule),
                ("--basis", args.basis),
                ("--displacement", args.displacement),
                ("--max-tolerance", args.max_tolerance),
                ("--rms-tolerance", args.rms_tolerance),
            ]
            if value is None
        ]
        if missing:
            parser.error("contextual validation summaries require " + ", ".join(missing))
        try:
            summary = build_validation_summary(
                analytic,
                reference,
                method=args.method,
                td_type=args.td_type,
                state=args.state,
                molecule=args.molecule,
                basis=args.basis,
                displacement=args.displacement,
                max_tolerance=args.max_tolerance,
                rms_tolerance=args.rms_tolerance,
                asymmetry_tolerance=args.asymmetry_tolerance,
                top_n=args.top,
            )
        except ValueError as exc:
            parser.error(str(exc))
    else:
        summary = compare_hessians(analytic, reference, top_n=args.top)
    summary["matrix_sources"] = {
        "analytic": str(args.analytic),
        "reference": str(args.reference),
    }
    summary["matrix_source_sha256"] = {
        "analytic": hashlib.sha256(args.analytic.read_bytes()).hexdigest(),
        "reference": hashlib.sha256(args.reference.read_bytes()).hexdigest(),
    }
    payload = summary_to_json(summary)
    if args.output:
        args.output.write_text(payload + "\n")
    else:
        print(payload)
    if "passed" in summary and not summary["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
