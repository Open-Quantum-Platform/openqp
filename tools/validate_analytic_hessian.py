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
import json
from pathlib import Path
from typing import Any

import numpy as np


def _as_float_matrix(name: str, value: Any) -> np.ndarray:
    matrix = np.asarray(value, dtype=float)
    if matrix.ndim != 2:
        raise ValueError(f"{name} must be a 2-D matrix, got shape {matrix.shape}")
    return matrix


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
    if analytic_matrix.shape != reference_matrix.shape:
        raise ValueError(
            "analytic and reference Hessians must have the same shape; "
            f"got {analytic_matrix.shape} and {reference_matrix.shape}"
        )

    diff = analytic_matrix - reference_matrix
    abs_diff = np.abs(diff)
    max_abs_diff = float(np.max(abs_diff)) if diff.size else 0.0
    rms_diff = float(np.sqrt(np.mean(diff * diff))) if diff.size else 0.0

    if analytic_matrix.shape[0] == analytic_matrix.shape[1] and analytic_matrix.size:
        max_asymmetry = float(np.max(np.abs(analytic_matrix - analytic_matrix.T)))
    else:
        max_asymmetry = 0.0

    ranked_indices = sorted(
        np.ndindex(diff.shape),
        key=lambda rc: (-float(abs_diff[rc]), int(rc[0]), int(rc[1])),
    )
    largest_components = []
    for row, col in ranked_indices[: max(0, min(int(top_n), diff.size))]:
        largest_components.append(
            {
                "row": int(row),
                "col": int(col),
                "computed_value": float(analytic_matrix[row, col]),
                "reference_value": float(reference_matrix[row, col]),
                "diff": float(diff[row, col]),
                "abs_diff": float(abs_diff[row, col]),
            }
        )

    return {
        "shape": [int(analytic_matrix.shape[0]), int(analytic_matrix.shape[1])],
        "max_abs_diff": max_abs_diff,
        "rms_diff": rms_diff,
        "max_asymmetry": max_asymmetry,
        "largest_components": largest_components,
    }


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
    top_n: int = 10,
) -> dict[str, Any]:
    """Compare Hessians and attach validation context plus pass/fail status.

    This stays deliberately dependency-light and file-format agnostic: callers can
    supply analytic Hessians from a native kernel and reference Hessians from
    central finite differences, while this helper records the scientific context
    needed to interpret the comparison without dumping full matrices.
    """

    summary = compare_hessians(analytic, reference, top_n=top_n)
    tolerances = {
        "max_abs_diff": float(max_tolerance),
        "rms_diff": float(rms_tolerance),
    }
    failed_metrics = []
    if summary["max_abs_diff"] > tolerances["max_abs_diff"]:
        failed_metrics.append("max_abs_diff")
    if summary["rms_diff"] > tolerances["rms_diff"]:
        failed_metrics.append("rms_diff")

    return {
        "method": method,
        "td_type": td_type,
        "state": int(state),
        "molecule": molecule,
        "basis": basis,
        "displacement": float(displacement),
        "tolerances": tolerances,
        "passed": not failed_metrics,
        "failed_metrics": failed_metrics,
        **summary,
    }


def summary_to_json(summary: dict[str, Any]) -> str:
    """Serialize a compact Hessian validation summary deterministically."""

    return json.dumps(summary, indent=2, sort_keys=True)


def _load_matrix(path: Path) -> np.ndarray:
    if path.suffix == ".npy":
        return np.load(path)
    return np.loadtxt(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare analytic and reference/finite-difference Hessian matrices."
    )
    parser.add_argument("analytic", type=Path, help="Analytic Hessian matrix (.npy or text)")
    parser.add_argument("reference", type=Path, help="Reference Hessian matrix (.npy or text)")
    parser.add_argument("--top", type=int, default=10, help="Largest components to report")
    parser.add_argument("--output", type=Path, help="Optional JSON summary output path")
    parser.add_argument("--method", help="Method label for contextual validation summaries")
    parser.add_argument("--td-type", default="none", help="TDHF type label, or 'none' for ground state")
    parser.add_argument("--state", type=int, help="State/root index used for the comparison")
    parser.add_argument("--molecule", help="Molecule/system label")
    parser.add_argument("--basis", help="Basis-set label")
    parser.add_argument("--displacement", type=float, help="Finite-difference displacement used for the reference Hessian")
    parser.add_argument("--max-tolerance", type=float, help="Pass/fail threshold for max_abs_diff")
    parser.add_argument("--rms-tolerance", type=float, help="Pass/fail threshold for rms_diff")
    args = parser.parse_args(argv)

    analytic = _load_matrix(args.analytic)
    reference = _load_matrix(args.reference)
    contextual_fields = [
        args.method,
        args.state,
        args.molecule,
        args.basis,
        args.displacement,
        args.max_tolerance,
        args.rms_tolerance,
    ]
    if any(value is not None for value in contextual_fields):
        missing = [
            name
            for name, value in [
                ("--method", args.method),
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
            top_n=args.top,
        )
    else:
        summary = compare_hessians(analytic, reference, top_n=args.top)
    payload = summary_to_json(summary)
    if args.output:
        args.output.write_text(payload + "\n")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
