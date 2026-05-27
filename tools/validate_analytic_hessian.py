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
    args = parser.parse_args(argv)

    summary = compare_hessians(_load_matrix(args.analytic), _load_matrix(args.reference), top_n=args.top)
    payload = summary_to_json(summary)
    if args.output:
        args.output.write_text(payload + "\n")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
