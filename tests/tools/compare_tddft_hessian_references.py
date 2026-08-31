#!/usr/bin/env python3
"""Compare OpenQP TDDFT Hessians with an authoritative GAMESS reference.

The GAMESS input is the punched ``$HESS``/normal-mode ``.dat`` file.  OpenQP
inputs are the analytic and numerical ``.hess.json`` files.  Matrix quantities
are in the units present in the files (normally Hartree/bohr**2); frequencies
are compared in cm**-1.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import numpy as np


_FLOAT_RE = re.compile(r"[+-]?(?:\d+\.\d*|\.\d+)(?:[DE][+-]?\d+)?", re.I)
_FREQUENCY_RE = re.compile(
    r"MODE\s+\d+\s+FREQUENCY\s*=\s*([+-]?(?:\d+\.\d*|\.\d+)(?:[DE][+-]?\d+)?)",
    re.I,
)


def load_gamess_dat(path: str | Path) -> dict[str, Any]:
    """Read the final square matrix in a GAMESS ``$HESS`` section."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    matches = list(re.finditer(r"(?ms)^\s*\$HESS\s*$.*?^\s*\$END\s*$", text))
    if not matches:
        raise ValueError(f"no $HESS section found in {path}")

    rows: dict[int, list[float]] = {}
    for line in matches[-1].group(0).splitlines()[1:]:
        prefix = re.match(r"^\s*(\d+)\s+([12])(?=[ +\-])", line)
        if prefix is None:
            continue
        row = int(prefix.group(1))
        values = [
            float(x.replace("D", "E").replace("d", "e"))
            for x in _FLOAT_RE.findall(line[prefix.end() :])
        ]
        rows.setdefault(row, []).extend(values)
    if not rows:
        raise ValueError(f"empty $HESS section in {path}")
    n = max(rows)
    if set(rows) != set(range(1, n + 1)) or any(len(rows[i]) != n for i in range(1, n + 1)):
        sizes = {i: len(v) for i, v in rows.items()}
        raise ValueError(f"malformed {n}x{n} $HESS section in {path}: row sizes {sizes}")

    frequencies = [
        float(value.replace("D", "E").replace("d", "e"))
        for value in _FREQUENCY_RE.findall(text)
    ]
    return {"hessian": np.asarray([rows[i] for i in range(1, n + 1)]), "frequencies": frequencies}


def load_openqp_json(path: str | Path) -> dict[str, Any]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    hessian = np.asarray(data["hessian"], dtype=float)
    if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
        raise ValueError(f"OpenQP Hessian in {path} is not square: {hessian.shape}")
    natom = len(data.get("atoms", [])) or hessian.shape[0] // 3
    return {
        "hessian": hessian,
        "frequencies": [float(x) for x in data.get("freqs", [])],
        "natom": natom,
    }


def matrix_metrics(left: np.ndarray, right: np.ndarray) -> dict[str, Any]:
    if left.shape != right.shape:
        raise ValueError(f"matrix shape mismatch: {left.shape} != {right.shape}")
    delta = np.asarray(left, dtype=float) - np.asarray(right, dtype=float)
    flat_left, flat_right = left.ravel(), right.ravel()
    if np.std(flat_left) == 0.0 or np.std(flat_right) == 0.0:
        correlation = 1.0 if np.array_equal(flat_left, flat_right) else None
    else:
        correlation = float(np.corrcoef(flat_left, flat_right)[0, 1])
    index = np.unravel_index(int(np.argmax(np.abs(delta))), delta.shape)
    return {
        "max_abs": float(np.max(np.abs(delta))),
        "rms": float(np.sqrt(np.mean(delta * delta))),
        "correlation": correlation,
        "max_abs_index_1based": [int(value + 1) for value in index],
        "signed_difference_at_max": float(delta[index]),
    }


def symmetry_metrics(hessian: np.ndarray) -> dict[str, float]:
    delta = hessian - hessian.T
    return {
        "max_abs": float(np.max(np.abs(delta))),
        "rms": float(np.sqrt(np.mean(delta * delta))),
    }


def translational_metrics(hessian: np.ndarray, natom: int) -> dict[str, Any]:
    if hessian.shape != (3 * natom, 3 * natom):
        raise ValueError(f"Hessian shape {hessian.shape} is inconsistent with {natom} atoms")
    right, left = [], []
    for axis in range(3):
        translation = np.zeros(3 * natom)
        translation[axis::3] = 1.0
        right.append(hessian @ translation)
        left.append(translation @ hessian)
    right_array, left_array = np.asarray(right), np.asarray(left)
    return {
        "right_max_abs": float(np.max(np.abs(right_array))),
        "right_rms": float(np.sqrt(np.mean(right_array * right_array))),
        "left_max_abs": float(np.max(np.abs(left_array))),
        "left_rms": float(np.sqrt(np.mean(left_array * left_array))),
        "right_max_abs_by_axis": [float(np.max(np.abs(row))) for row in right_array],
        "left_max_abs_by_axis": [float(np.max(np.abs(row))) for row in left_array],
    }


def frequency_metrics(openqp: list[float], gamess: list[float]) -> dict[str, Any] | None:
    """Compare physical-mode magnitudes; GAMESS punch output omits mode signs."""
    if not openqp or len(gamess) < len(openqp):
        return None
    oqp = np.sort(np.abs(np.asarray(openqp, dtype=float)))
    gms = np.sort(np.abs(np.asarray(gamess, dtype=float)))[-len(oqp) :]
    result = matrix_metrics(oqp, gms)
    result.update(
        {
            "openqp_abs_sorted_cm-1": oqp.tolist(),
            "gamess_largest_abs_sorted_cm-1": gms.tolist(),
            "note": "Magnitude-only comparison; GAMESS punched frequencies do not encode imaginary-mode signs.",
        }
    )
    return result


def compare_case(
    name: str,
    gamess_path: str | Path,
    analytic_path: str | Path,
    numerical_path: str | Path,
) -> dict[str, Any]:
    gamess = load_gamess_dat(gamess_path)
    analytic = load_openqp_json(analytic_path)
    numerical = load_openqp_json(numerical_path)
    natom = analytic["natom"]
    matrices = {
        "gamess": gamess["hessian"],
        "openqp_analytic": analytic["hessian"],
        "openqp_numerical": numerical["hessian"],
    }
    if numerical["natom"] != natom:
        raise ValueError("analytic and numerical OpenQP files have different atom counts")
    if any(matrix.shape != analytic["hessian"].shape for matrix in matrices.values()):
        raise ValueError(f"incompatible Hessian dimensions in case {name}")
    return {
        "name": name,
        "files": {
            "gamess": str(Path(gamess_path).resolve()),
            "openqp_analytic": str(Path(analytic_path).resolve()),
            "openqp_numerical": str(Path(numerical_path).resolve()),
        },
        "dimension": int(3 * natom),
        "comparisons": {
            "analytic_vs_numerical": matrix_metrics(matrices["openqp_analytic"], matrices["openqp_numerical"]),
            "analytic_vs_gamess": matrix_metrics(matrices["openqp_analytic"], matrices["gamess"]),
            "numerical_vs_gamess": matrix_metrics(matrices["openqp_numerical"], matrices["gamess"]),
        },
        "symmetry": {key: symmetry_metrics(value) for key, value in matrices.items()},
        "translation": {key: translational_metrics(value, natom) for key, value in matrices.items()},
        "frequencies": {
            "analytic_vs_gamess": frequency_metrics(analytic["frequencies"], gamess["frequencies"]),
            "numerical_vs_gamess": frequency_metrics(numerical["frequencies"], gamess["frequencies"]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        nargs=4,
        metavar=("NAME", "GAMESS_DAT", "OPENQP_ANALYTIC_JSON", "OPENQP_NUMERICAL_JSON"),
        required=True,
        help="comparison case; repeat for multiple functionals",
    )
    parser.add_argument("--output-json", type=Path, help="also write the JSON report to this path")
    parser.add_argument(
        "--fail-max-analytic-numerical",
        type=float,
        help="exit with status 1 if any analytic/numerical maximum error exceeds this value",
    )
    args = parser.parse_args()
    report = {"cases": [compare_case(*case) for case in args.case]}
    rendered = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    print(rendered, end="")
    if args.output_json:
        args.output_json.write_text(rendered, encoding="utf-8")
    if args.fail_max_analytic_numerical is not None:
        maximum = max(case["comparisons"]["analytic_vs_numerical"]["max_abs"] for case in report["cases"])
        if maximum > args.fail_max_analytic_numerical:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
