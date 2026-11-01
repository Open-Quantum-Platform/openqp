#!/usr/bin/env python3
"""Compare GAMESS TDDFT Hessian punch components with an OpenQP dump.

GAMESS writes ``$TDDXR``, ``$TDD1G``, and ``$TDD2G`` one coordinate
at a time: the vector following ``COORD K`` is column K of the Cartesian
matrix.  ``$TDH2XC`` instead contains lower-triangular atom-pair 3x3 blocks.
OpenQP's ``.components`` file prints ordinary matrix rows after a label.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import numpy as np


_FLOAT_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DE][+-]?\d+)?", re.I)
_COLUMN_LABEL = {"TDDXR": "DXR", "TDD1G": "D1G", "TDD2G": "D2G"}
_OPENQP_LABEL = {
    "TDDXR": "rowsxc",
    "TDD1G": "rows_one",
    "TDD2G": "rows_two",
    "TDH2XC": "hxc",
}


def _floats(text: str) -> list[float]:
    return [float(value.replace("D", "E").replace("d", "e")) for value in _FLOAT_RE.findall(text)]


def _section(text: str, tag: str) -> tuple[str, str]:
    match = re.search(rf"(?ms)^\s*\${re.escape(tag)}\b([^\n]*)\n(.*?)^\s*\$END\s*$", text)
    if match is None:
        raise ValueError(f"no ${tag} section found")
    return match.group(1), match.group(2)


def load_gamess_coordinate_component(path: str | Path, tag: str) -> np.ndarray:
    """Load a TDD component, interpreting each ``COORD K`` vector as a column."""
    if tag not in _COLUMN_LABEL:
        raise ValueError(f"unsupported coordinate component {tag!r}")
    header, body = _section(Path(path).read_text(encoding="utf-8", errors="replace"), tag)
    size_match = re.search(r"NXYZ\s*=\s*(\d+)", header, re.I)
    if size_match is None:
        raise ValueError(f"${tag} header has no NXYZ")
    nxyz = int(size_match.group(1))
    label = _COLUMN_LABEL[tag]
    pieces = re.split(rf"(?m)^\s*{label}\s+COORD\s+(\d+)\s*$", body)
    columns: dict[int, list[float]] = {}
    for index in range(1, len(pieces), 2):
        columns[int(pieces[index])] = _floats(pieces[index + 1])
    if set(columns) != set(range(1, nxyz + 1)):
        raise ValueError(f"${tag} has incomplete coordinate labels: {sorted(columns)}")
    if any(len(columns[k]) != nxyz for k in columns):
        raise ValueError(f"${tag} coordinate vectors are not length {nxyz}")
    return np.asarray([columns[k] for k in range(1, nxyz + 1)], dtype=float).T


def load_gamess_tdh2xc(path: str | Path) -> np.ndarray:
    """Expand GAMESS lower-triangular atom-pair 3x3 blocks into a square matrix."""
    header, body = _section(Path(path).read_text(encoding="utf-8", errors="replace"), "TDH2XC")
    nat_match = re.search(r"NAT\s*,\s*NAT2\s*=\s*(\d+)\s+(\d+)", header, re.I)
    if nat_match is None:
        raise ValueError("$TDH2XC header has no NAT,NAT2")
    nat, nat2 = map(int, nat_match.groups())
    if nat2 != nat * (nat + 1) // 2:
        raise ValueError("$TDH2XC NAT2 is inconsistent with NAT")
    values = _floats(body)
    if len(values) != 9 * nat2:
        raise ValueError(f"$TDH2XC has {len(values)} values; expected {9 * nat2}")
    matrix = np.zeros((3 * nat, 3 * nat), dtype=float)
    offset = 0
    for atom_i in range(nat):
        for atom_j in range(atom_i + 1):
            block = np.asarray(values[offset : offset + 9]).reshape(3, 3)
            offset += 9
            si, sj = slice(3 * atom_i, 3 * atom_i + 3), slice(3 * atom_j, 3 * atom_j + 3)
            matrix[si, sj] = block
            if atom_i != atom_j:
                matrix[sj, si] = block.T
    return matrix


def load_openqp_component(path: str | Path, label: str, size: int) -> np.ndarray:
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    matches = [index for index, line in enumerate(lines) if line.strip() == label]
    if not matches:
        raise ValueError(f"OpenQP component label {label!r} not found")
    start = matches[-1] + 1
    rows = [_floats(lines[start + row]) for row in range(size)]
    if any(len(row) != size for row in rows):
        raise ValueError(f"OpenQP component {label!r} is not {size}x{size}")
    return np.asarray(rows, dtype=float)


def matrix_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict[str, Any]:
    if reference.shape != candidate.shape:
        raise ValueError(f"component shapes differ: {reference.shape} != {candidate.shape}")
    delta = candidate - reference
    denominator = float(np.vdot(reference, reference))
    scale = float(np.vdot(reference, candidate) / denominator) if denominator else None
    fitted = candidate if scale is None else scale * reference
    fit_delta = fitted - candidate
    if np.std(reference) == 0.0 or np.std(candidate) == 0.0:
        correlation = 1.0 if np.array_equal(reference, candidate) else None
    else:
        correlation = float(np.corrcoef(reference.ravel(), candidate.ravel())[0, 1])
    index = np.unravel_index(int(np.argmax(np.abs(delta))), delta.shape)
    return {
        "max_abs": float(np.max(np.abs(delta))),
        "rms": float(np.sqrt(np.mean(delta * delta))),
        "correlation": correlation,
        "reference_norm": float(np.linalg.norm(reference)),
        "candidate_norm": float(np.linalg.norm(candidate)),
        "least_squares_scale": scale,
        "scaled_max_abs": float(np.max(np.abs(fit_delta))),
        "scaled_rms": float(np.sqrt(np.mean(fit_delta * fit_delta))),
        "max_abs_index_1based": [int(value + 1) for value in index],
        "candidate_minus_reference_at_max": float(delta[index]),
    }


def compare_components(gamess_path: str | Path, openqp_path: str | Path) -> dict[str, Any]:
    gamess = {
        tag: load_gamess_coordinate_component(gamess_path, tag)
        for tag in _COLUMN_LABEL
    }
    gamess["TDH2XC"] = load_gamess_tdh2xc(gamess_path)
    size = gamess["TDDXR"].shape[0]
    if any(matrix.shape != (size, size) for matrix in gamess.values()):
        raise ValueError("GAMESS component dimensions disagree")
    openqp = {
        tag: load_openqp_component(openqp_path, label, size)
        for tag, label in _OPENQP_LABEL.items()
    }
    direct = {tag: matrix_metrics(gamess[tag], openqp[tag]) for tag in gamess}
    transposed = {tag: matrix_metrics(gamess[tag].T, openqp[tag]) for tag in gamess}
    gamess_rows = gamess["TDD1G"] + gamess["TDD2G"] + gamess["TDDXR"]
    openqp_rows = openqp["TDD1G"] + openqp["TDD2G"] + openqp["TDDXR"]
    return {
        "files": {"gamess": str(Path(gamess_path).resolve()), "openqp": str(Path(openqp_path).resolve())},
        "dimension": size,
        "components": direct,
        "transposed_control": transposed,
        "combined_rows": matrix_metrics(gamess_rows, openqp_rows),
        "combined_rows_and_hxc": matrix_metrics(
            gamess_rows + gamess["TDH2XC"], openqp_rows + openqp["TDH2XC"]
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gamess", help="GAMESS punch/data file")
    parser.add_argument("openqp", help="OpenQP .components file")
    parser.add_argument("--output", help="optional JSON output file")
    args = parser.parse_args()
    result = compare_components(args.gamess, args.openqp)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
