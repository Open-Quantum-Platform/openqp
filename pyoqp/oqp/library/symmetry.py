"""Backend-free symmetry utilities for one-electron diagnostics.

This module currently provides metadata-only diagnostics used by the symmetry
planning gates. It does not change SCF/integral/response execution behavior.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping

import numpy as np


@dataclass(frozen=True)
class OneElectronBlockLeakSummary:
    """Normalized summary payload for one-electron block leakage checks."""

    max_off_block_abs: float
    max_off_block_indices: list[tuple[int, int]]
    off_block_element_count: int
    within_tolerance: bool
    orthogonality_ok: bool
    orthogonality_max_deviation: float
    status: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "max_off_block_abs": float(self.max_off_block_abs),
            "max_off_block_indices": [list(pair) for pair in self.max_off_block_indices],
            "off_block_element_count": int(self.off_block_element_count),
            "within_tolerance": bool(self.within_tolerance),
            "transform_orthogonality": {
                "ok": bool(self.orthogonality_ok),
                "max_deviation": float(self.orthogonality_max_deviation),
            },
            "status": self.status,
            "use_integral_symmetry": False,
            "use_response_symmetry": False,
        }


def _to_float_array(values: Any) -> np.ndarray:
    return np.asarray(values, dtype=float)


def _as_list(values: Any) -> list[Any]:
    if values is None:
        return []
    if isinstance(values, list):
        return values
    return list(values)


def one_electron_block_diagnostics(
    matrix: Any,
    symmetry_adapted_transform: Any,
    basis_labels: Iterable[str],
    tolerance: float = 1.0e-6,
) -> dict[str, Any]:
    """Compute off-block leakage diagnostics in a symmetry-adapted basis.

    Parameters
    ----------
    matrix:
        One-electron matrix (overlap or core-Hamiltonian like).
    symmetry_adapted_transform:
        AO-space transformation matrix to symmetry-adapted AO basis.
    basis_labels:
        Per-AO symmetry labels in the transformed basis.
    tolerance:
        Off-block absolute-coupling tolerance.
    """

    if tolerance <= 0:
        raise ValueError("tolerance must be positive")

    m = _to_float_array(matrix)
    if m.ndim != 2 or m.shape[0] != m.shape[1]:
        raise ValueError("matrix must be a square 2D array")

    label_list = [str(lbl) for lbl in _as_list(basis_labels)]
    if len(label_list) != m.shape[0]:
        raise ValueError("basis_labels length must match matrix dimension")

    u = _to_float_array(symmetry_adapted_transform)
    if u.shape != (m.shape[0], m.shape[0]):
        raise ValueError("symmetry_adapted_transform must be square with matching dimension")

    # Transform into symmetry-adapted AO basis.
    transformed = u.T @ m @ u

    # Identify off-block entries by label mismatch.
    labels = np.array(label_list)
    mismatch = labels[:, None] != labels[None, :]
    off_block_values = np.abs(transformed[mismatch])
    if off_block_values.size == 0:
        max_off_block = 0.0
        max_indices: list[tuple[int, int]] = []
        off_count = 0
    else:
        max_off_block = float(np.max(off_block_values))
        flat_index = int(np.argmax(off_block_values))
        # Recover one flattened index from compressed mismatch view to AO pair indices.
        mismatch_positions = np.argwhere(mismatch)
        pair = mismatch_positions[flat_index]
        max_indices = [(int(pair[0]), int(pair[1]))]
        off_count = int(mismatch.sum())

    orthogonality = u.T @ u
    orthogonality_max_deviation = float(np.max(np.abs(orthogonality - np.eye(u.shape[0]))))
    orthogonality_ok = orthogonality_max_deviation <= tolerance

    payload = OneElectronBlockLeakSummary(
        max_off_block_abs=max_off_block,
        max_off_block_indices=max_indices,
        off_block_element_count=off_count,
        within_tolerance=(max_off_block <= tolerance),
        orthogonality_ok=orthogonality_ok,
        orthogonality_max_deviation=orthogonality_max_deviation,
        status="diagnostic_only_no_reductions",
    ).as_dict()
    payload["transform_shape"] = tuple(int(x) for x in u.shape)

    return payload


def update_one_electron_block_diagnostics(
    symmetry_metadata: Mapping[str, Any] | None,
    matrices: Mapping[str, Any],
    symmetry_adapted_transform: Any,
    basis_labels: Iterable[str],
    tolerance: float = 1.0e-6,
) -> dict[str, Any]:
    """Attach named one-electron diagnostics under
    ``symmetry_metadata['one_electron_block_diagnostics']``.

    This is metadata-only and intentionally does not enable any symmetry
    acceleration behavior.
    """

    if symmetry_metadata is None:
        metadata: dict[str, Any] = {}
    else:
        metadata = dict(symmetry_metadata)

    diagnostics: dict[str, Any] = {}
    for key, matrix in matrices.items():
        diagnostics[str(key)] = one_electron_block_diagnostics(
            matrix, symmetry_adapted_transform, basis_labels,
            tolerance=tolerance,
        )

    metadata["one_electron_block_diagnostics"] = diagnostics
    metadata["use_integral_symmetry"] = False
    metadata["use_response_symmetry"] = False
    metadata["one_electron_block_diagnostics_status"] = "diagnostic_only_no_reductions"

    return metadata


def _cartesian_shell_size(l: int) -> int:
    """Return Cartesian basis-function count for a shell with angular momentum l."""

    if l < 0:
        raise ValueError("angular momentum must be non-negative")
    if l == 0:
        return 1
    if l == 1:
        return 3
    if l == 2:
        return 6
    raise ValueError("only s/p/d Cartesian shells are supported in this metadata-only scaffold")
