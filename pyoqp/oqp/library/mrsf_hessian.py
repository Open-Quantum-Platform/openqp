"""MRSF Hessian oracle helpers.

These helpers are backend-free diagnostics for Gate 3B. They do not enable an
MRSF Hessian by themselves; they provide the root-character checks that the
finite-difference oracle must pass before a displaced gradient/energy can be
used in a Hessian assembly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class MRSFRootTrackingError(RuntimeError):
    """Raised when a displaced MRSF calculation cannot be assigned safely."""


@dataclass(frozen=True)
class MRSFRootTrackingDiagnostic:
    requested_state: int
    matched_state: int
    overlap: float
    overlap_gap: float
    energy_gap: float
    root_flip: bool
    ok: bool


def physical_root_label(state: int) -> str:
    """Return the OpenQP/MRSF physical root label for a state index."""

    if state <= 0:
        return "high-spin reference"
    return f"physical S{state - 1}"


def _as_2d_vectors(vectors: np.ndarray | list[list[float]], *, name: str) -> np.ndarray:
    array = np.asarray(vectors, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"{name} must be a 2D array of state vectors")
    norms = np.linalg.norm(array, axis=1)
    if np.any(norms <= 0.0):
        raise ValueError(f"{name} contains a zero-norm state vector")
    return array / norms[:, None]


def assess_root_tracking(
    *,
    requested_state: int,
    reference_vectors: np.ndarray | list[list[float]],
    displaced_vectors: np.ndarray | list[list[float]],
    reference_energies: list[float] | np.ndarray,
    displaced_energies: list[float] | np.ndarray,
    min_overlap: float = 0.85,
    min_overlap_gap: float = 0.05,
    min_energy_gap: float = 1.0e-4,
) -> MRSFRootTrackingDiagnostic:
    """Assign a displaced MRSF root to the requested reference root.

    The requested state index is preserved. The displaced root is assigned by
    absolute vector overlap against the requested reference vector. A mismatch,
    weak/ambiguous overlap, or near-degenerate assigned root fails loudly.
    """

    if requested_state < 1:
        raise ValueError("MRSF Hessian root tracking requires a positive physical state index")

    ref_vec = _as_2d_vectors(reference_vectors, name="reference_vectors")
    disp_vec = _as_2d_vectors(displaced_vectors, name="displaced_vectors")
    ref_e = np.asarray(reference_energies, dtype=float)
    disp_e = np.asarray(displaced_energies, dtype=float)

    nstate = min(ref_vec.shape[0], disp_vec.shape[0], ref_e.size, disp_e.size)
    if requested_state >= nstate:
        raise ValueError(
            f"requested root {requested_state} is outside available roots 0..{nstate - 1}"
        )

    overlaps = np.abs(disp_vec[:nstate] @ ref_vec[requested_state])
    order = np.argsort(overlaps)[::-1]
    matched_state = int(order[0])
    best_overlap = float(overlaps[matched_state])
    second_overlap = float(overlaps[order[1]]) if nstate > 1 else 0.0
    overlap_gap = best_overlap - second_overlap

    assigned_energy_gaps = np.abs(disp_e[:nstate] - disp_e[matched_state])
    nonzero_gaps = assigned_energy_gaps[assigned_energy_gaps > 0.0]
    energy_gap = float(nonzero_gaps.min()) if nonzero_gaps.size else float("inf")

    diagnostic = MRSFRootTrackingDiagnostic(
        requested_state=requested_state,
        matched_state=matched_state,
        overlap=best_overlap,
        overlap_gap=overlap_gap,
        energy_gap=energy_gap,
        root_flip=(matched_state != requested_state),
        ok=(matched_state == requested_state and best_overlap >= min_overlap and overlap_gap >= min_overlap_gap and energy_gap >= min_energy_gap),
    )

    if diagnostic.root_flip:
        raise MRSFRootTrackingError(
            f"MRSF root flip detected: requested root {requested_state} "
            f"({physical_root_label(requested_state)}) matched root {matched_state} "
            f"with overlap {best_overlap:.6f}."
        )
    if best_overlap < min_overlap or overlap_gap < min_overlap_gap:
        raise MRSFRootTrackingError(
            f"MRSF ambiguous vector overlap for requested root {requested_state}: "
            f"best overlap {best_overlap:.6f}, overlap gap {overlap_gap:.6f}."
        )
    if energy_gap < min_energy_gap:
        raise MRSFRootTrackingError(
            f"MRSF ambiguous near-degenerate energy gap for requested root {requested_state}: "
            f"minimum displaced gap {energy_gap:.6e} Hartree."
        )

    return diagnostic
