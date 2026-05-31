"""Utilities for stable overlap-based numerical nonadiabatic couplings."""

from itertools import permutations

import numpy as np


def _assignment_max_abs(overlap):
    """Return a column order that maximizes diagonal absolute overlap.

    The state counts in OpenQP excited-state NAC jobs are normally small. Use an
    exact exhaustive assignment for small matrices and a deterministic greedy
    fallback for larger matrices to avoid adding a SciPy dependency in the core
    runtime path.
    """
    score = np.abs(np.asarray(overlap, dtype=float))
    nrow, ncol = score.shape
    if nrow != ncol:
        raise ValueError(f"state overlap must be square, got {score.shape}")
    n = nrow

    if n <= 8:
        best_order = None
        best_score = -np.inf
        for order in permutations(range(n)):
            total = float(sum(score[i, order[i]] for i in range(n)))
            if total > best_score:
                best_score = total
                best_order = order
        return np.array(best_order, dtype=int)

    remaining = set(range(n))
    order = np.empty(n, dtype=int)
    for i in range(n):
        j = max(remaining, key=lambda col: score[i, col])
        order[i] = j
        remaining.remove(j)
    return order


def closest_orthogonal(matrix):
    """Return the closest orthogonal matrix to ``matrix`` in Frobenius norm."""
    u, _, vt = np.linalg.svd(np.asarray(matrix, dtype=float), full_matrices=False)
    return u @ vt


def align_state_overlap(overlap, min_match=0.5):
    """Gauge-align a real state-overlap matrix for stable numerical NAC.

    Parameters
    ----------
    overlap
        Square matrix whose rows are previous/reference states and whose columns
        are current/displaced states.
    min_match
        Diagnostic threshold for the assigned absolute diagonal overlaps. Values
        below this indicate ambiguous root tracking or too-large nuclear steps.

    Returns
    -------
    aligned, diagnostics
        ``aligned`` has columns permuted and sign-corrected so assigned diagonal
        overlaps are positive. ``diagnostics`` records the permutation, signs,
        assigned overlaps, and an ambiguity flag.
    """
    raw = np.asarray(overlap, dtype=float)
    if raw.ndim != 2 or raw.shape[0] != raw.shape[1]:
        raise ValueError(f"state overlap must be a square matrix, got {raw.shape}")

    order = _assignment_max_abs(raw)
    aligned = raw[:, order].copy()
    diag = np.diag(aligned).copy()
    signs = np.where(diag < 0.0, -1.0, 1.0)
    # Treat exact zeros as +1 to avoid erasing a column.
    signs[signs == 0.0] = 1.0
    aligned *= signs.reshape((1, -1))
    assigned = np.diag(aligned).copy()

    diagnostics = {
        "permutation": order,
        "signs": signs,
        "assigned_overlaps": assigned,
        "min_abs_assigned_overlap": float(np.min(np.abs(assigned))) if assigned.size else 0.0,
        "max_offdiag_overlap": float(np.max(np.abs(aligned - np.diag(np.diag(aligned))))) if aligned.size else 0.0,
        "orthogonality_error_raw": float(np.linalg.norm(raw.T @ raw - np.eye(raw.shape[0]))),
        "ambiguous": bool(np.any(np.abs(assigned) < min_match)),
    }
    return aligned, diagnostics


def overlap_to_derivative_coupling(overlap, dt, method="polar"):
    """Convert a gauge-aligned overlap matrix to derivative couplings.

    ``method='antisym'`` is the Hammes-Schiffer-Tully small-step formula,
    ``(S-S^T)/(2*dt)``. ``method='polar'`` first projects the overlap onto the
    nearest orthogonal matrix before applying the same anti-symmetrization,
    suppressing symmetric/basis-following noise and arbitrary sign flips.
    """
    if dt == 0:
        raise ValueError("dt must be nonzero for numerical NAC")

    s = np.asarray(overlap, dtype=float)
    if method == "antisym":
        used = s
    elif method == "polar":
        used = closest_orthogonal(s)
    else:
        raise ValueError(f"unknown overlap-to-NAC method: {method}")

    dc = (used - used.T) / (2.0 * dt)
    np.fill_diagonal(dc, 0.0)
    return dc, used


def stable_nac_from_overlap(overlap, dt, method="polar", min_match=0.5):
    """Align a state-overlap matrix and compute a phase-stable NACME matrix."""
    aligned, diagnostics = align_state_overlap(overlap, min_match=min_match)
    dc, used = overlap_to_derivative_coupling(aligned, dt, method=method)
    diagnostics = dict(diagnostics)
    diagnostics.update(
        {
            "method": method,
            "orthogonality_error_aligned": float(np.linalg.norm(aligned.T @ aligned - np.eye(aligned.shape[0]))),
            "orthogonality_error_used": float(np.linalg.norm(used.T @ used - np.eye(used.shape[0]))),
            "antisymmetry_error": float(
                np.linalg.norm(dc + dc.T) / (np.linalg.norm(dc) + 1.0e-30)
            ),
        }
    )
    return dc, aligned, used, diagnostics


def symmetrize_derivative_coupling(dc):
    """Project derivative couplings onto the antisymmetric state-pair part."""
    arr = np.asarray(dc, dtype=float)
    out = 0.5 * (arr - np.swapaxes(arr, 0, 1))
    if out.ndim >= 2 and out.shape[0] == out.shape[1]:
        idx = np.arange(out.shape[0])
        out[idx, idx, ...] = 0.0
    return out


def verify_nac_conventions(dc, nac, atol=1.0e-6, label="", strict=True):
    """Check derivative/interstate NAC symmetry conventions.

    Returns the relative residuals ``(d_res, h_res)``. With ``strict=True`` an
    ``AssertionError`` is raised when either residual exceeds ``atol``; with
    ``strict=False`` callers can log diagnostics without killing production NAMD.
    """
    dc = np.asarray(dc, dtype=float)
    nac = np.asarray(nac, dtype=float)
    d_res = np.linalg.norm(dc + np.swapaxes(dc, 0, 1)) / (np.linalg.norm(dc) + 1.0e-30)
    h_res = np.linalg.norm(nac - np.swapaxes(nac, 0, 1)) / (np.linalg.norm(nac) + 1.0e-30)
    if strict:
        tag = f" ({label})" if label else ""
        assert d_res < atol, \
            f"NAC convention{tag}: derivative coupling d not antisymmetric, ||d+d^T||/||d||={d_res:.2e}"
        assert h_res < atol, \
            f"NAC convention{tag}: interstate coupling h not symmetric, ||h-h^T||/||h||={h_res:.2e}"
    return d_res, h_res
