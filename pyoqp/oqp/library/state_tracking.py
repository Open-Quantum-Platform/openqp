"""Thin Python boundary for OpenQP's native maximum-overlap tracker."""

from __future__ import annotations

import numpy as np

import oqp


def maximum_overlap_assignment(overlap):
    """Return global row-to-column correspondence and phase diagnostics.

    The numerical assignment is performed by the resident Fortran kernel.
    Returned columns are zero-based.  ``matched`` is the selected absolute
    overlap and ``margin`` is selected minus the best unselected value in the
    same row; a small or negative margin flags an ambiguous correspondence.
    """
    matrix = np.ascontiguousarray(overlap, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"tracking overlap must be square, got {matrix.shape}")
    n = matrix.shape[0]
    if n == 0:
        empty_i = np.empty(0, dtype=np.int32)
        empty_f = np.empty(0, dtype=np.float64)
        return empty_i, empty_f.copy(), empty_f.copy(), empty_f.copy()

    assignment = np.empty(n, dtype=np.int32)
    signs = np.empty(n, dtype=np.float64)
    matched = np.empty(n, dtype=np.float64)
    margins = np.empty(n, dtype=np.float64)
    info = oqp.lib.oqp_maximum_overlap_assignment(
        n,
        oqp.ffi.cast("double *", matrix.ctypes.data),
        oqp.ffi.cast("int *", assignment.ctypes.data),
        oqp.ffi.cast("double *", signs.ctypes.data),
        oqp.ffi.cast("double *", matched.ctypes.data),
        oqp.ffi.cast("double *", margins.ctypes.data),
    )
    if info != 0:
        raise ValueError(f"native maximum-overlap assignment failed (info={info})")
    return assignment, signs, matched, margins
