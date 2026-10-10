"""Helpers for presenting internal OpenQP arrays in JSON."""

import numpy as np


_TD_BVEC_KEYS = {
    "OQP::td_bvec_mo",
    "OQP::td_bvec_mo_old",
    "OQP::td_bvec_mo_s",
    "OQP::td_bvec_mo_t",
}


def json_array(key, value):
    """Return an array with axes laid out as documented for JSON output.

    TD response vectors originate in Fortran as ``(N_DRF, N_States)``.  The
    tag-array bridge exposes the contiguous buffer with that shape to NumPy,
    whose default C-order interpretation groups adjacent values by the second
    axis.  Rebuild these arrays from their state-major buffer so JSON rows are
    DRFs and columns are states.  This conversion is intentionally limited to
    serialization; the runtime/tag-array representation remains unchanged.
    """
    array = np.asarray(value)
    if key in _TD_BVEC_KEYS and array.ndim == 2:
        n_drf, n_states = array.shape
        array = array.ravel().reshape(n_states, n_drf).T
    return array.tolist()


def tag_array_from_json(key, value):
    """Restore the runtime tag-array view from documented JSON axes.

    This is the exact inverse of :func:`json_array` for TD response vectors.
    Without it, a ``guess=json`` restart compares the current state-major
    runtime buffer with a previous DRF-by-state JSON matrix as if their memory
    layouts were identical, which can invent a root exchange.
    """
    array = np.asarray(value)
    if key in _TD_BVEC_KEYS and array.ndim == 2:
        n_drf, n_states = array.shape
        array = np.ascontiguousarray(array.T).reshape(n_drf, n_states)
    return array
