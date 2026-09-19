"""Small convention helpers for overlap-based nonadiabatic couplings."""

from __future__ import annotations

import json
import os

import numpy as np


# Numerical-NAC scratch files predate the canonical PR #160 convention.  A
# sidecar marker prevents restart=True from silently mixing those historical
# files with current workers after a partial restart.
NUMERICAL_NAC_CACHE_SCHEMA_VERSION = 2
NUMERICAL_NAC_CACHE_CONVENTION = (
    "canonical-old-new__retained-column-normalized__"
    "s-minus-st-over-2step__central-root-order-v2"
)


def numerical_nac_cache_metadata(step, nstate, worker_index):
    """Return the complete identity of a reusable numerical-NAC worker file."""
    return {
        "schema_version": NUMERICAL_NAC_CACHE_SCHEMA_VERSION,
        "convention": NUMERICAL_NAC_CACHE_CONVENTION,
        "step": float(step),
        "nstate": int(nstate),
        "worker_index": int(worker_index),
    }


def load_numerical_nac_cache(dat, step, nstate, worker_index):
    """Load a current, complete numerical-NAC cache entry or return None."""
    marker = f"{dat}.meta.json"
    try:
        with open(marker, "r") as handle:
            metadata = json.load(handle)
        if metadata != numerical_nac_cache_metadata(step, nstate, worker_index):
            return None
        dcme = np.loadtxt(dat).reshape(-1)
        if dcme.size != nstate * nstate or not np.all(np.isfinite(dcme)):
            return None
        return dcme
    except (OSError, TypeError, ValueError):
        return None


def write_numerical_nac_cache_marker(dat, step, nstate, worker_index):
    """Atomically mark a successfully read worker file as current."""
    marker = f"{dat}.meta.json"
    temporary = f"{marker}.{os.getpid()}.tmp"
    with open(temporary, "w") as handle:
        json.dump(
            numerical_nac_cache_metadata(step, nstate, worker_index), handle,
            indent=2, sort_keys=True,
        )
        handle.write("\n")
    os.replace(temporary, marker)


def canonical_state_overlap(native_overlap):
    """Return ``S[i,j] = <old i | new j>`` from an OpenQP native tag.

    OpenQP's two-dimensional tag buffers expose Fortran column-major storage
    through a C-order NumPy view.  The resulting Python array is therefore the
    transpose of the Fortran state-overlap matrix.  The DFTB/xTB adapters
    deliberately return the same native tag layout.
    """
    native = np.asarray(native_overlap, dtype=float)
    if native.ndim != 2 or native.shape[0] != native.shape[1]:
        raise ValueError(
            f"state overlap must be a square matrix, got {native.shape}"
        )
    return np.array(native.T, copy=True)


def normalize_retained_state_overlap(state_overlap):
    """Normalize retained-state columns for a spatial finite-difference TLF.

    The published spatial TLF vector treats each displaced retained-state
    projection as a normalized electronic state before the +dx/-dx central
    difference.  Time-dependent propagation is deliberately different: its
    finite-window overlap remains raw so population outside the retained
    manifold is not folded back into it.
    """
    overlap = np.asarray(state_overlap, dtype=float)
    if overlap.ndim != 2 or overlap.shape[0] != overlap.shape[1]:
        raise ValueError(
            f"state overlap must be a square matrix, got {overlap.shape}"
        )
    norms = np.linalg.norm(overlap, axis=0)
    if np.any(norms <= np.finfo(float).tiny):
        raise ValueError("retained-state overlap contains a zero-norm column")
    return overlap / norms.reshape(1, -1)


def hst_derivative_coupling(state_overlap, step):
    """Hammes-Schiffer--Tully coupling from a canonical state overlap."""
    if step == 0:
        raise ValueError("the NAC time/displacement step must be nonzero")
    overlap = np.asarray(state_overlap, dtype=float)
    if overlap.ndim != 2 or overlap.shape[0] != overlap.shape[1]:
        raise ValueError(
            f"state overlap must be a square matrix, got {overlap.shape}"
        )
    coupling = (overlap - overlap.T) / (2.0 * step)
    np.fill_diagonal(coupling, 0.0)
    return coupling


def interstate_coupling(derivative_coupling, energies):
    """Return ``h_ij = (E_j - E_i) d_ij`` for a state-pair tensor."""
    derivative = np.asarray(derivative_coupling, dtype=float)
    energy = np.asarray(energies, dtype=float).reshape(-1)
    nstate = energy.size
    if derivative.ndim < 2 or derivative.shape[:2] != (nstate, nstate):
        raise ValueError(
            "derivative coupling leading dimensions must match the energies: "
            f"{derivative.shape} versus {nstate} states"
        )
    gap = energy.reshape((1, nstate)) - energy.reshape((nstate, 1))
    gap = gap.reshape((nstate, nstate) + (1,) * (derivative.ndim - 2))
    return derivative * gap
