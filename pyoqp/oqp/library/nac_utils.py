"""Small convention helpers for overlap-based nonadiabatic couplings."""

from __future__ import annotations

import numpy as np


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
