"""Grid ingredient export: n, sigma, tau and xi^alpha for several alpha at once.

Usage (after a converged SCF, e.g. ``runner.run()``)::

    from oqp.library.xi_ingredients import grid_ingredients
    g = grid_ingredients(runner.mol, alphas=[1.0, 0.5, 0.0, -0.5], p=[1, 1, 0, 0])
    g["xi"].shape      # (npts, 2*len(alphas)) : alpha-spin, beta-spin per alpha
    g["tau"], g["rho"], g["sigma"], g["xyzw"]

The ingredient is evaluated on the DFT grid of the molecule ([dftgrid]) with
the converged density matrices; the functional attached to the molecule is
not modified.
"""

import numpy as np

import oqp


def grid_ingredients(mol, alphas, p=None, scale=0, cutoff=0.0):
    """Return a dict of numpy arrays (points along the first axis)."""
    alphas = np.ascontiguousarray(np.asarray(alphas, dtype=np.float64).ravel())
    n = alphas.size
    if n == 0:
        raise ValueError("grid_ingredients: at least one alpha is required")
    if np.any(alphas > 1.0):
        raise ValueError("grid_ingredients: alpha must be <= 1")
    if p is None:
        ps = np.full(n, -1, dtype=np.int64)
    else:
        ps = np.ascontiguousarray(np.asarray(p, dtype=np.int64).ravel())
        if ps.size != n:
            raise ValueError("grid_ingredients: p must have one entry per alpha")
        if np.any((ps == 0) & (alphas > 0.0)):
            raise ValueError("grid_ingredients: p=0 requires alpha <= 0")
    ps = np.where(ps < 0, np.where(alphas > 0.0, 1, 0), ps).astype(np.int64)
    oqp.oqp_xi_grid_ingredients(
        mol,
        int(n),
        oqp.ffi.cast("double *", alphas.ctypes.data),
        oqp.ffi.cast("int64_t *", ps.ctypes.data),
        int(scale),
        float(cutoff),
    )
    data = mol.data
    out = {
        "alpha": np.array(data["OQP::XI_GRID_ALPHA"]).ravel(),
        "p": np.array(data["OQP::XI_GRID_P"]).ravel().astype(int),
        "xyzw": np.array(data["OQP::XI_GRID_XYZW"]).reshape(-1, 4),
        "rho": np.array(data["OQP::XI_GRID_RHO"]).reshape(-1, 2),
        "sigma": np.array(data["OQP::XI_GRID_SIGMA"]).reshape(-1, 3),
        "tau": np.array(data["OQP::XI_GRID_TAU"]).reshape(-1, 2),
        "xi": np.array(data["OQP::XI_GRID_XI"]).reshape(-1, 2 * n),
    }
    return out
