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

# Largest p - alpha for which the xi^alpha kernel is validated
# (source/dftlib/xi_alpha_kernel.F90, XI_C_MAX).
XI_C_MAX = 2.0


def grid_ingredients(mol, alphas, p=None, scale=0, cutoff=0.0):
    """Return a dict of numpy arrays (points along the first axis)."""
    alphas = np.ascontiguousarray(np.asarray(alphas, dtype=np.float64).ravel())
    n = alphas.size
    if n == 0:
        raise ValueError("grid_ingredients: at least one alpha is required")
    if not np.all(np.isfinite(alphas)):
        raise ValueError("grid_ingredients: alpha must be finite")
    if np.any(alphas > 1.0):
        raise ValueError("grid_ingredients: alpha must be <= 1")
    if p is None:
        ps = np.full(n, -1, dtype=np.int64)
    else:
        p_raw = np.asarray(p).ravel()
        if p_raw.dtype.kind not in "iu" and not np.all(np.equal(np.mod(p_raw, 1), 0)):
            raise ValueError("grid_ingredients: p entries must be integers")
        ps = np.ascontiguousarray(p_raw.astype(np.int64))
        if ps.size != n:
            raise ValueError("grid_ingredients: p must have one entry per alpha")
        if not np.all(np.isin(ps, (-1, 0, 1))):
            raise ValueError("grid_ingredients: p entries must be -1 (auto), 0 or 1")
        if np.any((ps == 0) & (alphas > 0.0)):
            raise ValueError("grid_ingredients: p=0 requires alpha <= 0")
    if isinstance(scale, bool) or scale not in (0, 1):
        raise ValueError("grid_ingredients: scale must be 0 or 1")
    cutoff = float(cutoff)
    if not np.isfinite(cutoff) or cutoff < 0.0:
        raise ValueError("grid_ingredients: cutoff must be a finite number >= 0 (bohr; 0 = exact)")
    ps = np.where(ps < 0, np.where(alphas > 0.0, 1, 0), ps).astype(np.int64)
    if np.any(ps - alphas > XI_C_MAX):
        raise ValueError("grid_ingredients: p - alpha must be <= 2 (validated kernel range; "
                         "alpha >= -2 with p=0, alpha >= -1 with p=1)")
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
