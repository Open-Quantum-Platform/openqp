"""Excited-state character descriptors for MRSF states.

* ``participation_ratio`` -- NTO participation ratio PR = (sum w)^2 / sum w^2.
* ``tozer_lambda`` -- Tozer's Lambda: NTO-weighted spatial overlap of the
  hole/particle moduli (Lambda in [0,1]; high = local, low = charge transfer).
* ``fragment_ct_matrix`` -- Plasser/Lischka Omega matrix from the Loewdin-
  orthogonalized spin-flip transition density, partitioned over atom fragments.
"""
import numpy as np
from scipy.linalg import sqrtm

__all__ = ["participation_ratio", "tozer_lambda", "fragment_ct_matrix"]


def participation_ratio(weights):
    """PR = (sum sigma^2)^2 / sum sigma^4 from NTO weights (= sigma^2)."""
    w = np.asarray(weights, dtype=float)
    num = w.sum() ** 2
    den = (w ** 2).sum()
    return float(num / den) if den > 0 else 0.0


def tozer_lambda(ao, nto_exc, grid_points, dV, weight_thresh=1.0e-8):
    """Tozer Lambda for an excitation, from its hole/particle NTOs.

    Lambda = sum_k w_k O_k / sum_k w_k, with O_k = \\int |phi_hole,k||phi_part,k|.
    The NTO orbitals are valence and smooth, so a moderate grid suffices."""
    w = nto_exc["weights"]
    keep = w > weight_thresh * (w.max() if w.size else 1.0)
    holes = nto_exc["holes_ao"][:, keep]
    parts = nto_exc["particles_ao"][:, keep]
    wk = w[keep]
    # evaluate NTOs on the grid in chunks
    npts = grid_points.shape[0]
    O = np.zeros(holes.shape[1])
    norm_h = np.zeros(holes.shape[1])
    norm_p = np.zeros(holes.shape[1])
    chunk = 300000
    for i in range(0, npts, chunk):
        chi = ao.eval_ao(grid_points[i:i + chunk])      # (npc, nbf)
        ph = chi @ holes                                 # (npc, k)
        pp = chi @ parts
        O += np.sum(np.abs(ph) * np.abs(pp), axis=0) * dV
        norm_h += np.sum(ph * ph, axis=0) * dV
        norm_p += np.sum(pp * pp, axis=0) * dV
    # guard against grid under-normalisation
    O = O / np.sqrt(norm_h * norm_p)
    lam = float(np.sum(wk * O) / np.sum(wk))
    return lam, {"O_k": O, "weights": wk, "norm_hole": norm_h, "norm_part": norm_p}


def fragment_ct_matrix(states, ao, n, fragments):
    """Fragment charge-transfer Omega matrix for excitation to root ``n``.

    ``fragments``: list of lists of 0-based atom indices.  Builds the AO
    spin-flip transition density T = C_occ X^{(n)} C_vir^T, Loewdin-orthogonalizes
    (Ttil = S^{1/2} T S^{1/2}), and sums Ttil_{mu,nu}^2 over fragment blocks.
    Row index = hole fragment, column index = particle fragment."""
    X = states.amplitude_matrix(n)                      # (noca, nvirb)
    C = states.C
    na, nb, nbf = states.na, states.nb, states.nbf
    T_ao = C[:, :na] @ X @ C[:, nb:nbf].T               # (nbf, nbf)
    Shalf = np.real(sqrtm(states.S))
    Ttil = Shalf @ T_ao @ Shalf
    P = Ttil ** 2                                       # element-wise CT weights

    nfrag = len(fragments)
    # map atom -> fragment
    atom_frag = {}
    for f, atoms in enumerate(fragments):
        for a in atoms:
            atom_frag[int(a)] = f
    ao_frag = np.array([atom_frag[int(a)] for a in ao.ao_atom])

    Omega = np.zeros((nfrag, nfrag))
    for A in range(nfrag):
        rows = ao_frag == A
        for B in range(nfrag):
            cols = ao_frag == B
            Omega[A, B] = P[np.ix_(rows, cols)].sum()
    total = P.sum()
    return {
        "Omega": Omega,
        "total": float(total),
        "hole_pop": Omega.sum(axis=1),       # per-fragment hole population
        "particle_pop": Omega.sum(axis=0),   # per-fragment particle population
        "ct_fraction": float(1.0 - np.trace(Omega) / total) if total > 0 else 0.0,
        "amplitude_norm": float((X ** 2).sum()),
    }
