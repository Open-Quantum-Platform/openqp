"""Natural transition orbitals (NTOs) for MRSF excited states.

Two complementary, well-defined objects (the MRSF S0-as-root subtlety means they
do *not* coincide as they would in closed-shell TDDFT):

* :func:`nto_excitation` -- SVD of the spin-adapted spin-flip amplitude matrix
  ``X^{(n)}`` of root *n* (alpha-occupied -> beta-virtual).  These are the
  hole/particle NTOs describing the *character* of state *n* relative to the
  high-spin reference; ``sum(sigma^2) == ||X^{(n)}||^2`` (~1 for a normalized
  single root).  Holes are alpha, particles are beta (spin-resolved).
* :func:`nto_transition` -- SVD of the state-interaction 1-TDM
  ``gamma^{i->j}`` (validated at GATE 2).  Truncated reconstruction reproduces
  the transition dipole; these are the natural orbitals of the genuine S_i->S_j
  transition density (which has occ-occ/vir-vir structure for MRSF).
"""
import numpy as np

__all__ = ["nto_excitation", "nto_transition"]


def nto_excitation(states, n, thresh=1.0e-6):
    """Hole/particle NTOs of root ``n`` from its spin-flip amplitude matrix.

    Returns a dict with spin-resolved AO-basis NTO coefficients, singular
    values and weights (sigma^2), plus the amplitude norm sum(sigma^2)."""
    X = states.amplitude_matrix(n)                  # (noca, nvirb)
    U, sig, Wt = np.linalg.svd(X, full_matrices=False)
    W = Wt.T
    C = states.C
    na, nb, nbf = states.na, states.nb, states.nbf
    C_occ = C[:, :na]                               # alpha-occupied
    C_vir = C[:, nb:nbf]                            # beta-virtual (spatial = alpha)
    holes_ao = C_occ @ U                            # (nbf, r) alpha holes
    parts_ao = C_vir @ W                            # (nbf, r) beta particles
    weights = sig ** 2
    keep = weights > thresh * weights.max() if weights.size else np.array([], bool)
    return {
        "state": n,
        "sigma": sig,
        "weights": weights,
        "sum_sigma2": float(weights.sum()),
        "holes_ao": holes_ao,
        "particles_ao": parts_ao,
        "hole_spin": "alpha",
        "particle_spin": "beta",
        "n_significant": int(keep.sum()),
        "amplitude_matrix": X,
    }


def nto_transition(states, i, j, thresh=1.0e-6):
    """NTOs of the state-interaction 1-TDM gamma^{i->j} (alpha-MO basis SVD).

    The left/right singular vectors are the particle/hole natural orbitals of
    the S_i->S_j transition density.  ``reconstruct_tdm`` rebuilds gamma from the
    leading pairs so the transition dipole can be re-derived (GATE 3 check)."""
    g = states.tdm_mo(i, j)                         # (nbf, nbf), alpha-MO basis
    # Particles = left vectors (target index), holes = right vectors (source).
    U, sig, Wt = np.linalg.svd(g)
    W = Wt.T
    C = states.C
    particles_ao = C @ U                            # (nbf, nbf)
    holes_ao = C @ W
    weights = sig ** 2
    keep = weights > thresh * weights.max() if weights.size else np.array([], bool)

    def reconstruct_tdm(rank=None):
        r = sig.size if rank is None else rank
        g_mo = (U[:, :r] * sig[:r]) @ Wt[:r, :]
        return g_mo

    return {
        "pair": (i, j),
        "sigma": sig,
        "weights": weights,
        "sum_sigma2": float(weights.sum()),
        "holes_ao": holes_ao,
        "particles_ao": particles_ao,
        "n_significant": int(keep.sum()),
        "reconstruct_tdm": reconstruct_tdm,
    }
