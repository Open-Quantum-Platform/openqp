"""Attachment / detachment densities for MRSF excited states (unrelaxed).

Given two roots ``i`` (reference state, default S0) and ``n``, the difference
density is ``Delta = gamma^n - gamma^i``.  Because both state 1-RDMs share the
same high-spin reference, the reference cancels and
``Delta = Delta gamma^n - Delta gamma^i`` (the two exposed traceless diagonal
blocks).  Diagonalizing the symmetric ``Delta`` (in the orthonormal MO basis)
and splitting by eigenvalue sign yields the detachment (negative) and
attachment (positive) densities of Head-Gordon et al.
"""
import numpy as np

__all__ = ["attachment_detachment"]


def attachment_detachment(states, n, ref=0):
    """Unrelaxed attachment/detachment analysis for the ref -> n transition.

    Returns A, D (AO-basis density matrices), the natural orbitals/occupations,
    and the promotion number n_promoted = Tr(A) = Tr(D)."""
    delta_mo = states.diff_density_mo(n) - states.diff_density_mo(ref)
    delta_mo = 0.5 * (delta_mo + delta_mo.T)          # symmetrize (numerical)
    w, V = np.linalg.eigh(delta_mo)                   # MO-basis, orthonormal

    pos = w > 0
    neg = w < 0
    A_mo = (V[:, pos] * w[pos]) @ V[:, pos].T         # attachment (positive)
    D_mo = (V[:, neg] * (-w[neg])) @ V[:, neg].T      # detachment (positive def)

    C = states.C
    A_ao = C @ A_mo @ C.T
    D_ao = C @ D_mo @ C.T

    n_attach = float(w[pos].sum())
    n_detach = float(-w[neg].sum())

    return {
        "state": n,
        "ref": ref,
        "delta_mo": delta_mo,
        "eigenvalues": w,
        "natural_orbitals_ao": C @ V,                 # attach/detach natural orbitals
        "A_mo": A_mo,
        "D_mo": D_mo,
        "A_ao": A_ao,
        "D_ao": D_ao,
        "n_attach": n_attach,
        "n_detach": n_detach,
        "n_promoted": 0.5 * (n_attach + n_detach),
    }
