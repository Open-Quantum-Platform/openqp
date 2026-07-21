"""Excited-state character descriptors for MRSF states.

* ``participation_ratio`` -- NTO participation ratio PR = (sum w)^2 / sum w^2.
* ``tozer_lambda`` -- Tozer's Lambda: NTO-weighted spatial overlap of the
  hole/particle moduli (Lambda in [0,1]; high = local, low = charge transfer).
* ``fragment_ct_matrix`` -- Plasser/Lischka Omega matrix from the Loewdin-
  orthogonalized physical-root state-interaction 1-TDM, partitioned over atom
  fragments.
"""
import numpy as np

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


def _lowdin_sqrt(overlap, threshold=1.0e-12):
    """Symmetric positive-semidefinite square root of an AO overlap matrix."""
    values, vectors = np.linalg.eigh(np.asarray(overlap, dtype=float))
    if values.size and values.min() < -threshold:
        raise ValueError("AO overlap matrix is not positive semidefinite")
    values = np.clip(values, 0.0, None)
    return (vectors * np.sqrt(values)) @ vectors.T


def _fragment_map(ao_atom, fragments):
    """Validate a disjoint, complete atom partition and return AO fragments."""
    atom_frag = {}
    for frag, atoms in enumerate(fragments):
        if not atoms:
            raise ValueError(f"fragment {frag} is empty")
        for atom in atoms:
            atom = int(atom)
            if atom in atom_frag:
                raise ValueError(f"atom {atom} occurs in more than one fragment")
            atom_frag[atom] = frag
    missing = sorted(set(np.asarray(ao_atom, dtype=int)) - set(atom_frag))
    if missing:
        raise ValueError(f"fragments do not assign AO-bearing atoms {missing}")
    return np.array([atom_frag[int(atom)] for atom in ao_atom], dtype=int)


def fragment_ct_matrix(states, ao, n, fragments, ref=0, norm_tolerance=1.0e-10):
    """Fragment charge-transfer matrix for the physical ``ref -> n`` transition.

    ``fragments`` is a disjoint partition of 0-based atom indices.  The genuine
    MRSF state-interaction 1-TDM is Loewdin-orthogonalized as
    ``Ttilde = S^(1/2) T_AO S^(1/2)`` and

    ``Omega[A, B] = sum_(mu in B, nu in A) |Ttilde[mu, nu]|^2``.

    Rows therefore identify the *hole* fragment and columns the *particle*
    fragment.  This function deliberately does not use ``amplitude_matrix``:
    that matrix describes a response root relative to the auxiliary high-spin
    determinant, whereas ``tdm_ao(ref, n)`` describes the requested pair of
    physical MRSF roots.

    ``total`` is the squared 1-TDM norm that normalizes LE/CT.  A symmetry-
    forbidden pair has no appreciable 1-TDM, so below ``norm_tolerance`` the
    fractions would be numerical noise dressed up as three decimal places:
    ``le_fraction``/``ct_fraction`` are then NaN and ``well_defined`` is False.

    .. note::
       This is a deliberate behaviour change from the pre-0.9 version, which
       analyzed ``amplitude_matrix(n)`` -- a response root relative to the
       auxiliary high-spin determinant, not a physical root pair.  The
       ``amplitude_norm`` result key belonged to that object and is gone; use
       :func:`~oqp.analysis.nto_excitation` for amplitude-based character.
    """
    ref = int(ref)
    n = int(n)
    if ref == n:
        raise ValueError("fragment CT analysis requires two different states")
    if not (0 <= ref < states.nstates and 0 <= n < states.nstates):
        raise IndexError("state index is outside the MRSF root range")

    T_ao = np.asarray(states.tdm_ao(ref, n), dtype=float)
    if T_ao.shape != (states.nbf, states.nbf):
        raise ValueError("state-interaction AO transition density has the wrong shape")
    Shalf = _lowdin_sqrt(states.S)
    Ttil = Shalf @ T_ao @ Shalf
    P = np.abs(Ttil) ** 2

    nfrag = len(fragments)
    ao_frag = _fragment_map(ao.ao_atom, fragments)

    Omega = np.zeros((nfrag, nfrag))
    for A in range(nfrag):
        hole_columns = ao_frag == A
        for B in range(nfrag):
            particle_rows = ao_frag == B
            Omega[A, B] = P[np.ix_(particle_rows, hole_columns)].sum()
    total = float(P.sum())
    well_defined = total > float(norm_tolerance)
    local = float(np.trace(Omega) / total) if well_defined else float("nan")
    return {
        "pair": (ref, n),
        "Omega": Omega,
        "total": total,
        "hole_pop": Omega.sum(axis=1),       # per-fragment hole population
        "particle_pop": Omega.sum(axis=0),   # per-fragment particle population
        "le_fraction": local,
        "ct_fraction": 1.0 - local,
        "well_defined": well_defined,
        "transition_density_orthogonalized": Ttil,
    }
