"""Physical-root MRSF state analysis.

Spatial LE/CT character and orbital character answer different questions and
are therefore reported separately.  LE/CT comes from a fragment-partitioned
state-interaction 1-TDM.  The ``n/pi/sigma -> pi*/sigma*`` fractions come from
the same physical-root NTO pairs after a Loewdin population analysis.
"""
from __future__ import annotations

import numpy as np

from .descriptors import _lowdin_sqrt, fragment_ct_matrix, participation_ratio
from .nto import nto_transition

__all__ = ["analyze_mrsf_transition", "infer_molecular_plane"]

_CHANNELS = (
    "n->pi*", "pi->pi*", "sigma->pi*",
    "n->sigma*", "pi->sigma*", "sigma->sigma*",
)
_PRETTY = {
    "n->pi*": "nπ*", "pi->pi*": "ππ*", "sigma->pi*": "σπ*",
    "n->sigma*": "nσ*", "pi->sigma*": "πσ*", "sigma->sigma*": "σσ*",
}


def _unit_vector(vector, name):
    vector = np.asarray(vector, dtype=float).reshape(3)
    norm = np.linalg.norm(vector)
    if not np.isfinite(norm) or norm <= 1.0e-12:
        raise ValueError(f"{name} must be a finite nonzero 3-vector")
    return vector / norm


def infer_molecular_plane(ao, atoms=None, relative_tolerance=0.12):
    """Return the least-squares molecular-plane normal for a planar molecule.

    Hydrogen atoms are omitted when atomic numbers are supplied.  Linear and
    appreciably non-planar geometries are rejected because their pi direction
    is not defined by a single global plane; callers may provide an explicit
    local ``plane_normal`` to :func:`analyze_mrsf_transition` instead.
    """
    coords = np.asarray(ao.coords, dtype=float).reshape(-1, 3)
    if atoms is not None:
        atoms = np.asarray(atoms, dtype=int).reshape(-1)
        if atoms.size != coords.shape[0]:
            raise ValueError("atomic numbers and coordinates have different lengths")
        selected = coords[atoms > 1]
        if selected.shape[0] < 3:
            selected = coords
    else:
        selected = coords
    if selected.shape[0] < 3:
        raise ValueError("at least three atoms or an explicit plane_normal are required")
    centered = selected - selected.mean(axis=0)
    covariance = centered.T @ centered / selected.shape[0]
    values, vectors = np.linalg.eigh(covariance)
    if values[1] <= 1.0e-12 * max(values[2], 1.0):
        raise ValueError("a linear geometry has no unique molecular plane")
    extent = np.sqrt(max(values[2], 1.0e-30))
    rms_out_of_plane = np.sqrt(max(values[0], 0.0))
    if rms_out_of_plane / extent > relative_tolerance:
        raise ValueError(
            "geometry is not described by one molecular plane; provide a local plane_normal")
    return _unit_vector(vectors[:, 0], "inferred plane normal")


def _p_shell_groups(ao):
    """Group Cartesian p AOs by shell, with the axis of each component.

    A p shell has to be projected on the plane normal *coherently*, so the
    grouping is what makes the pi fraction essentially independent of how the
    molecule happens to be oriented in the input frame.  Shells with l >= 2
    carry no pi weight here and are counted as in-plane; polarization d
    functions hold little population, so this only blurs the n/sigma split
    slightly.

    "Essentially" because a residue remains: q = S^(1/2)c is covariant under
    rotations only when the AO basis is, and normalized *Cartesian* d shells
    are not (the six components carry an extra s-like combination and mix
    non-orthogonally under rotation).  Since ``AOBasis`` accepts only Cartesian
    bases, rigidly rotating a molecule still moves the reported fractions a
    little through S^(1/2): measured at 5.4e-3 for formaldehyde MRSF/6-31G*,
    and below 1.5e-13 for the same molecule in unpolarized 6-31G.  That is far
    below the ``label_threshold`` margin, but do not read the third decimal of
    a fraction as physical when polarization functions are present.
    """
    shells = {}
    for mu, (shell, powers, _scale) in enumerate(ao.ao_index):
        if sum(powers) == 1:
            shells.setdefault(shell, []).append((int(np.argmax(powers)), mu))
    return [(np.array([mu for _axis, mu in members], dtype=int),
             np.array([axis for axis, _mu in members], dtype=int))
            for members in shells.values()]


def _nto_orbital_character(coefficients, shalf, ao, normal, hetero_atoms,
                           p_groups):
    # q = S^(1/2)c are coefficients in the symmetric orthonormal AO basis.
    q = shalf @ np.asarray(coefficients, dtype=float)
    population = q * q
    norm = population.sum()
    if norm <= 1.0e-20:
        return {"n": 0.0, "pi": 0.0, "sigma": 0.0}
    hetero_mask = np.isin(ao.ao_atom, np.asarray(hetero_atoms, dtype=int))
    pi_population = 0.0
    lone_pair = 0.0
    in_p_shell = np.zeros(ao.nbf, dtype=bool)
    for indices, axes in p_groups:
        shell_population = float(population[indices].sum())
        # (n.q)^2, not sum_a n_a^2 q_a^2: dropping the p_x/p_y/p_z cross terms
        # would make the pi fraction correct only when the plane normal is a
        # Cartesian axis, and underestimate it by up to 3x otherwise.
        pi_shell = min(float(normal[axes] @ q[indices]) ** 2, shell_population)
        pi_population += pi_shell
        if hetero_mask[indices[0]]:      # one shell sits on one atom
            lone_pair += shell_population - pi_shell
        in_p_shell[indices] = True
    lone_pair += float(population[hetero_mask & ~in_p_shell].sum())
    pi = pi_population / norm
    n = min(max(lone_pair / norm, 0.0), max(1.0 - pi, 0.0))
    sigma = max(1.0 - pi - n, 0.0)
    total = n + pi + sigma
    return {"n": n / total, "pi": pi / total, "sigma": sigma / total}


def _orbital_channel_analysis(nto, states, ao, normal, hetero_atoms,
                              weight_threshold, norm_tolerance):
    weights = np.asarray(nto["weights"], dtype=float)
    if weights.size == 0 or weights.sum() <= norm_tolerance:
        # sum(sigma^2) is the squared 1-TDM norm.  Normalizing it away would
        # turn the numerical noise of a forbidden transition into fractions
        # that look precise; report "no character" instead.
        return {"fractions": {key: float("nan") for key in _CHANNELS},
                "pairs": [], "classified": False,
                "reason": "the ref->n transition density has no appreciable norm"}
    keep = weights > weight_threshold * weights.max()
    if not np.any(keep):
        keep[np.argmax(weights)] = True
    selected = np.flatnonzero(keep)
    selected_weights = weights[selected]
    selected_weights /= selected_weights.sum()
    shalf = _lowdin_sqrt(states.S)
    p_groups = _p_shell_groups(ao)
    fractions = np.zeros((3, 2))  # hole n/pi/sigma x particle pi*/sigma*
    pairs = []
    for rank, pair_weight in zip(selected, selected_weights):
        hole = _nto_orbital_character(
            nto["holes_ao"][:, rank], shalf, ao, normal, hetero_atoms, p_groups)
        particle_hps = _nto_orbital_character(
            nto["particles_ao"][:, rank], shalf, ao, normal, (), p_groups)
        particle = {"pi*": particle_hps["pi"],
                    "sigma*": particle_hps["n"] + particle_hps["sigma"]}
        h = np.array([hole["n"], hole["pi"], hole["sigma"]])
        p = np.array([particle["pi*"], particle["sigma*"]])
        fractions += pair_weight * np.outer(h, p)
        pairs.append({"rank": int(rank + 1), "weight": float(pair_weight),
                      "hole": hole, "particle": particle})
    flat = fractions.ravel(order="F")
    return {"fractions": dict(zip(_CHANNELS, map(float, flat))), "pairs": pairs,
            "classified": True, "reason": None}


def analyze_mrsf_transition(states, ao, n, fragments, ref=0, plane_normal=None,
                            hetero_atoms=None, weight_threshold=1.0e-4,
                            label_threshold=0.55, norm_tolerance=1.0e-10):
    """Analyze a physical MRSF root pair as LE/CT and orbital-character fractions.

    State indices and atom indices are 0-based.  ``fragments`` defines the
    chemical units used for LE/CT; atoms are not silently treated as fragments.
    For a planar chromophore the plane is inferred.  Supply ``plane_normal``
    for a local chromophore in a non-planar molecule.  If no unique plane can
    be found, LE/CT and NTO results are still returned while orbital character
    is marked ``unclassified`` and its fractions are NaN.  The same happens
    when the ``ref -> n`` transition density has no appreciable norm, since
    every fraction here is normalized by it.

    The returned label is ``mixed`` unless the largest of the six orbital
    channels is at least ``label_threshold``.  The fractions remain the primary
    result and should be reported alongside any compact label.

    .. warning::
       ``n`` is *heteroatom in-plane population*, not a localized lone pair:
       a heteroatom s orbital and a polarized C-O sigma bond both land in it.
       The n/sigma split is therefore indicative only; ``pi`` (and hence the
       pi*/sigma* axis) is the well-defined quantity, being a projection on the
       molecular-plane normal.  Treat ``n->pi*`` vs ``sigma->pi*`` as a hint and
       confirm with the NTO pairs in ``orbital_character["pairs"]``.
    """
    ref = int(ref)
    n = int(n)
    atoms = np.asarray(states.mol.get_atoms(), dtype=int)
    plane_reason = None
    if plane_normal is None:
        try:
            normal = infer_molecular_plane(ao, atoms)
        except ValueError as exc:
            normal = None
            plane_reason = str(exc)
    else:
        normal = _unit_vector(plane_normal, "plane_normal")
    if hetero_atoms is None:
        hetero_atoms = np.flatnonzero((atoms != 1) & (atoms != 6))
    hetero_atoms = np.asarray(hetero_atoms, dtype=int)

    ct = fragment_ct_matrix(states, ao, n, fragments, ref=ref,
                            norm_tolerance=norm_tolerance)
    nto = nto_transition(states, ref, n)
    if normal is None:
        orbital = {"fractions": {key: float("nan") for key in _CHANNELS},
                   "pairs": [], "classified": False, "reason": plane_reason}
    else:
        orbital = _orbital_channel_analysis(
            nto, states, ao, normal, hetero_atoms, float(weight_threshold),
            float(norm_tolerance))
    if orbital["classified"]:
        best = max(orbital["fractions"], key=orbital["fractions"].get)
        best_fraction = orbital["fractions"][best]
        orbital["label_ascii"] = best if best_fraction >= label_threshold else "mixed"
        orbital["label"] = _PRETTY[best] if best_fraction >= label_threshold else "mixed"
    else:
        orbital["label_ascii"] = "unclassified"
        orbital["label"] = "unclassified"
    orbital["plane_normal"] = normal
    orbital["hetero_atoms"] = hetero_atoms
    orbital["method"] = "Loewdin NTO p-perpendicular/heteroatom in-plane population"

    return {
        "pair": (ref, n),
        "energy_gap": float(states.energies[n] - states.energies[ref]),
        "fragment_ct": ct,
        "le_fraction": ct["le_fraction"],
        "ct_fraction": ct["ct_fraction"],
        "nto": nto,
        "nto_participation_ratio": participation_ratio(nto["weights"]),
        "orbital_character": orbital,
    }
