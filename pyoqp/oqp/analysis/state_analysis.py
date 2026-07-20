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


def _ao_pi_projection(ao, normal):
    projection = np.zeros(ao.nbf)
    for mu, (_shell, powers, _scale) in enumerate(ao.ao_index):
        if sum(powers) == 1:
            axis = int(np.argmax(powers))
            projection[mu] = normal[axis] ** 2
    return projection


def _nto_orbital_character(coefficients, shalf, ao, normal, hetero_atoms):
    # q = S^(1/2)c are coefficients in the symmetric orthonormal AO basis.
    q = shalf @ np.asarray(coefficients, dtype=float)
    population = q * q
    norm = population.sum()
    if norm <= 1.0e-20:
        return {"n": 0.0, "pi": 0.0, "sigma": 0.0}
    population /= norm
    pi_projection = _ao_pi_projection(ao, normal)
    pi = float(np.dot(population, pi_projection))
    hetero_mask = np.isin(ao.ao_atom, np.asarray(hetero_atoms, dtype=int))
    in_plane = 1.0 - pi_projection
    lone_pair = float(np.dot(population * hetero_mask, in_plane))
    n = min(max(lone_pair, 0.0), max(1.0 - pi, 0.0))
    sigma = max(1.0 - pi - n, 0.0)
    total = n + pi + sigma
    return {"n": n / total, "pi": pi / total, "sigma": sigma / total}


def _orbital_channel_analysis(nto, states, ao, normal, hetero_atoms, weight_threshold):
    weights = np.asarray(nto["weights"], dtype=float)
    if weights.size == 0 or weights.sum() <= 0.0:
        return {"fractions": {key: 0.0 for key in _CHANNELS}, "pairs": []}
    keep = weights > weight_threshold * weights.max()
    if not np.any(keep):
        keep[np.argmax(weights)] = True
    selected = np.flatnonzero(keep)
    selected_weights = weights[selected]
    selected_weights /= selected_weights.sum()
    shalf = _lowdin_sqrt(states.S)
    fractions = np.zeros((3, 2))  # hole n/pi/sigma x particle pi*/sigma*
    pairs = []
    for rank, pair_weight in zip(selected, selected_weights):
        hole = _nto_orbital_character(
            nto["holes_ao"][:, rank], shalf, ao, normal, hetero_atoms)
        particle_hps = _nto_orbital_character(
            nto["particles_ao"][:, rank], shalf, ao, normal, ())
        particle = {"pi*": particle_hps["pi"],
                    "sigma*": particle_hps["n"] + particle_hps["sigma"]}
        h = np.array([hole["n"], hole["pi"], hole["sigma"]])
        p = np.array([particle["pi*"], particle["sigma*"]])
        fractions += pair_weight * np.outer(h, p)
        pairs.append({"rank": int(rank + 1), "weight": float(pair_weight),
                      "hole": hole, "particle": particle})
    flat = fractions.ravel(order="F")
    return {"fractions": dict(zip(_CHANNELS, map(float, flat))), "pairs": pairs}


def analyze_mrsf_transition(states, ao, n, fragments, ref=0, plane_normal=None,
                            hetero_atoms=None, weight_threshold=1.0e-4,
                            label_threshold=0.55):
    """Analyze a physical MRSF root pair as LE/CT and orbital-character fractions.

    State indices and atom indices are 0-based.  ``fragments`` defines the
    chemical units used for LE/CT; atoms are not silently treated as fragments.
    For a planar chromophore the plane is inferred.  Supply ``plane_normal``
    for a local chromophore in a non-planar molecule.  If no unique plane can
    be found, LE/CT and NTO results are still returned while orbital character
    is marked ``unclassified``.

    The returned label is ``mixed`` unless the largest of the six orbital
    channels is at least ``label_threshold``.  The fractions remain the primary
    result and should be reported alongside any compact label.
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

    ct = fragment_ct_matrix(states, ao, n, fragments, ref=ref)
    nto = nto_transition(states, ref, n)
    if normal is None:
        orbital = {"fractions": {key: 0.0 for key in _CHANNELS}, "pairs": []}
        orbital["label_ascii"] = "unclassified"
        orbital["label"] = "unclassified"
        orbital["classified"] = False
        orbital["reason"] = plane_reason
    else:
        orbital = _orbital_channel_analysis(
            nto, states, ao, normal, hetero_atoms, float(weight_threshold))
        best = max(orbital["fractions"], key=orbital["fractions"].get)
        best_fraction = orbital["fractions"][best]
        orbital["label_ascii"] = best if best_fraction >= label_threshold else "mixed"
        orbital["label"] = _PRETTY[best] if best_fraction >= label_threshold else "mixed"
        orbital["classified"] = True
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
