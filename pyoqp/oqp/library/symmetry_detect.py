"""Backend-free molecular point-group detection for symmetry metadata.

This module detects the full Schoenflies point group of a molecular geometry
and resolves the largest usable abelian subgroup from the D2h family
(C1, Ci, Cs, C2, C2v, C2h, D2, D2h) together with a standard orientation and
per-operation atom permutations.

Detection is geometry-only metadata: it does not change SCF/integral/response
execution behavior, and all symmetry reductions stay off.
"""

from __future__ import annotations

from typing import Any, MutableMapping

import numpy as np

# Operations of the D2h family expressed as sign matrices diag(sx, sy, sz)
# in the standard orientation.
SIGN_OPERATIONS: dict[str, tuple[int, int, int]] = {
    'E': (1, 1, 1),
    'C2z': (-1, -1, 1),
    'C2y': (-1, 1, -1),
    'C2x': (1, -1, -1),
    'i': (-1, -1, -1),
    'sxy': (1, 1, -1),
    'sxz': (1, -1, 1),
    'syz': (-1, 1, 1),
}

# Operation lists per abelian group (fixed order; character tables use the
# same op order).
ABELIAN_GROUP_OPS: dict[str, list[str]] = {
    'c1': ['E'],
    'ci': ['E', 'i'],
    'cs': ['E', 'sxy'],
    'c2': ['E', 'C2z'],
    'c2v': ['E', 'C2z', 'sxz', 'syz'],
    'c2h': ['E', 'C2z', 'i', 'sxy'],
    'd2': ['E', 'C2z', 'C2y', 'C2x'],
    'd2h': ['E', 'C2z', 'C2y', 'C2x', 'i', 'sxy', 'sxz', 'syz'],
}

# Mulliken-convention character tables, rows keyed by irrep label, columns
# following ABELIAN_GROUP_OPS order.
CHARACTER_TABLES: dict[str, dict[str, list[int]]] = {
    'c1': {'a': [1]},
    'ci': {'ag': [1, 1], 'au': [1, -1]},
    'cs': {"a'": [1, 1], "a''": [1, -1]},
    'c2': {'a': [1, 1], 'b': [1, -1]},
    'c2v': {
        'a1': [1, 1, 1, 1],
        'a2': [1, 1, -1, -1],
        'b1': [1, -1, 1, -1],
        'b2': [1, -1, -1, 1],
    },
    'c2h': {
        'ag': [1, 1, 1, 1],
        'bg': [1, -1, 1, -1],
        'au': [1, 1, -1, -1],
        'bu': [1, -1, -1, 1],
    },
    'd2': {
        'a': [1, 1, 1, 1],
        'b1': [1, 1, -1, -1],
        'b2': [1, -1, 1, -1],
        'b3': [1, -1, -1, 1],
    },
    'd2h': {
        'ag': [1, 1, 1, 1, 1, 1, 1, 1],
        'b1g': [1, 1, -1, -1, 1, 1, -1, -1],
        'b2g': [1, -1, 1, -1, 1, -1, 1, -1],
        'b3g': [1, -1, -1, 1, 1, -1, -1, 1],
        'au': [1, 1, 1, 1, -1, -1, -1, -1],
        'b1u': [1, 1, -1, -1, -1, -1, 1, 1],
        'b2u': [1, -1, 1, -1, -1, 1, -1, 1],
        'b3u': [1, -1, -1, 1, -1, 1, 1, -1],
    },
}

_MAX_PROPER_ORDER = 8


def _normalize(vec: np.ndarray) -> np.ndarray | None:
    norm = float(np.linalg.norm(vec))
    if norm < 1.0e-8:
        return None
    return vec / norm


def _rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Right-handed rotation by `angle` about unit vector `axis`."""
    x, y, z = axis
    c = np.cos(angle)
    s = np.sin(angle)
    cc = 1.0 - c
    return np.array([
        [c + x * x * cc, x * y * cc - z * s, x * z * cc + y * s],
        [y * x * cc + z * s, c + y * y * cc, y * z * cc - x * s],
        [z * x * cc - y * s, z * y * cc + x * s, c + z * z * cc],
    ])


def _reflection_matrix(normal: np.ndarray) -> np.ndarray:
    n = np.asarray(normal, dtype=float)
    return np.eye(3) - 2.0 * np.outer(n, n)


def _match_permutation(
    charges: np.ndarray,
    coords: np.ndarray,
    transformed: np.ndarray,
    tolerance: float,
) -> list[int] | None:
    """Return permutation p with transformed[i] == coords[p[i]] or None."""
    natom = coords.shape[0]
    used = np.zeros(natom, dtype=bool)
    permutation: list[int] = []
    for i in range(natom):
        dist = np.linalg.norm(coords - transformed[i], axis=1)
        dist[charges != charges[i]] = np.inf
        dist[used] = np.inf
        j = int(np.argmin(dist))
        if dist[j] > tolerance:
            return None
        used[j] = True
        permutation.append(j)
    return permutation


def _candidate_directions(charges: np.ndarray, coords: np.ndarray) -> list[np.ndarray]:
    """Axis/normal candidates: principal axes, atom vectors, same-element
    pair sums/differences, and pairwise cross products."""
    candidates: list[np.ndarray] = []

    inertia = _inertia_tensor(charges, coords)
    _, axes = np.linalg.eigh(inertia)
    for k in range(3):
        candidates.append(axes[:, k])

    vectors = [v for v in (coords[i] for i in range(coords.shape[0]))]
    for v in vectors:
        nv = _normalize(v)
        if nv is not None:
            candidates.append(nv)

    natom = coords.shape[0]
    for i in range(natom):
        for j in range(i + 1, natom):
            if charges[i] != charges[j]:
                continue
            for combo in (coords[i] + coords[j], coords[i] - coords[j]):
                nv = _normalize(combo)
                if nv is not None:
                    candidates.append(nv)

    for i in range(natom):
        for j in range(i + 1, natom):
            nv = _normalize(np.cross(coords[i], coords[j]))
            if nv is not None:
                candidates.append(nv)

    return _dedupe_directions(candidates)


def _dedupe_directions(candidates: list[np.ndarray]) -> list[np.ndarray]:
    unique: list[np.ndarray] = []
    for vec in candidates:
        if all(abs(float(np.dot(vec, u))) < 1.0 - 1.0e-6 for u in unique):
            unique.append(vec)
    return unique


def _inertia_tensor(charges: np.ndarray, coords: np.ndarray) -> np.ndarray:
    # Atomic numbers stand in for masses: symmetry-equivalent atoms share Z.
    r2 = np.sum(coords**2, axis=1)
    tensor = np.einsum('i,i,jk->jk', charges, r2, np.eye(3))
    tensor -= np.einsum('i,ij,ik->jk', charges, coords, coords)
    return tensor


def _is_symmetry_op(
    charges: np.ndarray,
    coords: np.ndarray,
    matrix: np.ndarray,
    tolerance: float,
) -> bool:
    transformed = coords @ matrix.T
    return _match_permutation(charges, coords, transformed, tolerance) is not None


def _proper_axis_order(
    charges: np.ndarray,
    coords: np.ndarray,
    axis: np.ndarray,
    tolerance: float,
) -> int:
    """Largest n <= _MAX_PROPER_ORDER with C_n about `axis` a symmetry op."""
    best = 1
    for order in range(2, _MAX_PROPER_ORDER + 1):
        rot = _rotation_matrix(axis, 2.0 * np.pi / order)
        if _is_symmetry_op(charges, coords, rot, tolerance):
            best = order
    return best


def _is_linear(coords: np.ndarray, tolerance: float) -> np.ndarray | None:
    """Return the molecular axis if all atoms are collinear, else None."""
    axis = None
    for i in range(coords.shape[0]):
        axis = _normalize(coords[i])
        if axis is not None:
            break
    if axis is None:  # all atoms at the origin-like degenerate case
        return np.array([0.0, 0.0, 1.0])
    for i in range(coords.shape[0]):
        residual = coords[i] - np.dot(coords[i], axis) * axis
        if np.linalg.norm(residual) > tolerance:
            return None
    return axis


def _perpendicular_vector(axis: np.ndarray) -> np.ndarray:
    trial = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(trial, axis))) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])
    perp = trial - np.dot(trial, axis) * axis
    result = _normalize(perp)
    assert result is not None
    return result


class _ElementSurvey:
    """Verified symmetry elements of a centered geometry."""

    def __init__(self, charges: np.ndarray, coords: np.ndarray, tolerance: float):
        self.charges = charges
        self.coords = coords
        self.tolerance = tolerance
        self.has_inversion = _is_symmetry_op(charges, coords, -np.eye(3), tolerance)
        self.proper_axes: list[tuple[np.ndarray, int]] = []
        self.mirror_normals: list[np.ndarray] = []

        for direction in _candidate_directions(charges, coords):
            order = _proper_axis_order(charges, coords, direction, tolerance)
            if order > 1:
                self.proper_axes.append((direction, order))
            if _is_symmetry_op(charges, coords, _reflection_matrix(direction), tolerance):
                self.mirror_normals.append(direction)

    @property
    def c2_axes(self) -> list[np.ndarray]:
        # A C2 about an axis exists only when its maximal order is even.
        return [axis for axis, order in self.proper_axes if order % 2 == 0]

    def principal(self) -> tuple[np.ndarray | None, int]:
        if not self.proper_axes:
            return None, 1
        axis, order = max(self.proper_axes, key=lambda item: item[1])
        return axis, order


def _schoenflies_from_survey(survey: _ElementSurvey) -> str:
    axis, order = survey.principal()

    high_order_axes = [o for _, o in survey.proper_axes if o >= 3]
    if len(high_order_axes) >= 2:
        # Cubic/icosahedral families.
        if any(o >= 5 for o in high_order_axes):
            return 'ih' if survey.has_inversion else 'i'
        if any(o == 4 for o in high_order_axes):
            return 'oh' if survey.has_inversion else 'o'
        if survey.has_inversion:
            return 'th'
        return 'td' if survey.mirror_normals else 't'

    if axis is None:
        if survey.mirror_normals:
            return 'cs'
        return 'ci' if survey.has_inversion else 'c1'

    tol_angle = 1.0e-6
    perpendicular_c2 = [
        a for a, o in survey.proper_axes
        if abs(float(np.dot(a, axis))) < 1.0e-3 and o % 2 == 0
    ]
    sigma_h = any(
        abs(float(np.dot(n, axis))) > 1.0 - tol_angle for n in survey.mirror_normals
    )
    sigma_v = [
        n for n in survey.mirror_normals if abs(float(np.dot(n, axis))) < 1.0e-3
    ]

    if len(perpendicular_c2) >= order:
        if sigma_h:
            return f'd{order}h'
        if len(sigma_v) >= order:
            return f'd{order}d'
        return f'd{order}'

    if sigma_h:
        return f'c{order}h'
    if len(sigma_v) >= order:
        return f'c{order}v'

    improper = _rotation_matrix(axis, np.pi / order) @ _reflection_matrix(axis)
    if _is_symmetry_op(survey.charges, survey.coords, improper, survey.tolerance):
        return f's{2 * order}'
    return f'c{order}'


def _orthonormal_frame(z: np.ndarray, x: np.ndarray | None = None) -> np.ndarray:
    """Rotation matrix with rows (x, y, z) forming a right-handed frame."""
    if x is None:
        x = _perpendicular_vector(z)
    x = x - np.dot(x, z) * z
    x_unit = _normalize(x)
    assert x_unit is not None
    y = np.cross(z, x_unit)
    return np.vstack([x_unit, y, z])


def _verify_group(
    charges: np.ndarray,
    coords: np.ndarray,
    rotation: np.ndarray,
    group: str,
    tolerance: float,
) -> list[dict[str, Any]] | None:
    """Verify all ops of `group` on the rotated geometry; return op payloads."""
    rotated = coords @ rotation.T
    operations: list[dict[str, Any]] = []
    for name in ABELIAN_GROUP_OPS[group]:
        matrix = np.diag(SIGN_OPERATIONS[name]).astype(float)
        permutation = _match_permutation(
            charges, rotated, rotated @ matrix.T, tolerance,
        )
        if permutation is None:
            return None
        operations.append({
            'name': name,
            'matrix': matrix.tolist(),
            'permutation': permutation,
        })
    return operations


def _abelian_resolution(
    survey: _ElementSurvey,
    principal_axis: np.ndarray | None,
) -> tuple[str, np.ndarray, list[dict[str, Any]]]:
    """Pick the largest verified D2h-family subgroup and standard orientation."""
    charges, coords, tol = survey.charges, survey.coords, survey.tolerance

    c2_axes = survey.c2_axes
    mirrors = survey.mirror_normals

    def axis_priority(axis: np.ndarray) -> float:
        if principal_axis is None:
            return 0.0
        return -abs(float(np.dot(axis, principal_axis)))

    candidates: list[tuple[str, np.ndarray]] = []

    # D2h / D2: three mutually perpendicular C2 axes.
    ordered_axes = sorted(c2_axes, key=axis_priority)
    for z_axis in ordered_axes:
        for x_axis in ordered_axes:
            if abs(float(np.dot(z_axis, x_axis))) > 1.0e-3:
                continue
            frame = _orthonormal_frame(z_axis, x_axis)
            candidates.append(('d2h', frame))
            candidates.append(('d2', frame))

    # C2h: C2 axis with parallel mirror normal and inversion.
    if survey.has_inversion:
        for z_axis in ordered_axes:
            if any(abs(float(np.dot(n, z_axis))) > 1.0 - 1.0e-6 for n in mirrors):
                candidates.append(('c2h', _orthonormal_frame(z_axis)))

    # C2v: C2 axis with two perpendicular mirror planes containing it.
    for z_axis in ordered_axes:
        containing = [n for n in mirrors if abs(float(np.dot(n, z_axis))) < 1.0e-3]
        for a in range(len(containing)):
            for b in range(len(containing)):
                if a == b:
                    continue
                if abs(float(np.dot(containing[a], containing[b]))) > 1.0e-3:
                    continue
                # sigma_xz has normal y, sigma_yz has normal x.
                frame = _orthonormal_frame(z_axis, containing[b])
                candidates.append(('c2v', frame))

    for z_axis in ordered_axes:
        candidates.append(('c2', _orthonormal_frame(z_axis)))

    for normal in mirrors:
        candidates.append(('cs', _orthonormal_frame(normal)))

    if survey.has_inversion:
        candidates.append(('ci', np.eye(3)))

    candidates.append(('c1', np.eye(3)))

    group_rank = {'d2h': 0, 'd2': 1, 'c2h': 2, 'c2v': 3, 'c2': 4, 'cs': 5, 'ci': 6, 'c1': 7}
    # Secondary key: prefer frames closest to the identity. This makes the
    # frame choice deterministic for geometries already in a valid standard
    # orientation (degenerate axis assignments, e.g. the three C2 axes of
    # d2h, would otherwise let repeated detections ping-pong between
    # equivalent frames and the reorientation iteration never converge).
    candidates.sort(key=lambda item: (
        group_rank[item[0]],
        float(np.abs(item[1] - np.eye(3)).max()),
    ))

    for group, rotation in candidates:
        operations = _verify_group(charges, coords, rotation, group, tol)
        if operations is not None:
            return group, rotation, operations

    # Unreachable: 'c1' always verifies.
    raise RuntimeError('abelian subgroup resolution failed')


def detect_point_group(
    atomic_numbers: Any,
    coordinates: Any,
    tolerance: float = 1.0e-5,
) -> dict[str, Any]:
    """Detect full point group and abelian (D2h-family) subgroup metadata.

    Parameters
    ----------
    atomic_numbers:
        Per-atom nuclear charges (used as symmetry-equivalence weights).
    coordinates:
        Cartesian coordinates, shape (natom, 3) (any consistent length unit).
    tolerance:
        Absolute displacement tolerance for accepting a symmetry operation.
    """
    if tolerance <= 0:
        raise ValueError('tolerance must be positive')

    charges = np.asarray(atomic_numbers, dtype=float).ravel()
    coords = np.asarray(coordinates, dtype=float).reshape(-1, 3).copy()
    if charges.size != coords.shape[0]:
        raise ValueError('atomic_numbers and coordinates disagree on atom count')
    if charges.size == 0:
        raise ValueError('geometry must contain at least one atom')

    center = np.einsum('i,ij->j', charges, coords) / float(np.sum(charges))
    coords -= center

    natom = charges.size

    if natom == 1:
        group = 'kh'
        rotation = np.eye(3)
        operations = _verify_group(charges, coords, rotation, 'd2h', tolerance)
        assert operations is not None
        subgroup = 'd2h'
    else:
        linear_axis = _is_linear(coords, tolerance)
        if linear_axis is not None:
            rotation = _orthonormal_frame(linear_axis)
            if _is_symmetry_op(charges, coords, -np.eye(3), tolerance):
                group, subgroup = 'dooh', 'd2h'
            else:
                group, subgroup = 'coov', 'c2v'
            operations = _verify_group(charges, coords, rotation, subgroup, tolerance)
            assert operations is not None
        else:
            survey = _ElementSurvey(charges, coords, tolerance)
            group = _schoenflies_from_survey(survey)
            principal_axis, _ = survey.principal()
            subgroup, rotation, operations = _abelian_resolution(survey, principal_axis)

    # Operation matrices conjugated back to the input frame: AO/MO data
    # produced by the backend lives in input coordinates, not the standard
    # orientation, so O_input = R^T O_standard R (atom permutations are
    # frame-independent).
    for op in operations:
        standard = np.asarray(op['matrix'])
        op['matrix_input_frame'] = (rotation.T @ standard @ rotation).tolist()

    return {
        'point_group': group,
        'abelian_subgroup': subgroup,
        'origin': center.tolist(),
        'orientation': rotation.tolist(),
        'operations': operations,
        'irreps': list(CHARACTER_TABLES[subgroup].keys()),
        'character_table': {
            irrep: list(row) for irrep, row in CHARACTER_TABLES[subgroup].items()
        },
        'tolerance': float(tolerance),
        'n_atoms': int(natom),
    }


def attach_detection_metadata(
    symmetry_metadata: MutableMapping[str, Any],
    atomic_numbers: Any,
    coordinates: Any,
) -> MutableMapping[str, Any]:
    """Run detection and record results in `symmetry_metadata` in place.

    Metadata-only: resolves 'auto' point-group/subgroup requests against the
    detected values and never enables symmetry reductions.
    """
    tolerance = float(symmetry_metadata.get('tolerance', 1.0e-5))
    detection = detect_point_group(atomic_numbers, coordinates, tolerance=tolerance)

    symmetry_metadata['detected_point_group'] = detection['point_group']
    symmetry_metadata['detected_subgroup'] = detection['abelian_subgroup']
    symmetry_metadata['detection'] = detection

    requested_group = str(symmetry_metadata.get('requested_point_group', 'auto')).lower()
    requested_subgroup = str(symmetry_metadata.get('requested_subgroup', 'auto')).lower()

    if requested_group == 'auto':
        symmetry_metadata['point_group'] = detection['point_group']
    if requested_subgroup == 'auto':
        symmetry_metadata['subgroup'] = detection['abelian_subgroup']

    symmetry_metadata['requested_matches_detected'] = (
        requested_group in ('auto', detection['point_group'])
        and requested_subgroup in ('auto', detection['abelian_subgroup'])
    )

    # Reduction flags are policy-controlled by the input checker (integral:
    # experimental opt-in; response: rejected); detection only defaults them.
    symmetry_metadata.setdefault('use_integral_symmetry', False)
    symmetry_metadata.setdefault('use_response_symmetry', False)

    return symmetry_metadata
