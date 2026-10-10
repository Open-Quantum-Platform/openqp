"""Tolerance-aware molecular point-group analysis and principal-axis alignment."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from math import cos, gcd, isfinite, pi, sin

import numpy as np

from .molden import SYMBOLS
from .structure_io import Atom, parse_xyz

_MASSES = [
    0.0, 1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999,
    18.998, 20.180, 22.990, 24.305, 26.982, 28.085, 30.974, 32.06,
    35.45, 39.948, 39.0983, 40.078, 44.956, 47.867, 50.942, 51.996,
    54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630,
    74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224,
    92.906, 95.95, 98.0, 101.07, 102.906, 106.42, 107.868, 112.414,
    114.818, 118.710, 121.760, 127.60, 126.904, 131.293, 132.905, 137.327,
    138.905, 140.116, 140.908, 144.242, 145.0, 150.36, 151.964, 157.25,
    158.925, 162.500, 164.930, 167.259, 168.934, 173.045, 174.967, 178.49,
    180.948, 183.84, 186.207, 190.23, 192.217, 195.084, 196.967, 200.592,
    204.38, 207.2, 208.980, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0,
    232.038, 231.036, 238.029, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0,
    252.0, 257.0, 258.0, 259.0, 266.0, 267.0, 268.0, 269.0, 270.0,
    277.0, 278.0, 281.0, 282.0, 285.0, 286.0, 289.0, 290.0, 293.0,
    294.0, 294.0,
]
ATOMIC_MASS = dict(zip(SYMBOLS, _MASSES))
MAX_MATCH_WORK = 100_000_000
MAX_ASSIGNMENT_EDGES = 10_000
MAX_ASSIGNMENT_VISITS = 100_000
MAX_XYZ_CHARACTERS = 1_000_000
MAX_ROTATION_SCREEN_WORK = 5_000_000


@dataclass
class Operation:
    label: str
    matrix: np.ndarray
    mapping: list[int]
    residual: float


def _unit(vector: np.ndarray) -> np.ndarray | None:
    length = float(np.linalg.norm(vector))
    return vector / length if length > 1.0e-10 else None


def _unique_axes(vectors: list[np.ndarray]) -> list[np.ndarray]:
    axes: list[np.ndarray] = []
    for vector in vectors:
        axis = _unit(vector)
        if axis is None or any(abs(float(np.dot(axis, other))) > 1 - 1.0e-6 for other in axes):
            continue
        axes.append(axis)
    return axes


def _rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    x, y, z = axis
    skew = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) * cos(angle) + (1 - cos(angle)) * np.outer(axis, axis) + sin(angle) * skew


def _rotation_orders(symbols: list[str], coordinates: np.ndarray, axis: np.ndarray,
                     tolerance: float) -> list[int]:
    """Orders compatible with atoms detectably moved by each candidate rotation."""
    return _operation_orders(symbols, coordinates, axis, tolerance, improper=False)


def _improper_orders(symbols: list[str], coordinates: np.ndarray, axis: np.ndarray,
                     tolerance: float) -> list[int]:
    """Orders compatible with displacement from rotation followed by reflection."""
    return _operation_orders(symbols, coordinates, axis, tolerance, improper=True)


def _operation_orders(symbols: list[str], coordinates: np.ndarray, axis: np.ndarray,
                      tolerance: float, *, improper: bool) -> list[int]:
    maximum = max(Counter(symbols).values(), default=1)
    reflection = np.eye(3) - 2 * np.outer(axis, axis)
    orders: list[int] = []
    for order in range(2, maximum + 1):
        matrix = _rotation(axis, 2 * pi / order)
        if improper:
            matrix = matrix @ reflection
        displacement = np.linalg.norm(coordinates @ matrix.T - coordinates, axis=1)
        counts = Counter(
            symbol for symbol, moved in zip(symbols, displacement) if moved > tolerance
        )
        multiplicity = 0
        for count in counts.values():
            multiplicity = gcd(multiplicity, count)
        if multiplicity and multiplicity % order == 0:
            orders.append(order)
        elif improper and counts:
            period = order if order % 2 == 0 else 2 * order
            divisors = [value for value in range(2, period + 1) if period % value == 0]

            def representable(count: int, allowed: tuple[int, ...] = tuple(divisors)) -> bool:
                reachable = [True] + [False] * count
                for total in range(1, count + 1):
                    reachable[total] = any(
                        divisor <= total and reachable[total - divisor] for divisor in allowed
                    )
                return reachable[count]

            if (all(representable(count) for count in counts.values())
                    and any(count >= order and representable(count - order)
                            for count in counts.values())):
                orders.append(order)
    return orders


def _distinct_rotation_axes(entries: list[tuple[np.ndarray, Operation]],
                            coordinates: np.ndarray, tolerance: float) -> list[np.ndarray]:
    """Cluster axes that cannot be distinguished at the coordinate tolerance."""
    radius = float(np.max(np.linalg.norm(coordinates, axis=1))) if len(coordinates) else 0.0
    if radius <= 1.0e-12:
        angular_cosine = -1.0
    else:
        # atan(tolerance / radius) remains meaningful when tolerance exceeds
        # the molecular radius, without collapsing all directions together.
        angular_cosine = radius / (radius ** 2 + tolerance ** 2) ** 0.5
    axes: list[np.ndarray] = []
    for axis, _operation in entries:
        if any(abs(float(np.dot(axis, other))) >= angular_cosine for other in axes):
            continue
        axes.append(axis)
    return axes


def _face_normal_axes(coordinates: np.ndarray) -> list[np.ndarray]:
    """Candidate axes through polygon centers, ranked by repeated local support."""
    if len(coordinates) < 4:
        return []
    clusters: dict[tuple[float, float, float], tuple[np.ndarray, int]] = {}
    neighbor_count = min(7, len(coordinates) - 1)
    for index, point in enumerate(coordinates):
        distance = np.linalg.norm(coordinates - point, axis=1)
        neighbors = [value for value in np.argsort(distance) if value != index][:neighbor_count]
        for first, second in combinations(neighbors, 2):
            axis = _unit(np.cross(coordinates[first] - point, coordinates[second] - point))
            if axis is None:
                continue
            first_nonzero = next((value for value in axis if abs(value) > 1.0e-10), 1.0)
            canonical = axis if first_nonzero > 0 else -axis
            key = tuple(float(round(value, 6)) for value in canonical)
            existing, support = clusters.get(key, (canonical, 0))
            clusters[key] = (existing, support + 1)
    ranked = sorted(clusters.items(), key=lambda item: (-item[1][1], item[0]))
    return [axis for _key, (axis, _support) in ranked[:240]]


def _assignment(distances: np.ndarray, tolerance: float) -> list[int] | None:
    """Find a same-element one-to-one assignment within the stated tolerance."""
    choices = [list(np.flatnonzero(row <= tolerance)) for row in distances]
    if sum(map(len, choices)) > MAX_ASSIGNMENT_EDGES:
        raise ValueError("symmetry assignment exceeded its work limit; use a tighter tolerance")
    order = sorted(range(len(choices)), key=lambda index: len(choices[index]))
    if any(not choices[index] for index in order):
        return None
    mapping = [-1] * len(choices)
    owner = [-1] * len(choices)
    visits = 0

    def augment(source: int, seen: set[int]) -> bool:
        nonlocal visits
        for target in sorted(choices[source], key=lambda index: distances[source, index]):
            visits += 1
            if visits > MAX_ASSIGNMENT_VISITS:
                raise ValueError("symmetry assignment exceeded its work limit; use a tighter tolerance")
            if target in seen:
                continue
            seen.add(target)
            if owner[target] < 0 or augment(owner[target], seen):
                owner[target] = source
                mapping[source] = target
                return True
        return False

    return mapping if all(augment(source, set()) for source in order) else None


def _match(symbols: list[str], coordinates: np.ndarray, matrix: np.ndarray,
           tolerance: float) -> tuple[list[int], float] | None:
    transformed = coordinates @ matrix.T
    mapping = [-1] * len(symbols)
    residual = 0.0
    for symbol in sorted(set(symbols)):
        indices = [index for index, value in enumerate(symbols) if value == symbol]
        source = transformed[indices]
        target = coordinates[indices]
        distances = np.linalg.norm(source[:, None, :] - target[None, :, :], axis=2)
        local = _assignment(distances, tolerance)
        if local is None:
            return None
        for row, column in enumerate(local):
            mapping[indices[row]] = indices[column]
            residual = max(residual, float(distances[row, column]))
    return mapping, residual


def _moves_coordinates(coordinates: np.ndarray, matrix: np.ndarray, tolerance: float) -> bool:
    """A numerically indistinguishable rotation is not a symmetry operation."""
    displacement = np.linalg.norm(coordinates @ matrix.T - coordinates, axis=1)
    return bool(len(displacement) and float(np.max(displacement)) > tolerance)


def _operation_powers(label: str, matrix: np.ndarray, matched: tuple[list[int], float],
                      order: int, coordinates: np.ndarray,
                      tolerance: float) -> list[Operation] | None:
    """Close a matched generator by composing its atom mapping."""
    generator_mapping, generator_residual = matched
    result = [Operation(label, matrix, generator_mapping, generator_residual)]
    mapping = generator_mapping
    for power in range(2, order):
        mapping = [generator_mapping[target] for target in mapping]
        powered = np.linalg.matrix_power(matrix, power)
        transformed = coordinates @ powered.T
        residual = max(
            float(np.linalg.norm(transformed[source] - coordinates[target]))
            for source, target in enumerate(mapping)
        )
        if residual > tolerance:
            return None
        result.append(Operation(f"{label}^{power}", powered, mapping.copy(), residual))
    return result


def _matrix_key(matrix: np.ndarray) -> bytes:
    rounded = np.round(matrix, 7)
    rounded[np.abs(rounded) < 1.0e-7] = 0.0
    return rounded.tobytes()


def _strict_xyz(xyz: str) -> list[Atom]:
    if len(xyz) > MAX_XYZ_CHARACTERS:
        raise ValueError(f"symmetry input is limited to {MAX_XYZ_CHARACTERS:,} characters")
    lines = xyz.splitlines()
    first = next((index for index, line in enumerate(lines) if line.strip()), len(lines))
    lines = lines[first:]
    while lines and not lines[-1].strip():
        lines.pop()
    if not lines:
        raise ValueError("no Cartesian coordinates were found")
    try:
        declared = int(lines[0].strip())
    except ValueError:
        declared = None
    if declared is not None:
        if declared > 300:
            raise ValueError("symmetry analysis supports at most 300 atoms")
        if declared < 1 or len(lines) < declared + 2:
            raise ValueError("XYZ atom count or coordinate rows are invalid")
        if any(line.strip() for line in lines[declared + 2:]):
            raise ValueError("symmetry analysis accepts one XYZ structure at a time")
        source = "\n".join(lines[:declared + 2])
        coordinate_rows = lines[2:declared + 2]
        frames = parse_xyz(source)
        expected = declared
    else:
        rows = [line for line in lines if line.strip()]
        if len(rows) > 300:
            raise ValueError("symmetry analysis supports at most 300 atoms")
        frames = parse_xyz("\n".join(rows))
        coordinate_rows = rows
        expected = len(rows)
    for row in coordinate_rows:
        fields = row.split(maxsplit=4)
        if len(fields) != 4 or fields[0] not in ATOMIC_MASS:
            raise ValueError("every coordinate row must use an exact element symbol and three numbers")
    if len(frames) != 1 or len(frames[0].atoms) != expected:
        raise ValueError("every coordinate row must contain a valid element and three numbers")
    atoms = frames[0].atoms
    if any(not isfinite(value) for atom in atoms for value in atom[1:]):
        raise ValueError("coordinates must be finite numbers")
    return atoms


def _principal_axes(atoms: list[Atom]) -> tuple[list[str], np.ndarray, np.ndarray]:
    symbols = [atom[0] for atom in atoms]
    coordinates = np.asarray([atom[1:] for atom in atoms], dtype=float)
    weights = np.asarray([ATOMIC_MASS[symbol] for symbol in symbols])
    if float(np.sum(weights)) <= 0:
        raise ValueError("symmetry analysis requires at least one atom with positive mass")
    center = np.average(coordinates, axis=0, weights=weights)
    centered = coordinates - center
    inertia = sum(
        weight * (float(np.dot(point, point)) * np.eye(3) - np.outer(point, point))
        for weight, point in zip(weights, centered)
    )
    _, axes = np.linalg.eigh(inertia)
    if np.linalg.det(axes) < 0:
        axes[:, -1] *= -1
    return symbols, centered @ axes, center


def _equivalent_atoms(size: int, operations: list[Operation]) -> list[list[int]]:
    parent = list(range(size))

    def root(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for operation in operations:
        for source, target in enumerate(operation.mapping):
            a, b = root(source), root(target)
            if a != b:
                parent[b] = a
    groups: dict[int, list[int]] = {}
    for index in range(size):
        groups.setdefault(root(index), []).append(index + 1)
    return sorted(groups.values(), key=lambda group: group[0])


def analyze(xyz: str, tolerance: float = 0.05) -> dict:
    atoms = _strict_xyz(xyz)
    if len(atoms) > 300:
        raise ValueError("symmetry analysis supports at most 300 atoms")
    symbols, coordinates, center = _principal_axes(atoms)
    single_atom = len(atoms) == 1
    if len(atoms) > 2:
        _left, _spread, directions = np.linalg.svd(coordinates, full_matrices=False)
        principal_line = directions[0]
        projection = np.outer(coordinates @ principal_line, principal_line)
        maximum_perpendicular = float(np.max(np.linalg.norm(coordinates - projection, axis=1)))
        axial_extent = float(np.max(np.abs(coordinates @ principal_line)))
    else:
        maximum_perpendicular = 0.0
        axial_extent = 0.0
    linear = len(atoms) == 2 or (
        len(atoms) > 2 and maximum_perpendicular <= tolerance
        and maximum_perpendicular <= 0.05 * max(axial_extent, 1.0e-12)
    )

    seed_axes = [np.eye(3)[index] for index in range(3)]
    atom_axes = _unique_axes([point for point in coordinates])[:80]
    basis_axes = _unique_axes([*seed_axes, *atom_axes])
    cross_axes = _unique_axes([
        np.cross(a, b) for a, b in combinations(basis_axes[:43], 2)
    ])[:160]
    bisector_axes = _unique_axes([
        vector
        for a, b in combinations(atom_axes[:40], 2)
        for vector in (a + b, a - b)
    ])[:160]
    face_axes = _face_normal_axes(coordinates)
    axes = _unique_axes([*basis_axes, *face_axes, *cross_axes, *bisector_axes])

    match_tests = 0
    max_match_tests = max(100, MAX_MATCH_WORK // (len(atoms) ** 2))

    def bounded_match(matrix: np.ndarray) -> tuple[list[int], float] | None:
        nonlocal match_tests
        if match_tests >= max_match_tests:
            raise ValueError(
                "symmetry search exceeded its work limit; use fewer atoms or a tighter tolerance"
            )
        match_tests += 1
        return _match(symbols, coordinates, matrix, tolerance)

    identity = Operation("E", np.eye(3), list(range(len(atoms))), 0.0)
    operations = [identity]
    inversion_match = bounded_match(-np.eye(3))
    inversion = inversion_match is not None
    if inversion_match:
        operations.append(Operation("i", -np.eye(3), *inversion_match))

    rotations: dict[int, list[tuple[np.ndarray, Operation]]] = {}
    improper: dict[int, list[tuple[np.ndarray, Operation]]] = {}
    rotation_screen_work = 0
    maximum_order = max(Counter(symbols).values(), default=1)
    for axis in axes:
        rotation_screen_work += 2 * len(atoms) * max(0, maximum_order - 1)
        if rotation_screen_work > MAX_ROTATION_SCREEN_WORK:
            raise ValueError(
                "symmetry rotation screening exceeded its work limit; "
                "use fewer atoms or a tighter tolerance"
            )
        proper_orders = set(_rotation_orders(symbols, coordinates, axis, tolerance))
        improper_orders = set(_improper_orders(symbols, coordinates, axis, tolerance))
        reflection = np.eye(3) - 2 * np.outer(axis, axis)
        for order in sorted(proper_orders | improper_orders):
            matrix = _rotation(axis, 2 * pi / order)
            matched = (bounded_match(matrix) if order in proper_orders
                       and _moves_coordinates(coordinates, matrix, tolerance) else None)
            if matched:
                powers = _operation_powers(
                    f"C{order}", matrix, matched, order, coordinates, tolerance,
                )
                if powers:
                    rotations.setdefault(order, []).append((axis, powers[0]))
                    operations.extend(powers)
            improper_matrix = matrix @ reflection
            improper_match = (bounded_match(improper_matrix) if order in improper_orders
                              and _moves_coordinates(coordinates, improper_matrix, tolerance)
                              else None)
            if improper_match:
                powers = _operation_powers(
                    f"S{order}", improper_matrix, improper_match,
                    2 * order if order % 2 else order, coordinates, tolerance,
                )
                if powers:
                    improper.setdefault(order, []).append((axis, powers[0]))
                    operations.extend(powers)

    mirrors: list[tuple[np.ndarray, Operation]] = []
    for normal in axes:
        matrix = np.eye(3) - 2 * np.outer(normal, normal)
        matched = bounded_match(matrix)
        if matched:
            operation = Operation("sigma", matrix, *matched)
            mirrors.append((normal, operation))
            operations.append(operation)

    if single_atom:
        point_group = "Kh"
    elif linear:
        point_group = "Dinfh" if inversion else "Cinfv"
    elif len(_distinct_rotation_axes(rotations.get(5, []), coordinates, tolerance)) >= 6:
        point_group = "Ih" if inversion else "I"
    elif len(_distinct_rotation_axes(rotations.get(4, []), coordinates, tolerance)) >= 3:
        point_group = "Oh" if inversion else "O"
    elif (len(_distinct_rotation_axes(rotations.get(3, []), coordinates, tolerance)) >= 4
          and len(_distinct_rotation_axes(rotations.get(2, []), coordinates, tolerance)) >= 3):
        point_group = "Th" if inversion else "Td" if mirrors else "T"
    else:
        principal_order = max(rotations, default=1)
        if principal_order == 1:
            point_group = "Ci" if inversion else "Cs" if mirrors else "C1"
        else:
            improper_order = max(improper, default=1)
            principal_axis = (
                improper[improper_order][0][0]
                if improper_order > principal_order else rotations[principal_order][0][0]
            )
            perpendicular_c2 = any(
                abs(float(np.dot(axis, principal_axis))) < 0.2
                for axis, _operation in rotations.get(2, [])
            )
            horizontal = any(
                abs(float(np.dot(normal, principal_axis))) > 0.98
                for normal, _operation in mirrors
            )
            vertical = any(
                abs(float(np.dot(normal, principal_axis))) < 0.2
                for normal, _operation in mirrors
            )
            family = "D" if perpendicular_c2 else "C"
            if (not perpendicular_c2 and not horizontal and not vertical
                    and improper_order > principal_order):
                point_group = f"S{improper_order}"
            else:
                suffix = (
                    "h" if horizontal else "d" if perpendicular_c2 and vertical
                    else "v" if vertical else ""
                )
                point_group = f"{family}{principal_order}{suffix}"

    unique_operations: dict[bytes, Operation] = {}
    for operation in operations:
        unique_operations.setdefault(_matrix_key(operation.matrix), operation)
    accepted = list(unique_operations.values())
    aligned = [
        [symbol, float(x), float(y), float(z)]
        for symbol, (x, y, z) in zip(symbols, coordinates)
    ]
    return {
        "point_group": point_group,
        "tolerance_angstrom": tolerance,
        "center_angstrom": center.tolist(),
        "max_deviation_angstrom": max(operation.residual for operation in accepted),
        "operations": sorted({operation.label for operation in accepted}),
        "operation_count": len(accepted),
        "equivalent_atoms": _equivalent_atoms(len(atoms), accepted),
        "aligned_atoms": aligned,
    }
