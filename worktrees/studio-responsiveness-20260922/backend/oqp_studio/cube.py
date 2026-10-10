"""Gaussian cube validation and pointwise arithmetic."""

from __future__ import annotations

import re
from dataclasses import dataclass
from io import StringIO, TextIOBase
from math import isfinite, prod

from .molden import BOHR_TO_ANGSTROM, SYMBOLS

MAX_GRID_VALUES = 2_000_000
MAX_CUBE_DATASETS = 10_000
MAX_CUBE_ATOMS = 10_000
MAX_HEADER_LINE = 16 * 1024
MAX_CUBE_HEADER = 4 * 1024 * 1024


@dataclass
class Cube:
    header: list[str]
    values: list[float]
    shape: tuple[int, int, int]
    datasets: int
    geometry: tuple[tuple[float, ...], ...]
    dataset_ids: tuple[int, ...]
    axis_units: tuple[bool, bool, bool]


def _number(token: str) -> float:
    value = float(token.replace("D", "E").replace("d", "e"))
    if not isfinite(value):
        raise ValueError("cube contains a non-finite number")
    return value


def _tokens(line: str):
    """Yield fields without allocating a list for a potentially huge line."""
    return (match.group() for match in re.finditer(r"\S+", line))


def _record(line: str, lengths: tuple[int, ...]) -> list[str]:
    """Read a fixed-width record without expanding all fields from malformed input."""
    maximum = max(lengths)
    values: list[str] = []
    for token in _tokens(line):
        if len(values) >= maximum:
            raise ValueError("cube record has too many fields")
        values.append(token)
    if len(values) not in lengths:
        raise ValueError("cube record has the wrong number of fields")
    return values


def parse(text: str, *, header_only: bool = False) -> Cube:
    stream = StringIO(text)
    lines: list[str] = []
    header_size = 0

    def header_line() -> str:
        nonlocal header_size
        line = stream.readline()
        if not line:
            raise ValueError("cube file is incomplete")
        header_size += len(line)
        if header_size > MAX_CUBE_HEADER:
            raise ValueError(
                f"cube header is limited to {MAX_CUBE_HEADER // (1024 * 1024)} MiB"
            )
        value = line.rstrip("\r\n")
        lines.append(value)
        return value

    try:
        for _index in range(6):
            header_line()
    except ValueError as exc:
        if "incomplete" in str(exc):
            raise ValueError("cube file is incomplete") from exc
        raise
    if len(lines) < 6:
        raise ValueError("cube file is incomplete")
    try:
        origin_record = _record(lines[2], (4, 5))
        atom_count = int(origin_record[0])
        atoms = abs(atom_count)
        if atoms > MAX_CUBE_ATOMS:
            raise ValueError
        axis_records = [_record(lines[index], (4,)) for index in (3, 4, 5)]
        voxel_counts = tuple(int(record[0]) for record in axis_records)
        shape = tuple(abs(count) for count in voxel_counts)
        if any(size == 0 for size in shape):
            raise ValueError
        if len({count < 0 for count in voxel_counts}) != 1:
            raise ValueError
    except (IndexError, ValueError) as exc:
        raise ValueError("cube header is invalid") from exc
    try:
        atom_records = [_record(header_line(), (5,)) for _index in range(atoms)]
    except ValueError as exc:
        if "header is limited" in str(exc):
            raise
        if "incomplete" in str(exc):
            raise ValueError("cube atom header is incomplete") from exc
        raise ValueError("cube atom header is invalid") from exc
    if any(len(record) != 5 for record in atom_records):
        raise ValueError("cube atom header is invalid")
    datasets = int(origin_record[4]) if len(origin_record) > 4 else 1
    if datasets < 1:
        raise ValueError("cube dataset count must be positive")
    voxel_points = prod(shape)
    if not header_only and datasets > MAX_CUBE_DATASETS:
        raise ValueError("cube dataset count exceeds the supported grid limit")
    if not header_only and voxel_points * datasets > MAX_GRID_VALUES:
        raise ValueError(f"cube grid is limited to {MAX_GRID_VALUES:,} values")
    dataset_id_values: tuple[int, ...] = ()
    if atom_count < 0:
        dataset_ids: int | None = None
        identifiers: list[int] = []
        while True:
            try:
                identifier_line = header_line()
            except ValueError as exc:
                if "incomplete" in str(exc):
                    raise ValueError("cube dataset identifier record is incomplete") from exc
                raise
            tokens = _tokens(identifier_line)
            for token in tokens:
                if dataset_ids is None:
                    try:
                        dataset_ids = int(token)
                    except ValueError as exc:
                        raise ValueError("cube dataset identifiers are invalid") from exc
                    maximum = min(MAX_CUBE_DATASETS, MAX_GRID_VALUES // voxel_points)
                    if dataset_ids < 1 or dataset_ids > maximum:
                        raise ValueError(
                            "cube dataset identifier count exceeds the supported grid limit"
                        )
                    continue
                if len(identifiers) >= dataset_ids:
                    raise ValueError("cube dataset identifier record is invalid")
                try:
                    identifiers.append(int(token))
                except ValueError as exc:
                    raise ValueError("cube dataset identifiers are invalid") from exc
            if dataset_ids is not None and len(identifiers) == dataset_ids:
                break
        if dataset_ids is None or len(identifiers) < dataset_ids:
            raise ValueError("cube dataset identifier record is incomplete")
        dataset_id_values = tuple(identifiers)
        if datasets != 1 and datasets != dataset_ids:
            raise ValueError("cube dataset counts disagree")
        datasets = dataset_ids
    values: list[float] = []
    if not header_only:
        expected = voxel_points * datasets
        if expected > MAX_GRID_VALUES:
            raise ValueError(f"cube grid is limited to {MAX_GRID_VALUES:,} values")
        try:
            for line in stream:
                for token in _tokens(line):
                    if len(values) >= expected:
                        raise ValueError(
                            f"cube grid has more than the expected {expected} values"
                        )
                    values.append(_number(token))
        except ValueError as exc:
            if "non-finite" in str(exc) or "more than the expected" in str(exc):
                raise
            raise ValueError("cube grid contains a nonnumeric value") from exc
        if len(values) != expected:
            raise ValueError(f"cube grid has {len(values)} values; expected {expected}")
    try:
        geometry = (
            tuple(_number(token) for token in origin_record[1:4]),
            *(tuple(_number(token) for token in record[1:4]) for record in axis_records),
            *(tuple(_number(token) for token in record) for record in atom_records),
        )
    except ValueError as exc:
        if "non-finite" in str(exc):
            raise
        raise ValueError("cube geometry header is invalid") from exc
    if any(not record[0].is_integer() for record in geometry[4:]):
        raise ValueError("cube atom header contains a nonintegral atomic number")
    axis_units = tuple(count < 0 for count in voxel_counts)
    return Cube(
        lines, values, shape, datasets, geometry, dataset_id_values, axis_units,
    )


def parse_header(stream: TextIOBase) -> Cube:
    """Read only the bounded Gaussian cube header from an open text stream."""
    total = 0

    def header_line() -> str:
        nonlocal total
        line = stream.readline(MAX_HEADER_LINE + 1)
        if len(line) > MAX_HEADER_LINE:
            raise ValueError(f"cube header lines are limited to {MAX_HEADER_LINE:,} characters")
        total += len(line)
        if total > MAX_CUBE_HEADER:
            raise ValueError(f"cube header is limited to {MAX_CUBE_HEADER // (1024 * 1024)} MiB")
        return line

    lines = [header_line() for _ in range(6)]
    if any(line == "" for line in lines):
        raise ValueError("cube file is incomplete")
    try:
        atom_count = int(_record(lines[2], (4, 5))[0])
    except (IndexError, ValueError) as exc:
        raise ValueError("cube header is invalid") from exc
    atoms = abs(atom_count)
    if atoms > MAX_CUBE_ATOMS:
        raise ValueError(f"cube geometry is limited to {MAX_CUBE_ATOMS:,} atoms")
    for _index in range(atoms):
        line = header_line()
        if not line:
            raise ValueError("cube atom header is incomplete")
        lines.append(line)
    if atom_count < 0:
        origin = _record(lines[2], (4, 5))
        origin[0] = str(atoms)
        lines[2] = " ".join(origin) + "\n"
    return parse("".join(lines), header_only=True)


def combine(left_text: str, right_text: str, operation: str) -> str:
    left = parse(left_text)
    right = parse(right_text)
    if left.shape != right.shape:
        raise ValueError("cube grids have different dimensions")
    if left.datasets != right.datasets:
        raise ValueError("cube grids have different dataset counts")
    if left.axis_units != right.axis_units:
        raise ValueError("cube grids use different coordinate units")
    if left.dataset_ids != right.dataset_ids:
        raise ValueError("cube grids have different dataset identifiers")
    if len(left.geometry) != len(right.geometry) or any(
        len(a) != len(b) or any(abs(x - y) > 1.0e-8 for x, y in zip(a, b))
        for a, b in zip(left.geometry, right.geometry)
    ):
        raise ValueError("cube grids use different origins, axes, or atoms")
    if operation == "difference":
        values = [a - b for a, b in zip(left.values, right.values)]
        label = "OQP Studio cube difference (left - right)"
    elif operation == "sum":
        values = [a + b for a, b in zip(left.values, right.values)]
        label = "OQP Studio cube sum (left + right)"
    else:
        raise ValueError("cube operation must be difference or sum")
    if any(not isfinite(value) for value in values):
        raise ValueError("cube arithmetic produced a non-finite value")
    output = [label, left.header[1], *left.header[2:]]
    output.extend(" ".join(f"{value:13.5E}" for value in values[index:index + 6])
                  for index in range(0, len(values), 6))
    return "\n".join(output) + "\n"


def geometry_xyz(cube: Cube) -> str:
    atom_records = cube.geometry[4:]
    if not atom_records:
        raise ValueError("cube file contains no atomic geometry")
    factor = 1.0 if all(cube.axis_units) else BOHR_TO_ANGSTROM
    rows: list[str] = []
    for record in atom_records:
        atomic_value = record[0]
        if not atomic_value.is_integer():
            raise ValueError(f"cube atom has nonintegral atomic number {record[0]}")
        atomic_number = int(atomic_value)
        if atomic_number < 1 or atomic_number >= len(SYMBOLS):
            raise ValueError(f"cube atom has unsupported atomic number {atomic_number}")
        x, y, z = (value * factor for value in record[2:5])
        rows.append(f"{SYMBOLS[atomic_number]} {x:.10f} {y:.10f} {z:.10f}")
    return f"{len(rows)}\ncube geometry\n" + "\n".join(rows) + "\n"
