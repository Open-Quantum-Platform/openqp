"""Dependency-light helpers for geomeTRIC NEB endpoint paths."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class NEBImage:
    """One NEB path image in user-facing Angstrom coordinates."""

    symbols: list[str]
    coordinates_angstrom: list[list[float]]


def _read_xyz(path: str | Path) -> NEBImage:
    xyz_path = Path(path)
    lines = xyz_path.read_text().splitlines()
    if not lines:
        raise ValueError(f"XYZ endpoint is empty: {xyz_path}")
    try:
        natom = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"XYZ endpoint first line must be an atom count: {xyz_path}") from exc
    atom_lines = lines[2:]
    if len(atom_lines) < natom:
        raise ValueError(f"XYZ endpoint has fewer atom lines than declared: {xyz_path}")

    symbols: list[str] = []
    coordinates: list[list[float]] = []
    for index, line in enumerate(atom_lines[:natom], start=1):
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"XYZ atom line {index} must contain symbol and x/y/z: {xyz_path}")
        symbols.append(fields[0])
        try:
            coordinates.append([float(fields[1]), float(fields[2]), float(fields[3])])
        except ValueError as exc:
            raise ValueError(f"XYZ atom line {index} contains non-numeric coordinates: {xyz_path}") from exc
    return NEBImage(symbols=symbols, coordinates_angstrom=coordinates)


def _interpolate_coordinates(
    reactant: Iterable[Iterable[float]],
    product: Iterable[Iterable[float]],
    fraction: float,
) -> list[list[float]]:
    return [
        [r_value + fraction * (p_value - r_value) for r_value, p_value in zip(r_atom, p_atom)]
        for r_atom, p_atom in zip(reactant, product)
    ]


def interpolate_xyz_endpoints(
    reactant_xyz: str | Path,
    product_xyz: str | Path,
    nimage: int,
) -> list[NEBImage]:
    """Return linearly interpolated NEB images including both endpoints.

    Input and output coordinates are Angstrom, matching OpenQP XYZ/system input
    conventions. Atom symbols/count/order must match exactly so each image can
    be evaluated on the same state-specific surface.
    """

    if nimage < 3:
        raise ValueError("NEB interpolation requires nimage >= 3")
    reactant = _read_xyz(reactant_xyz)
    product = _read_xyz(product_xyz)
    if reactant.symbols != product.symbols:
        raise ValueError("NEB endpoints must have identical atom symbols and order")

    denominator = nimage - 1
    images: list[NEBImage] = []
    for image_index in range(nimage):
        fraction = image_index / denominator
        images.append(
            NEBImage(
                symbols=list(reactant.symbols),
                coordinates_angstrom=_interpolate_coordinates(
                    reactant.coordinates_angstrom,
                    product.coordinates_angstrom,
                    fraction,
                ),
            )
        )
    return images
