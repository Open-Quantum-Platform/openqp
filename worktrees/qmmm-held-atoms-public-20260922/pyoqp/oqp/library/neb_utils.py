"""Dependency-light helpers for geomeTRIC NEB endpoint paths."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


BOHR_TO_ANGSTROM = 0.52917721092


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


def kabsch_align(
    reference: Iterable[Iterable[float]],
    mobile: Iterable[Iterable[float]],
) -> list[list[float]]:
    """Rigidly align ``mobile`` coordinates onto ``reference``.

    Coordinates may use any common unit.  The returned coordinates retain the
    reference centroid, and the fitted transform is constrained to a proper
    rotation so a reflection cannot change molecular handedness.
    """

    reference_array = np.asarray(reference, dtype=float)
    mobile_array = np.asarray(mobile, dtype=float)
    if reference_array.shape != mobile_array.shape:
        raise ValueError("Kabsch alignment requires coordinate arrays with the same shape")
    if reference_array.ndim != 2 or reference_array.shape[1] != 3:
        raise ValueError("Kabsch alignment requires coordinate arrays shaped (natom, 3)")
    if reference_array.shape[0] == 0:
        raise ValueError("Kabsch alignment requires at least one atom")

    reference_centroid = reference_array.mean(axis=0)
    mobile_centroid = mobile_array.mean(axis=0)
    reference_centered = reference_array - reference_centroid
    mobile_centered = mobile_array - mobile_centroid

    covariance = mobile_centered.T @ reference_centered
    left, _, right_t = np.linalg.svd(covariance)
    if np.linalg.det(left @ right_t) < 0.0:
        left[:, -1] *= -1.0
    rotation = left @ right_t
    aligned = mobile_centered @ rotation + reference_centroid
    return aligned.tolist()


def interpolate_xyz_endpoints(
    reactant_xyz: str | Path,
    product_xyz: str | Path,
    nimage: int,
    align: bool = False,
) -> list[NEBImage]:
    """Return linearly interpolated NEB images including both endpoints.

    Input and output coordinates are Angstrom, matching OpenQP XYZ/system input
    conventions. Atom symbols/count/order must match exactly so each image can
    be evaluated on the same state-specific surface.  With ``align=True``, the
    product endpoint is first rigidly superimposed on the reactant, removing
    overall translation and rotation before interpolation.
    """

    if nimage < 3:
        raise ValueError("NEB interpolation requires nimage >= 3")
    reactant = _read_xyz(reactant_xyz)
    product = _read_xyz(product_xyz)
    if reactant.symbols != product.symbols:
        raise ValueError("NEB endpoints must have identical atom symbols and order")

    product_coordinates = product.coordinates_angstrom
    if align:
        product_coordinates = kabsch_align(
            reactant.coordinates_angstrom,
            product_coordinates,
        )

    denominator = nimage - 1
    images: list[NEBImage] = []
    for image_index in range(nimage):
        fraction = image_index / denominator
        images.append(
            NEBImage(
                symbols=list(reactant.symbols),
                coordinates_angstrom=_interpolate_coordinates(
                    reactant.coordinates_angstrom,
                    product_coordinates,
                    fraction,
                ),
            )
        )
    return images


def write_neb_xyz(
    path: str | Path,
    symbols: Iterable[str],
    images_bohr: Iterable[Iterable[float]],
    energies: Iterable[float] | None = None,
) -> Path:
    """Write a final NEB path as a standard multi-frame XYZ file.

    ``images_bohr`` accepts flat ``3*natom`` arrays or ``(natom, 3)`` arrays.
    Optional energies are written in Hartree in each frame's comment line.
    """

    xyz_path = Path(path)
    xyz_path.parent.mkdir(parents=True, exist_ok=True)
    symbol_list = [str(symbol) for symbol in symbols]
    if not symbol_list:
        raise ValueError("NEB XYZ output requires at least one atom symbol")

    image_list = [np.asarray(image, dtype=float) for image in images_bohr]
    if not image_list:
        raise ValueError("NEB XYZ output requires at least one image")

    energy_list = None if energies is None else [float(energy) for energy in energies]
    if energy_list is not None and len(energy_list) != len(image_list):
        raise ValueError("NEB XYZ energies must match the number of images")

    lines: list[str] = []
    nimage = len(image_list)
    natom = len(symbol_list)
    for image_index, image in enumerate(image_list):
        try:
            coordinates = image.reshape((natom, 3)) * BOHR_TO_ANGSTROM
        except ValueError as exc:
            raise ValueError(
                "Each NEB image must contain exactly 3*natom coordinates"
            ) from exc

        lines.append(str(natom))
        comment = f"NEB image {image_index + 1}/{nimage}"
        if energy_list is not None:
            comment += f" energy_hartree={energy_list[image_index]:.16f}"
        lines.append(comment)
        for symbol, (x_coord, y_coord, z_coord) in zip(symbol_list, coordinates):
            lines.append(
                f"{symbol:<3s} {x_coord: .16f} {y_coord: .16f} {z_coord: .16f}"
            )

    xyz_path.write_text("\n".join(lines) + "\n")
    return xyz_path
