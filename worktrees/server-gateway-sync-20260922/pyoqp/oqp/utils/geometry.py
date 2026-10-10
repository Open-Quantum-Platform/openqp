"""Geometry helpers for OpenQP Python inputs."""

import os
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen

from oqp.utils.constants import ANGSTROM_TO_BOHR


_BOHR_UNITS = {"bohr", "b", "au", "a.u.", "atomic_unit", "atomic_units"}
_GEOMETRY_ALIASES = {
    "h2": "\nH 0.000000 0.000000 -0.370000\nH 0.000000 0.000000 0.370000",
    "hydrogen": "\nH 0.000000 0.000000 -0.370000\nH 0.000000 0.000000 0.370000",
    "water": "\nO 0.000000 0.000000 0.000000\nH 0.758602 0.000000 0.504284\nH -0.758602 0.000000 0.504284",
    "h2o": "\nO 0.000000 0.000000 0.000000\nH 0.758602 0.000000 0.504284\nH -0.758602 0.000000 0.504284",
    "methane": "\nC 0.000000 0.000000 0.000000\nH 0.629118 0.629118 0.629118\nH -0.629118 -0.629118 0.629118\nH -0.629118 0.629118 -0.629118\nH 0.629118 -0.629118 -0.629118",
    "ch4": "\nC 0.000000 0.000000 0.000000\nH 0.629118 0.629118 0.629118\nH -0.629118 -0.629118 0.629118\nH -0.629118 0.629118 -0.629118\nH 0.629118 -0.629118 -0.629118",
    "ammonia": "\nN 0.000000 0.000000 0.116000\nH 0.000000 0.939000 -0.271000\nH 0.813000 -0.469500 -0.271000\nH -0.813000 -0.469500 -0.271000",
    "nh3": "\nN 0.000000 0.000000 0.116000\nH 0.000000 0.939000 -0.271000\nH 0.813000 -0.469500 -0.271000\nH -0.813000 -0.469500 -0.271000",
    "co2": "\nO 0.000000 0.000000 -1.160000\nC 0.000000 0.000000 0.000000\nO 0.000000 0.000000 1.160000",
    "carbon dioxide": "\nO 0.000000 0.000000 -1.160000\nC 0.000000 0.000000 0.000000\nO 0.000000 0.000000 1.160000",
}


class GeometryLookupError(ValueError):
    """Raised when a named geometry cannot be resolved."""


def builtin_geometry(query):
    """Return a built-in approximate Angstrom geometry for common small molecules."""
    key = _geometry_key(query)
    try:
        return _GEOMETRY_ALIASES[key]
    except KeyError as exc:
        raise GeometryLookupError(f"No built-in geometry is available for '{query}'.") from exc


def pubchem_geometry(query, timeout=10):
    """Fetch a 3D conformer from PubChem PUG-REST and return OpenQP geometry text."""
    encoded = quote(str(query).strip(), safe="")
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded}/SDF?record_type=3d"
    )
    try:
        with urlopen(url, timeout=timeout) as response:
            sdf = response.read().decode("utf-8", errors="replace")
    except (HTTPError, URLError, TimeoutError, OSError) as exc:
        raise GeometryLookupError(f"PubChem geometry lookup failed for '{query}'.") from exc
    return geometry_from_sdf(sdf, query)


def get_geometry(query, source="auto", timeout=10):
    """
    Resolve a named molecular geometry to OpenQP inline geometry text.

    source="auto" tries OpenQP's small built-in table first, then PubChem.
    source="builtin" only uses the local table.
    source="pubchem" only uses PubChem PUG-REST.
    """
    mode = str(source).strip().lower()
    if mode in ("builtin", "local", "generator"):
        return builtin_geometry(query)
    if mode == "pubchem":
        return pubchem_geometry(query, timeout=timeout)
    if mode != "auto":
        raise ValueError("geometry source must be auto, builtin, or pubchem.")

    try:
        return builtin_geometry(query)
    except GeometryLookupError:
        return pubchem_geometry(query, timeout=timeout)


def geometry_from_sdf(sdf, query):
    """Convert a PubChem-style V2000 SDF string to OpenQP inline geometry text."""
    lines = sdf.splitlines()
    counts_index = None
    atom_count = None
    for idx, line in enumerate(lines):
        if "V2000" not in line and "V3000" not in line:
            continue
        fields = line.split()
        try:
            atom_count = int(fields[0])
        except (IndexError, ValueError) as exc:
            raise GeometryLookupError(f"Could not read PubChem atom count for '{query}'.") from exc
        counts_index = idx
        break

    if counts_index is None or atom_count is None:
        raise GeometryLookupError(f"PubChem did not return a readable SDF geometry for '{query}'.")

    rows = []
    for line in lines[counts_index + 1:counts_index + 1 + atom_count]:
        fields = line.split()
        if len(fields) < 4:
            raise GeometryLookupError(f"PubChem returned an incomplete atom row for '{query}'.")
        try:
            x = float(fields[0])
            y = float(fields[1])
            z = float(fields[2])
        except ValueError as exc:
            raise GeometryLookupError(f"PubChem returned nonnumeric coordinates for '{query}'.") from exc
        rows.append(f"{fields[3]} {x:.10g} {y:.10g} {z:.10g}")

    if not rows:
        raise GeometryLookupError(f"PubChem returned an empty geometry for '{query}'.")
    return "\n" + "\n".join(rows)


def normalize_system(system, unit="Angstrom"):
    """
    Normalize inline molecular geometries for OpenQP input.system.

    Coordinates are written in Angstrom because the OpenQP input reader converts
    the text geometry into its internal Bohr representation. Bohr inputs can be
    passed with unit="Bohr".
    """
    if isinstance(system, (list, tuple)):
        lines = [_normalize_atom_row(item, unit) for item in system]
    elif isinstance(system, str):
        s = system.strip()

        if os.path.exists(s):
            return s

        if "\n" in s:
            lines = [_normalize_system_line(ln, unit) for ln in s.split("\n") if ln.strip()]
        elif ";" in s or "," in s:
            sep = ";" if ";" in s else ","
            parts = [p.strip() for p in s.split(sep) if p.strip()]
            if not all(len(p.split()) >= 4 for p in parts):
                raise ValueError("Inline atom must be 'El x y z' per entry.")
            lines = [_normalize_system_line(part, unit) for part in parts]
        else:
            toks = s.split()
            if len(toks) >= 4:
                lines = [_normalize_system_line(s, unit)]
            else:
                return s
    else:
        raise TypeError("input.system must be str or list/tuple of atoms.")

    return "\n" + "\n".join(lines)


def _geometry_key(query):
    return " ".join(str(query).strip().lower().replace("_", " ").split())


def _is_bohr(unit):
    return str(unit).strip().lower() in _BOHR_UNITS


def _format_coord(value, unit):
    if not _is_bohr(unit):
        return str(value)
    coord = float(str(value).replace("D", "E").replace("d", "e"))
    return f"{coord * ANGSTROM_TO_BOHR:.12g}"


def _is_coord_sequence(value):
    if isinstance(value, (str, bytes)):
        return False
    try:
        iter(value)
        len(value)
    except TypeError:
        return False
    return True


def _normalize_atom_row(item, unit):
    if isinstance(item, str):
        return _normalize_system_line(item, unit)

    if not _is_coord_sequence(item) or len(item) < 2:
        raise ValueError("Each atom must be 'El x y z', (symbol, x, y, z), or (symbol, (x, y, z)).")

    element = item[0]
    if len(item) >= 4:
        coords = item[1:4]
        extra = item[4:]
    elif _is_coord_sequence(item[1]) and len(item[1]) >= 3:
        coords = item[1][:3]
        extra = item[2:]
    else:
        raise ValueError("Each atom must be 'El x y z', (symbol, x, y, z), or (symbol, (x, y, z)).")

    fields = [str(element), *(_format_coord(coord, unit) for coord in coords), *(str(x) for x in extra)]
    return " ".join(fields)


def _normalize_system_line(line, unit):
    tokens = line.split()
    if len(tokens) < 4:
        return line.strip()
    fields = [tokens[0], *(_format_coord(coord, unit) for coord in tokens[1:4]), *tokens[4:]]
    return " ".join(fields)
