"""Read molecular structures out of the file formats users already have.

One entry point, :func:`parse`, turns a file's bytes into frames of atoms.
Several formats hold more than one geometry — an optimization log, an XYZ
trajectory, a multi-model PDB, a NAMD record — so every result is a list of
frames even when it has length one.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field

from .molden import ATOMIC_NUMBER, SYMBOLS, parse_molden

Atom = tuple[str, float, float, float]

BOHR_TO_ANGSTROM = 0.529177210903


class UnsupportedFormat(ValueError):
    """The file's format is not one this reader understands."""


@dataclass
class Frame:
    atoms: list[Atom]
    label: str = ""


@dataclass
class Structure:
    format: str
    frames: list[Frame] = field(default_factory=list)


def _symbol(token: str) -> str | None:
    """Element symbol from a symbol or an atomic number, or None."""
    token = token.strip()
    if not token:
        return None
    if token.replace(".", "", 1).isdigit():
        number = round(float(token))
        return SYMBOLS[number] if 0 < number < len(SYMBOLS) else None
    symbol = token[:2].capitalize() if len(token) > 1 else token.upper()
    if symbol in ATOMIC_NUMBER:
        return symbol
    head = token[:1].upper()
    return head if head in ATOMIC_NUMBER else None


def _atoms_from_lines(lines: list[str]) -> list[Atom]:
    """Read 'symbol x y z' rows, ignoring anything that does not match."""
    atoms: list[Atom] = []
    for line in lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        symbol = _symbol(parts[0])
        if symbol is None:
            continue
        try:
            x, y, z = (float(parts[i]) for i in (1, 2, 3))
        except ValueError:
            continue
        atoms.append((symbol, x, y, z))
    return atoms


def parse_xyz(text: str) -> list[Frame]:
    """XYZ, including multi-frame trajectories."""
    lines = text.splitlines()
    frames: list[Frame] = []
    i = 0
    while i < len(lines):
        header = lines[i].strip()
        if not header:
            i += 1
            continue
        try:
            count = int(header)
        except ValueError:
            break
        comment = lines[i + 1].strip() if i + 1 < len(lines) else ""
        atoms = _atoms_from_lines(lines[i + 2 : i + 2 + count])
        if atoms:
            frames.append(Frame(atoms, comment or f"frame {len(frames) + 1}"))
        i += 2 + count
    if not frames:  # headerless coordinate list
        atoms = _atoms_from_lines(lines)
        if atoms:
            frames.append(Frame(atoms))
    return frames


def parse_pdb(text: str) -> list[Frame]:
    """PDB, honouring MODEL records so NMR-style ensembles come through."""
    frames: list[Frame] = []
    atoms: list[Atom] = []
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            symbol = _symbol(line[76:78]) or _symbol(line[12:16].strip())
            if symbol is None:
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            atoms.append((symbol, x, y, z))
        elif line.startswith("ENDMDL") and atoms:
            frames.append(Frame(atoms, f"model {len(frames) + 1}"))
            atoms = []
    if atoms:
        frames.append(Frame(atoms))
    return frames


def parse_oqp(text: str) -> list[Frame]:
    """`.oqp` route input: the geometry lives in a geom=… value."""
    block = re.search(r'geom(?:etry)?\s*=\s*"""(.*?)"""', text, re.DOTALL | re.IGNORECASE)
    if block:
        return [Frame(_atoms_from_lines(block.group(1).splitlines()))]
    inline = re.search(r'geom(?:etry)?\s*=\s*"([^"]*)"', text, re.IGNORECASE)
    if inline:
        value = inline.group(1)
        if "\\n" in value:  # escaped one-line spelling
            return [Frame(_atoms_from_lines(value.replace("\\n", "\n").splitlines()))]
        raise UnsupportedFormat(
            f"this .oqp file points at a separate geometry file ({value}); open that file instead"
        )
    raise UnsupportedFormat("no geom=… value found in this .oqp file")


def parse_inp(text: str) -> list[Frame]:
    """Sectioned `.inp` input: coordinates follow `system=` in `[input]`."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if re.match(r"\s*system\s*=", line, re.IGNORECASE):
            tail = line.split("=", 1)[1]
            body = [tail, *lines[i + 1 :]]
            collected: list[str] = []
            for row in body:
                if re.match(r"\s*\[", row) or re.match(r"\s*\w+\s*=", row):
                    break
                collected.append(row)
            atoms = _atoms_from_lines(collected)
            if atoms:
                return [Frame(atoms)]
    raise UnsupportedFormat("no [input] system= geometry found in this .inp file")


def parse_openqp_log(text: str) -> list[Frame]:
    """Every Cartesian block in an OpenQP log, so an optimization plays back."""
    lines = text.splitlines()
    frames: list[Frame] = []
    step = 0
    i = 0
    while i < len(lines):
        marker = re.search(r"Geometry Optimization Step\s+(\d+)", lines[i])
        if marker:
            step = int(marker.group(1))
        if "Cartesian Coordinate in Angstrom" in lines[i]:
            i += 1
            while i < len(lines) and not re.match(r"\s*-{5,}", lines[i]):
                i += 1
            i += 1
            atoms: list[Atom] = []
            while i < len(lines):
                row = re.match(
                    r"\s*\d+\s+(\d+(?:\.\d+)?)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
                    lines[i],
                )
                if not row:
                    break
                symbol = _symbol(row.group(1))
                if symbol:
                    atoms.append(
                        (symbol, float(row.group(2)), float(row.group(3)), float(row.group(4)))
                    )
                i += 1
            if atoms:
                frames.append(Frame(atoms, f"step {step or len(frames) + 1}"))
            continue
        i += 1
    if not frames:
        raise UnsupportedFormat("no Cartesian coordinate block found in this log")
    return frames


def parse_openqp_json(text: str) -> list[Frame]:
    """OpenQP JSON results, including the geometry a Hessian run records."""
    data = json.loads(text)

    def atoms_from(source) -> list[Atom]:
        if not isinstance(source, dict):
            return []
        elements = (
            source.get("elements")
            or source.get("symbols")
            or source.get("atomic_numbers")
            or source.get("atoms")
        )
        coords = (
            source.get("coord")
            or source.get("coordinates")
            or source.get("geometry")
            or source.get("xyz")
        )
        if elements is None or coords is None:
            return []
        flat = coords
        if flat and isinstance(flat[0], (list, tuple)):
            flat = [value for row in flat for value in row]
        # OpenQP writes its own results in Bohr under "coord"; other producers
        # label the unit, and Angstrom is the usual default elsewhere.
        default_unit = "bohr" if "coord" in source else "angstrom"
        unit = str(source.get("units") or source.get("unit") or default_unit).lower()
        scale = BOHR_TO_ANGSTROM if unit.startswith(("bohr", "au")) else 1.0
        atoms: list[Atom] = []
        for index, element in enumerate(elements):
            symbol = _symbol(str(element))
            if symbol is None or 3 * index + 2 >= len(flat):
                continue
            x, y, z = (float(flat[3 * index + k]) * scale for k in range(3))
            atoms.append((symbol, x, y, z))
        return atoms

    for source in (data, data.get("molecule"), data.get("input"), data.get("results")):
        atoms = atoms_from(source)
        if atoms:
            return [Frame(atoms)]
    raise UnsupportedFormat("no geometry found in this JSON file")


def parse_namd_trajectory(path: str) -> list[Frame]:
    """OpenQP's packed NAMD record, read through OpenQP's own reader."""
    try:
        from oqp.library.namd import read_namd_trajectory
    except ImportError as exc:
        raise UnsupportedFormat(
            "reading .namd.trj needs OpenQP installed in this environment "
            "(the packed record is read with oqp.library.namd)"
        ) from exc

    metadata, trj = read_namd_trajectory(path)
    coordinates = trj["coordinates_angstrom"] if "coordinates_angstrom" in trj else trj["coordinates"]
    elements = metadata.get("elements") or metadata.get("symbols") or []
    times = trj["time_fs"] if "time_fs" in trj else range(len(coordinates))
    frames: list[Frame] = []
    for index, geometry in enumerate(coordinates):
        atoms: list[Atom] = []
        for atom_index, position in enumerate(geometry):
            symbol = _symbol(str(elements[atom_index])) if atom_index < len(elements) else "C"
            atoms.append((symbol or "C", float(position[0]), float(position[1]), float(position[2])))
        if atoms:
            frames.append(Frame(atoms, f"{float(times[index]):.1f} fs"))
    if not frames:
        raise UnsupportedFormat("the trajectory holds no geometries")
    return frames


def parse_with_rdkit(text: str, suffix: str) -> list[Frame]:
    """Chemistry formats RDKit already reads: MOL/SDF, MOL2, CDXML, SMILES."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise UnsupportedFormat(f"reading {suffix} needs RDKit on the backend") from exc

    if suffix in {".mol", ".sdf", ".sd"}:
        mol = Chem.MolFromMolBlock(text, removeHs=False)
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2Block(text, removeHs=False)
    elif suffix in {".cdxml", ".cdx"}:
        reader = getattr(Chem, "MolFromCDXML", None)
        if reader is None or suffix == ".cdx":
            raise UnsupportedFormat(
                "ChemDraw binary .cdx is not readable; save it as CDXML "
                "(File ▸ Save As ▸ ChemDraw XML) and open that"
            )
        mols = reader(text)
        mol = mols[0] if mols else None
    else:  # SMILES
        mol = Chem.MolFromSmiles(text.strip().splitlines()[0] if text.strip() else "")
    if mol is None:
        raise UnsupportedFormat(f"could not read this {suffix} file")

    if mol.GetNumConformers() == 0 or suffix in {".cdxml", ".smi", ".smiles", ".txt"}:
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise UnsupportedFormat("could not generate 3D coordinates for this structure")
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:  # noqa: BLE001, S110 — unoptimized coordinates still usable
            pass
    conformer = mol.GetConformer()
    atoms = [
        (atom.GetSymbol(), *(round(v, 6) for v in conformer.GetAtomPosition(i)))
        for i, atom in enumerate(mol.GetAtoms())
    ]
    return [Frame(atoms)]


# Suffixes are matched longest-first so ".hess.json" wins over ".json".
_TEXT_READERS: list[tuple[tuple[str, ...], str, callable]] = [
    ((".namd.restart.oqp", ".oqp"), "oqp", parse_oqp),
    ((".inp",), "inp", parse_inp),
    ((".xyz",), "xyz", parse_xyz),
    ((".pdb", ".ent"), "pdb", parse_pdb),
    ((".molden", ".freq.molden"), "molden", None),
    ((".hess.json", ".json"), "json", parse_openqp_json),
    ((".log", ".out", ".txt"), "log", parse_openqp_log),
]

_RDKIT_SUFFIXES = {".mol", ".sdf", ".sd", ".mol2", ".cdxml", ".cdx", ".smi", ".smiles"}


def parse(filename: str, data: bytes, path: str | None = None) -> Structure:
    """Read `data` according to `filename`'s extension."""
    name = filename.lower()

    if name.endswith((".namd.trj", ".trj")):
        if path is None:
            raise UnsupportedFormat("the packed trajectory must be read from a file on disk")
        return Structure("namd.trj", parse_namd_trajectory(path))

    if name.endswith((".molden", ".freq.molden")):
        molden = parse_molden(data.decode(errors="replace"))
        if not molden.atoms:
            raise UnsupportedFormat("no atoms found in this Molden file")
        return Structure("molden", [Frame(list(molden.atoms))])

    text = data.decode(errors="replace")

    for suffix in sorted(_RDKIT_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix):
            return Structure(suffix.lstrip("."), parse_with_rdkit(text, suffix))

    for suffixes, label, reader in _TEXT_READERS:
        if reader and any(name.endswith(s) for s in suffixes):
            frames = reader(text)
            if frames:
                return Structure(label, frames)

    # Unknown extension: try the formats that identify themselves from content.
    for label, reader in (("xyz", parse_xyz), ("pdb", parse_pdb), ("log", parse_openqp_log)):
        try:
            frames = reader(text)
        except (UnsupportedFormat, ValueError):
            continue
        if frames:
            return Structure(label, frames)
    raise UnsupportedFormat(f"could not tell what format {filename} is")
