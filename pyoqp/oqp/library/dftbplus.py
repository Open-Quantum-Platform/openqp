"""External DFTB+ backend helpers for OpenQP.

This module intentionally keeps the DFTB+ dependency optional.  Fixture tests
exercise parsing and input generation without requiring the ``dftb+`` executable
or Slater-Koster parameter files to be installed on CI/developer machines.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile
from typing import Iterable, Sequence

BOHR_TO_ANGSTROM = 0.529177210903
_SYMBOLS = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
    19: "K", 20: "Ca", 35: "Br", 53: "I",
}


class DFTBPlusError(RuntimeError):
    """Raised when the optional external DFTB+ backend cannot complete."""


@dataclass
class DFTBPlusResult:
    energy: float | None = None
    gradient: list[list[float]] | None = None
    stdout: str = ""
    stderr: str = ""
    workdir: str | None = None


def _read_text(path: str | os.PathLike[str]) -> str:
    return Path(path).read_text(encoding="utf-8")


def parse_detailed_out(path: str | os.PathLike[str]) -> DFTBPlusResult:
    """Parse the total energy from a DFTB+ ``detailed.out`` file."""

    text = _read_text(path)
    energy = None
    for pattern in (
        r"Total\s+Energy\s*:\s*([-+0-9.Ee]+)\s*H",
        r"Total\s+Mermin\s+free\s+energy\s*:\s*([-+0-9.Ee]+)\s*H",
    ):
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            energy = float(match.group(1))
            break
    if energy is None:
        raise DFTBPlusError(f"Could not parse total energy from {path}")
    return DFTBPlusResult(energy=energy)


def _parse_tag_header(header: str) -> tuple[str, list[int]]:
    # Example: forces:real:2:3,3
    fields = header.strip().split(":")
    name = fields[0].strip()
    shape = []
    if len(fields) >= 4 and fields[3]:
        shape = [int(x) for x in fields[3].split(",") if x]
    return name, shape


def parse_results_tag(path: str | os.PathLike[str]) -> DFTBPlusResult:
    """Parse energy and gradients from DFTB+ ``results.tag``.

    DFTB+ stores forces in Hartree/Bohr. OpenQP gradients are dE/dR, so this
    parser returns ``gradient = -forces``.
    """

    lines = _read_text(path).splitlines()
    energy = None
    gradient = None
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line or ":" not in line:
            index += 1
            continue
        name, shape = _parse_tag_header(line)
        index += 1
        count = 1
        if shape:
            count = 1
            for dim in shape:
                count *= dim
        values: list[float] = []
        while index < len(lines) and len(values) < count:
            next_line = lines[index].strip()
            if next_line:
                values.extend(float(token) for token in next_line.split())
            index += 1
        if name in {"total_energy", "mermin_energy", "extrapolated0_energy", "forcerelated_energy"} and values:
            if energy is None or name == "total_energy":
                energy = values[0]
        elif name == "forces" and values:
            if not shape or len(shape) != 2 or 3 not in shape:
                raise DFTBPlusError(f"Unexpected forces shape in {path}: {shape}")
            if shape[1] == 3:
                natom = shape[0]
                gradient = [[-values[3 * i + j] for j in range(3)] for i in range(natom)]
            else:
                natom = shape[1]
                gradient = [[-values[3 * i + j] for j in range(3)] for i in range(natom)]
    if energy is None:
        raise DFTBPlusError(f"Could not parse total_energy from {path}")
    return DFTBPlusResult(energy=energy, gradient=gradient)


def _symbols_for_atoms(atoms: Sequence[int]) -> list[str]:
    symbols = []
    for atom in atoms:
        atomic_number = int(atom)
        try:
            symbols.append(_SYMBOLS[atomic_number])
        except KeyError as exc:
            raise DFTBPlusError(f"No DFTB+ symbol mapping for atomic number {atomic_number}") from exc
    return symbols


def _coords_rows(coords_bohr: Sequence[float]) -> list[list[float]]:
    if len(coords_bohr) % 3:
        raise DFTBPlusError("Coordinate array length must be a multiple of 3")
    return [
        [float(coords_bohr[3 * i + j]) * BOHR_TO_ANGSTROM for j in range(3)]
        for i in range(len(coords_bohr) // 3)
    ]


_DEFAULT_MAX_ANGULAR_MOMENTUM = {
    "H": "s", "He": "s",
    "Li": "p", "Be": "p", "B": "p", "C": "p", "N": "p", "O": "p", "F": "p", "Ne": "p",
    "Na": "p", "Mg": "p", "Al": "p", "Si": "p", "P": "p", "S": "p", "Cl": "p", "Ar": "p",
    "K": "p", "Ca": "p", "Br": "p", "I": "p",
    "Cr": "d",
}


def _max_angular_momentum_block(symbols: Sequence[str], config: dict) -> str:
    dftb = _dftb_config(config)
    overrides = dftb.get("max_angular_momentum", {}) or {}
    lines = ["  MaxAngularMomentum = {"]
    for symbol in dict.fromkeys(symbols):
        value = overrides.get(symbol, _DEFAULT_MAX_ANGULAR_MOMENTUM.get(symbol))
        if value is None:
            raise DFTBPlusError(f"No DFTB+ MaxAngularMomentum default for atom type {symbol}")
        lines.append(f'    {symbol} = "{value}"')
    lines.append("  }")
    return "\n".join(lines)


def _dftb_config(config: dict) -> dict:
    return config.get("dftb", {}) if config else {}


def write_dftbplus_input(
    workdir: str | os.PathLike[str],
    atoms: Sequence[int],
    coords_bohr: Sequence[float],
    config: dict,
    *,
    gradient: bool = False,
) -> Path:
    """Write ``dftb_in.hsd`` and ``geo.gen`` for an external DFTB+ run."""

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    dftb = _dftb_config(config)
    sk_path = str(dftb.get("sk_path", "")).strip()
    if sk_path and not sk_path.endswith(os.sep):
        sk_path += os.sep
    scc = bool(dftb.get("scc", True))
    max_scc = int(dftb.get("max_scc_iterations", 100))
    charge = int(config.get("input", {}).get("charge", 0))
    symbols = _symbols_for_atoms(atoms)
    unique_symbols = list(dict.fromkeys(symbols))
    max_l_block = _max_angular_momentum_block(unique_symbols, config)
    rows = _coords_rows(coords_bohr)
    if len(rows) != len(symbols):
        raise DFTBPlusError("Number of atoms does not match coordinate rows")

    geo_path = workdir / "geo.gen"
    with geo_path.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(symbols)} C\n")
        handle.write(" ".join(unique_symbols) + "\n")
        for idx, (symbol, xyz) in enumerate(zip(symbols, rows), start=1):
            type_index = unique_symbols.index(symbol) + 1
            handle.write(f"{idx:5d} {type_index:3d} {xyz[0]:18.10f} {xyz[1]:18.10f} {xyz[2]:18.10f}\n")

    hsd = f'''Geometry = GenFormat {{
  <<< "geo.gen"
}}
# AtomTypes = {" ".join(unique_symbols)}

Hamiltonian = DFTB {{
  SCC = {"Yes" if scc else "No"}
  MaxSCCIterations = {max_scc}
  Charge = {charge}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{sk_path}"
    Separator = "-"
    Suffix = ".skf"
  }}
{max_l_block}
}}

Options {{
  WriteResultsTag = Yes
}}
'''
    if gradient:
        hsd += "\nDriver = ConjugateGradient { MaxSteps = 0 }\n"
    input_path = workdir / "dftb_in.hsd"
    input_path.write_text(hsd, encoding="utf-8")
    return input_path


class DFTBPlusRunner:
    """Small shell-out wrapper around the optional ``dftb+`` executable."""

    def __init__(self, config: dict):
        self.config = config or {}
        dftb = _dftb_config(self.config)
        self.executable = str(dftb.get("executable", "dftb+") or "dftb+")
        self.sk_path = str(dftb.get("sk_path", "") or "")
        self.keep_workdir = bool(dftb.get("keep_workdir", False))
        self.timeout = int(dftb.get("timeout", 300) or 300)

    def _check_prerequisites(self) -> None:
        exe = self.executable
        if os.path.sep in exe:
            if not (Path(exe).exists() and os.access(exe, os.X_OK)):
                raise DFTBPlusError(f"DFTB+ executable not found or not executable: {exe}")
        elif shutil.which(exe) is None:
            raise DFTBPlusError(f"DFTB+ executable not found on PATH: {exe}")
        if not self.sk_path:
            raise DFTBPlusError("DFTB+ parameter directory is not configured; set [dftb] sk_path")
        if not Path(self.sk_path).is_dir():
            raise DFTBPlusError(f"DFTB+ parameter directory not found: {self.sk_path}")

    def run(self, atoms: Sequence[int], coords_bohr: Sequence[float], *, gradient: bool = False) -> DFTBPlusResult:
        self._check_prerequisites()
        if self.keep_workdir:
            workdir = Path(tempfile.mkdtemp(prefix="openqp-dftbplus-"))
            return self._run_in_workdir(workdir, atoms, coords_bohr, gradient=gradient)

        with tempfile.TemporaryDirectory(prefix="openqp-dftbplus-") as tmp:
            return self._run_in_workdir(Path(tmp), atoms, coords_bohr, gradient=gradient)

    def _run_in_workdir(self, workdir: Path, atoms: Sequence[int], coords_bohr: Sequence[float], *, gradient: bool) -> DFTBPlusResult:
        write_dftbplus_input(workdir, atoms, coords_bohr, self.config, gradient=gradient)
        proc = subprocess.run(
            [self.executable],
            cwd=workdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=self.timeout,
            check=False,
        )
        if proc.returncode != 0:
            raise DFTBPlusError(
                f"DFTB+ exited with status {proc.returncode}: {proc.stderr.strip() or proc.stdout.strip()}"
            )
        result_path = workdir / "results.tag"
        detail_path = workdir / "detailed.out"
        if result_path.exists():
            result = parse_results_tag(result_path)
        elif detail_path.exists():
            result = parse_detailed_out(detail_path)
        else:
            raise DFTBPlusError("DFTB+ completed but neither results.tag nor detailed.out was written")
        result.stdout = proc.stdout
        result.stderr = proc.stderr
        result.workdir = str(workdir)
        return result


def run_openqp_molecule(mol, *, gradient: bool = False) -> DFTBPlusResult:
    """Execute DFTB+ for an OpenQP Molecule-like object and populate results."""

    result = DFTBPlusRunner(mol.config).run(mol.get_atoms(), mol.get_system(), gradient=gradient)
    if result.energy is not None:
        mol.data._data.mol_energy.energy = result.energy
        if hasattr(mol, "energies"):
            mol.energies = [result.energy]
    if gradient and result.gradient is not None:
        flat = [component for row in result.gradient for component in row]
        try:
            import numpy as np
            from oqp import ffi

            grad_array = np.asarray(flat, dtype=np.float64)
            ffi.memmove(mol.data._data.grad, grad_array, grad_array.nbytes)
        except Exception:
            if hasattr(mol, "grads"):
                mol.grads = [flat]
    return result
