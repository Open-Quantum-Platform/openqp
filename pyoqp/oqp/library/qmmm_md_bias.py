"""ODP umbrella and conservative restraints for ground-state QM/MM MD.

The OpenMM-integrated QM/MM MD driver (``md(...) qmmm(...)``) adds these
terms to the QM/MM energy and force before handing them to the integrator.
They have the same meaning as in QM/MM NAMD:

- ``odp(...)``: one-dimensional projection umbrella (native kernel); CV atom
  labels are 1-based over all PDB atoms.
- ``droplet(...)``: spherical wall on solvent molecules (native kernel).
- ``solute_com(...)``: harmonic restraint of a centre of mass to a point
  (native kernel).

Everything here is in atomic units; the caller converts to and from OpenMM
units.  ``series_record`` and ``series_metadata`` define the per-step record
that ``oqp.library.odp.odp_wham`` reads back.
"""

from __future__ import annotations

import hashlib
import json
import os
import re

import numpy as np

import oqp
from oqp.library.odp import odp_from_config

BOHR_TO_ANGSTROM = 0.52917721090
ANGSTROM_TO_BOHR = 1.0/BOHR_TO_ANGSTROM
KCALMOLANG2_TO_HARTREEBOHR2 = (1.0/627.5094740631)/ANGSTROM_TO_BOHR**2
SERIES_FORMAT = "openqp-qmmm-md-odp-series-1"


def _as_bool(value):
    return value is True or str(value).strip().lower() in {
        "true", "t", "1", "yes", "on", ".true."}


def _vector3(value, label):
    if isinstance(value, str):
        value = [item for item in value.replace(",", " ").split() if item]
    array = np.asarray(value, dtype=np.float64).reshape(-1)
    if array.size != 3 or not np.all(np.isfinite(array)):
        raise ValueError(f"{label} requires three finite values")
    return array


def _index_list(spec, natom, label):
    """0-based atom indices of a list such as "0-3,7 9:11"; every index must
    lie in [0, natom) and every range must ascend (no index from the end)."""
    out = []
    for token in str(spec or "").replace(",", " ").split():
        match = re.fullmatch(r"(\d+)(?:[-:](\d+))?", token)
        if match is None:
            raise ValueError(
                f"{label} atoms: {token!r} is not a 0-based index or range")
        first = int(match.group(1))
        last = int(match.group(2)) if match.group(2) is not None else first
        if last < first:
            raise ValueError(f"{label} atoms: range {token!r} descends")
        if last >= natom:
            raise ValueError(
                f"{label} atoms: index {last} is outside 0..{natom - 1}")
        out.extend(range(first, last + 1))
    return out


def requested_sections(sections):
    """Names of the bias sections that are enabled in a config mapping."""
    return [name for name in ("odp", "droplet", "solute_com")
            if _as_bool(sections.get(name, {}).get("enabled", False))]


class QMMMMDBias:
    """ODP umbrella plus optional droplet and solute-COM restraints."""

    def __init__(self, sections, topology, qm_atoms, masses, periodic):
        self.odp = odp_from_config(sections)
        natom = int(topology.getNumAtoms())
        self.natom = natom
        masses = np.asarray(masses, dtype=np.float64)
        # Virtual sites (TIP4P/TIP5P) carry no mass: they belong to no
        # restraint group, and the native kernels, which require positive
        # masses for every atom, get a placeholder that is never used.
        self._massless = masses <= 0.0
        self.masses = np.ascontiguousarray(np.where(self._massless, 1.0, masses))
        if self.odp is not None:
            self.odp.validate_atom_count(natom)

        droplet = sections.get("droplet", {})
        self.droplet = _as_bool(droplet.get("enabled", False))
        self._droplet_groups = None
        if self.droplet:
            if periodic:
                raise ValueError(
                    "[droplet] finite spherical containment requires "
                    "qmmm(cutoff=NoCutoff)")
            self.droplet_settings = {
                "center_angstrom": _vector3(
                    droplet.get("center", (0.0, 0.0, 0.0)), "[droplet] center").tolist(),
                "radius_angstrom": float(droplet.get("radius", 20.0)),
                "buffer_angstrom": float(droplet.get("buffer", 1.0)),
                "force_constant_kcal_mol_angstrom2": float(
                    droplet.get("force_constant", 10.0)),
                "max_penetration_angstrom": float(
                    droplet.get("max_penetration", 10.0)),
                "target": str(droplet.get("target", "water_com")).strip().lower(),
                "atoms_zero_based": str(droplet.get("atoms", "") or "").strip(),
            }
            names = droplet.get("water_resnames", "hoh,wat,sol,tip3,tip3p")
            if isinstance(names, str):
                names = names.replace(",", " ").split()
            names = tuple(str(n).strip().lower() for n in names if str(n).strip())
            self.droplet_settings["water_resnames"] = list(names)
            s = self.droplet_settings
            if (s["radius_angstrom"] <= 0.0 or s["buffer_angstrom"] < 0.0
                    or s["force_constant_kcal_mol_angstrom2"] <= 0.0
                    or s["max_penetration_angstrom"] < 0.0):
                raise ValueError(
                    "[droplet] radius/force_constant must be positive and "
                    "buffer/max_penetration non-negative")
            target = s["target"].replace("-", "_")
            if target in ("com", "molecule_com", "water_molecule_com"):
                target = "water_com"
            if target in ("o", "water_oxygen"):
                target = "oxygen"
            if target not in ("water_com", "oxygen", "atoms"):
                raise ValueError("[droplet] target must be water_com, oxygen, or atoms")
            s["target"] = target
            groups = np.zeros(natom, dtype=np.int64)
            if target == "atoms":
                for k, atom in enumerate(sorted(set(_index_list(s["atoms_zero_based"], natom,
                                                                          "[droplet]")))):
                    groups[atom] = k + 1
            else:
                group = 0
                for residue in topology.residues():
                    if str(residue.name).strip().lower() not in names:
                        continue
                    atoms = list(residue.atoms())
                    atoms = [a for a in atoms if not self._massless[int(a.index)]]
                    if target == "oxygen":
                        atoms = [a for a in atoms if getattr(
                            getattr(a, "element", None), "symbol", "") == "O"][:1]
                    if not atoms:
                        continue
                    group += 1
                    for atom in atoms:
                        groups[int(atom.index)] = group
            if int(groups.max()) == 0:
                raise ValueError("[droplet] selected no solvent molecules")
            self._droplet_groups = groups
            self._droplet_center = (np.asarray(s["center_angstrom"])
                                    * ANGSTROM_TO_BOHR)

        solute = sections.get("solute_com", {})
        self.solute_com = _as_bool(solute.get("enabled", False))
        self._solute_selected = None
        if self.solute_com:
            spec = str(solute.get("atoms", "") or "").strip()
            atoms = (_index_list(spec, natom, "[solute_com]") if spec
                     else [int(i) for i in qm_atoms])
            selected = np.zeros(natom, dtype=np.int64)
            selected[sorted(set(atoms))] = 1
            selected[self._massless] = 0
            self.solute_settings = {
                "center_angstrom": _vector3(
                    solute.get("center", (0.0, 0.0, 0.0)),
                    "[solute_com] center").tolist(),
                "force_constant_kcal_mol_angstrom2": float(
                    solute.get("force_constant", 5.0)),
                "atoms_zero_based": str(solute.get("atoms", "") or "").strip(),
            }
            if not self.solute_settings["force_constant_kcal_mol_angstrom2"] > 0.0:
                raise ValueError("[solute_com] force_constant must be positive")
            self._solute_selected = selected
            self._solute_center = (np.asarray(self.solute_settings["center_angstrom"])
                                   * ANGSTROM_TO_BOHR)
        self.last = None

    @property
    def enabled(self):
        return self.odp is not None or self.droplet or self.solute_com

    def evaluate(self, xyz_bohr, box_bohr=None):
        """Energy (Hartree), force (Hartree/bohr) and the per-step record."""
        coords = np.ascontiguousarray(xyz_bohr, dtype=np.float64).reshape((-1, 3))
        force = np.zeros_like(coords)
        record = {"odp_energy_hartree": 0.0, "droplet_energy_hartree": 0.0,
                  "solute_com_energy_hartree": 0.0}
        if self.odp is not None:
            odp = self.odp.evaluate(coords, box=box_bohr)
            force += odp["force"]
            record.update({
                "odp_energy_hartree": odp["energy"],
                "odp_window": self.odp.window,
                "odp_xi": odp["xi"],
                "odp_cv_raw": odp["cv_raw"],
                "odp_cv_scaled": odp["cv_scaled"],
                "odp_cv_perpendicular": odp["cv_perpendicular"],
                "odp_perpendicular_norm": odp["perpendicular_norm"],
                "odp_bias_parallel_hartree": odp["energy_parallel"],
                "odp_bias_perpendicular_hartree": odp["energy_perpendicular"],
                "odp_bias_hartree": odp["energy"],
            })
        if self.droplet:
            s = self.droplet_settings
            out = np.zeros_like(coords)
            energy = np.zeros(1)
            penetration = np.zeros(1)
            active = np.zeros(1, dtype=np.int64)
            status = int(oqp.oqp_namd_droplet_boundary(
                len(coords), int(self._droplet_groups.max()),
                oqp.ffi.cast("double *", coords.ctypes.data),
                oqp.ffi.cast("double *", self.masses.ctypes.data),
                oqp.ffi.cast("int64_t *", self._droplet_groups.ctypes.data),
                oqp.ffi.cast("double *", self._droplet_center.ctypes.data),
                s["radius_angstrom"]*ANGSTROM_TO_BOHR,
                s["buffer_angstrom"]*ANGSTROM_TO_BOHR,
                s["force_constant_kcal_mol_angstrom2"]*KCALMOLANG2_TO_HARTREEBOHR2,
                s["max_penetration_angstrom"]*ANGSTROM_TO_BOHR,
                oqp.ffi.cast("double *", energy.ctypes.data),
                oqp.ffi.cast("double *", out.ctypes.data),
                oqp.ffi.cast("double *", penetration.ctypes.data),
                oqp.ffi.cast("int64_t *", active.ctypes.data)))
            if status == 1:
                raise RuntimeError(
                    "droplet boundary failsafe: maximum penetration "
                    f"{penetration[0]*BOHR_TO_ANGSTROM:.6f} angstrom exceeds "
                    f"{s['max_penetration_angstrom']:.6f} angstrom")
            if status != 0:
                raise RuntimeError(
                    f"native droplet boundary rejected coordinates (status={status})")
            force += out
            record["droplet_energy_hartree"] = float(energy[0])
        if self.solute_com:
            out = np.zeros_like(coords)
            energy = np.zeros(1)
            displacement = np.zeros(1)
            status = int(oqp.oqp_namd_com_restraint(
                len(coords), oqp.ffi.cast("double *", coords.ctypes.data),
                oqp.ffi.cast("double *", self.masses.ctypes.data),
                oqp.ffi.cast("int64_t *", self._solute_selected.ctypes.data),
                oqp.ffi.cast("double *", self._solute_center.ctypes.data),
                self.solute_settings["force_constant_kcal_mol_angstrom2"]
                * KCALMOLANG2_TO_HARTREEBOHR2,
                oqp.ffi.cast("double *", energy.ctypes.data),
                oqp.ffi.cast("double *", out.ctypes.data),
                oqp.ffi.cast("double *", displacement.ctypes.data)))
            if status != 0:
                raise RuntimeError(
                    f"native solute COM restraint rejected coordinates (status={status})")
            force += out
            record["solute_com_energy_hartree"] = float(energy[0])
        energy = (record["odp_energy_hartree"] + record["droplet_energy_hartree"]
                  + record["solute_com_energy_hartree"])
        self.last = record
        return energy, force, record

    def settings(self):
        out = {"odp": (self.odp.provenance() if self.odp is not None
                       else {"enabled": False})}
        out["droplet"] = (dict(self.droplet_settings, enabled=True)
                          if self.droplet else {"enabled": False})
        out["solute_com"] = (dict(self.solute_settings, enabled=True)
                             if self.solute_com else {"enabled": False})
        return out


def require_orthorhombic(cell, names):
    """Refuse a cell with off-diagonal lattice components: the CVs and the
    restraints image each Cartesian axis on its own, which is the minimum
    image only in a rectangular cell."""
    cell = np.asarray(cell, dtype=np.float64).reshape(3, 3)
    if np.any(np.abs(cell - np.diag(np.diag(cell))) > 1.0e-12):
        raise ValueError(
            "%s in periodic QM/MM MD requires an orthorhombic cell; this cell "
            "has off-diagonal lattice components" % ", ".join(names))


def frozen_coordinates_digest(positions_nm, frozen_atoms):
    """SHA-256 of the frozen atoms' positions rounded to 1e-5 nm."""
    digest = hashlib.sha256()
    xyz = np.asarray(positions_nm, dtype=np.float64)
    for atom in sorted(int(i) for i in frozen_atoms):
        x, y, z = xyz[atom]
        digest.update(f"{atom}:{x:.5f},{y:.5f},{z:.5f};".encode())
    return digest.hexdigest()


def require_whole_group(topology, positions, box, atoms, label):
    """Refuse a group whose bonds cross the cell as given: a bond longer
    than half a box length along an axis means the group starts split
    between periodic images."""
    natom = int(topology.getNumAtoms())
    if isinstance(atoms, str):
        selected = set(_index_list(atoms, natom, label))
    else:
        selected = {int(i) for i in atoms}
    positions = np.asarray(positions, dtype=np.float64)
    box = np.asarray(box, dtype=np.float64).reshape(3)
    for bond in topology.bonds():
        i, j = int(bond[0].index), int(bond[1].index)
        if i in selected and j in selected and np.any(
                np.abs(positions[i] - positions[j]) > 0.5*box):
            raise ValueError(
                f"{label} atoms {i} and {j} are bonded across the periodic "
                "cell in the input geometry; write the solute whole")


def _forcefield_path(name):
    """The file OpenMM's ForceField loads for this name: a path, or a file
    shipped in OpenMM's data directories."""
    if os.path.isfile(name):
        return name
    from openmm.app import forcefield as openmm_forcefield
    for directory in openmm_forcefield._getDataDirectories():
        candidate = os.path.join(directory, name)
        if os.path.isfile(candidate):
            return candidate
    raise FileNotFoundError(f"force field file {name!r} not found")


def mm_hamiltonian_identity(forcefield_files, options):
    """Content digests of the force-field files and the QM/MM Hamiltonian
    options, so windows sampled with different potentials are not combined.

    The files are recorded by content alone, in loading order (the given
    files, then their includes): a renamed copy is the same force field.
    An include is named inside its parent, whose digest therefore also
    fixes the dependency structure."""
    import xml.etree.ElementTree as ElementTree

    files = []
    seen = set()
    # follow <Include file="..."> transitively, resolved as OpenMM's
    # ForceField.loadFile resolves it (next to the including file when that
    # exists there, otherwise through the data directories)
    queue = [_forcefield_path(str(name)) for name in forcefield_files]
    while queue:
        path = queue.pop(0)
        real = os.path.realpath(path)
        if real in seen:
            continue
        seen.add(real)
        with open(real, "rb") as stream:
            data = stream.read()
        files.append(hashlib.sha256(data).hexdigest())
        for included in ElementTree.fromstring(data).findall("Include"):
            name = included.attrib["file"]
            joined = os.path.join(os.path.dirname(path), name)
            queue.append(joined if os.path.isfile(joined)
                         else _forcefield_path(name))
    return {"forcefield_files": files,
            "options": {k: options[k] for k in sorted(options)}}


def _file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def qm_surface_identity(config, input_file=None):
    """Every QM option that changes the energy surface, as NAMD records it
    (input method/basis/library/ispher/..., scf, dftgrid, pcm, ...), with
    ``file:`` basis definitions replaced by the SHA-256 of their contents."""
    from oqp.library.namd import _electronic_config_identity

    identity = _electronic_config_identity(config)
    method = str(identity.get("input", {}).get("method", "")).lower()
    if isinstance(identity.get(method), dict):
        # backend artifacts are content-addressed by tight_binding_identity;
        # their path spellings would only make identical runs differ
        for key in ("parameter_path", "library_path", "executable"):
            identity[method].pop(key, None)
    directory = os.path.dirname(input_file) if input_file else ""
    for section, key in (("input", "basis"), ("scf", "init_basis")):
        value = identity.get(section, {}).get(key)
        if isinstance(value, str) and value.startswith("file:"):
            path = os.path.join(directory, value[len("file:"):])
            # by content alone: a copy under another name is the same basis
            identity[section][key] = (
                {"sha256": _file_digest(path)} if os.path.isfile(path)
                else {"missing": value})
    return identity


def espf_environment_identity():
    """The ESPF_* environment controls that change electrostatic QM/MM
    energies and forces, as the NAMD identity records them."""
    from oqp.library.namd import NAMD

    return NAMD._espf_environment_identity()


def tight_binding_identity(mol):
    """Content-addressed DFTB/xTB Hamiltonian (parameter directory, native
    library or probe executable, resolved like the calculation resolves
    them), through the NAMD implementation; None for other methods."""
    import types

    from oqp.library.namd import NAMD

    config = getattr(mol, "config", None) or {}
    method = str(config.get("input", {}).get("method", "")).lower()
    if method not in ("dftb", "xtb"):
        return None
    view = types.SimpleNamespace(mol=mol)
    # every NAMD method the resolution calls on ``self`` (the configuration
    # view first: QM/MM MD does not switch response manifolds, so it is
    # the molecule's configuration as given)
    for name in ("_signature_config", "_effective_tight_binding_settings",
                 "_resolved_tight_binding_artifacts", "_tight_binding_identity"):
        setattr(view, name, types.MethodType(getattr(NAMD, name), view))
    identity = view._tight_binding_identity()
    # copies of the same artifacts in other directories are the same
    # Hamiltonian: keep the content digests, not the resolved paths
    settings = {}
    for key, value in identity["settings"].items():
        if isinstance(value, dict) and "sha256" in value:
            value = {"sha256": value["sha256"]}
        elif isinstance(value, dict) and "resolved" in value:
            value = {"missing": os.path.basename(str(value["resolved"]))}
        settings[key] = value
    return {"section": identity["section"], "settings": settings}


def topology_identity(topology, qm_atoms, electronic):
    """Hash of what makes windows comparable: atoms, bonds, QM region, model."""
    digest = hashlib.sha256()
    for atom in topology.atoms():
        element = getattr(getattr(atom, "element", None), "symbol", "")
        digest.update(f"{atom.residue.name}|{atom.name}|{element};".encode())
    for bond in topology.bonds():
        digest.update(f"{bond[0].index}-{bond[1].index};".encode())
    digest.update(json.dumps(sorted(int(i) for i in qm_atoms)).encode())
    digest.update(json.dumps(electronic, sort_keys=True).encode())
    return {"kind": "qmmm_md", "natom": int(topology.getNumAtoms()),
            "sha256": digest.hexdigest(), "electronic": electronic}
