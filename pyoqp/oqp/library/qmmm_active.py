"""Active / frozen atom selection shared by the QM/MM drivers.

``[qmmm] active_atoms`` names the atoms a QM/MM run is allowed to move, and
``frozen_atoms`` takes atoms back out of that set.  The spelling follows ORCA's
``%qmmm ActiveAtoms {0:5 16 21:30} end`` -- indices count from 0 and a range may
be written ``first:last`` -- and our own ``0-5,16,21-30`` form is accepted too,
so a selection can be copied from either side without editing.  ``name:`` groups
select by PDB atom name (``name:P,OP1,O5'``, a nucleic-acid backbone in one
line), ``active_radius`` adds whole MM residues near the QM region (ORCA's
``ActiveCore_Extension``), and ``active_from_pdb`` reads the selection from the
B-factor column of ``[qmmm] pdb_file``, 1 = active, as ORCA's
``Use_Active_InfoFromPDB`` does.

Defaults keep every existing deck running: an optimisation with no selection
moves the QM region only, dynamics with no selection propagates every atom.

Freezing changes what moves, never the physics: a frozen atom keeps its charge,
its embedding field and its force contribution.  Each driver holds it in its own
way -- the optimiser leaves the coordinates out of the optimisation vector, the
NAMD drivers do not integrate those rows, and ``runtype=md`` gives the particle
zero mass, which is how OpenMM holds an atom in place.
"""
import os

import numpy as np

__all__ = ["parse_atom_selection", "active_from_pdb_file", "residues_within_radius",
           "resolve_active_set", "selection_requested", "freeze_constrained_partners",
           "SELECTION_KEYS"]

#: the ``[qmmm]`` keys this module reads (ORCA's names)
SELECTION_KEYS = ("active_atoms", "frozen_atoms", "active_radius", "active_from_pdb")

_BREAK = str.maketrans({"{": " ", "}": " ", ";": " ", ",": " "})


def _truthy(value):
    return (value is True) or (str(value or "").strip().lower() in ("1", "true", "yes", "on", "t"))


def _radius_of(cfg):
    raw = cfg.get("active_radius", 0.0)
    text = "" if raw is None else str(raw).strip()
    if text.lower() in ("", "none"):
        return 0.0
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"[qmmm] active_radius must be a distance in angstrom, got {raw!r}.")


def selection_requested(qmmm_cfg):
    """True when the deck asks for a selection at all; otherwise each driver
    keeps its own default (the QM region, or every atom for dynamics)."""
    cfg = qmmm_cfg or {}
    return bool(str(cfg.get("active_atoms", "") or "").strip()
                or str(cfg.get("frozen_atoms", "") or "").strip()
                or _truthy(cfg.get("active_from_pdb", False))
                or _radius_of(cfg) > 0.0)


def _groups(text):
    """Split a selection into groups, keeping ``name:a,b,c`` groups whole."""
    out = []
    for chunk in str(text).translate(str.maketrans({"{": " ", "}": " ", ";": " "})).split():
        if chunk.lower().startswith("name:"):
            out.append(chunk)
        else:
            out.extend(part for part in chunk.translate(_BREAK).split() if part)
    return out


def parse_atom_selection(spec, atoms, key):
    """``0:5 16 21:30`` / ``0-5,16,21-30`` / ``name:P,OP1`` -> 0-based indices.

    ``atoms`` is the OpenMM topology atom list.  Virtual sites (element None)
    are never selected: OpenMM places them from their parents, so they are not
    independent coordinates.
    """
    text = str(spec or "").strip()
    if not text:
        return set()
    by_name = {}
    for atom in atoms:
        if atom.element is not None:
            by_name.setdefault(str(atom.name).strip().lower(), []).append(int(atom.index))
    natom, out = len(atoms), set()
    for group in _groups(text):
        if group.lower().startswith("name:"):
            names = [n for n in group[5:].split(",") if n.strip()]
            if not names:
                raise ValueError(f"[qmmm] {key}: 'name:' needs at least one atom name, got {group!r}.")
            for name in names:
                hit = by_name.get(name.strip().lower())
                if not hit:
                    raise ValueError(f"[qmmm] {key}: no atom is named {name.strip()!r} in [qmmm] pdb_file "
                                     "(names are matched case-insensitively).")
                out.update(hit)
            continue
        parts = group.replace(":", "-").split("-")
        try:
            if len(parts) == 1:
                lo = hi = int(parts[0])
            elif len(parts) == 2:
                lo, hi = int(parts[0]), int(parts[1])
            else:
                raise ValueError
        except ValueError:
            raise ValueError(f"[qmmm] {key}: {group!r} is neither a 0-based index, a range "
                             "(first-last or ORCA's first:last), nor a 'name:...' group.")
        if lo > hi:
            raise ValueError(f"[qmmm] {key}: range {group!r} runs backwards.")
        if lo < 0 or hi >= natom:
            raise ValueError(f"[qmmm] {key}: {group!r} is outside the 0-based range 0..{natom - 1} "
                             "of [qmmm] pdb_file.")
        out.update(i for i in range(lo, hi + 1) if atoms[i].element is not None)
    return out


def active_from_pdb_file(path, atoms, key="active_from_pdb"):
    """The B-factor column of ``[qmmm] pdb_file``: 1 = active, 0 = frozen.

    OpenMM's PDBFile drops that column, so the file is read again through the
    parser OpenMM uses internally (pdbstructure), which keeps it.
    """
    from openmm import unit as openmm_unit
    from openmm.app.internal.pdbstructure import PdbStructure
    if not path or not os.path.isfile(path):
        raise ValueError(f"[qmmm] {key}=true needs a readable [qmmm] pdb_file, got {path!r}.")
    with open(path) as handle:
        structure = PdbStructure(handle)

    def flag(atom):
        # pdbstructure hands the B-factor back as a Quantity in angstrom^2, not
        # a bare number, and float() refuses a Quantity outright
        value = getattr(atom, "temperature_factor", 0.0)
        if openmm_unit.is_quantity(value):
            value = value.value_in_unit(openmm_unit.angstrom ** 2)
        return float(value or 0.0)

    flags = [flag(atom)
             for chain in structure.iter_chains()
             for residue in chain.iter_residues()
             for atom in residue.iter_atoms()]
    if len(flags) != len(atoms):
        raise ValueError(f"[qmmm] {key}=true read {len(flags)} atoms from {os.path.basename(path)} but the "
                         f"topology has {len(atoms)}; it must be the file [qmmm] pdb_file names.")
    picked = {i for i, b in enumerate(flags) if b >= 0.5 and atoms[i].element is not None}
    if not picked:
        raise ValueError(f"[qmmm] {key}=true selected no atom: the B-factor column of "
                         f"{os.path.basename(path)} is 0 everywhere (1 marks an active atom).")
    return picked


def residues_within_radius(radius, positions_ang, qm_atoms, topology, box_ang=None):
    """QM atoms plus every MM residue with an atom within ``radius`` of one.

    The residue carrying the QM atoms is selected atom by atom (a covalent cut
    runs through it); every other residue is taken whole, which keeps waters and
    side chains intact.
    """
    qm = set(int(i) for i in qm_atoms)
    if radius <= 0.0:
        return set(qm)
    X = np.asarray(positions_ang, dtype=float)
    qx = X[sorted(qm)]
    box = None if box_ang is None else np.asarray(box_ang, dtype=float)

    def dist(idx):
        d = X[idx][:, None, :] - qx[None, :, :]
        if box is not None:
            d = d - box * np.round(d / box)
        return np.linalg.norm(d, axis=2).min(axis=1)     # per atom: nearest QM atom

    out = set(qm)
    for residue in topology.residues():
        # virtual sites (element None, e.g. the TIP4P M site) are not independent
        # coordinates: OpenMM places them from their parents
        idx = [a.index for a in residue.atoms() if a.element is not None]
        if not idx:
            continue
        if qm.intersection(idx):
            mm_here = [i for i in idx if i not in qm]
            if mm_here:
                out.update(int(i) for i, r in zip(mm_here, dist(mm_here)) if r <= radius)
        elif dist(idx).min() <= radius:
            out.update(idx)
    return out


def resolve_active_set(qmmm_cfg, topology, positions_ang, qm_atoms,
                       box_ang=None, default_all=False, pdb_path=None, radius=None):
    """Return ``(active_indices, frozen_indices)`` for one driver.

    ``default_all`` is what the driver does when the deck asks for nothing:
    dynamics propagates every atom, an optimisation moves the QM region only.
    The QM atoms are always active and may not be frozen.
    """
    cfg = dict(qmmm_cfg or {})
    if radius is not None:
        cfg["active_radius"] = radius
    atoms = list(topology.atoms())
    real = {i for i, a in enumerate(atoms) if a.element is not None}
    qm = set(int(i) for i in qm_atoms)

    frozen = parse_atom_selection(cfg.get("frozen_atoms", ""), atoms, "frozen_atoms")
    clash = sorted(frozen.intersection(qm))
    if clash:
        raise ValueError(f"[qmmm] frozen_atoms names QM atoms {clash} (0-based). The QM region is what "
                         "the calculation moves; freeze MM atoms only, or take those atoms out of "
                         "[qmmm] qm_atoms.")
    radius = _radius_of(cfg)
    if not np.isfinite(radius) or radius < 0.0:
        raise ValueError("[qmmm] active_radius must be a finite distance >= 0 angstrom, "
                         f"got {cfg.get('active_radius')!r}.")
    explicit = str(cfg.get("active_atoms", "") or "").strip()
    from_pdb = _truthy(cfg.get("active_from_pdb", False))

    if explicit or from_pdb or radius > 0.0:
        base = residues_within_radius(radius, positions_ang, qm, topology, box_ang)
        base |= parse_atom_selection(explicit, atoms, "active_atoms")
        if from_pdb:
            base |= active_from_pdb_file(pdb_path, atoms)
    else:
        # nothing selected: frozen_atoms alone still means "everything else moves"
        base = set(real) if default_all else set(qm)

    active = ((base - frozen) | qm) & real
    return np.array(sorted(active), dtype=int), frozen


def freeze_constrained_partners(pairs, active, frozen):
    """Hold the constrained partner of a frozen atom, and let a constrained
    partner of a moving atom move.

    A rigid-water constraint between a moving atom and a fixed one is not
    something SHAKE/RATTLE can satisfy -- it would drag the fixed atom.  So a
    pair that straddles the boundary is resolved: if either atom was explicitly
    frozen the whole pair freezes, otherwise the whole pair moves.  Returns
    ``(active, frozen)`` as sets.
    """
    active = set(int(i) for i in active)
    frozen = set(int(i) for i in frozen)
    changed = True
    while changed:
        changed = False
        for p1, p2 in pairs:
            p1, p2 = int(p1), int(p2)
            if (p1 in active) == (p2 in active):
                continue
            if frozen.intersection((p1, p2)):
                for p in (p1, p2):
                    if p in active:
                        active.discard(p)
                        frozen.add(p)
                        changed = True
            else:
                active.update((p1, p2))
                changed = True
    return active, frozen
