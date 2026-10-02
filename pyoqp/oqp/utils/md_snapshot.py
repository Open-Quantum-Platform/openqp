"""Phase-space snapshots of a full MD system: one small ``.npz`` file.

A snapshot is the complete nuclear state of a (QM + MM) system at one instant:
positions, velocities and, for a periodic system, the cell.  The same file
serves three purposes, which is why it is a format of its own rather than a
detail of one driver:

* the checkpoint ground-state QM/MM MD continues from (``md(restart=true)``);
* the starting point handed from an equilibration run to a production run --
  QM/MM MD or QM/MM surface hopping -- with ``md(snapshot="...")``;
* the hand-over from a classical (pure-MM) equilibration, written by a few
  lines of user script with :func:`write_snapshot`.

Units are OpenMM's, because every producer and consumer holds the system in an
OpenMM topology: nm, nm/ps, dalton, ps.  This module needs neither OpenMM nor
the OpenQP runtime.
"""
import os
import re
import tempfile

import numpy as np

SNAPSHOT_FORMAT = "openqp-md-snapshot-1"

#: 1 atomic unit of velocity in nm/ps (bohr / atomic time unit).
AU_VELOCITY_TO_NM_PER_PS = 2187.69126364
NM_TO_BOHR = 1.0 / 0.052917721067


def write_snapshot(path, *, positions_nm, velocities_nm_ps, masses_dalton,
                   box_nm=None, step=0, time_ps=0.0,
                   integrator_velocities_nm_ps=None, timestep_ps=None,
                   qm_atoms=None, source="user"):
    """Write one snapshot atomically (a reader never sees a partial file).

    ``velocities_nm_ps`` are the velocities AT the positions.  A leapfrog
    integrator keeps its own velocities half a step behind; a driver that
    wants to continue such a run exactly stores those as well, in
    ``integrator_velocities_nm_ps``, together with the ``timestep_ps`` they
    belong to: the half-step offset depends on the time step, so they are
    only meaningful to a run that uses the same one.  ``box_nm`` is the three orthorhombic
    cell lengths, or None for a non-periodic system.
    """
    positions = np.asarray(positions_nm, dtype=float)
    velocities = np.asarray(velocities_nm_ps, dtype=float)
    masses = np.asarray(masses_dalton, dtype=float).reshape(-1)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("snapshot positions must have shape (natom, 3)")
    natom = positions.shape[0]
    if velocities.shape != (natom, 3):
        raise ValueError("snapshot velocities must have shape (natom, 3)")
    if masses.shape != (natom,):
        raise ValueError("snapshot masses must have one entry per atom")
    payload = {
        "format": np.array(SNAPSHOT_FORMAT),
        "source": np.array(str(source)),
        "positions_nm": positions,
        "velocities_nm_ps": velocities,
        "masses_dalton": masses,
        "box_nm": (np.zeros(0) if box_nm is None
                   else np.asarray(box_nm, dtype=float).reshape(3)),
        "step": np.array(int(step), dtype=np.int64),
        "time_ps": np.array(float(time_ps)),
        "qm_atoms": (np.zeros(0, dtype=np.int64) if qm_atoms is None
                     else np.asarray(sorted(int(i) for i in qm_atoms), dtype=np.int64)),
    }
    if integrator_velocities_nm_ps is not None:
        leapfrog = np.asarray(integrator_velocities_nm_ps, dtype=float)
        if leapfrog.shape != (natom, 3):
            raise ValueError("integrator velocities must have shape (natom, 3)")
        payload["integrator_velocities_nm_ps"] = leapfrog
        if timestep_ps is not None:
            if not float(timestep_ps) > 0.0:
                raise ValueError("snapshot timestep_ps must be positive")
            payload["timestep_ps"] = np.array(float(timestep_ps))
    for key in ("positions_nm", "velocities_nm_ps", "masses_dalton", "box_nm"):
        if not np.all(np.isfinite(payload[key])):
            raise ValueError(f"snapshot {key} contains NaN or infinity")

    directory = os.path.dirname(os.path.abspath(path)) or "."
    os.makedirs(directory, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=".snapshot-", suffix=".npz", dir=directory)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            np.savez(stream, **payload)
        # mkstemp creates the file owner-only; give it the usual permissions
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(temporary, 0o666 & ~umask)
        os.replace(temporary, path)
    except BaseException:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise


def read_snapshot(path):
    """Read and validate a snapshot; returns a plain dict of arrays/scalars.

    ``box_nm`` is None for a non-periodic system and
    ``integrator_velocities_nm_ps`` is None when the writer did not store it.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"MD snapshot not found: {path}")
    try:
        with np.load(path, allow_pickle=False) as data:
            raw = {key: data[key] for key in data.files}
    except Exception as error:
        raise ValueError(f"{path} is not a readable MD snapshot: {error}") from error
    if str(raw.get("format", "")) != SNAPSHOT_FORMAT:
        raise ValueError(
            f"{path} is not an OpenQP MD snapshot (format "
            f"{str(raw.get('format', 'missing'))!r}, expected {SNAPSHOT_FORMAT!r}). "
            "A surface-hopping checkpoint is continued with namd(restart=true), "
            "not with snapshot=.")
    positions = np.asarray(raw["positions_nm"], dtype=float)
    velocities = np.asarray(raw["velocities_nm_ps"], dtype=float)
    masses = np.asarray(raw["masses_dalton"], dtype=float).reshape(-1)
    natom = positions.shape[0]
    if (positions.ndim != 2 or positions.shape[1] != 3
            or velocities.shape != (natom, 3) or masses.shape != (natom,)):
        raise ValueError(f"{path}: inconsistent snapshot array shapes")
    if not (np.all(np.isfinite(positions)) and np.all(np.isfinite(velocities))):
        raise ValueError(f"{path}: snapshot contains NaN or infinity")
    box = np.asarray(raw.get("box_nm", np.zeros(0)), dtype=float).reshape(-1)
    leapfrog = raw.get("integrator_velocities_nm_ps")
    return {
        "path": str(path),
        "source": str(raw.get("source", "")),
        "natom": natom,
        "positions_nm": positions,
        "velocities_nm_ps": velocities,
        "integrator_velocities_nm_ps": (
            None if leapfrog is None else np.asarray(leapfrog, dtype=float)),
        "timestep_ps": (float(raw["timestep_ps"]) if "timestep_ps" in raw else None),
        "masses_dalton": masses,
        "box_nm": box if box.size == 3 else None,
        "step": int(raw.get("step", 0)),
        "time_ps": float(raw.get("time_ps", 0.0)),
        "qm_atoms": np.asarray(raw.get("qm_atoms", np.zeros(0)), dtype=np.int64),
    }


def check_snapshot_matches(snapshot, masses_dalton, *, label="snapshot"):
    """Refuse a snapshot that belongs to a different system.

    Atom count and per-atom masses identify the system well enough to catch a
    snapshot of another topology or another atom order, which would otherwise
    start a trajectory from scrambled coordinates.  Massless rows (virtual
    sites, atoms a run holds fixed) are compared by count only.
    """
    masses = np.asarray(masses_dalton, dtype=float).reshape(-1)
    if snapshot["natom"] != masses.size:
        raise ValueError(
            f"{label} {snapshot['path']} has {snapshot['natom']} atoms, the "
            f"system has {masses.size}")
    stored = snapshot["masses_dalton"]
    both = (masses > 0.0) & (stored > 0.0)
    if not np.allclose(stored[both], masses[both], rtol=0.0, atol=1.0e-3):
        worst = int(np.argmax(np.abs(np.where(both, stored - masses, 0.0))))
        raise ValueError(
            f"{label} {snapshot['path']} does not describe this system: atom "
            f"{worst} has mass {stored[worst]:.4f} there and {masses[worst]:.4f} "
            "here (different topology or atom order)")


def numbered_snapshot_path(restart_file, step):
    """``run.restart.npz`` -> ``run.snapshot.00000100.npz`` for step 100."""
    base = str(restart_file)
    if base.lower().endswith(".npz"):
        base = base[:-4]
    if base.lower().endswith(".restart"):
        base = base[:-8]
    return f"{base}.snapshot.{int(step):08d}.npz"


def _directory_entry(path):
    """Where a write to ``path`` lands: the directory resolved through links,
    the last component not -- ``os.replace`` replaces a link standing there."""
    path = os.path.abspath(str(path))
    return os.path.join(os.path.realpath(os.path.dirname(path)),
                        os.path.basename(path))


def _link_chain(path, limit=40):
    """Every directory entry a name passes through, then its final referent."""
    current = _directory_entry(path)
    seen = set()
    while len(seen) < limit:
        yield current
        if current in seen or not os.path.islink(current):
            break
        seen.add(current)
        current = _directory_entry(
            os.path.join(os.path.dirname(current), os.readlink(current)))
    yield os.path.realpath(str(path))


def numbered_snapshot_step(restart_file, path):
    """The step whose numbered snapshot of ``restart_file`` ``path`` reaches.

    Every entry on the link chain of ``path`` is compared, not only its final
    referent: a snapshot is written by replacing the directory entry
    (``os.replace``), so a link standing at the destination is itself the file
    that is replaced, and a name that reaches the destination through it is
    redirected to the new snapshot afterwards.  Names are compared without
    regard to case, because on a case-insensitive filesystem (the macOS
    default) two spellings name one file; refusing a case-only variant on a
    case-sensitive one costs nothing.  Returns None when ``path`` never reaches
    one of the numbered snapshots this checkpoint name generates.
    """
    for entry in _link_chain(path):
        match = re.search(r"\.snapshot\.(\d+)\.npz$", os.path.basename(entry),
                          re.IGNORECASE)
        if match is None:
            continue
        step = int(match.group(1))
        destination = _directory_entry(numbered_snapshot_path(restart_file, step))
        if destination.casefold() == entry.casefold():
            return step
    return None
