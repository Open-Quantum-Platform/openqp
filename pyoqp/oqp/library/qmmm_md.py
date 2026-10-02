import openmm.app as app
import openmm as mm
import openmm.unit as unit
import numpy as np
import os
import re
import time
from copy import deepcopy
import sys
from oqp.library.qmmm_active import (
    freeze_constrained_partners, held_atoms, resolve_active_set,
    selection_requested)
from oqp.library.qmmm_driver import OpenQpQMMM, read_xyz, is_periodic_method
from oqp.utils.md_controls import (
    continuation_seed,
    explicit_md_control_keys,
    merge_explicit_md_controls,
    openmm_random_seed,
    resolve_output_name,
    sectioned_md_keys,
)
from oqp.utils.md_snapshot import (
    check_snapshot_matches,
    numbered_snapshot_path,
    read_snapshot,
    write_snapshot,
)


def _copy_virtual_site(site):
    """A new OpenMM virtual site with the same definition (a System owns the
    sites it is given, so sys0's cannot be shared)."""
    name = type(site).__name__
    p = [site.getParticle(k) for k in range(site.getNumParticles())]
    if name == "TwoParticleAverageSite":
        return mm.TwoParticleAverageSite(p[0], p[1], site.getWeight(0), site.getWeight(1))
    if name == "ThreeParticleAverageSite":
        return mm.ThreeParticleAverageSite(p[0], p[1], p[2], site.getWeight(0), site.getWeight(1), site.getWeight(2))
    if name == "OutOfPlaneSite":
        return mm.OutOfPlaneSite(p[0], p[1], p[2], site.getWeight12(), site.getWeight13(), site.getWeightCross())
    raise NotImplementedError(f"virtual site type {name} is not supported by the QM/MM MD driver")


def _rigid_water_constraints(forcefield, topology, qm_atoms):
    """(i, j, distance) water constraints of an OpenMM rigidWater system of
    this topology, QM atoms excluded.  Built without a cutoff: constraints do
    not depend on the nonbonded treatment, and a periodic reference would tie
    this lookup to OpenMM's default 1 nm cutoff, which a box shorter than
    2 nm cannot hold."""
    ref = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff,
                                  constraints=None, rigidWater=True)
    qm = set(int(i) for i in qm_atoms)
    out = []
    for k in range(ref.getNumConstraints()):
        p1, p2, dist = ref.getConstraintParameters(k)
        if p1 in qm or p2 in qm:
            continue
        out.append((p1, p2, dist))
    return out


def _to_kJmol(energy):
    """Energy from the force backend (a Quantity or a bare float already in
    kJ/mol) as a float in kJ/mol."""
    if unit.is_quantity(energy):
        return float(energy.value_in_unit(unit.kilojoules_per_mole))
    return float(energy)


class _PDBTrajectoryReporter:
    """Multi-model PDB trajectory that is complete on disk after every run().

    ``app.PDBReporter`` writes its END record only when the object is
    garbage-collected.  A driver object kept alive after its run -- the usual
    situation in a script that then restarts the trajectory -- would write
    that record later, at its old file offset, into the middle of the
    continued file.  This reporter writes the same header, models and footer
    with OpenMM's own writers, but closes the file at the end of run() and
    takes the END record back off if the driver is run again.

    Every model is preceded by a ``REMARK OPENQP STEP n`` record.  A file need
    not start at step 0 (a restart that names a new trajectory file starts it
    at the checkpoint), so the step of a frame cannot be inferred from its
    position in the file.

    ``keep_through_step=None`` starts a new file; an integer continues an
    existing one, keeping its header and the models up to that step (frames
    written after the checkpoint being continued are dropped).  Models of a
    file without step records are taken to start at step 0."""

    STEP_RECORD = "REMARK OPENQP STEP"

    def __init__(self, path, interval, keep_through_step=None):
        self._path = path
        self._interval = int(interval)
        self._topology = None
        if keep_through_step is None:
            self._out = open(path, "w")
            self._next_model = 1                 # OpenMM numbers models from 1
            self._header_written = False
            return
        with open(path, "r") as stream:
            lines = stream.readlines()
        kept, models, recorded = [], 0, None
        for line in lines:
            if line.startswith(self.STEP_RECORD):
                recorded = int(line[len(self.STEP_RECORD):].split()[0])
                continue                         # re-emitted with its model
            if line.startswith("MODEL"):
                step = models * self._interval if recorded is None else recorded
                if step > keep_through_step:
                    break
                if recorded is not None:
                    kept.append(f"{self.STEP_RECORD} {recorded}\n")
                recorded = None
                models += 1
            if line.startswith("END") and not line.startswith("ENDMDL"):
                continue
            if line.startswith("CONECT"):
                continue
            kept.append(line)
        with open(path, "w") as stream:
            stream.writelines(kept)
        self._out = open(path, "a")
        self._next_model = models + 1
        self._header_written = True
        self.kept_models = models

    def describeNextReport(self, simulation):
        steps = self._interval - simulation.currentStep % self._interval
        return (steps, True, False, False, False, None)

    def _reopen(self):
        """Take the END record back off and append again."""
        with open(self._path, "r") as stream:
            lines = [line for line in stream
                     if not (line.startswith("CONECT")
                             or (line.startswith("END") and not line.startswith("ENDMDL")))]
        with open(self._path, "w") as stream:
            stream.writelines(lines)
        self._out = open(self._path, "a")

    def report(self, simulation, state):
        self._topology = simulation.topology
        if self._out.closed:
            self._reopen()
        if not self._header_written:
            app.PDBFile.writeHeader(simulation.topology, self._out)
            self._header_written = True
        self._out.write(f"{self.STEP_RECORD} {int(simulation.currentStep)}\n")
        app.PDBFile.writeModel(simulation.topology, state.getPositions(),
                               self._out, self._next_model)
        self._next_model += 1
        self._out.flush()

    def close_appended(self):
        if self._out is not None and not self._out.closed:
            if self._topology is not None:
                app.PDBFile.writeFooter(self._topology, self._out)
            self._out.close()


# ======================================================================
#  INI parser
# ======================================================================

def parse_ini_to_config(filepath):
    """
    Parse an INI-style configuration file into a flat dictionary
    with ``section.key`` keys.

    Numeric values are auto-converted to int or float.
    Entries with empty values are skipped.
    """
    config = {}
    section = None
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith(";"):
                continue
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1]
                continue
            if "=" not in line or section is None:
                continue
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            if not value:
                continue
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass
            config[f"{section}.{key}"] = value
    return config


# ======================================================================
#  Config helpers
# ======================================================================

_CUTOFF_MAP = {
    "pme":                app.PME,
    "nocutoff":           app.NoCutoff,
    "cutoffnonperiodic":  app.CutoffNonPeriodic,
    "cutoffperiodic":     app.CutoffPeriodic,
    "ewald":              app.Ewald,
}

_TRAJ_REPORTERS = {
    "pdb": _PDBTrajectoryReporter,
    "dcd": app.DCDReporter,
}

AU_VELOCITY_TO_NM_PER_PS = 2187.69126364

_VALID_ENSEMBLES = ("nve", "nvt")   # npt: see the message in __init__


def _parse_int_list(value):
    """
    Convert *value* to a list of ints.

    Accepts:
      - a list / ndarray  -> returned as-is (cast to int)
      - an int            -> [value]
      - a string          -> ints separated by commas and/or whitespace, each
        item optionally a ``start-end`` range
        e.g. ``"0,1,2"``  ``"0 1 2"`` (the Python API's form)  ``"0-2"``  ``"0-3, 8 9"``
    """
    import re as _re
    if isinstance(value, (list, np.ndarray)):
        return [int(v) for v in value]
    if isinstance(value, (int, np.integer)):
        return [int(value)]
    out = []
    for item in _re.split(r"[,\s]+", str(value).strip()):
        if not item:
            continue
        if "-" in item[1:]:                       # a range; a leading '-' would be a sign
            a, b = item.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(item))
    return out


def _parse_str_list(value):
    """Force-field file list -> list of strings.

    Accepts a list, a comma-separated string (the legacy form), or a
    whitespace-separated string (the form the NAMD driver has always taken).
    A single path that contains spaces is kept whole when it names an
    existing file, so ``/data/my forcefield.xml`` still works.
    """
    import re as _re
    if isinstance(value, list):
        return value
    text = str(value).strip()
    if not text:
        return []
    if "," in text:
        return [s.strip() for s in text.split(",") if s.strip()]
    if os.path.exists(text):
        return [text]
    return [s for s in _re.split(r"\s+", text) if s]


def _resolve_cutoff(value):
    """Map a string (or pass-through object) to an OpenMM cutoff method."""
    if value is None:
        return app.PME
    if isinstance(value, str):
        key = value.strip().lower().replace("_", "").replace(".", "")
        if key not in _CUTOFF_MAP:
            raise ValueError(
                f"Unknown cutoff '{value}'. "
                f"Choices: {', '.join(_CUTOFF_MAP)}"
            )
        return _CUTOFF_MAP[key]
    return value


def _extract_qmmm_config(oqp_cfg=None, mol=None):
    """
    Return a plain dict of QM/MM parameters gathered from either
    a flat ``oqp_cfg`` dict (``"qmmm.key"``) **or** from
    ``mol.config["qmmm"]``.

    Also returns the *remaining* oqp_cfg entries (QM settings)
    stripped of the ``qmmm.*`` keys, or ``None`` when in mol-mode.
    """
    if mol is not None:
        qmmm = dict(mol.config.get("qmmm", {}))
        return qmmm, None

    qmmm = {}
    qm_cfg = {}
    for k, v in oqp_cfg.items():
        if k.startswith("qmmm."):
            qmmm[k[5:]] = v
        else:
            qm_cfg[k] = v
    return qmmm, qm_cfg


# ======================================================================
#  QMMM_MD
# ======================================================================

class QMMM_MD:
    """
    QM/MM Molecular Dynamics in the NVE and NVT ensembles.

    There is no NPT: the OpenMM system that integrates the motion carries only
    a force that is linear in the coordinates, so a barostat would have to
    evaluate the full QM/MM energy at every trial cell, and a QM/MM run is far
    too short to equilibrate a density.  The cell is equilibrated classically
    and handed over with ``md(snapshot=...)``.

    All parameters are read from the configuration. Provide **exactly one** of:

    *   ``oqp_cfg``  - a ``dict`` with flat ``"section.key"`` entries
        **or** a path to an INI file. QM/MM settings live under the
        ``[qmmm]`` section / ``"qmmm.*"`` keys.

    *   ``mol`` - a pre-built OpenQP mol object whose
        ``mol.config["qmmm"]`` contains the same keys.

    Recognised ``[qmmm]`` keys
    --------------------------
    pdb_file           : str            (required)
    forcefield_files   : str            comma-separated list
    qm_atoms           : str or list    ``"0,1,2"`` or ``"0-2"``
    cutoff             : str            PME | NoCutoff | Ewald | ...
    embedding          : str            mechanical | electrostatic
    n_steps            : int            default 1000
    timestep           : float (fs)     default 1.0
    temperature        : float (K)      default 300.0
    ensemble           : str            nve | nvt         (default nve; npt rejected)
    friction           : float (ps^-1)  default 1.0      (NVT/NPT only)
    pressure           : float (bar)    default 1.0      (NPT only)
    barostat_interval  : int            default 25       (NPT only)
    trajectory_format  : str            pdb | dcd        (default pdb)
    trajectory_file    : str            default qmmm_trajectory.<format>
    log_file           : str            default qmmm_trajectory.dat
    report_interval    : int            default 1
    energy_file        : str            default total_energy.npz
    qm_atoms_xyz       : str            optional XYZ file
    qm_list            : str or list    optional index mapping

    Energies are sampled at the START of every step, at the positions the
    QM/MM force was just computed for: ``E_pot`` is the QM/MM energy returned
    by the force backend at those positions (not OpenMM's first-order
    extrapolation ``E(r0) - F.(r - r0)``, which is what the linear
    CustomExternalForce evaluates to once the atoms have moved), and ``E_kin``
    is the time-centred kinetic energy at the same instant.  Row ``step = s``
    belongs to the geometry after ``s`` steps; the trajectory file starts with
    the initial geometry as frame 0 and ``run()`` samples the geometry after
    the last step too, so rows 0..n_steps and frames 0..n_steps pair up.

    ``rigidwater`` (true unless the deck sets it, the NAMD driver's behaviour)
    puts the MM water bond/angle constraints of an OpenMM rigidWater system on
    the MD system, so OpenMM's Verlet applies SHAKE/RATTLE to them; QM atoms
    are never constrained.  The command line hands this driver the deck itself
    (config mode), so an omitted key is seen as omitted; the schema default
    stays false because it also feeds the NAMD restart identity and the
    legacy qmmm.py builders.  Without it the stiff O-H stretch is integrated
    explicitly and at 0.5 fs the total energy of a solvated box fluctuates by
    ~5 kJ/mol per 1000 atoms (the same figure a pure-MM run gives).

    Saved observables (in ``energy_file`` as ``.npz``)
    --------------------------------------------------
    step          : MD step index
    time_ps       : simulation time (ps)
    E_pot         : potential energy (kJ/mol)
    E_kin         : kinetic energy (kJ/mol)
    E_tot         : total energy (kJ/mol)
    temperature   : instantaneous T from kinetic energy (K)
    volume_nm3    : box volume (nm^3), only populated for NPT
    """

    def __init__(self, oqp_cfg=None, mol=None):

        # ------ validate ---------------------------------------------------
        if oqp_cfg is None and mol is None:
            raise ValueError("Either 'oqp_cfg' or 'mol' must be provided.")
        if oqp_cfg is not None and mol is not None:
            raise ValueError(
                "'oqp_cfg' and 'mol' are mutually exclusive - provide only one."
            )

        # ------ resolve oqp_cfg from file if needed -----------------------
        self._deck_path = getattr(mol, "input_file", None) if mol is not None else None
        if isinstance(oqp_cfg, str):
            if not os.path.isfile(oqp_cfg):
                raise FileNotFoundError(f"Config file not found: {oqp_cfg}")
            self._deck_path = oqp_cfg
            oqp_cfg = parse_ini_to_config(oqp_cfg)

        # ------ one definition of the QM region is enough -----------------
        # (mol mode was reconciled by Molecule.load_config)
        if isinstance(oqp_cfg, dict):
            from oqp.utils.qm_selection import reconcile_qm_selection
            _system, _qm_atoms = reconcile_qm_selection(
                oqp_cfg.get("input.system"), oqp_cfg.get("qmmm.qm_atoms"),
                oqp_cfg.get("qmmm.pdb_file"))
            oqp_cfg = dict(oqp_cfg)
            if _system is not None:
                oqp_cfg["input.system"] = _system
            if _qm_atoms is not None:
                oqp_cfg["qmmm.qm_atoms"] = _qm_atoms

        # ------ split qmmm.* settings from QM settings --------------------
        qmmm_cfg, qm_cfg = _extract_qmmm_config(oqp_cfg=oqp_cfg, mol=mol)

        # Nuclear propagation belongs to [md] for both all-QM and QM/MM
        # dynamics.  Retain the older qmmm.n_steps/timestep/... spellings as
        # explicit overrides while allowing the concise md(...) qmmm(...)
        # composition to drive this established OpenMM integrator.
        if mol is not None:
            md_cfg = dict(mol.config.get("md", {}))
            # mol.config is filled with schema defaults; only the deck says
            # which [md] keys a sectioned input really set
            written_md_keys = sectioned_md_keys(getattr(mol, "input_file", None))
        else:
            md_cfg = {
                str(key).split(".", 1)[1]: value
                for key, value in (qm_cfg or {}).items()
                if str(key).startswith("md.")
            }
            # a deck read here, or a dict handed in, holds no defaults:
            # whatever [md] key is present was written by the user
            written_md_keys = set(md_cfg)
        explicit_md_controls = explicit_md_control_keys(md_cfg, written_md_keys)
        if explicit_md_controls and not str(md_cfg.get("common_control_keys", "") or "").strip():
            md_cfg["common_control_keys"] = ",".join(sorted(explicit_md_controls))
        qmmm_cfg = merge_explicit_md_controls(qmmm_cfg, md_cfg)

        # ------ continuation and hand-over ---------------------------------
        # restart=true continues THIS run from restart_file: same step counter,
        # outputs appended.  snapshot=<file> starts a NEW run from a stored
        # phase-space point (an equilibration run of this driver, or a
        # classical one).  Both replace the PDB coordinates, the cell and the
        # initial velocities.
        _true = ("1", "true", "yes", "on", "t")
        self.restart = str(md_cfg.get("restart", False)).strip().lower() in _true
        self.restart_file = (str(md_cfg.get("restart_file", "") or "").strip()
                             or "qmmm_md.restart.npz")
        self.snapshot_file = str(md_cfg.get("snapshot", "") or "").strip()
        self.restart_interval = int(md_cfg.get("restart_interval", 10) or 0)
        self.snapshot_interval = int(md_cfg.get("snapshot_interval", 0) or 0)
        if self.restart_interval < 0 or self.snapshot_interval < 0:
            raise ValueError(
                "[md] restart_interval and snapshot_interval must be >= 0")
        if self.restart and self.snapshot_file:
            raise ValueError(
                "[md] restart continues a run and snapshot starts a new one; "
                "give one of them, not both")
        # A restart is the same deck plus restart=true, so its velocity= line
        # is simply superseded by the checkpoint.  A snapshot starts a new run,
        # where two sources of velocities would be a contradiction.
        if self.snapshot_file and "velocity" in explicit_md_controls:
            raise ValueError(
                "[md] velocity cannot be combined with snapshot: the stored "
                "velocities are the starting velocities")
        self._start_state = None
        if self.restart:
            self._start_state = read_snapshot(self.restart_file)
        elif self.snapshot_file:
            self._start_state = read_snapshot(self.snapshot_file)

        # ------ extract QM/MM parameters with defaults --------------------
        pdb_file = qmmm_cfg.get("pdb_file")
        if pdb_file is None:
            raise ValueError("'qmmm.pdb_file' is required in the configuration.")

        self._pdb_path = pdb_file
        self.pdb = app.PDBFile(pdb_file)
        if self._start_state is not None:
            self._apply_start_state()
        # [qmmm] active_atoms / frozen_atoms / active_radius / active_from_pdb.
        # Resolved in _build_md_system, once the OpenMM system exists.
        self._selection_cfg = qmmm_cfg

        ff_files = _parse_str_list(qmmm_cfg.get("forcefield_files", ""))
        if not ff_files:
            raise ValueError(
                "'qmmm.forcefield_files' is required in the configuration."
            )
        self.forcefield = app.ForceField(*ff_files)

        qm_atoms_raw = qmmm_cfg.get("qm_atoms")
        if qm_atoms_raw is None:
            raise ValueError("'qmmm.qm_atoms' is required in the configuration.")
        self.qm_atoms = np.array(_parse_int_list(qm_atoms_raw))

        # Fall back to the documented [qmmm] cutoff default rather than a
        # second, different default here: a concise qmmm(...) call leaves
        # unset keys out of the lowered config, so the old "PME" fallback
        # asked OpenMM for periodic boundaries on a non-periodic PDB and
        # aborted ("Requested periodic boundary conditions for a Topology
        # that does not specify periodic box dimensions").  A periodic box
        # with no cutoff named is reported, since the schema default treats
        # the system as non-periodic.
        from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
        schema_cutoff = OQP_CONFIG_SCHEMA["qmmm"]["cutoff"]["default"]
        cutoff_named = qmmm_cfg.get("cutoff") not in (None, "")
        self.cutoff    = _resolve_cutoff(
            qmmm_cfg.get("cutoff") if cutoff_named else schema_cutoff)
        if not cutoff_named and not is_periodic_method(self.cutoff):
            try:
                periodic_box = self.pdb.topology.getPeriodicBoxVectors()
            except Exception:
                periodic_box = None
            if periodic_box is not None:
                print("[QM/MM] the PDB carries a periodic box but [qmmm] "
                      "cutoff is unset; using the %s default. Set "
                      "cutoff=PME for periodic electrostatics." % schema_cutoff)
        self.embedding = str(qmmm_cfg.get("embedding", "electrostatic"))
        self.frontier_scheme = str(qmmm_cfg.get("frontier_scheme", "none"))
        _et = qmmm_cfg.get("ewald_tol", None)
        self.ewald_tol = None if _et in (None, "", "none", "None") else float(_et)
        self.lj_switch = str(qmmm_cfg.get("lj_switch", "false")).strip().lower() in ("1", "true", "yes", "on")
        self.h_lj = str(qmmm_cfg.get("h_lj", "false")).strip().lower() in ("1", "true", "yes", "on")
        _w = qmmm_cfg.get("mm_charge_width", None)
        self.mm_charge_width = None if _w in (None, "", "none", "None", 0, 0.0, "0") else float(_w)
        self.n_steps   = int(qmmm_cfg.get("n_steps", 1000))
        self.timestep  = float(qmmm_cfg.get("timestep", 1.0)) * unit.femtoseconds
        self.temperature = float(qmmm_cfg.get("temperature", 300.0)) * unit.kelvin
        self.initial_temperature = float(qmmm_cfg.get(
            "initial_temperature", qmmm_cfg.get("temperature", 300.0)
        )) * unit.kelvin
        self.mo_reuse = str(md_cfg.get("mo_reuse", True)).strip().lower() in (
            "1", "true", "yes", "on", "t"
        )
        self.velocity_source = (
            str(md_cfg.get("velocity", "maxwell"))
            if "velocity" in explicit_md_controls else "maxwell"
        )
        if explicit_md_controls.intersection({"seed", "rng_stream"}):
            self.random_seed = openmm_random_seed(
                md_cfg.get("seed", 0), md_cfg.get("rng_stream", 1))
        else:
            self.random_seed = None

        # ------ ensemble settings -----------------------------------------
        self.ensemble = str(qmmm_cfg.get("ensemble", "nve")).lower()
        if self.ensemble == "npt":
            raise NotImplementedError(
                "ensemble=npt is not available: a QM/MM barostat would need "
                "the QM/MM energy at each barostat trial box, and a QM/MM run "
                "is too short to equilibrate a density anyway. Equilibrate "
                "the cell classically (pure MM, NPT) and start from it with "
                "md(snapshot=...).")
        if self.ensemble not in _VALID_ENSEMBLES:
            raise ValueError(
                f"Unknown ensemble '{self.ensemble}'. "
                f"Choices: {', '.join(_VALID_ENSEMBLES)}"
            )
        self.friction = float(qmmm_cfg.get("friction", 1.0)) / unit.picosecond


        # ------ trajectory format -----------------------------------------
        fmt = str(qmmm_cfg.get("trajectory_format", "pdb")).lower()
        if fmt not in _TRAJ_REPORTERS:
            raise ValueError(
                f"Unknown trajectory_format '{fmt}'. "
                f"Choices: {', '.join(_TRAJ_REPORTERS)}"
            )
        self.trajectory_format = fmt

        # resolve_output_name, not get(key, fallback): the [qmmm] output-name
        # schema defaults are empty strings, so a materialised config carries
        # the key with an empty value and the fallback never applies.
        default_traj = f"qmmm_trajectory.{self.trajectory_format}"
        self.trajectory_file = resolve_output_name(
            qmmm_cfg, "trajectory_file", default_traj)
        self.log_file        = resolve_output_name(
            qmmm_cfg, "log_file", "qmmm_trajectory.dat")
        self.report_interval = int(qmmm_cfg.get("report_interval", 1))
        if self.report_interval <= 0:
            # md(trajectory_interval=0) means "automatic": about every 10 fs
            self.report_interval = max(1, int(round(
                10.0 / self.timestep.value_in_unit(unit.femtoseconds))))
        self.rigidwater = str(qmmm_cfg.get("rigidwater", True)).strip().lower() in (
            "1", "true", "yes", "on")
        self.energy_file     = resolve_output_name(
            qmmm_cfg, "energy_file", "total_energy.npz")
        # everything else this run reads, or that is already being written
        if mol is not None:
            geometry = mol.config.get("input", {}).get("system")
        else:
            geometry = (oqp_cfg or {}).get("input.system")
        self._borrowed_files = [
            ("[md] velocity", os.path.expanduser(self.velocity_source)),
            ("[qmmm] qm_atoms_xyz", qmmm_cfg.get("qm_atoms_xyz")),
            ("[input] system", geometry),
            ("OpenQP log", getattr(mol, "log", None)),
        ]
        self._validate_output_paths(ff_files)

        # ------ optional XYZ override for QM positions --------------------
        qm_atoms_xyz = qmmm_cfg.get("qm_atoms_xyz")
        qm_list_raw  = qmmm_cfg.get("qm_list")
        qm_list = _parse_int_list(qm_list_raw) if qm_list_raw is not None else None

        if qm_atoms_xyz is not None:
            if self._start_state is not None:
                raise ValueError(
                    "[qmmm] qm_atoms_xyz cannot be combined with [md] restart "
                    "or snapshot: the stored positions are the starting positions")
            self._apply_xyz_positions(str(qm_atoms_xyz), qm_list)

        if self.restart and self._start_state["step"] >= self.n_steps:
            raise ValueError(
                f"[md] restart: {self.restart_file} is at step "
                f"{self._start_state['step']} and nstep is {self.n_steps}; nstep is "
                "the total length of the trajectory, so raise it to continue")

        # ------ store QM config / mol for the driver ----------------------
        # The outer runtype=md only selects this ground-state QM/MM MD driver
        # (see pyoqp dispatch); the QM subsystem itself runs a single-point
        # energy+gradient each step. Force the internal QM runtype to 'energy'
        # so the QM engine's input check accepts a config that arrived with
        # runtype=md (config mode; harmless when the key is absent).
        if isinstance(qm_cfg, dict):
            for _k in list(qm_cfg):
                if str(_k).split('.')[-1].strip().lower() == 'runtype':
                    qm_cfg[_k] = 'energy'
            # That rewrite hides the dynamics from Molecule.get_config, which
            # would otherwise default [scf] verbose to 0 and stop the SCF from
            # printing one MO coefficient table per step (see
            # Molecule._quiet_orbitals_in_dynamics).  This IS a dynamics run,
            # so apply the same default here; an explicit verbose >= 2 in the
            # deck still prints, and verbose = 0 was already silent.
            _vkeys = [k for k in qm_cfg
                      if str(k).split('.')[-1].strip().lower() == 'verbose'
                      and str(k).split('.')[0].strip().lower() in ('scf', 'qm_cfg')]
            if not _vkeys:
                qm_cfg['scf.verbose'] = '0'
            elif all(str(qm_cfg[k]).strip() == '1' for k in _vkeys):
                for k in _vkeys:
                    qm_cfg[k] = '0'
        self.oqp_cfg = qm_cfg
        self.mol     = mol

        # ------ internal state --------------------------------------------
        self.oqp_driver = None
        self.mm_systems = None
        self.simulation_md = None
        self.system_md = None
        self.qmmm_ext = None

        self._traj_data = {
            "step":        [],
            "time_ps":     [],
            "E_pot":       [],
            "E_kin":       [],
            "E_tot":       [],
            "temperature": [],
            "volume_nm3":  [],
        }

    # ------------------------------------------------------------------
    #  Ensemble label
    # ------------------------------------------------------------------

    def _ensemble_label(self):
        return {
            "nve": "NVE (Verlet, lagged QM/MM forces - expect energy drift)",
            "nvt": "NVT (Langevin)",
        }[self.ensemble]

    # ------------------------------------------------------------------
    #  XYZ override
    # ------------------------------------------------------------------

    def _validate_output_paths(self, forcefield_files):
        """No output may land on another output or on an input.

        Checked once, after every destination is resolved and before any file
        is opened.  The inputs are borrowed: a snapshot is typically the start
        of many trajectories, and writing a checkpoint, a trajectory, a log or
        an energy table over it (or over the PDB or a force-field file) would
        destroy it for all of them."""
        outputs = {
            "trajectory_file": self.trajectory_file,
            "log_file": self.log_file,
            "energy_file": self._energy_npz_path(),
            "restart_file": self.restart_file,
        }
        resolved = {label: os.path.realpath(path) for label, path in outputs.items()}
        owner = {}
        for label, path in resolved.items():
            if path in owner:
                raise ValueError(
                    f"QM/MM MD outputs must be distinct files: {owner[path]} and "
                    f"{label} both name {outputs[label]!r}")
            owner[path] = label

        inputs = [("[qmmm] pdb_file", self._pdb_path)]
        if self.snapshot_file:
            inputs.append(("[md] snapshot", self.snapshot_file))
        # every force-field file, not only the last one named
        inputs.extend(("[qmmm] forcefield file", str(name))
                      for name in forcefield_files if os.path.isfile(str(name)))
        inputs.append(("input deck", getattr(self, "_deck_path", None)))
        inputs.extend(getattr(self, "_borrowed_files", ()))
        # keywords such as velocity=maxwell or system="pdb 1-4" name no file
        inputs = [(label, str(path)) for label, path in inputs
                  if path and os.path.isfile(str(path))]
        for input_label, input_path in inputs:
            hit = owner.get(os.path.realpath(input_path))
            if hit is not None:
                raise ValueError(
                    f"{hit}={outputs[hit]!r} would overwrite the input "
                    f"{input_label} ({input_path!r}); give the output another name")

        if self.snapshot_file and self.snapshot_interval:
            # a numbered snapshot of this run's own family is an output too
            number = re.search(r"\.snapshot\.(\d+)\.npz$", self.snapshot_file)
            if number is not None and os.path.realpath(numbered_snapshot_path(
                    self.restart_file, int(number.group(1)))) == os.path.realpath(
                        self.snapshot_file):
                raise ValueError(
                    f"[md] snapshot={self.snapshot_file!r} is one of the numbered "
                    "snapshots this run would write, which would overwrite the "
                    "starting point. Give restart_file a different name.")

    def _apply_start_state(self):
        """Replace the PDB coordinates (and the cell of a periodic system) by
        those of the restart checkpoint or the snapshot, before any OpenMM
        context is built from them."""
        state = self._start_state
        natom = self.pdb.topology.getNumAtoms()
        if state["natom"] != natom:
            raise ValueError(
                f"{state['path']} has {state['natom']} atoms, the PDB has {natom}")
        self.pdb.positions = unit.Quantity(
            [mm.Vec3(*row) for row in state["positions_nm"]], unit.nanometer)
        if state["box_nm"] is not None:
            lx, ly, lz = (float(v) for v in state["box_nm"])
            self.pdb.topology.setPeriodicBoxVectors(unit.Quantity(
                [mm.Vec3(lx, 0.0, 0.0), mm.Vec3(0.0, ly, 0.0), mm.Vec3(0.0, 0.0, lz)],
                unit.nanometer))

    def _apply_xyz_positions(self, xyz_path, qm_list=None):
        """Overwrite PDB positions for QM atoms from an XYZ file."""
        symbols, xyz_coords = read_xyz(xyz_path)

        if qm_list is None:
            qm_list = list(range(len(self.qm_atoms)))
        qm_list = np.asarray(qm_list, dtype=int)

        if len(qm_list) != len(self.qm_atoms):
            raise ValueError(
                f"qm_list length ({len(qm_list)}) != "
                f"qm_atoms length ({len(self.qm_atoms)})"
            )
        if np.any(qm_list >= len(symbols)):
            raise ValueError(
                f"qm_list index >= atoms in XYZ file ({len(symbols)})"
            )

        self.pdb.positions = deepcopy(self.pdb.positions)
        ang_to_nm = 0.1
        for k, pdb_idx in enumerate(self.qm_atoms):
            xyz_idx = qm_list[k]
            x, y, z = xyz_coords[xyz_idx] * ang_to_nm
            self.pdb.positions[pdb_idx] = mm.Vec3(x, y, z) * unit.nanometer

    # ------------------------------------------------------------------
    #  Setup
    # ------------------------------------------------------------------

    def _build_oqp_driver(self):
        self.oqp_driver = OpenQpQMMM(
            positions=self.pdb.positions,
            topology=self.pdb.topology,
            forcefield=self.forcefield,
            qm_atoms=self.qm_atoms,
            oqp_cfg=self.oqp_cfg,
            mol=self.mol,
            Cutoff=self.cutoff,
            Embedding=self.embedding,
            frontier_scheme=self.frontier_scheme,
            ewald_tol=self.ewald_tol,
            lj_switch=self.lj_switch,
            h_lj=self.h_lj,
            mm_charge_width=self.mm_charge_width,
        )
        # The first electronic evaluation must build a fresh guess.  The
        # helper below enables reuse only after that evaluation succeeds.
        self.oqp_driver._reuse_orbitals = False
        self.mm_systems = self.oqp_driver.mm_systems

    def _compute_qmmm_force(self, positions):
        """Evaluate one force and enable orbital reuse only after success."""
        try:
            result = self.oqp_driver.compute_force(
                positions, self.pdb.topology, self.mm_systems, self.qm_atoms
            )
        except Exception:
            # A failed SCF may leave unconverged orbitals resident.  A caller
            # that elects to resume must rebuild a fresh guess.
            self.oqp_driver._reuse_orbitals = False
            raise
        self.oqp_driver._reuse_orbitals = self.mo_reuse
        return result

    def _resolve_frozen_atoms(self):
        """The atoms OpenMM holds fixed for this run, from the ``[qmmm]``
        selection keys (see ``oqp.library.qmmm_active`` for the syntax).

        Empty when the deck asks for nothing, so an existing run propagates
        every atom as it always did.  A rigid-water constraint may not tie a
        moving atom to a fixed one, so constrained partners are frozen together.
        """
        cfg = getattr(self, "_selection_cfg", {}) or {}
        if not selection_requested(cfg):
            return set()
        box = self.oqp_driver._box_lengths_bohr()
        active, frozen = resolve_active_set(
            cfg,
            self.pdb.topology,
            self.pdb.positions.value_in_unit(unit.angstrom),
            self.qm_atoms,
            box_ang=None if box is None else [b * 0.52917721067 for b in box],   # bohr -> angstrom
            default_all=True,          # dynamics propagates everything unless asked
            pdb_path=getattr(self, "_pdb_path", None),
        )
        if self.rigidwater:
            pairs = [(p1, p2) for p1, p2, _ in _rigid_water_constraints(
                self.forcefield, self.pdb.topology, self.qm_atoms)]
            active, frozen = freeze_constrained_partners(pairs, active, frozen)
        # Everything outside the active set is held, not merely what
        # frozen_atoms named: an active_atoms / active_radius / active_from_pdb
        # selection holds every atom it did not select.
        held = held_atoms(self.pdb.topology, active)
        natom = self.pdb.topology.getNumAtoms()
        print(f"[QM/MM MD] active atoms: {natom - len(held)} of {natom} propagated, "
              f"{len(held)} held fixed (their charges and forces still act)")
        return held

    def _build_md_system(self):
        sys0 = self.mm_systems["sys0"]
        self.system_md = mm.System()
        for i in range(sys0.getNumParticles()):
            self.system_md.addParticle(sys0.getParticleMass(i))
        # virtual sites (e.g. the TIP4P M site) must follow their parents, as in
        # sys0; OpenMM then also moves forces applied to them onto the parents
        for i in range(sys0.getNumParticles()):
            if sys0.isVirtualSite(i):
                self.system_md.setVirtualSite(i, _copy_virtual_site(sys0.getVirtualSite(i)))

        # [qmmm] active_atoms / frozen_atoms: OpenMM holds an atom in place by
        # giving it zero mass, and that is what "frozen" means for runtype=md.
        # The atom keeps its charge, its embedding field and its force
        # contribution -- it simply does not move.  With no selection every atom
        # moves, exactly as before.
        self.frozen_atoms = self._resolve_frozen_atoms()
        for i in sorted(self.frozen_atoms):
            self.system_md.setParticleMass(i, 0.0)
        # The integration system starts as an empty mm.System(), whose default
        # box is OpenMM's 2 nm cube, not the cell of this topology.  Give it
        # the real cell so the logged volume and DCD frames describe it.
        _cell = self.pdb.topology.getPeriodicBoxVectors()
        if _cell is not None:
            self.system_md.setDefaultPeriodicBoxVectors(*_cell)

        # MM rigid-water constraints (O-H, O-H, H-H per TIP3P water), as the
        # NAMD driver's _build_constraints: taken from a rigidWater system of
        # the same topology, QM atoms excluded.  The MM forces still come from
        # the flexible sys0, whose water bond/angle terms vanish at the
        # constrained geometry.
        self.n_constraints = 0
        if self.rigidwater:
            for p1, p2, dist in _rigid_water_constraints(self.forcefield, self.pdb.topology, self.qm_atoms):
                if p1 in self.frozen_atoms or p2 in self.frozen_atoms:
                    # OpenMM rejects a constraint on a massless particle, and a
                    # held water has nothing to constrain: both atoms are fixed
                    # (_resolve_frozen_atoms freezes constrained partners together).
                    continue
                self.system_md.addConstraint(p1, p2, dist)
                self.n_constraints += 1
            print(f"[QM/MM MD] rigid water: {self.n_constraints} MM constraints applied "
                  f"(SHAKE/RATTLE in the integrator); QM atoms unconstrained")

        self.qmmm_ext = mm.CustomExternalForce(
            "-grad_x*x - grad_y*y - grad_z*z + qmmm_energy - ecorr"
        )
        self.system_md.addForce(self.qmmm_ext)
        self.qmmm_ext.addPerParticleParameter("grad_x")
        self.qmmm_ext.addPerParticleParameter("grad_y")
        self.qmmm_ext.addPerParticleParameter("grad_z")

        qmmm_energy, qmmm_force = self._compute_qmmm_force(self.pdb.positions)
        n_particles = self.system_md.getNumParticles()
        self._qmmm_energy_kJ = _to_kJmol(qmmm_energy)
        qmmm_energy = qmmm_energy / n_particles
        self.qmmm_ext.addGlobalParameter("qmmm_energy", qmmm_energy)

        ecorr = 0.0 * unit.kilojoules_per_mole
        for i in range(n_particles):
            self.qmmm_ext.addParticle(i, qmmm_force[i])
            for d in range(3):
                ecorr -= (
                    self.pdb.positions[i][d]
                    * qmmm_force[i][d]
                    / unit.nanometer
                    * unit.kilojoules_per_mole
                )
        ecorr /= n_particles
        self.qmmm_ext.addGlobalParameter("ecorr", ecorr)

        self._force_xyz = None     # positions (nm) the installed force belongs to
        self._force_kj_nm = np.asarray(qmmm_force, dtype=float)
        self._sys0_masses = np.array([
            sys0.getParticleMass(i).value_in_unit(unit.dalton)
            for i in range(n_particles)])
        if self._start_state is not None:
            check_snapshot_matches(
                self._start_state, self._sys0_masses,
                label="[md] restart_file" if self.restart else "[md] snapshot")

    def _build_integrator(self):
        if self.ensemble == "nve":
            return mm.VerletIntegrator(self.timestep)
        # NVT and NPT both use Langevin for temperature control
        integrator = mm.LangevinMiddleIntegrator(
            self.temperature, self.friction, self.timestep
        )
        if self.random_seed is not None:
            seed = self.random_seed
            if self.restart:
                # a new integrator with the same seed would replay the noise
                # of the first segment; give the continuation its own stream
                seed = continuation_seed(seed, self._start_state["step"])
            self.integrator_seed = seed
            integrator.setRandomNumberSeed(seed)
        return integrator

    def _set_initial_velocities(self):
        source = self.velocity_source.strip().lower()
        context = self.simulation_md.context
        if source in {"zero", "none", "0"}:
            velocities = np.zeros((self.system_md.getNumParticles(), 3))
            context.setVelocities(
                velocities * (unit.nanometer / unit.picosecond))
            return
        if source in {"maxwell", "boltzmann", "random"}:
            if self.random_seed is None:
                context.setVelocitiesToTemperature(self.initial_temperature)
            else:
                context.setVelocitiesToTemperature(
                    self.initial_temperature, self.random_seed)
            return
        path = os.path.abspath(os.path.expanduser(self.velocity_source))
        if not os.path.isfile(path):
            raise ValueError(
                f"[md] velocity={self.velocity_source!r} is not zero/maxwell "
                "or a readable file"
            )
        velocities = np.loadtxt(path, dtype=float).reshape(
            (self.system_md.getNumParticles(), 3))
        if not np.all(np.isfinite(velocities)):
            raise ValueError("[md] velocity file contains NaN or infinity")
        context.setVelocities(
            velocities * AU_VELOCITY_TO_NM_PER_PS
            * (unit.nanometer / unit.picosecond))

    def _build_simulation(self):
        integrator = self._build_integrator()
        self.simulation_md = app.Simulation(
            self.pdb.topology, self.system_md, integrator
        )
        context = self.simulation_md.context
        context.setPositions(self.pdb.positions)
        start = self._start_state
        if start is None:
            self._set_initial_velocities()
            if self.velocity_source.strip().lower() not in {
                    "maxwell", "boltzmann", "random"}:
                # velocity=zero or a velocity file states the velocities AT
                # the starting positions.  (A Maxwell-Boltzmann draw is a
                # random sample either way and is left as OpenMM draws it.)
                given = np.asarray(context.getState(getVelocities=True).getVelocities(
                    asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond))
                self._install_onstep_velocities(given)
        else:
            # A restart resumes the integrator's own (leapfrog) velocities, so
            # an NVE trajectory continues as if it had never stopped; a
            # snapshot supplies the velocities at its positions.
            velocities = start["velocities_nm_ps"]
            # The integrator's own (half-step-behind) velocities continue a
            # run exactly, but only with the time step they were written for:
            # the offset is half a step of THAT step.  With another time step
            # (or a checkpoint that does not say), start from the velocities
            # at the positions and let the new step define the offset.
            dt_ps = self.timestep.value_in_unit(unit.picoseconds)
            exact = (self.restart
                     and start["integrator_velocities_nm_ps"] is not None
                     and start.get("timestep_ps") is not None
                     and abs(start["timestep_ps"] - dt_ps) <= 1.0e-12 * dt_ps)
            if exact:
                velocities = start["integrator_velocities_nm_ps"]
            elif self.restart and start.get("timestep_ps") is not None:
                print(f"[QM/MM MD] restart: time step changed from "
                      f"{start['timestep_ps'] * 1000.0:g} fs to {dt_ps * 1000.0:g} fs; "
                      "continuing from the on-step velocities")
            from_snapshot = not exact
            if from_snapshot:
                self._install_onstep_velocities(velocities)
            else:
                velocities = np.array(velocities, dtype=float)
                for i in range(self.system_md.getNumParticles()):
                    if self.system_md.getParticleMass(i).value_in_unit(unit.dalton) <= 0.0:
                        velocities[i] = 0.0      # held atoms and virtual sites
                context.setVelocities(velocities * (unit.nanometer / unit.picosecond))
        step0 = 0
        if self.restart:
            step0 = int(start["step"])
            self.simulation_md.currentStep = step0
            context.setTime(float(start["time_ps"]) * unit.picoseconds)

        resumed = (self.restart and os.path.isfile(self.trajectory_file)
                   and self._resume_trajectory(step0))
        if not resumed:
            TrajReporter = _TRAJ_REPORTERS[self.trajectory_format]
            self.simulation_md.reporters.append(
                TrajReporter(self.trajectory_file, self.report_interval)
            )
            # frame 0 = the starting geometry, so that trajectory frame s and
            # energy row s (sampled before step s+1) describe the same structure
            state0 = context.getState(getPositions=True)
            for rep in self.simulation_md.reporters:
                rep.report(self.simulation_md, state0)

        # Energies are written by the driver itself (see ``_report_energies``):
        # OpenMM's StateDataReporter would report the potential of the linear
        # CustomExternalForce at the post-step positions, which is only a
        # first-order extrapolation of the QM/MM energy.
        self._log_columns = ['"Time (ps)"', '"Potential Energy (kJ/mole)"',
                             '"Kinetic Energy (kJ/mole)"',
                             '"Total Energy (kJ/mole)"', '"Temperature (K)"']
        if self.ensemble == "npt":
            self._log_columns.append('"Box Volume (nm^3)"')
        self._log_columns.append('"Speed (ns/day)"')
        if self.restart:
            # the energy table and the text log are independent files: a log
            # that was removed or renamed must not cost the table its history
            self._resume_energy_table(start)
        if self.restart and os.path.isfile(self.log_file):
            self._resume_text_log(start)
            self._log_handle = open(self.log_file, "a")
        else:
            self._log_handle = open(self.log_file, "w")
            self._log_handle.write("#" + ",".join(self._log_columns) + "\n")
            self._log_handle.flush()
        sys.stdout.write("#" + ",".join(['"Step"'] + self._log_columns) + "\n")
        sys.stdout.flush()
        self._wall_t0 = None
        if start is not None:
            print(f"[QM/MM MD] {'restart' if self.restart else 'snapshot'}: "
                  f"positions, velocities"
                  f"{'' if start['box_nm'] is None else ' and cell'} from "
                  f"{start['path']}"
                  + (f"; continuing at step {step0} of {self.n_steps}"
                     if self.restart else ""))

    def _install_onstep_velocities(self, velocities_nm_ps):
        """Hand velocities that are valid AT the current positions to the
        integrator.

        OpenMM's leapfrog integrators keep their velocities half a step behind
        the positions.  Installed unshifted, on-position velocities give the
        first step an extra half kick (v dt + a dt^2 instead of
        v dt + a dt^2 / 2) and a kinetic energy that is not the one asked for
        -- velocity=zero started with a nonzero kinetic energy.  Step them
        back by half a kick of the force just evaluated at these positions,
        then project onto the constraints."""
        context = self.simulation_md.context
        velocities = np.array(velocities_nm_ps, dtype=float)
        masses = np.array([
            self.system_md.getParticleMass(i).value_in_unit(unit.dalton)
            for i in range(self.system_md.getNumParticles())])
        moving = masses > 0.0
        velocities[~moving] = 0.0                # held atoms and virtual sites
        dt_ps = self.timestep.value_in_unit(unit.picoseconds)
        velocities[moving] -= (0.5 * dt_ps * self._force_kj_nm[moving]
                               / masses[moving, None])
        context.setVelocities(velocities * (unit.nanometer / unit.picosecond))
        context.applyVelocityConstraints(1.0e-8)

    # ------------------------------------------------------------------
    #  Restart: appending to the outputs of the run being continued
    # ------------------------------------------------------------------

    def _resume_trajectory(self, step0):
        """Append to the existing trajectory instead of replacing it.

        Frames written after the checkpoint (a run that stopped between two
        checkpoints) are dropped, so frames and steps stay aligned.  The file
        need not start at step 0: a restart that named a new trajectory file
        started it at that checkpoint, so the step of each frame is read from
        the file (the DCD header, the PDB step records), never assumed.

        Returns False when the file holds no frame at or before the
        checkpoint; the caller then starts it afresh, as for a missing file."""
        if self.trajectory_format == "dcd":
            with open(self.trajectory_file, "rb") as stream:
                header = stream.read(20)
            if len(header) < 20 or header[4:8] != b"CORD":
                raise ValueError(
                    f"[md] restart: {self.trajectory_file} is not a DCD "
                    "trajectory; name a new trajectory_file.")
            frames, first, interval = (
                int(v) for v in np.frombuffer(header[8:20], dtype="<i4"))
            interval = interval or self.report_interval
            # the first frame is written where the file was started, the
            # others on multiples of the interval
            keep = (1 + step0 // interval - first // interval
                    if step0 >= first else 0)
            if frames > keep:
                raise ValueError(
                    f"[md] restart: {self.trajectory_file} starts at step {first} "
                    f"and holds {frames} frames, but the checkpoint at step "
                    f"{step0} corresponds to {keep}; a DCD file cannot be "
                    "shortened in place. Restart from the final checkpoint, or "
                    "name a new trajectory_file.")
            if frames == 0:
                return False
            self.simulation_md.reporters.append(app.DCDReporter(
                self.trajectory_file, self.report_interval, append=True))
            return True
        reporter = _PDBTrajectoryReporter(
            self.trajectory_file, self.report_interval, step0)
        if reporter.kept_models == 0:
            reporter._out.close()
            return False
        self.simulation_md.reporters.append(reporter)
        return True

    def _resume_text_log(self, start):
        """Drop log rows written after the checkpoint, so the continued run
        appends rows instead of repeating or reordering them."""
        t0 = float(start["time_ps"])
        with open(self.log_file, "r") as stream:
            lines = stream.readlines()
        kept = []
        for line in lines:
            if line.startswith("#") or not line.strip():
                kept.append(line)
                continue
            try:
                if float(line.split(",")[0]) <= t0 + 1.0e-9:
                    kept.append(line)
            except ValueError:
                kept.append(line)
        with open(self.log_file, "w") as stream:
            stream.writelines(kept)

    def _resume_energy_table(self, start):
        """Reload the energy rows up to the checkpoint; the continued run then
        rewrites the table with the earlier rows still in it."""
        step0 = int(start["step"])
        energy_npz = self._energy_npz_path()
        if os.path.isfile(energy_npz):
            with np.load(energy_npz) as data:
                if all(key in data.files for key in self._traj_data):
                    mask = np.asarray(data["step"]) <= step0
                    for key in self._traj_data:
                        self._traj_data[key] = list(np.asarray(data[key])[mask])
        if self._traj_data["step"] and int(self._traj_data["step"][-1]) == step0:
            # the row of the checkpoint step is already on file, and the force
            # installed while building the system belongs to these positions
            xyz = np.asarray(self.simulation_md.context.getState(
                getPositions=True).getPositions(asNumpy=True).value_in_unit(
                    unit.nanometer))
            self._sampled_step, self._sampled_xyz = step0, xyz
            self._force_xyz = xyz

    # ------------------------------------------------------------------
    #  Checkpoints and snapshots
    # ------------------------------------------------------------------

    def _write_snapshot_file(self, path):
        """Positions, velocities and cell of the current state.

        The integrator keeps its velocities half a step behind the positions.
        Both are stored: those, for an exact continuation, and the velocities
        AT the positions (shifted by half a step of the current force), which
        is what a velocity-Verlet run such as surface hopping starts from."""
        context = self.simulation_md.context
        state = context.getState(getPositions=True, getVelocities=True)
        positions = np.asarray(state.getPositions(asNumpy=True).value_in_unit(
            unit.nanometer))
        leapfrog = np.asarray(state.getVelocities(asNumpy=True).value_in_unit(
            unit.nanometer / unit.picosecond))
        masses = np.array([
            self.system_md.getParticleMass(i).value_in_unit(unit.dalton)
            for i in range(self.system_md.getNumParticles())])
        dt_ps = self.timestep.value_in_unit(unit.picoseconds)
        onstep = leapfrog.copy()
        moving = masses > 0.0
        if self._force_is_current(positions):
            onstep[moving] += (0.5 * dt_ps * self._force_kj_nm[moving]
                               / masses[moving, None])
            if self.system_md.getNumConstraints():
                # The half-step kick uses the unconstrained force, so it adds
                # velocity along the rigid-water bonds.  Let OpenMM project it
                # out, then put the integrator's own velocities back untouched.
                quantity = unit.nanometer / unit.picosecond
                context.setVelocities(onstep * quantity)
                context.applyVelocityConstraints(1.0e-8)
                onstep = np.asarray(context.getState(getVelocities=True).getVelocities(
                    asNumpy=True).value_in_unit(quantity))
                context.setVelocities(state.getVelocities())
        box = (self._box_lengths_nm()
               if self.pdb.topology.getPeriodicBoxVectors() is not None else None)
        write_snapshot(
            path, positions_nm=positions, velocities_nm_ps=onstep,
            integrator_velocities_nm_ps=leapfrog, timestep_ps=dt_ps,
            masses_dalton=self._sys0_masses, box_nm=box,
            step=self.simulation_md.currentStep,
            time_ps=state.getTime().value_in_unit(unit.picoseconds),
            qm_atoms=self.qm_atoms, source="qmmm_md")

    def _write_periodic_snapshots(self, step_idx):
        if self.restart_interval and step_idx % self.restart_interval == 0:
            self._write_snapshot_file(self.restart_file)
        if (self.snapshot_interval and step_idx > 0
                and step_idx % self.snapshot_interval == 0):
            self._write_snapshot_file(
                numbered_snapshot_path(self.restart_file, step_idx))

    def setup(self):
        """Full setup: build driver, MD system, and simulation context."""
        self._build_oqp_driver()
        self._build_md_system()
        self._build_simulation()

    # ------------------------------------------------------------------
    #  MD loop helpers
    # ------------------------------------------------------------------

    def _update_qmmm_force(self, positions):
        qmmm_energy, qmmm_force = self._compute_qmmm_force(positions)
        self._install_qmmm_force(positions, qmmm_energy, qmmm_force)

    def _install_qmmm_force(self, positions, qmmm_energy, qmmm_force):
        """Hand an evaluated QM/MM energy and force to the integrator."""
        n_particles = self.system_md.getNumParticles()
        self._force_kj_nm = np.asarray(qmmm_force, dtype=float)
        self._qmmm_energy_kJ = _to_kJmol(qmmm_energy)
        qmmm_energy = qmmm_energy / n_particles
        self.simulation_md.context.setParameter("qmmm_energy", qmmm_energy)

        ecorr = 0.0 * unit.kilojoules_per_mole
        for i in range(n_particles):
            self.qmmm_ext.setParticleParameters(i, i, qmmm_force[i])
            for d in range(3):
                ecorr -= (
                    positions[i][d]
                    * qmmm_force[i][d]
                    / unit.nanometer
                    * unit.kilojoules_per_mole
                )
        ecorr /= n_particles
        self.simulation_md.context.setParameter("ecorr", ecorr)
        self.qmmm_ext.updateParametersInContext(self.simulation_md.context)

    def _instantaneous_temperature(self, E_kin_kJmol):
        """Compute T from kinetic energy: T = 2 * E_kin / (dof * k_B)."""
        # massless particles (virtual sites) carry no kinetic degrees of freedom
        massive = sum(1 for i in range(self.system_md.getNumParticles())
                      if self.system_md.getParticleMass(i).value_in_unit(unit.dalton) > 0.0)
        dof = 3 * massive - self.n_constraints
        if dof <= 0:
            return 0.0
        kB_kJ = unit.MOLAR_GAS_CONSTANT_R.value_in_unit(
            unit.kilojoules_per_mole / unit.kelvin
        )
        return (2.0 * E_kin_kJmol) / (dof * kB_kJ)

    # ------------------------------------------------------------------
    #  Force cache and cell
    # ------------------------------------------------------------------

    def _force_is_current(self, xyz_nm):
        cached = getattr(self, "_force_xyz", None)
        return cached is not None and np.array_equal(cached, xyz_nm)

    def _box_lengths_nm(self):
        vectors = self.pdb.topology.getPeriodicBoxVectors()
        return np.array([vectors[i][i].value_in_unit(unit.nanometer)
                         for i in range(3)], dtype=float)

    # ------------------------------------------------------------------
    #  Single-step interface
    # ------------------------------------------------------------------

    def step(self):
        """
        Perform one QM/MM MD step.

        Returns
        -------
        E_tot : float   Total energy (kJ/mol).
        """
        if self.simulation_md is None:
            self.setup()

        E_tot = self._sample_energy()

        self.simulation_md.step(1)

        state_md = self.simulation_md.context.getState(getPositions=True)
        pos0 = state_md.getPositions()

        sim0 = self.mm_systems["sim0"]
        sim0.context.setPositions(pos0)
        if is_periodic_method(self.cutoff):
            self.mm_systems["simew"].context.setPositions(pos0)
            self.mm_systems["simor"].context.setPositions(pos0)

        return E_tot

    def _sample_energy(self):
        """Refresh the QM/MM force at the current positions and record the
        Hamiltonian there (one energy row).  Returns E_tot (kJ/mol)."""
        sim0 = self.mm_systems["sim0"]

        # PR #205 review (M1c): update the QM/MM force at the CURRENT positions
        # BEFORE integrating, so the integrator applies a force consistent with the
        # positions it acts on. The previous order (step first, update after) left
        # the QM force one step stale and broke energy conservation.
        state_pre = self.simulation_md.context.getState(getPositions=True)
        pos_pre = state_pre.getPositions()
        xyz_now = np.asarray(state_pre.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        if (getattr(self, "_sampled_step", None) == self.simulation_md.currentStep
                and getattr(self, "_sampled_xyz", None) is not None
                and np.array_equal(self._sampled_xyz, xyz_now)):
            # a continued run() starts where the previous run() sampled its
            # last row: the force is already current and the row is written
            return self._traj_data["E_tot"][-1]
        if not self._force_is_current(xyz_now):
            # (a restart arrives with the force of its positions installed)
            sim0.context.setPositions(pos_pre)
            if is_periodic_method(self.cutoff):
                self.mm_systems["simew"].context.setPositions(pos_pre)
                self.mm_systems["simor"].context.setPositions(pos_pre)
            self._update_qmmm_force(pos_pre)
            self._force_xyz = xyz_now

        # Sample the Hamiltonian HERE, at the positions the force was computed
        # for: the linear term of the external force cancels identically at
        # pos_pre, so the potential is the QM/MM energy itself, and OpenMM's
        # kinetic energy is time-centred with the forces now in the context.
        state_energy = self.simulation_md.context.getState(getEnergy=True)
        E_pot = self._qmmm_energy_kJ
        E_kin = state_energy.getKineticEnergy().value_in_unit(
            unit.kilojoules_per_mole
        )
        E_tot = E_pot + E_kin
        T_inst = self._instantaneous_temperature(E_kin)

        # Box volume (only meaningful for periodic systems / NPT)
        if self.ensemble == "npt":
            box = state_pre.getPeriodicBoxVectors()
            vol = (box[0][0] * box[1][1] * box[2][2]).value_in_unit(
                unit.nanometer ** 3
            )
        else:
            vol = np.nan

        step_idx = self.simulation_md.currentStep
        t_ps = (step_idx * self.timestep).value_in_unit(unit.picoseconds)

        self._traj_data["step"].append(step_idx)
        self._traj_data["time_ps"].append(t_ps)
        self._traj_data["E_pot"].append(E_pot)
        self._traj_data["E_kin"].append(E_kin)
        self._traj_data["E_tot"].append(E_tot)
        self._traj_data["temperature"].append(T_inst)
        self._traj_data["volume_nm3"].append(vol)
        self._report_energies(step_idx, t_ps, E_pot, E_kin, E_tot, T_inst, vol)
        self._sampled_step, self._sampled_xyz = step_idx, xyz_now
        self._write_periodic_snapshots(step_idx)
        return E_tot

    def _report_energies(self, step_idx, t_ps, E_pot, E_kin, E_tot, T_inst, vol):
        """Write one energy row to ``log_file`` and to stdout (every
        ``report_interval`` steps, step 0 included)."""
        if step_idx % self.report_interval != 0:
            return
        if self._log_handle is None or self._log_handle.closed:
            self._log_handle = open(self.log_file, "a")     # a continued run() appends
        now = time.time()
        if self._wall_t0 is None:
            self._wall_t0 = (now, t_ps)
            speed = 0.0
        else:
            elapsed = now - self._wall_t0[0]
            speed = ((t_ps - self._wall_t0[1]) / 1000.0 * 86400.0 / elapsed
                     if elapsed > 0 else 0.0)
        cols = [f"{t_ps:.8g}", f"{E_pot:.14g}", f"{E_kin:.14g}",
                f"{E_tot:.14g}", f"{T_inst:.8g}"]
        if self.ensemble == "npt":
            cols.append(f"{vol:.8g}")
        cols.append(f"{speed:.3g}")
        self._log_handle.write(",".join(cols) + "\n")
        self._log_handle.flush()
        sys.stdout.write(",".join([str(step_idx)] + cols) + "\n")
        sys.stdout.flush()

    # ------------------------------------------------------------------
    #  Persistence
    # ------------------------------------------------------------------

    def _energy_npz_path(self):
        """The file the energy table is really written to.  ``np.savez``
        appends ``.npz`` to any other name, so ``energy_file="energies"`` lands
        in ``energies.npz``; a restart has to look for it under that name, or
        it would reload nothing and then overwrite the earlier rows."""
        base, ext = os.path.splitext(self.energy_file)
        if ext.lower() == ".npy":
            return base + ".npz"
        if ext.lower() != ".npz":
            return self.energy_file + ".npz"
        return self.energy_file

    def _save_traj_data(self):
        """Persist all collected per-step observables to a single .npz file."""
        # Written through an open handle: given a name, np.savez appends
        # ".npz" by its own (case-sensitive) rule, and "energies.NPZ" would
        # land in "energies.NPZ.npz" -- a file a restart never looks for.
        with open(self._energy_npz_path(), "wb") as stream:
            np.savez(
                stream,
                **{k: np.asarray(v) for k, v in self._traj_data.items()},
            )

    # ------------------------------------------------------------------
    #  Run all steps
    # ------------------------------------------------------------------

    def run(self):
        """
        Run the full simulation in the configured ensemble.

        Returns
        -------
        dict of np.ndarray
            Keys: step, time_ps, E_pot, E_kin, E_tot, temperature, volume_nm3.
        """
        try:
            if self.simulation_md is None:
                self.setup()

            print(f"\n\nStarting {self._ensemble_label()} dynamics:\n")

            n_todo = self.n_steps
            if (getattr(self, "restart", False)
                    and not getattr(self, "_restart_consumed", False)):
                # nstep is the total length; a restart runs what is left of it
                n_todo = self.n_steps - int(self._start_state["step"])
                self._restart_consumed = True
            for step_i in range(n_todo):
                self.step()

                # Persist on the same cadence as the trajectory reporters
                if (step_i + 1) % self.report_interval == 0:
                    self._save_traj_data()

            # the geometry after the last step gets its energy row too
            self._sample_energy()
            # Final save (covers n_steps not a multiple of report_interval)
            self._save_traj_data()
            # the end of the run is always a checkpoint
            self._write_snapshot_file(self.restart_file)
        finally:
            if getattr(self, "_log_handle", None) is not None:
                self._log_handle.close()
                self._log_handle = None
            for reporter in getattr(getattr(self, "simulation_md", None),
                                    "reporters", None) or ():
                finish = getattr(reporter, "close_appended", None)
                if finish is not None:
                    finish()                 # END record of a continued PDB

        return {k: np.asarray(v) for k, v in self._traj_data.items()}


# ======================================================================
#  Example usage
# ======================================================================
if __name__ == "__main__":

    # ---- Option A: single INI file contains everything -------------------
    #
    #   [input]
    #   functional = bhhlyp
    #   basis      = sto-3g
    #   method     = tdhf
    #
    #   [scf]
    #   type  = rohf
    #   maxit = 100
    #
    #   [tdhf]
    #   type   = mrsf
    #   nstate = 6
    #
    #   [properties]
    #   export    = true
    #   nac       = nacme
    #   back_door = true
    #   grad      = 5
    #
    #   [qmmm]
    #   pdb_file          = water_dimer.pdb
    #   forcefield_files  = tip3p.xml
    #   qm_atoms          = 0-2
    #   cutoff            = PME
    #   embedding         = electrostatic
    #   ensemble          = nvt          ; nve | nvt | npt
    #   friction          = 1.0          ; ps^-1
    #   pressure          = 1.0          ; bar (NPT only)
    #   barostat_interval = 25           ; steps (NPT only)
    #   trajectory_format = dcd
    #   n_steps           = 100
    #   timestep          = 1.0
    #   temperature       = 300.0
    #
    md = QMMM_MD(oqp_cfg="run.inp")
    data = md.run()

    # Post-run analysis:
    #   data = np.load("total_energy.npz")
    #   print(data["temperature"].mean(), data["E_tot"].std())
    #   if data["volume_nm3"][0] == data["volume_nm3"][0]:  # not NaN
    #       print("density-related volume:", data["volume_nm3"].mean())
