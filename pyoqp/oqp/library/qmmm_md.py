import openmm.app as app
import openmm as mm
import openmm.unit as unit
import numpy as np
import os
from copy import deepcopy
from sys import stdout
from oqp.library.qmmm_driver import OpenQpQMMM, read_xyz


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
    "pdb": app.PDBReporter,
    "dcd": app.DCDReporter,
}

_VALID_ENSEMBLES = ("nve", "nvt", "npt")


def _parse_int_list(value):
    """
    Convert *value* to a list of ints.

    Accepts:
      - a list / ndarray  -> returned as-is (cast to int)
      - an int            -> [value]
      - a string          -> comma-separated ints, or ``start-end`` range
        e.g. ``"0,1,2"``  or  ``"0-2"``  or  ``"0, 1, 2"``
    """
    if isinstance(value, (list, np.ndarray)):
        return [int(v) for v in value]
    if isinstance(value, (int, np.integer)):
        return [int(value)]
    value = str(value).strip()
    if "-" in value and "," not in value:
        parts = value.split("-")
        return list(range(int(parts[0]), int(parts[1]) + 1))
    return [int(v) for v in value.split(",")]


def _parse_str_list(value):
    """Comma-separated string or list -> list of stripped strings."""
    if isinstance(value, list):
        return value
    return [s.strip() for s in str(value).split(",") if s.strip()]


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
    QM/MM Molecular Dynamics with NVE / NVT / NPT support.

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
    ensemble           : str            nve | nvt | npt   (default nve)
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
        if isinstance(oqp_cfg, str):
            if not os.path.isfile(oqp_cfg):
                raise FileNotFoundError(f"Config file not found: {oqp_cfg}")
            oqp_cfg = parse_ini_to_config(oqp_cfg)

        # ------ split qmmm.* settings from QM settings --------------------
        qmmm_cfg, qm_cfg = _extract_qmmm_config(oqp_cfg=oqp_cfg, mol=mol)

        # ------ extract QM/MM parameters with defaults --------------------
        pdb_file = qmmm_cfg.get("pdb_file")
        if pdb_file is None:
            raise ValueError("'qmmm.pdb_file' is required in the configuration.")

        self.pdb = app.PDBFile(pdb_file)

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

        self.cutoff    = _resolve_cutoff(qmmm_cfg.get("cutoff", "PME"))
        self.embedding = str(qmmm_cfg.get("embedding", "electrostatic"))
        self.n_steps   = int(qmmm_cfg.get("n_steps", 1000))
        self.timestep  = float(qmmm_cfg.get("timestep", 1.0)) * unit.femtoseconds
        self.temperature = float(qmmm_cfg.get("temperature", 300.0)) * unit.kelvin

        # ------ ensemble settings -----------------------------------------
        self.ensemble = str(qmmm_cfg.get("ensemble", "nve")).lower()
        if self.ensemble not in _VALID_ENSEMBLES:
            raise ValueError(
                f"Unknown ensemble '{self.ensemble}'. "
                f"Choices: {', '.join(_VALID_ENSEMBLES)}"
            )
        self.friction = float(qmmm_cfg.get("friction", 1.0)) / unit.picosecond
        self.pressure = float(qmmm_cfg.get("pressure", 1.0)) * unit.bar
        self.barostat_interval = int(qmmm_cfg.get("barostat_interval", 25))

        if self.ensemble == "npt" and self.cutoff is app.NoCutoff:
            raise ValueError(
                "NPT requires a periodic cutoff method "
                "(PME / Ewald / CutoffPeriodic)."
            )

        # ------ trajectory format -----------------------------------------
        fmt = str(qmmm_cfg.get("trajectory_format", "pdb")).lower()
        if fmt not in _TRAJ_REPORTERS:
            raise ValueError(
                f"Unknown trajectory_format '{fmt}'. "
                f"Choices: {', '.join(_TRAJ_REPORTERS)}"
            )
        self.trajectory_format = fmt

        default_traj = f"qmmm_trajectory.{self.trajectory_format}"
        self.trajectory_file = str(qmmm_cfg.get("trajectory_file", default_traj))
        self.log_file        = str(qmmm_cfg.get("log_file",
                                                "qmmm_trajectory.dat"))
        self.report_interval = int(qmmm_cfg.get("report_interval", 1))
        self.energy_file     = str(qmmm_cfg.get("energy_file",
                                                "total_energy.npz"))

        # ------ optional XYZ override for QM positions --------------------
        qm_atoms_xyz = qmmm_cfg.get("qm_atoms_xyz")
        qm_list_raw  = qmmm_cfg.get("qm_list")
        qm_list = _parse_int_list(qm_list_raw) if qm_list_raw is not None else None

        if qm_atoms_xyz is not None:
            self._apply_xyz_positions(str(qm_atoms_xyz), qm_list)

        # ------ store QM config / mol for the driver ----------------------
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
            "npt": "NPT (Langevin + MonteCarloBarostat)",
        }[self.ensemble]

    # ------------------------------------------------------------------
    #  XYZ override
    # ------------------------------------------------------------------

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
        )
        self.mm_systems = self.oqp_driver.mm_systems

    def _build_md_system(self):
        sys0 = self.mm_systems["sys0"]
        self.system_md = mm.System()
        for i in range(sys0.getNumParticles()):
            self.system_md.addParticle(sys0.getParticleMass(i))

        self.qmmm_ext = mm.CustomExternalForce(
            "-grad_x*x - grad_y*y - grad_z*z + qmmm_energy - ecorr"
        )
        self.system_md.addForce(self.qmmm_ext)
        self.qmmm_ext.addPerParticleParameter("grad_x")
        self.qmmm_ext.addPerParticleParameter("grad_y")
        self.qmmm_ext.addPerParticleParameter("grad_z")

        qmmm_energy, qmmm_force = self.oqp_driver.compute_force(
            self.pdb.positions, self.pdb.topology, self.mm_systems, self.qm_atoms
        )
        n_particles = self.system_md.getNumParticles()
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

        # Add barostat for NPT
        if self.ensemble == "npt":
            barostat = mm.MonteCarloBarostat(
                self.pressure, self.temperature, self.barostat_interval
            )
            self.system_md.addForce(barostat)

    def _build_integrator(self):
        if self.ensemble == "nve":
            return mm.VerletIntegrator(self.timestep)
        # NVT and NPT both use Langevin for temperature control
        return mm.LangevinMiddleIntegrator(
            self.temperature, self.friction, self.timestep
        )

    def _build_simulation(self):
        integrator = self._build_integrator()
        self.simulation_md = app.Simulation(
            self.pdb.topology, self.system_md, integrator
        )
        self.simulation_md.context.setPositions(self.pdb.positions)
        self.simulation_md.context.setVelocitiesToTemperature(self.temperature)

        TrajReporter = _TRAJ_REPORTERS[self.trajectory_format]
        self.simulation_md.reporters.append(
            TrajReporter(self.trajectory_file, self.report_interval)
        )

        report_kwargs = dict(
            time=True, potentialEnergy=True,
            kineticEnergy=True, totalEnergy=True,
            temperature=True, speed=True,
        )
        if self.ensemble == "npt":
            report_kwargs["volume"] = True
            report_kwargs["density"] = True

        self.simulation_md.reporters.append(
            app.StateDataReporter(
                self.log_file, self.report_interval, **report_kwargs
            )
        )
        self.simulation_md.reporters.append(
            app.StateDataReporter(
                stdout, self.report_interval,
                step=True, **report_kwargs,
            )
        )

    def setup(self):
        """Full setup: build driver, MD system, and simulation context."""
        self._build_oqp_driver()
        self._build_md_system()
        self._build_simulation()

    # ------------------------------------------------------------------
    #  MD loop helpers
    # ------------------------------------------------------------------

    def _update_qmmm_force(self, positions):
        qmmm_energy, qmmm_force = self.oqp_driver.compute_force(
            positions, self.pdb.topology, self.mm_systems, self.qm_atoms
        )
        n_particles = self.system_md.getNumParticles()
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
        dof = 3 * self.system_md.getNumParticles()
        if dof <= 0:
            return 0.0
        kB_kJ = unit.MOLAR_GAS_CONSTANT_R.value_in_unit(
            unit.kilojoules_per_mole / unit.kelvin
        )
        return (2.0 * E_kin_kJmol) / (dof * kB_kJ)

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

        sim0 = self.mm_systems["sim0"]
        self.simulation_md.step(1)

        state_md = self.simulation_md.context.getState(getPositions=True)
        pos0 = state_md.getPositions()

        sim0.context.setPositions(pos0)
        if self.cutoff is not app.NoCutoff:
            self.mm_systems["simew"].context.setPositions(pos0)
            self.mm_systems["simor"].context.setPositions(pos0)

        self._update_qmmm_force(pos0)

        state_energy = self.simulation_md.context.getState(getEnergy=True)
        E_pot = state_energy.getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole
        )
        E_kin = state_energy.getKineticEnergy().value_in_unit(
            unit.kilojoules_per_mole
        )
        E_tot = E_pot + E_kin
        T_inst = self._instantaneous_temperature(E_kin)

        # Box volume (only meaningful for periodic systems / NPT)
        if self.ensemble == "npt":
            box = self.simulation_md.context.getState().getPeriodicBoxVectors()
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

        return E_tot

    # ------------------------------------------------------------------
    #  Persistence
    # ------------------------------------------------------------------

    def _save_traj_data(self):
        """Persist all collected per-step observables to a single .npz file."""
        base, ext = os.path.splitext(self.energy_file)
        out = base + ".npz" if ext == ".npy" else self.energy_file
        np.savez(
            out,
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
        if self.simulation_md is None:
            self.setup()

        print(f"\n\nStarting {self._ensemble_label()} dynamics:\n")

        for step_i in range(self.n_steps):
            self.step()

            # Persist on the same cadence as the trajectory reporters
            if (step_i + 1) % self.report_interval == 0:
                self._save_traj_data()

        # Final save (covers n_steps not a multiple of report_interval)
        self._save_traj_data()

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
