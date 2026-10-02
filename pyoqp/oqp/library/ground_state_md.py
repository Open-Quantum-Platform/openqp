"""Gas-phase ground-state Born--Oppenheimer molecular dynamics.

The driver evaluates the ground-state energy and analytic gradient at every
geometry and propagates the nuclei in atomic units with velocity Verlet.  It
shares the legacy ``[md]`` section with NAMD so that the concise input surface
can put common nuclear-propagation controls in ``md(...)``.
"""

from __future__ import annotations

import csv
import os
from datetime import date

import numpy as np

from oqp.library.single_point import Gradient, LastStep, SinglePoint
from oqp.periodic_table import ELEMENTS_NAME
from oqp.utils.file_utils import dump_log


FS_TO_AU = 41.341374575751
KB_HARTREE = 3.166811563e-6
AMU_TO_AU = 1822.888486209
BOHR_TO_ANGSTROM = 0.529177210903
_MISSING = object()


def _as_bool(value):
    return value is True or str(value).strip().lower() in {
        "1", "true", "yes", "on", "t",
    }


def _resolved_seed(value):
    seed = int(value)
    return seed if seed else int(date.today().strftime("%Y%m%d"))


class GroundStateMD:
    """All-QM ground-state BOMD using analytic gradients."""

    def __init__(self, mol):
        self.mol = mol
        manager = getattr(mol, "mpi_manager", None)
        if manager is not None and bool(getattr(manager, "use_mpi", False)):
            raise ValueError(
                "[md] ground-state all-QM MD currently requires a single "
                "MPI rank; rerun with --nompi"
            )
        md = mol.config["md"]
        self.natom = int(mol.data["natom"])
        self.mass = np.asarray(mol.get_mass(), dtype=float) * AMU_TO_AU
        if self.mass.shape != (self.natom,) or np.any(self.mass <= 0.0):
            raise ValueError("ground-state MD requires one positive mass per atom")

        self.nstep = int(md.get("nstep", 100))
        self.dt_fs = float(md.get("dt", 0.5))
        self.dt = self.dt_fs * FS_TO_AU
        if self.nstep < 0:
            raise ValueError("[md] nstep must be non-negative")
        if not np.isfinite(self.dt_fs) or self.dt_fs <= 0.0:
            raise ValueError("[md] dt must be finite and positive")

        self.temperature = float(md.get(
            "thermostat_temperature", md.get("init_temp", 300.0)))
        self.init_temp = float(md.get("init_temp", self.temperature))
        self.friction = float(md.get("thermostat_friction", 1.0))
        self.thermostat = str(md.get("thermostat", "off")).strip().lower()
        legacy_ensemble = str(md.get("ensemble", "nve")).strip().lower()
        if self.thermostat not in {"off", "langevin"}:
            raise ValueError("[md] thermostat must be off or langevin")
        if legacy_ensemble not in {"nve", "nvt"}:
            raise ValueError("[md] ensemble must be nve or nvt")
        if legacy_ensemble == "nvt" and self.thermostat == "off":
            self.thermostat = "langevin"
        if not np.isfinite(self.init_temp) or self.init_temp < 0.0:
            raise ValueError("[md] initial temperature must be finite and non-negative")
        if not np.isfinite(self.temperature) or self.temperature < 0.0:
            raise ValueError("[md] thermostat temperature must be finite and non-negative")
        if self.thermostat == "langevin" and (
                not np.isfinite(self.friction) or self.friction <= 0.0):
            raise ValueError("[md] friction must be positive for thermostat=langevin")

        self.seed = _resolved_seed(md.get("seed", 0))
        self.rng_stream = int(md.get("rng_stream", 1))
        if self.rng_stream < 0:
            raise ValueError("[md] rng_stream must be non-negative")
        mixed_seed = (
            (self.seed & ((1 << 64) - 1))
            ^ ((self.rng_stream * 0x9E3779B97F4A7C15) & ((1 << 64) - 1))
        )
        self.rng = np.random.Generator(np.random.PCG64(mixed_seed))
        self.velocity_source = str(md.get("velocity", "maxwell"))
        self.velocity = self._initial_velocity()

        # Checkpointing is a surface-hopping feature: this driver writes no
        # checkpoint and reads none, so an accepted-but-ignored restart request
        # would make fresh velocities, start at step zero and then truncate the
        # trajectory and energy files of the run it was asked to continue.
        # Refuse here, before the first electronic evaluation and before any
        # output file is opened.
        if _as_bool(md.get("restart", False)):
            raise ValueError(
                "[md] restart is not available for ground-state MD: this "
                "driver has no checkpoint to continue from. Remove restart, "
                "or run the trajectory from the start.")
        if str(md.get("restart_file", "") or "").strip():
            raise ValueError(
                "[md] restart_file is not available for ground-state MD: this "
                "driver writes no checkpoint. Remove restart_file; "
                "trajectory_file and energy_file carry the outputs.")

        # Output cadence of the xyz trajectory; the csv energy table keeps
        # every step.  0 is the documented automatic cadence, about 10 fs.
        self.trajectory_interval = int(md.get("trajectory_interval", 1) or 0)
        if self.trajectory_interval < 0:
            raise ValueError("[md] trajectory_interval must be >= 0")
        if self.trajectory_interval == 0:
            self.trajectory_interval = max(1, int(round(10.0 / self.dt_fs)))

        if (str(md.get("snapshot", "") or "").strip()
                or int(md.get("snapshot_interval", 0) or 0)):
            raise ValueError(
                "[md] snapshot / snapshot_interval belong to QM/MM dynamics "
                "(they store every atom of the embedded system); gas-phase "
                "ground-state MD takes its geometry from the input and its "
                "velocities from velocity=.")

        prefix = os.path.splitext(os.path.abspath(mol.log))[0]
        trajectory = str(md.get("trajectory_file", "") or "").strip()
        self.trajectory_file = trajectory or prefix + ".md.xyz"
        energy_file = str(md.get("energy_file", "") or "").strip()
        self.energy_file = energy_file or prefix + ".md.csv"
        if os.path.realpath(self.trajectory_file) == os.path.realpath(self.energy_file):
            raise ValueError("[md] trajectory_file and energy_file must be different")
        self._protect_inputs_from_outputs()

        self._thermostat_exchange = 0.0
        self._thermostat_exchange_cumulative = 0.0

    def _protect_inputs_from_outputs(self):
        """run() opens both outputs for writing; neither may be a file this
        calculation reads (geometry, velocities, the input deck) or its log."""
        mol = self.mol
        inputs = {"the log file": getattr(mol, "log", None)}
        for attribute in ("input_file", "oqp_input_source"):
            inputs["the input deck (%s)" % attribute] = getattr(mol, attribute, None)
        source = self.velocity_source.strip()
        if source.lower() not in {"zero", "none", "0", "maxwell", "boltzmann", "random"}:
            inputs["[md] velocity"] = os.path.expanduser(source)
        for key in ("system", "system2"):
            geometry = getattr(mol, "config", {}).get("input", {}).get(key, "")
            if isinstance(geometry, str) and geometry.strip() and "\n" not in geometry.strip():
                candidate = geometry.strip().split()[0]
                if os.path.isfile(candidate):
                    inputs["[input] %s" % key] = candidate
        outputs = {"trajectory_file": self.trajectory_file,
                   "energy_file": self.energy_file}
        for label, path in inputs.items():
            if not path or not isinstance(path, str):
                continue
            for name, output in outputs.items():
                if os.path.realpath(output) == os.path.realpath(path):
                    raise ValueError(
                        f"[md] {name}={output!r} would overwrite {label}; "
                        "give the output another name")

    def _remove_com_velocity(self, velocity):
        velocity = np.asarray(velocity, dtype=float).reshape((self.natom, 3))
        return velocity - np.sum(
            self.mass[:, None] * velocity, axis=0
        ) / np.sum(self.mass)

    def _initial_velocity(self):
        source = self.velocity_source.strip().lower()
        if source in {"zero", "none", "0"}:
            return np.zeros((self.natom, 3), dtype=float)
        if source in {"maxwell", "boltzmann", "random"}:
            sigma = np.sqrt(KB_HARTREE * self.init_temp / self.mass)
            velocity = self.rng.normal(size=(self.natom, 3)) * sigma[:, None]
            return self._remove_com_velocity(velocity)
        path = os.path.abspath(os.path.expanduser(self.velocity_source))
        if not os.path.isfile(path):
            raise ValueError(
                f"[md] velocity={self.velocity_source!r} is not zero/maxwell "
                "or a readable file"
            )
        velocity = np.loadtxt(path, dtype=float).reshape((self.natom, 3))
        if not np.all(np.isfinite(velocity)):
            raise ValueError("[md] velocity file contains NaN or infinity")
        return self._remove_com_velocity(velocity)

    def _electronic(self, *, continuation):
        mol = self.mol
        guess = mol.config["guess"]
        properties = mol.config["properties"]
        original_guess = guess.get("type", _MISSING)
        original_grad = properties.get("grad", _MISSING)
        try:
            reuse_orbitals = (
                continuation
                and _as_bool(mol.config["md"].get("mo_reuse", True))
            )
            if reuse_orbitals:
                guess["type"] = "previous"
            properties["grad"] = [0]
            # A resident production-basis orbital vector cannot be consumed by
            # an initial-basis SCF.  Continuation therefore starts directly in
            # the production basis; fresh evaluations retain the configured
            # initial-SCF convergence aid.
            SinglePoint(mol).energy(do_init_scf=not reuse_orbitals)
            Gradient(mol).gradient()
            LastStep(mol).compute(mol, grad_list=[0])
            energy = float(np.asarray(mol.energies, dtype=float).reshape(-1)[0])
            gradient = np.asarray(mol.grads[0], dtype=float).reshape((self.natom, 3))
            if not np.isfinite(energy) or not np.all(np.isfinite(gradient)):
                raise RuntimeError(
                    "ground-state MD received a non-finite energy or gradient"
                )
            return energy, gradient
        finally:
            if original_guess is _MISSING:
                guess.pop("type", None)
            else:
                guess["type"] = original_guess
            if original_grad is _MISSING:
                properties.pop("grad", None)
            else:
                properties["grad"] = original_grad

    def _apply_langevin(self):
        self._thermostat_exchange = 0.0
        if self.thermostat == "off":
            return
        kinetic_before = self._kinetic_energy()
        decay = np.exp(-self.friction * self.dt_fs / 1000.0)
        sigma = np.sqrt(
            (1.0 - decay * decay) * KB_HARTREE * self.temperature / self.mass
        )
        self.velocity = (
            decay * self.velocity
            + self.rng.normal(size=(self.natom, 3)) * sigma[:, None]
        )
        self.velocity = self._remove_com_velocity(self.velocity)
        self._thermostat_exchange = self._kinetic_energy() - kinetic_before
        self._thermostat_exchange_cumulative += self._thermostat_exchange

    def _kinetic_energy(self):
        return float(0.5 * np.sum(self.mass[:, None] * self.velocity ** 2))

    def _instantaneous_temperature(self):
        dof = max(1, 3 * self.natom - 3)
        return 2.0 * self._kinetic_energy() / (dof * KB_HARTREE)

    def _write_xyz(self, stream, coordinates, step, total_energy):
        atoms = np.asarray(self.mol.get_atoms(), dtype=int).reshape(-1)
        stream.write(f"{self.natom}\n")
        stream.write(
            f"step={step} time_fs={step * self.dt_fs:.8f} "
            f"energy_hartree={total_energy:.12f}\n"
        )
        xyz = np.asarray(coordinates) * BOHR_TO_ANGSTROM
        for atomic_number, row in zip(atoms, xyz):
            symbol = ELEMENTS_NAME[int(atomic_number)].strip()
            stream.write(
                f"{symbol:2s} {row[0]:18.10f} {row[1]:18.10f} {row[2]:18.10f}\n"
            )

    def _record(self, xyz_stream, csv_writer, step, coordinates, potential):
        kinetic = self._kinetic_energy()
        total = potential + kinetic
        temperature = self._instantaneous_temperature()
        if step % getattr(self, "trajectory_interval", 1) == 0:
            self._write_xyz(xyz_stream, coordinates, step, total)
        csv_writer.writerow([
            step, f"{step * self.dt_fs:.10f}", f"{potential:.16g}",
            f"{kinetic:.16g}", f"{total:.16g}", f"{temperature:.10f}",
            f"{self._thermostat_exchange:.16g}",
            f"{self._thermostat_exchange_cumulative:.16g}",
        ])
        dump_log(
            self.mol,
            title=(
                f"MD step {step:6d}  t={step * self.dt_fs:9.3f} fs  "
                f"E_tot={total:.10f}  E_pot={potential:.10f}  "
                f"E_kin={kinetic:.10f}  T={temperature:.3f} K  "
                f"dQ_therm={self._thermostat_exchange:+.3e}"
            ),
        )

    def run(self):
        mol = self.mol
        dump_log(mol, title="PyOQP: Ground-State Born-Oppenheimer Molecular Dynamics")
        coordinates = np.asarray(mol.get_system(), dtype=float).reshape((self.natom, 3))
        potential, gradient = self._electronic(continuation=False)
        acceleration = -gradient / self.mass[:, None]

        os.makedirs(os.path.dirname(os.path.abspath(self.trajectory_file)), exist_ok=True)
        os.makedirs(os.path.dirname(os.path.abspath(self.energy_file)), exist_ok=True)
        with open(self.trajectory_file, "w", encoding="utf-8") as xyz_stream, open(
                self.energy_file, "w", newline="", encoding="utf-8") as energy_stream:
            writer = csv.writer(energy_stream)
            writer.writerow([
                "step", "time_fs", "potential_hartree", "kinetic_hartree",
                "total_hartree", "temperature_kelvin",
                "thermostat_exchange_hartree",
                "thermostat_exchange_cumulative_hartree",
            ])
            self._record(xyz_stream, writer, 0, coordinates, potential)

            for step in range(1, self.nstep + 1):
                coordinates = (
                    coordinates + self.velocity * self.dt
                    + 0.5 * acceleration * self.dt * self.dt
                )
                mol.update_system(coordinates.reshape(-1))
                potential, gradient = self._electronic(continuation=True)
                acceleration_new = -gradient / self.mass[:, None]
                self.velocity += 0.5 * (acceleration + acceleration_new) * self.dt
                self._apply_langevin()
                acceleration = acceleration_new
                self._record(xyz_stream, writer, step, coordinates, potential)

        mol.md_velocity = self.velocity.copy()
        dump_log(
            mol,
            title=(
                "PyOQP: ground-state MD trajectory complete; "
                f"trajectory={self.trajectory_file}; energies={self.energy_file}"
            ),
        )
