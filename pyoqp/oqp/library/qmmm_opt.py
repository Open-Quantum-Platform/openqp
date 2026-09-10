"""Geometry optimisation on the QM/MM potential (``runtype = optimize`` with
``qmmm_flag = True``).

The all-QM optimiser moves the QM coordinates and asks the gas-phase engine for
energy and gradient; it knows nothing about MM atoms, the embedding field or
the link atoms.  This driver instead minimises the QM/MM energy returned by
``OpenQpQMMM.compute_force`` (embedded SCF with the ESPF field, the ESPF and
link-atom gradient terms, the classical MM forces) with respect to a chosen
set of MOVABLE atoms: the QM region plus, optionally, every MM residue with an
atom within ``[optimize] qmmm_radius`` angstrom of a QM atom.  Everything else
is held fixed, which is the usual practice for a solvated system and is what
keeps the problem well posed (the classical solvent has no minimum worth
finding).

State selection follows the dynamics drivers: ``[properties] grad`` names the
TDHF/MRSF root (1 = the lowest root); ground-state HF/DFT needs nothing.

Outputs: the optimised full-system PDB (``[optimize] qmmm_output``, default
``<project>_opt.pdb``), one log line per iteration, and the usual mol data
(energies, geometry) for the QM fragment.
"""
import os
import numpy as np
from openmm import app, unit
from oqp.library.oqp_engine import OQPEngine

from oqp.library.qmmm_driver import OpenQpQMMM
from oqp.library.qmmm_md import (
    _extract_qmmm_config, _parse_int_list, _parse_str_list, _resolve_cutoff)
from oqp.utils.file_utils import dump_log

HARTREE_TO_KJMOL = 2625.499639
BOHR_TO_NM = 0.052917721067
FORCE_KJMOLNM_TO_HABOHR = 1.0 / 49614.75


class QMMM_Opt:
    """Minimise the QM/MM energy over the movable atoms with the native
    trust-radius RFO/BFGS engine that drives the all-QM optimizer."""

    def __init__(self, mol):
        self.mol = mol
        qmmm_cfg, qm_cfg = _extract_qmmm_config(mol=mol)
        opt = mol.config.get("optimize", {})

        pdb_file = qmmm_cfg.get("pdb_file")
        if pdb_file is None:
            raise ValueError("'qmmm.pdb_file' is required for a QM/MM optimisation.")
        self.pdb = app.PDBFile(self._resolve_aux_file(pdb_file))
        ff_files = _parse_str_list(qmmm_cfg.get("forcefield_files", ""))
        if not ff_files:
            raise ValueError("'qmmm.forcefield_files' is required for a QM/MM optimisation.")
        self.forcefield = app.ForceField(*ff_files)
        qm_raw = qmmm_cfg.get("qm_atoms")
        if qm_raw is None:
            raise ValueError("'qmmm.qm_atoms' is required for a QM/MM optimisation.")
        self.qm_atoms = np.array(sorted(_parse_int_list(qm_raw)), dtype=int)

        cutoff = _resolve_cutoff(qmmm_cfg.get("cutoff", "PME"))
        _et = qmmm_cfg.get("ewald_tol", None)
        _w = qmmm_cfg.get("mm_charge_width", None)
        _flag = lambda k: str(qmmm_cfg.get(k, "false")).strip().lower() in ("1", "true", "yes", "on")
        self.driver = OpenQpQMMM(
            positions=self.pdb.positions, topology=self.pdb.topology,
            forcefield=self.forcefield, qm_atoms=self.qm_atoms, mol=mol,
            Cutoff=cutoff, Embedding=str(qmmm_cfg.get("embedding", "electrostatic")),
            frontier_scheme=str(qmmm_cfg.get("frontier_scheme", "none")),
            ewald_tol=None if _et in (None, "", "none", "None") else float(_et),
            lj_switch=_flag("lj_switch"), h_lj=_flag("h_lj"),
            mm_charge_width=None if _w in (None, "", "none", "None", 0, 0.0, "0") else float(_w),
        )

        # ---- state: [optimize] istate, as for the all-QM optimizer (0 = the
        # ground state of an HF/DFT run, n = the n-th TDHF/MRSF root); the
        # driver reads the root it differentiates from [properties] grad.
        self.istate = int(opt.get("istate", 0))
        mol.config.setdefault("properties", {})["grad"] = [self.istate]

        # ---- movable set: QM atoms + whole MM residues within qmmm_radius ----
        self.radius = float(opt.get("qmmm_radius", 0.0))
        self.movable = self._movable_atoms(self.radius)
        self.maxit = max(1, int(opt.get("maxit", 30)))    # the checker rejects < 1; never skip the first evaluation
        # [optimize] init_scf=true asks for a fresh initial SCF at every geometry
        # (the all-QM optimizer's policy); otherwise the converged orbitals of
        # the previous step are the guess.
        self.init_scf = bool(opt.get("init_scf", False))
        # the same five-part test as the native all-QM optimizer
        # (_native_metrics_converged): energy change, rms/max step, rms/max gradient
        self.rmsd_grad = float(opt.get("rmsd_grad", 1e-4))     # Hartree/bohr
        self.max_grad = float(opt.get("max_grad", 3e-4))
        self.rmsd_step = float(opt.get("rmsd_step", 1e-3))     # bohr
        self.max_step = float(opt.get("max_step", 2e-3))
        self.energy_shift = float(opt.get("energy_shift", 1e-6))   # Hartree
        project = os.path.splitext(os.path.basename(str(mol.config["input"].get("system", "qmmm")).split()[0]))[0]
        self.output = str(opt.get("qmmm_output", "") or f"{getattr(mol, 'project_name', project)}_opt.pdb")
        self.history = []

    # ------------------------------------------------------------------ #
    def _resolve_aux_file(self, name):
        """A relative [qmmm] path not found in the working directory is looked
        up next to the input deck, the rule the NAMD driver applies, so a deck
        can be run from any directory (e.g. by ``openqp --run_tests``)."""
        value = str(name or "")
        input_file = getattr(self.mol, "input_file", None)
        if value and input_file and not os.path.isabs(value) and not os.path.exists(value):
            candidate = os.path.join(os.path.dirname(os.path.abspath(input_file)), value)
            if os.path.exists(candidate):
                return candidate
        return value

    def _movable_atoms(self, radius):
        qm = set(int(i) for i in self.qm_atoms)
        if radius <= 0.0:
            return np.array(sorted(qm), dtype=int)
        X = np.array(self.pdb.positions.value_in_unit(unit.angstrom))
        box = self.driver._box_lengths_bohr()
        box_ang = None if box is None else np.asarray(box) / 1.8897259886
        qx = X[sorted(qm)]
        movable = set(qm)

        def dist(idx):
            d = X[idx][:, None, :] - qx[None, :, :]
            if box_ang is not None:
                d -= box_ang * np.round(d / box_ang)
            return np.linalg.norm(d, axis=2).min(axis=1)     # per atom: nearest QM atom

        for res in self.pdb.topology.residues():
            idx = [a.index for a in res.atoms()]
            if qm.intersection(idx):
                # A covalent cut inside a residue: its MM atoms (the link host
                # and beyond) are selected one by one, otherwise the residue
                # that carries the QM atoms could never move at all.
                mm_here = [i for i in idx if i not in qm]
                if mm_here:
                    movable.update(int(i) for i, r in zip(mm_here, dist(mm_here)) if r <= radius)
            elif dist(idx).min() <= radius:
                movable.update(idx)          # whole residue: keeps waters and side chains intact
        return np.array(sorted(movable), dtype=int)

    # ------------------------------------------------------------------ #
    def _energy_force(self, positions_nm):
        """QM/MM energy (Hartree) and force on every atom (Hartree/bohr)."""
        pos = unit.Quantity(np.asarray(positions_nm, dtype=float), unit.nanometer)
        # compute_force takes the QM geometry from ``positions`` but the MM
        # energy and forces from the OpenMM contexts, which the caller owns
        # (the MD driver moves them with its integrator).  Move every context
        # first, or the classical part is evaluated at the starting geometry
        # while the QM atoms walk away from it.
        for key, sim in self.driver.mm_systems.items():
            if hasattr(sim, "context"):
                sim.context.setPositions(pos)
        e_q, f_q = self.driver.compute_force(pos, self.pdb.topology, self.driver.mm_systems, self.qm_atoms)
        e = e_q.value_in_unit(unit.kilojoule_per_mole) / HARTREE_TO_KJMOL
        f = (f_q.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
             if hasattr(f_q, "value_in_unit") else np.asarray(f_q)) * FORCE_KJMOLNM_TO_HABOHR
        return float(e), np.asarray(f, dtype=float)

    def optimize(self):
        mol = self.mol
        X0 = np.array(self.pdb.positions.value_in_unit(unit.nanometer))
        mv = self.movable
        at = list(self.pdb.topology.atoms())
        symbols = [int(at[i].element.atomic_number) for i in mv]   # the engine takes atomic numbers
        opt = mol.config.get("optimize", {})
        # the engine controls live in the [oqp] section (concise opt(trust=...)
        # is lowered there too), as for the all-QM native optimizer
        eng = mol.config.get("oqp", {})
        trust = float(eng.get("trust", 0.2))
        trust_max = float(eng.get("trust_max", 0.5))
        dump_log(mol, title=(f"PyOQP: QM/MM geometry optimisation: {len(self.qm_atoms)} QM atoms, "
                             f"{len(mv)} movable atoms (radius {self.radius:.1f} A), "
                             f"{len(X0) - len(mv)} fixed; state {self.istate}; native RFO/BFGS, "
                             f"Cartesian, trust {trust:.2f} (max {trust_max:.2f}) bohr, maxit {self.maxit}"))
        it = [0]

        def energy_gradient(x_bohr):
            X = X0.copy()
            X[mv] = np.asarray(x_bohr, dtype=float).reshape(-1, 3) * BOHR_TO_NM
            e, f = self._energy_force(X)
            self.driver._reuse_orbitals = not self.init_scf   # later steps start from these orbitals
            g = -f[mv].reshape(-1)
            it[0] += 1
            rms = float(np.sqrt(np.mean(g * g))); mx = float(np.abs(g).max())
            if self.history:
                dx = np.asarray(x_bohr, dtype=float) - self.history[-1]["x"]
                rms_step, max_step = float(np.sqrt(np.mean(dx * dx))), float(np.abs(dx).max())
                de = e - self.history[-1]["e"]
            else:
                rms_step = max_step = de = float("inf")
            self.history.append({"x": np.array(x_bohr, dtype=float), "e": e, "rms": rms, "max": mx,
                                 "rms_step": rms_step, "max_step": max_step, "de": de})
            dump_log(mol, title=(f"PyOQP: QM/MM optimisation step {it[0]}: E = {e:.10f} Hartree, "
                                 f"dE {de:+.2e}, rms/max grad {rms:.2e}/{mx:.2e} Hartree/bohr, "
                                 f"rms/max step {rms_step:.2e}/{max_step:.2e} bohr"), section="")
            return e, g

        def converged_at(h):
            return (h["max"] <= self.max_grad and h["rms"] <= self.rmsd_grad
                    and h["max_step"] <= self.max_step and h["rms_step"] <= self.rmsd_step
                    and abs(h["de"]) <= self.energy_shift)

        def on_converged():
            if converged_at(self.history[-1]):
                raise StopIteration

        x0 = (X0[mv] / BOHR_TO_NM).reshape(-1)
        engine = OQPEngine(symbols, x0, mode="min", trust=trust, trust_max=trust_max,
                           maxiter=self.maxit, coordsys="cartesian")
        try:
            engine.run(energy_gradient, on_converged=on_converged)
        finally:
            self.driver._reuse_orbitals = False
        best = min(self.history, key=lambda h: h["e"])
        last = self.history[-1]
        converged = converged_at(last)
        final = last if converged else best
        X = X0.copy(); X[mv] = final["x"].reshape(-1, 3) * BOHR_TO_NM
        if final is not last:                          # leave mol/driver on the geometry we report
            self._energy_force(X)
        with open(self.output, "w") as fh:
            app.PDBFile.writeFile(self.pdb.topology, unit.Quantity(X, unit.nanometer), fh, keepIds=True)
        dump_log(mol, title=(f"PyOQP: QM/MM optimisation {'converged' if converged else 'NOT converged'} "
                             f"after {it[0]} evaluations: E = {final['e']:.10f} Hartree, "
                             f"rms grad {final['rms']:.2e}, max grad {final['max']:.2e} Hartree/bohr; "
                             f"full-system geometry written to {self.output}"))
        if not converged:
            dump_log(mol, title=(f"PyOQP: QM/MM optimisation reached maxit={self.maxit} without meeting "
                                 f"max_grad={self.max_grad:.0e} rmsd_grad={self.rmsd_grad:.0e} "
                                 f"max_step={self.max_step:.0e} rmsd_step={self.rmsd_step:.0e} "
                                 f"energy_shift={self.energy_shift:.0e}; the lowest-energy geometry "
                                 f"visited is written"), section="")
        self.converged = converged
        self.energy = final["e"]
        self.positions_nm = X
        return final["e"], X


def run_qmmm_optimization(mol):
    opt = QMMM_Opt(mol)
    opt.optimize()
    if mol.config["guess"].get("save_mol"):
        mol.save_data()
    return opt
