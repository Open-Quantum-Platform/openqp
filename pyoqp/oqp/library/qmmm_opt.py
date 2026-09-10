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
import math
import os
import numpy as np
from openmm import app, unit
from oqp.library.oqp_engine import OQPEngine
from oqp.library.liboqp import _OQPRunner

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

        # the parsed schema supplies "" for an omitted key, so blank == missing
        pdb_file = str(qmmm_cfg.get("pdb_file") or "").strip()
        if not pdb_file:
            raise ValueError("'qmmm.pdb_file' is required for a QM/MM optimisation.")
        self.pdb = app.PDBFile(self._resolve_aux_file(pdb_file))
        ff_files = self._forcefield_paths(qmmm_cfg.get("forcefield_files", ""))
        if not ff_files:
            raise ValueError("'qmmm.forcefield_files' is required for a QM/MM optimisation.")
        self.forcefield = app.ForceField(*ff_files)
        qm_raw = qmmm_cfg.get("qm_atoms")
        if qm_raw is None or (isinstance(qm_raw, str) and not qm_raw.strip()):
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
        self._validate_qm_molecule_layout()

        # ---- state: [optimize] istate, as for the all-QM optimizer (0 = the
        # ground state of an HF/DFT run, n = the n-th TDHF/MRSF root); the
        # driver reads the root it differentiates from [properties] grad.
        self.istate = int(opt.get("istate", 0))
        self._validate_istate(self.istate, mol.config["input"].get("method", "hf"))
        self._reject_swapmo(mol.config.get("guess", {}).get("swapmo", ""))
        self._reject_continue_geom(mol.config.get("guess", {}).get("continue_geom", False))
        mol.config.setdefault("properties", {})["grad"] = [self.istate]

        # ---- movable set: QM atoms + whole MM residues within qmmm_radius ----
        self.radius = float(opt.get("qmmm_radius", 0.0))
        if not math.isfinite(self.radius) or self.radius < 0.0:
            raise ValueError("'optimize.qmmm_radius' must be a finite distance >= 0 angstrom, "
                             f"got {opt.get('qmmm_radius')!r}.")
        self.movable = self._movable_atoms(self.radius)
        # [qmmm] rigidwater / constraints, held on the movable MM atoms (may
        # add constrained partners to the movable set)
        self.frozen_pairs = self._constraint_pairs(qmmm_cfg)
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
        self._last_gradient = None      # full-system gradient of the last evaluation

    # ------------------------------------------------------------------ #
    @staticmethod
    def _resolve_coordsys(value):
        """[oqp] coordsys for a QM/MM optimisation.  'auto' means Cartesian.
        Cartesian and TRIC keep the collective translations and rotations of
        the movable atoms, which are real degrees of freedom against the fixed
        MM atoms and the periodic cell; DLC/RIC drop them (the engine's
        isolated-molecule assumption), so they are rejected."""
        cs = str(value or "auto").strip().lower()
        if cs == "auto":
            return "cartesian"
        if cs in ("cart", "cartesian", "tric"):
            return cs
        raise ValueError(f"[oqp] coordsys={value!r} is not available for a QM/MM optimisation: DLC/RIC "
                         "remove the collective translations and rotations of the movable atoms, which "
                         "move against the fixed MM atoms; use auto (Cartesian), cartesian or tric.")

    @staticmethod
    def _reject_continue_geom(value):
        """[guess] continue_geom restores a saved QM-fragment geometry, but this
        driver starts from [qmmm] pdb_file, so the restart would silently rerun
        from the original structure."""
        on = value if isinstance(value, bool) else str(value or "").strip().lower() in ("1", "true", "yes", "on", "t")
        if on:
            raise ValueError("[guess] continue_geom is not used by the QM/MM optimisation, which starts from "
                             "[qmmm] pdb_file; to continue, use the optimised full-system PDB "
                             "([optimize] qmmm_output) as the new pdb_file.")

    @staticmethod
    def _reject_swapmo(value):
        """[guess] swapmo is applied by SinglePoint.reference(); the embedded
        SCF of this driver builds its guesses directly, so a requested
        non-Aufbau occupation would be silently lost."""
        given = (len(value) > 0) if isinstance(value, (list, tuple)) else bool(str(value or "").strip())
        if given:
            raise ValueError("[guess] swapmo is not applied by the QM/MM optimisation's embedded SCF; "
                             "remove it, or optimise without qmmm_flag.")

    @staticmethod
    def _validate_istate(istate, method):
        """[optimize] istate: >= 0, and >= 1 for TDHF/MRSF (a response root,
        1 = the lowest); a negative root would index the state arrays from the end."""
        tdhf = str(method or "").strip().lower() == "tdhf"
        if istate < 0 or (tdhf and istate < 1):
            raise ValueError(f"[optimize] istate={istate} is not a valid state for a QM/MM optimisation: "
                             + ("use >= 1 (1 = the lowest MRSF/TDHF root)." if tdhf else "use >= 0."))

    def _validate_qm_molecule_layout(self):
        """The QM Molecule was built from [input] system before this driver
        existed.  It must be the [qmmm] selection in topology order followed by
        one hydrogen per cut bond (the NAMD driver's check), or every
        evaluation would put this driver's coordinates on the wrong atoms, or
        past the end of the native coordinate buffer."""
        nlink = len(getattr(self.driver, "link_atoms", None) or [])
        z_mol = np.asarray(self.mol.get_atoms2("charge"), dtype=float).reshape(-1)
        z_top = {a.index: (0 if a.element is None else a.element.atomic_number)
                 for a in self.pdb.topology.atoms()}
        if any(int(i) not in z_top for i in self.qm_atoms):
            raise ValueError(f"QM/MM optimisation: [qmmm] qm_atoms reaches beyond the "
                             f"{len(z_top)} atoms of [qmmm] pdb_file.")
        z_expected = [z_top[int(i)] for i in self.qm_atoms] + [1] * nlink
        if z_mol.size != len(z_expected) or any(abs(z_mol[k] - z_expected[k]) > 0.5
                                               for k in range(z_mol.size)):
            raise ValueError(
                "QM/MM optimisation: the QM molecule built from [input] system does not match "
                f"[qmmm] pdb_file / qm_atoms: molecule Z={z_mol.astype(int).tolist()}, expected "
                f"{z_expected} ({len(self.qm_atoms)} QM atoms in topology order + {nlink} link "
                "hydrogen(s)). '[input] system = file.pdb ...' indices are 1-based while "
                "[qmmm] qm_atoms are 0-based.")

    def _unwrap_constrained(self, X_nm, pairs):
        """Coordinates (nm) with every group of constrained atoms made whole
        under the periodic cell: from the lowest index of each connected
        group, each partner is placed at the minimum-image vector from the
        atom it was reached from, so the engine's frozen distances start from
        bond lengths even for atoms stored on opposite faces of the cell.
        Without a cell the coordinates are returned unchanged."""
        box = self.driver._box_lengths_bohr()
        if box is None or not pairs:
            return X_nm
        L = np.asarray(box, dtype=float) * BOHR_TO_NM
        adj = {}
        for i, j in pairs:
            adj.setdefault(int(i), []).append(int(j))
            adj.setdefault(int(j), []).append(int(i))
        X = np.array(X_nm, dtype=float)
        seen = set()
        for root in sorted(adj):
            if root in seen:
                continue
            seen.add(root)
            stack = [root]
            while stack:
                a = stack.pop()
                for b in adj[a]:
                    if b in seen:
                        continue
                    d = X[b] - X[a]
                    X[b] = X[a] + d - L * np.round(d / L)
                    seen.add(b)
                    stack.append(b)
        return X

    def _constraint_pairs(self, qmmm_cfg):
        """[qmmm] rigidwater / constraints for the movable MM atoms, as
        topology index pairs the engine holds at their starting distance
        (the native optimizer's frozen distances: tangent projection plus a
        SHAKE-like correction).  QM atoms are never constrained, as in the
        MD and NAMD drivers.  A constrained partner of a movable atom joins
        the movable set, so no constraint ties a moving atom to a fixed one."""
        rigid = str(qmmm_cfg.get("rigidwater", False)).strip().lower() in ("1", "true", "yes", "on")
        name = str(qmmm_cfg.get("constraints", "None") or "None").strip().lower()
        kinds = {"none": None, "": None, "hbonds": app.HBonds, "allbonds": app.AllBonds, "hangles": app.HAngles}
        if name not in kinds:
            raise ValueError("[qmmm] constraints must be None, HBonds, AllBonds or HAngles, "
                             f"got {qmmm_cfg.get('constraints')!r}.")
        if not rigid and kinds[name] is None:
            return []
        ref = self.forcefield.createSystem(self.pdb.topology, nonbondedMethod=app.NoCutoff,
                                          constraints=kinds[name], rigidWater=rigid)
        qm = set(int(i) for i in self.qm_atoms)
        allpairs = []
        for k in range(ref.getNumConstraints()):
            p1, p2, _ = ref.getConstraintParameters(k)
            if int(p1) in qm or int(p2) in qm:
                continue
            allpairs.append((int(p1), int(p2)))
        movable = set(int(i) for i in self.movable)
        changed = True
        while changed:
            changed = False
            for p1, p2 in allpairs:
                if (p1 in movable) != (p2 in movable):
                    movable.update((p1, p2))
                    changed = True
        self.movable = np.array(sorted(movable), dtype=int)
        return [(p1, p2) for p1, p2 in allpairs if p1 in movable]

    def _forcefield_paths(self, raw):
        """[qmmm] forcefield_files -> list of paths.  The unsplit value is
        resolved against the deck first, so a single deck-relative file whose
        name contains spaces stays one file; only then is it split into a
        list (commas or whitespace), each entry resolved the same way."""
        text = str(raw or "").strip()
        if not text:
            return []
        whole = self._resolve_aux_file(text)
        if os.path.isfile(whole):
            return [whole]
        return [self._resolve_aux_file(f) for f in _parse_str_list(text)]

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
            # virtual sites (element None, e.g. the TIP4P M site) are not
            # independent coordinates: OpenMM places them from their parents
            idx = [a.index for a in res.atoms() if a.element is not None]
            if not idx:
                continue
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
    def _with_virtual_sites(self, positions_nm):
        """Positions (nm) with every virtual site placed from its parent atoms
        by OpenMM, since the optimiser moves only real atoms."""
        X = np.array(positions_nm, dtype=float)
        vs = getattr(self, "_virtual_sites", None)
        if vs is None:
            sys0 = self.driver.mm_systems.get("sys0")
            vs = ([i for i in range(sys0.getNumParticles()) if sys0.isVirtualSite(i)]
                  if sys0 is not None else [])
            self._virtual_sites = vs
        if vs:
            ctx = self.driver.mm_systems["sim0"].context
            ctx.setPositions(unit.Quantity(X, unit.nanometer))
            ctx.computeVirtualSites()
            P = np.asarray(ctx.getState(getPositions=True).getPositions(asNumpy=True)
                           .value_in_unit(unit.nanometer))
            X[vs] = P[vs]
        return X

    def _energy_force(self, positions_nm):
        """QM/MM energy (Hartree) and force on every atom (Hartree/bohr)."""
        pos = unit.Quantity(self._with_virtual_sites(positions_nm), unit.nanometer)
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
        f = np.asarray(f, dtype=float)
        self._last_gradient = -f
        return float(e), f

    def optimize(self):
        mol = self.mol
        X0 = np.array(self.pdb.positions.value_in_unit(unit.nanometer))
        mv = self.movable
        at = list(self.pdb.topology.atoms())
        symbols = [int(at[i].element.atomic_number) for i in mv]   # the engine takes atomic numbers
        pos = {int(a): k for k, a in enumerate(mv)}
        engine_pairs = [(pos[i] + 1, pos[j] + 1) for i, j in self.frozen_pairs]   # engine numbering is 1-based
        opt = mol.config.get("optimize", {})
        # the engine controls live in the [oqp] section (concise opt(trust=...)
        # is lowered there too), as for the all-QM native optimizer
        eng = mol.config.get("oqp", {})
        trust = float(eng.get("trust", 0.2))
        trust_max = float(eng.get("trust_max", 0.5))
        # the native optimizer's recovery stage: an unconverged search is
        # restarted from its lowest-energy geometry with a small trust radius,
        # a fresh model Hessian and a budget of max(maxit - used, recovery_maxit)
        auto_recovery = bool(eng.get("auto_recovery", True))
        recovery_maxit = int(eng.get("recovery_maxit", 30))
        recovery_trust = float(eng.get("recovery_trust", 0.02))
        # [oqp] coordsys: 'auto' means Cartesian here (the movable set can be
        # several disconnected fragments, and Cartesian coordinates are safe
        # for that); an explicit choice is passed to the engine as requested.
        self.coordsys = self._resolve_coordsys(eng.get("coordsys", "auto"))
        dump_log(mol, title=(f"PyOQP: QM/MM geometry optimisation: {len(self.qm_atoms)} QM atoms, "
                             f"{len(mv)} movable atoms (radius {self.radius:.1f} A, "
                             f"{len(self.frozen_pairs)} constrained distances), "
                             f"{len(X0) - len(mv)} fixed; state {self.istate}; native RFO/BFGS, "
                             f"{self.coordsys} coordinates, trust {trust:.2f} (max {trust_max:.2f}) bohr, "
                             f"maxit {self.maxit}"))
        it = [0]
        electronic_failure = [None]
        self._prev_eval = None          # the point the next step and energy change are measured from

        def energy_gradient(x_bohr):
            X = X0.copy()
            X[mv] = np.asarray(x_bohr, dtype=float).reshape(-1, 3) * BOHR_TO_NM
            # the driver and mol move to this geometry even if the solve fails
            self._attempted_x = np.array(x_bohr, dtype=float)
            self._last_attempt_ok = False
            try:
                e, f = self._energy_force(X)
            except Exception as error:
                # as the native optimizer: an SCF/response solve that does not
                # converge at a trial geometry rejects that geometry; with no
                # evaluated geometry to fall back on it is a real failure
                if not self.history or not _OQPRunner._is_electronic_nonconvergence(error):
                    raise
                electronic_failure[0] = error
                dump_log(mol, title=("PyOQP: QM/MM optimisation: electronic solver did not converge "
                                     "at the trial geometry; rejecting it.\n   %s" % error))
                raise StopIteration from error
            self.driver._reuse_orbitals = not self.init_scf   # later steps start from these orbitals
            self._last_attempt_ok = True
            g = -f[mv].reshape(-1)
            if engine_pairs:
                # as the native optimizer: the constrained gradient drives the
                # search and the convergence test
                g = engine._project_constraint_tangent(g, np.asarray(x_bohr, dtype=float))
            it[0] += 1
            rms = float(np.sqrt(np.mean(g * g))); mx = float(np.abs(g).max())
            prev = self._prev_eval
            if prev is not None:
                dx = np.asarray(x_bohr, dtype=float) - prev["x"]
                rms_step, max_step = float(np.sqrt(np.mean(dx * dx))), float(np.abs(dx).max())
                de = e - prev["e"]
            else:
                rms_step = max_step = de = float("inf")
            self.history.append({"x": np.array(x_bohr, dtype=float), "e": e, "rms": rms, "max": mx,
                                 "rms_step": rms_step, "max_step": max_step, "de": de, "id": it[0]})
            self._prev_eval = self.history[-1]
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

        if self.frozen_pairs:
            X0 = self._unwrap_constrained(X0, self.frozen_pairs)
        x0 = (X0[mv] / BOHR_TO_NM).reshape(-1)
        engine = OQPEngine(symbols, x0, mode="min", trust=trust, trust_max=trust_max,
                           maxiter=self.maxit, coordsys=self.coordsys, frozen_distances=engine_pairs)
        recovered = False
        try:
            engine.run(energy_gradient, on_converged=on_converged)
            if not converged_at(self.history[-1]) and auto_recovery:
                best = min(self.history, key=lambda h: h["e"])
                steps = max(self.maxit - it[0], recovery_maxit)
                r_trust = min(recovery_trust, trust_max)
                r_trust_max = min(trust_max, max(r_trust, 2.5 * r_trust))
                reason = ("electronic non-convergence at a trial geometry"
                          if electronic_failure[0] is not None else "the iteration limit")
                dump_log(mol, title=(f"PyOQP: QM/MM optimisation recovery selected after {reason}. "
                                     f"Restarting the lowest-energy geometry (E = {best['e']:.10f} "
                                     f"Hartree) with {self.coordsys} coordinates, trust={r_trust:.3f}, "
                                     f"and a fresh model Hessian for up to {steps} steps"))
                recovered = True
                self._prev_eval = best      # the restarted search measures its first step from here
                engine = OQPEngine(symbols, best["x"].copy(), mode="min", trust=r_trust,
                                   trust_max=r_trust_max, maxiter=steps, coordsys=self.coordsys,
                                   frozen_distances=engine_pairs)
                engine.run(energy_gradient, on_converged=on_converged)
        finally:
            self.driver._reuse_orbitals = False
        best = min(self.history, key=lambda h: h["e"])
        last = self.history[-1]
        converged = converged_at(last)
        final = last if converged else best
        X = X0.copy(); X[mv] = final["x"].reshape(-1, 3) * BOHR_TO_NM
        held_is_final = bool(getattr(self, "_last_attempt_ok", False)) and last.get("id") == final.get("id")
        if not held_is_final:
            # the driver and mol hold the last attempted evaluation, which is not
            # the one we report (a later trial, a rejected one, or another
            # evaluation at the same coordinates, as a recovery restart makes):
            # bring them back to the reported geometry, and report what this
            # evaluation returns, so energy, metrics and gradient belong together
            e_re, f_re = self._energy_force(X)
            g_re = -f_re[mv].reshape(-1)
            if engine_pairs:
                g_re = engine._project_constraint_tangent(g_re, final["x"])
            final = dict(final, e=float(e_re), rms=float(np.sqrt(np.mean(g_re * g_re))),
                         max=float(np.abs(g_re).max()))
        X = self._with_virtual_sites(X)
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
        # Publish the result the way Runner.results() reads it: index istate of
        # mol.energies holds the QM/MM total energy (Hartree) of the reported
        # geometry (the embedded QM state energy + nuclear-MM + MM terms, the
        # objective that was minimised); the QM-fragment records themselves
        # are left as the last evaluation set them.
        energies = [float("nan")] * self.istate + [float(final["e"])]
        mol.energies = energies
        # and the gradient that belongs to that energy, shaped like the atoms
        # and coordinates Runner.results() and the saved JSON publish (the QM
        # molecule): row k is the derivative of the QM/MM objective with
        # respect to QM atom k (link-atom forces already carried to their
        # hosts); the link hydrogens are not independent coordinates and get
        # zero rows.  The full-system gradient (every PDB atom) of the same
        # evaluation is kept as self.gradient_full.
        self.gradient_full = np.array(self._last_gradient, dtype=float)
        nmol = np.asarray(mol.get_atoms2("charge")).reshape(-1).size
        grad = np.zeros((nmol, 3))
        grad[:len(self.qm_atoms)] = self.gradient_full[self.qm_atoms]
        mol.grads = [np.full_like(grad, float("nan"))] * self.istate + [grad]
        mol.qmmm_optimization = {
            "converged": bool(converged), "energy_hartree": float(final["e"]),
            "evaluations": int(it[0]), "recovery": bool(recovered),
            "electronic_failure": electronic_failure[0] is not None,
            "rms_grad": float(final["rms"]), "max_grad": float(final["max"]),
            "constraints": len(self.frozen_pairs),
            "movable_atoms": [int(i) for i in mv], "output": self.output,
            "grad_fragment": grad.tolist(),
        }
        return final["e"], X


def run_qmmm_optimization(mol):
    opt = QMMM_Opt(mol)
    opt.optimize()
    if mol.config["guess"].get("save_mol"):
        mol.save_data()
    return opt
