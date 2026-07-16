"""OQP geometry optimizer backend for OpenQP.

A NumPy/SciPy-only alternative to the external ``geomeTRIC`` driver.  The
electronic-structure work, convergence test and logging are reused verbatim from
:class:`oqp.library.libscipy.StateSpecificOpt`; this module only supplies the
step-determination loop through :class:`oqp.library.oqp_engine.OQPEngine`
(redundant internal coordinates + restricted-step RFO / P-RFO with model-Hessian
BFGS/Bofill updates).

Selected with ``[optimize] lib=oqp``.  Supported runtypes:

* ``optimize`` -> :class:`OQPOpt`   (state-specific minimum)
* ``ts``       -> :class:`OQPTSOpt` (transition state, eigenvector following)

Options are read from the ``[oqp]`` input section (see
``oqp.molecule.oqpdata``).
"""

from __future__ import annotations

import numpy as np

import os

from oqp.library.libscipy import Optimizer, StateSpecificOpt, MECIOpt, MECPOpt
from oqp.library.baeka import BaekAState
from oqp.library.oqp_engine import OQPEngine, parse_frozen_distance_spec
from oqp.library.oqp_neb import NEB
from oqp.library.neb_utils import _read_xyz, kabsch_align, write_neb_xyz
from oqp.periodic_table import ELEMENTS_NAME, SYMBOL_MAP
from oqp.utils.file_utils import dump_log, dump_data

ANGSTROM_TO_BOHR = 1.0 / 0.52917721092


class _OQPRunner:
    """Mixin that drives a PyOQP ``one_step`` objective with OQPEngine."""

    mode = "min"

    def _oqp_config(self):
        cfg = self.mol.config.get("oqp", {})
        return {
            "coordsys": cfg.get("coordsys", "auto"),
            "trust": float(cfg.get("trust", 0.2)),
            "trust_max": float(cfg.get("trust_max", 0.5)),
            "frozen_distances": parse_frozen_distance_spec(cfg.get("freeze", "")),
            "follow_mode": int(cfg.get("follow", 0)),
        }

    def _initial_hessian(self):
        """Return an optional Cartesian Hessian for the native engine."""
        self._initial_hessian_gradient = None
        return None

    def optimize(self):
        opts = self._oqp_config()
        atoms = np.asarray(self.mol.get_atoms(), dtype=int).reshape(-1)
        x0 = np.asarray(self.pre_coord, dtype=float).reshape(-1)

        initial_hessian = self._initial_hessian()
        engine = OQPEngine(
            atoms, x0,
            mode=self.mode,
            trust=opts["trust"],
            trust_max=opts["trust_max"],
            follow_mode=opts["follow_mode"],
            coordsys=opts["coordsys"],
            # OQPEngine's own loop is a backstop; the OQP convergence test
            # (which raises StopIteration at maxit) governs termination.
            maxiter=self.maxit,
            initial_hessian=initial_hessian,
            initial_gradient=getattr(self, "_initial_hessian_gradient", None),
            masses=np.asarray(self.mol.get_mass(), dtype=float).reshape(-1),
            project_global_rigid_modes=not bool(
                self.mol.config.get("input", {}).get("qmmm_flag", False)
            ),
            frozen_distances=opts["frozen_distances"],
        )
        dump_log(
            self.mol,
            title="PyOQP: OQP optimizer [mode=%s, coordsys=%s, trust=%.3f%s]"
            % (
                self.mode,
                engine.coordsys,
                opts["trust"],
                ", frozen distances=%s" % len(engine.frozen_distances)
                if engine.frozen_distances else "",
            ),
        )

        def energy_gradient(coords):
            # one_step updates metrics/pre_coord/pre_energy and returns
            # (energy [Hartree], gradient [Hartree/Bohr]) for the active state.
            energy, gradient = self.one_step(
                np.asarray(coords, dtype=float).reshape(-1)
            )
            if engine.frozen_distances:
                gradient = engine._project_constraint_tangent(gradient, coords)
                self.metrics["rmsd_grad"] = float(
                    np.sqrt(np.mean(np.asarray(gradient) ** 2))
                )
                self.metrics["max_grad"] = float(
                    np.max(np.abs(np.asarray(gradient)))
                )
            return energy, gradient

        try:
            engine.run(energy_gradient, on_converged=self.check_convergence)
        except StopIteration:
            pass


class OQPOpt(_OQPRunner, StateSpecificOpt):
    """State-specific minimum via redundant internals + restricted-step RFO."""

    mode = "min"

    def __init__(self, mol):
        super().__init__(mol)


class OQPTSOpt(_OQPRunner, StateSpecificOpt):
    """Transition-state search via partitioned RFO (eigenvector following)."""

    mode = "ts"

    def __init__(self, mol):
        super().__init__(mol)

    def _initial_hessian(self):
        policy = str(
            self.mol.config.get("oqp", {}).get("init_hessian", "model")
        ).strip().lower()
        if policy == "model":
            return None
        if policy not in {"numerical", "analytical"}:
            raise ValueError(
                "oqp.init_hessian must be model, numerical, or analytical"
            )

        from oqp.library.single_point import Hessian

        dump_log(
            self.mol,
            title="PyOQP: Native TS initial Hessian [%s]" % policy,
        )
        self.mol.update_system(np.asarray(self.pre_coord, dtype=float).reshape(-1))
        energies = self.sp.energy(do_init_scf=True)
        self.grad.grads = [self.istate]
        grads = self.grad.gradient()
        energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
        self.mol.energies = energies
        self.mol.grads = grads
        self._initial_hessian_gradient = np.asarray(
            grads[self.istate], dtype=float
        ).reshape(-1).copy()

        hess_cfg = self.mol.config["hess"]
        previous_type = hess_cfg.get("type", "numerical")
        previous_state = hess_cfg.get("state", 0)
        previous_read = hess_cfg.get("read", False)
        previous_restart = hess_cfg.get("restart", False)
        try:
            hess_cfg["type"] = policy
            hess_cfg["state"] = self.istate
            # ``hessian=numerical|analytical`` requests a fresh Hessian for the
            # current TS.  A stale legacy hess.read flag must not substitute a
            # cached matrix, and restart files from old displaced geometries
            # must not be assembled into the new numerical Hessian.
            hess_cfg["read"] = False
            hess_cfg["restart"] = False
            hessian = Hessian(self.mol).hessian(analysis=False)
            if hessian is None:
                raise RuntimeError("Native TS initial Hessian calculation failed")
            hessian = np.asarray(hessian, dtype=float)
        finally:
            hess_cfg["type"] = previous_type
            hess_cfg["state"] = previous_state
            hess_cfg["read"] = previous_read
            hess_cfg["restart"] = previous_restart

        expected = (3 * self.natom, 3 * self.natom)
        if hessian.shape != expected:
            raise ValueError(
                "Native TS Hessian has shape %s; expected %s"
                % (hessian.shape, expected)
            )
        return hessian


class OQPMECIOpt(_OQPRunner, MECIOpt):
    """MECI search: minimize the penalty/UBP objective with the oqp engine.

    Reuses ``MECIOpt.one_step`` (which returns the penalty objective and its
    gradient) and ``MECIOpt.check_convergence`` (which adds the energy-gap
    criterion) verbatim -- only the step determination is oqp.
    """

    mode = "min"

    def __init__(self, mol):
        MECIOpt.__init__(self, mol)


class OQPMECPOpt(_OQPRunner, MECPOpt):
    """MECP search: minimize the gap-penalty objective with the oqp engine."""

    mode = "min"

    def __init__(self, mol):
        MECPOpt.__init__(self, mol)


class OQPBaekAOpt(_OQPRunner, Optimizer):
    r"""Generalized Baek adaptive-penalty CI search for ``N >= 2`` states.

    The selected roots are ordered and only adjacent energy gaps are penalized,
    removing the redundant outer gap of an all-pairs formulation.  ``sigma`` is
    increased additively by ``pen_delta`` on ordinary iterations and by the next
    entry of ``pen_jump`` when the projected objective is locally stationary but
    the outer-state span still exceeds ``energy_gap``.

    Reference: Y. S. Baek, S. Lee, M. Filatov, and C. H. Choi,
    J. Phys. Chem. A 125, 1994--2006 (2021),
    https://doi.org/10.1021/acs.jpca.0c11294.
    """

    mode = "min"

    def __init__(self, mol, states=None):
        Optimizer.__init__(self, mol)
        cfg = mol.config["optimize"]

        if states is None:
            states = cfg.get("states", [])
        if isinstance(states, str):
            states = [item for item in states.replace(",", " ").split() if item]
        states = [int(state) for state in states]
        if not states:
            states = [self.istate, self.jstate]
        if len(states) < 2:
            raise ValueError("BaekA requires at least two states")
        if any(state < 1 for state in states):
            raise ValueError("BaekA states must be positive response roots")
        if any(b <= a for a, b in zip(states, states[1:])):
            raise ValueError("BaekA states must be strictly increasing")
        if any(b != a + 1 for a, b in zip(states, states[1:])):
            raise ValueError("BaekA requires consecutive states")
        self.states = states

        alpha = float(cfg.get("pen_alpha", 0.0))
        if alpha == 0.0:
            # The global legacy schema must retain pen_alpha=0 for the older
            # two-state MECI path.  For BaekA, zero is the method-default
            # sentinel and maps to the established energy-valued smoothing
            # parameter; negative values remain invalid.
            alpha = 0.02
        if alpha < 0.0:
            raise ValueError("BaekA pen_alpha must be zero or a positive energy")

        self.alpha = alpha
        self.delta_beta = float(cfg.get("pen_delta", 0.025))
        self.jump_schedule = cfg.get(
            "pen_jump", "10,10,25,25,100,100,1000,1000,3000"
        )
        self.weights = float(cfg.get("gap_weight", 1.0))
        if not np.isfinite(self.weights) or self.weights != 1.0:
            raise ValueError(
                "BaekA gap_weight is fixed at 1.0; use pen_sigma for the "
                "published penalty multiplier"
            )
        self.baeka = BaekAState(
            sigma=float(cfg.get("pen_sigma", 1.0)),
            alpha=self.alpha,
            delta_beta=self.delta_beta,
            beta_schedule=self.jump_schedule,
            gap_tol=self.energy_gap,
            tol_f=self.energy_shift,
            tol_g=self.rmsd_grad,
            gap_weight=self.weights,
        )
        self._baeka_evaluation = None
        self.metrics.update({
            "meci_search": "baeka",
            "states": self.states.copy(),
            "state_count": len(self.states),
            "alpha": self.alpha,
            "delta_beta": self.delta_beta,
            "jump_schedule": self.baeka.beta_schedule,
            "tol_f": self.energy_shift,
            "tol_g": self.rmsd_grad,
        })

    def optimize(self):
        """Run BaekA with same-sigma-rebased native RFO/BFGS steps."""

        opts = self._oqp_config()
        atoms = np.asarray(self.mol.get_atoms(), dtype=int).reshape(-1)
        x0 = np.asarray(self.pre_coord, dtype=float).reshape(-1)
        engine = OQPEngine(
            atoms, x0,
            mode=self.mode,
            trust=opts["trust"],
            trust_max=opts["trust_max"],
            follow_mode=opts["follow_mode"],
            coordsys=opts["coordsys"],
            maxiter=self.maxit,
            masses=np.asarray(self.mol.get_mass(), dtype=float).reshape(-1),
            project_global_rigid_modes=not bool(
                self.mol.config.get("input", {}).get("qmmm_flag", False)
            ),
        )
        dump_log(
            self.mol,
            title="PyOQP: OQP BaekA optimizer [states=%s, coordsys=%s, trust=%.3f]"
            % ("/".join(map(str, self.states)), engine.coordsys, opts["trust"]),
        )

        previous_objective_data = [None]

        def energy_gradient(coords):
            energy, gradient = self.one_step(
                np.asarray(coords, dtype=float).reshape(-1)
            )
            current_data = self._baeka_evaluation.data
            previous_data = previous_objective_data[0]
            if previous_data is not None:
                # Recombine the previous geometry at the exact sigma of the
                # objective returned for this geometry.  OQPEngine then forms
                # its BFGS secant and trust ratio from one common objective
                # without discarding accumulated curvature.
                previous_energy = (
                    previous_data.average_energy
                    + current_data.sigma * previous_data.penalty
                )
                previous_gradient = (
                    previous_data.average_gradient
                    + current_data.sigma * previous_data.penalty_gradient
                )
                engine.rebase_previous_objective(
                    previous_energy, previous_gradient
                )
            previous_objective_data[0] = current_data
            return energy, gradient

        engine.run(energy_gradient, on_converged=self.check_convergence)

    def one_step(self, coordinates):
        self.itr += 1
        dump_log(self.mol, title="PyOQP: Geometry Optimization Step %s" % self.itr)
        do_init_scf = True if self.itr == 1 else self.init_scf

        self.mol.update_system(coordinates)
        energies = self.sp.energy(do_init_scf=do_init_scf)
        self.grad.grads = self.states.copy()
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads
        energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
        self.mol.energies = energies
        self.mol.grads = grads

        selected_energies = np.asarray(
            [energies[state] for state in self.states], dtype=float
        )
        selected_gradients = np.asarray(
            [np.asarray(grads[state], dtype=float).reshape(-1)
             for state in self.states],
            dtype=float,
        )
        evaluation = self.baeka.evaluate(selected_energies, selected_gradients)
        self._baeka_evaluation = evaluation
        data = evaluation.data
        tested_data = evaluation.tested_data

        previous_coord = np.asarray(self.pre_coord, dtype=float).reshape(-1)
        current_coord = np.asarray(coordinates, dtype=float).reshape(-1)
        displacement = current_coord - previous_coord
        rmsd_step = float(np.sqrt(np.mean(displacement * displacement)))
        max_step = float(np.max(np.abs(displacement)))
        rmsd_grad = float(np.sqrt(np.mean(data.gradient * data.gradient)))
        max_grad = float(np.max(np.abs(data.gradient)))

        self.metrics.update({
            "itr": self.itr,
            "sigma": tested_data.sigma,
            "active_sigma": data.sigma,
            "next_sigma": evaluation.next_sigma,
            "effective_sigma": data.effective_sigma,
            "de": evaluation.delta_f,
            "gap": tested_data.span,
            "gaps": tested_data.gaps.copy(),
            "objective": data.objective,
            "tested_objective": tested_data.objective,
            "penalty": tested_data.penalty,
            "projector_rank": tested_data.projector_rank,
            "parallel_grad": tested_data.parallel_scaled_norm,
            "perpendicular_grad": tested_data.perpendicular_norm,
            "stationary": evaluation.local_stationary,
            "gap_converged": evaluation.gap_converged,
            "converged": evaluation.converged,
            "action": evaluation.action,
            "jump": evaluation.jump,
            "jump_index": evaluation.jump_index,
            "rmsd_step": rmsd_step,
            "max_step": max_step,
            "rmsd_grad": rmsd_grad,
            "max_grad": max_grad,
        })
        dump_data(
            self.mol,
            (
                self.itr, self.atoms, current_coord, data.objective,
                evaluation.delta_f, tested_data.span, tested_data.gaps,
                tested_data.sigma, data.sigma, rmsd_step, max_step,
                tested_data.parallel_scaled_norm,
                tested_data.perpendicular_norm, evaluation.action,
            ),
            title="BAEKA",
            fpath=self.mol.log_path,
        )
        self.pre_energy = data.objective
        self.pre_coord = current_coord.copy()
        return data.objective, data.gradient

    def check_convergence(self):
        dump_log(
            self.mol,
            title="Geometry Optimization Convergence %s" % self.itr,
            section="baeka",
            info=self.metrics,
        )
        if self._baeka_evaluation.converged:
            dump_log(self.mol, title="PyOQP: BaekA Geometry Optimization Has Converged")
            raise StopIteration
        if self.itr >= self.maxit:
            dump_log(
                self.mol,
                title="PyOQP: BaekA Geometry Optimization Has Not Converged. "
                      "Reached The Maximum Iteration",
            )
            raise StopIteration


class OQPTCIOpt(_OQPRunner, MECIOpt):
    r"""Legacy three-state CI search retained for input compatibility.

    This is the historical OpenQP TCI prototype: it penalizes the two
    consecutive gaps and multiplies ``sigma`` by ``pen_incre`` every step.
    The published additive Baek adaptive algorithm is independently selected
    with ``runtype=meci`` and ``meci_search=baeka``.
    """

    mode = "min"

    def __init__(self, mol):
        MECIOpt.__init__(self, mol)

    def one_step(self, coordinates):
        if not (self.istate < self.jstate < self.kstate):
            raise ValueError(
                "TCI requires istate < jstate < kstate, got "
                f"{self.istate}/{self.jstate}/{self.kstate}")

        self.itr += 1
        dump_log(self.mol, title="PyOQP: Geometry Optimization Step %s" % self.itr)
        do_init_scf = True if self.itr == 1 else self.init_scf

        self.mol.update_system(coordinates)
        energies = self.sp.energy(do_init_scf=do_init_scf)

        self.grad.grads = [self.istate, self.jstate, self.kstate]
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads
        energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
        self.mol.energies = energies
        self.mol.grads = grads

        return self._tci_penalty(coordinates, energies, grads)

    def _tci_penalty(self, coordinates, energies, grads):
        ei = energies[self.istate]
        ej = energies[self.jstate]
        ek = energies[self.kstate]
        gi = grads[self.istate].reshape(-1)
        gj = grads[self.jstate].reshape(-1)
        gk = grads[self.kstate].reshape(-1)

        g_ji = ej - ei
        g_kj = ek - ej
        dg_ji = gj - gi
        dg_kj = gk - gj

        if self.alpha == 0:
            scale = max(np.mean(dg_ji ** 2) ** 0.5,
                        np.mean(dg_kj ** 2) ** 0.5)
            alpha = scale if scale > 1.0e-12 else 1.0e-6
        else:
            alpha = self.alpha

        self.sigma *= self.incre

        def pen(d):
            return d * d / (d + alpha)

        def dpen(d):
            return (d * d + 2.0 * alpha * d) / (d + alpha) ** 2

        avg_e = (ei + ej + ek) / 3.0
        f = avg_e + self.weights * self.sigma * (pen(g_ji) + pen(g_kj))
        df_avg = (gi + gj + gk) / 3.0
        df_pen = dpen(g_ji) * dg_ji + dpen(g_kj) * dg_kj
        df = df_avg + self.weights * self.sigma * df_pen

        max_gap = max(abs(g_ji), abs(g_kj))
        de = f - self.pre_energy
        rmsd_step = np.mean((coordinates - self.pre_coord) ** 2) ** 0.5
        max_step = np.amax(np.abs(coordinates - self.pre_coord))
        rmsd_grad = np.mean(df ** 2) ** 0.5
        max_grad = np.amax(np.abs(df))

        self.metrics['itr'] = self.itr
        self.metrics['sigma'] = self.sigma
        self.metrics['de'] = de
        self.metrics['gap'] = max_gap
        self.metrics['rmsd_step'] = rmsd_step
        self.metrics['max_step'] = max_step
        self.metrics['rmsd_grad'] = rmsd_grad
        self.metrics['max_grad'] = max_grad

        self.pre_energy = f
        self.pre_coord = coordinates.copy()
        dump_data(
            self.mol,
            (self.itr, self.atoms, coordinates, f, de, max_gap, rmsd_step,
             max_step, rmsd_grad, max_grad, g_ji, g_kj),
            title='TCI',
            fpath=self.mol.log_path,
        )
        return f, df


class OQPNEBOpt(StateSpecificOpt):
    """Nudged elastic band reaction path via the oqp FIRE band optimizer.

    Reactant is the ``[input] system`` geometry; product is read from the
    ``[neb] product`` XYZ endpoint.  Images are linearly interpolated and the
    interior images optimized with improved-tangent NEB + optional climbing
    image (see :mod:`oqp.library.oqp_neb`).  Energies/gradients per image use
    the state-specific (``istate``) surface via the reused OQP machinery.
    """

    def __init__(self, mol):
        StateSpecificOpt.__init__(self, mol)
        neb_cfg = mol.config.get("neb", {})
        nat_cfg = mol.config.get("oqp", {})
        self.nimage = int(neb_cfg.get("nimage", 5))
        self.product = neb_cfg.get("product", "")
        self.k_spring = float(nat_cfg.get("spring", 0.05))
        self.climbing = bool(nat_cfg.get("climb", True))
        self.neb_fmax = float(nat_cfg.get("fmax", 2.0e-3))
        self.neb_frms = float(nat_cfg.get("frms", self.neb_fmax))
        # Relax-then-climb threshold, FIRE step, and per-image step cap.
        self.climb_fmax = float(nat_cfg.get("climb_fmax", 0.05))
        self.neb_dt = float(nat_cfg.get("neb_dt", 0.5))
        self.maxmove = float(nat_cfg.get("maxmove", 0.2))
        self.align = bool(nat_cfg.get("align", True))
        # Pre-optimize the endpoints to minima before building the band
        # (a relaxed reactant/product makes the climbing image far more stable).
        self.opt_ends = bool(nat_cfg.get("opt_ends", True))
        self.end_fmax = float(nat_cfg.get("end_fmax", 1.0e-3))
        self.neb_output = str(nat_cfg.get("neb_output", "")).strip()

    def _resolve_product(self):
        if os.path.isabs(self.product):
            return self.product
        candidates = [self.product]
        input_file = getattr(self.mol, "input_file", "")
        if input_file:
            candidates.append(os.path.join(os.path.dirname(input_file), self.product))
        for c in candidates:
            if os.path.exists(os.path.abspath(c)):
                return os.path.abspath(c)
        return os.path.abspath(candidates[-1])

    def _energy_gradient(self, coords):
        # Fresh state-specific energy+gradient for one image (Hartree, H/Bohr).
        self.mol.update_system(coords)
        energies = self.sp.energy(do_init_scf=True)
        self.grad.grads = [self.istate]
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads
        energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
        self.mol.energies = energies
        self.mol.grads = grads
        return float(energies[self.istate]), grads[self.istate].reshape(-1)

    def _relax_endpoint(self, x0, label):
        """Minimize an endpoint to a nearby minimum with the oqp engine."""
        atoms = np.asarray(self.mol.get_atoms(), dtype=int).reshape(-1)
        engine = OQPEngine(atoms, x0, mode="min", maxiter=self.maxit)
        last = {
            "x": np.asarray(x0, dtype=float).reshape(-1).copy(),
            "g": None,
            "converged": False,
        }

        def eg(x):
            e, g = self._energy_gradient(x)
            last["x"] = np.asarray(x, dtype=float).reshape(-1).copy()
            last["g"] = g
            return e, g

        def conv():
            if last["g"] is not None and np.max(np.abs(last["g"])) < self.end_fmax:
                last["converged"] = True
                raise StopIteration

        engine.run(eg, on_converged=conv)
        gmax = float(np.max(np.abs(last["g"]))) if last["g"] is not None else float("nan")
        if last["converged"]:
            dump_log(self.mol, title="PyOQP: NEB endpoint '%s' relaxed (gmax=%.3e)"
                     % (label, gmax))
        else:
            dump_log(
                self.mol,
                title=("PyOQP: WARNING NEB endpoint '%s' did not reach end_fmax "
                       "%.3e within %d iterations (gmax=%.3e); continuing with "
                       "the last evaluated geometry")
                % (label, self.end_fmax, self.maxit, gmax),
            )
        # OQPEngine may take an unevaluated trial step after the final allowed
        # evaluation.  Only return the last endpoint whose energy/gradient were
        # actually computed.
        return last["x"]

    def optimize(self):
        dump_log(self.mol, title="PyOQP: OQP NEB [%d images, climbing=%s, k=%.3f, opt_ends=%s]"
                 % (self.nimage, self.climbing, self.k_spring, self.opt_ends))

        # Endpoints: reactant from the molecule (Bohr), product from XYZ (Ang).
        reactant = np.asarray(self.pre_coord, dtype=float).reshape(-1)
        prod_img = _read_xyz(self._resolve_product())

        # The product must describe the SAME atoms in the SAME order as the
        # reactant -- matching atom count alone would silently interpolate a path
        # between mismatched atoms (e.g. HCN vs CNH) and corrupt the energies.
        reactant_z = [int(z) for z in np.asarray(self.mol.get_atoms(), dtype=int).reshape(-1)]
        try:
            product_z = [int(SYMBOL_MAP[s]) for s in prod_img.symbols]
        except KeyError as exc:
            raise ValueError(
                "NEB product endpoint has an unrecognized element symbol: %s" % exc)
        if product_z != reactant_z:
            raise ValueError(
                "NEB product endpoint atoms must match the reactant in element "
                "and order (got %d product atoms vs %d reactant atoms)"
                % (len(product_z), len(reactant_z)))

        product = np.asarray(prod_img.coordinates_angstrom, dtype=float).reshape(-1) \
            * ANGSTROM_TO_BOHR

        if self.opt_ends:
            reactant = self._relax_endpoint(reactant, "reactant")
            product = self._relax_endpoint(product, "product")

        if self.align:
            product = np.asarray(
                kabsch_align(reactant.reshape(-1, 3), product.reshape(-1, 3)),
                dtype=float,
            ).reshape(-1)

        # Linear interpolation of all images (Bohr).
        fractions = np.linspace(0.0, 1.0, self.nimage)
        images = [reactant + f * (product - reactant) for f in fractions]

        neb = NEB(images, k_spring=self.k_spring, climbing=self.climbing,
                  climb_fmax=self.climb_fmax)

        def log_iter(it, fmax, energies):
            rel = np.array(energies) - energies[0]
            self.metrics["itr"] = it + 1
            self.metrics["max_grad"] = fmax
            dump_log(self.mol,
                     title="NEB Iteration %d  fmax=%.3e  Emax(rel)=%.6f"
                     % (it + 1, fmax, float(np.max(rel))))

        result = neb.run(self._energy_gradient, fmax_tol=self.neb_fmax,
                         maxiter=self.maxit, dt=self.neb_dt,
                         maxmove=self.maxmove, on_iteration=log_iter,
                         frms_tol=self.neb_frms)

        # Report and dump the final path.
        e0 = result["energies"][0]
        for idx, (img, en) in enumerate(zip(neb.images, result["energies"])):
            dump_data(self.mol,
                      (idx, self.atoms, img.reshape(-1, 3), en, en - e0),
                      title="NEB", fpath=self.mol.log_path)
        if result["converged"]:
            dump_log(self.mol, title="PyOQP: NEB Has Converged (fmax=%.3e, frms=%.3e in %d iters)"
                     % (result["fmax"], result["frms"], result["iters"]))
        else:
            dump_log(self.mol, title="PyOQP: NEB Reached Max Iterations (fmax=%.3e, frms=%.3e)"
                     % (result["fmax"], result["frms"]))

        output = self.neb_output
        if not output:
            project = getattr(self.mol, "project_name", "openqp") or "openqp"
            output = os.path.join(getattr(self.mol, "log_path", os.getcwd()),
                                  "%s_neb.xyz" % project)
        elif not os.path.isabs(output):
            output = os.path.join(getattr(self.mol, "log_path", os.getcwd()), output)
        symbols = [ELEMENTS_NAME[int(z)].strip() for z in reactant_z]
        write_neb_xyz(output, symbols, result["images"], result["energies"])
        dump_log(self.mol, title="PyOQP: Native NEB path written to %s" % output)

        # Preserve the synchronized final band for programmatic consumers and
        # leave the molecule at the highest-energy interior image, the natural
        # starting point for a subsequent native TS search.
        self.mol.neb_images = [image.copy() for image in result["images"]]
        self.mol.neb_energies = list(result["energies"])
        self.mol.neb_gradients = [gradient.copy() for gradient in result["gradients"]]
        top = 1 + int(np.argmax(result["energies"][1:-1]))
        self.pre_coord = result["images"][top].copy()
        self.neb_result = result

        # Downstream compute_scf_prop and JSON serialization require the full
        # wavefunction, density, state spectrum, and gradients to match the
        # geometry left on mol. Re-evaluate that context at the top image, but
        # keep the synchronized band/result/XYZ immutable: their values and
        # force metrics belong to the NEB snapshot above.
        context_energy, _ = self._energy_gradient(self.pre_coord)
        self.pre_energy = float(context_energy)
        dump_log(
            self.mol,
            title=("PyOQP: Native NEB final electronic context synchronized at "
                   "image %d (band E=%.12f, context E=%.12f)")
            % (top + 1, result["energies"][top], context_energy),
        )


class OQPIRCOpt(StateSpecificOpt):
    """Intrinsic reaction coordinate by the oqp Gonzalez-Schlegel IRC.

    Computes the Hessian at the starting transition state, takes the
    imaginary-frequency mode as the initial (mass-weighted) direction, and
    traces the mass-weighted steepest-descent path forward or backward (see
    :mod:`oqp.library.oqp_irc`).
    """

    def __init__(self, mol):
        StateSpecificOpt.__init__(self, mol)
        nat_cfg = mol.config.get("oqp", {})
        geo_cfg = mol.config.get("geometric", {})
        self.irc_step = float(nat_cfg.get("irc_step", 0.10))
        direction = nat_cfg.get("irc_direction",
                                geo_cfg.get("irc_direction", "forward"))
        direction = str(direction).strip().lower()
        if direction == "reverse":
            direction = "backward"
        if direction not in {"forward", "backward"}:
            raise ValueError("oqp.irc_direction must be forward or backward")
        self.irc_sign = -1 if direction == "backward" else 1
        self.path_gtol = float(nat_cfg.get("path_gtol", 1.0e-4))

    def _energy_gradient(self, coords):
        self.mol.update_system(coords)
        energies = self.sp.energy(do_init_scf=True)
        self.grad.grads = [self.istate]
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads
        energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
        return float(energies[self.istate]), grads[self.istate].reshape(-1)

    def optimize(self):
        from oqp.library.single_point import Hessian
        from oqp.library.oqp_irc import IRC, imaginary_mode

        dump_log(self.mol, title="PyOQP: OQP IRC [direction=%s, step=%.3f]"
                 % ("backward" if self.irc_sign < 0 else "forward", self.irc_step))

        # Hessian at the transition state -> imaginary mode.  The Hessian class
        # reads mol.energies, so converge the SCF/reference at the TS first.
        self.mol.update_system(np.asarray(self.pre_coord, dtype=float).reshape(-1))
        self.mol.energies = self.sp.energy(do_init_scf=True)
        hess_cfg = self.mol.config["hess"]
        previous_state = hess_cfg.get("state", 0)
        previous_restart = hess_cfg.get("restart", False)
        if previous_restart:
            raise ValueError(
                "Native IRC does not accept hess.restart=true because numerical "
                "Hessian restart gradients do not yet carry geometry/model signatures; "
                "set restart=false"
            )
        try:
            hess_cfg["state"] = self.istate
            # A numerical IRC Hessian must never assemble filename-only restart
            # gradients left by another geometry/model. A versioned hess.read
            # sidecar remains allowed and is validated by Molecule.read_freqs.
            hess_cfg["restart"] = False
            hess = Hessian(self.mol).hessian(analysis=False)
        finally:
            hess_cfg["state"] = previous_state
            hess_cfg["restart"] = previous_restart
        if hess is None:
            raise RuntimeError("Native IRC Hessian calculation failed")
        hess = np.asarray(hess, dtype=float)
        mass = np.asarray(self.mol.get_mass(), dtype=float).reshape(-1)
        x_ts = np.asarray(self.pre_coord, dtype=float).reshape(-1)
        mode, curv = imaginary_mode(mass, hess, coordinates=x_ts)
        irc = IRC(mass, x_ts, mode, step=self.irc_step, sign=self.irc_sign)

        def log_point(k, energy, gmax):
            dump_log(self.mol, title="IRC Point %d  E=%.8f  gmax=%.3e"
                     % (k + 1, energy, gmax))

        result = irc.run(self._energy_gradient, max_points=self.maxit,
                         gtol=self.path_gtol,
                         on_point=log_point)

        e0 = result["points"][0]["energy"]
        for idx, p in enumerate(result["points"]):
            dump_data(self.mol,
                      (idx, self.atoms, p["x"].reshape(-1, 3),
                       p["energy"], p["energy"] - e0),
                      title="IRC", fpath=self.mol.log_path)
        if result["converged"]:
            if result.get("reason") == "gradient":
                title = "PyOQP: IRC Reached Gradient Convergence (%d points)" % result["npoint"]
            else:
                title = "PyOQP: IRC Bracketed A Minimum Basin (%d points)" % result["npoint"]
            dump_log(self.mol, title=title)
        else:
            dump_log(self.mol, title="PyOQP: IRC Stopped At Max Points (%d)"
                     % result["npoint"])
        self.irc_result = result


class OQPMEPOpt(StateSpecificOpt):
    """Minimum energy path: mass-weighted steepest descent from the start point.

    Reuses the Gonzalez-Schlegel constrained-hypersphere integrator
    (:mod:`oqp.library.oqp_irc`) but, unlike IRC, begins at an arbitrary
    (non-stationary) geometry along the steepest-descent direction
    ``-g`` instead of a transition-state imaginary mode, tracing the path down
    to the nearest minimum.  Matches the semantics of the mass-weighted
    constrained-sphere ``MEP``/``ConstrainOpt`` driver.
    """

    def __init__(self, mol):
        StateSpecificOpt.__init__(self, mol)
        nat_cfg = mol.config.get("oqp", {})
        self.mep_step = float(nat_cfg.get("mep_step",
                                          nat_cfg.get("irc_step", 0.10)))
        self.path_gtol = float(nat_cfg.get("path_gtol", 1.0e-4))

    def _energy_gradient(self, coords):
        self.mol.update_system(coords)
        energies = self.sp.energy(do_init_scf=True)
        self.grad.grads = [self.istate]
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads
        energies, grads = self.ls.compute(self.mol, grad_list=self.grad.grads)
        return float(energies[self.istate]), grads[self.istate].reshape(-1)

    def optimize(self):
        from oqp.library.oqp_irc import IRC

        dump_log(self.mol, title="PyOQP: OQP MEP [step=%.3f]" % self.mep_step)

        x0 = np.asarray(self.pre_coord, dtype=float).reshape(-1)
        mass = np.asarray(self.mol.get_mass(), dtype=float).reshape(-1)
        sqm = np.sqrt(np.repeat(mass, 3))

        # Initial direction: mass-weighted steepest descent at the start point.
        _, g0 = self._energy_gradient(x0)
        g_mw = np.asarray(g0, dtype=float).reshape(-1) / sqm
        gnorm = np.linalg.norm(g_mw)
        if gnorm < 1.0e-8:
            dump_log(self.mol, title="PyOQP: MEP start is already stationary; nothing to do")
            return
        direction = -g_mw / gnorm

        mep = IRC(mass, x0, direction, step=self.mep_step, sign=1)

        def log_point(k, energy, gmax):
            dump_log(self.mol, title="MEP Point %d  E=%.8f  gmax=%.3e"
                     % (k + 1, energy, gmax))

        result = mep.run(self._energy_gradient, max_points=self.maxit,
                         gtol=self.path_gtol,
                         on_point=log_point)

        e0 = result["points"][0]["energy"]
        for idx, p in enumerate(result["points"]):
            dump_data(self.mol,
                      (idx, self.atoms, p["x"].reshape(-1, 3),
                       p["energy"], p["energy"] - e0),
                      title="MEP", fpath=self.mol.log_path)
        if result["converged"]:
            if result.get("reason") == "gradient":
                title = "PyOQP: MEP Reached Gradient Convergence (%d points)" % result["npoint"]
            else:
                title = "PyOQP: MEP Bracketed A Minimum Basin (%d points)" % result["npoint"]
            dump_log(self.mol, title=title)
        else:
            dump_log(self.mol, title="PyOQP: MEP Stopped At Max Points (%d)"
                     % result["npoint"])
        self.mep_result = result
