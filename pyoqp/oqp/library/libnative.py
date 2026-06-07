"""Native geometry optimizer backend for OpenQP.

A NumPy/SciPy-only alternative to the external ``geomeTRIC`` driver.  The
electronic-structure work, convergence test and logging are reused verbatim from
:class:`oqp.library.libscipy.StateSpecificOpt`; this module only supplies the
step-determination loop through :class:`oqp.library.native_engine.NativeEngine`
(redundant internal coordinates + restricted-step RFO / P-RFO with model-Hessian
BFGS/Bofill updates).

Selected with ``[optimize] lib=native``.  Supported runtypes:

* ``optimize`` -> :class:`NativeOpt`   (state-specific minimum)
* ``ts``       -> :class:`NativeTSOpt` (transition state, eigenvector following)

Options are read from the ``[native]`` input section (see
``oqp.molecule.oqpdata``).
"""

from __future__ import annotations

import numpy as np

import os

from oqp.library.libscipy import StateSpecificOpt, MECIOpt, MECPOpt
from oqp.library.native_engine import NativeEngine
from oqp.library.native_neb import NEB
from oqp.library.neb_utils import _read_xyz
from oqp.utils.file_utils import dump_log, dump_data

ANGSTROM_TO_BOHR = 1.0 / 0.52917721092


class _NativeRunner:
    """Mixin that drives a PyOQP ``one_step`` objective with NativeEngine."""

    mode = "min"

    def _native_config(self):
        cfg = self.mol.config.get("native", {})
        return {
            "coordsys": cfg.get("coordsys", "auto"),
            "trust": float(cfg.get("trust", 0.2)),
            "trust_max": float(cfg.get("trust_max", 0.5)),
            "follow_mode": int(cfg.get("follow", 0)),
        }

    def optimize(self):
        opts = self._native_config()
        atoms = np.asarray(self.mol.get_atoms(), dtype=int).reshape(-1)
        x0 = np.asarray(self.pre_coord, dtype=float).reshape(-1)

        engine = NativeEngine(
            atoms, x0,
            mode=self.mode,
            trust=opts["trust"],
            trust_max=opts["trust_max"],
            follow_mode=opts["follow_mode"],
            coordsys=opts["coordsys"],
            # NativeEngine's own loop is a backstop; the OQP convergence test
            # (which raises StopIteration at maxit) governs termination.
            maxiter=self.maxit,
        )
        dump_log(
            self.mol,
            title="PyOQP: Native optimizer [mode=%s, coordsys=%s, trust=%.3f]"
            % (self.mode, engine.coordsys, opts["trust"]),
        )

        def energy_gradient(coords):
            # one_step updates metrics/pre_coord/pre_energy and returns
            # (energy [Hartree], gradient [Hartree/Bohr]) for the active state.
            return self.one_step(np.asarray(coords, dtype=float).reshape(-1))

        try:
            engine.run(energy_gradient, on_converged=self.check_convergence)
        except StopIteration:
            pass


class NativeOpt(_NativeRunner, StateSpecificOpt):
    """State-specific minimum via redundant internals + restricted-step RFO."""

    mode = "min"

    def __init__(self, mol):
        super().__init__(mol)


class NativeTSOpt(_NativeRunner, StateSpecificOpt):
    """Transition-state search via partitioned RFO (eigenvector following)."""

    mode = "ts"

    def __init__(self, mol):
        super().__init__(mol)


class NativeMECIOpt(_NativeRunner, MECIOpt):
    """MECI search: minimize the penalty/UBP objective with the native engine.

    Reuses ``MECIOpt.one_step`` (which returns the penalty objective and its
    gradient) and ``MECIOpt.check_convergence`` (which adds the energy-gap
    criterion) verbatim -- only the step determination is native.
    """

    mode = "min"

    def __init__(self, mol):
        MECIOpt.__init__(self, mol)


class NativeMECPOpt(_NativeRunner, MECPOpt):
    """MECP search: minimize the gap-penalty objective with the native engine."""

    mode = "min"

    def __init__(self, mol):
        MECPOpt.__init__(self, mol)


class NativeTCIOpt(_NativeRunner, MECIOpt):
    r"""Three-state conical intersection (TCI) search by adaptive penalty.

    Generalizes the Levine-Martinez two-state penalty to three states
    (i < j < k) by penalizing the two consecutive gaps, following the adaptive
    penalty-function method for three-state CIs with MRSF-TDDFT
    (Lee, Back, et al., J. Phys. Chem. A 2021, 125, 9580):

        F  = (E_i + E_j + E_k) / 3 + sigma * [ P(g_ji) + P(g_kj) ]
        P(d)  = d^2 / (d + alpha)
        P'(d) = (d^2 + 2 alpha d) / (d + alpha)^2
        dF = (G_i + G_j + G_k) / 3
             + sigma * [ P'(g_ji) (G_j - G_i) + P'(g_kj) (G_k - G_j) ]

    with g_ji = E_j - E_i, g_kj = E_k - E_j.  ``sigma`` grows by ``pen_incre``
    each step (the adaptive tightening); ``alpha`` defaults to the RMS of the
    larger gap-gradient when ``pen_alpha == 0``.  Convergence requires BOTH
    consecutive gaps below ``energy_gap``.  Reuses ``MECIOpt`` infrastructure;
    only the objective is three-state.
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


class NativeNEBOpt(StateSpecificOpt):
    """Nudged elastic band reaction path via the native FIRE band optimizer.

    Reactant is the ``[input] system`` geometry; product is read from the
    ``[neb] product`` XYZ endpoint.  Images are linearly interpolated and the
    interior images optimized with improved-tangent NEB + optional climbing
    image (see :mod:`oqp.library.native_neb`).  Energies/gradients per image use
    the state-specific (``istate``) surface via the reused OQP machinery.
    """

    def __init__(self, mol):
        StateSpecificOpt.__init__(self, mol)
        neb_cfg = mol.config.get("neb", {})
        nat_cfg = mol.config.get("native", {})
        self.nimage = int(neb_cfg.get("nimage", 5))
        self.product = neb_cfg.get("product", "")
        self.k_spring = float(nat_cfg.get("spring", 0.05))
        self.climbing = bool(nat_cfg.get("climb", True))
        self.neb_fmax = float(nat_cfg.get("fmax", 2.0e-3))

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
        return float(energies[self.istate]), grads[self.istate].reshape(-1)

    def optimize(self):
        dump_log(self.mol, title="PyOQP: Native NEB [%d images, climbing=%s, k=%.3f]"
                 % (self.nimage, self.climbing, self.k_spring))

        # Endpoints: reactant from the molecule (Bohr), product from XYZ (Ang).
        reactant = np.asarray(self.pre_coord, dtype=float).reshape(-1)
        prod_img = _read_xyz(self._resolve_product())
        product = np.asarray(prod_img.coordinates_angstrom, dtype=float).reshape(-1) \
            * ANGSTROM_TO_BOHR
        if product.shape != reactant.shape:
            raise ValueError("NEB product endpoint atom count must match the reactant")

        # Linear interpolation of all images (Bohr).
        fractions = np.linspace(0.0, 1.0, self.nimage)
        images = [reactant + f * (product - reactant) for f in fractions]

        neb = NEB(images, k_spring=self.k_spring, climbing=self.climbing)

        def log_iter(it, fmax, energies):
            rel = np.array(energies) - energies[0]
            self.metrics["itr"] = it + 1
            self.metrics["max_grad"] = fmax
            dump_log(self.mol,
                     title="NEB Iteration %d  fmax=%.3e  Emax(rel)=%.6f"
                     % (it + 1, fmax, float(np.max(rel))))

        result = neb.run(self._energy_gradient, fmax_tol=self.neb_fmax,
                         maxiter=self.maxit, on_iteration=log_iter)

        # Report and dump the final path.
        e0 = result["energies"][0]
        for idx, (img, en) in enumerate(zip(neb.images, result["energies"])):
            dump_data(self.mol,
                      (idx, self.atoms, img.reshape(-1, 3), en, en - e0),
                      title="NEB", fpath=self.mol.log_path)
        if result["converged"]:
            dump_log(self.mol, title="PyOQP: NEB Has Converged (fmax=%.3e in %d iters)"
                     % (result["fmax"], result["iters"]))
        else:
            dump_log(self.mol, title="PyOQP: NEB Reached Max Iterations (fmax=%.3e)"
                     % result["fmax"])
        self.neb_result = result
