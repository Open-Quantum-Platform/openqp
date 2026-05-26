"""geomeTRIC optimizer interface for OpenQP."""

from __future__ import annotations

import os

import numpy as np

from oqp.library.libscipy import MECIOpt, MECPOpt, StateSpecificOpt
from oqp.utils.file_utils import dump_log

BOHR_TO_ANGSTROM = 0.52917721092
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM


class GeometricEngine:
    """geomeTRIC custom engine that delegates energy/gradient calls to OpenQP."""

    def __init__(self, optimizer, molecule):
        self.optimizer = optimizer
        self.M = molecule
        self.stored_calcs = {}

    def calc(self, coords, dirname, read_data=False, copydir=None):
        coord_hash = hash(np.asarray(coords, dtype=float).tobytes())
        if coord_hash not in self.stored_calcs:
            self.stored_calcs[coord_hash] = self.calc_new(coords, dirname)
        return self.stored_calcs[coord_hash]

    def clearCalcs(self):
        self.stored_calcs = {}

    def detect_dft(self):
        return False

    def save_guess_files(self, dirname):
        """No-op hook used by geomeTRIC engines that persist SCF guesses."""

    def load_guess_files(self, dirname):
        """No-op hook used by geomeTRIC engines that restore SCF guesses."""

    def calc_new(self, coords, dirname):
        energy, gradient = self.optimizer.one_step(np.asarray(coords, dtype=float).reshape(-1))
        return {
            "energy": float(energy),
            "gradient": np.asarray(gradient, dtype=float).reshape(-1),
        }


class _GeometricRunner:
    """Mixin that runs a PyOQP optimizer objective through geomeTRIC."""

    def __init__(self, mol):
        self.geometric_config = mol.config.get("geometric", {})
        self.coordsys = self.geometric_config.get("coordsys", "tric")
        self.trust = self.geometric_config.get("trust", 0.1)
        self.tmax = self.geometric_config.get("tmax", 0.3)
        self.convergence_set = self.geometric_config.get("convergence_set", "GAU")
        self.prefix = self.geometric_config.get("prefix", "geometric")
        self.hessian = self.geometric_config.get("hessian", "never")

    def _optimizer_keywords(self):
        keywords = {}
        constraints_file = self.geometric_config.get("constraints_file", "")
        if constraints_file:
            if not os.path.isabs(constraints_file):
                mol = getattr(self, "mol", None)
                input_file = getattr(mol, "input_file", "")
                if input_file:
                    constraints_file = os.path.join(os.path.dirname(input_file), constraints_file)
            keywords["constraints"] = constraints_file
            keywords["enforce"] = self.geometric_config.get("enforce", 0.0)
            keywords["conmethod"] = self.geometric_config.get("conmethod", 0)
        return keywords

    def optimize(self):
        try:
            from geometric.molecule import Elements, Molecule
            from geometric.optimize import run_optimizer
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "geomeTRIC is required for [optimize] lib=geometric. "
                "Install it with `pip install geometric`."
            ) from exc

        geometric_molecule = Molecule()
        atoms = np.asarray(self.mol.get_atoms(), dtype=int).reshape(-1)
        coords_bohr = np.asarray(self.pre_coord, dtype=float).reshape((self.natom, 3))
        geometric_molecule.elem = [Elements[int(atom)] for atom in atoms]
        geometric_molecule.xyzs = [coords_bohr * BOHR_TO_ANGSTROM]
        geometric_molecule.comms = ["OpenQP geomeTRIC optimization"]
        geometric_molecule.build_topology()

        engine = GeometricEngine(self, geometric_molecule)
        log_path = getattr(self.mol, "log_path", os.getcwd())
        prefix = os.path.join(log_path, self.prefix)

        dump_log(self.mol, title="PyOQP: geomeTRIC Optimizer Started")
        optimizer_keywords = self._optimizer_keywords()
        progress = run_optimizer(
            customengine=engine,
            prefix=prefix,
            coordsys=self.coordsys,
            maxiter=self.maxit,
            trust=self.trust,
            tmax=self.tmax,
            hessian=self.hessian,
            convergence_set=self.convergence_set,
            convergence_energy=self.energy_shift,
            convergence_grms=self.rmsd_grad,
            convergence_gmax=self.max_grad,
            convergence_drms=self.rmsd_step * BOHR_TO_ANGSTROM,
            convergence_dmax=self.max_step * BOHR_TO_ANGSTROM,
            **optimizer_keywords,
        )

        if getattr(progress, "xyzs", None):
            final_coords = np.asarray(progress.xyzs[-1], dtype=float).reshape((self.natom, 3)) * ANGSTROM_TO_BOHR
            self.mol.update_system(final_coords)
            self.pre_coord = final_coords.reshape(-1)

        dump_log(self.mol, title="PyOQP: geomeTRIC Optimizer Finished")
        return progress


class GeometricOpt(_GeometricRunner, StateSpecificOpt):
    """State-specific geometry optimization driven by geomeTRIC."""

    def __init__(self, mol):
        StateSpecificOpt.__init__(self, mol)
        _GeometricRunner.__init__(self, mol)


class GeometricMECIOpt(_GeometricRunner, MECIOpt):
    """MECI penalty/UBP objective minimized by geomeTRIC."""

    def __init__(self, mol):
        MECIOpt.__init__(self, mol)
        _GeometricRunner.__init__(self, mol)


class GeometricMECPOpt(_GeometricRunner, MECPOpt):
    """MECP objective minimized by geomeTRIC."""

    def __init__(self, mol):
        MECPOpt.__init__(self, mol)
        _GeometricRunner.__init__(self, mol)


class GeometricTSOpt(_GeometricRunner, StateSpecificOpt):
    """Transition-state optimization driven by geomeTRIC."""

    def __init__(self, mol):
        StateSpecificOpt.__init__(self, mol)
        _GeometricRunner.__init__(self, mol)
        if self.hessian == "never":
            self.hessian = "first"

    def _optimizer_keywords(self):
        return {"transition": True}


class GeometricIRCOpt(_GeometricRunner, StateSpecificOpt):
    """Intrinsic reaction coordinate calculation driven by geomeTRIC."""

    def __init__(self, mol):
        StateSpecificOpt.__init__(self, mol)
        _GeometricRunner.__init__(self, mol)
        self.irc_direction = self.geometric_config.get("irc_direction", "forward").lower()
        if self.irc_direction == "reverse":
            self.irc_direction = "backward"
        if self.irc_direction not in {"forward", "backward"}:
            raise ValueError("geometric.irc_direction must be forward or backward")
        if self.hessian == "never":
            self.hessian = "first"
        if self.tmax > self.trust:
            self.tmax = self.trust

    def _optimizer_keywords(self):
        return {"irc": True, "irc_direction": self.irc_direction}


class GeometricNEBOpt:
    """Dependency-light scaffold for geomeTRIC NEB path calculations."""

    def __init__(self, mol):
        self.mol = mol
        self.neb_config = mol.config.get("neb", {})
        self.nimage = int(self.neb_config.get("nimage", 5))

    def plan_image_directories(self):
        """Return isolated per-image working directories for a NEB path."""
        log_path = getattr(self.mol, "log_path", os.getcwd())
        neb_root = os.path.join(log_path, "neb")
        return [
            os.path.join(neb_root, f"image_{image_index:03d}")
            for image_index in range(self.nimage)
        ]

    def optimize(self):
        raise NotImplementedError(
            "geomeTRIC NEB runner wiring is not implemented yet; "
            "the current scaffold only validates input and plans isolated image directories."
        )
