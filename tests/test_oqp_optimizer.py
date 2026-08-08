"""Tests for the oqp geometry optimizer (coordinates, engine, wiring).

The numerical core (oqp_coords, oqp_engine) is pure NumPy/SciPy and is
tested directly.  The dispatcher/validator wiring is checked with lightweight
stubs so the compiled OQP backend is not required, mirroring
``test_geometric_optimizer.py``.
"""
import importlib.util
import sys
import tempfile
import types
import unittest
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "pyoqp" / "oqp" / "library"
EXAMPLES_OPT = ROOT / "examples" / "OPT"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _load_oqp_modules():
    """Load oqp_coords and oqp_engine without the oqp backend."""
    package_names = ("oqp", "oqp.library")
    module_names = (
        "oqp.library.oqp_coords", "oqp.library.oqp_engine",
        "oqp.library.oqp_neb", "oqp.library.oqp_irc", "oqp.library.baeka",
    )
    saved = {name: sys.modules.get(name) for name in package_names + module_names}
    try:
        for pkg in package_names:
            module = types.ModuleType(pkg)
            module.__path__ = []
            sys.modules[pkg] = module
        nc = _load("oqp.library.oqp_coords", LIB / "oqp_coords.py")
        ne = _load("oqp.library.oqp_engine", LIB / "oqp_engine.py")
        nb = _load("oqp.library.oqp_neb", LIB / "oqp_neb.py")
        ni = _load("oqp.library.oqp_irc", LIB / "oqp_irc.py")
        ba = _load("oqp.library.baeka", LIB / "baeka.py")
        return nc, ne, nb, ni, ba
    finally:
        # Do not leave the backend-free package stubs in global import state:
        # later test modules may legitimately import the real Python package.
        for name, module in saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


NC, NE, NB, NI, BA = _load_oqp_modules()
NEB = NB.NEB
IRC = NI.IRC


def _restore_modules(saved):
    for name, module in saved.items():
        if module is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = module


def _load_baeka_file_utils():
    """Load the real BAEKA artifact writer with dependency-light stubs."""

    dependency_names = (
        "oqp", "oqp.molden", "oqp.molden.moldenwriter",
        "oqp.periodic_table", "oqp.runtime", "oqp.utils",
        "oqp.utils.constants", "oqp.utils.mpi_utils", "oqp.utils.qmmm",
        "oqp.utils.state_labels", "oqp.utils.dftb_trace",
    )
    saved = {name: sys.modules.get(name) for name in dependency_names}
    try:
        oqp = types.ModuleType("oqp")
        oqp.__path__ = []
        sys.modules["oqp"] = oqp

        molden = types.ModuleType("oqp.molden")
        molden.__path__ = []
        sys.modules["oqp.molden"] = molden
        moldenwriter = types.ModuleType("oqp.molden.moldenwriter")
        moldenwriter.write_frequency = lambda *args, **kwargs: ""
        sys.modules["oqp.molden.moldenwriter"] = moldenwriter

        periodic_table = types.ModuleType("oqp.periodic_table")
        periodic_table.SYMBOL_MAP = {1: "H"}
        periodic_table.ELEMENTS_NAME = {"H": "H"}
        sys.modules["oqp.periodic_table"] = periodic_table

        runtime = types.ModuleType("oqp.runtime")
        runtime.basis_search_paths = lambda: []
        sys.modules["oqp.runtime"] = runtime

        oqp_utils = types.ModuleType("oqp.utils")
        oqp_utils.__path__ = []
        sys.modules["oqp.utils"] = oqp_utils
        constants = types.ModuleType("oqp.utils.constants")
        constants.ANGSTROM_TO_BOHR = 1.0
        sys.modules["oqp.utils.constants"] = constants
        mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
        mpi_utils.mpi_dump = lambda function: function
        sys.modules["oqp.utils.mpi_utils"] = mpi_utils
        qmmm = types.ModuleType("oqp.utils.qmmm")
        qmmm.gradient_qmmm = lambda *args, **kwargs: None
        sys.modules["oqp.utils.qmmm"] = qmmm
        dftb_trace = types.ModuleType("oqp.utils.dftb_trace")
        dftb_trace.final_energy_annotation = lambda *args, **kwargs: ""
        dftb_trace.final_energy_header = lambda *args, **kwargs: ""
        dftb_trace.final_energy_label = lambda *args, **kwargs: ""
        sys.modules["oqp.utils.dftb_trace"] = dftb_trace
        state_labels = types.ModuleType("oqp.utils.state_labels")
        state_labels.format_calculation_request = lambda *args, **kwargs: ""
        state_labels.format_dftb_settings = lambda *args, **kwargs: ""
        state_labels.is_mrsf = lambda *args, **kwargs: False
        state_labels.public_method_name = lambda *args, **kwargs: ""
        state_labels.public_state_label = lambda *args, **kwargs: ""
        state_labels.spin_name = lambda *args, **kwargs: ""
        sys.modules["oqp.utils.state_labels"] = state_labels

        return _load(
            "_oqp_file_utils_baeka_runtime_under_test",
            ROOT / "pyoqp" / "oqp" / "utils" / "file_utils.py",
        )
    finally:
        _restore_modules(saved)


def _load_baeka_runtime(fake_engine, optimizer_base, dump_log, dump_data):
    """Load OQPBaekAOpt around fakes while retaining its production callback."""

    dependency_names = (
        "oqp", "oqp.library", "oqp.library.libscipy",
        "oqp.library.baeka", "oqp.library.oqp_engine",
        "oqp.library.oqp_neb", "oqp.library.neb_utils",
        "oqp.periodic_table", "oqp.utils", "oqp.utils.file_utils",
    )
    saved = {name: sys.modules.get(name) for name in dependency_names}
    try:
        oqp = types.ModuleType("oqp")
        oqp.__path__ = []
        sys.modules["oqp"] = oqp
        library = types.ModuleType("oqp.library")
        library.__path__ = []
        sys.modules["oqp.library"] = library
        oqp_utils = types.ModuleType("oqp.utils")
        oqp_utils.__path__ = []
        sys.modules["oqp.utils"] = oqp_utils

        libscipy = types.ModuleType("oqp.library.libscipy")
        libscipy.Optimizer = optimizer_base
        libscipy.StateSpecificOpt = type("StateSpecificOpt", (optimizer_base,), {})
        libscipy.MECIOpt = type("MECIOpt", (optimizer_base,), {})
        libscipy.MECPOpt = type("MECPOpt", (optimizer_base,), {})
        sys.modules["oqp.library.libscipy"] = libscipy
        sys.modules["oqp.library.baeka"] = BA

        oqp_engine = types.ModuleType("oqp.library.oqp_engine")
        oqp_engine.OQPEngine = fake_engine
        oqp_engine.parse_frozen_distance_spec = lambda value: []
        sys.modules["oqp.library.oqp_engine"] = oqp_engine
        oqp_neb = types.ModuleType("oqp.library.oqp_neb")
        oqp_neb.NEB = type("NEB", (), {})
        sys.modules["oqp.library.oqp_neb"] = oqp_neb
        neb_utils = types.ModuleType("oqp.library.neb_utils")
        neb_utils._read_xyz = lambda path: None
        neb_utils.kabsch_align = lambda reference, mobile: mobile
        neb_utils.write_neb_xyz = lambda *args, **kwargs: None
        sys.modules["oqp.library.neb_utils"] = neb_utils

        periodic_table = types.ModuleType("oqp.periodic_table")
        periodic_table.ELEMENTS_NAME = {"H": "H"}
        periodic_table.SYMBOL_MAP = {1: "H"}
        sys.modules["oqp.periodic_table"] = periodic_table
        file_utils = types.ModuleType("oqp.utils.file_utils")
        file_utils.dump_log = dump_log
        file_utils.dump_data = dump_data
        sys.modules["oqp.utils.file_utils"] = file_utils

        return _load(
            "_oqp_liboqp_baeka_runtime_under_test",
            LIB / "liboqp.py",
        )
    finally:
        _restore_modules(saved)


# --------------------------------------------------------------------------- #
class TestNativeRecoveryLadder(unittest.TestCase):
    def test_electronic_failure_rejects_geometry_and_walks_shared_ladder(self):
        logs = []

        class FakeEngine:
            instances = []

            def __init__(
                    self, atoms, x0, coordsys="auto", trust=0.2,
                    trust_max=0.5, frozen_distances=None, **kwargs):
                self.x = np.asarray(x0, dtype=float).reshape(-1)
                self.coordsys = str(coordsys).upper()
                self.trust = float(trust)
                self.trust_max = float(trust_max)
                self.frozen_distances = list(frozen_distances or [])
                self.__class__.instances.append(self)

            def run(self, energy_gradient, on_converged=None):
                try:
                    energy_gradient(self.x)
                except StopIteration:
                    pass
                return self.x

        class FakeOptimizer:
            def __init__(self, mol):
                cfg = mol.config["optimize"]
                self.mol = mol
                self.maxit = int(cfg["maxit"])
                self.itr = 0
                self.pre_energy = 0.0
                self.pre_coord = mol.get_system().copy()
                self.metrics = {}
                self.energy_shift = 1.0e-6
                self.energy_gap = 1.0e-5
                self.rmsd_step = 1.0e-3
                self.max_step = 2.0e-3
                self.rmsd_grad = 1.0e-4
                self.max_grad = 3.0e-4

        runtime = _load_baeka_runtime(
            FakeEngine,
            FakeOptimizer,
            lambda mol, title=None, **kwargs: logs.append(title or ""),
            lambda *args, **kwargs: None,
        )

        class FakeMolecule:
            def __init__(self):
                self.config = {
                    "input": {"qmmm_flag": False, "runtype": "optimize"},
                    "optimize": {"maxit": 30},
                    "oqp": {
                        "coordsys": "auto",
                        "trust": 0.2,
                        "trust_max": 0.5,
                        "auto_recovery": True,
                        "recovery_maxit": 2,
                        "recovery_trust": 0.02,
                        "freeze": "",
                        "follow": 0,
                    },
                }
                self.coordinates = np.zeros((30, 3))

            @staticmethod
            def get_atoms():
                return np.ones(30, dtype=int)

            @staticmethod
            def get_mass():
                return np.ones(30)

            def get_system(self):
                return self.coordinates.copy()

            def update_system(self, coordinates):
                self.coordinates = np.asarray(coordinates, dtype=float).reshape(30, 3)

        optimizer = runtime.OQPOpt(FakeMolecule())

        improving_oscillation = [
            {
                "rmsd_grad": value,
                "rmsd_step": step,
                "de": (-1.0 if index % 2 else 1.0) * 1.0e-4,
            }
            for index, (value, step) in enumerate([
                (0.0049, 0.0277), (0.0034, 0.0187),
                (0.0027, 0.0093), (0.0045, 0.0075),
                (0.0030, 0.0059), (0.0037, 0.0081),
                (0.0023, 0.0050), (0.0052, 0.0081),
            ])
        ]
        self.assertFalse(
            optimizer._native_progress_stalled(improving_oscillation)
        )
        flat_high_gradient = [
            {
                "rmsd_grad": 0.01,
                "rmsd_step": 0.001,
                "de": (-1.0 if index % 2 else 1.0) * 1.0e-4,
            }
            for index in range(8)
        ]
        self.assertTrue(
            optimizer._native_progress_stalled(flat_high_gradient)
        )
        self.assertTrue(optimizer._is_electronic_nonconvergence(
            RuntimeError(
                "openqp-dftb LC density-matrix CPHF did not converge; "
                "residual 5.1"
            )
        ))
        self.assertFalse(optimizer._is_electronic_nonconvergence(
            RuntimeError("openqp-dftb parameter set is incomplete")
        ))

        def fail_electronics(_coordinates):
            optimizer.itr += 1
            raise RuntimeError("SCF did not converge at this geometry")

        optimizer.one_step = fail_electronics
        optimizer.check_convergence = lambda: None
        optimizer.optimize()

        self.assertEqual(
            [engine.coordsys for engine in FakeEngine.instances],
            ["DLC", "CART", "DLC"],
        )
        self.assertEqual(
            [engine.trust for engine in FakeEngine.instances],
            [0.02, 0.02, 0.01],
        )
        self.assertTrue(any("Rejecting that geometry" in log for log in logs))
        self.assertTrue(any("Last electronic failure" in log for log in logs))


class TestBaekAKernel(unittest.TestCase):
    def test_three_state_objective_uses_only_adjacent_gaps(self):
        energies = np.array([1.0, 1.2, 1.5])
        gradients = np.array([[0.1, -0.2], [0.4, 0.3], [-0.1, 0.7]])
        sigma = 2.0
        alpha = 0.02
        result = BA.baeka_objective(
            energies, gradients, sigma=sigma, alpha=alpha
        )

        adjacent = (0.2 ** 2 / (0.2 + alpha)
                    + 0.3 ** 2 / (0.3 + alpha))
        expected = np.mean(energies) + sigma * adjacent
        self.assertAlmostEqual(result.objective, expected, places=14)
        self.assertTrue(np.allclose(result.gaps, [0.2, 0.3]))
        # The redundant outer gap (E2-E0) must not appear.
        with_outer = expected + sigma * 0.5 ** 2 / (0.5 + alpha)
        self.assertNotAlmostEqual(result.objective, with_outer, places=8)

    def test_generalized_objective_has_n_minus_one_gaps(self):
        energies = np.array([0.0, 0.1, 0.25, 0.7])
        gradients = np.arange(12.0).reshape(4, 3) / 10.0
        result = BA.baeka_objective(
            energies, gradients, sigma=1.0, alpha=0.02
        )
        self.assertEqual(result.gaps.shape, (3,))
        self.assertEqual(result.gap_gradients.shape, (3, 3))
        self.assertAlmostEqual(result.span, 0.7)

    def test_objective_gradient_matches_finite_difference(self):
        def surfaces(x):
            x0, x1 = x
            energies = np.array([
                0.10 * x0 * x0 + 0.05 * x1 * x1,
                0.40 + 0.20 * x0 * x0 + 0.10 * x1 * x1,
                0.90 + 0.05 * x0 * x0 + 0.30 * x1 * x1,
            ])
            gradients = np.array([
                [0.20 * x0, 0.10 * x1],
                [0.40 * x0, 0.20 * x1],
                [0.10 * x0, 0.60 * x1],
            ])
            return energies, gradients

        x = np.array([0.31, -0.22])
        energies, gradients = surfaces(x)
        analytic = BA.baeka_objective(
            energies, gradients, sigma=1.7, alpha=0.02
        ).gradient
        numeric = np.zeros(2)
        step = 1.0e-6
        for coordinate in range(2):
            xp = x.copy(); xp[coordinate] += step
            xm = x.copy(); xm[coordinate] -= step
            ep, gp = surfaces(xp)
            em, gm = surfaces(xm)
            fp = BA.baeka_objective(
                ep, gp, sigma=1.7, alpha=0.02
            ).objective
            fm = BA.baeka_objective(
                em, gm, sigma=1.7, alpha=0.02
            ).objective
            numeric[coordinate] = (fp - fm) / (2.0 * step)
        self.assertLess(np.max(np.abs(analytic - numeric)), 1.0e-9)

    def test_projector_retains_constraints_when_sum_cancels(self):
        gap_gradients = np.array([[1.0, 0.0], [-1.0, 0.0]])
        self.assertTrue(np.array_equal(np.sum(gap_gradients, axis=0), [0.0, 0.0]))
        projector, rank, _ = BA.penalty_projector(gap_gradients)
        self.assertEqual(rank, 1)
        self.assertTrue(np.allclose(projector, projector.T))
        self.assertTrue(np.allclose(projector @ projector, projector))
        self.assertTrue(np.allclose(projector @ np.array([2.0, 3.0]), [2.0, 0.0]))

    def test_sigma_updates_are_additive_and_same_sigma_delta_is_used(self):
        state = BA.BaekAState(
            sigma=1.0,
            alpha=0.02,
            delta_beta=0.025,
            beta_schedule="10,25",
            gap_tol=1.0e-4,
            tol_f=1.0e-8,
            tol_g=1.0e-8,
        )
        energies = np.array([0.0, 0.2])
        gradients = np.zeros((2, 3))

        first = state.evaluate(energies, gradients)
        self.assertEqual(first.action, "step")
        self.assertAlmostEqual(first.data.sigma, 1.0)
        self.assertAlmostEqual(first.next_sigma, 1.025)

        second = state.evaluate(energies, gradients)
        # Recombining the previous point at current sigma makes delta_F zero;
        # comparing F(sigma=1.0) directly with F(sigma=1.025) would not.
        self.assertAlmostEqual(second.delta_f, 0.0, places=15)
        self.assertTrue(second.local_stationary)
        self.assertEqual(second.action, "jump")
        self.assertAlmostEqual(second.tested_sigma, 1.025)
        self.assertAlmostEqual(second.jump, 10.0)
        self.assertAlmostEqual(second.tested_data.sigma, 1.025)
        self.assertAlmostEqual(second.data.sigma, 11.025)
        self.assertAlmostEqual(second.next_sigma, 11.025)

    def test_tight_gap_converges_without_large_jump(self):
        state = BA.BaekAState(
            sigma=1.0,
            alpha=0.02,
            delta_beta=0.025,
            beta_schedule=[10.0],
            gap_tol=1.0e-4,
            tol_f=1.0e-8,
            tol_g=1.0e-8,
        )
        energies = np.array([0.0, 5.0e-5])
        gradients = np.zeros((2, 3))
        state.evaluate(energies, gradients)
        result = state.evaluate(energies, gradients)
        self.assertTrue(result.converged)
        self.assertEqual(result.action, "converged")
        self.assertIsNone(result.jump)

    def test_invalid_energy_smoothing_parameter_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "alpha"):
            BA.baeka_objective(
                [0.0, 0.1], np.zeros((2, 3)), sigma=1.0, alpha=0.0
            )
        with self.assertRaisesRegex(ValueError, "pen_jump"):
            BA.parse_jump_schedule("10,,25")
        with self.assertRaisesRegex(ValueError, "gap_weight is fixed at 1.0"):
            BA.baeka_objective(
                [0.0, 0.1], np.zeros((2, 3)), sigma=1.0, alpha=0.02,
                gap_weight=2.0,
            )

    def test_runtime_rebases_current_sigma_and_appends_baeka_artifacts(self):
        file_utils = _load_baeka_file_utils()
        convergence_logs = []
        artifact_calls = []

        def dump_log(mol, title=None, section=None, info=None):
            if section == "baeka":
                convergence_logs.append(dict(info))

        def dump_data(mol, data, title=None, fpath="."):
            if title == "BAEKA":
                artifact_calls.append(data)
            return file_utils.dump_data(
                mol, data, title=title, fpath=fpath
            )

        class FakeEngine:
            instances = []

            def __init__(self, atoms, x0, coordsys="auto", **kwargs):
                self.x = np.asarray(x0, dtype=float).reshape(-1)
                self.coordsys = "cart" if coordsys == "auto" else coordsys
                self.rebases = []
                self.evaluations = []
                self.__class__.instances.append(self)

            def rebase_previous_objective(self, energy, gradient):
                self.rebases.append(
                    (float(energy), np.asarray(gradient, dtype=float).copy())
                )
                return True

            def run(self, energy_gradient, on_converged=None):
                coordinates = (self.x.copy(), self.x + np.array([0.1, 0.0, 0.0]))
                for point in coordinates:
                    self.evaluations.append(energy_gradient(point))
                    try:
                        if on_converged is not None:
                            on_converged()
                    except StopIteration:
                        break
                return coordinates[-1]

        class FakeOptimizer:
            def __init__(self, mol):
                cfg = mol.config["optimize"]
                self.mol = mol
                self.maxit = int(cfg["maxit"])
                self.rmsd_grad = float(cfg["rmsd_grad"])
                self.istate = int(cfg["istate"])
                self.jstate = int(cfg["jstate"])
                self.energy_shift = float(cfg["energy_shift"])
                self.energy_gap = float(cfg["energy_gap"])
                self.init_scf = bool(cfg["init_scf"])
                self.natom = 1
                self.atoms = np.asarray(mol.get_atoms(), dtype=int).reshape(1, 1)
                self.sp = mol.sp
                self.grad = mol.grad
                self.ls = mol.ls
                self.itr = 0
                self.pre_energy = 0.0
                self.pre_coord = mol.get_system().copy()
                self.metrics = {}

        runtime = _load_baeka_runtime(
            FakeEngine, FakeOptimizer, dump_log, dump_data
        )

        class FakeSinglePoint:
            def __init__(self, mol):
                self.mol = mol
                self.init_scf_calls = []

            def energy(self, do_init_scf=False):
                self.init_scf_calls.append(bool(do_init_scf))
                return self.mol.surface_energies.copy()

        class FakeGradient:
            def __init__(self, mol):
                self.mol = mol
                self.grads = []

            def gradient(self):
                return self.mol.surface_gradients.copy()

        class FakeLastStep:
            @staticmethod
            def compute(mol, grad_list=None):
                return mol.energies, mol.grads

        class FakeMolecule:
            def __init__(self, log_path):
                self.log_path = str(log_path)
                self.config = {
                    "input": {"qmmm_flag": False},
                    "tdhf": {"nstate": 3},
                    "oqp": {
                        "coordsys": "cart", "trust": 0.2,
                        "trust_max": 0.5, "follow": 0,
                        "auto_recovery": False,
                    },
                    "optimize": {
                        "lib": "oqp", "optimizer": "bfgs", "maxit": 2,
                        "istate": 1, "jstate": 2, "kstate": 3,
                        "states": [1, 2, 3], "energy_shift": 1.0e-12,
                        "energy_gap": 1.0e-4, "rmsd_grad": 1.0,
                        "init_scf": False, "pen_sigma": 1.0,
                        "pen_alpha": 0.02, "pen_delta": 0.025,
                        "pen_jump": "10,25", "gap_weight": 1.0,
                    },
                }
                self.data = {"natom": 1}
                self.coordinates = np.zeros(3)
                self.surface_energies = np.array([-1.0, 0.0, 0.2, 0.4])
                self.surface_gradients = np.array([
                    [[0.0, 0.0, 0.0]],
                    [[0.0, 0.02, -0.01]],
                    [[0.1, 0.02, -0.01]],
                    [[0.2, 0.02, -0.01]],
                ])
                self.sp = FakeSinglePoint(self)
                self.grad = FakeGradient(self)
                self.ls = FakeLastStep()

            @staticmethod
            def get_atoms():
                return np.array([1])

            def get_system(self):
                return self.coordinates.copy()

            @staticmethod
            def get_mass():
                return np.array([1.0])

            def update_system(self, coordinates):
                self.coordinates = np.asarray(coordinates, dtype=float).copy()

        with tempfile.TemporaryDirectory() as tmpdir:
            mol = FakeMolecule(tmpdir)
            optimizer = runtime.OQPBaekAOpt(mol)
            optimizer.optimize()

            self.assertEqual(mol.sp.init_scf_calls, [True, False])
            self.assertEqual(len(convergence_logs), 2)
            first, second = convergence_logs
            self.assertAlmostEqual(first["sigma"], 1.0)
            self.assertAlmostEqual(first["active_sigma"], 1.0)
            self.assertAlmostEqual(first["next_sigma"], 1.025)
            self.assertAlmostEqual(
                first["next_sigma"], first["sigma"] + 0.025
            )
            self.assertEqual(first["action"], "step")
            self.assertAlmostEqual(second["sigma"], first["next_sigma"])
            self.assertAlmostEqual(second["active_sigma"], 11.025)
            self.assertAlmostEqual(
                second["active_sigma"], second["sigma"] + 10.0
            )
            self.assertAlmostEqual(second["next_sigma"], second["active_sigma"])
            self.assertEqual(second["action"], "jump")

            engine = FakeEngine.instances[-1]
            self.assertEqual(len(engine.rebases), 1)
            selected_energies = mol.surface_energies[[1, 2, 3]]
            selected_gradients = mol.surface_gradients[[1, 2, 3]].reshape(3, -1)
            expected_rebase = BA.baeka_objective(
                selected_energies, selected_gradients,
                sigma=second["active_sigma"], alpha=0.02,
            )
            rebased_energy, rebased_gradient = engine.rebases[0]
            self.assertAlmostEqual(rebased_energy, expected_rebase.objective)
            self.assertTrue(np.allclose(rebased_gradient, expected_rebase.gradient))
            tested_objective = BA.baeka_objective(
                selected_energies, selected_gradients,
                sigma=second["sigma"], alpha=0.02,
            )
            self.assertNotAlmostEqual(rebased_energy, tested_objective.objective)

            self.assertEqual(len(artifact_calls), 2)
            self.assertEqual([row[-1] for row in artifact_calls], ["step", "jump"])
            status_lines = (Path(tmpdir) / "opt_status.txt").read_text().splitlines()
            status_rows = []
            for line in status_lines:
                fields = line.split()
                if fields and fields[0].isdigit():
                    status_rows.append(fields)
            self.assertEqual(len(status_rows), 2)
            self.assertEqual(
                [float(row[4]) for row in status_rows], [1.0, 1.025]
            )
            self.assertEqual(
                [float(row[5]) for row in status_rows], [1.0, 11.025]
            )
            self.assertEqual([row[10] for row in status_rows], ["step", "jump"])
            geometry_lines = (Path(tmpdir) / "opt_geom.xyz").read_text().splitlines()
            self.assertEqual(geometry_lines.count("1"), 2)


class TestNativeFrozenDistances(unittest.TestCase):
    def test_parser_and_engine_preserve_the_initial_distance(self):
        self.assertEqual(
            NE.parse_frozen_distance_spec("distance(1,2);r(2,3)"),
            [(1, 2), (2, 3)],
        )
        x0 = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        target = np.array([[0.2, 1.0, 0.0], [1.7, 1.0, 0.0]])
        engine = NE.OQPEngine(
            [1, 1], x0, coordsys="cart", trust=0.15, maxiter=10,
            frozen_distances=[(1, 2)],
        )
        seen = []

        def energy_gradient(x):
            xyz = np.asarray(x).reshape(-1, 3)
            seen.append(np.linalg.norm(xyz[0] - xyz[1]))
            displacement = xyz - target
            return 0.5 * float(np.sum(displacement ** 2)), displacement.reshape(-1)

        result = engine.run(energy_gradient).reshape(-1, 3)
        self.assertTrue(seen)
        self.assertLess(max(abs(distance - 2.0) for distance in seen), 1.0e-10)
        self.assertAlmostEqual(np.linalg.norm(result[0] - result[1]), 2.0, places=10)
        self.assertGreater(np.mean(result[:, 1]), 0.5)

    def test_constraint_normal_is_removed_from_gradient(self):
        engine = NE.OQPEngine(
            [1, 1], [[0, 0, 0], [2, 0, 0]], coordsys="cart",
            frozen_distances=[(1, 2)],
        )
        normal = engine._constraint_jacobian(engine.x)[0]
        projected = engine._project_constraint_tangent(normal)
        self.assertLess(np.linalg.norm(projected), 1.0e-12)


# --------------------------------------------------------------------------- #
class TestInternalCoordinates(unittest.TestCase):
    def _fd_brow(self, prim, x):
        x = x.reshape(-1, 3)
        ana = np.zeros_like(x)
        for atom, d in prim.derivatives(x):
            ana[atom] += d
        num = np.zeros_like(x)
        h = 1.0e-6
        for a in range(x.shape[0]):
            for c in range(3):
                xp = x.copy(); xp[a, c] += h
                xm = x.copy(); xm[a, c] -= h
                dv = prim.value(xp) - prim.value(xm)
                dv = (dv + np.pi) % (2.0 * np.pi) - np.pi
                num[a, c] = dv / (2.0 * h)
        return np.max(np.abs(ana - num))

    def test_primitive_derivatives_match_fd(self):
        rng = np.random.default_rng(0)
        for seed in range(4):
            x = rng.normal(size=(4, 3)) * 1.5
            self.assertLess(self._fd_brow(NC.Bond(0, 1), x), 1e-6)
            self.assertLess(self._fd_brow(NC.Angle(0, 1, 2), x), 1e-6)
            self.assertLess(self._fd_brow(NC.Dihedral(0, 1, 2, 3), x), 1e-6)

    def test_collapsed_primitive_is_finite_and_warning_free(self):
        x = np.zeros((3, 3))
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            bond = NC.Bond(0, 1)
            angle = NC.Angle(0, 1, 2)
            bond_deriv = bond.derivatives(x)
            angle_value = angle.value(x)
            angle_deriv = angle.derivatives(x)
        # The invalid derivative row is intentional: it forces an existing
        # internal-coordinate engine into bounded Cartesian recovery instead
        # of silently accepting a zero gradient direction.
        self.assertFalse(np.all(np.isfinite([d for _, d in bond_deriv])))
        self.assertTrue(np.isfinite(angle_value))
        self.assertFalse(np.all(np.isfinite([d for _, d in angle_deriv])))

    def test_bmatrix_matches_fd_water_and_ethane(self):
        cases = [
            ([8, 1, 1], np.array([[0, 0, 0.12], [0, 1.43, -0.96],
                                  [0, -1.43, -0.96]], float)),
            ([6, 6, 1, 1, 1, 1, 1, 1],
             np.array([[0, 0, 1.45], [0, 0, -1.45], [1.7, 0, 2.2],
                       [-0.85, 1.47, 2.2], [-0.85, -1.47, 2.2],
                       [-1.7, 0, -2.2], [0.85, 1.47, -2.2],
                       [0.85, -1.47, -2.2]], float)),
        ]
        for atoms, x in cases:
            ic = NC.build_coordinates(atoms, x)
            self.assertIsInstance(ic, NC.DelocalizedInternalCoordinates)
            B = ic.b_matrix(x)
            h = 1.0e-6
            Bn = np.zeros_like(B)
            xf = x.reshape(-1)
            for k in range(len(xf)):
                xp = xf.copy(); xp[k] += h
                xm = xf.copy(); xm[k] -= h
                Bn[:, k] = ic.q_displacement(ic.q(xp), ic.q(xm)) / (2.0 * h)
            self.assertLess(np.max(np.abs(B - Bn)), 1e-6)

    def test_linear_molecule_falls_back_to_cartesian(self):
        # CO2 is linear: neither RIC nor (full-rank) TRIC can span it.
        for cs in ("ric", "tric", "dlc"):
            co2 = NC.build_coordinates([6, 8, 8],
                                       np.array([[0, 0, 0], [0, 0, 2.2],
                                                 [0, 0, -2.2]], float), coordsys=cs)
            self.assertIsInstance(co2, NC.CartesianCoordinates)


class TestTRICandDLC(unittest.TestCase):
    WATER_AT = [8, 1, 1]
    WATER_X = np.array([[0, 0, 0.12], [0, 1.43, -0.96], [0, -1.43, -0.96]], float)

    def _bmatrix_fd(self, ic, x):
        B = ic.b_matrix(x)
        Bn = np.zeros_like(B)
        xf = x.reshape(-1)
        h = 1.0e-6
        for k in range(len(xf)):
            xp = xf.copy(); xp[k] += h
            xm = xf.copy(); xm[k] -= h
            Bn[:, k] = ic.q_displacement(ic.q(xp.reshape(-1, 3)),
                                         ic.q(xm.reshape(-1, 3))) / (2.0 * h)
        return np.max(np.abs(B - Bn))

    def test_rotation_primitive_bmatrix(self):
        rng = np.random.default_rng(0)
        ref = rng.normal(size=(4, 3))
        rot = NC.RotationComponent([0, 1, 2, 3], 1, ref)
        x = rng.normal(size=(4, 3))
        ana = np.zeros((4, 3))
        for a, d in rot.derivatives(x):
            ana[a] = d
        num = np.zeros((4, 3))
        h = 1.0e-6
        for a in range(4):
            for c in range(3):
                xp = x.copy(); xp[a, c] += h
                xm = x.copy(); xm[a, c] -= h
                num[a, c] = (rot.value(xp) - rot.value(xm)) / (2.0 * h)
        self.assertLess(np.max(np.abs(ana - num)), 1e-6)

    def test_tric_full_rank_and_bmatrix(self):
        ic = NC.build_coordinates(self.WATER_AT, self.WATER_X, coordsys="tric")
        self.assertIsInstance(ic, NC.RedundantInternalCoordinates)
        # 2 bonds + 1 angle + 3 translations + 3 rotations
        self.assertEqual(len(ic.primitives), 9)
        self.assertLess(self._bmatrix_fd(ic, self.WATER_X), 1e-6)

    def test_dlc_nonredundant_and_bmatrix(self):
        ic = NC.build_coordinates(self.WATER_AT, self.WATER_X, coordsys="dlc")
        self.assertIsInstance(ic, NC.DelocalizedInternalCoordinates)
        self.assertEqual(len(ic.q(self.WATER_X)), 3)  # 3N-6 active coords
        self.assertLess(self._bmatrix_fd(ic, self.WATER_X), 1e-6)

    def test_auto_and_omitted_coordsys_select_dlc(self):
        for coordsys in (None, "auto"):
            ic = NC.build_coordinates(
                self.WATER_AT, self.WATER_X, coordsys=coordsys)
            self.assertIsInstance(ic, NC.DelocalizedInternalCoordinates)
            self.assertEqual(len(ic.q(self.WATER_X)), 3)

        engine = NE.OQPEngine(self.WATER_AT, self.WATER_X)
        self.assertIsInstance(engine.coords, NC.DelocalizedInternalCoordinates)
        self.assertEqual(engine.coordsys, "DLC")

    def test_water_dimer_dlc_spans_all_twelve_vibrational_modes(self):
        second = self.WATER_X + np.array([0.25, 0.10, 5.50])
        geometry = np.vstack([self.WATER_X, second])
        atoms = self.WATER_AT + self.WATER_AT

        _, _, fragments = NC._covalent_graph(atoms, geometry)
        self.assertEqual(len(fragments), 2)

        ic = NC.build_coordinates(atoms, geometry)
        self.assertIsInstance(ic, NC.DelocalizedInternalCoordinates)
        self.assertEqual(len(ic.q(geometry)), 12)
        _, rank = ic._g_inverse(ic.b_matrix(geometry))
        self.assertEqual(rank, 12)

    def test_rank_deficient_dimer_primitives_are_not_silently_accepted(self):
        second = self.WATER_X + np.array([0.25, 0.10, 5.50])
        geometry = np.vstack([self.WATER_X, second])
        atoms = self.WATER_AT + self.WATER_AT
        bonds, neighbours, fragments = NC._covalent_graph(atoms, geometry)
        self.assertEqual(len(fragments), 2)
        intrafragment = NC.RedundantInternalCoordinates(
            NC._make_internals(bonds, neighbours, geometry), len(atoms))

        with self.assertRaisesRegex(ValueError, "rank-deficient DLC basis"):
            NC.DelocalizedInternalCoordinates.from_ric(
                intrafragment, geometry)

    def test_ill_conditioned_dimer_gets_direct_interfragment_distance(self):
        second = self.WATER_X + np.array([0.25, 0.10, 5.50])
        geometry = np.vstack([self.WATER_X, second])
        atoms = self.WATER_AT + self.WATER_AT
        bonds, neighbours, fragments = NC._covalent_graph(atoms, geometry)
        bonds = NC._connect_fragments(
            bonds, neighbours, geometry, len(atoms))
        bonds, neighbours = NC._augment_internals(
            atoms, geometry, bonds, neighbours, len(atoms))
        primitives = NC._make_internals(bonds, neighbours, geometry)
        _, quality_before = NC._primitive_metric_quality(
            primitives, len(atoms), geometry)

        conditioned = NC._augment_interfragment_distances(
            primitives, fragments, geometry, len(atoms), quality_tol=0.03)
        rank_after, quality_after = NC._primitive_metric_quality(
            conditioned, len(atoms), geometry)

        self.assertGreater(len(conditioned), len(primitives))
        self.assertEqual(rank_after, 12)
        self.assertGreater(quality_after, quality_before)

    def test_dlc_linear_algebra_is_warning_free(self):
        """Coordinate products must not leak handled Accelerate FP flags."""
        ic = NC.build_coordinates(self.WATER_AT, self.WATER_X, coordsys="dlc")
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            q = ic.q(self.WATER_X)
            b = ic.b_matrix(self.WATER_X)
            gq = ic.grad_to_q(self.WATER_X, np.ones(9))
            hessian = ic.guess_hessian(self.WATER_X)
        self.assertTrue(np.all(np.isfinite(q)))
        self.assertTrue(np.all(np.isfinite(b)))
        self.assertTrue(np.all(np.isfinite(gq)))
        self.assertTrue(np.all(np.isfinite(hessian)))

    def test_dlc_g_inverse_rejects_overflowing_metric(self):
        """An overflowing B B^T metric must request coordinate fallback."""
        ic = NC.build_coordinates(self.WATER_AT, self.WATER_X, coordsys="dlc")
        b = np.full_like(ic.b_matrix(self.WATER_X), 1.0e308)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            g_inverse, rank = ic._g_inverse(b)
        self.assertEqual(rank, 0)
        self.assertTrue(np.array_equal(g_inverse, np.zeros_like(g_inverse)))

    def test_dlc_g_inverse_rejects_overflowing_symmetrization(self):
        """Finite G must remain safe when G + G.T would overflow."""
        ic = NC.build_coordinates(self.WATER_AT, self.WATER_X, coordsys="dlc")
        b = np.full_like(ic.b_matrix(self.WATER_X), 9.0e153)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            g_inverse, rank = ic._g_inverse(b)
        self.assertEqual(rank, 0)
        self.assertTrue(np.array_equal(g_inverse, np.zeros_like(g_inverse)))

    def test_dlc_rank_zero_gradient_requests_recovery(self):
        ic = NC.build_coordinates(self.WATER_AT, self.WATER_X, coordsys="dlc")
        b = np.full_like(ic.b_matrix(self.WATER_X), 9.0e153)
        ic.b_matrix = lambda _x: b
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            gradient = ic.grad_to_q(self.WATER_X, np.ones(9))
        self.assertTrue(np.all(np.isnan(gradient)))

    def test_dlc_partial_rank_gradient_requests_recovery(self):
        ic = NC.build_coordinates(self.WATER_AT, self.WATER_X, coordsys="dlc")
        ncoord = len(ic.q(self.WATER_X))
        original_inverse = ic._g_inverse
        try:
            ic._g_inverse = lambda _b: (np.eye(ncoord), ncoord - 1)
            gradient = ic.grad_to_q(self.WATER_X, np.ones(9))
        finally:
            ic._g_inverse = original_inverse
        self.assertTrue(np.all(np.isnan(gradient)))

    def test_safe_matmul_suppresses_handled_overflow(self):
        left = np.full((2, 2), 1.0e308)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            product = NC._safe_matmul(left, left)
        self.assertFalse(np.all(np.isfinite(product)))

    def test_internal_back_transform_rejects_overflow_without_poisoning_geometry(self):
        for coordsys in ("ric", "dlc"):
            with self.subTest(coordsys=coordsys):
                ic = NC.build_coordinates(
                    self.WATER_AT, self.WATER_X, coordsys=coordsys
                )
                huge_step = np.full(len(ic.q(self.WATER_X)), 1.0e308)
                with warnings.catch_warnings():
                    warnings.simplefilter("error", RuntimeWarning)
                    x_new, ok = ic.back_transform(self.WATER_X, huge_step)
                self.assertFalse(ok)
                self.assertTrue(np.all(np.isfinite(x_new)))
                self.assertTrue(np.allclose(x_new, self.WATER_X.reshape(-1)))

    def test_multifragment_tric(self):
        # Two well-separated triangular C3 fragments (each non-collinear, so
        # rotation is well defined): TRIC must span the full 3N = 18.
        tri = np.array([[0, 0, 0], [2.0, 0, 0], [1.0, 1.6, 0.0]], float)
        x = np.vstack([tri, tri + np.array([0.0, 0.0, 12.0])])
        ic = NC.build_coordinates([6, 6, 6, 6, 6, 6], x, coordsys="tric")
        self.assertIsInstance(ic, NC.RedundantInternalCoordinates)
        _, rank = ic._g_inverse(ic.b_matrix(x))
        self.assertEqual(rank, 18)

    def test_convergence_each_coordsys(self):
        def eg(xflat):
            x = xflat.reshape(-1, 3)
            pairs = [(0, 1), (0, 2), (1, 2)]
            E = 0.0
            g = np.zeros_like(x)
            for (i, j) in pairs:
                d = x[i] - x[j]
                r = np.linalg.norm(d)
                ex = np.exp(-(r - 2.0))
                E += 0.15 * (1 - ex) ** 2
                u = d / r
                f = 2 * 0.15 * (1 - ex) * ex
                g[i] += f * u
                g[j] -= f * u
            return E, g.reshape(-1)
        rng = np.random.default_rng(1)
        x0 = np.array([[0, 0, 0], [2.3, 0, 0], [1.0, 1.6, 0.0]], float) \
            + rng.normal(scale=0.1, size=(3, 3))
        for cs in ("ric", "dlc", "tric"):
            eng = NE.OQPEngine([6, 6, 6], x0.copy(), mode="min",
                                   maxiter=100, coordsys=cs)

            def stop():
                _, g = eg(eng.x)
                if np.max(np.abs(g)) < 1e-5:
                    raise StopIteration
            eng.run(eg, on_converged=stop)
            _, g = eg(eng.x)
            self.assertLess(np.max(np.abs(g)), 1e-5, f"coordsys={cs} failed")


# --------------------------------------------------------------------------- #
def _morse(pairs, r_e=2.0, De=0.15, a=1.0):
    def eg(xflat):
        x = xflat.reshape(-1, 3)
        E = 0.0
        g = np.zeros_like(x)
        for (i, j) in pairs:
            d = x[i] - x[j]
            r = np.linalg.norm(d)
            ex = np.exp(-a * (r - r_e))
            E += De * (1.0 - ex) ** 2
            dEdr = 2.0 * De * (1.0 - ex) * (a * ex)
            u = d / r
            g[i] += dEdr * u
            g[j] -= dEdr * u
        return E, g.reshape(-1)
    return eg


class TestOQPEngineConverges(unittest.TestCase):
    def test_changed_objective_rebases_previous_point_without_losing_history(self):
        engine = NE.OQPEngine(
            [1], [0.0, 0.0, 0.0], coordsys="cart",
            project_global_rigid_modes=False,
        )
        self.assertFalse(
            engine.rebase_previous_objective(0.0, np.zeros(3))
        )

        initial_x = engine.x.copy()
        engine._take_step(1.0, np.array([0.20, -0.10, 0.05]))
        hessian_before = engine.H.copy()
        previous_dq = engine._prev["dq"].copy()

        rebased_gradient = np.array([-0.30, 0.25, 0.10])
        self.assertTrue(
            engine.rebase_previous_objective(2.5, rebased_gradient)
        )
        expected_gradient_q = engine.coords.grad_to_q(
            initial_x, rebased_gradient
        )
        expected_prediction = float(
            expected_gradient_q @ previous_dq
            + 0.5 * previous_dq @ hessian_before @ previous_dq
        )

        self.assertTrue(np.array_equal(engine.H, hessian_before))
        self.assertTrue(np.array_equal(engine._prev["x"], initial_x))
        self.assertAlmostEqual(engine._prev["e"], 2.5)
        self.assertTrue(np.allclose(engine._prev["g_q"], expected_gradient_q))
        self.assertAlmostEqual(engine._prev["pred"], expected_prediction)

    def _optimize(self, atoms, x0, pairs, gtol=1e-4, maxiter=100):
        eg = _morse(pairs)
        eng = NE.OQPEngine(atoms, x0, mode="min", maxiter=maxiter)

        def stop():
            _, g = eg(eng.x)
            if np.max(np.abs(g)) < gtol:
                raise StopIteration

        eng.run(eg, on_converged=stop)
        _, g = eg(eng.x)
        return eng, np.max(np.abs(g))

    def test_morse_triangle(self):
        rng = np.random.default_rng(1)
        x0 = np.array([[0, 0, 0], [2.3, 0, 0], [1.0, 1.6, 0.0]], float)
        x0 = x0 + rng.normal(scale=0.1, size=(3, 3))
        eng, gmax = self._optimize([6, 6, 6], x0,
                                   [(0, 1), (0, 2), (1, 2)])
        self.assertLess(gmax, 1e-4)
        x = eng.x.reshape(-1, 3)
        bl = [np.linalg.norm(x[i] - x[j]) for (i, j) in [(0, 1), (0, 2), (1, 2)]]
        self.assertTrue(np.allclose(bl, 2.0, atol=1e-2))

    def test_morse_chain(self):
        rng = np.random.default_rng(2)
        x0 = np.array([[i * 2.2, 0, 0] for i in range(5)], float)
        x0 = x0 + rng.normal(scale=0.15, size=(5, 3))
        _, gmax = self._optimize([6] * 5, x0, [(i, i + 1) for i in range(4)])
        self.assertLess(gmax, 1e-4)

    def test_methane_like(self):
        # Tetrahedral CH4-like start (well separated H, so coordinate
        # connectivity matches the Morse bond list), distorted then relaxed.
        rng = np.random.default_rng(3)
        tet = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
                       float)
        tet = tet / np.linalg.norm(tet[0]) * 2.05
        x0 = np.vstack([[0, 0, 0], tet]) + rng.normal(scale=0.12, size=(5, 3))
        eng, gmax = self._optimize([6, 1, 1, 1, 1], x0,
                                   [(0, j) for j in range(1, 5)], maxiter=120)
        self.assertLess(gmax, 1e-4)
        x = eng.x.reshape(-1, 3)
        bl = [np.linalg.norm(x[0] - x[j]) for j in range(1, 5)]
        self.assertTrue(np.allclose(bl, 2.0, atol=1e-2))


class TestOQPEngineInitialHessian(unittest.TestCase):
    WATER_AT = [8, 1, 1]
    WATER_X = np.array([[0, 0, 0.12], [0, 1.43, -0.96],
                        [0, -1.43, -0.96]], float)

    def test_finite_failed_back_transform_preserves_internal_model(self):
        eng = NE.OQPEngine(self.WATER_AT, self.WATER_X, coordsys="dlc",
                           trust=0.1)
        original = eng.x.copy()
        eng.coords.back_transform = lambda x, dq: (np.asarray(x).copy(), False)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            x_new = eng._take_step(0.0, np.linspace(-0.1, 0.1, 9))
        self.assertTrue(np.all(np.isfinite(x_new)))
        self.assertFalse(np.allclose(x_new, original))
        self.assertIsNotNone(eng._prev)
        self.assertEqual(eng.nonfinite_step_rejections, 0)
        self.assertIsInstance(eng.coords, NC.DelocalizedInternalCoordinates)
        self.assertNotIn("CART(nonfinite recovery)", eng.coordsys)

    def test_nonfinite_hessian_switches_to_cartesian_recovery(self):
        eng = NE.OQPEngine(self.WATER_AT, self.WATER_X, coordsys="dlc",
                           trust=0.1)
        eng.H.fill(np.inf)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            x_new = eng._take_step(0.0, np.linspace(-0.1, 0.1, 9))
        self.assertTrue(np.all(np.isfinite(x_new)))
        self.assertEqual(eng.nonfinite_step_rejections, 1)
        self.assertIsInstance(eng.coords, NC.CartesianCoordinates)

    def test_collapsed_internal_geometry_switches_to_cartesian_recovery(self):
        eng = NE.OQPEngine(self.WATER_AT, self.WATER_X.copy(), coordsys="dlc",
                           trust=0.1)
        eng.x[:] = 0.0
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            x_new = eng._take_step(0.0, np.ones(9))
        self.assertTrue(np.all(np.isfinite(x_new)))
        self.assertEqual(eng.nonfinite_step_rejections, 1)
        self.assertIsInstance(eng.coords, NC.CartesianCoordinates)
        self.assertFalse(np.allclose(x_new, np.zeros(9)))

    def test_trust_restriction_scales_huge_finite_step_without_overflow(self):
        eng = NE.OQPEngine([1], [0.0, 0.0, 0.0], coordsys="cart",
                           trust=0.1)
        huge_step = np.array([1.0e200, -5.0e199, 2.5e199])
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            restricted = eng._restrict_to_trust(np.eye(3), huge_step)
        self.assertTrue(np.all(np.isfinite(restricted)))
        self.assertGreater(np.max(np.abs(restricted)), 0.0)
        self.assertLessEqual(NC._scaled_rms(restricted), eng.trust)

    def test_ts_cartesian_recovery_preserves_eigenvector_following(self):
        h_cart = np.diag([-0.5, 0.4, 0.8])
        eng = NE.OQPEngine(
            [1], [0.0, 0.0, 0.0], mode="ts", coordsys="cart",
            trust=0.1, initial_hessian=h_cart,
        )
        eng.coords.back_transform = lambda x, dq: (np.asarray(x).copy(), False)
        gradient = np.array([0.1, 0.2, 0.3])
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            x_new = eng._take_step(0.0, gradient)
        displacement = x_new - eng.x
        self.assertGreater(displacement[0] * gradient[0], 0.0)
        self.assertLess(np.dot(displacement[1:], gradient[1:]), 0.0)
        self.assertTrue(np.all(np.isfinite(x_new)))

    def test_meci_objective_rebase_drops_nonfinite_secant_history(self):
        eng = NE.OQPEngine(self.WATER_AT, self.WATER_X, coordsys="dlc")
        eng._prev = {
            "x": eng.x.copy(), "dq": np.ones(3), "g_q": np.zeros(3),
            "e": 0.0, "pred": 0.0, "cart_step": 0.01,
        }
        eng.H.fill(np.inf)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            rebased = eng.rebase_previous_objective(
                0.0, np.linspace(-0.1, 0.1, 9)
            )
        self.assertFalse(rebased)
        self.assertIsNone(eng._prev)

    def test_omitted_hessian_keeps_existing_model(self):
        eng = NE.OQPEngine(self.WATER_AT, self.WATER_X, coordsys="dlc")

        expected = eng.coords.guess_hessian(eng.x)

        self.assertTrue(np.array_equal(eng.H, expected))

    def test_cartesian_hessian_regularizes_only_global_external_subspace(self):
        h_cart = np.diag(np.linspace(-0.4, 0.8, 9))

        eng = NE.OQPEngine(
            self.WATER_AT, self.WATER_X, coordsys="cart",
            initial_hessian=h_cart, project_global_rigid_modes=True,
        )
        external = eng._global_external_mass_weighted_basis()
        p_external = external @ external.T
        p_internal = np.eye(9) - p_external
        model = np.eye(9) * 0.5
        expected = (
            p_internal @ h_cart @ p_internal
            + p_external @ model @ p_external
        )

        self.assertTrue(np.allclose(eng.H, expected, atol=1e-12))

    def test_cartesian_hessian_uses_b_pseudoinverse_in_working_coordinates(self):
        h_cart = np.diag(np.linspace(0.2, 1.0, 9))

        eng = NE.OQPEngine(
            self.WATER_AT, self.WATER_X, coordsys="dlc",
            initial_hessian=h_cart,
        )
        bmat = eng.coords.b_matrix(eng.x)
        g_inverse, _ = eng.coords._g_inverse(bmat)
        b_pinv = bmat.T @ g_inverse
        expected = b_pinv.T @ h_cart @ b_pinv

        self.assertTrue(np.allclose(eng.H, expected, atol=1e-12))
        self.assertTrue(np.allclose(eng.H, eng.H.T, atol=1e-14))

    def test_nonstationary_gradient_adds_internal_coordinate_curvature(self):
        h_cart = np.diag(np.linspace(0.2, 1.0, 9))
        gradient = np.linspace(-0.03, 0.04, 9)

        eng = NE.OQPEngine(
            self.WATER_AT, self.WATER_X, coordsys="dlc",
            initial_hessian=h_cart, initial_gradient=gradient,
        )
        bmat = eng.coords.b_matrix(eng.x)
        g_inverse, _ = eng.coords._g_inverse(bmat)
        b_pinv = bmat.T @ g_inverse
        g_work = eng.coords.grad_to_q(eng.x, gradient)
        correction = np.zeros_like(h_cart)
        step = 1.0e-4
        for column in range(9):
            xp = eng.x.copy(); xp[column] += step
            xm = eng.x.copy(); xm[column] -= step
            dbdx = (eng.coords.b_matrix(xp) - eng.coords.b_matrix(xm)) / (2 * step)
            correction[:, column] = g_work @ dbdx
        correction = 0.5 * (correction + correction.T)
        expected = b_pinv.T @ (h_cart - correction) @ b_pinv
        uncorrected = b_pinv.T @ h_cart @ b_pinv

        self.assertTrue(np.allclose(eng.H, expected, atol=1e-10))
        self.assertFalse(np.allclose(eng.H, uncorrected, atol=1e-8))

    def test_redundant_null_space_retains_model_curvature(self):
        # Four mutually bonded atoms generate more primitive internals than the
        # 3N-6 physical internal-coordinate rank.
        atoms = [6, 6, 6, 6]
        xyz = np.array([[0.0, 0.0, 0.0], [2.0, 0.1, 0.0],
                        [0.9, 1.7, 0.2], [0.8, 0.6, 1.8]])
        h_cart = np.eye(12) * 0.5

        eng = NE.OQPEngine(atoms, xyz, coordsys="ric",
                           initial_hessian=h_cart)
        bmat = eng.coords.b_matrix(eng.x)
        g_inverse, _ = eng.coords._g_inverse(bmat)
        b_pinv = bmat.T @ g_inverse
        null_projector = np.eye(bmat.shape[0]) - bmat @ b_pinv
        model = eng.coords.guess_hessian(eng.x)
        expected = (b_pinv.T @ h_cart @ b_pinv
                    + null_projector.T @ model @ null_projector)

        self.assertGreater(bmat.shape[0], np.linalg.matrix_rank(bmat))
        self.assertTrue(np.allclose(eng.H, expected, atol=1e-12))

    def test_injected_hessian_uses_the_optimizer_active_subspace_cutoff(self):
        class NearNullCoordinates:
            def b_matrix(self, _x):
                return np.diag([1.0, 1.0e-4, 1.0])

            def _g_inverse(self, bmat):
                # Mirrors the RIC relative 1e-6 eigenvalue cutoff: the middle
                # G eigenvalue is 1e-8 and is deliberately inactive.
                return np.diag([1.0, 0.0, 1.0]), 2

        eng = NE.OQPEngine([1], np.zeros(3), coordsys="cart")
        eng.coords = NearNullCoordinates()
        eng._global_external_mass_weighted_basis = lambda: np.empty((3, 0))
        model = np.eye(3) * 0.5

        transformed = eng._transform_initial_hessian(np.eye(3), model)

        self.assertTrue(np.allclose(transformed, np.diag([1.0, 0.5, 1.0])))

    def test_real_tric_ts_hessian_does_not_follow_external_zero_mode_noise(self):
        mass = np.array([15.999, 1.008, 1.008])
        xyz = np.array([[0.0, 0.0, 0.1], [0.0, 1.4, -0.8],
                        [0.0, -1.4, -0.8]])
        vibrational = NI._vibrational_basis(mass, xyz)
        external_projector = np.eye(9) - vibrational @ vibrational.T
        mw_hessian = (
            vibrational @ np.diag([-1.0e-6, 0.5, 0.8]) @ vibrational.T
            - 1.0e-5 * external_projector
        )
        mr = np.repeat(mass, 3)
        h_cart = mw_hessian * np.sqrt(np.outer(mr, mr))

        eng = NE.OQPEngine(
            [8, 1, 1], xyz, mode="ts", coordsys="tric",
            initial_hessian=h_cart, initial_gradient=np.zeros(9), masses=mass,
            project_global_rigid_modes=True,
        )
        values, vectors = np.linalg.eigh(eng.H)
        bmat = eng.coords.b_matrix(eng.x)
        g_inverse, _ = eng.coords._g_inverse(bmat)
        displacement = bmat.T @ (g_inverse @ vectors[:, 0])
        q_mw = displacement * np.sqrt(mr)
        q_mw /= np.linalg.norm(q_mw)

        self.assertEqual(np.count_nonzero(values < 0.0), 1)
        self.assertLess(float(q_mw @ external_projector @ q_mw), 1.0e-6)

    def test_external_regularization_does_not_erase_internal_ric_or_dlc_mode(self):
        mass = np.array([1.0, 35.0])
        xyz = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]])
        vibrational = NI._vibrational_basis(mass, xyz)
        mr = np.repeat(mass, 3)
        h_cart = (-0.5 * (vibrational @ vibrational.T)) * np.sqrt(
            np.outer(mr, mr)
        )

        for coordsys in ("ric", "dlc", "tric"):
            with self.subTest(coordsys=coordsys):
                eng = NE.OQPEngine(
                    [1, 17], xyz, mode="ts", coordsys=coordsys,
                    initial_hessian=h_cart, initial_gradient=np.zeros(6),
                    masses=mass, project_global_rigid_modes=True,
                )
                self.assertLess(float(np.linalg.eigvalsh(eng.H)[0]), -0.1)
        self.assertEqual(eng.coordsys, "TRIC->RIC(fallback)")

    def test_initial_hessian_validation(self):
        valid = np.eye(9)
        invalid = (
            (np.eye(8), "shape"),
            (np.full((9, 9), np.nan), "finite"),
            (valid + np.triu(np.ones((9, 9)), 1), "symmetric"),
        )

        for hessian, message in invalid:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    NE.OQPEngine(
                        self.WATER_AT, self.WATER_X, coordsys="cart",
                        initial_hessian=hessian,
                    )

        with self.assertRaisesRegex(ValueError, "initial_gradient must have shape"):
            NE.OQPEngine(
                self.WATER_AT, self.WATER_X, coordsys="dlc",
                initial_hessian=np.eye(9), initial_gradient=np.zeros(8),
            )

    def test_ts_step_follows_injected_negative_cartesian_mode(self):
        h_cart = np.diag([-0.5, 0.4, 0.8])
        eng = NE.OQPEngine(
            [1], [0.0, 0.0, 0.0], mode="ts", coordsys="cart",
            initial_hessian=h_cart,
        )

        step = eng._rfo_step(eng.H, np.array([0.1, 0.1, 0.1]))

        self.assertEqual(int(np.argmax(np.abs(eng._followed))), 0)
        self.assertNotEqual(step[0], 0.0)
        self.assertTrue(np.array_equal(eng.H, h_cart))

    def test_ts_follow_mode_has_a_clear_runtime_bound_error(self):
        eng = NE.OQPEngine(
            [1], [0.0, 0.0, 0.0], mode="ts", coordsys="cart",
            follow_mode=3,
        )

        with self.assertRaisesRegex(ValueError, "out of range for 3"):
            eng._rfo_step(eng.H, np.array([0.1, 0.1, 0.1]))


# --------------------------------------------------------------------------- #
class TestDispatchAndValidation(unittest.TestCase):
    _IMPORT_STATE = (
        "oqp", "oqp.utils", "oqp.utils.mpi_utils", "oqp.utils.input_checker",
        "oqp.utils.perf_levels",
    )

    def setUp(self):
        # The dispatcher tests load input_checker around lightweight module
        # stubs.  Preserve the real package imports so they cannot leak into
        # later integration tests (notably QM/MM, which imports perf_levels).
        self._saved_import_state = {
            name: sys.modules.get(name) for name in self._IMPORT_STATE
        }
        utils = types.ModuleType("oqp.utils")
        utils.__path__ = []
        sys.modules["oqp.utils"] = utils

    def tearDown(self):
        _restore_modules(self._saved_import_state)

    def test_input_checker_accepts_oqp(self):
        # input_checker only needs a tiny mpi_utils stub.
        utils = sys.modules.setdefault("oqp.utils",
                                       types.ModuleType("oqp.utils"))
        utils.__path__ = []
        mpi = types.ModuleType("oqp.utils.mpi_utils")
        mpi.MPIManager = type("MPIManager", (), {})
        sys.modules["oqp.utils.mpi_utils"] = mpi
        ic = _load("oqp.utils.input_checker",
                   ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py")
        self.assertIn("oqp", ic.OPT_LIBS)

        report = ic.CheckReport()
        cfg = {"input": {"runtype": "optimize", "method": "hf"},
               "optimize": {"lib": "oqp", "istate": 0}}
        ic._check_optimize(cfg, report)
        self.assertEqual(report.errors, [],
                         "oqp/optimize should validate clean")

        # Omit optimize.lib entirely: every workflow implemented by the native
        # backend must inherit oqp rather than requiring boilerplate input.
        for runtype in ("optimize", "ts", "meci", "mecp", "tci",
                        "neb", "irc", "mep"):
            native_default = ic.CheckReport()
            method = "tdhf" if runtype in {"meci", "mecp", "tci"} else "hf"
            config = {
                "input": {"runtype": runtype, "method": method},
                "optimize": {
                    "istate": 1 if method == "tdhf" else 0,
                    "jstate": 2,
                    "kstate": 3,
                    "imult": 1,
                    "jmult": 3,
                },
            }
            ic._check_optimize(config, native_default)
            lib_errors = [
                diagnostic for diagnostic in native_default.errors
                if diagnostic.path == "optimize.lib"
            ]
            self.assertEqual(
                lib_errors, [],
                f"omitted optimize.lib should select oqp for {runtype}")

        # mep/meci are supported on oqp; a non-optimization runtype is not.
        report2 = ic.CheckReport()
        cfg2 = {"input": {"runtype": "energy", "method": "hf"},
                "optimize": {"lib": "oqp", "istate": 0}}
        ic._check_optimize(cfg2, report2)
        msgs = " ".join(d.message for d in report2.diagnostics)
        self.assertIn("oqp optimizer currently supports", msgs)

        for rt, extra in (("mep", {}), ("meci", {"jstate": 2})):
            rep = ic.CheckReport()
            cfg3 = {"input": {"runtype": rt,
                              "method": "tdhf" if rt == "meci" else "hf"},
                    "optimize": {"lib": "oqp", "istate": 1, **extra}}
            ic._check_optimize(cfg3, rep)
            lib_errs = [d for d in rep.errors if d.path == "optimize.lib"]
            self.assertEqual(lib_errs, [], f"oqp/{rt} should be allowed")

    def test_neb_product_resolves_relative_to_input_dir(self):
        utils = sys.modules.setdefault("oqp.utils",
                                       types.ModuleType("oqp.utils"))
        utils.__path__ = []
        mpi = types.ModuleType("oqp.utils.mpi_utils")
        mpi.MPIManager = type("MPIManager", (), {})
        sys.modules["oqp.utils.mpi_utils"] = mpi
        ic = _load("oqp.utils.input_checker",
                   ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py")
        ex_dir = str(ROOT / "examples" / "OPT")
        cfg = {"input": {"runtype": "neb", "method": "hf"},
               "optimize": {"istate": 0, "lib": "oqp"},
               "neb": {"product": "HCN_RHF-DFT_NEB_OQP_product.xyz", "nimage": 7}}

        # Without input_dir the product (stored beside the .inp) is not found...
        no_dir = ic.CheckReport()
        ic._check_neb(cfg, no_dir)
        self.assertTrue(any(d.path == "neb.product" for d in no_dir.errors))

        # ...but it resolves when the input directory is supplied.
        with_dir = ic.CheckReport()
        ic._check_neb(cfg, with_dir, input_dir=ex_dir)
        self.assertEqual([d for d in with_dir.errors if d.path == "neb.product"], [])

    def test_tci_validation(self):
        utils = sys.modules.setdefault("oqp.utils",
                                       types.ModuleType("oqp.utils"))
        utils.__path__ = []
        mpi = types.ModuleType("oqp.utils.mpi_utils")
        mpi.MPIManager = type("MPIManager", (), {})
        sys.modules["oqp.utils.mpi_utils"] = mpi
        ic = _load("oqp.utils.input_checker",
                   ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py")
        self.assertIn("tci", ic.SUPPORTED_RUNTYPES)

        # valid: oqp, tdhf, strictly increasing states
        ok = ic.CheckReport()
        ic._check_optimize({"input": {"runtype": "tci", "method": "tdhf"},
                            "optimize": {"lib": "oqp", "istate": 1,
                                         "jstate": 2, "kstate": 3}}, ok)
        self.assertEqual([d for d in ok.errors if d.path in
                          ("optimize.lib", "optimize.kstate")], [])

        # invalid ordering
        bad = ic.CheckReport()
        ic._check_optimize({"input": {"runtype": "tci", "method": "tdhf"},
                            "optimize": {"lib": "oqp", "istate": 1,
                                         "jstate": 2, "kstate": 2}}, bad)
        self.assertTrue(any(d.path == "optimize.kstate" for d in bad.errors))

        # tci only on oqp
        geo = ic.CheckReport()
        ic._check_optimize({"input": {"runtype": "tci", "method": "tdhf"},
                            "optimize": {"lib": "geometric", "istate": 1,
                                         "jstate": 2, "kstate": 3}}, geo)
        self.assertTrue(any("only through the oqp" in d.message
                            for d in geo.errors))


class TestOQPNEB(unittest.TestCase):
    """Analytic CI-NEB on a 1-atom double well: V=(x^2-a^2)^2 + y^2 + z^2.

    Minima at (+/-a,0,0); saddle at (0,0,0) with E=a^4.  The climbing image
    must land on the saddle and the bowed initial path must straighten.
    """

    def test_climbing_image_finds_saddle(self):
        a = 1.5

        def eg(xf):
            x, y, z = xf
            E = (x * x - a * a) ** 2 + y * y + z * z
            g = np.array([2 * (x * x - a * a) * 2 * x, 2 * y, 2 * z])
            return E, g

        M = 9
        react = np.array([-a, 0, 0.0])
        prod = np.array([a, 0, 0.0])
        images = [react + (i / (M - 1)) * (prod - react)
                  + np.array([0, 0.6 * np.sin(np.pi * i / (M - 1)), 0])
                  for i in range(M)]
        neb = NEB(images, k_spring=0.5, climbing=True, climb_fmax=0.1)
        res = neb.run(eg, fmax_tol=1e-3, maxiter=800, dt=0.05, dt_max=0.2)
        self.assertTrue(res["converged"])
        top = int(np.argmax(res["energies"]))
        self.assertLess(abs(neb.images[top][0]), 0.1)          # saddle at x=0
        self.assertLess(abs(res["energies"][top] - a ** 4), 0.05)  # E = a^4
        # bowed-out y-displacement relaxed back toward zero
        ymax = max(abs(im.reshape(3)[1]) for im in neb.images)
        self.assertLess(ymax, 0.05)

    def test_maxmove_caps_step(self):
        # A huge force must not move an image more than maxmove in one step.
        def eg(xf):
            return float(xf @ xf), 100.0 * xf  # enormous gradient
        images = [np.array([0., 0, 0]), np.array([1., 0, 0]),
                  np.array([2., 0, 0])]
        neb = NEB(images, k_spring=0.1, climbing=False)
        x_before = neb.images[1].copy()
        neb.run(eg, fmax_tol=0.0, maxiter=1, dt=0.5, maxmove=0.2)
        moved = np.linalg.norm(neb.images[1] - x_before)
        self.assertLessEqual(moved, 0.2 + 1e-9)

    def test_climbing_convergence_uses_climbing_force(self):
        # Straight symmetric double-well path (g_perp == 0 everywhere) with an
        # EVEN number of images: no image starts on the saddle.  Convergence must
        # NOT be declared until the climbing image actually reaches the saddle --
        # judging it by g_perp alone would converge immediately off-saddle.
        a = 1.0

        def eg(xf):
            x, y, z = xf
            return (x * x - a * a) ** 2 + y * y + z * z, \
                np.array([4 * x * (x * x - a * a), 2 * y, 2 * z])

        M = 8  # even -> highest image is off the saddle (x = -+1/7)
        images = [np.array([-a + i * (2 * a) / (M - 1), 0.0, 0.0])
                  for i in range(M)]
        neb = NEB(images, k_spring=0.5, climbing=True, climb_fmax=1.0)
        res = neb.run(eg, fmax_tol=1e-3, maxiter=600, dt=0.05, dt_max=0.2)
        top = int(np.argmax(res["energies"]))
        self.assertTrue(res["converged"])
        # the climbing image must have climbed to the saddle (x = 0, E = a**4)
        self.assertLess(abs(neb.images[top][0]), 0.05)
        self.assertLess(abs(res["energies"][top] - a ** 4), 0.02)


class TestOQPIRC(unittest.TestCase):
    """Gonzalez-Schlegel IRC on V=(x^2-1)^2 + y^2 + z^2 (1 atom).

    Saddle at the origin (imaginary mode along x); IRC forward/backward must
    descend monotonically into the +/-x minima.
    """

    def _eg(self, xf):
        x, y, z = xf
        return (x * x - 1) ** 2 + y * y + z * z, \
            np.array([4 * x * (x * x - 1), 2 * y, 2 * z])

    def test_imaginary_mode(self):
        mode, curv = NI.imaginary_mode([1.0], np.diag([-4.0, 2.0, 2.0]))
        self.assertLess(curv, 0.0)
        self.assertAlmostEqual(abs(mode[0]), 1.0, places=6)

    def test_imaginary_mode_rejects_minimum_and_higher_order_saddle(self):
        with self.assertRaisesRegex(ValueError, "found 0"):
            NI.imaginary_mode([1.0], np.diag([1.0, 2.0, 3.0]))
        with self.assertRaisesRegex(ValueError, "found 2"):
            NI.imaginary_mode([1.0], np.diag([-4.0, -1.0, 2.0]))

    def test_imaginary_mode_projects_translation_and_rotation_noise(self):
        mass = np.array([16.0, 1.0, 1.0])
        xyz = np.array([[0.0, 0.0, 0.1], [0.0, 1.4, -0.8],
                        [0.0, -1.4, -0.8]])
        basis = NI._vibrational_basis(mass, xyz)
        self.assertEqual(basis.shape, (9, 3))
        vibrational_curvature = np.diag([-0.4, 0.7, 1.1])
        internal = basis @ vibrational_curvature @ basis.T
        external_projector = np.eye(9) - basis @ basis.T
        # Deliberately negative external-mode noise would look like seven
        # imaginary modes without projection.
        mw_hessian = internal - 1.0e-5 * external_projector
        mr = np.repeat(mass, 3)
        hessian = mw_hessian * np.sqrt(np.outer(mr, mr))

        mode, curvature = NI.imaginary_mode(mass, hessian, coordinates=xyz)

        self.assertAlmostEqual(curvature, -0.4, places=10)
        self.assertAlmostEqual(np.linalg.norm(mode), 1.0, places=12)

    def test_irc_descends_both_directions(self):
        mode, _ = NI.imaginary_mode([1.0], np.diag([-4.0, 2.0, 2.0]))
        for sign, target in ((+1, 1.0), (-1, -1.0)):
            res = IRC([1.0], [0, 0, 0.0], mode, step=0.08, sign=sign).run(
                self._eg, max_points=80, gtol=1e-4)
            self.assertTrue(res["converged"])
            self.assertLess(abs(res["points"][-1]["x"][0] - target), 0.1)
            es = [p["energy"] for p in res["points"]]
            self.assertTrue(all(es[i] >= es[i + 1] - 1e-9
                                for i in range(len(es) - 1)))
            if res["reason"] == "minimum_bracket":
                # The uphill fixed-step overshoot is discarded and the lower
                # accepted side of the bracket is restored.
                self.assertLessEqual(abs(res["points"][-1]["x"][0]), 1.0)

    def test_last_allowed_point_is_the_last_evaluated_geometry(self):
        evaluated = []

        def tracked_eg(xf):
            evaluated.append(np.asarray(xf, dtype=float).copy())
            return self._eg(xf)

        mode, _ = NI.imaginary_mode([1.0], np.diag([-4.0, 2.0, 2.0]))
        res = IRC([1.0], [0, 0, 0.0], mode, step=0.08).run(
            tracked_eg, max_points=1, gtol=1.0e-12,
        )

        self.assertFalse(res["converged"])
        self.assertEqual(len(evaluated), 1)
        self.assertTrue(np.array_equal(evaluated[-1], res["points"][-1]["x"]))

    def test_minimum_bracket_restores_callback_to_lower_point(self):
        evaluated = []

        def tracked_eg(xf):
            evaluated.append(np.asarray(xf, dtype=float).copy())
            return self._eg(xf)

        mode, _ = NI.imaginary_mode([1.0], np.diag([-4.0, 2.0, 2.0]))
        res = IRC([1.0], [0, 0, 0.0], mode, step=0.08).run(
            tracked_eg, max_points=80, gtol=1.0e-4,
        )

        self.assertEqual(res["reason"], "minimum_bracket")
        self.assertTrue(np.array_equal(evaluated[-1], res["points"][-1]["x"]))
        self.assertLessEqual(abs(res["points"][-1]["x"][0]), 1.0)


class TestOQPMEP(unittest.TestCase):
    """MEP = steepest descent from a non-stationary point to a minimum.

    Same V=(x^2-1)^2+y^2+z^2; starting on the slope, the path must descend
    monotonically into the nearby minimum.
    """

    def test_mep_descends_to_minimum(self):
        def eg(xf):
            x, y, z = xf
            return (x * x - 1) ** 2 + y * y + z * z, \
                np.array([4 * x * (x * x - 1), 2 * y, 2 * z])

        x0 = np.array([0.3, 0.5, 0.0])
        mass = np.array([1.0])
        sqm = np.sqrt(np.repeat(mass, 3))
        _, g0 = eg(x0)
        gmw = g0 / sqm
        direction = -gmw / np.linalg.norm(gmw)
        res = IRC(mass, x0, direction, step=0.08, sign=1).run(
            eg, max_points=120, gtol=1e-4)
        self.assertTrue(res["converged"])
        self.assertLess(abs(res["points"][-1]["x"][0] - 1.0), 0.1)
        es = [p["energy"] for p in res["points"]]
        self.assertTrue(all(es[i] >= es[i + 1] - 1e-9
                            for i in range(len(es) - 2)))


class TestOQPExample(unittest.TestCase):
    def test_example_exists_and_capped(self):
        path = EXAMPLES_OPT / "H2O_RHF-DFT_OPTIMIZE_OQP.inp"
        self.assertTrue(path.exists(), "oqp example input missing")
        text = path.read_text()
        self.assertIn("lib=oqp", text)
        # honor the repo example maxit cap
        for line in text.splitlines():
            if line.strip().startswith("maxit"):
                self.assertLessEqual(int(line.split("=")[1]), 10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
