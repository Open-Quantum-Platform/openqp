"""Tests for the MECP/MECI crossing objectives in ``libscipy``.

The objectives are pure NumPy, so they are exercised directly on a two-state
model whose crossing seam is known analytically:

    E_lo = 0.5 k (x^2 + y^2)
    E_up = E_lo + a (x - x0)          -> seam is the line x = x0

The exact MECP is (x0, 0) with a zero gap.  Following the effective gradient
downhill must reach it for the converging objectives, while the legacy
fixed-weight quadratic penalty settles at a residual gap of -0.495/(1 + 0.18 w),
which is derived by solving its stationarity condition for this model.

The compiled OQP backend is not required; the module's backend imports are
stubbed as in ``test_geometric_optimizer.py``.
"""
import importlib.util
import sys
import types
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

K, A, X0 = 1.0, 0.3, 1.5


def _load(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# Names this module injects while importing libscipy.  They are removed again
# immediately afterwards: pytest shares one interpreter across test modules, so
# leaving a stub (or an emptied __path__) behind would break every later import
# of the real oqp.utils.* and oqp.library.* submodules.
_STUBBED = (
    "oqp",
    "oqp.library",
    "oqp.library.single_point",
    "oqp.library.baeka",
    "oqp.utils",
    "oqp.utils.file_utils",
    "oqp.utils.state_labels",
    "oqp.utils.qmmm",
    "oqp.utils.mpi_utils",
)


def _stub(name, is_package=False):
    module = types.ModuleType(name)
    if is_package:
        module.__path__ = []
    sys.modules[name] = module
    return module


def _load_isolated(path, extra_module=None):
    """Import a pyoqp module against stubs, then restore sys.modules."""
    saved = {name: sys.modules.get(name) for name in _STUBBED}

    _stub("oqp")
    _stub("oqp.library", is_package=True)
    _stub("oqp.utils", is_package=True)

    single_point = _stub("oqp.library.single_point")
    for name in ("SinglePoint", "Gradient", "LastStep"):
        setattr(single_point, name, type(name, (), {}))

    file_utils = _stub("oqp.utils.file_utils")
    file_utils.dump_log = lambda *args, **kwargs: None
    file_utils.dump_data = lambda *args, **kwargs: None

    state_labels = _stub("oqp.utils.state_labels")
    state_labels.is_mrsf = lambda config: False
    state_labels.public_state_label = lambda *args, **kwargs: ""

    _stub("oqp.utils.qmmm")
    mpi_utils = _stub("oqp.utils.mpi_utils")
    mpi_utils.MPIManager = type("MPIManager", (), {})

    try:
        # baeka is pure NumPy, so the real kernel is used
        baeka = _load("oqp.library.baeka", "pyoqp/oqp/library/baeka.py")
        loaded = _load("oqp.library.libscipy", path)
        checker = (_load(extra_module[0], extra_module[1])
                   if extra_module else None)
    finally:
        for name, module in saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module
        sys.modules.pop("oqp.library.libscipy", None)
        sys.modules.pop("oqp.utils.input_checker", None)

    return loaded, baeka, checker


libscipy, baeka_module, input_checker = _load_isolated(
    "pyoqp/oqp/library/libscipy.py",
    extra_module=("oqp.utils.input_checker", "pyoqp/oqp/utils/input_checker.py"),
)
BaekAState = baeka_module.BaekAState


def model(coordinates):
    """Return (E_lo, G_lo, E_up, G_up) of the two-state model."""
    x, y = coordinates
    e_lo = 0.5 * K * (x * x + y * y)
    g_lo = np.array([K * x, K * y])
    e_up = e_lo + A * (x - X0)
    g_up = g_lo + np.array([A, 0.0])
    return e_lo, g_lo, e_up, g_up


def exact_gap(coordinates):
    return A * (coordinates[0] - X0)


def make_mecp(search, **attributes):
    """Build a MECPOpt bound to the model without touching the backend."""
    opt = object.__new__(libscipy.MECPOpt)
    opt.istate, opt.jstate, opt.nstate = 1, 1, 2
    opt.mecp_search = search
    opt.weights = 1.0
    opt.sigma = 1.0
    opt.incre = 1.0
    opt.alpha = 0.02
    opt.mu = 10.0
    opt.lagrange = 0.0
    opt.mu_scaled = None
    opt.metrics = {}
    opt.energy_gap = 1.0e-4
    opt.energy_shift = 1.0e-6
    opt.rmsd_grad = 1.0e-4
    for key, value in attributes.items():
        setattr(opt, key, value)
    # keep the objectives free of file I/O and optimizer bookkeeping
    def _record(coordinates, f, df, gap_e, gap_term):
        opt.metrics["gap"] = gap_e
        return f, df

    opt.record = _record
    return opt


def descend(objective, coordinates, steps=4000, step=0.05):
    """Follow the effective gradient downhill and return the final point."""
    coordinates = np.array(coordinates, dtype=float)
    for _ in range(steps):
        e_lo, g_lo, e_up, g_up = model(coordinates)
        # index layout expected by MECPOpt: state j sits after nstate entries
        energies = np.array([0.0, e_lo, 0.0, e_up])
        grads = np.array([np.zeros(2), g_lo, np.zeros(2), g_up])
        _, df = objective(coordinates, energies, grads)
        coordinates = coordinates - step * df
    return coordinates


class MECPObjectiveTest(unittest.TestCase):
    def test_auglag_reaches_the_exact_crossing(self):
        opt = make_mecp("auglag")
        final = descend(opt.auglag, [0.0, 0.6])
        self.assertLess(abs(exact_gap(final)), 1.0e-10)
        self.assertAlmostEqual(final[0], X0, places=8)
        self.assertAlmostEqual(final[1], 0.0, places=8)

    def test_auglag_reduces_to_the_bearpark_projection(self):
        """At pen_sigma = 1 the multiplier form is the published projection."""
        self.assertNotIn("ubp", libscipy.MECPOpt.__dict__)
        opt = make_mecp("auglag")
        coordinates = np.array([0.3, 0.4])
        e_lo, g_lo, e_up, g_up = model(coordinates)
        energies = np.array([0.0, e_lo, 0.0, e_up])
        grads = np.array([np.zeros(2), g_lo, np.zeros(2), g_up])
        # gap_sigma = 1 must reproduce the Bearpark projection 2 dE x once the
        # multiplier has been seeded from this geometry
        opt.mu = 1.0
        opt.auglag(coordinates, energies, grads)
        opt.mu_scaled = 2.0 * opt.mu / np.linalg.norm(g_up - g_lo)
        _, df = opt.auglag(coordinates, energies, grads)
        gap_g = g_up - g_lo
        x_hat = gap_g / np.linalg.norm(gap_g)
        mean_g = 0.5 * (g_lo + g_up)
        reference = 2.0 * (e_up - e_lo) * x_hat + (
            mean_g - np.dot(mean_g, x_hat) * x_hat
        )
        np.testing.assert_allclose(df, reference, atol=1.0e-12)

    def test_auglag_gradient_matches_a_finite_difference_of_its_value(self):
        """The multipliers are constants within one evaluation.

        If lam and mu were recomputed from the geometry while f is being
        differentiated, df would not be the gradient of the reported f and the
        optimizer's line searches and Hessian updates would see an
        inconsistent pair.
        """
        def value_and_grad(coordinates, lagrange, mu_scaled):
            opt = make_mecp("auglag", lagrange=lagrange, mu_scaled=mu_scaled)
            e_lo, g_lo, e_up, g_up = model(coordinates)
            energies = np.array([0.0, e_lo, 0.0, e_up])
            grads = np.array([np.zeros(2), g_lo, np.zeros(2), g_up])
            return opt.auglag(coordinates, energies, grads)

        # freeze the multipliers at the values a previous evaluation would leave
        seed = make_mecp("auglag")
        centre = np.array([0.3, 0.4])
        e_lo, g_lo, e_up, g_up = model(centre)
        seed.auglag(centre, np.array([0.0, e_lo, 0.0, e_up]),
                    np.array([np.zeros(2), g_lo, np.zeros(2), g_up]))
        lagrange, mu_scaled = seed.lagrange, seed.mu_scaled

        _, df = value_and_grad(centre, lagrange, mu_scaled)
        step = 1.0e-6
        for axis in (0, 1):
            plus, minus = centre.copy(), centre.copy()
            plus[axis] += step
            minus[axis] -= step
            f_plus, _ = value_and_grad(plus, lagrange, mu_scaled)
            f_minus, _ = value_and_grad(minus, lagrange, mu_scaled)
            numerical = (f_plus - f_minus) / (2.0 * step)
            self.assertAlmostEqual(df[axis], numerical, places=6)

    def test_quad_settles_at_a_residual_gap(self):
        """The legacy objective cannot meet a tight energy_gap criterion."""
        for weight in (1.0, 10.0, 100.0):
            opt = make_mecp("quad", weights=weight)
            final = descend(opt.quad, [0.0, 0.6])
            expected = -0.495 / (1.0 + 0.18 * weight)
            self.assertAlmostEqual(exact_gap(final), expected, places=8)
            self.assertGreater(abs(exact_gap(final)), 1.0e-4)

    def test_search_names_are_registered(self):
        checker = input_checker
        self.assertIn("auglag", checker.MECP_SEARCH)
        self.assertIn("quad", checker.MECP_SEARCH)
        self.assertNotIn("ubp", checker.MECP_SEARCH)
        self.assertNotIn("baeka", checker.MECP_SEARCH)
        self.assertNotIn("nonsense", checker.MECP_SEARCH)
        self.assertIn("auglag", checker.MECI_SEARCH)


class MECPSQPTest(unittest.TestCase):
    """The SQP driver on the same model, with a curved seam.

    A straight seam makes the constraint linear and SQP exact in one step,
    which would flatter it, so the gap is given a quadratic term in y.
    """

    CURVE = 0.05  # must stay below A/(2 X0) or (X0, 0) is not the seam minimum

    def curved(self, coordinates):
        x, y = coordinates
        e_lo = 0.5 * K * (x * x + y * y)
        g_lo = np.array([K * x, K * y])
        c = A * (x - X0) + self.CURVE * y * y
        g_c = np.array([A, 2.0 * self.CURVE * y])
        return e_lo, g_lo, e_lo + c, g_lo + g_c

    def curved_gap(self, coordinates):
        return A * (coordinates[0] - X0) + self.CURVE * coordinates[1] ** 2

    def make_driver(self):
        opt = object.__new__(libscipy.SQPMECPOpt)
        opt.istate, opt.jstate, opt.nstate = 1, 1, 2
        opt.mecp_search = "sqp"
        opt.lagrange = 0.0
        opt.metrics = {}
        opt.trust, opt.trust_max, opt.trust_min = 0.5, 1.0, 1.0e-6
        return opt

    def test_kkt_step_solves_the_multiplier_and_step_together(self):
        """The multiplier is a result of the step equations, not a formula."""
        opt = self.make_driver()
        hess = np.eye(2)
        mean_g = np.array([0.4, -0.2])
        gap_g = np.array([A, 0.1])
        c = 0.03
        step, lagrange = opt.kkt_step(hess, mean_g, gap_g, c)
        # the KKT rows must hold exactly
        np.testing.assert_allclose(hess @ step + mean_g + lagrange * gap_g,
                                   np.zeros(2), atol=1.0e-12)
        self.assertAlmostEqual(gap_g @ step + c, 0.0, places=12)

    def test_singular_system_is_reported_not_guessed(self):
        opt = self.make_driver()
        # identical state gradients leave the branching direction undefined
        self.assertIsNone(
            opt.kkt_step(np.zeros((2, 2)), np.array([1.0, 0.0]),
                         np.zeros(2), 0.1)
        )

    def test_driver_reaches_the_exact_crossing_on_a_curved_seam(self):
        """Drive the KKT loop directly, without the electronic-structure layer."""
        opt = self.make_driver()
        coordinates = np.array([0.0, 0.6])
        hess = np.eye(2)
        previous = None
        for _ in range(60):
            e_lo, g_lo, e_up, g_up = self.curved(coordinates)
            c, gap_g = e_up - e_lo, g_up - g_lo
            mean_g = 0.5 * (g_lo + g_up)
            if previous is not None:
                step_old, grad_old = previous
                dgrad = mean_g + opt.lagrange * gap_g - grad_old
                curvature = step_old @ hess @ step_old
                if curvature > 0 and step_old @ dgrad > 1.0e-14:
                    hess = (hess + np.outer(dgrad, dgrad) / (step_old @ dgrad)
                            - np.outer(hess @ step_old, hess @ step_old) / curvature)
            step, lagrange = opt.kkt_step(hess, mean_g, gap_g, c)
            opt.lagrange = lagrange
            length = np.linalg.norm(step)
            if length > opt.trust:
                step = step * (opt.trust / length)
            previous = (step, mean_g + lagrange * gap_g)
            coordinates = coordinates + step
        self.assertLess(abs(self.curved_gap(coordinates)), 1.0e-9)
        self.assertAlmostEqual(coordinates[0], X0, places=6)
        self.assertAlmostEqual(coordinates[1], 0.0, places=6)

    def test_sqp_needs_no_gap_sigma(self):
        """The driver never reads the auglag penalty scale.

        Checked on the parsed tree rather than the text, so the prose in the
        docstring explaining that gap_sigma is not needed does not count.
        """
        import ast

        source = (ROOT / "pyoqp/oqp/library/libscipy.py").read_text()
        tree = ast.parse(source)
        driver = next(node for node in tree.body
                      if isinstance(node, ast.ClassDef)
                      and node.name == "SQPMECPOpt")
        for node in ast.walk(driver):
            if isinstance(node, ast.Constant) and node.value == "gap_sigma":
                self.fail("SQPMECPOpt reads gap_sigma")
            if isinstance(node, ast.Attribute) and node.attr == "mu":
                self.fail("SQPMECPOpt reads the auglag penalty parameter")


class MECIObjectiveTest(unittest.TestCase):
    def make_meci(self):
        opt = object.__new__(libscipy.MECIOpt)
        opt.istate, opt.jstate = 0, 1
        opt.meci_search = "auglag"
        opt.weights = 1.0
        opt.mu = 10.0
        opt.lagrange = 0.0
        opt.mu_scaled = None
        opt.x = np.zeros(0)
        opt.y = np.zeros(0)
        opt.metrics = {}
        opt.itr = 0
        opt.pre_energy = 0.0
        opt.pre_coord = np.zeros(2)
        opt.atoms = np.zeros((1, 1))
        opt.mol = types.SimpleNamespace(log_path=".")
        return opt

    def descend(self, objective, opt, steps=4000, step=0.05):
        coordinates = np.array([0.0, 0.6])
        for _ in range(steps):
            e_lo, g_lo, e_up, g_up = model(coordinates)
            _, df = objective(coordinates, np.array([e_lo, e_up]),
                              np.array([g_lo, g_up]))
            coordinates = coordinates - step * df
        return coordinates

    def test_penalty_settles_at_a_residual_gap(self):
        """Why auto no longer resolves to the plain penalty.

        Its stationary point sits at a finite gap, so on its own it cannot
        reach energy_gap; only the BaekA escalation used to rescue it.
        """
        opt = self.make_meci()
        opt.sigma, opt.alpha, opt.incre = 1.0, 0.0, 1.0
        final = self.descend(opt.penalty, opt, steps=2000)
        self.assertGreater(abs(exact_gap(final)), 1.0e-2)

    def test_auglag_reaches_the_exact_crossing(self):
        opt = self.make_meci()
        coordinates = np.array([0.0, 0.6])
        for _ in range(4000):
            e_lo, g_lo, e_up, g_up = model(coordinates)
            energies = np.array([e_lo, e_up])
            grads = np.array([g_lo, g_up])
            _, df = opt.auglag(coordinates, energies, grads)
            coordinates = coordinates - 0.05 * df
        self.assertLess(abs(exact_gap(coordinates)), 1.0e-8)
        self.assertAlmostEqual(coordinates[0], X0, places=6)


if __name__ == "__main__":
    unittest.main()
