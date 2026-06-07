"""Tests for the native geometry optimizer (coordinates, engine, wiring).

The numerical core (native_coords, native_engine) is pure NumPy/SciPy and is
tested directly.  The dispatcher/validator wiring is checked with lightweight
stubs so the compiled OQP backend is not required, mirroring
``test_geometric_optimizer.py``.
"""
import importlib.util
import sys
import types
import unittest
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


def _load_native_modules():
    """Load native_coords and native_engine without the oqp backend."""
    for pkg in ("oqp", "oqp.library"):
        m = sys.modules.get(pkg) or types.ModuleType(pkg)
        m.__path__ = []
        sys.modules[pkg] = m
    nc = _load("oqp.library.native_coords", LIB / "native_coords.py")
    ne = _load("oqp.library.native_engine", LIB / "native_engine.py")
    nb = _load("oqp.library.native_neb", LIB / "native_neb.py")
    return nc, ne, nb


NC, NE, NB = _load_native_modules()
NEB = NB.NEB


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
            self.assertIsInstance(ic, NC.RedundantInternalCoordinates)
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
        co2 = NC.build_coordinates([6, 8, 8],
                                   np.array([[0, 0, 0], [0, 0, 2.2],
                                             [0, 0, -2.2]], float))
        self.assertIsInstance(co2, NC.CartesianCoordinates)


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


class TestNativeEngineConverges(unittest.TestCase):
    def _optimize(self, atoms, x0, pairs, gtol=1e-4, maxiter=100):
        eg = _morse(pairs)
        eng = NE.NativeEngine(atoms, x0, mode="min", maxiter=maxiter)

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


# --------------------------------------------------------------------------- #
class TestDispatchAndValidation(unittest.TestCase):
    def test_input_checker_accepts_native(self):
        # input_checker only needs a tiny mpi_utils stub.
        utils = sys.modules.setdefault("oqp.utils",
                                       types.ModuleType("oqp.utils"))
        utils.__path__ = []
        mpi = types.ModuleType("oqp.utils.mpi_utils")
        mpi.MPIManager = type("MPIManager", (), {})
        sys.modules["oqp.utils.mpi_utils"] = mpi
        ic = _load("oqp.utils.input_checker",
                   ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py")
        self.assertIn("native", ic.OPT_LIBS)

        report = ic.CheckReport()
        cfg = {"input": {"runtype": "optimize", "method": "hf"},
               "optimize": {"lib": "native", "istate": 0}}
        ic._check_optimize(cfg, report)
        self.assertEqual(report.errors, [],
                         "native/optimize should validate clean")

        # meci is supported on native; irc is not.
        report2 = ic.CheckReport()
        cfg2 = {"input": {"runtype": "irc", "method": "hf"},
                "optimize": {"lib": "native", "istate": 0}}
        ic._check_optimize(cfg2, report2)
        msgs = " ".join(d.message for d in report2.diagnostics)
        self.assertIn("native optimizer currently supports", msgs)

        report3 = ic.CheckReport()
        cfg3 = {"input": {"runtype": "meci", "method": "tdhf"},
                "optimize": {"lib": "native", "istate": 1, "jstate": 2}}
        ic._check_optimize(cfg3, report3)
        lib_errs = [d for d in report3.errors if d.path == "optimize.lib"]
        self.assertEqual(lib_errs, [], "native/meci should be allowed")

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

        # valid: native, tdhf, strictly increasing states
        ok = ic.CheckReport()
        ic._check_optimize({"input": {"runtype": "tci", "method": "tdhf"},
                            "optimize": {"lib": "native", "istate": 1,
                                         "jstate": 2, "kstate": 3}}, ok)
        self.assertEqual([d for d in ok.errors if d.path in
                          ("optimize.lib", "optimize.kstate")], [])

        # invalid ordering
        bad = ic.CheckReport()
        ic._check_optimize({"input": {"runtype": "tci", "method": "tdhf"},
                            "optimize": {"lib": "native", "istate": 1,
                                         "jstate": 2, "kstate": 2}}, bad)
        self.assertTrue(any(d.path == "optimize.kstate" for d in bad.errors))

        # tci only on native
        geo = ic.CheckReport()
        ic._check_optimize({"input": {"runtype": "tci", "method": "tdhf"},
                            "optimize": {"lib": "geometric", "istate": 1,
                                         "jstate": 2, "kstate": 3}}, geo)
        self.assertTrue(any("only through the native" in d.message
                            for d in geo.errors))


class TestNativeNEB(unittest.TestCase):
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


class TestNativeExample(unittest.TestCase):
    def test_example_exists_and_capped(self):
        path = EXAMPLES_OPT / "H2O_RHF-DFT_OPTIMIZE_NATIVE.inp"
        self.assertTrue(path.exists(), "native example input missing")
        text = path.read_text()
        self.assertIn("lib=native", text)
        # honor the repo example maxit cap
        for line in text.splitlines():
            if line.strip().startswith("maxit"):
                self.assertLessEqual(int(line.split("=")[1]), 10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
