"""Tests for the builtin geometry optimizer (coordinates, engine, wiring).

The numerical core (builtin_coords, builtin_engine) is pure NumPy/SciPy and is
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


def _load_builtin_modules():
    """Load builtin_coords and builtin_engine without the oqp backend."""
    for pkg in ("oqp", "oqp.library"):
        m = sys.modules.get(pkg) or types.ModuleType(pkg)
        m.__path__ = []
        sys.modules[pkg] = m
    nc = _load("oqp.library.builtin_coords", LIB / "builtin_coords.py")
    ne = _load("oqp.library.builtin_engine", LIB / "builtin_engine.py")
    nb = _load("oqp.library.builtin_neb", LIB / "builtin_neb.py")
    ni = _load("oqp.library.builtin_irc", LIB / "builtin_irc.py")
    return nc, ne, nb, ni


NC, NE, NB, NI = _load_builtin_modules()
NEB = NB.NEB
IRC = NI.IRC


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


class TestBuiltinEngineConverges(unittest.TestCase):
    def _optimize(self, atoms, x0, pairs, gtol=1e-4, maxiter=100):
        eg = _morse(pairs)
        eng = NE.BuiltinEngine(atoms, x0, mode="min", maxiter=maxiter)

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
    def test_input_checker_accepts_builtin(self):
        # input_checker only needs a tiny mpi_utils stub.
        utils = sys.modules.setdefault("oqp.utils",
                                       types.ModuleType("oqp.utils"))
        utils.__path__ = []
        mpi = types.ModuleType("oqp.utils.mpi_utils")
        mpi.MPIManager = type("MPIManager", (), {})
        sys.modules["oqp.utils.mpi_utils"] = mpi
        ic = _load("oqp.utils.input_checker",
                   ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py")
        self.assertIn("builtin", ic.OPT_LIBS)

        report = ic.CheckReport()
        cfg = {"input": {"runtype": "optimize", "method": "hf"},
               "optimize": {"lib": "builtin", "istate": 0}}
        ic._check_optimize(cfg, report)
        self.assertEqual(report.errors, [],
                         "builtin/optimize should validate clean")

        # mep/meci are supported on builtin; a non-optimization runtype is not.
        report2 = ic.CheckReport()
        cfg2 = {"input": {"runtype": "energy", "method": "hf"},
                "optimize": {"lib": "builtin", "istate": 0}}
        ic._check_optimize(cfg2, report2)
        msgs = " ".join(d.message for d in report2.diagnostics)
        self.assertIn("builtin optimizer currently supports", msgs)

        for rt, extra in (("mep", {}), ("meci", {"jstate": 2})):
            rep = ic.CheckReport()
            cfg3 = {"input": {"runtype": rt,
                              "method": "tdhf" if rt == "meci" else "hf"},
                    "optimize": {"lib": "builtin", "istate": 1, **extra}}
            ic._check_optimize(cfg3, rep)
            lib_errs = [d for d in rep.errors if d.path == "optimize.lib"]
            self.assertEqual(lib_errs, [], f"builtin/{rt} should be allowed")

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

        # valid: builtin, tdhf, strictly increasing states
        ok = ic.CheckReport()
        ic._check_optimize({"input": {"runtype": "tci", "method": "tdhf"},
                            "optimize": {"lib": "builtin", "istate": 1,
                                         "jstate": 2, "kstate": 3}}, ok)
        self.assertEqual([d for d in ok.errors if d.path in
                          ("optimize.lib", "optimize.kstate")], [])

        # invalid ordering
        bad = ic.CheckReport()
        ic._check_optimize({"input": {"runtype": "tci", "method": "tdhf"},
                            "optimize": {"lib": "builtin", "istate": 1,
                                         "jstate": 2, "kstate": 2}}, bad)
        self.assertTrue(any(d.path == "optimize.kstate" for d in bad.errors))

        # tci only on builtin
        geo = ic.CheckReport()
        ic._check_optimize({"input": {"runtype": "tci", "method": "tdhf"},
                            "optimize": {"lib": "geometric", "istate": 1,
                                         "jstate": 2, "kstate": 3}}, geo)
        self.assertTrue(any("only through the builtin" in d.message
                            for d in geo.errors))


class TestBuiltinNEB(unittest.TestCase):
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


class TestBuiltinIRC(unittest.TestCase):
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

    def test_irc_descends_both_directions(self):
        mode, _ = NI.imaginary_mode([1.0], np.diag([-4.0, 2.0, 2.0]))
        for sign, target in ((+1, 1.0), (-1, -1.0)):
            res = IRC([1.0], [0, 0, 0.0], mode, step=0.08, sign=sign).run(
                self._eg, max_points=80, gtol=1e-4)
            self.assertTrue(res["converged"])
            self.assertLess(abs(res["points"][-1]["x"][0] - target), 0.1)
            es = [p["energy"] for p in res["points"]]
            # monotonic descent up to the final basin-overshoot point
            self.assertTrue(all(es[i] >= es[i + 1] - 1e-9
                                for i in range(len(es) - 2)))


class TestBuiltinMEP(unittest.TestCase):
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


class TestBuiltinExample(unittest.TestCase):
    def test_example_exists_and_capped(self):
        path = EXAMPLES_OPT / "H2O_RHF-DFT_OPTIMIZE_BUILTIN.inp"
        self.assertTrue(path.exists(), "builtin example input missing")
        text = path.read_text()
        self.assertIn("lib=builtin", text)
        # honor the repo example maxit cap
        for line in text.splitlines():
            if line.strip().startswith("maxit"):
                self.assertLessEqual(int(line.split("=")[1]), 10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
