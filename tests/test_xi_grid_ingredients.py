"""Grid ingredient export (n, sigma, tau, xi^alpha for several alpha).

Checks, after a converged H2O M06-2X SCF, that the exported xi^1 equals tau,
xi^0 (p=0) equals n/2, the weights integrate the density to the electron
count, and a mid alpha is positive and distinct.  Skipped without a compiled
OpenQP runtime.
"""

import os
import unittest
from pathlib import Path

INPUT = """[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g
method=hf
functional=m06-2x
d4=False
[guess]
type=huckel
[scf]
multiplicity=1
type=rhf
conv=1.0e-8
[dftgrid]
rad_npts=48
ang_npts=110
pruned=
"""


def _runtime_available():
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return hasattr(oqp, "oqp_xi_grid_ingredients")
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class XiGridIngredientsTest(unittest.TestCase):
    def test_export_matches_tau_and_density(self):
        import numpy as np
        from oqp.pyoqp import Runner
        from oqp.library.xi_ingredients import grid_ingredients

        workdir = Path("/tmp/oqp_xi_grid_test")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "h2o.inp"
        inp.write_text(INPUT)
        runner = Runner(project="h2o_xi_grid", input_file=str(inp), log=str(workdir / "h2o.log"))
        runner.run()
        g = grid_ingredients(runner.mol, alphas=[1.0, 0.5, 0.0, -0.5], p=[1, 1, 0, 0])
        npts = g["xyzw"].shape[0]
        self.assertGreater(npts, 1000)
        for key, ncol in (("rho", 2), ("sigma", 3), ("tau", 2), ("xi", 8)):
            self.assertEqual(g[key].shape, (npts, ncol), key)
        w = g["xyzw"][:, 3]
        nelec = float(np.sum(w * (g["rho"][:, 0] + g["rho"][:, 1])))
        self.assertAlmostEqual(nelec, 10.0, places=4)
        # xi^1 == tau, xi^0- == n/2
        self.assertLess(np.max(np.abs(g["xi"][:, 0:2] - g["tau"])), 1e-10)
        self.assertLess(np.max(np.abs(g["xi"][:, 4:6] - 0.5 * g["rho"])), 1e-10)
        # fractional orders: non-negative, finite, different from tau
        self.assertTrue(np.all(np.isfinite(g["xi"])))
        self.assertGreaterEqual(float(np.min(g["xi"])), -1e-14)
        self.assertGreater(np.max(np.abs(g["xi"][:, 2] - g["tau"][:, 0])), 1e-3)
        self.assertEqual(list(g["p"]), [1, 1, 0, 0])


if __name__ == "__main__":
    unittest.main()
