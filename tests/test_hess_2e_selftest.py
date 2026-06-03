"""Integration test for the analytic two-electron (ERI) Hessian skeleton.

Drives the ``grd2_hess_selftest`` bind(C) harness, which finite-differences the
production analytic 2e gradient (grd2::grd2_driver) and compares it to the
analytic per-quartet 2e second-derivative driver (grd2::grd2_hess_driver),
both contracted with the same fixed converged RHF density. The identity
d/dR [ sum P P dX/dR ] = sum P P d2X/dR2 holds for a fixed density P, so the
test isolates the ERI second-derivative integral contraction.

Skipped unless the compiled OpenQP runtime and its Python dependencies are
importable (i.e. a built tree with OPENQP_ROOT set).
"""

import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELFTEST_OUT = Path("/tmp/grd2_hess_selftest.out")

INPUT = """[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g*
method=hf
[guess]
type=huckel
[scf]
multiplicity=1
type=rhf
"""


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class Hess2eSelfTest(unittest.TestCase):
    def test_two_electron_hessian_matches_finite_difference(self):
        import oqp
        from oqp.pyoqp import Runner

        workdir = Path("/tmp/oqp_hess2e_test")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "h2o.inp"
        inp.write_text(INPUT)
        log = workdir / "h2o.log"

        if SELFTEST_OUT.exists():
            SELFTEST_OUT.unlink()

        runner = Runner(project="h2o_hess2e", input_file=str(inp), log=str(log))
        runner.run()
        oqp.grd2_hess_selftest(runner.mol)

        self.assertTrue(SELFTEST_OUT.exists(), "self-test produced no output file")
        result = SELFTEST_OUT.read_text()
        self.assertIn("GRD2_HESS_SELFTEST PASS", result,
                      "analytic 2e Hessian disagrees with finite difference:\n" + result)


if __name__ == "__main__":
    unittest.main()
