"""Integration test for the native 1e (overlap + kinetic) Hessian terms.

Drives the ``hess1_selftest`` bind(C) harness, which finite-differences the
production gradient routines (grd1::grad_ee_overlap / grad_ee_kinetic) and
compares them to the analytic second-derivative drivers
(grd1::hess_ee_overlap / hess_ee_kinetic) on a real basis.

The test is skipped unless the compiled OpenQP runtime and its Python
dependencies are importable (i.e. a built tree with OPENQP_ROOT set), so it is
a no-op in environments that only run the dependency-light unit tests.
"""

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]

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
    if not (ROOT / "lib").glob("liboqp.*"):
        pass
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class Hess1eSelfTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="oqp_hess1e_")
        self.addCleanup(tmp.cleanup)
        self.tmpdir = Path(tmp.name)
        # The native self-test writes its report to $OQP_SELFTEST_DIR
        # (source/io/selftest_output.F90); keep it private to this test.
        env = mock.patch.dict(os.environ, {"OQP_SELFTEST_DIR": str(self.tmpdir)})
        env.start()
        self.addCleanup(env.stop)

    def test_overlap_and_kinetic_hessian_match_finite_difference(self):
        import oqp
        from oqp.pyoqp import Runner

        workdir = self.tmpdir
        inp = workdir / "h2o.inp"
        inp.write_text(INPUT)
        log = workdir / "h2o.log"

        runner = Runner(project="h2o_hess1e", input_file=str(inp), log=str(log))
        runner.run()
        oqp.hess1_selftest(runner.mol)
        report = self.tmpdir / "hess1_selftest.out"

        self.assertTrue(report.exists(), "self-test produced no output file")
        result = report.read_text()
        self.assertIn("HESS1E_SELFTEST PASS", result, msg=result)


if __name__ == "__main__":
    unittest.main()
