"""Integration test for the open-shell (UHF/ROHF) analytic Hessian skeleton.

Drives the ``hess_skel_open_selftest`` bind(C) harness, the open-shell analog
of ``hess_skel_selftest``. It builds the analytic fixed-density Hessian skeleton
(1e + 2e, excluding nuclear repulsion) from the open-shell density
(OQP_DM_A + OQP_DM_B), the open-shell energy-weighted density returned by
grd1::eijden, and the open-shell two-electron contraction
(grd2_uhf_compute_data_t), then compares it to a central finite difference of
the open-shell frozen-density gradient. The identity
d/dR [ sum M dX/dR ] = sum M d2X/dR2 holds for fixed W and spin densities, so a
PASS validates the open-shell skeleton assembly (signs, W, total-vs-spin
density routing) for both UHF and ROHF -- the prerequisite for the open-shell
CPHF response term.

The test is skipped unless the compiled OpenQP runtime is importable.
"""

import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELFTEST_OUT = Path("/tmp/hess_skel_open_selftest.out")

# OH radical (doublet), small enough to finite-difference quickly while still
# exercising d functions (6-31g* on oxygen) and a genuine open shell.
INPUT_TMPL = """[input]
system=
   8   0.000000000   0.000000000   0.000000000
   1   0.000000000   0.000000000   0.969000000
charge=0
runtype=energy
basis=6-31g*
method=hf
[guess]
type=huckel
[scf]
multiplicity=2
type={scftype}
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
class HessSkelOpenSelfTest(unittest.TestCase):
    def _run_case(self, scftype):
        import oqp
        from oqp.pyoqp import Runner

        workdir = Path(f"/tmp/oqp_hess_skel_open_{scftype}")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "oh.inp"
        inp.write_text(INPUT_TMPL.format(scftype=scftype))
        log = workdir / "oh.log"

        if SELFTEST_OUT.exists():
            SELFTEST_OUT.unlink()

        runner = Runner(project=f"oh_{scftype}", input_file=str(inp), log=str(log))
        runner.run()
        oqp.hess_skel_open_selftest(runner.mol)

        self.assertTrue(SELFTEST_OUT.exists(), "self-test produced no output file")
        result = SELFTEST_OUT.read_text()
        self.assertIn("HESS_SKEL_OPEN_SELFTEST PASS", result, msg=result)

    def test_uhf_skeleton_matches_finite_difference(self):
        self._run_case("uhf")

    def test_rohf_skeleton_matches_finite_difference(self):
        self._run_case("rohf")


if __name__ == "__main__":
    unittest.main()
