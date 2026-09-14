"""Native TRAH: the descent loop taken at floating-point precision is capped.

Near an ROHF solution whose orbital gradient cannot reach the requested
tolerance at this precision, the loop used to repeat forever without output
(24 DNA thymine NAMD trajectories hung for up to 67 hours at full CPU).
"""
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TRAH_CORE = ROOT / "source" / "trah_core.F90"

DECK = """[input]
system=
   8   0.000000000   0.000000000   0.000000000
   1   0.000000000  -0.757000000   0.587000000
   1   0.000000000   0.757000000   0.587000000
charge=0
runtype=energy
basis=6-31g*
functional=bhhlyp
method=hf

[guess]
type=huckel

[scf]
multiplicity=3
type=rohf
converger_type=trah
conv=1e-14
"""


def _runtime_available():
    try:
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


class TestTrahRefinementLoopSource(unittest.TestCase):
    def test_loop_carries_an_iteration_cap(self):
        src = TRAH_CORE.read_text()
        self.assertIn("integer,  parameter :: max_fp_refine = 8", src)
        self.assertRegex(
            src,
            r"do while \(gnorm > par%conv_tol \.and\. snorm > 0\.0_dp \.and\. pred <= pred_floor &\s*\n"
            r"\s*\.and\. n_fp < max_fp_refine\)")


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestUnreachableGradientToleranceTerminates(unittest.TestCase):
    def test_rohf_triplet_trah_with_unreachable_tolerance_finishes(self):
        with tempfile.TemporaryDirectory() as tmp:
            deck = Path(tmp) / "h2o_trah_tight.inp"
            deck.write_text(DECK)
            env = dict(os.environ, OMP_NUM_THREADS="2")
            try:
                proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", deck.name], cwd=tmp, env=env,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=240)
            except subprocess.TimeoutExpired:
                self.fail("TRAH did not finish within 240 s with an unreachable gradient tolerance")
            log = (Path(tmp) / "h2o_trah_tight.log").read_text(errors="ignore")
        self.assertEqual(proc.returncode, 0, proc.stdout.decode(errors="ignore")[-2000:])
        m = re.search(r"CONVERGED \(FP precision, (\d+) refinement steps\)", log)
        self.assertIsNotNone(m, "TRAH did not report convergence at floating-point precision")
        self.assertLessEqual(int(m.group(1)), 8)


if __name__ == "__main__":
    unittest.main()
