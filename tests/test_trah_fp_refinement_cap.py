"""Native TRAH: the descent loop taken at floating-point precision terminates.

Near an ROHF solution whose orbital gradient cannot reach the requested
tolerance at this precision, the loop used to repeat forever without output
(24 DNA thymine NAMD trajectories hung for up to 67 hours at full CPU).  It now
stops once |g| has not dropped by a quarter for a fixed number of steps, while
a refinement that is still converging -- even slowly -- runs to the requested
tolerance.
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
trh_impl=native
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
    def test_loop_stops_on_stagnation_and_is_bounded(self):
        src = TRAH_CORE.read_text()
        self.assertIn("integer,  parameter :: fp_stall_steps = 16", src)
        self.assertIn("real(dp), parameter :: fp_progress = 0.75_dp", src)
        self.assertIn("integer,  parameter :: max_fp_refine = 100", src)
        self.assertRegex(
            src,
            r"do while \(gnorm > par%conv_tol \.and\. snorm > 0\.0_dp \.and\. pred <= pred_floor &\s*\n"
            r"\s*\.and\. n_stall < fp_stall_steps \.and\. n_fp < max_fp_refine\)")
        self.assertRegex(src, r"max_fp_refine\)\s*\n\s*n_fp = n_fp \+ 1\n")
        # only a step without progress (a quarter drop since the last one) counts towards the stall
        self.assertRegex(src, r"if \(gnorm < fp_progress\*g_ref\) then\s*\n\s*g_ref   = gnorm\s*\n\s*n_stall = 0\s*\n"
                              r"\s*else\s*\n\s*n_stall = n_stall \+ 1")

    def test_block_bound_reached_while_improving_returns_to_the_macro_loop(self):
        """A block that hits its step bound while |g| still converges cycles into
        the nmac-bounded macro loop; a stagnant block stops."""
        src = TRAH_CORE.read_text()
        self.assertRegex(src, r"if \(n_fp >= max_fp_refine \.and\. n_stall < fp_stall_steps \.and\. macro < par%nmac\) then"
                              r"[\s\S]{0,400}?\bcycle\b")


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestTightGradientToleranceTerminates(unittest.TestCase):
    """Tight molecular targets may converge after full-Fock refinement. Success
    requires the measured residual; an unreachable target must still fail."""

    def _run(self, conv):
        with tempfile.TemporaryDirectory() as tmp:
            deck = Path(tmp) / "h2o_trah_tight.inp"
            deck.write_text(DECK.replace("conv=1e-14", f"conv={conv}"))
            env = dict(os.environ, OMP_NUM_THREADS="2")
            try:
                proc = subprocess.run([sys.executable, "-m", "oqp.pyoqp", deck.name], cwd=tmp, env=env,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=240)
            except subprocess.TimeoutExpired:
                self.fail(f"TRAH did not finish within 240 s at conv={conv}")
            return proc, (Path(tmp) / "h2o_trah_tight.log").read_text(errors="ignore")

    def test_tight_tolerance_uses_the_fresh_fock_residual(self):
        proc, log = self._run("1e-14")
        diagnostic = proc.stdout.decode(errors="ignore")[-2000:] + "\n" + log[-4000:]
        residuals = re.findall(r"final fresh-Fock residual\s*=\s*([0-9.E+-]+)", log)
        if proc.returncode == 0:
            self.assertTrue(residuals, diagnostic)
            self.assertLess(float(residuals[-1]), 1e-14, diagnostic)
            self.assertIn("SCF convergence achieved", log)
        else:
            # Whether 1e-14 is reachable depends on the numerical libraries.
            # A platform that cannot reach it must report failure honestly.
            self.assertIn("SCF did not converge: TRAH failed the requested criterion", log, diagnostic)
            self.assertNotIn("SCF convergence achieved", log, diagnostic)
            self.assertIn("SCF energy is not converged", log, diagnostic)

    def test_unreachable_tolerance_stops_without_claiming_convergence(self):
        proc, log = self._run("1e-30")
        diagnostic = proc.stdout.decode(errors="ignore")[-2000:] + "\n" + log[-4000:]
        self.assertNotEqual(proc.returncode, 0, diagnostic)
        self.assertIn("SCF did not converge: TRAH failed the requested criterion", log, diagnostic)
        self.assertIn("SCF energy is not converged", log, diagnostic)
        self.assertNotIn("SCF convergence achieved", log, diagnostic)
        self.assertNotIn("CONVERGED (FP precision", log, diagnostic)
        # Blocks may continue while making progress, but each block remains
        # bounded and the subprocess timeout detects an unbounded refinement.
        stopped = re.findall(r"\s(\S+)\s+refinement stopped after (\d+) steps above conv", log)
        self.assertTrue(stopped, diagnostic)
        for residual, steps in stopped:
            self.assertGreater(float(residual), 1e-30)
            self.assertLessEqual(int(steps), 100)


if __name__ == "__main__":
    unittest.main()
