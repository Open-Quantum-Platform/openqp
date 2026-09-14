"""Native TRAH: the descent loop taken at floating-point precision terminates.

Near an ROHF solution whose orbital gradient cannot reach the requested
tolerance at this precision, the loop used to repeat forever without output
(24 DNA thymine NAMD trajectories hung for up to 67 hours at full CPU).  It now
stops once |g| has not halved for a fixed number of steps, while a refinement
that is still converging -- slowly, halving every few steps -- runs to the
requested tolerance.
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
    def test_loop_stops_on_stagnation_and_is_bounded(self):
        src = TRAH_CORE.read_text()
        self.assertIn("integer,  parameter :: fp_stall_steps = 16", src)
        self.assertIn("integer,  parameter :: max_fp_refine = 100", src)
        self.assertRegex(
            src,
            r"do while \(gnorm > par%conv_tol \.and\. snorm > 0\.0_dp \.and\. pred <= pred_floor &\s*\n"
            r"\s*\.and\. n_stall < fp_stall_steps \.and\. n_fp < max_fp_refine\)")
        self.assertRegex(src, r"max_fp_refine\)\s*\n\s*n_fp = n_fp \+ 1\n")
        # only a step that does not halve |g| since its last halving counts towards the stall
        self.assertRegex(src, r"if \(gnorm < 0\.5_dp\*g_ref\) then\s*\n\s*g_ref   = gnorm\s*\n\s*n_stall = 0\s*\n"
                              r"\s*else\s*\n\s*n_stall = n_stall \+ 1")

    def test_block_bound_reached_while_improving_returns_to_the_macro_loop(self):
        """A block that hits its step bound while |g| still halves cycles into
        the nmac-bounded macro loop; a stagnant block stops."""
        src = TRAH_CORE.read_text()
        self.assertRegex(src, r"if \(n_fp >= max_fp_refine \.and\. n_stall < fp_stall_steps \.and\. macro < par%nmac\) then"
                              r"[\s\S]{0,400}?\bcycle\b")


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class TestUnreachableGradientToleranceTerminates(unittest.TestCase):
    """The refinement loop stops, and above the requested tolerance the core no
    longer claims convergence: the log states where it stopped.  The SCF driver
    then applies its own |g| < 1e-4 acceptance and finishes the calculation."""

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

    def test_unreachable_tolerance_stops_without_claiming_convergence(self):
        proc, log = self._run("1e-14")
        self.assertEqual(proc.returncode, 0, proc.stdout.decode(errors="ignore")[-2000:])
        m = re.search(r"\s(\S+)\s+refinement stopped after (\d+) steps above conv", log)
        self.assertIsNotNone(m, "TRAH did not report where the refinement stopped")
        n_steps, g_stop = int(m.group(2)), float(m.group(1))
        self.assertGreater(g_stop, 1e-14)
        # |g| halves every few steps for well over eight steps before it reaches
        # the noise floor, so the refinement must run past any fixed eight-step cap
        self.assertGreater(n_steps, 8)
        entry = re.findall(r"^\s+\d+\s+-?\d+\.\d+\s+(\d\.\d+E[-+]\d+)\s+[-\d.]+\s+[\d.]+\s+\d+\s+acc", log, re.M)
        self.assertTrue(entry, "no accepted TRAH macroiteration in the log")
        self.assertLess(g_stop, 1e-2*float(entry[-1]))
        self.assertNotIn("CONVERGED (FP precision", log)
        self.assertIn("SCF convergence achieved", log)        # the SCF driver's own acceptance


if __name__ == "__main__":
    unittest.main()
