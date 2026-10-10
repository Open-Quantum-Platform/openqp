"""CPHF solves print one convergence summary at the default log level.

An analytic Hessian solves CPHF for every nuclear displacement, and each solve
used to write an initial residual, a residual per iteration, a completion line
per right-hand side and a DFT XC timing line per iteration: 152 of the 688
lines of examples/HESS/H2O_B3LYP5_RPA_ANA_HESS.  Those move to verbose >= 2.
"""
import os
import re
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DECK = """[input]
system=
   1   0.000000000   0.000000000   0.000000000
   1   0.000000000   0.000000000   0.740000000
charge=0
method=hf
functional=bhhlyp
basis=6-31g
runtype=hess
{verbose}
[guess]
type=huckel
save_mol=false
[scf]
type=rhf
multiplicity=1
conv=1.0e-8
[hess]
type=analytical
state=0
clean=True
"""

OPEN_SHELL_DECK = """[input]
system=
   O   0.000000000   0.000000000   0.000000000
   H   0.000000000   0.000000000   0.970000000
charge=0
method=hf
basis=6-31g
runtype=hess
[guess]
type=huckel
save_mol=false
[scf]
type={reference}
multiplicity=2
conv=1.0e-8
[hess]
type=analytical
state=0
clean=True
"""

SUMMARY = re.compile(r"converged +\d+ of +\d+ right-hand sides in +\d+ - +\d+ iterations")


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime unavailable")
class CphfLogLevel(unittest.TestCase):
    def _run(self, deck):
        from oqp.pyoqp import Runner
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "case.inp"
            inp.write_text(deck)
            log = Path(tmp) / "case.log"
            # Native Hessian kernels write fort.6 in the working directory; run
            # there so parallel tests cannot read each other's capture.
            os.chdir(tmp)
            try:
                Runner(project="case", input_file=str(inp), log=str(log),
                       silent=1, usempi=False).run()
            finally:
                os.chdir(cwd)
            return log.read_text(errors="replace")

    def _log(self, verbose=None):
        return self._run(DECK.format(verbose="" if verbose is None else f"verbose={verbose}"))

    def test_default_level_prints_one_summary_per_solve(self):
        text = self._log()
        self.assertGreater(len(SUMMARY.findall(text)), 0)
        self.assertNotIn("CPHF ITER RHS", text)
        self.assertNotIn("INITIAL CPHF ERROR", text)
        self.assertIsNone(re.search(r"CPHF RHS +\d+ completed", text))
        self.assertNotIn("DFT XC integration time", text)

    def test_detailed_level_keeps_the_residuals(self):
        import oqp
        text = self._log(2)
        self.assertIn("CPHF ITER RHS", text)
        self.assertIsNotNone(re.search(r"CPHF RHS +\d+ completed", text))
        if bool(oqp.lib.oqp_have_openmp()):
            # the timing line is only compiled into OpenMP builds
            self.assertIn("DFT XC integration time", text)

    def test_open_shell_solvers_report_their_summary(self):
        """The UHF and ROHF solvers write through the captured fort.6 unit, so
        the summary is lost unless it is flushed before the capture is read."""
        for reference in ("uhf", "rohf"):
            text = self._run(OPEN_SHELL_DECK.format(reference=reference))
            self.assertEqual(len(SUMMARY.findall(text)), 1, reference)
            self.assertIsNone(re.search(r"(UHF|ROHF) CPHF RHS +\d+ completed", text), reference)


if __name__ == "__main__":
    unittest.main()
