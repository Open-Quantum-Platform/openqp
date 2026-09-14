"""CASSCF follows the global log level: verbose = 0 drops the macroiteration table.

The CASSCF log is written from Python (``CASSCF._write_log``), so the native
print gates never saw it and a quiet run still listed every macroiteration.
"""
import os
import re
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DECK = ROOT / "examples" / "WF_methods" / "H2O_CASSCF_CAS44.inp"


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return True
    except Exception:
        return False


@unittest.skipUnless(DECK.exists() and _runtime_available(), "compiled OpenQP runtime unavailable")
class CasscfLogLevel(unittest.TestCase):
    def _log(self, verbose=None):
        from oqp.pyoqp import Runner
        text = DECK.read_text()
        if verbose is not None:
            text = re.sub(r"^(\[input\][^\n]*\n)", rf"\g<1>verbose={verbose}\n", text,
                          count=1, flags=re.M)
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "cas.inp"
            inp.write_text(text)
            log = Path(tmp) / "cas.log"
            os.chdir(tmp)
            try:
                Runner(project="cas", input_file=str(inp), log=str(log),
                       silent=1, usempi=False).run()
            finally:
                os.chdir(cwd)
            return log.read_text(errors="replace")

    def test_quiet_level_keeps_only_the_summary(self):
        normal = self._log()
        quiet = self._log(0)
        self.assertIn("--- macro iterations ---", normal)      # control: the table exists
        self.assertNotIn("--- macro iterations ---", quiet)
        for text in (normal, quiet):
            self.assertIn("PyOQP CASSCF converged:", text)
            self.assertIn("PyOQP CASSCF macro iterations:", text)


if __name__ == "__main__":
    unittest.main()
