"""End-to-end regression decks for the native EDIIS/ADIIS QP path."""

import os
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples" / "SCF"


def _runtime_available():
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp
        from oqp.utils.oqp_tester import OQPTester  # noqa: F401

        return hasattr(oqp.lib, "oqp_simplex_qp_solve")
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime not available")
class EDIISADIISRegressionTests(unittest.TestCase):
    def _run_deck(self, name):
        from oqp.utils.mpi_utils import MPIManager
        from oqp.utils.oqp_tester import OQPTester

        with tempfile.TemporaryDirectory(prefix="openqp_simplex_qp_") as tmp:
            tester = OQPTester(
                output_dir=str(Path(tmp) / "results"),
                total_cpus=1,
                omp_threads=1,
                mpi_manager=MPIManager(),
            )
            tester.run(str(EXAMPLES / name))
            self.assertEqual(len(tester.results), 1, tester.results)
            result = tester.results[0]
            self.assertEqual(
                result["status"], "PASSED",
                f"{name} regression failed:\n{result['message']}",
            )

    def test_ediis_reference_deck(self):
        self._run_deck("h2o_rhf_6-31g_pbe_ediis.oqp")

    def test_adiis_reference_deck(self):
        self._run_deck("h2o_rhf_6-31g_pbe_adiis.oqp")


if __name__ == "__main__":
    unittest.main()
