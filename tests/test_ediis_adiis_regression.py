"""End-to-end regression decks for the native EDIIS/ADIIS QP path."""

import os
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples" / "SCF"
NATIVE_TESTS_REQUIRED = os.getenv("OQP_REQUIRE_NATIVE_TESTS") == "1"
_RUNTIME_ERROR = ""


def _runtime_available():
    global _RUNTIME_ERROR
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp
        from oqp.utils.oqp_tester import OQPTester  # noqa: F401

        missing = [
            name for name in (
                "oqp_simplex_qp_solve",
                "oqp_simplex_qp_solve_avoid",
            )
            if not hasattr(oqp.lib, name)
        ]
        if missing:
            _RUNTIME_ERROR = f"missing native symbols: {', '.join(missing)}"
            return False
        return True
    except Exception as exc:
        _RUNTIME_ERROR = f"{type(exc).__name__}: {exc}"
        return False


RUNTIME_AVAILABLE = _runtime_available()


class NativeEDIISADIISBuildGate(unittest.TestCase):
    @unittest.skipUnless(
        NATIVE_TESTS_REQUIRED,
        "source-only run; set OQP_REQUIRE_NATIVE_TESTS=1 after building OpenQP",
    )
    def test_required_runtime_symbols_and_tester_import(self):
        self.assertTrue(
            RUNTIME_AVAILABLE,
            f"compiled EDIIS/A-DIIS runtime is required: {_RUNTIME_ERROR}",
        )


@unittest.skipUnless(RUNTIME_AVAILABLE, "compiled OpenQP runtime not available")
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
