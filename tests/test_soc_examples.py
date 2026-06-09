import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOC_EXAMPLES = ROOT / "examples" / "SOC"


class TestSOCExamples(unittest.TestCase):

    def _run_example(self, inp_file):
        from oqp.utils.oqp_tester import OQPTester
        from oqp.utils.mpi_utils import MPIManager
        mpi_manager = MPIManager()
        tester = OQPTester(
            output_dir='openqp_soc_test_tmp',
            omp_threads=2,
            mpi_manager=mpi_manager,
        )
        tester.run(str(inp_file))
        result = tester.results[0]
        self.assertEqual(
            result['status'], 'PASSED',
            f"{inp_file.name} regression failed:\n{result['message']}"
        )

    def test_soc_example_files_present(self):
        for name in (
            'H2O_BHHLYP_SOC.inp', 'H2O_BHHLYP_SOC.json',
            'CH3Br-BHHLYP-SOC.inp', 'CH3Br-BHHLYP-SOC.json',
        ):
            self.assertTrue((SOC_EXAMPLES / name).is_file(), f"missing: {name}")

    def test_h2o_soc(self):
        self._run_example(SOC_EXAMPLES / "H2O_BHHLYP_SOC.inp")

    def test_ch3br_soc(self):
        self._run_example(SOC_EXAMPLES / "CH3Br-BHHLYP-SOC.inp")


if __name__ == "__main__":
    unittest.main()
