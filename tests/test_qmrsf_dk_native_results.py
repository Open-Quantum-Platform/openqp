"""Unit coverage for the compact OpenQP-native QMRSF-DK result handoff."""

import tempfile
import unittest
from pathlib import Path
import importlib.util


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "pyoqp" / "oqp" / "library" / "qmrsf_results.py"
NATIVE_SOURCE = ROOT / "source" / "modules" / "qmrsf_dk_paper_native.F90"
ENTRY_SOURCE = ROOT / "source" / "modules" / "tdhf_qmrsf_dk.F90"
SPEC = importlib.util.spec_from_file_location("qmrsf_results_under_test", MODULE_PATH)
QMRSF_RESULTS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(QMRSF_RESULTS)
HARTREE_TO_EV = QMRSF_RESULTS.HARTREE_TO_EV
build_qmrsf_dk_results = QMRSF_RESULTS.build_qmrsf_dk_results
format_qmrsf_dk_log_table = QMRSF_RESULTS.format_qmrsf_dk_log_table
parse_qmrsf_dk_dump = QMRSF_RESULTS.parse_qmrsf_dk_dump


class QmrsfDkNativeResultsTests(unittest.TestCase):
    def test_native_entry_point_and_paper_exchange_seam(self):
        native = NATIVE_SOURCE.read_text().lower()
        entry = ENTRY_SOURCE.read_text().lower()
        self.assertIn("call qmrsf_dk_paper_native(inf)", entry)
        self.assertNotIn("import pyscf", native)
        self.assertIn("call native_exchange_eri(eri_act,eri_k)", native)
        self.assertIn("hdet = hdet - (1.0_dp-c_h)*hkdet", native)
        self.assertIn("if (p==q .or. r==s) eri_k(p,q,r,s)=0.0_dp", native)

    def _write_dump(self, directory):
        path = Path(directory) / "qmrsf_dk_full_live.dat"
        singlets = [-1.20 + 0.01 * i for i in range(20)]
        triplets = [-1.10 + 0.01 * i for i in range(15)]
        path.write_text(
            "QMRSF_DK_NATIVE_V1 4 20 15 1\n"
            "0.5000000000000000 0.5000000000000000\n"
            "-1.0000000000000000\n"
            + " ".join("%.16f" % value for value in singlets) + "\n"
            + " ".join("%.16f" % value for value in triplets) + "\n"
            "-1.0000000000000000\n"
            "1.0e-15 2.0e-14\n"
        )
        return path

    def test_parse_and_paper_energy_law(self):
        with tempfile.TemporaryDirectory() as tmp:
            dump = parse_qmrsf_dk_dump(self._write_dump(tmp))

        self.assertEqual(dump["format"], "native_v1")
        self.assertEqual((dump["nsinglet"], dump["ntriplet"], dump["nquintet"]),
                         (20, 15, 1))

        reference = -10.0
        results = build_qmrsf_dk_results(dump, reference)
        self.assertEqual(len(results["states"]), 20)
        self.assertEqual(len(results["triplet_states"]), 15)
        self.assertEqual(len(results["quintet_states"]), 1)
        self.assertAlmostEqual(results["states"][0]["E_DK"], -10.2, places=12)
        self.assertAlmostEqual(results["states"][1]["exc_DK_eV"],
                               0.01 * HARTREE_TO_EV, places=10)
        self.assertAlmostEqual(results["quintet_states"][0]["E_DK"], reference,
                               places=12)
        self.assertTrue(results["diagnostics"]["pass"])
        self.assertAlmostEqual(
            results["diagnostics"]["discarded_cross_spin_block_max"], 2.0e-14
        )

        log = format_qmrsf_dk_log_table(results, max_states=2)
        self.assertIn("OpenQP-native paper matrix", log)
        self.assertIn("20 singlets, 15 triplets, 1 quintet", log)

    def test_rejects_inconsistent_quintet_anchor(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write_dump(tmp)
            text = path.read_text().replace(
                "-1.0000000000000000\n1.0e-15",
                "-0.9000000000000000\n1.0e-15",
            )
            path.write_text(text)
            with self.assertRaisesRegex(ValueError, "quintet anchor"):
                parse_qmrsf_dk_dump(path)


if __name__ == "__main__":
    unittest.main()
