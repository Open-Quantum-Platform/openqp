"""Static and reference-fixture checks for the NMR GIAO validation path.

These tests deliberately do not claim OpenQP GIAO is implemented. They guard the
benchmark-driven workflow: a distinct PySCF GIAO oracle exists, PySCF is not used
as the OpenQP implementation path, the benchmark driver records OpenQP GIAO as a
gated status until native integrals are connected, and the conservative Overleaf
claim boundary remains intact.
"""
from __future__ import annotations

import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "nmr" / "pyscf_giao_reference.json"
GENERATOR = ROOT / "tests" / "fixtures" / "nmr" / "generate_pyscf_giao_reference.py"
DRIVER = ROOT / "scripts" / "nmr_giao_benchmark_matrix.py"


class NMRGIAOBenchmarkPathTests(unittest.TestCase):
    def test_giao_reference_fixture_is_distinct_from_cgo_oracle(self):
        self.assertTrue(FIXTURE.exists(), "Generate pyscf_giao_reference.json before validating GIAO")
        data = json.loads(FIXTURE.read_text())
        self.assertEqual(data["gauge_formulation"], "giao")
        self.assertIsNone(data["gauge_origin"])
        self.assertIn("gauge_orig=None", data["provenance"])
        self.assertIn("HF", data["results"])
        self.assertIn("pbe", data["results"])
        hf = data["results"]["HF"]
        self.assertIn("sigma_total_coupled_tensor", hf)
        self.assertEqual(len(hf["sigma_total_coupled_tensor"]), 3)
        # This value is intentionally different from the CGO coupled HF oxygen
        # total shielding (~286.79 ppm in the CGO fixture), proving this is not
        # a CGO fixture copy.
        self.assertGreater(abs(hf["sigma_total_coupled"][0] - 286.79), 50.0)

    def test_generator_uses_pyscf_as_oracle_not_openqp_implementation(self):
        text = GENERATOR.read_text()
        self.assertIn("gauge_orig=None", text)
        self.assertIn("use as oracle only", text)
        self.assertNotIn("oqp.nmr_shielding", text)
        self.assertNotIn("nmr_gauge=giao", text)

    def test_benchmark_driver_records_gated_openqp_giao_without_cgo_fallback(self):
        text = DRIVER.read_text()
        self.assertIn("OpenQP GIAO is expected to remain gated", text)
        self.assertIn('"status": "gated"', text)
        self.assertIn("not yet implemented", text)
        self.assertNotIn("nmr_gauge=cgo for GIAO", text)

    def test_benchmark_results_capture_reference_and_gated_status(self):
        path = ROOT / "tests" / "fixtures" / "nmr" / "benchmark_results" / "nmr_giao_benchmark_results.json"
        self.assertTrue(path.exists(), "Run scripts/nmr_giao_benchmark_matrix.py before reporting benchmarks")
        data = json.loads(path.read_text())
        self.assertEqual(data["status"], "reference_ready_openqp_giao_gated")
        gauges = {(row["backend"], row["gauge"]) for row in data["rows"]}
        self.assertIn(("pyscf", "giao"), gauges)
        self.assertIn(("pyscf", "cgo"), gauges)
        self.assertIn(("openqp", "giao"), gauges)
        self.assertTrue(any(row["backend"] == "pyscf" and row["gauge"] == "giao" and row["status"] == "ok" for row in data["rows"]))
        self.assertTrue(any(row["backend"] == "openqp" and row["gauge"] == "giao" and row["status"] in {"gated", "skipped"} for row in data["rows"]))


if __name__ == "__main__":
    unittest.main()
