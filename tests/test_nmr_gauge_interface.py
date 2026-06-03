"""Interface and claim-boundary tests for NMR gauge handling.

CGO is the validated baseline and must remain the default.  GIAO is exposed as
an explicit, gated development option until full gauge-including integrals and
benchmarks are implemented.
"""

import unittest
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class NMRGaugeInterfaceTests(unittest.TestCase):
    def test_config_schema_has_cgo_default_and_giao_option(self):
        text = (ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py").read_text()
        self.assertIn("'nmr_gauge'", text)
        self.assertIn("'default': 'cgo'", text)

    def test_input_checker_accepts_only_cgo_or_giao(self):
        text = (ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py").read_text()
        self.assertIn("NMR_GAUGES", text)
        self.assertIn('"cgo"', text)
        self.assertIn('"giao"', text)
        self.assertIn("properties.nmr_gauge", text)

    def test_manuscript_uses_conservative_giao_wording(self):
        main = Path("/Users/cheolhochoi/Dropbox/Documents/overleaf/hf-dft-nmr/main.tex")
        if not main.exists():
            self.skipTest("Overleaf hf-dft-nmr manuscript not present")
        text = main.read_text()
        self.assertIn("common gauge-origin formulation", text)
        self.assertRegex(
            text,
            r"gauge-including atomic orbitals (will be introduced|are being introduced)",
        )
        self.assertNotIn("implementation now supports gauge-including atomic orbitals", text)

    def test_giao_benchmark_matrix_tracks_required_comparisons(self):
        matrix = json.loads((ROOT / "tests" / "fixtures" / "nmr" / "giao_benchmark_matrix.json").read_text())
        self.assertEqual(matrix["baseline"], "cgo")
        self.assertEqual(matrix["candidate"], "giao")
        for required in (
            "isotropic_shielding_ppm",
            "tensor_components_ppm",
            "gauge_origin_translation_delta_ppm",
            "wall_time_seconds",
            "peak_memory_mb",
            "trusted_reference_delta_ppm",
        ):
            self.assertIn(required, matrix["metrics"])
        self.assertGreaterEqual(len(matrix["small_regression_systems"]), 3)
        self.assertGreaterEqual(len(matrix["origin_dependence_systems"]), 2)


if __name__ == "__main__":
    unittest.main()
