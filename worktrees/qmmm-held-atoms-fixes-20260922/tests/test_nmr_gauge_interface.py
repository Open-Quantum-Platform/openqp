"""Interface and claim-boundary tests for NMR gauge handling.

CGO is the validated default; GIAO is the gauge-origin-independent option.
"""

import unittest
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


if __name__ == "__main__":
    unittest.main()
