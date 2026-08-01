"""Static guard against call-history dependent MRSF NAC amplitudes."""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


class MRSFNACTDPIsolationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = (
            ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
        ).read_text()
        cls.body = source.split(
            "subroutine mrsf_nac_amp(infos, only_istate, only_jstate)", 1
        )[1].split(
            "end subroutine mrsf_nac_amp", 1
        )[0]

    def test_amplitude_engine_does_not_read_diagonal_gradient_scratch(self):
        self.assertNotIn("OQP_td_p", self.body)
        self.assertNotIn('"OQP::td_p"', self.body)
        self.assertNotIn("OQP_NAC_AMP_NOP2", self.body)
        self.assertNotIn("has_records", self.body)

    def test_amplitude_engine_owns_a_zero_p2_channel(self):
        self.assertIn("pIJ(nbf,nbf,2)", self.body)
        self.assertIn("source=0.0_dp", self.body)
        self.assertIn("p2 = pIJ", self.body)
        assignments = re.findall(r"(?mi)^\s*pIJ(?:\([^\n]*\))?\s*=", self.body)
        self.assertEqual(assignments, [])


if __name__ == "__main__":
    unittest.main()
