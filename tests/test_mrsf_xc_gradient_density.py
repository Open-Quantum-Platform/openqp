import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MRSF_GRADIENT = ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
MRSF_Z_VECTOR = ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90"
XC_GRADIENT = ROOT / "source" / "dftlib" / "dft_gridint_tdxc_grad.F90"


def _call_block(source: str, callee: str) -> str:
    start = re.search(rf"call\s+{callee}\s*\(", source, re.IGNORECASE)
    if not start:
        raise AssertionError(f"Could not find call to {callee}")
    tail = source[start.start():]
    end = re.search(r"\n\s*infos\s*=\s*infos\s*\)", tail, re.IGNORECASE)
    if not end:
        raise AssertionError(f"Could not find end of {callee} keyword call")
    return tail[: end.end()]


class MrsfXcGradientDensityTests(unittest.TestCase):
    def test_utddft_optional_xa_xb_is_the_fxc_gate(self):
        source = XC_GRADIENT.read_text()
        self.assertRegex(source, r"doFxc\s*=\s*present\(xa\)")
        self.assertRegex(source, r"if \(doFxc\) nxcder\s*=\s*3")
        self.assertRegex(source, r"dat%xa\s*=>\s*xa")
        self.assertRegex(source, r"dat%xb\s*=>\s*xb")

    def test_mrsf_z_vector_preserves_legacy_channels_and_adds_spin_resolved_xc_candidate(self):
        source = MRSF_Z_VECTOR.read_text()
        self.assertRegex(source, r"reserve_data\(OQP_td_mrsf_density,\s*TA_TYPE_REAL64,\s*nbf\*nbf\*9,\s*\(/\s*9,\s*nbf,\s*nbf\s*/\)")
        self.assertRegex(source, r"call\s+mrsfcbc\(infos,\s*mo_a,\s*mo_a,\s*wrk1,\s*fmrst1\(1,:,:,:\)\)")
        self.assertRegex(source, r"td_mrsf_den\(1:7,:,:\)\s*=\s*fmrst1\(1,1:7,:,:\)")
        self.assertRegex(source, r"call\s+umrsfcbc\(infos,\s*mo_a,\s*mo_a,\s*wrk1,\s*umrsf_xc_den\)")
        self.assertRegex(source, r"td_mrsf_den\(8,:,:\)\s*=\s*umrsf_xc_den\(1,:,:\)\s*\+\s*umrsf_xc_den\(3,:,:\)\s*\+\s*&?\s*umrsf_xc_den\(5,:,:\)\s*\+\s*umrsf_xc_den\(7,:,:\)")
        self.assertRegex(source, r"td_mrsf_den\(9,:,:\)\s*=\s*umrsf_xc_den\(2,:,:\)\s*\+\s*umrsf_xc_den\(4,:,:\)\s*\+\s*&?\s*umrsf_xc_den\(6,:,:\)\s*\+\s*umrsf_xc_den\(8,:,:\)")
        self.assertNotRegex(source, r"td_mrsf_den\(8,:,:\)\s*=\s*td_abxc")
        self.assertNotRegex(source, r"td_mrsf_den\(9,:,:\)\s*=\s*td_abxc")

    def test_mrsf_gradient_passes_only_candidate_channels_to_xc_xa_xb(self):
        source = MRSF_GRADIENT.read_text()
        block = _call_block(source, "utddft_xc_gradient")
        self.assertRegex(block, r"pa\s*=\s*p\(:,:,1:1\)")
        self.assertRegex(block, r"pb\s*=\s*p\(:,:,2:2\)")
        self.assertRegex(block, r"xa\s*=\s*spc\(:,:,8:8\)")
        self.assertRegex(block, r"xb\s*=\s*spc\(:,:,9:9\)")
        self.assertNotRegex(block, r"xa\s*=\s*td_abxc")
        self.assertNotRegex(block, r"xb\s*=\s*td_abxc")
        self.assertNotRegex(block, r"xa\s*=\s*v\(:,:,1:1\)")
        self.assertNotRegex(block, r"xb\s*=\s*v\(:,:,2:2\)")


if __name__ == "__main__":
    unittest.main()
