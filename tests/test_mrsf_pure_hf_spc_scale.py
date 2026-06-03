"""Regression guard for the pure-HF MRSF-TDDFT analytic gradient.

Background
----------
The MRSF spin-pair-coupling scales (spc_coco/ovov/coov) default to the response
HF-exchange scale infos%tddft%HFscale, which in turn defaults to the reference
infos%dft%HFscale. For a pure-HF reference (functional=hf, no DFT functional)
infos%dft%HFscale is left at its -1.0 sentinel, so without intervention the
spin-pair coupling inherits -1.0.

The MRSF energy tolerates this (its fmrst2 rescale is skipped because
spc == HFscale either way), but the MRSF gradient kernel
(grd2_mrsf_compute_data_t_get_density) multiplies the spin-pair-coupling
density terms by spc directly. With spc = -1.0 instead of the correct +1.0,
the coupling gradient terms came out with the wrong sign/magnitude, so the
pure-HF MRSF analytic gradient disagreed with finite differences by ~1e-2
(both singlet and triplet) even though the excitation energy was correct.

Fix: tdhf_mrsf_energy resolves the response HF-exchange scale to 1.0 for a
pure-HF (non-DFT) reference, so the spin-pair coupling defaults to +1.0.

Validation at the time of the fix (H2O/6-31G*, central FD): mrsf-s/mrsf-t hf
analytic-vs-FD max|diff| dropped from ~1.0e-2 to ~5e-6; all DFT-functional MRSF
gradients and the MRSF energy were unchanged.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENERGY_SRC = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"


class MrsfPureHfSpcScaleTests(unittest.TestCase):
    def test_response_hfscale_resolved_for_pure_hf(self):
        """Non-DFT (pure HF) reference must resolve tddft%HFscale to 1.0
        before the spin-pair-coupling defaults are assigned, so the coupling
        does not inherit the -1.0 sentinel that breaks the MRSF gradient."""
        text = ENERGY_SRC.read_text()

        resolve = re.search(
            r"\.not\.\s*dft\s*\)\s*infos%tddft%HFscale\s*=\s*1\.0_dp",
            text,
        )
        self.assertIsNotNone(
            resolve,
            "tdhf_mrsf_energy must set infos%tddft%HFscale = 1.0_dp for a "
            "non-DFT (pure HF) reference so the spin-pair coupling resolves "
            "to +1.0 instead of the -1.0 sentinel.",
        )

        # It must come before the spc_coco default assignment.
        spc_default = text.find("infos%tddft%spc_coco==-1.0_dp")
        self.assertNotEqual(spc_default, -1, "spc_coco default not found")
        self.assertLess(
            resolve.start(),
            spc_default,
            "The HFscale resolution must precede the spin-pair-coupling "
            "defaults so spc_coco/ovov/coov pick up the resolved 1.0.",
        )


if __name__ == "__main__":
    unittest.main()
