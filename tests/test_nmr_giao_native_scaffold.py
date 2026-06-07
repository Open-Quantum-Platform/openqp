"""Source guards for the first native OpenQP GIAO NMR building block.

These tests do not claim a validated OpenQP GIAO shielding implementation.  They
only require that native London/GIAO overlap-derivative code exists separately
from the CGO path while the user-facing nmr_gauge=giao gate remains closed.
"""
from __future__ import annotations

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INT1 = ROOT / "source" / "integrals" / "int1.F90"
PRIM = ROOT / "source" / "integrals" / "mod_1e_primitives.F90"
RUNFUNC = ROOT / "pyoqp" / "oqp" / "library" / "runfunc.py"
STATUS = ROOT / "source" / "modules" / "NMR_SHIELDING_STATUS.md"
BUILD_NOTES = ROOT / "source" / "modules" / "NMR_BUILD_NOTES.md"


class NMRGIAONativeScaffoldTests(unittest.TestCase):
    def test_native_giao_overlap_derivative_symbols_exist(self):
        int1 = INT1.read_text()
        prim = PRIM.read_text()
        self.assertIn("public giao_overlap_derivative", int1)
        self.assertIn("subroutine giao_overlap_derivative", int1.lower())
        self.assertIn("giao_overlap_deriv_ints", int1)
        self.assertIn("comp_giao_overlap_deriv_prim", prim)
        self.assertIn("S10_a(mu,nu)", prim)
        self.assertIn("(R_mu - R_nu) x <mu|r|nu>", prim)

    def test_native_giao_h10_core_symbols_exist(self):
        int1 = INT1.read_text()
        prim = PRIM.read_text()
        debug = (ROOT / "source" / "modules" / "nmr_giao_debug.F90").read_text()
        self.assertIn("public giao_h10_core", int1)
        self.assertIn("subroutine giao_h10_core", int1.lower())
        self.assertIn("giao_h10_core_ints", int1)
        self.assertIn("int1_giao_h10_core", int1)
        self.assertIn("comp_giao_h10_core_prim", prim)
        self.assertIn("h10_core = - int1e_ignuc(asym) - int1e_igkin", prim)
        self.assertIn("h10_onee = -0.5*int1e_giao_irjxp", int1)
        self.assertIn("GIAO two-electron Fock derivative", prim)
        self.assertIn("public nmr_giao_h10_twoe_debug", debug)
        self.assertIn("subroutine nmr_giao_h10_twoe_debug", debug.lower())
        self.assertIn("GIAO_H10_TWOE_DEBUG_FULL", debug)

    def test_native_giao_piece_is_not_wired_as_cgo_fallback(self):
        runfunc = RUNFUNC.read_text()
        self.assertIn("NotImplementedError", runfunc)
        self.assertNotIn("nmr_gauge=cgo for GIAO", runfunc)
        self.assertNotIn("giao_overlap_derivative", runfunc)

    def test_build_notes_keep_the_abi_caution(self):
        # The detailed build recipes are environment-specific and may evolve;
        # only require that the ABI/build caution note ships with the NMR code.
        notes = BUILD_NOTES.read_text()
        self.assertIn(
            "Do not solve future ABI/build issues by hiding them behind global compiler",
            notes,
        )


if __name__ == "__main__":
    unittest.main()
