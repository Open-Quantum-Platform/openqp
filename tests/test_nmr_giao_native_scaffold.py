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

    def test_native_giao_piece_is_not_wired_as_cgo_fallback(self):
        runfunc = RUNFUNC.read_text()
        self.assertIn("NotImplementedError", runfunc)
        self.assertNotIn("nmr_gauge=cgo for GIAO", runfunc)
        self.assertNotIn("giao_overlap_derivative", runfunc)


if __name__ == "__main__":
    unittest.main()
