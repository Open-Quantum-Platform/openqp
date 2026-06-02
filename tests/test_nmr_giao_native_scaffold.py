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
        self.assertIn("public giao_h10_core", int1)
        self.assertIn("subroutine giao_h10_core", int1.lower())
        self.assertIn("giao_h10_core_ints", int1)
        self.assertIn("int1_giao_h10_core", int1)
        self.assertIn("comp_giao_h10_core_prim", prim)
        self.assertIn("h10_core = - int1e_ignuc(asym) - int1e_igkin", prim)
        self.assertIn("h10_onee = -0.5*int1e_giao_irjxp", int1)
        self.assertIn("GIAO two-electron Fock derivative", prim)

    def test_native_giao_piece_is_not_wired_as_cgo_fallback(self):
        runfunc = RUNFUNC.read_text()
        self.assertIn("NotImplementedError", runfunc)
        self.assertNotIn("nmr_gauge=cgo for GIAO", runfunc)
        self.assertNotIn("giao_overlap_derivative", runfunc)

    def test_status_ledger_keeps_native_giao_checkpoint_explicit(self):
        status = STATUS.read_text()
        self.assertIn("GIAO development checkpoint", status)
        self.assertIn("GIAO overlap magnetic derivative `S10`", status)
        self.assertIn("Native one-electron building block implemented", status)
        self.assertIn("GIAO first-order core Hamiltonian `h10`", status)
        self.assertIn("int1.F90::giao_h10_core", status)
        self.assertIn("mod_1e_primitives.F90::comp_giao_h10_core_prim", status)
        self.assertIn("make_h10(..., gauge_orig=None)", status)
        self.assertIn("not in the production CGO shielding path", status)
        self.assertIn("first-order core-Hamiltonian magnetic derivative", status)
        self.assertIn("kinetic magnetic-derivative contribution", status)
        self.assertIn("nuclear-attraction magnetic-derivative", status)
        self.assertIn("real-valued storage for the real coefficient of an imaginary", status)
        self.assertIn("Component ordering is Cartesian `(x, y, z)`", status)
        self.assertIn("h10` alone still does not enable native", status)
        self.assertIn("Native one-electron building block implemented; validation in progress", status)
        self.assertIn("h10_onee = -0.5*int1e_giao_irjxp", status)
        self.assertIn("Not implemented; no CGO fallback allowed", status)
        self.assertIn("no CGO fallback allowed", status)
        self.assertIn("S10_a(mu,nu) = 0.5", status)
        self.assertIn("intentionally not", status)

    def test_build_notes_record_verified_blas_recipes(self):
        notes = BUILD_NOTES.read_text()
        for required in (
            "Homebrew GCC/GFortran and native Apple Accelerate",
            "env -u BLAS_LIBRARIES -u LAPACK_LIBRARIES",
            "CC=/opt/homebrew/bin/gcc-15",
            "-DLINALG_LIB=auto -DLINALG_LIB_INT64=OFF -DBLA_VENDOR=Apple",
            "/home/cheolhochoi.guest/venvs/openqp311/bin/python -m pip install . -v",
            "-DLINALG_LIB=OpenBLAS -DLINALG_LIB_INT64=OFF",
            "BLA_SIZEOF_INTEGER=4",
            "INTEGER_SIZE=4",
            "/opt/homebrew`, not `/usr/local`",
            "Do not solve future ABI/build issues by hiding them behind global compiler",
        ):
            self.assertIn(required, notes)


if __name__ == "__main__":
    unittest.main()
