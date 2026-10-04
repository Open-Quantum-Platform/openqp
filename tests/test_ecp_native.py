"""Native ECP integrals (source/ecp_native.F90) against independent references.

Drives the ``oqp_ecp_native_selftest`` bind(C) harness
(tests/fortran/ecp_native_selftest.F90): an O p primitive and a Br s and p
primitive with the LANL2DZ Br ECP at the HOBr O-Br separation.

- An off-centre matrix element is compared with an mpmath (local channel) and
  adaptive-quadrature (semi-local channels) evaluation that does not use the
  Bessel expansion.  libecpint 1.0.7, even with
  external/fix_libecpint_accuracy.py, is 8e-10 away from this value.
- Two on-centre elements are compared with their Gamma-function closed forms.
- An f-g element on fluorine next to the def2 iodine ECP is compared with an
  independent evaluation; libecpint 1.0.7 returned -4.0e-5 instead of -2.70e-4
  for it (g functions on an atom without an ECP).
- The first and second nuclear derivatives are compared with five-point
  finite differences of the value and of the first derivative, and the first
  derivatives must sum to zero over the atoms.
"""

import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NATIVE_TESTS_REQUIRED = os.getenv("OQP_REQUIRE_NATIVE_TESTS") == "1"
_RUNTIME_ERROR = ""


def _runtime_available():
    global _RUNTIME_ERROR
    try:
        os.environ.setdefault("OPENQP_ROOT", str(ROOT))
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp

        if not hasattr(oqp.lib, "oqp_ecp_native_selftest"):
            _RUNTIME_ERROR = "missing native symbol: oqp_ecp_native_selftest"
            return False
        return True
    except Exception as exc:
        _RUNTIME_ERROR = f"{type(exc).__name__}: {exc}"
        return False


RUNTIME_AVAILABLE = _runtime_available()


class NativeEcpBuildGate(unittest.TestCase):
    @unittest.skipUnless(
        NATIVE_TESTS_REQUIRED,
        "source-only run; set OQP_REQUIRE_NATIVE_TESTS=1 after building OpenQP",
    )
    def test_required_runtime_and_symbols(self):
        self.assertTrue(
            RUNTIME_AVAILABLE,
            f"compiled OpenQP runtime is required: {_RUNTIME_ERROR}",
        )


@unittest.skipUnless(RUNTIME_AVAILABLE, f"compiled OpenQP runtime unavailable: {_RUNTIME_ERROR}")
class NativeEcpSelfTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import oqp

        err = oqp.ffi.new("double[9]")
        oqp.lib.oqp_ecp_native_selftest(err)
        cls.err = [err[i] for i in range(9)]

    def test_off_centre_element_matches_independent_quadrature(self):
        self.assertLess(self.err[0], 1e-12)

    def test_on_centre_elements_match_closed_forms(self):
        self.assertLess(self.err[1], 1e-13)
        self.assertLess(self.err[2], 1e-13)

    def test_first_derivatives_match_finite_differences(self):
        self.assertLess(self.err[3], 1e-10)

    def test_second_derivatives_match_finite_differences(self):
        self.assertLess(self.err[4], 1e-10)

    def test_first_derivatives_are_translationally_invariant(self):
        self.assertLess(self.err[5], 1e-14)

    def test_f_g_element_matches_independent_quadrature(self):
        self.assertLess(self.err[6], 1e-13)

    def test_f_g_derivatives_match_finite_differences(self):
        self.assertLess(self.err[7], 1e-9)
        self.assertLess(self.err[8], 1e-7)


if __name__ == "__main__":
    unittest.main()
