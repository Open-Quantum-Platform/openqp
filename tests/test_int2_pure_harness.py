"""Call the CTest ERI checks across the shared-library boundary."""
import ctypes
import os
import unittest


class PureEriHarnessTests(unittest.TestCase):
    def test_rys_and_rotated_axis_reference_cases(self):
        try:
            from oqp import runtime
            root, source = runtime.resolve_oqp_root()
            library = ctypes.CDLL(runtime.library_path(root, source))
            checks = [getattr(library, name) for name in
                      ("oqp_test_int2_rys_pure", "oqp_test_int2_rotaxis_pure")]
        except (ImportError, OSError, AttributeError, RuntimeError) as exc:
            if os.getenv("OQP_REQUIRE_NATIVE_TESTS") == "1":
                self.fail(f"compiled ERI checks are required: {exc}")
            self.skipTest(str(exc))
        for check in checks:
            check.argtypes = []
            check.restype = None
            check()
