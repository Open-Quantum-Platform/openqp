"""Run the CTest ERI checks in isolated shared-library subprocesses."""
import ctypes
import os
import subprocess
import sys
import unittest


class PureEriHarnessTests(unittest.TestCase):
    def test_rys_and_rotated_axis_reference_cases(self):
        names = ("oqp_test_int2_rys_pure", "oqp_test_int2_rotaxis_pure")
        try:
            from oqp import runtime
            root, source = runtime.resolve_oqp_root()
            library_path = runtime.library_path(root, source)
            library = ctypes.CDLL(library_path)
            for name in names:
                getattr(library, name)
        except (ImportError, OSError, AttributeError, RuntimeError) as exc:
            if os.getenv("OQP_REQUIRE_NATIVE_TESTS") == "1":
                self.fail(f"compiled ERI checks are required: {exc}")
            self.skipTest(str(exc))
        script = """
import ctypes, sys
library = ctypes.CDLL(sys.argv[1])
check = getattr(library, sys.argv[2])
check.argtypes = []
check.restype = None
check()
"""
        for name in names:
            with self.subTest(check=name):
                # Fortran ERROR STOP terminates its process. Keep pytest alive
                # and retain the native diagnostic as an ordinary test failure.
                result = subprocess.run(
                    [sys.executable, "-c", script, str(library_path), name],
                    capture_output=True, text=True, timeout=120)
                self.assertEqual(result.returncode, 0,
                                 f"{name} exited {result.returncode}:\n"
                                 + result.stdout + result.stderr)
