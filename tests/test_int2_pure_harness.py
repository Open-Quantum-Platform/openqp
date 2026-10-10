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
from oqp import runtime
# Initialize each child's DLL search directories and package preloads.
runtime.library_path(sys.argv[3])
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
                # The rys check alone takes ~10-20 s on one core; under
                # `pytest -n auto` it shares a small CI runner with workers
                # running multithreaded OpenQP jobs, and 120 s was exceeded
                # repeatedly on ubuntu-24.04 x86 (MPI OFF) without any change
                # to the ERI code.  The limit only guards against a hang.
                result = subprocess.run(
                    [sys.executable, "-c", script, str(library_path), name, str(root)],
                    capture_output=True, text=True, timeout=600)
                self.assertEqual(result.returncode, 0,
                                 f"{name} exited {result.returncode}:\n"
                                 + result.stdout + result.stderr)
