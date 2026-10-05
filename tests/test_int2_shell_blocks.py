"""Native signed-matrix comparison for shell-block Fock and response contractions."""
import ctypes
import os
import unittest
from unittest.mock import patch


class ShellBlockContractions(unittest.TestCase):
    def test_compact_inactive_images(self):
        self._check("int2_td_shell_images_selftest", 16)

    def test_all_shell_permutations_and_spin_channels(self):
        # Exercise the driver's final layout resolution, both cutoffs and
        # 1/2/4 thread images, plus the forced direct/mixed, CAM and FP32 cases.
        # The native check asserts the expected configuration for every kind.
        for layout in (None, "legacy", "shell", "unknown", "x" * 17):
            with self.subTest(layout=layout), patch.dict(os.environ):
                if layout is None:
                    os.environ.pop("OQP_INT2_LAYOUT", None)
                else:
                    os.environ["OQP_INT2_LAYOUT"] = layout
                self._check("int2_shell_block_selftest", 44000)

    def _check(self, symbol, expected_count):
        try:
            import oqp.runtime as runtime
            root, source = runtime.resolve_oqp_root()
            lib = ctypes.CDLL(runtime.library_path(root, source))
            check = getattr(lib, symbol)
        except (ImportError, OSError, AttributeError, RuntimeError) as exc:
            if os.getenv("OQP_REQUIRE_NATIVE_TESTS") == "1":
                self.fail(f"compiled shell-block regression is required: {exc}")
            self.skipTest(str(exc))
        check.argtypes = [ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_int),
                          ctypes.POINTER(ctypes.c_int)]
        check.restype = None
        for _ in range(2):  # Includes repeated allocation, cleanup, and re-entry.
            error, count, failed = ctypes.c_double(), ctypes.c_int(), ctypes.c_int()
            check(ctypes.byref(error), ctypes.byref(count), ctypes.byref(failed))
            self.assertEqual(count.value, expected_count)
            self.assertEqual(failed.value, 0)
            self.assertLessEqual(error.value, 2e-10)


if __name__ == "__main__":
    unittest.main()
