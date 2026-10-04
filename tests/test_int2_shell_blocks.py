"""Native signed-matrix comparison for shell-block Fock and response contractions."""
import ctypes
import os
import unittest


class ShellBlockContractions(unittest.TestCase):
    def test_all_shell_permutations_and_spin_channels(self):
        try:
            import oqp.runtime as runtime
            root, source = runtime.resolve_oqp_root()
            lib = ctypes.CDLL(runtime.library_path(root, source))
            check = lib.int2_shell_block_selftest
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
            self.assertEqual(count.value, 17600)
            self.assertEqual(failed.value, 0)
            self.assertLessEqual(error.value, 2e-10)


if __name__ == "__main__":
    unittest.main()
