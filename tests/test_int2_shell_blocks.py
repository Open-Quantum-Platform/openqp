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
        # The last 16 checks pin zero-scaled J/K digestion (RHF and UROHF;
        # legacy update and forced shell consume_shell; replicated and atomic
        # accumulators; one and four images) by linearity against the full path.
        for layout in (None, "legacy", "shell", "unknown", "x" * 17):
            with self.subTest(layout=layout), patch.dict(os.environ):
                if layout is None:
                    os.environ.pop("OQP_INT2_LAYOUT", None)
                else:
                    os.environ["OQP_INT2_LAYOUT"] = layout
                self._check("int2_shell_block_selftest", 44016)

    def test_zero_scale_checks_under_forced_accumulator_mode(self):
        # parallel_start sizes the Fock accumulator from OQP_FOCK_ATOMIC; the
        # zero-scale checks then select each mode themselves and must resize
        # the accumulator to match (one shared copy vs one per image).
        for atomic in ("0", "1"):
            with self.subTest(OQP_FOCK_ATOMIC=atomic), patch.dict(os.environ):
                os.environ.pop("OQP_INT2_LAYOUT", None)
                os.environ["OQP_FOCK_ATOMIC"] = atomic
                self._check("int2_shell_block_selftest", 44016)

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
