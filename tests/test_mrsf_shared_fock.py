"""Exercise the actual Fock kernels and reduction, with an independent signed oracle.

Run with the clean installed package's Python, at OMP_NUM_THREADS=1 and 4.
For a compiled MPI package also run with mpiexec -n 2 and
OQP_SHARED_REQUIRE_MPI=1; a serial library fails that mode, rather than skipping.
"""
import ctypes
import os
import unittest


class SharedMRSFFockTests(unittest.TestCase):
    def test_signed_channels_fp32_cam_resize_and_cleanup(self):
        import oqp.runtime as runtime

        library = ctypes.CDLL(runtime.library_path(*runtime.resolve_oqp_root()))
        check = library.int2_mrsf_shared_selftest
        check.argtypes = [ctypes.c_int, ctypes.c_bool, ctypes.c_int, ctypes.c_int,
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        check.restype = None
        use_mpi = os.getenv("OQP_SHARED_REQUIRE_MPI") == "1"
        comm, ranks = 0, 1
        if use_mpi:
            from mpi4py import MPI
            comm, ranks = MPI.COMM_WORLD.py2f(), MPI.COMM_WORLD.size
            self.assertGreaterEqual(ranks, 2, "multirank gate needs at least two ranks")
        threads = int(os.environ["OMP_NUM_THREADS"])
        self.assertIn(threads, (1, 4))
        for _ in range(2):
            error, checks, failures = ctypes.c_double(), ctypes.c_int(), ctypes.c_int()
            check(comm, use_mpi, ranks, threads, ctypes.byref(error),
                  ctypes.byref(checks), ctypes.byref(failures))
            self.assertEqual(checks.value, 24)
            self.assertEqual(failures.value, 0)
            self.assertEqual(error.value, 0.0)


if __name__ == "__main__":
    unittest.main()
