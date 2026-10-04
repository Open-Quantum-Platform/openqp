"""Native response allocation reuse, inactive images and cleanup regression."""
import ctypes
import unittest

class ResponseWorkspaceTest(unittest.TestCase):
    def test_response_workspace(self):
        import oqp
        from oqp.runtime import library_path
        lib = ctypes.CDLL(str(library_path(oqp.oqp_root, oqp.suffix)))
        func = lib.response_workspace_selftest
        func.argtypes = [ctypes.POINTER(ctypes.c_int)]
        func.restype = None
        failed = ctypes.c_int()
        func(ctypes.byref(failed))
        self.assertEqual(failed.value, 0)
