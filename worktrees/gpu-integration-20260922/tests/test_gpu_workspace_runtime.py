"""METC-C1a workspace runtime tests.

These compile ``source/gpu_workspace_runtime.c`` in its host-fallback mode (no
``OQP_CUDA_ENABLE``) and exercise the real C ABI via ctypes: a contiguous arena
allocator with handle-based pointer lookup by byte offset.  The arena sizes and
the d3/f3 byte offsets are driven by a real METC workspace manifest from the
Python ``GpuWorkspaceManager`` so the source runtime and the Python control
plane agree end to end.

The METC contraction wrapper is deliberately NOT exercised here -- C1a only
proves the runtime allocation/ptr/release bridge.  If no C compiler is
available the whole module is skipped.
"""

import ctypes
import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNTIME_C = ROOT / "src" / "metc" / "gpu_workspace_runtime.c"


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _find_cc():
    for cc in ("cc", "gcc", "clang"):
        path = shutil.which(cc)
        if path:
            return path
    return None


def _compile_runtime():
    """Compile the runtime as a host-fallback shared library; return its path."""

    cc = _find_cc()
    if cc is None:
        return None
    tmpdir = Path(tempfile.mkdtemp(prefix="oqp-ws-rt-"))
    lib = tmpdir / "libgpu_ws_runtime.so"
    proc = subprocess.run(
        [cc, "-shared", "-fPIC", "-O0", "-o", str(lib), str(RUNTIME_C)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise AssertionError(
            "failed to compile gpu_workspace_runtime.c:\n" + proc.stderr
        )
    return lib


def _bind(lib_path):
    lib = ctypes.CDLL(str(lib_path))
    lib.oqp_gpu_ws_acquire.restype = ctypes.c_int
    lib.oqp_gpu_ws_acquire.argtypes = [ctypes.c_int, ctypes.c_int64]
    lib.oqp_gpu_ws_validate.restype = ctypes.c_int
    lib.oqp_gpu_ws_validate.argtypes = [ctypes.c_int, ctypes.c_int64, ctypes.c_int64]
    lib.oqp_gpu_ws_ptr.restype = ctypes.c_void_p
    lib.oqp_gpu_ws_ptr.argtypes = [ctypes.c_int, ctypes.c_int64]
    lib.oqp_gpu_ws_total_bytes.restype = ctypes.c_int64
    lib.oqp_gpu_ws_total_bytes.argtypes = [ctypes.c_int]
    lib.oqp_gpu_ws_release.restype = ctypes.c_int
    lib.oqp_gpu_ws_release.argtypes = [ctypes.c_int]
    return lib


@unittest.skipIf(_find_cc() is None, "no C compiler available")
class WorkspaceRuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib = _bind(_compile_runtime())
        cls.ws = load_module("gpu_workspace_runtime_ws", "python/openqp_gpu/gpu_workspace.py")
        cls.buffers = load_module(
            "gpu_metc_buffers_runtime", "python/openqp_gpu/gpu_metc_buffers.py"
        )

    def _metc_manifest(self, *, nthreads, policy=None):
        plan = self.buffers.PersistentMetcBufferPlan.from_problem(
            nbf=7, nf=5, nmatrix=3, max_integrals=11
        )
        manager = self.ws.GpuWorkspaceManager()
        if policy is None:
            alloc = manager.from_metc_plan(plan, nthreads=nthreads)
        else:
            alloc = manager.from_metc_plan(plan, f3_policy=policy, nthreads=nthreads)
        return plan, self.ws.workspace_manifest(alloc)

    def _row(self, manifest, name):
        return next(r for r in manifest["rows"] if r["name"] == name)

    def test_acquire_returns_positive_handle_and_total_bytes(self):
        _, man = self._metc_manifest(nthreads=1)
        handle = self.lib.oqp_gpu_ws_acquire(man["target_code"], man["total_bytes"])
        try:
            self.assertGreater(handle, 0)
            self.assertEqual(self.lib.oqp_gpu_ws_total_bytes(handle), man["total_bytes"])
        finally:
            self.lib.oqp_gpu_ws_release(handle)

    def test_acquire_rejects_bad_args(self):
        self.assertEqual(self.lib.oqp_gpu_ws_acquire(0, 0), 0)      # zero bytes
        self.assertEqual(self.lib.oqp_gpu_ws_acquire(0, -8), 0)     # negative bytes
        self.assertEqual(self.lib.oqp_gpu_ws_acquire(9, 64), 0)     # unknown target

    def test_ptr_offset_zero_is_non_null(self):
        _, man = self._metc_manifest(nthreads=1)
        handle = self.lib.oqp_gpu_ws_acquire(man["target_code"], man["total_bytes"])
        try:
            self.assertIsNotNone(self.lib.oqp_gpu_ws_ptr(handle, 0))
        finally:
            self.lib.oqp_gpu_ws_release(handle)

    def test_d3_and_f3_pointers_are_distinct(self):
        _, man = self._metc_manifest(nthreads=2)
        d3 = self._row(man, "density")   # METC d3 input
        f3 = self._row(man, "fock")      # METC f3 accumulator
        handle = self.lib.oqp_gpu_ws_acquire(man["target_code"], man["total_bytes"])
        try:
            p_d3 = self.lib.oqp_gpu_ws_ptr(handle, d3["offset"])
            p_f3 = self.lib.oqp_gpu_ws_ptr(handle, f3["offset"])
            self.assertIsNotNone(p_d3)
            self.assertIsNotNone(p_f3)
            self.assertNotEqual(p_d3, p_f3)
            # The pointer delta equals the manifest offset delta.
            self.assertEqual(p_f3 - p_d3, f3["offset"] - d3["offset"])
            # Each row validates within the arena.
            self.assertEqual(self.lib.oqp_gpu_ws_validate(handle, d3["offset"], d3["bytes"]), 0)
            self.assertEqual(self.lib.oqp_gpu_ws_validate(handle, f3["offset"], f3["bytes"]), 0)
        finally:
            self.lib.oqp_gpu_ws_release(handle)

    def test_invalid_offsets_and_rows_are_rejected(self):
        _, man = self._metc_manifest(nthreads=1)
        total = man["total_bytes"]
        handle = self.lib.oqp_gpu_ws_acquire(man["target_code"], total)
        try:
            # Offset at/after end -> null pointer.
            self.assertIsNone(self.lib.oqp_gpu_ws_ptr(handle, total))
            self.assertIsNone(self.lib.oqp_gpu_ws_ptr(handle, total + 8))
            # A row whose offset+bytes exceeds the arena fails validation.
            self.assertNotEqual(self.lib.oqp_gpu_ws_validate(handle, total - 8, 16), 0)
            self.assertNotEqual(self.lib.oqp_gpu_ws_validate(handle, -1, 8), 0)
        finally:
            self.lib.oqp_gpu_ws_release(handle)

    def test_invalid_handle_is_rejected(self):
        self.assertIsNone(self.lib.oqp_gpu_ws_ptr(999, 0))
        self.assertEqual(self.lib.oqp_gpu_ws_total_bytes(999), -1)
        self.assertNotEqual(self.lib.oqp_gpu_ws_validate(999, 0, 8), 0)

    def test_release_invalidates_handle_and_double_release_is_safe(self):
        _, man = self._metc_manifest(nthreads=1)
        handle = self.lib.oqp_gpu_ws_acquire(man["target_code"], man["total_bytes"])
        self.assertEqual(self.lib.oqp_gpu_ws_release(handle), 0)
        # Released handle returns no pointer and no size.
        self.assertIsNone(self.lib.oqp_gpu_ws_ptr(handle, 0))
        self.assertEqual(self.lib.oqp_gpu_ws_total_bytes(handle), -1)
        # Double release is safe (nonzero error, no crash).
        self.assertNotEqual(self.lib.oqp_gpu_ws_release(handle), 0)

    def test_per_thread_vs_single_atomic_arena_sizes(self):
        plan_pt, man_pt = self._metc_manifest(
            nthreads=8, policy=self.ws.F3AccumulatorPolicy.PER_THREAD
        )
        _, man_sa = self._metc_manifest(
            nthreads=8, policy=self.ws.F3AccumulatorPolicy.SINGLE_ATOMIC
        )
        fock_pt = self._row(man_pt, "fock")["bytes"]
        fock_sa = self._row(man_sa, "fock")["bytes"]
        base = plan_pt.bytes_for("fock")
        self.assertEqual(fock_pt, base * 8)   # PER_THREAD: nthreads x base
        self.assertEqual(fock_sa, base)       # SINGLE_ATOMIC: one base accumulator
        # The runtime honors the differing arena totals.
        h_pt = self.lib.oqp_gpu_ws_acquire(man_pt["target_code"], man_pt["total_bytes"])
        h_sa = self.lib.oqp_gpu_ws_acquire(man_sa["target_code"], man_sa["total_bytes"])
        try:
            self.assertEqual(self.lib.oqp_gpu_ws_total_bytes(h_pt), man_pt["total_bytes"])
            self.assertEqual(self.lib.oqp_gpu_ws_total_bytes(h_sa), man_sa["total_bytes"])
            self.assertGreater(man_pt["total_bytes"], man_sa["total_bytes"])
        finally:
            self.lib.oqp_gpu_ws_release(h_pt)
            self.lib.oqp_gpu_ws_release(h_sa)

    def test_validate_rejects_int64_overflow(self):
        h = self.lib.oqp_gpu_ws_acquire(0, 64)
        try:
            self.assertNotEqual(self.lib.oqp_gpu_ws_validate(h, 2**63-4, 8), 0)
        finally:
            self.lib.oqp_gpu_ws_release(h)

    def test_session_rejects_overflow_shape(self):
        fn = self.lib.oqp_gpu_metc_session_begin
        fn.argtypes = [ctypes.c_int] * 5
        fn.restype = ctypes.c_int
        self.assertEqual(fn(2**30, 11, 100, 1, 10), 0)

    def test_host_only_contraction_reports_unavailable(self):
        begin = self.lib.oqp_gpu_metc_session_begin
        begin.argtypes = [ctypes.c_int] * 5
        begin.restype = ctypes.c_int
        session = begin(1, 7, 1, 1, 1)
        self.assertGreater(session, 0)
        try:
            fn = self.lib.oqp_gpu_metc_session_contract
            fn.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                           ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int,
                           ctypes.c_double, ctypes.c_double, ctypes.c_int]
            fn.restype = ctypes.c_int
            self.assertEqual(fn(session, 0, (ctypes.c_int*4)(1,1,1,1),
                                (ctypes.c_double*1)(.5), 1, 1, 1., 1., 0), 9)
        finally:
            self.lib.oqp_gpu_metc_session_end(session)


if __name__ == "__main__":
    unittest.main()
