"""Native ECP calls cap and restore BLAS/OpenMP settings on repeated entry.

The fixture supplies only thread controls, including a shared OpenMP setter.
All ECP arithmetic still uses the installed native library's numerical BLAS.
Each width runs in a fresh process so dynamic symbol lookup is independent.
"""
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

CHILD = r'''
import ctypes, json, sys
fixture = ctypes.CDLL(sys.argv[1], mode=ctypes.RTLD_GLOBAL)
try:
    import oqp
    from oqp import runtime
    root, suffix = runtime.resolve_oqp_root()
    native = ctypes.CDLL(str(runtime.library_path(root, suffix)))
    if not hasattr(oqp.lib, 'oqp_ecp_selftest'):
        raise RuntimeError('Missing native ECP selftest')
except Exception:
    if sys.argv[3] == '1':
        raise
    sys.exit(77)
native.oqp_have_openmp.restype = ctypes.c_int
if not native.oqp_have_openmp():
    sys.exit(77)
set_omp = native.omp_set_num_threads
set_omp.argtypes = [ctypes.c_int]
set_omp.restype = None
get_omp = native.omp_get_max_threads
get_omp.restype = ctypes.c_int
native.oqp_blas_thread_count.restype = ctypes.c_int64
fixture.ecp_fixture_callback.argtypes = [ctypes.c_void_p]
fixture.ecp_fixture_callback(ctypes.cast(set_omp, ctypes.c_void_p))
width = int(sys.argv[2])
set_omp(width)
before = (native.oqp_blas_thread_count(), get_omp())
assert before == (4, width), before
previous = None
previous_caps = 0
for repeat in range(2):
    error = oqp.ffi.new('double[10]')
    oqp.lib.oqp_ecp_selftest(error)
    after = (native.oqp_blas_thread_count(), get_omp())
    assert after == before, (repeat, before, after)
    caps = fixture.ecp_fixture_caps()
    calls = fixture.ecp_fixture_calls()
    if width > 1:
        assert caps > previous_caps, (repeat, caps, previous_caps)
        assert calls == 2 * caps, (calls, caps)
    else:
        assert calls == 0, calls
    errors = list(error)
    if previous is not None:
        assert errors == previous, (previous, errors)
    previous, previous_caps = errors, caps
print(json.dumps({'before': before, 'after': after, 'caps': caps, 'calls': calls}))
'''


@unittest.skipIf(os.name == 'nt', 'Windows BLAS control is explicitly a no-op')
class NativeEcpThreadControlTests(unittest.TestCase):
    def test_repeated_serial_and_parallel_calls_restore_thread_settings(self):
        required = os.getenv('OQP_REQUIRE_NATIVE_TESTS') == '1'
        compiler = shutil.which('cc') or shutil.which('gcc')
        if compiler is None or importlib.util.find_spec('oqp') is None:
            if required:
                self.fail('Native OpenQP and a C compiler are required')
            self.skipTest('Native OpenQP and a C compiler are unavailable')
        with tempfile.TemporaryDirectory(prefix='oqp-ecp-threads-') as directory:
            fixture = Path(directory) / ('controls.dylib' if sys.platform == 'darwin' else 'controls.so')
            source = Path(__file__).with_name('ecp_thread_control_fixture.c')
            flags = ['-dynamiclib'] if sys.platform == 'darwin' else ['-shared', '-fPIC']
            built = subprocess.run([compiler, *flags, str(source), '-o', str(fixture)],
                                   capture_output=True, text=True, timeout=30)
            self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
            for width in (1, 4):
                with self.subTest(omp_threads=width):
                    env = os.environ.copy()
                    env['OMP_NUM_THREADS'] = str(width)
                    env.pop('LD_PRELOAD', None)
                    result = subprocess.run([sys.executable, '-c', CHILD, str(fixture), str(width),
                                             '1' if required else '0'],
                                            env=env, capture_output=True, text=True, timeout=180)
                    if result.returncode == 77:
                        self.skipTest('Native OpenMP ECP runtime is unavailable')
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
