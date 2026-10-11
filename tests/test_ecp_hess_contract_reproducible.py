"""The threaded ECP Hessian contraction does not depend on thread scheduling.

``ecp_hess_contract`` (source/ecp.F90) splits shell pairs over OpenMP
threads.  It used ``schedule(dynamic)`` with ``reduction(+:hloc)``: which
thread summed which pairs, and the order in which OpenMP combined the private
copies, followed the order in which threads happened to run.  err[9] of
``oqp_ecp_selftest`` therefore changed in the last bits with scheduling: two
calls in one GitHub Actions process gave 5.55e-17 and 8.33e-17.  With four
threads on TPX (GCC 13) the old code gives 0x1p-54 when the threads share one
CPU and 0x1.8p-54 when they have two, every time; so pinning reproduces what
CI hit by chance.

For a fixed thread count every call must now give the same bits, wherever the
threads run.  On Linux each child is pinned to one CPU and to two CPUs before
OpenMP starts; elsewhere only repeated unpinned calls are compared.
"""
import importlib.util
import json
import os
import subprocess
import sys
import unittest

REPEATS = 10
THREADS = 4

CHILD = r'''
import ctypes, json, os, sys
repeats, threads, required, cpus = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3] == '1', sys.argv[4]
if cpus != 'none':
    # before OpenMP creates its threads, so that they all inherit the mask
    os.sched_setaffinity(0, [int(c) for c in cpus.split(',')])
try:
    import oqp
    from oqp import runtime
    root, suffix = runtime.resolve_oqp_root()
    native = ctypes.CDLL(str(runtime.library_path(root, suffix)))
    if not hasattr(oqp.lib, 'oqp_ecp_selftest'):
        raise RuntimeError('Missing native ECP selftest')
except Exception:
    if required:
        raise
    sys.exit(77)
native.oqp_have_openmp.restype = ctypes.c_int
if not native.oqp_have_openmp():
    if required:
        raise RuntimeError('OpenQP was built without OpenMP')
    sys.exit(77)
native.omp_get_max_threads.restype = ctypes.c_int
assert native.omp_get_max_threads() == threads, native.omp_get_max_threads()
vectors = set()
for _ in range(repeats):
    error = oqp.ffi.new('double[10]')
    oqp.lib.oqp_ecp_selftest(error)
    vectors.add(tuple(error))
print('RESULT ' + json.dumps(sorted(vectors)))
'''


class EcpHessContractReproducibleTests(unittest.TestCase):
    def _run(self, cpus, required):
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(THREADS)
        env['OMP_DYNAMIC'] = 'false'
        for name in ('OMP_PROC_BIND', 'OMP_PLACES', 'GOMP_CPU_AFFINITY', 'KMP_AFFINITY'):
            env.pop(name, None)
        result = subprocess.run([sys.executable, '-c', CHILD, str(REPEATS), str(THREADS),
                                 '1' if required else '0', cpus],
                                env=env, capture_output=True, text=True, timeout=600)
        if result.returncode == 77:
            self.skipTest('Native OpenMP ECP runtime is unavailable')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        line = [l for l in result.stdout.splitlines() if l.startswith('RESULT ')][-1]
        return [tuple(v) for v in json.loads(line[len('RESULT '):])]

    def test_four_thread_selftest_is_bitwise_identical_however_scheduled(self):
        required = os.getenv('OQP_REQUIRE_NATIVE_TESTS') == '1'
        if importlib.util.find_spec('oqp') is None:
            if required:
                self.fail('Native OpenQP is required')
            self.skipTest('Native OpenQP is unavailable')
        if hasattr(os, 'sched_getaffinity') and len(os.sched_getaffinity(0)) >= 2:
            mine = sorted(os.sched_getaffinity(0))
            placements = [str(mine[0]), f'{mine[0]},{mine[1]}']
        else:
            placements = ['none']
        vectors = {}
        for cpus in placements:
            with self.subTest(cpus=cpus):
                vectors[cpus] = self._run(cpus, required)
                self.assertEqual(len(vectors[cpus]), 1,
                                 f'{REPEATS} calls on CPUs {cpus} gave {vectors[cpus]}')
                # the contraction still agrees with the stored second derivatives
                self.assertLess(vectors[cpus][0][9], 1e-13)
        distinct = {v for vs in vectors.values() for v in vs}
        self.assertEqual(len(distinct), 1,
                         f'{THREADS}-thread error vectors depend on CPU placement: {vectors}')


if __name__ == '__main__':
    unittest.main()
