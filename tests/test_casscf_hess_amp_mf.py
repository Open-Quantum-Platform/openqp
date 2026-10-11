"""Equivalence pin: dense vs matrix-free (MF) CASSCF Hessian amplitude kernel.

Both paths compute the same

    amp[k] = _apply_active_operator(f_der[k], g_der[k], stack, wmat) @ vecs

but the MF path replaces the materialised [nact,nact,ndet,ndet] excitation
stack with an on-the-fly sparse walk.  The two must agree to machine precision.
"""
import numpy as np
import pytest
import sys

from oqp.library.casscf_hessian import (
    _apply_active_operator,
    _excitation_matrices_python,
    _hess_amp_backend,
    _lib_hess_amp,
)
from oqp.library.fci import _as_f64c, _as_i64c

def _mf_wmat(nact, dets, ci):
    backend = _hess_amp_backend()
    assert backend is not None, "Native Hessian backend is required"
    lib, ffi = backend
    assert hasattr(lib, 'casscf_exc_stack_apply_wmat')
    buffers = [_as_i64c(dets), _as_i64c(np.sort(dets)),
               _as_i64c(np.argsort(dets).astype(np.int64))]
    ci = _as_f64c(ci)
    out = np.full((nact, nact, len(dets)), np.nan)
    status = lib.casscf_exc_stack_apply_wmat(
        nact, len(dets), *(ffi.cast('int64_t *', x.ctypes.data) for x in buffers),
        ffi.cast('double *', ci.ctypes.data), ffi.cast('double *', out.ctypes.data))
    assert status == 0
    return out


@pytest.mark.parametrize('nact,nalpha,nbeta',
                         [(3, 1, 1), (4, 2, 2), (5, 3, 2), (6, 3, 3)])
def test_matrix_free_w_matches_excitation_operator(nact, nalpha, nbeta):
    dets, stack = _excitation_matrices_python(nact, nalpha, nbeta)
    rng = np.random.default_rng(470 + nact)
    # Repeat with distinct vectors in one process to expose stale output.
    for _ in range(2):
        ci = rng.normal(size=len(dets))
        ci /= np.linalg.norm(ci)
        actual = _mf_wmat(nact, dets, ci)
        expected = np.tensordot(stack, ci, axes=([3], [0]))
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)


def _case(nact, nalpha, nbeta, npar, seed):
    dets, stack = _excitation_matrices_python(nact, nalpha, nbeta)
    ndet = len(dets)
    rng = np.random.default_rng(seed)
    f_der = rng.standard_normal((npar, nact, nact))
    g_der = rng.standard_normal((npar,) + (nact,) * 4)
    ci = rng.standard_normal(ndet)
    ci /= np.linalg.norm(ci)
    wmat = _mf_wmat(nact, dets, ci)
    reference_w = np.tensordot(stack, ci, axes=([3], [0]))
    np.testing.assert_allclose(wmat, reference_w, rtol=0, atol=1e-14)
    vecs = np.linalg.qr(rng.standard_normal((ndet, ndet)))[0]
    return dets, stack, f_der, g_der, wmat, vecs, ndet, ci


def _mf_hess_amp(dets, f_der, g_der, wmat, vecs):
    """Call casscf_hess_amp_mf via the CFFI backend directly."""
    backend = _hess_amp_backend()
    assert backend is not None, "Native Hessian backend is required"
    lib, ffi = backend
    assert hasattr(lib, "casscf_hess_amp_mf")
    nact = int(f_der.shape[1])
    ndet = int(vecs.shape[0])
    npar = int(f_der.shape[0])

    det_arr = _as_i64c(np.ascontiguousarray(dets, dtype=np.int64))
    skeys  = _as_i64c(np.ascontiguousarray(np.sort(dets), dtype=np.int64))
    sperm  = _as_i64c(np.ascontiguousarray(np.argsort(dets).astype(np.int64), dtype=np.int64))
    fd = _as_f64c(f_der)
    gd = _as_f64c(g_der)
    wm = _as_f64c(wmat)
    vc = _as_f64c(vecs)

    amp = np.zeros((npar, ndet), dtype=np.float64)
    status = lib.casscf_hess_amp_mf(
        nact, ndet, npar,
        ffi.cast("int64_t *", det_arr.ctypes.data),
        ffi.cast("int64_t *", skeys.ctypes.data),
        ffi.cast("int64_t *", sperm.ctypes.data),
        ffi.cast("double *", fd.ctypes.data),
        ffi.cast("double *", gd.ctypes.data),
        ffi.cast("double *", wm.ctypes.data),
        ffi.cast("double *", vc.ctypes.data),
        ffi.cast("double *", amp.ctypes.data),
    )
    assert status == 0, f"casscf_hess_amp_mf returned {status}"
    return amp


@pytest.mark.parametrize(
    "nact,nalpha,nbeta,npar",
    [
        (4, 2, 2, 15),
        (5, 3, 2, 23),
        (6, 3, 3, 31),
        (4, 2, 1, 9),   # unequal alpha/beta
        (3, 1, 1, 5),   # smallest useful active space
    ],
)
def test_mf_amp_matches_dense(nact, nalpha, nbeta, npar):
    dets, stack, f_der, g_der, wmat, vecs, ndet, ci = _case(
        nact, nalpha, nbeta, npar, seed=nact * 101 + npar
    )

    # Dense path
    dense = _lib_hess_amp(stack, f_der, g_der, wmat, vecs)
    assert dense is not None
    assert dense.shape == (npar, ndet)

    # MF path
    mf = _mf_hess_amp(dets, f_der, g_der, wmat, vecs)
    assert mf.shape == (npar, ndet)

    # Must agree to machine precision
    # Independent Python operator, using E_tu @ CI rather than the native W.
    reference_w = np.tensordot(stack, ci, axes=([3], [0]))
    reference = np.array([
        _apply_active_operator(f, g, stack, reference_w) @ vecs
        for f, g in zip(f_der, g_der)
    ])
    scale = max(1.0, np.abs(reference).max())
    np.testing.assert_allclose(mf, dense, rtol=0, atol=1e-13 * scale)
    np.testing.assert_allclose(mf, reference, rtol=0, atol=1e-13 * scale)
    # Scheduling must not change the result in a process with repeated calls.
    np.testing.assert_array_equal(mf, _mf_hess_amp(dets, f_der, g_der, wmat, vecs))

@pytest.mark.parametrize('nact,ndet', [(0, 1), (32, 1), (3, 0), (3, -1)])
def test_invalid_w_dimensions_return_status_without_writing(nact, ndet):
    lib, ffi = _hess_amp_backend()
    keys = np.zeros(1, dtype=np.int64)
    values = np.full(32 * 32, 123.0)
    kp = ffi.cast('int64_t *', keys.ctypes.data)
    vp = ffi.cast('double *', values.ctypes.data)
    status = lib.casscf_exc_stack_apply_wmat(nact, ndet, kp, kp, kp, vp, vp)
    assert status == -1
    np.testing.assert_array_equal(values, np.full_like(values, 123.0))
    status = lib.casscf_hess_amp_mf(nact, ndet, 1, kp, kp, kp, vp, vp, vp, vp, vp)
    assert status == -1
    np.testing.assert_array_equal(values, np.full_like(values, 123.0))


@pytest.mark.parametrize('position', [-1, 1])
def test_one_based_or_negative_permutation_is_rejected(position):
    lib, ffi = _hess_amp_backend()
    keys = np.zeros(1, dtype=np.int64)
    perm = np.array([position], dtype=np.int64)
    values = np.full(16, 123.0)
    kp = ffi.cast('int64_t *', keys.ctypes.data)
    pp = ffi.cast('int64_t *', perm.ctypes.data)
    vp = ffi.cast('double *', values.ctypes.data)
    assert lib.casscf_exc_stack_apply_wmat(1, 1, kp, kp, pp, vp, vp) == -1
    assert lib.casscf_hess_amp_mf(1, 1, 1, kp, kp, pp, vp, vp, vp, vp, vp) == -1
    np.testing.assert_array_equal(values, np.full_like(values, 123.0))


@pytest.mark.parametrize('nact,ndet,npar', [
    (1, 1, 0), (1, 1, -1),
    (1, np.iinfo(np.int64).max, 2),
    (31, np.iinfo(np.int64).max, 1),
])
def test_amplitude_extent_rejected_before_reading_arrays(nact, ndet, npar):
    """A huge declared extent must fail before even reading the permutation."""
    lib, ffi = _hess_amp_backend()
    keys = np.zeros(1, dtype=np.int64)
    values = np.full(1, 123.0)
    kp = ffi.cast('int64_t *', keys.ctypes.data)
    vp = ffi.cast('double *', values.ctypes.data)
    assert lib.casscf_hess_amp_mf(
        nact, ndet, npar, kp, kp, kp, vp, vp, vp, vp, vp) == -1
    np.testing.assert_array_equal(values, [123.0])


def test_w_extent_rejected_before_reading_arrays():
    lib, ffi = _hess_amp_backend()
    keys = np.zeros(1, dtype=np.int64)
    values = np.full(1, 123.0)
    kp = ffi.cast('int64_t *', keys.ctypes.data)
    vp = ffi.cast('double *', values.ctypes.data)
    assert lib.casscf_exc_stack_apply_wmat(
        31, np.iinfo(np.int64).max, kp, kp, kp, vp, vp) == -1
    np.testing.assert_array_equal(values, [123.0])


@pytest.mark.skipif(sys.platform != 'linux', reason='glibc allocation injection')
def test_allocation_failure_leaves_output_and_allows_recovery(tmp_path):
    """Fail each of the three actual Fortran allocations, then call again."""
    import os
    import subprocess
    from pathlib import Path
    source = tmp_path / 'fault.c'
    source.write_text(r'''
#include <stddef.h>
#include <stdint.h>
extern void *__libc_malloc(size_t);
extern void *__libc_calloc(size_t, size_t);
static size_t target;
static int failures;
void arm(size_t bytes) { target=bytes; failures=0; }
int disarm(void) { target=0; return failures; }
void *malloc(size_t n) {
  if (target && n==target) { target=0; ++failures; return NULL; }
  return __libc_malloc(n);
}
void *calloc(size_t n, size_t width) {
  if (width && n<=SIZE_MAX/width && target && n*width==target) {
    target=0; ++failures; return NULL;
  }
  return __libc_calloc(n,width);
}
''')
    shim = tmp_path / 'fault.so'
    subprocess.run([os.environ.get('CC', 'gcc'), '-shared', '-fPIC',
                    str(source), '-o', str(shim)], check=True)
    child = tmp_path / 'allocation_test.py'
    child.write_text(r'''
import ctypes, os, runpy
import numpy as np
t=runpy.run_path(os.environ['CASSCF_TEST_SOURCE'])
dets,stack,f,g,w,v,nd,ci=t['_case'](3,1,1,7,509)
lib,ffi=t['_hess_amp_backend']()
keys=np.asarray(dets,dtype=np.int64)
sorted_keys=np.sort(keys)
perm=np.argsort(keys).astype(np.int64)
amp=np.full((7,nd),123.0)
buffers=[keys,sorted_keys,perm,f,g,w,v,amp]
pointers=[ffi.cast('int64_t *' if i<3 else 'double *',b.ctypes.data)
          for i,b in enumerate(buffers)]
fault=ctypes.CDLL(os.environ['LD_PRELOAD'])
fault.arm.argtypes=[ctypes.c_size_t]
fault.disarm.restype=ctypes.c_int
reference=t['_mf_hess_amp'](dets,f,g,w,v)
for size in [nd*7*8, nd*3*3*7*8, 3*3*7*8]:
    amp.fill(123.0)
    fault.arm(size)
    status=lib.casscf_hess_amp_mf(3,nd,7,*pointers)
    failures=fault.disarm()
    assert failures==1, (size,failures)
    assert status==-1, (size,status)
    np.testing.assert_array_equal(amp,np.full_like(amp,123.0))
    assert lib.casscf_hess_amp_mf(3,nd,7,*pointers)==0
    np.testing.assert_array_equal(amp,reference)
print('three allocation failures and recoveries passed')
''')
    env = dict(os.environ, LD_PRELOAD=str(shim),
               CASSCF_TEST_SOURCE=str(Path(__file__).resolve()))
    subprocess.run([sys.executable, str(child)], env=env, check=True,
                   capture_output=True, text=True)


def test_w_encoding_limit_and_reordered_determinants():
    det = np.array([(1 << 30) | (1 << 61)], dtype=np.int64)
    actual = _mf_wmat(31, det, np.array([0.7]))
    expected = np.zeros_like(actual)
    expected[30, 30, 0] = 1.4
    np.testing.assert_array_equal(actual, expected)
    # At the bit-width boundary, the sole occupied spatial orbital has two
    # electrons: its one-/two-body coefficients are 2 and 1 respectively.
    f = np.zeros((1, 31, 31))
    g = np.zeros((1, 31, 31, 31, 31))
    f[0, 30, 30] = 0.3
    g[0, 30, 30, 30, 30] = -0.2
    amp = _mf_hess_amp(det, f, g, actual, np.ones((1, 1)))
    np.testing.assert_allclose(amp, [[0.7 * (2 * 0.3 - 0.2)]], rtol=0, atol=1e-14)
    dets, stack = _excitation_matrices_python(4, 2, 1)
    order = np.random.default_rng(15).permutation(len(dets))
    ci = np.random.default_rng(16).normal(size=len(dets))
    expected = np.tensordot(stack[:, :, order][:, :, :, order], ci, axes=([3], [0]))
    np.testing.assert_allclose(_mf_wmat(4, np.asarray(dets)[order], ci),
                               expected, rtol=0, atol=1e-14)


def test_amplitude_chunk_boundary_and_partial_last_chunk():
    """Cross the real 128 MiB x budget with 36 determinants, not a large CAS."""
    nact = 4
    dets, stack = _excitation_matrices_python(nact, 2, 2)
    ndet = len(dets)
    chunk_pairs = 16777216 // (nact * nact * ndet)
    npar = chunk_pairs + 1
    rng = np.random.default_rng(718)
    ci = rng.normal(size=ndet)
    ci /= np.linalg.norm(ci)
    wmat = _mf_wmat(nact, dets, ci)
    f = rng.normal(size=(npar, nact, nact))
    g = rng.normal(size=(npar,) + (nact,) * 4)
    vecs = np.linalg.qr(rng.normal(size=(ndet, ndet)))[0]
    actual = _mf_hess_amp(dets, f, g, wmat, vecs)
    dense = _lib_hess_amp(stack, f, g, wmat, vecs)
    np.testing.assert_allclose(actual, dense, rtol=0, atol=1e-13 * np.abs(dense).max())
    reference_w = np.tensordot(stack, ci, axes=([3], [0]))
    for k in [0, chunk_pairs - 1, chunk_pairs]:
        reference = _apply_active_operator(f[k], g[k], stack, reference_w) @ vecs
        np.testing.assert_allclose(actual[k], reference, rtol=0,
                                   atol=1e-13 * max(1.0, np.abs(reference).max()))


def test_direct_products_are_bit_identical_across_thread_counts():
    import os
    import subprocess
    import sys
    from pathlib import Path
    code = """
import hashlib, runpy, sys
t = runpy.run_path(sys.argv[1])
dets, stack, f, g, w, v, nd, ci = t['_case'](6, 3, 3, 31, 637)
a = t['_mf_hess_amp'](dets, f, g, w, v)
print(hashlib.sha256(w.tobytes() + a.tobytes()).hexdigest())
"""
    digests = []
    for threads in [1, 4]:
        env = dict(os.environ, OMP_NUM_THREADS=str(threads),
                   VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        result = subprocess.run([sys.executable, '-c', code, str(Path(__file__).resolve())],
                                env=env, text=True, capture_output=True, check=True)
        digests.append(result.stdout.strip().splitlines()[-1])
    assert digests[0] == digests[1]
