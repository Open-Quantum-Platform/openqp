"""Equivalence pin: dense vs matrix-free (MF) CASSCF Hessian amplitude kernel.

Both paths compute the same

    amp[k] = _apply_active_operator(f_der[k], g_der[k], stack, wmat) @ vecs

but the MF path replaces the materialised [nact,nact,ndet,ndet] excitation
stack with an on-the-fly sparse walk.  The two must agree to machine precision.
"""
import numpy as np
import pytest

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
