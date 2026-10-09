"""Equivalence pin: dense vs matrix-free (MF) CASSCF Hessian amplitude kernel.

Both paths compute the same

    amp[k] = _apply_active_operator(f_der[k], g_der[k], stack, wmat) @ vecs

but the MF path replaces the materialised [nact,nact,ndet,ndet] excitation
stack with an on-the-fly sparse walk.  The two must agree to machine precision.
"""
import numpy as np
import pytest

from oqp.library.casscf_hessian import (
    _excitation_matrices,
    _hess_amp_backend,
    _lib_hess_amp,
)
from oqp.library.fci import _as_f64c, _as_i64c

pytestmark = pytest.mark.skipif(
    _hess_amp_backend() is None,
    reason="liboqp casscf_hess_amp unavailable",
)


def _check_mf_backend():
    """True when the MF engine symbol is present."""
    backend = _hess_amp_backend()
    if backend is None:
        return False
    lib, _ffi = backend
    return hasattr(lib, "casscf_hess_amp_mf")


def _case(nact, nalpha, nbeta, npar, seed):
    dets, stack = _excitation_matrices(nact, nalpha, nbeta)
    ndet = len(dets)
    rng = np.random.default_rng(seed)
    f_der = rng.standard_normal((npar, nact, nact))
    g_der = rng.standard_normal((npar,) + (nact,) * 4)
    wmat = rng.standard_normal((nact, nact, ndet))
    vecs = np.linalg.qr(rng.standard_normal((ndet, ndet)))[0]
    return dets, stack, f_der, g_der, wmat, vecs, ndet


def _mf_hess_amp(dets, f_der, g_der, wmat, vecs):
    """Call casscf_hess_amp_mf via the CFFI backend directly."""
    lib, ffi = _hess_amp_backend()
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
    lib.casscf_hess_amp_mf(
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
    return amp


@pytest.mark.skipif(not _check_mf_backend(), reason="MF engine casscf_hess_amp_mf unavailable")
@pytest.mark.parametrize(
    "nact,nalpha,nbeta,npar",
    [
        (4, 2, 2, 15),
        (5, 3, 2, 23),
        (6, 3, 3, 31),
        (8, 4, 4, 45),
        (4, 2, 1, 9),   # unequal alpha/beta
        (3, 1, 1, 5),   # smallest useful active space
    ],
)
def test_mf_amp_matches_dense(nact, nalpha, nbeta, npar):
    dets, stack, f_der, g_der, wmat, vecs, ndet = _case(
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
    max_abs = max(np.abs(dense).max(), np.abs(mf).max())
    tol = 1e-11 * max_abs if max_abs > 1e-100 else 1e-14
    np.testing.assert_allclose(mf, dense, rtol=0, atol=tol)
