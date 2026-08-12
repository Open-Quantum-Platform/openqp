"""Equivalence pin for the Fortran CI-relaxation amplitude engine
(casscf_hess_kernel.F90).

The engine replaces the per-pair loop in ``casscf_hessian.py``

    amp[k] = _apply_active_operator(f_der[k], g_der[k], stack, wmat) @ vecs

with the same contractions batched over the rotation pair index, so the two
are algebraically identical and are held to round-off here.  The Python loop
stays in the tree as this pin and as the fallback.
"""
import numpy as np
import pytest

from oqp.library.casscf_hessian import (
    _apply_active_operator,
    _excitation_matrices,
    _hess_amp_backend,
    _lib_hess_amp,
)

pytestmark = pytest.mark.skipif(
    _hess_amp_backend() is None,
    reason="liboqp casscf_hess_amp unavailable",
)


def _case(nact, nalpha, nbeta, npar, seed):
    dets, stack = _excitation_matrices(nact, nalpha, nbeta)
    ndet = len(dets)
    rng = np.random.default_rng(seed)
    f_der = rng.standard_normal((npar, nact, nact))
    g_der = rng.standard_normal((npar,) + (nact,) * 4)
    wmat = rng.standard_normal((nact, nact, ndet))
    # an orthonormal eigenbasis, as _symmetric_eigh would return
    vecs = np.linalg.qr(rng.standard_normal((ndet, ndet)))[0]
    return stack, f_der, g_der, wmat, vecs, ndet


def _reference(stack, f_der, g_der, wmat, vecs, npar, ndet):
    amp = np.empty((npar, ndet))
    for k in range(npar):
        amp[k] = _apply_active_operator(f_der[k], g_der[k], stack, wmat) @ vecs
    return amp


@pytest.mark.parametrize(
    "nact,nalpha,nbeta,npar",
    [
        (4, 2, 2, 15),
        (5, 3, 2, 23),
        (6, 3, 3, 31),
        (4, 2, 1, 9),   # unequal alpha/beta
        (3, 1, 1, 5),   # smallest useful active space
        (4, 2, 2, 1),   # single pair: the degenerate chunk shape
    ],
)
def test_fortran_amp_matches_python_loop(nact, nalpha, nbeta, npar):
    stack, f_der, g_der, wmat, vecs, ndet = _case(
        nact, nalpha, nbeta, npar, seed=nact * 101 + npar)
    ref = _reference(stack, f_der, g_der, wmat, vecs, npar, ndet)

    got = _lib_hess_amp(stack, f_der, g_der, wmat, vecs)
    assert got is not None
    assert got.shape == (npar, ndet)
    np.testing.assert_allclose(got, ref, rtol=0, atol=1e-11 * np.abs(ref).max())


def test_amp_is_linear_in_the_derivative_integrals():
    """The engine batches pairs into shared GEMMs; each pair must stay independent.

    A pair-index slip (wrong chunk offset, wrong leading dimension) would still
    produce a plausible matrix, so this checks the property that pins the
    batching: row k must depend only on (f_der[k], g_der[k]).  Building the
    amplitudes one pair at a time has to reproduce the batched result exactly.
    """
    nact, nalpha, nbeta, npar = 5, 3, 2, 12
    stack, f_der, g_der, wmat, vecs, ndet = _case(nact, nalpha, nbeta, npar, seed=99)

    batched = _lib_hess_amp(stack, f_der, g_der, wmat, vecs)
    for k in range(npar):
        single = _lib_hess_amp(stack, f_der[k:k + 1], g_der[k:k + 1], wmat, vecs)
        np.testing.assert_allclose(
            single[0], batched[k], rtol=0, atol=1e-11 * np.abs(batched[k]).max())
