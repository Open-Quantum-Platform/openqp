"""Equivalence pin for the Fortran CI-relaxation accumulation.

``casscf_hess_relax`` builds the response weights ``2 coef_j / (E_I - eps_j)``
and accumulates ``sum_j factor_j amp[:,j] amp[:,j]^T``.  It replaces an
ndet-long Python loop plus a masked NumPy product, so the pin holds it to
round-off against that loop -- and, just as importantly, pins the two
degeneracy branches, where the right behaviour is to refuse rather than to
return a plausible number.
"""
import numpy as np
import pytest

from oqp.library.casscf_hessian import _COUPLING_NOISE, _hess_amp_backend, _lib_hess_relax

pytestmark = pytest.mark.skipif(
    _hess_amp_backend() is None, reason="liboqp casscf_hess_relax unavailable")


def _python_reference(npar, ndet, iavg, ovl, weights, eps, e_i, amp, hess, tol):
    w_i = float(weights[iavg])
    factors = np.zeros(ndet)
    for j in range(ndet):
        if ovl[iavg, j] > 0.5:
            continue
        m = int(np.argmax(ovl[:, j]))
        coef = 0.5 * (w_i - float(weights[m])) if ovl[m, j] > 0.5 else w_i
        if coef == 0.0:
            continue
        denom = e_i - float(eps[j])
        if abs(denom) < tol:
            if float(np.max(np.abs(amp[:, j]))) * abs(coef) < _COUPLING_NOISE:
                continue
            raise ValueError("degenerate")
        factors[j] = 2.0 * coef / denom
    nz = factors != 0.0
    if np.any(nz):
        hess += (amp[:, nz] * factors[nz]) @ amp[:, nz].T
    return hess


def _case(npar, ndet, navg, seed):
    rng = np.random.default_rng(seed)
    # Overlaps must stay below the 0.5 "this eigenstate IS the averaged root"
    # threshold everywhere except the one eigenstate each averaged root owns --
    # otherwise every state looks like a reference and the branches under test
    # are never reached.
    ovl = rng.random((navg, ndet)) * 0.4
    ovl[:, :navg] = 0.0
    for a in range(navg):
        ovl[a, a] = 0.9          # each averaged root owns one eigenstate
    weights = np.full(navg, 1.0 / navg)
    eps = np.sort(rng.standard_normal(ndet))
    amp = rng.standard_normal((npar, ndet))
    return ovl, weights, eps, amp


@pytest.mark.parametrize("npar,ndet,navg", [(12, 40, 1), (20, 60, 2), (31, 120, 3)])
def test_relax_matches_the_python_loop(npar, ndet, navg):
    ovl, weights, eps, amp = _case(npar, ndet, navg, seed=npar)
    e_i = float(eps[0]) + 1.7
    for iavg in range(navg):
        got = np.zeros((npar, npar))
        ref = np.zeros((npar, npar))
        assert _lib_hess_relax(npar, ndet, iavg, ovl, weights, eps, e_i, amp, got, 1e-8)
        _python_reference(npar, ndet, iavg, ovl, weights, eps, e_i, amp, ref, 1e-8)
        np.testing.assert_allclose(got, ref, rtol=0, atol=1e-11 * np.abs(ref).max())


def test_equal_weight_state_average_cancels_within_the_set():
    """The (w_I - w_J)/2 branch must vanish exactly for equal weights.

    This is the invariance of equal-weight state averaging; if the engine got
    that coefficient wrong it would still produce a symmetric, plausible
    matrix, so it is pinned directly: an eigenstate owned by another averaged
    root must contribute nothing.
    """
    npar, ndet, navg = 10, 30, 2
    ovl, weights, eps, amp = _case(npar, ndet, navg, seed=3)
    e_i = float(eps[0]) + 2.0
    got = np.zeros((npar, npar))
    _lib_hess_relax(npar, ndet, 0, ovl, weights, eps, e_i, amp, got, 1e-8)

    # zero out the partner's amplitudes: with equal weights it contributed 0,
    # so the result must be unchanged.
    amp2 = amp.copy()
    amp2[:, 1] = 0.0
    got2 = np.zeros((npar, npar))
    _lib_hess_relax(npar, ndet, 0, ovl, weights, eps, e_i, amp2, got2, 1e-8)
    np.testing.assert_allclose(got2, got, rtol=0, atol=1e-11 * np.abs(got).max())


def test_degenerate_root_with_real_coupling_is_refused():
    npar, ndet, navg = 8, 20, 1
    ovl, weights, eps, amp = _case(npar, ndet, navg, seed=11)
    e_i = float(eps[5])                     # exactly degenerate with eigenstate 5
    hess = np.zeros((npar, npar))
    with pytest.raises(ValueError, match="root degeneracy"):
        _lib_hess_relax(npar, ndet, 0, ovl, weights, eps, e_i, amp, hess, 1e-8)


def test_degenerate_but_uncoupled_partner_is_skipped():
    """Spin/symmetry-forbidden partners sit at zero coupling and must not raise."""
    npar, ndet, navg = 8, 20, 1
    ovl, weights, eps, amp = _case(npar, ndet, navg, seed=13)
    amp[:, 5] = 0.0                         # forbidden partner: exact zero coupling
    e_i = float(eps[5])
    hess = np.zeros((npar, npar))
    ref = np.zeros((npar, npar))
    assert _lib_hess_relax(npar, ndet, 0, ovl, weights, eps, e_i, amp, hess, 1e-8)
    _python_reference(npar, ndet, 0, ovl, weights, eps, e_i, amp, ref, 1e-8)
    np.testing.assert_allclose(hess, ref, rtol=0, atol=1e-11 * max(1.0, np.abs(ref).max()))
