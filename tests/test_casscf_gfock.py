"""Equivalence pin for the Fortran CASSCF Fock/gradient engine (casscf_kernel.F90).

The engine replaces the Python assembly in ``casscf.py``

    D, G = _full_rdms(gamma, Gamma, ...)          # explicit nbf^4 2-RDM
    F    = _generalized_fock(D, G, h1e, eri)      # O(nbf^5) einsum
    g    = 2 (F[q,p] - F[p,q])  over _nonredundant_pairs

with an O(nbf^4) J/K build that never forms G.  The two are algebraically
identical, so these tests hold them to round-off rather than to a tolerance
that would hide a real discrepancy.

The Python path stays in the tree precisely so it can serve as this pin.
"""
import numpy as np
import pytest

from oqp.library.casscf import (
    _averaged_rdm2,
    _full_rdms,
    _generalized_fock,
    _gfock_backend,
    _lib_gfock_grad,
    _nonredundant_pairs,
)

pytestmark = pytest.mark.skipif(
    _gfock_backend() is None,
    reason="liboqp casscf_gfock_grad unavailable",
)


def _symmetric_eri(nbf, rng):
    """Random ERIs carrying the full 8-fold permutational symmetry."""
    a = rng.standard_normal((nbf, nbf, nbf, nbf))
    a = a + a.transpose(1, 0, 2, 3)
    a = a + a.transpose(0, 1, 3, 2)
    return a + a.transpose(2, 3, 0, 1)


def _random_case(nbf, ncore, nact, nroot, seed):
    rng = np.random.default_rng(seed)
    h1e = rng.standard_normal((nbf, nbf))
    h1e = h1e + h1e.T
    eri = _symmetric_eri(nbf, rng)
    weights = rng.random(nroot)
    weights /= weights.sum()
    gammas = np.zeros((nroot, nact, nact))
    Gammas = np.zeros((nroot,) + (nact,) * 4)
    for r in range(nroot):
        g = rng.standard_normal((nact, nact))
        gammas[r] = g + g.T
        G = rng.standard_normal((nact,) * 4)
        G = G + G.transpose(1, 0, 2, 3)
        G = G + G.transpose(0, 1, 3, 2)
        Gammas[r] = G + G.transpose(2, 3, 0, 1)
    return h1e, eri, weights, gammas, Gammas


def _reference(h1e, eri, weights, gammas, Gammas, ncore, nact, nbf):
    """The Python Fock/gradient path, verbatim."""
    D = np.zeros((nbf, nbf))
    for w, gamma, Gamma in zip(weights, gammas, Gammas):
        D += w * _full_rdms(gamma, Gamma, ncore, nact, nbf)[0]
    G = _averaged_rdm2(gammas, Gammas, weights, ncore, nact, nbf)
    F = _generalized_fock(D, G, h1e, eri)
    pairs = _nonredundant_pairs(ncore, nact, nbf)
    grad = np.array([2.0 * (F[q, p] - F[p, q]) for (p, q) in pairs])
    return F, grad


# (nbf, ncore, nact, nroot) -- the last three are the degenerate shapes that
# an index slip in the kernel would survive on the generic case:
#   no inactive block, no virtual block, active space spanning the remainder.
@pytest.mark.parametrize(
    "nbf,ncore,nact,nroot",
    [
        (10, 2, 4, 1),
        (14, 3, 5, 1),
        (10, 2, 4, 3),
        (14, 4, 6, 2),
        (12, 0, 5, 1),
        (12, 4, 8, 1),
        (9, 3, 6, 1),
    ],
)
def test_fortran_gfock_matches_python(nbf, ncore, nact, nroot):
    h1e, eri, weights, gammas, Gammas = _random_case(nbf, ncore, nact, nroot, seed=nbf * 31 + nroot)
    F_ref, g_ref = _reference(h1e, eri, weights, gammas, Gammas, ncore, nact, nbf)

    npair = len(_nonredundant_pairs(ncore, nact, nbf))
    built = _lib_gfock_grad(nbf, ncore, nact, weights, gammas, Gammas, h1e, eri, npair)
    assert built is not None
    F_lib, g_lib = built

    assert F_lib.shape == (nbf, nbf)
    assert g_lib.shape == (npair,)
    np.testing.assert_allclose(F_lib, F_ref, rtol=0, atol=1e-11 * max(1.0, np.abs(F_ref).max()))
    np.testing.assert_allclose(g_lib, g_ref, rtol=0, atol=1e-11 * max(1.0, np.abs(g_ref).max()))


def _averaged_build(nbf, ncore, nact, weights, gammas, Gammas, h1e, eri, npair):
    """One build from the weight-averaged RDMs, instead of one build per root."""
    avg_gamma = np.einsum("r,rij->ij", weights, gammas)[None]
    avg_Gamma = np.einsum("r,rijkl->ijkl", weights, Gammas)[None]
    return _lib_gfock_grad(nbf, ncore, nact, np.array([1.0]), avg_gamma,
                           avg_Gamma, h1e, eri, npair)[0]


def test_state_average_loop_matches_averaged_rdms_when_normalized():
    """Averaging the RDMs first is equivalent -- but only for normalized weights.

    ``G_sep`` is quadratic in D, which suggests the per-root loop cannot be
    collapsed.  It can: the term quadratic in the root-dependent (active) part
    of D has all four indices active, and that is exactly the block the dG
    correction overwrites, so it never reaches F.  What does *not* survive
    collapsing is an unnormalized weight set -- the core-core term picks up
    ``sum_I w_I`` from the per-root loop and ``(sum_I w_I)^2`` from the
    averaged build.

    CASSCF always normalizes, so both agree in production; the kernel keeps the
    per-root loop so it does not silently depend on that, and this test pins
    both halves of the statement.
    """
    nbf, ncore, nact, nroot = 12, 2, 4, 2
    h1e, eri, weights, gammas, Gammas = _random_case(nbf, ncore, nact, nroot, seed=7)
    npair = len(_nonredundant_pairs(ncore, nact, nbf))

    weights = weights / weights.sum()
    F_lib, _g = _lib_gfock_grad(nbf, ncore, nact, weights, gammas, Gammas, h1e, eri, npair)
    F_ref, _g_ref = _reference(h1e, eri, weights, gammas, Gammas, ncore, nact, nbf)
    np.testing.assert_allclose(F_lib, F_ref, rtol=0, atol=1e-11 * np.abs(F_ref).max())

    F_avg = _averaged_build(nbf, ncore, nact, weights, gammas, Gammas, h1e, eri, npair)
    np.testing.assert_allclose(F_avg, F_lib, rtol=0, atol=1e-11 * np.abs(F_lib).max())

    # ... and the per-root loop is what keeps unnormalized weights correct.
    off = weights * 1.7
    F_off, _ = _lib_gfock_grad(nbf, ncore, nact, off, gammas, Gammas, h1e, eri, npair)
    F_off_ref, _ = _reference(h1e, eri, off, gammas, Gammas, ncore, nact, nbf)
    np.testing.assert_allclose(F_off, F_off_ref, rtol=0, atol=1e-11 * np.abs(F_off_ref).max())
    assert np.abs(_averaged_build(nbf, ncore, nact, off, gammas, Gammas, h1e, eri, npair)
                  - F_off).max() > 1e-6 * np.abs(F_off).max()
