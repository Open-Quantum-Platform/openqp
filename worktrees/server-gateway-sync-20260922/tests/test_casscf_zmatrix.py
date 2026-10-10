"""Pin for the symmetry argument behind ``casscf_hessian._z_matrix``.

``_z_matrix`` evaluates one four-index contraction where the reference form
``_z_matrix_general`` evaluates four, because all four are the same number
whenever both operands carry the pairwise index symmetry documented on
``_z_matrix``.  These tests pin both halves of that claim: that the operands
really do carry the symmetry, and that the two forms then agree to round-off.
"""
import numpy as np
import pytest

from oqp.library.casscf import _full_rdms
from oqp.library.casscf_hessian import (
    _one_index_derivative_g,
    _one_index_derivative_h,
    _z_matrix,
    _z_matrix_general,
)


def _eri(nbf, rng):
    """Random ERIs carrying the full eight-fold permutational symmetry."""
    a = rng.standard_normal((nbf, nbf, nbf, nbf))
    a = a + a.transpose(1, 0, 2, 3)
    a = a + a.transpose(0, 1, 3, 2)
    return a + a.transpose(2, 3, 0, 1)


def _rdms(nbf, ncore, nact, rng):
    gamma = rng.standard_normal((nact, nact))
    gamma = gamma + gamma.T
    Gamma = rng.standard_normal((nact,) * 4)
    Gamma = Gamma + Gamma.transpose(1, 0, 3, 2)   # G[pqrs] = G[qpsr]
    Gamma = Gamma + Gamma.transpose(2, 3, 0, 1)   # G[pqrs] = G[rspq]
    return _full_rdms(gamma, Gamma, ncore, nact, nbf)


def test_one_index_derivative_inherits_the_eri_symmetry():
    """The step the docstring of the old code got wrong: T *is* symmetric.

    The one-index derivative applies the same generator to all four slots, so
    it preserves every permutational symmetry of the ERIs it is built from.
    """
    rng = np.random.default_rng(0)
    eri = _eri(9, rng)
    T = _one_index_derivative_g(eri, 2, 6)
    scale = np.abs(T).max()
    for axes in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
        assert np.abs(T - T.transpose(axes)).max() <= 1e-12 * scale


def test_full_rdm2_carries_the_pair_symmetry():
    """G[pqrs] = G[qpsr] = G[rspq] -- but NOT G[pqrs] = G[qprs]."""
    rng = np.random.default_rng(1)
    _D, G = _rdms(9, 2, 4, rng)
    scale = np.abs(G).max()
    assert np.abs(G - G.transpose(1, 0, 3, 2)).max() <= 1e-12 * scale
    assert np.abs(G - G.transpose(2, 3, 0, 1)).max() <= 1e-12 * scale
    # the symmetry the argument deliberately does not rely on
    assert np.abs(G - G.transpose(1, 0, 2, 3)).max() > 1e-6 * scale


@pytest.mark.parametrize("nbf,ncore,nact", [(9, 2, 4), (8, 0, 4), (10, 3, 5)])
def test_z_matrix_matches_the_general_form_on_derivative_integrals(nbf, ncore, nact):
    rng = np.random.default_rng(nbf * 17 + nact)
    h1e = rng.standard_normal((nbf, nbf))
    h1e = h1e + h1e.T
    eri = _eri(nbf, rng)
    D, G = _rdms(nbf, ncore, nact, rng)

    for (p, q) in ((0, nbf - 1), (ncore, ncore + 1), (nact, nbf - 2)):
        t = _one_index_derivative_h(h1e, p, q)
        T = _one_index_derivative_g(eri, p, q)
        fast = _z_matrix(D, G, t, T)
        ref = _z_matrix_general(D, G, t, T)
        np.testing.assert_allclose(
            fast, ref, rtol=0, atol=1e-10 * max(1.0, np.abs(ref).max()))


def test_z_matrix_matches_on_the_undifferentiated_integrals():
    """The gradient call site: Z(h, g) must also agree with the general form."""
    rng = np.random.default_rng(5)
    nbf, ncore, nact = 9, 2, 4
    h1e = rng.standard_normal((nbf, nbf))
    h1e = h1e + h1e.T
    eri = _eri(nbf, rng)
    D, G = _rdms(nbf, ncore, nact, rng)
    np.testing.assert_allclose(
        _z_matrix(D, G, h1e, eri), _z_matrix_general(D, G, h1e, eri),
        rtol=0, atol=1e-10 * np.abs(_z_matrix_general(D, G, h1e, eri)).max())
