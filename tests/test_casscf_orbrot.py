"""Equivalence pins for the CASSCF-side Fortran engines added with the native
wavefunction stack:

  * ``casscf_orbrot.F90``  -- C <- C exp(K) (kappa build + matrix exponential)
  * ``casscf_ah.F90``      -- AH model step, lowest-mode step, DIIS coefficients
  * ``casscf_kernel.F90``  -- the canonicalization mean-field Fock
  * ``rdm_kernel.F90``     -- the spin-summed spatial 1-/2-RDMs

Each engine replaces a NumPy/scipy expression that stays in the tree precisely
so it can serve as this pin.  The engines are algebraically identical to their
references, so these tests hold them to round-off rather than to a tolerance
that would hide a real discrepancy -- and the two spatial RDMs are held to
BIT equality, because they are accumulated in the reference's order on purpose
(a reassociated spin sum was observed to move a near-degenerate H4 MCQDPT2
excited state in its 8th decimal).

The augmented-Hessian step is checked over randomized cases that include the
regime the optimizer actually converges into -- small gradient, one deep
negative curvature -- because that is where the rejected arrowhead shortcut
failed, and where a future attempt would fail again.
"""
import numpy as np
import pytest
from scipy.linalg import expm

from oqp.library import casscf_convergers as CC
from oqp.library.casscf import (
    CASSCF,
    _effective_fock_backend,
    _kappa_matrix,
    _lib_effective_fock,
    _nonredundant_pairs,
    _orbital_rotate,
    _rotate_backend,
)
from oqp.library.rdm import (
    _finalize,
    _lib_rdm1_spatial,
    _lib_rdm2_spatial,
    determinant_basis,
    make_rdm1_spatial,
    make_rdm1_spinorb,
    make_rdm2_spatial,
    make_rdm2_spinorb,
)


def _symmetric_eri(nbf, rng):
    """Random ERIs carrying the full 8-fold permutational symmetry."""
    a = rng.standard_normal((nbf, nbf, nbf, nbf))
    a = a + a.transpose(1, 0, 2, 3)
    a = a + a.transpose(0, 1, 3, 2)
    return a + a.transpose(2, 3, 0, 1)


# --------------------------------------------------------------- orbital rotation
rotate_pin = pytest.mark.skipif(
    _rotate_backend() is None, reason="liboqp casscf_orbital_rotate unavailable")


@rotate_pin
@pytest.mark.parametrize("nbf,ncore,nact", [(10, 2, 4), (24, 2, 6), (31, 3, 6)])
@pytest.mark.parametrize("scale", [1.0e-8, 1.0e-4, 0.2, 1.0])
def test_orbital_rotation_matches_expm(nbf, ncore, nact, scale):
    rng = np.random.default_rng(hash((nbf, ncore, nact)) % 2**32)
    pairs = _nonredundant_pairs(ncore, nact, nbf)
    vec = scale * rng.standard_normal(len(pairs))
    C = rng.standard_normal((nbf, nbf))
    ref = C @ expm(_kappa_matrix(vec, pairs, nbf))
    got = _orbital_rotate(C, vec, pairs, nbf)
    assert np.allclose(got, ref, rtol=0.0, atol=1.0e-11 * max(1.0, np.abs(ref).max()))


@rotate_pin
def test_orbital_rotation_accumulates_repeated_pairs():
    """``_kappa_matrix`` uses ``+=``/``-=``; the engine must do the same."""
    rng = np.random.default_rng(5)
    pairs = [(3, 1), (3, 1), (2, 0)]
    vec = np.array([0.3, -0.11, 0.7])
    C = rng.standard_normal((5, 5))
    ref = C @ expm(_kappa_matrix(vec, pairs, 5))
    assert np.allclose(_orbital_rotate(C, vec, pairs, 5), ref, rtol=0.0, atol=1.0e-13)


@rotate_pin
def test_rotation_of_orthonormal_orbitals_stays_orthonormal():
    """exp(antisymmetric) is orthogonal, so the rotation preserves the metric."""
    rng = np.random.default_rng(9)
    nbf, ncore, nact = 20, 2, 5
    pairs = _nonredundant_pairs(ncore, nact, nbf)
    C = np.linalg.qr(rng.standard_normal((nbf, nbf)))[0]
    Cn = _orbital_rotate(C, 0.4 * rng.standard_normal(len(pairs)), pairs, nbf)
    assert np.allclose(Cn.T @ Cn, np.eye(nbf), rtol=0.0, atol=1.0e-12)


# ------------------------------------------------------------------- AH model step
ah_pin = pytest.mark.skipif(
    CC._ah_backend("casscf_ah_model_step") is None,
    reason="liboqp casscf_ah_model_step unavailable")


def _hessian_case(npar, rng, deep_negative):
    H = rng.standard_normal((npar, npar))
    H = 0.5 * (H + H.T)
    if deep_negative:
        w, U = np.linalg.eigh(H)
        w[0] = -abs(w[0]) - 0.5
        H = U @ np.diag(w) @ U.T
    w, U = np.linalg.eigh(H)
    return w, U


@ah_pin
@pytest.mark.parametrize("npar", [6, 20, 60])
@pytest.mark.parametrize("deep_negative", [False, True])
@pytest.mark.parametrize("trust", [1.0e-3, 0.05, 0.2, 1.0])
def test_ah_model_step_matches_reference(npar, deep_negative, trust):
    rng = np.random.default_rng(hash((npar, deep_negative)) % 2**32)
    for _ in range(4):
        w, U = _hessian_case(npar, rng, deep_negative)
        # sweep the gradient scale: the interesting regime is a small gradient
        # against one deep negative mode.
        grad = 10.0 ** rng.uniform(-9, 0) * rng.standard_normal(npar)
        ref = CC._ah_model_step_reference(grad, w, U, trust, 32)
        got = CC._ah_model_step(grad, w, U, trust, 32)

        # the no-reference-component verdict must agree exactly: disagreeing
        # here is what the rejected arrowhead solver did wrong
        assert (ref[0] is None) == (got[0] is None)
        assert ref[3] == got[3], "microiteration count diverged"
        assert got[1] == pytest.approx(ref[1], rel=1.0e-11, abs=1.0e-13)
        if ref[0] is None:
            continue
        scale = max(1.0, np.abs(ref[0]).max())
        assert np.allclose(got[0], ref[0], rtol=0.0, atol=1.0e-7 * scale)
        assert got[2] == pytest.approx(ref[2], rel=1.0e-6, abs=1.0e-14)


@ah_pin
@pytest.mark.parametrize("npar", [4, 25, 60])
def test_lowest_mode_step_matches_reference(npar):
    rng = np.random.default_rng(npar)
    for _ in range(10):
        w, U = _hessian_case(npar, rng, deep_negative=True)
        grad = rng.standard_normal(npar)
        for trust in (1.0e-4, 0.3):
            ref = CC._lowest_mode_step_reference(grad, w, U, trust)
            got = CC._lowest_mode_step(grad, w, U, trust)
            assert np.allclose(got[0], ref[0], rtol=0.0, atol=1.0e-14)
            # pred is a difference of two cancelling terms, so it is pinned
            # relatively loosely on purpose
            assert got[1] == pytest.approx(ref[1], rel=1.0e-9, abs=1.0e-14)


# ---------------------------------------------------------------------- DIIS
diis_pin = pytest.mark.skipif(
    CC._ah_backend("casscf_diis_coeffs") is None,
    reason="liboqp casscf_diis_coeffs unavailable")


@diis_pin
@pytest.mark.parametrize("npar", [5, 40, 140])
@pytest.mark.parametrize("nvec", [2, 3, 5, 8])
def test_diis_coefficients_match_reference(npar, nvec):
    rng = np.random.default_rng(hash((npar, nvec)) % 2**32)
    for trial in range(8):
        gs = [rng.standard_normal(npar) * 10.0 ** rng.uniform(-8, 0)
              for _ in range(nvec)]
        if trial % 4 == 2 and nvec >= 3:
            gs[-1] = gs[-2] * (1.0 + 1.0e-13)        # near-dependent -> drop path
        if trial % 4 == 3:
            gs = [g * 10.0 ** (-2.0 * k) for k, g in enumerate(gs)]   # converging
        ref = CC._diis_coefficients_reference(gs)
        got = CC._diis_coefficients(gs)
        assert (ref is None) == (got is None)
        if ref is None:
            continue
        assert len(got) == len(ref), "different number of vectors retained"
        assert np.allclose(got, ref, rtol=1.0e-8, atol=1.0e-10)


# ------------------------------------------------------------- effective Fock
@pytest.mark.skipif(_effective_fock_backend() is None,
                    reason="liboqp casscf_effective_fock unavailable")
@pytest.mark.parametrize("nbf", [4, 9, 14])
def test_effective_fock_matches_einsum_reference(nbf):
    rng = np.random.default_rng(nbf)
    h1e = rng.standard_normal((nbf, nbf))
    h1e = h1e + h1e.T
    eri = _symmetric_eri(nbf, rng)
    D = rng.standard_normal((nbf, nbf))
    D = D + D.T
    ref = CASSCF._effective_fock_reference(h1e, eri, D)
    got = _lib_effective_fock(h1e, eri, D)
    assert np.allclose(got, ref, rtol=0.0, atol=1.0e-11 * max(1.0, np.abs(ref).max()))
    # and through the dispatching entry point the driver actually calls
    assert np.allclose(CASSCF._effective_fock(h1e, eri, D, 0, 0), ref,
                       rtol=0.0, atol=1.0e-11 * max(1.0, np.abs(ref).max()))


# ------------------------------------------------------------- spatial RDMs
@pytest.mark.skipif(_lib_rdm1_spatial(np.ones(1), [3], 1) is None,
                    reason="liboqp rdm1_spatial unavailable")
@pytest.mark.parametrize("norb,nelec", [(2, (1, 1)), (4, (2, 2)), (5, (3, 2)),
                                        (6, (3, 3))])
def test_spatial_rdms_are_bit_identical_to_spinorb_collapse(norb, nelec):
    """The engines must reproduce the spin-orbital build plus the NumPy spin
    sum EXACTLY, not just to round-off.

    Bit equality is the contract: these RDMs feed CASSCF, CASPT2/NEVPT2 and
    MCQDPT2, and a reassociated sum here is a last-digit change in a printed
    excited-state energy.
    """
    rng = np.random.default_rng(hash((norb, nelec)) % 2**32)
    dets = determinant_basis(norb, nelec)
    for _ in range(3):
        c = rng.standard_normal(len(dets))
        c /= np.linalg.norm(c)

        spin1 = make_rdm1_spinorb(c, dets, 2 * norb)
        ref1 = _finalize(spin1[:norb, :norb] + spin1[norb:, norb:])
        assert np.array_equal(make_rdm1_spatial(c, dets, norb), ref1)

        spin2 = make_rdm2_spinorb(c, dets, 2 * norb)
        ref2 = np.zeros((norb,) * 4)
        for so in (0, norb):
            for to in (0, norb):
                blk = spin2[so:so + norb, to:to + norb, so:so + norb, to:to + norb]
                ref2 += np.transpose(blk, (0, 2, 1, 3))
        assert np.array_equal(make_rdm2_spatial(c, dets, norb), _finalize(ref2))


@pytest.mark.skipif(_lib_rdm2_spatial(np.ones(1), [3], 1) is None,
                    reason="liboqp rdm2_spatial unavailable")
def test_spatial_rdm2_contracts_to_the_active_energy():
    """0.5 * sum (pq|rs) D2[p,q,r,s] + sum h[p,q] D1[p,q] is the CI energy.

    An independent check on the engine's index convention: a transposed or
    mis-summed spin block would still be symmetric-looking but would not
    reproduce <psi|H|psi>.
    """
    rng = np.random.default_rng(3)
    norb, nelec = 4, (2, 2)
    dets = determinant_basis(norb, nelec)
    h1e = rng.standard_normal((norb, norb))
    h1e = h1e + h1e.T
    eri = _symmetric_eri(norb, rng)
    c = rng.standard_normal(len(dets))
    c /= np.linalg.norm(c)

    d1 = make_rdm1_spatial(c, dets, norb)
    d2 = make_rdm2_spatial(c, dets, norb)
    energy = float(np.einsum("pq,pq->", h1e, d1) + 0.5 * np.einsum("pqrs,pqrs->", eri, d2))

    # explicit <psi|H|psi> from the spin-orbital RDMs
    s1 = make_rdm1_spinorb(c, dets, 2 * norb)
    s2 = make_rdm2_spinorb(c, dets, 2 * norb)
    hs = np.zeros((2 * norb, 2 * norb))
    for off in (0, norb):
        hs[off:off + norb, off:off + norb] = h1e
    gs = np.zeros((2 * norb,) * 4)
    for a in (0, norb):
        for b in (0, norb):
            gs[a:a + norb, a:a + norb, b:b + norb, b:b + norb] = eri
    ref = float(np.einsum("pq,pq->", hs, s1)
                + 0.5 * np.einsum("prqs,pqrs->", gs, s2))
    assert energy == pytest.approx(ref, rel=1.0e-12, abs=1.0e-12)
