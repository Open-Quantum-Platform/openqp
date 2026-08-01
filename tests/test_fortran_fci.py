"""Equivalence tests for the Fortran determinant-CI engine (fci_hamiltonian.F90).

The Python enumeration/eigensolvers remain the numerical pin; every Fortran
entry point must reproduce them:

* fci_dense_hamiltonian == the Python _iter_hamiltonian_elements build (1e-12);
* fci_hamiltonian_matvec == sparse-Hamiltonian .dot (1e-12), including blocks;
* fci_hamiltonian_diag == the dense diagonal;
* oqp_dsyevd eigenvalues == numpy eigh, eigenvectors reassemble the matrix;
* end-to-end CASCI energies are unchanged with the engine active (the whole
  existing WF suite is the broader pin).

All tests skip when liboqp lacks the symbols (pre-rebuild checkouts).
"""
import numpy as np
import pytest


def _backend_available() -> bool:
    try:
        import oqp
        from oqp.utils import file_utils  # noqa: F401
    except Exception:
        return False
    return bool(getattr(oqp, "BACKEND_AVAILABLE", False)) and hasattr(oqp, "fci_ao_integrals")


def _engine_available() -> bool:
    if not _backend_available():
        return False
    from oqp import lib
    return all(hasattr(lib, s) for s in
               ("oqp_dsyevd", "fci_dense_hamiltonian",
                "fci_hamiltonian_matvec", "fci_hamiltonian_diag"))


def _random_problem(norb, nalpha, nbeta, seed):
    from oqp.library.fci import _determinants, _spin_orbital_integrals

    rng = np.random.default_rng(seed)
    h1e = rng.standard_normal((norb, norb))
    h1e = 0.5 * (h1e + h1e.T)
    eri = rng.standard_normal((norb,) * 4) * 0.1
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)
    dets = _determinants(norb, (nalpha, nbeta))
    det_index = {d: i for i, d in enumerate(dets)}
    hspin, gspin = _spin_orbital_integrals(h1e, eri)
    return dets, det_index, hspin, gspin, 2 * norb


def _python_dense(dets, det_index, hspin, gspin, nspin, cutoff):
    from oqp.library.fci import _iter_hamiltonian_elements

    ndet = len(dets)
    h = np.zeros((ndet, ndet))
    for row, col, value in _iter_hamiltonian_elements(
            dets, det_index, hspin, gspin, nspin, cutoff):
        h[row, col] += value
    return 0.5 * (h + h.T)


@pytest.mark.parametrize("norb,na,nb,seed", [(4, 2, 2, 3), (5, 3, 2, 7), (6, 3, 3, 11)])
def test_dense_hamiltonian_matches_python(norb, na, nb, seed):
    if not _engine_available():
        pytest.skip("liboqp lacks the fci_hamiltonian engine (rebuild the core)")
    from oqp.library.fci import _lib_dense_hamiltonian

    dets, det_index, hspin, gspin, nspin = _random_problem(norb, na, nb, seed)
    h_py = _python_dense(dets, det_index, hspin, gspin, nspin, 0.0)
    h_f = _lib_dense_hamiltonian(dets, hspin, gspin, nspin, 0.0)
    assert h_f is not None
    np.testing.assert_allclose(h_f, h_py, atol=1.0e-12)


def test_matvec_and_diag_match_dense():
    if not _engine_available():
        pytest.skip("liboqp lacks the fci_hamiltonian engine (rebuild the core)")
    from oqp.library.fci import _LibHamiltonianOperator

    dets, det_index, hspin, gspin, nspin = _random_problem(5, 3, 2, 23)
    h_py = _python_dense(dets, det_index, hspin, gspin, nspin, 0.0)
    op = _LibHamiltonianOperator(dets, hspin, gspin, nspin, 0.0)

    np.testing.assert_allclose(op.diagonal(), np.diag(h_py), atol=1.0e-12)

    rng = np.random.default_rng(5)
    x1 = rng.standard_normal(len(dets))
    np.testing.assert_allclose(op.dot(x1), h_py @ x1, atol=1.0e-12)
    xb = rng.standard_normal((len(dets), 3))
    np.testing.assert_allclose(op.dot(xb), h_py @ xb, atol=1.0e-12)

    np.testing.assert_allclose(op.toarray(), h_py, atol=1.0e-12)


def test_oqp_dsyevd_matches_numpy():
    if not _engine_available():
        pytest.skip("liboqp lacks the fci_hamiltonian engine (rebuild the core)")
    from oqp.library.fci import _lib_dsyevd

    rng = np.random.default_rng(17)
    for n in (5, 60, 200):
        a = rng.standard_normal((n, n))
        a = 0.5 * (a + a.T)
        out = _lib_dsyevd(a)
        assert out is not None
        w, v = out
        np.testing.assert_allclose(w, np.linalg.eigvalsh(a), atol=1.0e-9)
        # eigenvectors reassemble the matrix (sign/degeneracy agnostic)
        np.testing.assert_allclose(v @ np.diag(w) @ v.T, a, atol=1.0e-9)
        np.testing.assert_allclose(v.T @ v, np.eye(n), atol=1.0e-9)


def test_davidson_lib_operator_matches_dense_casci():
    """solve_fci davidson (lib matvec operator) == dense solver energies."""
    if not _engine_available():
        pytest.skip("liboqp lacks the fci_hamiltonian engine (rebuild the core)")
    from oqp.library.fci import solve_fci

    rng = np.random.default_rng(29)
    norb = 6
    h1e = rng.standard_normal((norb, norb))
    h1e = 0.5 * (h1e + h1e.T)
    eri = rng.standard_normal((norb,) * 4) * 0.05
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)

    e_dense, _v1 = solve_fci(h1e, eri, (3, 3), nroot=3, solver="dense",
                             max_det=100000, max_memory=2048)
    e_dav, _v2 = solve_fci(h1e, eri, (3, 3), nroot=3, solver="davidson",
                           max_det=100000, max_memory=2048)
    np.testing.assert_allclose(e_dav, e_dense, atol=1.0e-8)


def _mo_engine_available() -> bool:
    if not _backend_available():
        return False
    from oqp import lib
    return all(hasattr(lib, s) for s in
               ("mo_transform_h1e", "mo_transform_eri",
                "fci_spin_orbital_integrals", "fci_fold_core"))


@pytest.mark.parametrize("nbf", [1, 2, 5, 13, 24])
def test_mo_transform_matches_the_einsum_reference(nbf):
    """mo_transform.F90 reorders the same four quarter transformations, so it
    agrees with the einsum to summation-order rounding, not bit for bit."""
    if not _mo_engine_available():
        pytest.skip("liboqp lacks the MO transformation engine (rebuild the core)")
    from oqp.library.fci import _transform_integrals, _transform_integrals_reference

    rng = np.random.default_rng(101 + nbf)
    hcore = rng.standard_normal((nbf, nbf))
    hcore = 0.5 * (hcore + hcore.T)
    eri = rng.standard_normal((nbf,) * 4)
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)
    coeff = rng.standard_normal((nbf, nbf))

    h_f, e_f = _transform_integrals(hcore, eri, coeff)
    h_p, e_p = _transform_integrals_reference(hcore, eri, coeff)

    scale = max(float(np.max(np.abs(e_p))), 1.0)
    np.testing.assert_allclose(h_f, h_p, atol=1.0e-11 * scale)
    np.testing.assert_allclose(e_f, e_p, atol=1.0e-11 * scale)


@pytest.mark.parametrize("norb", [1, 2, 3, 4, 6, 8, 10])
def test_spin_orbital_engine_is_bit_identical(norb):
    """The spin expansion is a permuted copy with no arithmetic in it, so the
    engine, the numpy build and the element loop must agree exactly."""
    if not _mo_engine_available():
        pytest.skip("liboqp lacks the fci_setup engine (rebuild the core)")
    from oqp.library.fci import (
        _spin_orbital_integrals,
        _spin_orbital_integrals_numpy,
        _spin_orbital_integrals_reference,
    )

    rng = np.random.default_rng(7 * norb)
    h1e = rng.standard_normal((norb, norb))
    h1e = 0.5 * (h1e + h1e.T)
    eri = rng.standard_normal((norb,) * 4)

    h_ref, g_ref = _spin_orbital_integrals_reference(h1e, eri)
    for build in (_spin_orbital_integrals, _spin_orbital_integrals_numpy):
        h_got, g_got = build(h1e, eri)
        assert np.array_equal(h_got, h_ref), build.__name__
        assert np.array_equal(g_got, g_ref), build.__name__


@pytest.mark.parametrize(
    ("norb", "active", "core"),
    [
        (6, [1, 2, 3, 4], [0]),
        (10, [2, 3, 4, 5, 6, 7], [0, 1]),
        (5, [0, 1, 2, 3, 4], []),
        (9, [1, 3, 4, 7], [0, 2]),          # non-sequential explicit selection
        (12, [4, 5, 6, 7], [0, 1, 2, 3]),
    ],
)
def test_fold_core_engine_is_bit_identical(norb, active, core):
    """The fold is a sum of ERI elements in a fixed order; the engine
    reproduces that order, so the result is exact, not merely close."""
    if not _mo_engine_available():
        pytest.skip("liboqp lacks the fci_setup engine (rebuild the core)")
    from oqp.library.fci import _fold_core, _fold_core_reference

    rng = np.random.default_rng(31 + norb)
    h1e = rng.standard_normal((norb, norb))
    h1e = 0.5 * (h1e + h1e.T)
    eri = rng.standard_normal((norb,) * 4)
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)

    h_ref, e_ref = _fold_core_reference(h1e, eri, active, core, ecore=-12.3456789)
    h_got, e_got = _fold_core(h1e, eri, active, core, ecore=-12.3456789)

    assert np.array_equal(h_got, h_ref)
    assert e_got == e_ref


@pytest.mark.parametrize(
    ("norb", "na", "nb"),
    [(2, 2, 0), (4, 2, 2), (5, 3, 2), (6, 3, 3), (6, 4, 2)],
)
def test_spin_square_engine_matches_compute_s2(norb, na, nb):
    """fci_spin_square shares one sort of the CI-ordered determinant list
    across the roots and scales by the norm at the end rather than normalizing
    up front, so it matches the per-root Python walk to rounding."""
    if not _mo_engine_available():
        pytest.skip("liboqp lacks the fci_setup engine (rebuild the core)")
    from oqp import lib
    if not hasattr(lib, "fci_spin_square"):
        pytest.skip("liboqp lacks fci_spin_square (rebuild the core)")
    from oqp.library.fci import _determinants, _lib_spin_square, compute_s2

    dets = _determinants(norb, (na, nb))
    rng = np.random.default_rng(3 * norb + na)
    nvec = min(4, len(dets))
    vectors = rng.standard_normal((len(dets), nvec))

    got = _lib_spin_square(vectors, dets, norb, (na, nb))
    assert got is not None
    index = {d: i for i, d in enumerate(dets)}
    ref = np.array([
        compute_s2(vectors[:, r], dets, norb, (na, nb), det_index=index)
        for r in range(nvec)
    ])
    np.testing.assert_allclose(got, ref, atol=1.0e-12)


def test_spin_square_engine_agrees_on_real_ci_roots():
    """On genuine CI eigenvectors the engine must reproduce the integer
    multiplicities, not merely a close <S^2>."""
    if not _engine_available():
        pytest.skip("liboqp lacks the fci_hamiltonian engine (rebuild the core)")
    from oqp.library.fci import _determinants, fci_spin_diagnostics, solve_fci

    rng = np.random.default_rng(77)
    norb = 4
    h1e = rng.standard_normal((norb, norb))
    h1e = 0.5 * (h1e + h1e.T)
    eri = rng.standard_normal((norb,) * 4) * 0.2
    eri = eri + eri.transpose(1, 0, 2, 3)
    eri = eri + eri.transpose(0, 1, 3, 2)
    eri = eri + eri.transpose(2, 3, 0, 1)

    _e, vecs = solve_fci(h1e, eri, (2, 2), nroot=6, solver="dense",
                         max_det=10000, max_memory=512)
    dets = _determinants(norb, (2, 2))
    s2, mult = fci_spin_diagnostics(vecs, dets, norb, (2, 2))

    index = {d: i for i, d in enumerate(dets)}
    from oqp.library.fci import compute_s2
    ref = np.array([compute_s2(vecs[:, r], dets, norb, (2, 2), det_index=index)
                    for r in range(vecs.shape[1])])
    np.testing.assert_allclose(s2, ref, atol=1.0e-11)
    # every root is a spin eigenstate, so <S^2> is s(s+1) for integer 2s
    np.testing.assert_allclose(
        s2, 0.25 * (mult ** 2 - 1.0), atol=1.0e-8)
