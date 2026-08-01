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
