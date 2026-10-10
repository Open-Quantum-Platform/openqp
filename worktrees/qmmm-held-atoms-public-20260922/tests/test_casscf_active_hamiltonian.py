"""Pin for routing the analytic Hessian's active operator through liboqp.

``casscf_hessian._active_hamiltonian`` builds the determinant-basis matrix of
``sum f_tu E_tu + 1/2 sum g_tuvw (E_tu E_vw - d_uv E_tw)`` with the same Fortran
engine the CI solver uses, instead of the ``nact^2 ndet^3`` tensordot in
``_active_operator_matrix``.  The claim being pinned is that these are the same
matrix -- the active operator of the folded integrals IS the CASCI Hamiltonian
of those integrals -- so the NumPy form stays as the reference.
"""
import numpy as np
import pytest

from oqp.library.casscf_hessian import (
    _active_hamiltonian,
    _active_operator_matrix,
    _excitation_matrices,
)


def _folded(nact, seed):
    """Random folded active integrals with the symmetries the real ones carry."""
    rng = np.random.default_rng(seed)
    f = rng.standard_normal((nact, nact))
    f = f + f.T
    g = rng.standard_normal((nact,) * 4)
    g = g + g.transpose(1, 0, 2, 3)
    g = g + g.transpose(0, 1, 3, 2)
    g = g + g.transpose(2, 3, 0, 1)
    return f, g


@pytest.mark.parametrize(
    "nact,nalpha,nbeta",
    [(4, 2, 2), (5, 3, 2), (6, 3, 3), (4, 2, 1), (3, 1, 1), (5, 2, 2)],
)
def test_active_hamiltonian_matches_the_tensordot_reference(nact, nalpha, nbeta):
    dets, stack = _excitation_matrices(nact, nalpha, nbeta)
    f, g = _folded(nact, seed=nact * 41 + nalpha)

    ref = _active_operator_matrix(f, g, stack)
    got = _active_hamiltonian(f, g, stack, dets, nact)

    assert got.shape == ref.shape
    np.testing.assert_allclose(got, ref, rtol=0, atol=1e-10 * np.abs(ref).max())


def test_active_hamiltonian_is_symmetric_and_preserves_the_spectrum():
    """What the caller actually consumes: the eigenpairs and <c|H|c>."""
    nact, nalpha, nbeta = 6, 3, 3
    dets, stack = _excitation_matrices(nact, nalpha, nbeta)
    f, g = _folded(nact, seed=7)

    ref = _active_operator_matrix(f, g, stack)
    got = _active_hamiltonian(f, g, stack, dets, nact)

    assert np.abs(got - got.T).max() <= 1e-10 * np.abs(got).max()

    w_ref = np.linalg.eigvalsh(0.5 * (ref + ref.T))
    w_got = np.linalg.eigvalsh(0.5 * (got + got.T))
    np.testing.assert_allclose(w_got, w_ref, rtol=0, atol=1e-10 * max(1.0, abs(w_ref).max()))

    rng = np.random.default_rng(11)
    c = rng.standard_normal(len(dets))
    c /= np.linalg.norm(c)
    assert abs(float(c @ (got @ c)) - float(c @ (ref @ c))) <= 1e-10 * max(1.0, abs(float(c @ (ref @ c))))
