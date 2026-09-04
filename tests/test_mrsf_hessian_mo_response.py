from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_mo_response.F90"


def test_mo_fock_and_orbital_energy_derivatives():
    rng = np.random.default_rng(150184111)
    c, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    dc = rng.normal(size=(6, 6))
    f = rng.normal(size=(6, 6))
    f = 0.5 * (f + f.T)
    df = rng.normal(size=(6, 6))
    df = 0.5 * (df + df.T)
    derivative = dc.T @ f @ c + c.T @ df @ c + c.T @ f @ dc
    step = 1.0e-5
    numerical = (
        (c + step * dc).T @ (f + step * df) @ (c + step * dc)
        - (c - step * dc).T @ (f - step * df) @ (c - step * dc)
    ) / (2.0 * step)
    np.testing.assert_allclose(derivative, numerical, atol=2e-8)
    np.testing.assert_allclose(np.diag(derivative), np.diag(numerical), atol=2e-8)


def test_fortran_uses_shared_moving_basis_derivative():
    text = SOURCE.read_text()
    assert "transform_mo_operator_derivative" in text
    assert "orbital_energy_derivative" in text
    assert "Hiroya" in text and "Nakata" in text


def test_common_rohf_canonical_rotation_preserves_metric_and_diagonalizes_sectors():
    rng = np.random.default_rng(149104101)
    nbf = 7
    c, _ = np.linalg.qr(rng.normal(size=(nbf, nbf)))
    overlap_derivative = rng.normal(size=(nbf, nbf))
    overlap_derivative = 0.5 * (overlap_derivative + overlap_derivative.T)
    connection = -0.5 * overlap_derivative
    energies = np.array([-3.0, -1.8, -0.7, -0.2, 0.6, 1.5, 2.8])
    dfeff = rng.normal(size=(nbf, nbf))
    dfeff = 0.5 * (dfeff + dfeff.T)
    rotation = np.zeros((nbf, nbf))
    for first, last in ((0, 2), (2, 4), (4, 7)):
        for q in range(first + 1, last):
            for p in range(first, q):
                rotation[p, q] = -dfeff[p, q] / (energies[p] - energies[q])
                rotation[q, p] = -rotation[p, q]
    completed = connection + rotation
    np.testing.assert_allclose(completed + completed.T, -overlap_derivative,
                               atol=2e-15)
    transformed = dfeff + np.diag(energies) @ rotation \
        - rotation @ np.diag(energies)
    for first, last in ((0, 2), (2, 4), (4, 7)):
        block = transformed[first:last, first:last].copy()
        np.fill_diagonal(block, 0.0)
        np.testing.assert_allclose(block, 0.0, atol=2e-15)
    # The AO derivative changes only by a common orbital rotation.
    assert np.linalg.norm(c @ rotation) > 0.0


def test_common_rohf_canonicalizer_fails_closed_at_unresolved_degeneracy():
    text = SOURCE.read_text()
    assert "abs(gap)<=tolerance" in text
    assert "abs(off_diagonal)>10.0_dp*tolerance" in text
    assert "status=-3" in text
