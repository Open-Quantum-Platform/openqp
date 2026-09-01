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
