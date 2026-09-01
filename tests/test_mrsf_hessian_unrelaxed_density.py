from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_unrelaxed_density.F90"


def test_separated_polarization_is_exact_for_density_polynomials():
    rng = np.random.default_rng(20260901)
    c = rng.normal(size=(5, 5))
    dc = rng.normal(size=(5, 5))
    x = rng.normal(size=(5, 5))
    dx = rng.normal(size=(5, 5))

    def linear_state_map(cmat, xmat):
        return cmat @ xmat @ cmat.T

    analytic_linear = (
        dc @ x @ c.T + c @ x @ dc.T + c @ dx @ c.T
    )
    polarized_linear = 0.5 * (
        linear_state_map(c + dc, x) - linear_state_map(c - dc, x)
    ) + linear_state_map(c, dx)
    np.testing.assert_allclose(polarized_linear, analytic_linear, atol=2e-13)

    def quadratic_state_map(cmat, xmat):
        y = cmat @ xmat
        return y @ y.T

    analytic_quadratic = (
        dc @ x @ x.T @ c.T
        + c @ x @ x.T @ dc.T
        + c @ dx @ x.T @ c.T
        + c @ x @ dx.T @ c.T
    )
    polarized_quadratic = 0.5 * (
        quadratic_state_map(c + dc, x) - quadratic_state_map(c - dc, x)
        + quadratic_state_map(c, x + dx) - quadratic_state_map(c, x - dx)
    )
    np.testing.assert_allclose(polarized_quadratic, analytic_quadratic, atol=2e-13)


def test_fortran_path_is_spin_adapted_and_contains_no_geometry_fd():
    text = SOURCE.read_text()
    assert "call mrsfxvec" in text
    assert "call mrsfcbc" in text
    assert "call sfdmat" in text
    assert "build_mrsf_density_orbital_derivative" in text
    forbidden = ("type(determinant", "slater_", "displaced_geometry", "finite_difference")
    lowered = text.lower()
    for token in forbidden:
        assert token not in lowered
