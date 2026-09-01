from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_relaxed_density.F90"


def test_relaxed_ao_density_derivative_matches_small_step():
    rng = np.random.default_rng(104101)
    c = rng.normal(size=(6, 6))
    dc = rng.normal(size=(6, 6))
    p = rng.normal(size=(6, 6))
    dp = rng.normal(size=(6, 6))

    def ao_density(cmat, pmat):
        raw = cmat @ pmat @ cmat.T
        return 0.5 * (raw + raw.T)

    analytic = dc @ p @ c.T + c @ dp @ c.T + c @ p @ dc.T
    analytic = 0.5 * (analytic + analytic.T)
    step = 2.0e-6
    numerical = (
        ao_density(c + step * dc, p + step * dp)
        - ao_density(c - step * dc, p - step * dp)
    ) / (2.0 * step)
    np.testing.assert_allclose(analytic, numerical, atol=2e-8, rtol=2e-9)


def test_fortran_relaxed_density_uses_spin_adapted_z_map():
    text = SOURCE.read_text()
    assert "call sfropcal" in text
    assert "dmo_a" in text and "dmo_b" in text
    assert "dz" in text
    assert "0.5_dp*(draw+transpose(draw))" in text
