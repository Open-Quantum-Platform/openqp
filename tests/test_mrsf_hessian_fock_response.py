"""Independent algebra for the ROHF spin-Fock nuclear derivative."""

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FOCK_RESPONSE = ROOT / "source/modules/tdhf_mrsf_hessian_fock_response.F90"
FOCK_DERIV = ROOT / "source/modules/fock_deriv.F90"


def _symmetrized_eri(rng, nbf):
    eri = rng.normal(size=(nbf, nbf, nbf, nbf))
    for permutation in (
        (1, 0, 2, 3),
        (0, 1, 3, 2),
        (2, 3, 0, 1),
    ):
        eri = eri + eri.transpose(permutation)
    return eri / 8.0


def _spin_fock(hcore, eri, pa, pb, exchange_scale, spin_density):
    total = pa + pb
    coulomb = np.einsum("uvls,ls->uv", eri, total, optimize=True)
    exchange = np.einsum("ulvs,ls->uv", eri, spin_density, optimize=True)
    return hcore + coulomb - exchange_scale * exchange


def test_spin_fock_total_derivative_matches_central_difference():
    rng = np.random.default_rng(149104101)
    nbf = 4
    h = rng.normal(size=(nbf, nbf)); h = 0.5 * (h + h.T)
    dh = rng.normal(size=(nbf, nbf)); dh = 0.5 * (dh + dh.T)
    eri = _symmetrized_eri(rng, nbf)
    deri = _symmetrized_eri(rng, nbf)
    pa = rng.normal(size=(nbf, nbf)); pa = 0.5 * (pa + pa.T)
    pb = rng.normal(size=(nbf, nbf)); pb = 0.5 * (pb + pb.T)
    dpa = rng.normal(size=(nbf, nbf)); dpa = 0.5 * (dpa + dpa.T)
    dpb = rng.normal(size=(nbf, nbf)); dpb = 0.5 * (dpb + dpb.T)
    scale = 0.54

    explicit_coulomb = np.einsum(
        "uvls,ls->uv", deri, pa + pb, optimize=True
    )
    explicit_exchange_a = np.einsum("ulvs,ls->uv", deri, pa, optimize=True)
    response_coulomb = np.einsum(
        "uvls,ls->uv", eri, dpa + dpb, optimize=True
    )
    response_exchange_a = np.einsum("ulvs,ls->uv", eri, dpa, optimize=True)
    analytic = (
        dh + explicit_coulomb - scale * explicit_exchange_a
        + response_coulomb - scale * response_exchange_a
    )

    step = 1.0e-5
    plus = _spin_fock(
        h + step * dh, eri + step * deri,
        pa + step * dpa, pb + step * dpb, scale, pa + step * dpa,
    )
    minus = _spin_fock(
        h - step * dh, eri - step * deri,
        pa - step * dpa, pb - step * dpb, scale, pa - step * dpa,
    )
    np.testing.assert_allclose((plus - minus) / (2.0 * step), analytic, atol=2e-10)


def test_orbital_density_derivative_matches_finite_difference():
    rng = np.random.default_rng(20260902)
    mo, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    dmo = rng.normal(size=(6, 6))
    nocc = 3
    analytic = (
        dmo[:, :nocc] @ mo[:, :nocc].T
        + mo[:, :nocc] @ dmo[:, :nocc].T
    )
    step = 1.0e-5
    plus = (mo[:, :nocc] + step * dmo[:, :nocc])
    minus = (mo[:, :nocc] - step * dmo[:, :nocc])
    finite_difference = (
        plus @ plus.T - minus @ minus.T
    ) / (2.0 * step)
    np.testing.assert_allclose(finite_difference, analytic, atol=2e-11)


def test_production_builder_has_explicit_and_relaxed_density_terms():
    response = "".join(FOCK_RESPONSE.read_text().lower().split())
    derivative = "".join(FOCK_DERIV.read_text().lower().split())
    assert "callfock_deriv_matrix_os" in response
    assert "callfock_jk" in response
    assert "subroutinebuild_orbital_density_derivatives" in response
    assert "hcore_derivative" in response
    assert "if(dft.and.infos%dft%cam_flag)thenstatus=-3" in response
    assert ".not.present(dvxc_a)" in response
    assert "subroutinefock_deriv_matrix_os" in derivative
    assert "fmat(mu,nu,:,:)=gx" in derivative
    assert "fmat(mu,nu,:,:)=2.0_dp*gx" not in derivative.split(
        "subroutinefock_deriv_matrix_os", 1
    )[1].split("endsubroutinefock_deriv_matrix_os", 1)[0]
