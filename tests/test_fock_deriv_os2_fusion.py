"""Static and algebraic guards for the fused open-shell spin contraction."""

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FOCK_DERIV = ROOT / "source" / "modules" / "fock_deriv.F90"
GRADIENT = ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
INTERCHANGE = ROOT / "source" / "modules" / "mrsf_nac_interchange.F90"


def _body(source, name):
    return source.split(f"subroutine {name}", 1)[1].split(
        f"end subroutine {name}", 1
    )[0]


def _legacy_term(pcoul, pexch, mmat, hfscale, i, j, k, l):
    coulomb = 4.0 * (
        mmat[i, j] * pcoul[k, l] + mmat[k, l] * pcoul[i, j]
    )
    exchange = 2.0 * hfscale * (
        mmat[i, k] * pexch[j, l]
        + mmat[i, l] * pexch[j, k]
        + pexch[i, k] * mmat[j, l]
        + pexch[i, l] * mmat[j, k]
    )
    return coulomb - exchange


def _fused_term(pcoul, pa, pb, ma, mb, hfscale, i, j, k, l):
    msum = ma + mb
    coulomb = 4.0 * (
        msum[i, j] * pcoul[k, l] + msum[k, l] * pcoul[i, j]
    )
    exchange = 2.0 * hfscale * (
        ma[i, k] * pa[j, l]
        + ma[i, l] * pa[j, k]
        + pa[i, k] * ma[j, l]
        + pa[i, l] * ma[j, k]
        + mb[i, k] * pb[j, l]
        + mb[i, l] * pb[j, k]
        + pb[i, k] * mb[j, l]
        + pb[i, l] * mb[j, k]
    )
    return coulomb - exchange


def test_fused_density_is_termwise_sum_of_two_legacy_spin_products():
    rng = np.random.default_rng(90210)
    nbf = 7
    pa = rng.normal(size=(nbf, nbf))
    pb = rng.normal(size=(nbf, nbf))
    ma = rng.normal(size=(nbf, nbf))
    mb = rng.normal(size=(nbf, nbf))
    pcoul = pa + pb
    for hfscale in (0.0, 0.2, 1.0):
        for _ in range(50):
            i, j, k, l = rng.integers(0, nbf, size=4)
            legacy = _legacy_term(pcoul, pa, ma, hfscale, i, j, k, l)
            legacy += _legacy_term(pcoul, pb, mb, hfscale, i, j, k, l)
            fused = _fused_term(pcoul, pa, pb, ma, mb, hfscale, i, j, k, l)
            np.testing.assert_allclose(fused, legacy, rtol=2e-15, atol=2e-14)


def test_nac_esum_and_hf_adjoint_use_one_fused_derivative_pass():
    fock = FOCK_DERIV.read_text()
    gradient = _body(GRADIENT.read_text(), "mrsf_nac_esum(infos, istate, jstate)")
    hf = _body(INTERCHANGE.read_text(), "mrsf_nac_rohf_hf_adjoint(infos)")
    assert "public :: fock_deriv_contract_os2" in fock
    assert "type, extends(grd2_compute_data_t) :: grd2_fockprobe_os2_data_t" in fock
    for body in (gradient, hf):
        assert body.count("call fock_deriv_contract_os2(") == 1
        assert "call fock_deriv_contract_os(" not in body

