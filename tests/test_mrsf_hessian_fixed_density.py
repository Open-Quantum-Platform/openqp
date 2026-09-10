"""Regression gates for the fixed relaxed-density MRSF Hessian skeleton."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_fixed_density.F90"
GRADIENT = ROOT / "source/modules/tdhf_sf_gradient.F90"


def _compact(path):
    return "".join(path.read_text().lower().split())


def test_fixed_density_hessian_differentiates_every_gradient_integral_family():
    source = _compact(SOURCE)
    for call in (
        "callhess_nn",
        "callhess_ee_overlap",
        "callhess_ee_kinetic",
        "callhess_en",
        "calladd_ecphess",
        "callgrd2_hess_driver",
    ):
        assert call in source
    assert "grd2_mrsf_compute_data_t" in source
    assert "oqp_td_mrsf_density" in source
    assert "h_fixed=0.5_dp*(h_fixed+transpose(h_fixed))" in source


def test_overlap_and_one_particle_density_conventions_match_mrsf_gradient():
    source = _compact(SOURCE)
    gradient = _compact(GRADIENT)
    assert "overlap_density=overlap_density+2.0_dp*wao" in source
    assert "dens=dens+2*w" in gradient
    assert "one_density=dmat_a+dmat_b+td_p(:,1)+td_p(:,2)" in source


def test_fixed_density_layer_does_not_claim_missing_xc_or_response_rows():
    text = SOURCE.read_text().lower()
    assert "semilocal xc quadrature derivatives" in text
    assert "first-order" in text
    for forbidden in ("slater", "fock_space", "fci_hamiltonian"):
        assert forbidden not in text
