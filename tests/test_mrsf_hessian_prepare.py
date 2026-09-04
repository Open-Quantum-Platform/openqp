from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_prepare.F90"
RESPONSE = ROOT / "source/modules/tdhf_hessian_response.F90"


def test_prepare_layer_connects_complete_first_response():
    text = SOURCE.read_text()
    for call in (
        "build_rohf_nuclear_response",
        "build_orbital_density_derivatives",
        "mrsf_xc_fock_total_derivative",
        "solve_mrsf_first_nuclear_response",
    ):
        assert call in text
    for result in ("dmo_a", "dmo_b", "dfock_a", "dfock_b", "dax", "dx"):
        assert result in text
    assert "require_two_somo=.true." in text


def test_prepare_layer_fails_closed_outside_initial_scope():
    compact = "".join(SOURCE.read_text().lower().split())
    for condition in (
        "infos%control%scftype/=3",
        "nocca-noccb/=2",
        "infos%tddft%umrsf",
        "infos%dft%cam_flag",
        "infos%functional%needtau",
        "infos%functional%needlapl",
    ):
        assert condition in compact
    assert "hiroya" in compact and "nakata" in compact


def test_eigenpair_criterion_respects_double_precision_and_solver_resolution():
    compact = "".join(RESPONSE.read_text().lower().split()).replace("&", "")
    threshold = (
        "eigenpair_tol=max(1.0e-8_dp,100.0_dp*solve_tol,"
        "sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(omega)))"
    )
    assert compact.count(threshold) == 2
    assert "residual_max>eigenpair_tol" in compact
