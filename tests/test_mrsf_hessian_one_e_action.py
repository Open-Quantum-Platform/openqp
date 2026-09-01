"""Independent checks for the MRSF one-electron operator derivative."""

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_one_e_action.F90"
TOTAL_ACTION = ROOT / "source/modules/tdhf_mrsf_hessian_operator_derivative.F90"
TWO_E_ACTION = ROOT / "source/modules/tdhf_mrsf_hessian_two_e_action.F90"
STATE_RESPONSE = ROOT / "source/modules/tdhf_mrsf_hessian_state_response.F90"


def _analytic_mo_derivative(c, dc, f, df):
    return dc.T @ f @ c + c.T @ df @ c + c.T @ f @ dc


def test_moving_basis_mo_operator_derivative_matches_finite_difference():
    rng = np.random.default_rng(20260901)
    c, _ = np.linalg.qr(rng.normal(size=(5, 5)))
    dc = rng.normal(size=(5, 5))
    f = rng.normal(size=(5, 5))
    f = 0.5 * (f + f.T)
    df = rng.normal(size=(5, 5))
    df = 0.5 * (df + df.T)
    expected = _analytic_mo_derivative(c, dc, f, df)

    errors = []
    for step in (2.0e-3, 1.0e-3, 5.0e-4, 2.5e-4):
        plus = (c + step * dc).T @ (f + step * df) @ (c + step * dc)
        minus = (c - step * dc).T @ (f - step * df) @ (c - step * dc)
        finite_difference = (plus - minus) / (2.0 * step)
        errors.append(np.max(np.abs(finite_difference - expected)))

    assert errors[1] < errors[0] / 3.9
    assert errors[2] < errors[1] / 3.9
    assert errors[3] < errors[2] / 3.9
    assert errors[3] < 1.0e-6


def test_production_action_uses_linear_mrsfesum_boundary():
    compact = "".join(SOURCE.read_text().lower().split())
    assert "calltransform_mo_operator_derivative" in compact
    assert "callmrsfesum" in compact
    assert "dfock_a_ao" in compact
    assert "dfock_b_ao" in compact
    for forbidden in ("slater", "fock_space", "fci_hamiltonian"):
        assert forbidden not in compact


def test_total_operator_derivative_combines_one_and_seven_density_actions():
    compact = "".join(TOTAL_ACTION.read_text().lower().split())
    assert "callprepare_mrsf_response_scaling" in compact
    assert "callbuild_mrsf_one_e_operator_derivative_action" in compact
    assert "callbuild_mrsf_two_e_operator_derivative_action" in compact
    assert "dax=one_e+two_e" in compact


def test_range_separated_derivative_eri_action_fails_closed():
    compact = "".join(TWO_E_ACTION.read_text().lower().split())
    assert "if(infos%control%hamilton==20.and.infos%dft%cam_flag)then" in compact
    assert "status=-2" in compact


def test_first_nuclear_response_connects_orbitals_fock_operator_and_amplitudes():
    compact = "".join(STATE_RESPONSE.read_text().lower().split())
    for call in (
        "callbuild_orbital_density_derivatives",
        "callbuild_rohf_total_fock_derivatives",
        "callbuild_mrsf_operator_derivative_action",
        "callsolve_mrsf_tda_amplitude_derivatives",
    ):
        assert call in compact
    assert "numericalnucleardisplacement" in compact
