"""Regression gates for MRSF Z/W Hessian intermediates."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source/modules/tdhf_mrsf_hessian_intermediates.F90"


def test_apb_explicit_eri_derivative_uses_symmetric_density_convention():
    text = SOURCE.read_text()
    assert "int_apb=.true." in text
    assert "dfa(:,:,coordinate)=2.0_dp*explicit_a" in text
    assert "dfb(:,:,coordinate)=2.0_dp*explicit_b" in text
    assert "ordinary open-shell Fock derivative is therefore doubled" in text


def test_reconstructed_gradient_intermediates_fail_closed():
    text = SOURCE.read_text()
    assert "max(baseline_error_a,baseline_error_b)>1.0e-8_dp" in text
    assert "data%z_rhs_hxa-(stored_hxa-2.0_dp*" in text
    assert "data%ppija-stored_ppija" in text
    assert "data%ppijb-stored_ppijb" in text


def test_no_unconditional_debug_files_or_diagnostic_stderr():
    text = SOURCE.read_text()
    assert "/private/tmp" not in text
    assert "error_unit" not in text
