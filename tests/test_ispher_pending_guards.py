from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_ecp_spherical_path_has_clear_guard():
    text = (ROOT / "source/ecp.F90").read_text()
    assert "subroutine guard_spherical_ecp" in text
    assert "ispher=false for ECP basis sets" in text
    assert text.count("call guard_spherical_ecp(basis)") >= 4


def test_hessian_spherical_path_has_python_and_native_guards():
    py = (ROOT / "pyoqp/oqp/library/single_point.py").read_text()
    f90 = (ROOT / "source/modules/hf_hessian.F90").read_text()
    assert "def _spherical_ao_active" in py
    assert "Analytic HF/DFT Hessian with spherical-harmonic AO dimensions" in py
    assert "Analytic HF/DFT Hessian with spherical-harmonic AO dimensions" in f90


def test_ekt_uses_full_matrix_space_no_ispher_guard():
    text = (ROOT / "source/modules/tdhf_mrsf_ekt.F90").read_text()
    assert "OQP_mrsf_ekt_density_mo" in text
    assert "HARMONIC_ACTIVE" not in text


def test_mrsf_gradient_spherical_path_is_not_guarded():
    text = (ROOT / "source/modules/tdhf_mrsf_gradient.F90").read_text()
    assert "grd2_mrsf_build_cart" in text
    assert "call gcomp%build_cart(basis)" in text
    assert "MRSF-TDDFT gradients with spherical-harmonic AO dimensions" not in text
