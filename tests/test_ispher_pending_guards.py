from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_ecp_spherical_path_uses_cartesian_to_spherical_transform():
    text = (ROOT / "source/ecp.F90").read_text()
    c2s = (ROOT / "source/integrals/cart2sph.F90").read_text()
    assert "subroutine transform_ecp_matrix" in text
    assert "cart2sph_mat_unit" in text
    assert "ECP integrals with spherical-harmonic AO dimensions" not in text
    assert "subroutine cart2sph_mat_unit" in c2s


def test_ground_state_hessian_spherical_path_is_not_guarded():
    py = (ROOT / "pyoqp/oqp/library/single_point.py").read_text()
    f90 = (ROOT / "source/modules/hf_hessian.F90").read_text()
    grd1 = (ROOT / "source/integrals/grd1.F90").read_text()
    assert "def _spherical_ao_active" in py
    assert "Analytic HF/DFT Hessian with spherical-harmonic AO dimensions" not in py
    assert "Analytic HF/DFT Hessian with spherical-harmonic AO dimensions" not in f90
    assert "reduce_der1_shell_block" in grd1
    assert "call gcomp%build_cart(basis)" in f90
    assert "TDDFT analytic Hessian is not implemented yet" in py
    assert "SF-TDDFT analytic Hessian is not implemented yet" in py


def test_ekt_uses_full_matrix_space_no_ispher_guard():
    text = (ROOT / "source/modules/tdhf_mrsf_ekt.F90").read_text()
    assert "OQP_mrsf_ekt_density_mo" in text
    assert "HARMONIC_ACTIVE" not in text


def test_mrsf_gradient_spherical_path_is_not_guarded():
    text = (ROOT / "source/modules/tdhf_mrsf_gradient.F90").read_text()
    assert "grd2_mrsf_build_cart" in text
    assert "call gcomp%build_cart(basis)" in text
    assert "MRSF-TDDFT gradients with spherical-harmonic AO dimensions" not in text


def test_pcm_ddx_sources_are_present_after_main_sync():
    solvent = (ROOT / "source/solvent_pcm.F90").read_text()
    schema = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()
    checker = (ROOT / "pyoqp/oqp/utils/input_checker.py").read_text()
    assert "subroutine add_pcm_reaction_field" in solvent
    assert "oqp_ddx_pcm_solve" in solvent
    assert "'pcm':" in schema
    assert "single-point energies only" in checker
