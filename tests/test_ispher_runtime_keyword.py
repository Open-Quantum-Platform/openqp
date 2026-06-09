from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_ispher_input_schema_defaults_enabled():
    text = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()

    assert "'ispher': {'type': bool, 'default': 'True'}" in text
    assert "oqp_set_harmonic_active" in text


def test_ispher_cffi_binding_is_declared():
    header = (ROOT / "include/oqp.h").read_text()
    constants = (ROOT / "source/constants.F90").read_text()

    assert "void oqp_set_harmonic_active(bool flag);" in header
    assert 'bind(C, name="oqp_set_harmonic_active")' in constants
    assert "logical :: HARMONIC_ACTIVE = .true." in constants


def test_apply_basis_reports_runtime_cartesian_override():
    text = (ROOT / "source/modules/apply_basis.F90").read_text()

    assert "use constants, only: HARMONIC_ACTIVE, NUM_CART_BF" in text
    assert "if (HARMONIC_ACTIVE .and. nsphsh > 0) then" in text
    assert "AO angular type: Cartesian (6d/10f/15g)" in text
