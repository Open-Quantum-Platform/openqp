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


def test_auto_spherical_basis_is_limited_to_implemented_transforms():
    set_basis = (ROOT / "pyoqp/oqp/library/set_basis.py").read_text()
    cart2sph = (ROOT / "source/integrals/cart2sph.F90").read_text()

    assert "ang_mom > 4" in set_basis
    assert "input.ispher=True only supports pure spherical shells through g" in set_basis
    assert "shell_harmonic = 1 if shell_is_spherical and ang_mom <= 4 else 0" in set_basis
    assert "pure spherical transforms are implemented only through g shells" in cart2sph


def test_cartesian_override_keeps_cartesian_normalization():
    text = (ROOT / "source/basis_tools.F90").read_text()

    assert "if (HARMONIC_ACTIVE .and. basis%harmonic(i) == 1) then" in text
