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


def test_auto_spherical_basis_uses_all_supported_transforms():
    set_basis = (ROOT / "pyoqp/oqp/library/set_basis.py").read_text()
    cart2sph = (ROOT / "source/integrals/cart2sph.F90").read_text()
    symmetry = (ROOT / "pyoqp/oqp/library/symmetry.py").read_text()

    assert "shell_harmonic = 1 if shell_is_spherical else 0" in set_basis
    assert "input.ispher=True only supports pure spherical shells through g" not in set_basis
    assert "c2s_build(l)" in cart2sph
    assert "do l = 2, BAS_MXANG" in cart2sph
    assert "pure spherical transforms are implemented only through g shells" not in cart2sph
    assert "5: [(5, 0, 0)" in symmetry
    assert "6: [(6, 0, 0)" in symmetry


def test_cartesian_override_keeps_cartesian_normalization():
    text = (ROOT / "source/basis_tools.F90").read_text()

    assert "if (HARMONIC_ACTIVE .and. basis%harmonic(i) == 1) then" in text


def test_two_e_backends_reduce_pure_blocks_before_consumers():
    int2 = (ROOT / "source/integrals/int2.F90").read_text()
    rot = (ROOT / "source/integrals/int_rotaxis.F90").read_text()
    rys = (ROOT / "source/integrals/int_rys.F90").read_text()

    assert "public genr22_reduce_pure" in rot
    assert "subroutine genr22_reduce_pure" in rot
    assert "public :: int2_rys_reduce_pure" in rys
    assert "subroutine int2_rys_reduce_pure" in rys
    assert "call genr22_reduce_pure(basis, eri_data%ids, eri_data%flips" in int2
    assert "call int2_rys_reduce_pure(basis, eri_data%gdat, eri_data%ints, nbf)" in int2
    assert "call genr22_reduce_pure(basis, shell_ids, flips, ints, nbf)" in int2
    assert "call int2_rys_reduce_pure(basis, gdat, ints, nbf)" in int2
    assert "call normalize_ints(nbf, am(flips), pints)" in int2
    assert "nbf = NUM_CART_BF(am(flips))" in int2
