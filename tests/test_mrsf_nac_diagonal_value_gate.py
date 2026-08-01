"""Pure orchestration contracts for the full Lee diagonal-value gate."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
GATE_PATH = ROOT / "tools" / "nac_lagrangian" / "diagonal_value_gate.py"
SPEC = importlib.util.spec_from_file_location("diagonal_value_gate", GATE_PATH)
assert SPEC is not None and SPEC.loader is not None
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)


def test_duplicate_slots_preserve_the_selected_amplitude():
    raw = np.arange(15.0).reshape(3, 5)
    duplicated, amplitude = GATE.duplicate_state_slots(raw, 3, 5, 0, 2)
    np.testing.assert_array_equal(amplitude, np.arange(5.0))
    np.testing.assert_array_equal(duplicated.reshape(3, 5)[0], amplitude)
    np.testing.assert_array_equal(duplicated.reshape(3, 5)[2], amplitude)


def test_duplicate_slots_reject_the_literal_diagonal():
    with pytest.raises(ValueError, match="must differ"):
        GATE.duplicate_state_slots(np.zeros((2, 3)), 2, 3, 1, 1)


def test_gate_requires_native_and_legacy_multiplier_identity():
    source = GATE_PATH.read_text()
    for token in (
        'mol.data["OQP::nac_zvec_solution"]',
        'mol.data["OQP::nac_rohf_solution"]',
        "native + 0.5*legacy",
        "zvector_convention_closure",
        "pair - excitation",
    ):
        assert token in source
    fortran = (
        ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90"
    ).read_text()
    assert "OQP::nac_zvec_solution" in fortran
    assert "zeta_native + xk_legacy/2 = 0" in fortran


def test_rohf_multiplier_comparison_is_flat_blockwise_and_sign_aware():
    slices = GATE.rohf_rotation_slices(nbf=6, noca=3, nocb=1)
    assert list(slices) == ["socc-docc", "virt-docc", "virt-socc"]
    assert slices["socc-docc"] == slice(0, 2)
    assert slices["virt-docc"] == slice(2, 5)
    assert slices["virt-socc"] == slice(5, 11)

    xk = np.arange(1.0, 12.0)
    zeta = -0.5*xk
    np.testing.assert_array_equal(
        GATE.zvector_convention_closure(zeta, xk, 6, 3, 1),
        np.zeros_like(xk),
    )


def test_rohf_multiplier_comparison_does_not_hide_shape_or_length_errors():
    with pytest.raises(RuntimeError, match="one-dimensional"):
        GATE.require_rohf_vector(
            np.zeros((1, 11)), "test", nbf=6, noca=3, nocb=1
        )
    with pytest.raises(RuntimeError, match="expected 11"):
        GATE.require_rohf_vector(
            np.zeros(10), "test", nbf=6, noca=3, nocb=1
        )
    bad = np.zeros(11)
    bad[5] = np.nan
    with pytest.raises(RuntimeError, match="non-finite"):
        GATE.require_rohf_vector(bad, "test", nbf=6, noca=3, nocb=1)


def test_fortran_export_is_one_dimensional_and_precedes_density_build():
    source = (
        ROOT / "source" / "modules" / "tdhf_mrsf_z_vector.F90"
    ).read_text()
    audit = source.split("! Diagonal Lee-limit audit:", 1)[1].split(
        "call build_mrsf_relaxed_density_and_w()", 1
    )[0]
    assert "ta_type_real64, lzdim, (/ lzdim /)" in audit
    assert "xk_dump = xk" in audit
    assert "flat lzdim vector in the SD/DV/SV loop ordering" in audit
    assert "sfropcal inserts xk/2" in audit
    assert "rohf_unpack_trial inserts zeta directly" in audit


def test_tight_input_selects_a_deterministic_legacy_zvector_root():
    source = (
        ROOT / "tools" / "nac_lagrangian" /
        "H2O_energy_tlf0_tight_zv_analytic.inp"
    ).read_text()
    for setting in (
        "runtype=energy",
        "conv=1.0e-10",
        "target=1",
        "z_solver=2",
        "zvconv=1.0e-10",
        "maxit_zv=500",
    ):
        assert setting in source
