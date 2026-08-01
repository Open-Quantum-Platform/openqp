"""Pure contracts for the resident exact-tlf gamma validation gate."""

from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
GATE_DIR = ROOT / "tools" / "nac_lagrangian"
sys.path.insert(0, str(GATE_DIR))
import gamma_gate as GATE  # noqa: E402
import nac_formula_kernel as FORMULA  # noqa: E402


def test_resident_tagarray_decoder_preserves_fortran_pair_and_orbital_order():
    nstate, nbf = 3, 4
    flat = np.empty(nstate*nstate*nbf*nbf)
    for istate in range(nstate):
        for jstate in range(nstate):
            for q in range(nbf):
                for p in range(nbf):
                    offset = (
                        (istate + jstate*nstate)*nbf*nbf + p + q*nbf
                    )
                    flat[offset] = 1000*istate + 100*jstate + 10*p + q

    decoded = GATE.decode_resident_gamma(flat, nstate, nbf)
    for istate in range(nstate):
        for jstate in range(nstate):
            for p in range(nbf):
                for q in range(nbf):
                    assert decoded[istate, jstate, p, q] == (
                        1000*istate + 100*jstate + 10*p + q
                    )


def test_all_ordered_state_directions_and_same_space_blocks_are_explicit():
    assert GATE.ordered_state_pairs(3) == [
        (1, 0), (2, 0), (0, 1), (2, 1), (0, 2), (1, 2)
    ]
    blocks = GATE.generator_blocks(nocb=2, noca=4, nbf=6)
    assert tuple(blocks) == ("dd", "ds", "dv", "ss", "sv", "vv")
    assert blocks["dd"] == ((1, 0),)
    assert blocks["ss"] == ((3, 2),)
    assert blocks["vv"] == ((5, 4),)


def _independent_ordered_tensors():
    nstate, nbf = 2, 6
    resident = np.zeros((nstate, nstate, nbf, nbf))
    # Deliberately unrelated reverse-state matrices: state antisymmetry must
    # neither be assumed nor manufactured by the gate.
    for p, q, value in ((2, 1, 1.5), (4, 0, -0.25), (5, 4, 0.75)):
        resident[0, 1, p, q] = value
        resident[0, 1, q, p] = -value
    for p, q, value in ((2, 1, -0.4), (4, 0, 2.0), (5, 4, -1.25)):
        resident[1, 0, p, q] = value
        resident[1, 0, q, p] = -value
    cofactor = resident.copy()
    derivative = np.zeros_like(resident)
    for istate, jstate in GATE.ordered_state_pairs(nstate):
        for q in range(nbf - 1):
            for p in range(q + 1, nbf):
                derivative[istate, jstate, p, q] = (
                    2.0*resident[istate, jstate, p, q]
                )
    return resident, cofactor, derivative


def test_analysis_does_not_force_state_or_orbital_antisymmetry():
    resident, cofactor, derivative = _independent_ordered_tensors()
    report = GATE.analyze_gamma(
        resident, cofactor, derivative, nocb=1, noca=3
    )
    for key in (
        "cofactor_pair_max_abs",
        "resident_generator_pair_max_abs",
        "cofactor_generator_pair_max_abs",
        "orbital_antisym_pair_max_abs",
    ):
        np.testing.assert_array_equal(report[key], 0.0)

    # Break only reverse pair (1,0), in the same-space ss block.  A state
    # projection would leak/cancel this error; the raw ordered comparison must
    # localize it to exactly that direction and block.
    resident[1, 0, 2, 1] += 0.125
    broken = GATE.analyze_gamma(
        resident, cofactor, derivative, nocb=1, noca=3
    )
    pair_index = list(map(tuple, broken["ordered_pairs"])).index((1, 0))
    other_index = list(map(tuple, broken["ordered_pairs"])).index((0, 1))
    ss_index = list(broken["block_names"]).index("ss")
    assert broken["cofactor_pair_max_abs"][pair_index] == pytest.approx(0.125)
    assert broken["cofactor_pair_max_abs"][other_index] == 0.0
    assert broken["cofactor_block_max_abs"][pair_index, ss_index] == pytest.approx(0.125)
    assert broken["orbital_antisym_pair_max_abs"][pair_index] == pytest.approx(0.125)


def test_generator_coordinate_and_production_separation_contracts():
    generator = FORMULA.antisymmetric_generator(5, p=4, q=1)
    assert generator[4, 1] == 1.0
    assert generator[1, 4] == -1.0
    np.testing.assert_array_equal(generator + generator.T, 0.0)
    with pytest.raises(ValueError, match="q < p"):
        FORMULA.antisymmetric_generator(5, p=1, q=4)
    rotation = FORMULA.exact_generator_rotation(5, p=4, q=1, theta=0.37)
    np.testing.assert_allclose(rotation.T @ rotation, np.eye(5), atol=1.0e-15)
    assert np.linalg.det(rotation) == pytest.approx(1.0)

    gate_source = (GATE_DIR / "gamma_gate.py").read_text()
    assert "oqp.mrsf_nac_metric_data(mol)" in gate_source
    assert "FORMULA.gamma_closed(context)" in gate_source
    assert "FORMULA.generator_derivative_sweep" in gate_source
    assert "no state antisymmetry is imposed" in gate_source
    production = (
        ROOT / "pyoqp" / "oqp" / "library" / "nac_analytic.py"
    ).read_text()
    assert "gamma_gate" not in production
    assert "nac_formula_kernel" not in production


def test_cofactor_and_exact_generator_oracles_agree_without_compiled_openqp():
    rng = np.random.default_rng(2)
    context = {
        "nstate": 2,
        "noca": 3,
        "nocb": 1,
        "nbf": 5,
        "nvirb": 4,
        "nij": 12,
        "noc": 2,
        "RS": 1.0/np.sqrt(2.0),
        "Xt": [rng.normal(size=(3, 4)) for _ in range(2)],
        "genmask": np.ones((3, 4)),
    }
    context["genmask"][1:3, 0:2] = 0.0
    cofactor = FORMULA.gamma_closed(context)
    derivative = FORMULA.generator_derivative_sweep(context)
    errors = []
    for q in range(context["nbf"] - 1):
        for p in range(q + 1, context["nbf"]):
            errors.append(np.max(np.abs(
                2.0*cofactor[:, :, p, q] - derivative[:, :, p, q]
            )))
    assert max(errors) <= 2.0e-10
