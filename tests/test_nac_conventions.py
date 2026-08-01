"""Convention regressions for overlap-based numerical NAC assembly."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest


MODULE = (
    Path(__file__).resolve().parents[1]
    / "pyoqp" / "oqp" / "library" / "nac_utils.py"
)
SPEC = importlib.util.spec_from_file_location("nac_utils_under_test", MODULE)
assert SPEC is not None and SPEC.loader is not None
NAC_UTILS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NAC_UTILS)


def test_native_tag_transpose_and_hst_factor_recover_signed_coupling():
    step = 2.5e-4
    expected = np.array(
        [[0.0, 0.7, -0.2], [-0.7, 0.0, 0.4], [0.2, -0.4, 0.0]]
    )
    canonical = np.eye(3) + step * expected
    native_tag = canonical.T

    overlap = NAC_UTILS.canonical_state_overlap(native_tag)
    coupling = NAC_UTILS.hst_derivative_coupling(overlap, step)

    assert np.allclose(overlap, canonical)
    assert np.allclose(coupling, expected)
    # Omitting the storage-boundary transpose gives the genuine sign defect
    # that PR #160 must prevent, not a removable factor-of-two convention.
    wrong_sign = NAC_UTILS.hst_derivative_coupling(native_tag, step)
    assert np.allclose(wrong_sign, -expected)


def test_interstate_coupling_uses_ej_minus_ei_for_vectors():
    energies = np.array([0.10, 0.25, 0.70])
    derivative = np.zeros((3, 3, 2, 3))
    derivative[0, 1, 0, 2] = 0.8
    derivative[1, 0, 0, 2] = -0.8

    interstate = NAC_UTILS.interstate_coupling(derivative, energies)

    assert interstate[0, 1, 0, 2] == pytest.approx((0.25 - 0.10) * 0.8)
    assert interstate[1, 0, 0, 2] == pytest.approx(interstate[0, 1, 0, 2])
    assert np.allclose(interstate, np.swapaxes(interstate, 0, 1))


def test_state_gauge_covariance_is_joint_not_pairwise():
    energies = np.array([0.10, 0.25, 0.70])
    derivative = np.array(
        [[0.0, 0.7, -0.2], [-0.7, 0.0, 0.4], [0.2, -0.4, 0.0]]
    )
    signs = np.array([1.0, -1.0, -1.0])
    gauge = signs[:, None] * signs[None, :]

    interstate = NAC_UTILS.interstate_coupling(derivative, energies)
    gauged = NAC_UTILS.interstate_coupling(derivative * gauge, energies)

    assert np.allclose(gauged, interstate * gauge)
    assert np.prod([gauge[0, 1], gauge[0, 2], gauge[1, 2]]) == 1.0


@pytest.mark.parametrize("shape", [(2, 3), (2,), (2, 2, 1)])
def test_canonical_overlap_rejects_non_square_matrices(shape):
    with pytest.raises(ValueError, match="square matrix"):
        NAC_UTILS.canonical_state_overlap(np.zeros(shape))
