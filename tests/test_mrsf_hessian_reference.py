"""Spin-adapted finite-model checks for the MRSF-TDDFT Hessian algebra.

The production and reference representations in this file are the physical
CO/OV/CV/OO response amplitudes. The two M_S parent contributions enter only
through their normalized exchange-symmetry combinations. The seven response
density channels are kept explicit so that the ordinary channel and the six
spin-pairing channels cannot be mixed accidentally.
"""

from __future__ import annotations

import math
import unittest

import numpy as np


_TOL = 1.0e-12
_SEVEN_DENSITY_TOPOLOGY = (
    ("CO", "+1"),
    ("OV", "+1"),
    ("CV", "+1"),
    ("OO", "both"),
    ("CV", "-1"),
    ("OV", "-1"),
    ("CO", "-1"),
)

# Minimal two-SOMO response space: C, O1, O2, V. The packed order is the
# order consumed by iatogen/mrsfcbc when nbf=4, nocca=3, and noccb=1.
_PACKED_LABELS = (
    "CO1", "L", "G", "CO2", "D", "R", "CV", "OV1", "OV2",
)
_SINGLET_POSITIONS = (0, 1, 2, 3, 4, 6, 7, 8)
_TRIPLET_POSITIONS = (0, 1, 3, 6, 7, 8)


def _selection(positions):
    transform = np.zeros((len(_PACKED_LABELS), len(positions)))
    for column, position in enumerate(positions):
        transform[position, column] = 1.0
    return transform


def _sector(label):
    if label.startswith("CO"):
        return "CO"
    if label.startswith("OV"):
        return "OV"
    if label == "CV":
        return "CV"
    return "OO"


def _symmetric_coefficients(size, seed=83):
    """Return A0, A0', and A0'' in a physical spin-adapted basis."""

    rng = np.random.default_rng(seed + size)
    result = []
    for order, scale in enumerate((0.32, 0.075, 0.028)):
        raw = rng.normal(scale=scale, size=(size, size))
        matrix = 0.5 * (raw + raw.T)
        if order == 0:
            matrix += np.diag(np.linspace(-0.7, 1.4, size))
        result.append(matrix)
    return result


def _spc_coefficients(labels):
    """Return C_SPC, C_SPC', and C_SPC'' in the allowed sectors only."""

    result = []
    scales = ((0.11, -0.07, 0.05), (0.035, 0.021, -0.017),
              (-0.012, 0.009, 0.006))
    for coco, ovov, coov in scales:
        matrix = np.zeros((len(labels), len(labels)))
        for i, left in enumerate(labels):
            for j, right in enumerate(labels):
                left_sector = _sector(left)
                right_sector = _sector(right)
                if left_sector == right_sector == "CO":
                    matrix[i, j] = coco / (1.0 + abs(i - j))
                elif left_sector == right_sector == "OV":
                    matrix[i, j] = ovov / (1.0 + abs(i - j))
                elif {left_sector, right_sector} == {"CO", "OV"}:
                    matrix[i, j] = coov / (1.0 + abs(i - j))
        result.append(matrix)
    return result


def _raw_seven_density_action(a0, c_spc, vector, multiplicity):
    """Apply the OpenQP raw-channel convention in a spin-adapted basis."""

    raw_special = -(c_spc @ vector)
    raw_multiplier = {1: 1.0, 3: -1.0}[multiplicity]
    return a0 @ vector + raw_multiplier * raw_special


def _physical_action(a0, c_spc, vector, multiplicity):
    physical_sign = {1: -1.0, 3: 1.0}[multiplicity]
    return (a0 + physical_sign * c_spc) @ vector


def _eigen_response(a0, a1, a2, root=0):
    values, vectors = np.linalg.eigh(a0)
    value = values[root]
    vector = vectors[:, root]
    derivative = float(vector @ a1 @ vector)
    vector_derivative = np.zeros_like(vector)
    for other in range(len(values)):
        if other == root:
            continue
        gap = value - values[other]
        if abs(gap) < 1.0e-8:
            raise ValueError("isolated-state response requested for a degenerate root")
        coupling = float(vectors[:, other] @ a1 @ vector)
        vector_derivative += vectors[:, other] * (coupling / gap)
    second = float(vector @ a2 @ vector + 2.0 * vector @ a1 @ vector_derivative)
    residual = np.max(np.abs(
        (a0 - value * np.eye(a0.shape[0])) @ vector_derivative
        + (a1 - derivative * np.eye(a0.shape[0])) @ vector
    ))
    return value, vector, derivative, second, residual


def _tracked_value(matrix, reference):
    values, vectors = np.linalg.eigh(matrix)
    overlaps = np.abs(vectors.T @ reference)
    root = int(np.argmax(overlaps))
    return float(values[root]), float(overlaps[root])


class TestMRSFSpinAdaptedHessianReference(unittest.TestCase):
    def test_seven_density_topology_is_complete_and_ordered(self):
        self.assertEqual(
            _SEVEN_DENSITY_TOPOLOGY,
            (("CO", "+1"), ("OV", "+1"), ("CV", "+1"),
             ("OO", "both"),
             ("CV", "-1"), ("OV", "-1"), ("CO", "-1")),
        )
        self.assertEqual(len(_SEVEN_DENSITY_TOPOLOGY), 7)

    def test_physical_packed_spaces_keep_each_oo_amplitude_once(self):
        singlet = _selection(_SINGLET_POSITIONS)
        triplet = _selection(_TRIPLET_POSITIONS)
        np.testing.assert_allclose(singlet.T @ singlet, np.eye(8),
                                   atol=_TOL, rtol=0.0)
        np.testing.assert_allclose(triplet.T @ triplet, np.eye(6),
                                   atol=_TOL, rtol=0.0)
        self.assertEqual(tuple(_PACKED_LABELS[i] for i in _SINGLET_POSITIONS),
                         ("CO1", "L", "G", "CO2", "D", "CV", "OV1", "OV2"))
        self.assertEqual(tuple(_PACKED_LABELS[i] for i in _TRIPLET_POSITIONS),
                         ("CO1", "L", "CO2", "CV", "OV1", "OV2"))
        self.assertTrue(np.all(singlet[5] == 0.0))
        self.assertTrue(np.all(triplet[[2, 4, 5]] == 0.0))

    def test_parent_exchange_gives_authoritative_singlet_triplet_phases(self):
        plus_parent = np.array([1.0, 0.0])
        minus_parent = np.array([0.0, 1.0])
        singlet = (plus_parent - minus_parent) / math.sqrt(2.0)
        triplet = (plus_parent + minus_parent) / math.sqrt(2.0)
        exchange = np.array([[0.0, 1.0], [1.0, 0.0]])
        np.testing.assert_allclose(exchange @ singlet, -singlet,
                                   atol=_TOL, rtol=0.0)
        np.testing.assert_allclose(exchange @ triplet, triplet,
                                   atol=_TOL, rtol=0.0)
        self.assertAlmostEqual(float(singlet @ triplet), 0.0, places=14)
        self.assertAlmostEqual(float(singlet @ singlet), 1.0, places=14)
        self.assertAlmostEqual(float(triplet @ triplet), 1.0, places=14)

        # Removing either parent contribution cannot preserve normalization or
        # the required parent-exchange symmetry.
        one_parent = plus_parent / math.sqrt(2.0)
        self.assertAlmostEqual(float(one_parent @ one_parent), 0.5, places=14)
        self.assertGreater(np.linalg.norm(exchange @ one_parent - one_parent), 0.5)

    def test_raw_seven_density_sign_matches_physical_spc_sign_columnwise(self):
        for multiplicity, positions in ((1, _SINGLET_POSITIONS),
                                        (3, _TRIPLET_POSITIONS)):
            labels = tuple(_PACKED_LABELS[i] for i in positions)
            a0 = _symmetric_coefficients(len(labels))[0]
            c_spc = _spc_coefficients(labels)[0]
            for column in range(len(labels)):
                unit = np.eye(len(labels))[:, column]
                np.testing.assert_allclose(
                    _raw_seven_density_action(a0, c_spc, unit, multiplicity),
                    _physical_action(a0, c_spc, unit, multiplicity),
                    atol=_TOL,
                    rtol=0.0,
                )

    def test_spc_is_confined_to_coco_coov_and_ovov(self):
        for positions in (_SINGLET_POSITIONS, _TRIPLET_POSITIONS):
            labels = tuple(_PACKED_LABELS[i] for i in positions)
            for matrix in _spc_coefficients(labels):
                for i, left in enumerate(labels):
                    for j, right in enumerate(labels):
                        sectors = {_sector(left), _sector(right)}
                        allowed = sectors in ({"CO"}, {"OV"}, {"CO", "OV"})
                        if not allowed:
                            self.assertEqual(matrix[i, j], 0.0)

    def test_spin_adapted_eigenvalue_hessian_matches_finite_difference(self):
        for multiplicity, positions in ((1, _SINGLET_POSITIONS),
                                        (3, _TRIPLET_POSITIONS)):
            labels = tuple(_PACKED_LABELS[i] for i in positions)
            a_coeff = _symmetric_coefficients(len(labels))
            c_coeff = _spc_coefficients(labels)
            sign = {1: -1.0, 3: 1.0}[multiplicity]
            a0, a1, a2 = [a + sign * c for a, c in zip(a_coeff, c_coeff)]
            value, vector, _first, second, residual = _eigen_response(a0, a1, a2)
            self.assertLess(residual, 1.0e-10)
            errors = []
            for step in (8.0e-3, 4.0e-3, 2.0e-3):
                plus = a0 + step * a1 + 0.5 * step * step * a2
                minus = a0 - step * a1 + 0.5 * step * step * a2
                eplus, overlap_plus = _tracked_value(plus, vector)
                eminus, overlap_minus = _tracked_value(minus, vector)
                self.assertGreater(min(overlap_plus, overlap_minus), 0.99)
                finite_difference = (eplus - 2.0 * value + eminus) / step**2
                errors.append(abs(finite_difference - second))
            self.assertLess(errors[-1], 2.0e-6, (multiplicity, errors))
            self.assertLess(errors[-1], errors[0] / 3.0)

    def test_mrsf_zero_has_no_singlet_triplet_spc_splitting(self):
        size = len(_TRIPLET_POSITIONS)
        a0 = _symmetric_coefficients(size)[0]
        vector = np.linspace(-0.4, 0.7, size)
        zero = np.zeros_like(a0)
        np.testing.assert_allclose(
            _raw_seven_density_action(a0, zero, vector, 1),
            _raw_seven_density_action(a0, zero, vector, 3),
            atol=_TOL,
            rtol=0.0,
        )


if __name__ == "__main__":
    unittest.main()
