"""Independent second-quantized audit of an MRSF nuclear Hessian.

The production method is *not* represented or solved in a determinant space;
it uses the spin-adapted CO/OV/CV/OO response amplitudes and seven transition
densities.  This deliberately isolated tiny-system audit expands those
spin-adapted functions only to compare a dense second-quantized Hamiltonian
with separate Slater--Condon expressions.  No OpenQP response routine imports
or consumes this representation.
"""

from __future__ import annotations

import itertools
import math
import unittest

import numpy as np


_TOL = 1.0e-11


def _parity(det: int, orbital: int) -> int:
    return -1 if (det & ((1 << orbital) - 1)).bit_count() % 2 else 1


def _annihilate(det: int, orbital: int):
    if not det & (1 << orbital):
        return None
    return det ^ (1 << orbital), _parity(det, orbital)


def _create(det: int, orbital: int):
    if det & (1 << orbital):
        return None
    return det | (1 << orbital), _parity(det, orbital)


def _apply(det: int, operations):
    phase = 1
    for operation, orbital in operations:
        result = (_annihilate if operation == "ann" else _create)(det, orbital)
        if result is None:
            return None
        det, local_phase = result
        phase *= local_phase
    return det, phase


def _determinants(nspin: int, nelec: int):
    return [sum(1 << p for p in occupied)
            for occupied in itertools.combinations(range(nspin), nelec)]


def _spin_orbital_integrals(h_spatial, eri_chemist):
    """Return h[p,r] and (pq|rs) in a spin-orbital physicist convention."""

    nspatial = h_spatial.shape[0]
    nspin = 2 * nspatial
    h = np.zeros((nspin, nspin))
    g = np.zeros((nspin, nspin, nspin, nspin))
    for p, q, r, s in itertools.product(range(nspin), repeat=4):
        ps, qs, rs, ss = (p % nspatial, q % nspatial,
                          r % nspatial, s % nspatial)
        pspin, qspin = p // nspatial, q // nspatial
        rspin, sspin = r // nspatial, s // nspatial
        if pspin == rspin and qspin == sspin:
            g[p, q, r, s] = eri_chemist[ps, rs, qs, ss]
    for p, r in itertools.product(range(nspin), repeat=2):
        if p // nspatial == r // nspatial:
            h[p, r] = h_spatial[p % nspatial, r % nspatial]
    return h, g


def _dense_second_quantized_hamiltonian(determinants, h, g):
    """Construct H by direct fermionic operator application."""

    index = {det: position for position, det in enumerate(determinants)}
    nspin = h.shape[0]
    matrix = np.zeros((len(determinants), len(determinants)))
    for ket_index, ket in enumerate(determinants):
        for p, r in itertools.product(range(nspin), repeat=2):
            result = _apply(ket, (("ann", r), ("cre", p)))
            if result is not None:
                bra, phase = result
                matrix[index[bra], ket_index] += phase * h[p, r]
        for p, q, r, s in itertools.product(range(nspin), repeat=4):
            result = _apply(
                ket,
                (("ann", r), ("ann", s), ("cre", q), ("cre", p)),
            )
            if result is not None:
                bra, phase = result
                matrix[index[bra], ket_index] += 0.5 * phase * g[p, q, r, s]
    return matrix


def _slater_condon_element(bra, ket, h, g):
    """Independent Slater--Condon expression for <bra|H|ket>."""

    nspin = h.shape[0]
    occupied_bra = [p for p in range(nspin) if bra & (1 << p)]
    occupied_ket = [p for p in range(nspin) if ket & (1 << p)]
    removed = [p for p in occupied_ket if p not in occupied_bra]
    added = [p for p in occupied_bra if p not in occupied_ket]
    if len(removed) != len(added) or len(removed) > 2:
        return 0.0
    if not removed:
        value = sum(h[i, i] for i in occupied_ket)
        value += 0.5 * sum(
            g[i, j, i, j] - g[i, j, j, i]
            for i in occupied_ket for j in occupied_ket
        )
        return value
    if len(removed) == 1:
        i, a = removed[0], added[0]
        result = _apply(ket, (("ann", i), ("cre", a)))
        if result is None or result[0] != bra:
            raise AssertionError("invalid single-excitation phase")
        common = [p for p in occupied_ket if p != i]
        value = h[a, i] + sum(
            g[a, j, i, j] - g[a, j, j, i] for j in common
        )
        return result[1] * value
    i, j = sorted(removed)
    a, b = sorted(added)
    result = _apply(
        ket,
        (("ann", i), ("ann", j), ("cre", b), ("cre", a)),
    )
    if result is None or result[0] != bra:
        raise AssertionError("invalid double-excitation phase")
    return result[1] * (g[a, b, i, j] - g[a, b, j, i])


def _slater_condon_sigma(determinants, h, g, vector):
    sigma = np.zeros_like(vector)
    for i, bra in enumerate(determinants):
        sigma[i] = sum(
            _slater_condon_element(bra, ket, h, g) * vector[j]
            for j, ket in enumerate(determinants)
        )
    return sigma


def _symmetrize_chemist(raw):
    """Impose the eight real two-electron-integral permutations."""

    result = np.zeros_like(raw)
    for p, q, r, s in itertools.product(range(raw.shape[0]), repeat=4):
        result[p, q, r, s] = np.mean([
            raw[p, q, r, s], raw[q, p, r, s],
            raw[p, q, s, r], raw[q, p, s, r],
            raw[r, s, p, q], raw[s, r, p, q],
            raw[r, s, q, p], raw[s, r, q, p],
        ])
    return result


def _integral_coefficients(seed=71):
    rng = np.random.default_rng(seed)
    h_coeff = []
    eri_coeff = []
    for order, scale in enumerate((0.45, 0.17, 0.08)):
        h_raw = rng.normal(scale=scale, size=(4, 4))
        h_coeff.append(0.5 * (h_raw + h_raw.T))
        eri_raw = rng.normal(scale=0.4 * scale, size=(4, 4, 4, 4))
        eri_coeff.append(_symmetrize_chemist(eri_raw))
    # Keep the lowest state separated from the rest for state-specific checks.
    h_coeff[0] += np.diag((-1.8, -0.35, 0.22, 1.35))
    return h_coeff, eri_coeff


def _spin_squared(determinants, nspatial):
    index = {det: position for position, det in enumerate(determinants)}
    n = len(determinants)
    s_plus = np.zeros((n, n))
    s_minus = np.zeros((n, n))
    sz = np.zeros((n, n))
    for ket_index, ket in enumerate(determinants):
        nalpha = sum(bool(ket & (1 << p)) for p in range(nspatial))
        nbeta = sum(bool(ket & (1 << (nspatial + p)))
                    for p in range(nspatial))
        sz[ket_index, ket_index] = 0.5 * (nalpha - nbeta)
        for p in range(nspatial):
            result = _apply(
                ket,
                (("ann", nspatial + p), ("cre", p)),
            )
            if result is not None:
                s_plus[index[result[0]], ket_index] += result[1]
            result = _apply(
                ket,
                (("ann", p), ("cre", nspatial + p)),
            )
            if result is not None:
                s_minus[index[result[0]], ket_index] += result[1]
    return sz @ sz + 0.5 * (s_plus @ s_minus + s_minus @ s_plus)


def _mrsf_spaces():
    """Enumerate the two-parent model and return singlet/triplet CSFs."""

    nspatial = 4
    closed, open1, open2, virtual = range(nspatial)
    alpha = lambda p: p
    beta = lambda p: nspatial + p
    plus_parent = sum(1 << p for p in (
        alpha(closed), beta(closed), alpha(open1), alpha(open2)))
    minus_parent = sum(1 << p for p in (
        alpha(closed), beta(closed), beta(open1), beta(open2)))
    determinants = _determinants(2 * nspatial, 4)
    index = {det: position for position, det in enumerate(determinants)}

    def response(parent, hole, particle):
        result = _apply(parent, (("ann", hole), ("cre", particle)))
        if result is None:
            raise AssertionError("invalid parent response")
        vector = np.zeros(len(determinants))
        vector[index[result[0]]] = result[1]
        return vector, result[0]

    slots = [(p, q) for p in (closed, open1, open2)
             for q in (open1, open2, virtual)]
    plus = {}
    minus = {}
    generated = set()
    for p, q in slots:
        plus[p, q], det_plus = response(
            plus_parent, alpha(p), beta(q))
        minus[p, q], det_minus = response(
            minus_parent, beta(p), alpha(q))
        generated.update((det_plus, det_minus))

    invroot2 = 1.0 / math.sqrt(2.0)
    singlet_labels = ["G", "D", "OOS", "CO1", "CO2", "OV1", "OV2", "CV"]
    singlet_vectors = [
        plus[open2, open1],
        plus[open1, open2],
        invroot2 * (plus[open1, open1] - plus[open2, open2]),
        invroot2 * (plus[closed, open1] - minus[closed, open1]),
        invroot2 * (plus[closed, open2] - minus[closed, open2]),
        invroot2 * (plus[open1, virtual] - minus[open1, virtual]),
        invroot2 * (plus[open2, virtual] - minus[open2, virtual]),
        invroot2 * (plus[closed, virtual] - minus[closed, virtual]),
    ]
    triplet_labels = ["OOT", "CO1", "CO2", "OV1", "OV2", "CV"]
    triplet_vectors = [
        invroot2 * (plus[open1, open1] + plus[open2, open2]),
        invroot2 * (plus[closed, open1] + minus[closed, open1]),
        invroot2 * (plus[closed, open2] + minus[closed, open2]),
        invroot2 * (plus[open1, virtual] + minus[open1, virtual]),
        invroot2 * (plus[open2, virtual] + minus[open2, virtual]),
        invroot2 * (plus[closed, virtual] + minus[closed, virtual]),
    ]
    topology = (
        ("CO", "+1"), ("OV", "+1"), ("CV", "+1"),
        ("OO", "both"),
        ("CV", "-1"), ("OV", "-1"), ("CO", "-1"),
    )
    return {
        "determinants": determinants,
        "generated_determinants": generated,
        "topology": topology,
        "singlet": (singlet_labels, np.column_stack(singlet_vectors)),
        "triplet": (triplet_labels, np.column_stack(triplet_vectors)),
        "plus": plus,
        "minus": minus,
    }


def _spc_coefficients(labels):
    """Return C, C', and C'' with only CO--CO, CO--OV, and OV--OV terms."""

    matrices = [np.zeros((len(labels), len(labels))) for _ in range(3)]
    co = [i for i, label in enumerate(labels) if label.startswith("CO")]
    ov = [i for i, label in enumerate(labels) if label.startswith("OV")]
    scales = ((0.11, -0.07, 0.05), (0.035, 0.021, -0.017),
              (-0.012, 0.009, 0.006))
    for matrix, (coco, ovov, coov) in zip(matrices, scales):
        for local_i, i in enumerate(co):
            for local_j, j in enumerate(co):
                matrix[i, j] = coco * (1.0 if local_i == local_j else -0.4)
        for local_i, i in enumerate(ov):
            for local_j, j in enumerate(ov):
                matrix[i, j] = ovov * (1.0 if local_i == local_j else 0.3)
        for local_i, i in enumerate(co):
            for local_j, j in enumerate(ov):
                value = coov * (1.0 if local_i == local_j else -0.6)
                matrix[i, j] = value
                matrix[j, i] = value
    return matrices


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


class TestMRSFHessianReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = _mrsf_spaces()
        cls.h_coeff, cls.eri_coeff = _integral_coefficients()
        cls.hso = []
        cls.gso = []
        cls.full_hamiltonian = []
        for h, eri in zip(cls.h_coeff, cls.eri_coeff):
            hso, gso = _spin_orbital_integrals(h, eri)
            cls.hso.append(hso)
            cls.gso.append(gso)
            cls.full_hamiltonian.append(_dense_second_quantized_hamiltonian(
                cls.model["determinants"], hso, gso))

    def test_parent_topology_rank_and_phases(self):
        self.assertEqual(
            self.model["topology"],
            (("CO", "+1"), ("OV", "+1"), ("CV", "+1"),
             ("OO", "both"),
             ("CV", "-1"), ("OV", "-1"), ("CO", "-1")),
        )
        self.assertEqual(len(self.model["generated_determinants"]), 14)
        for spin, dimension in (("singlet", 8), ("triplet", 6)):
            labels, transform = self.model[spin]
            self.assertEqual(len(labels), dimension)
            self.assertEqual(np.linalg.matrix_rank(transform, tol=1.0e-12), dimension)
            np.testing.assert_allclose(
                transform.T @ transform, np.eye(dimension), atol=_TOL, rtol=0.0)

        # A missing second parent leaves the CV singlet with norm 1/sqrt(2).
        cv_mutation = self.model["plus"][0, 3] / math.sqrt(2.0)
        self.assertGreater(abs(np.linalg.norm(cv_mutation) - 1.0), 0.2)

    def test_spin_separation_and_documented_type_iv_boundary(self):
        s2 = _spin_squared(self.model["determinants"], 4)
        singlet_labels, singlet = self.model["singlet"]
        triplet_labels, triplet = self.model["triplet"]
        for column, label in enumerate(singlet_labels):
            residual = np.max(np.abs(s2 @ singlet[:, column]))
            if label == "CV":
                # Four of the six type-IV spin complements are absent in the
                # founding method, so this two-determinant CV singlet is not an
                # exact S^2 eigenfunction in the finite determinant model.
                self.assertGreater(residual, 0.5)
            else:
                self.assertLess(residual, _TOL, label)
        for column, label in enumerate(triplet_labels):
            residual = np.max(np.abs(
                s2 @ triplet[:, column] - 2.0 * triplet[:, column]))
            self.assertLess(residual, _TOL, label)

        # Changing the authoritative OOS L-R amplitude fold to L+R changes
        # the physical function from a singlet to a triplet.
        oos_mutation = (
            self.model["plus"][1, 1] + self.model["plus"][2, 2]
        ) / math.sqrt(2.0)
        np.testing.assert_allclose(s2 @ oos_mutation, 2.0 * oos_mutation,
                                   atol=_TOL, rtol=0.0)

    def test_columnwise_slater_condon_sigma(self):
        dense = self.full_hamiltonian[0]
        determinants = self.model["determinants"]
        for spin in ("singlet", "triplet"):
            _labels, transform = self.model[spin]
            projected = transform.T @ dense @ transform
            for column in range(transform.shape[1]):
                sigma_det = _slater_condon_sigma(
                    determinants, self.hso[0], self.gso[0], transform[:, column])
                sigma = transform.T @ sigma_det
                np.testing.assert_allclose(
                    sigma, projected[:, column], atol=1.0e-10, rtol=0.0)

    def test_singlet_and_triplet_eigenvalue_hessians(self):
        h0, h1, h2 = self.full_hamiltonian
        # Nakata et al., founding SI Eqs. S8.9--S8.10: the physical triplet
        # response contains A+C and the physical singlet response A-C.
        for spin, spc_sign in (("singlet", -1.0), ("triplet", 1.0)):
            labels, transform = self.model[spin]
            c0, c1, c2 = _spc_coefficients(labels)
            a0 = transform.T @ h0 @ transform + spc_sign * c0
            a1 = transform.T @ h1 @ transform + spc_sign * c1
            a2 = transform.T @ h2 @ transform + spc_sign * c2
            _value, vector, first, second, residual = _eigen_response(a0, a1, a2)
            self.assertLess(residual, 1.0e-10)
            first_errors = []
            second_errors = []
            for step in (8.0e-3, 4.0e-3, 2.0e-3):
                plus = a0 + step * a1 + 0.5 * step * step * a2
                minus = a0 - step * a1 + 0.5 * step * step * a2
                eplus, overlap_plus = _tracked_value(plus, vector)
                eminus, overlap_minus = _tracked_value(minus, vector)
                self.assertGreater(min(overlap_plus, overlap_minus), 0.99)
                e0 = float(vector @ a0 @ vector)
                first_fd = (eplus - eminus) / (2.0 * step)
                second_fd = (eplus - 2.0 * e0 + eminus) / (step * step)
                first_errors.append(abs(first_fd - first))
                second_errors.append(abs(second_fd - second))
            self.assertLess(first_errors[-1], 2.0e-7, (spin, first_errors))
            self.assertLess(second_errors[-1], 2.0e-6, (spin, second_errors))
            self.assertLess(first_errors[-1], first_errors[0] / 3.0)
            self.assertLess(second_errors[-1], second_errors[0] / 3.0)

    def test_spc_has_only_the_three_authoritative_contractions(self):
        for spin in ("singlet", "triplet"):
            labels, _transform = self.model[spin]
            for matrix in _spc_coefficients(labels):
                for i, left in enumerate(labels):
                    for j, right in enumerate(labels):
                        allowed = (
                            (left.startswith("CO") and right.startswith("CO"))
                            or (left.startswith("OV") and right.startswith("OV"))
                            or (left.startswith("CO") and right.startswith("OV"))
                            or (left.startswith("OV") and right.startswith("CO"))
                        )
                        if not allowed:
                            self.assertEqual(matrix[i, j], 0.0)


if __name__ == "__main__":
    unittest.main()
