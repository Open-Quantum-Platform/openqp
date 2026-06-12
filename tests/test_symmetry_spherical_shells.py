import importlib.util
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import unittest


ROOT = Path(__file__).resolve().parents[1]
SUPPORTED_PURE_L = (2, 3, 4, 5, 6)


def load_module(name, relpath):
    spec = importlib.util.spec_from_file_location(name, str(ROOT / relpath))
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def random_rotation(seed):
    rng = np.random.RandomState(seed)
    q, _ = np.linalg.qr(rng.randn(3, 3))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


WATER_CHARGES = [8, 1, 1]
WATER_COORDS = np.array(
    [
        [0.0, 0.0, 0.1173],
        [0.0, 0.7572, -0.4692],
        [0.0, -0.7572, -0.4692],
    ]
)
# O: s, p, pure-d ; H: s each -> 1+3+5+1+1 = 11 AOs.
WATER_SPH_SHELLS = [(0, 0, False), (0, 1, False), (0, 2, True), (1, 0), (2, 0)]
WATER_SPH_NAO = 11


class TestSolidHarmonics(unittest.TestCase):
    def setUp(self):
        self.symmetry = load_module(
            "openqp_symmetry_spherical_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )

    def test_metric_orthonormality(self):
        for l in SUPPORTED_PURE_L:
            b, metric = self.symmetry._spherical_basis(l)
            deviation = np.max(np.abs(b @ metric @ b.T - np.eye(2 * l + 1)))
            self.assertLess(float(deviation), 1.0e-12, f"l={l}")

    def test_high_l_shell_sizes_match_openqp_h_i_limits(self):
        expected = {
            5: (21, 11),
            6: (28, 13),
        }
        for l, (ncart, nsph) in expected.items():
            with self.subTest(l=l):
                self.assertEqual(self.symmetry._cartesian_shell_size(l), ncart)
                self.assertEqual(self.symmetry._shell_size(l, pure=True), nsph)
                self.assertEqual(len(self.symmetry._CART_MONOMIALS[l]), ncart)

    def test_d_shell_matches_textbook_forms(self):
        b, _ = self.symmetry._spherical_basis(2)
        monomials = self.symmetry._CART_MONOMIALS[2]
        norms = self.symmetry._monomial_norms(2)
        index = {m: i for i, m in enumerate(monomials)}

        def mono_coeffs(row):
            return b[row] * norms

        # m = -2 ~ xy ; m = -1 ~ yz ; m = +1 ~ xz : single monomials.
        for row, mono in [(0, (1, 1, 0)), (1, (0, 1, 1)), (3, (1, 0, 1))]:
            coeffs = mono_coeffs(row)
            self.assertAlmostEqual(abs(coeffs[index[mono]]), 1.0, places=12)
            self.assertAlmostEqual(
                float(np.sum(np.abs(coeffs))), abs(coeffs[index[mono]]), places=12)
        # m = 0 ~ 2z^2 - x^2 - y^2 (ratios -1 : -1 : 2).
        c0 = mono_coeffs(2)
        self.assertAlmostEqual(
            c0[index[(0, 0, 2)]] / c0[index[(2, 0, 0)]], -2.0, places=12)
        self.assertAlmostEqual(
            c0[index[(2, 0, 0)]], c0[index[(0, 2, 0)]], places=12)
        # m = +2 ~ x^2 - y^2.
        c2 = mono_coeffs(4)
        self.assertAlmostEqual(
            c2[index[(2, 0, 0)]], -c2[index[(0, 2, 0)]], places=12)

    def test_rotation_blocks_orthogonal_and_representation(self):
        m1 = random_rotation(1)
        m2 = random_rotation(2)
        for l in (1, *SUPPORTED_PURE_L):
            b1 = self.symmetry._shell_block_spherical(l, m1)
            b2 = self.symmetry._shell_block_spherical(l, m2)
            b12 = self.symmetry._shell_block_spherical(l, m1 @ m2)
            n = 2 * l + 1
            self.assertLess(float(np.max(np.abs(b1 @ b1.T - np.eye(n)))), 1.0e-12)
            self.assertLess(float(np.max(np.abs(b12 - b1 @ b2))), 1.0e-10)

    def test_sign_operations_are_diagonal(self):
        for l in (1, *SUPPORTED_PURE_L):
            for signs in [(-1, 1, 1), (1, -1, 1), (1, 1, -1), (-1, -1, -1)]:
                comp = self.symmetry._component_signs_any(l, True, signs)
                self.assertEqual(len(comp), 2 * l + 1)
                self.assertTrue(all(c in (-1, 1) for c in comp))


class TestSphericalSALCsAndLabels(unittest.TestCase):
    def setUp(self):
        self.symmetry = load_module(
            "openqp_symmetry_sph_salc_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )
        self.detect = load_module(
            "openqp_symmetry_detect_sph_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )
        self.detection = self.detect.detect_point_group(
            WATER_CHARGES, WATER_COORDS, tolerance=1.0e-6
        )

    def test_salc_transform_complete_and_orthogonal(self):
        u, labels = self.symmetry.build_symmetry_adapted_transform(
            WATER_SPH_SHELLS, self.detection["operations"],
            self.detection["character_table"],
        )
        self.assertEqual(u.shape, (WATER_SPH_NAO, WATER_SPH_NAO))
        deviation = np.max(np.abs(u.T @ u - np.eye(WATER_SPH_NAO)))
        self.assertLess(float(deviation), 1.0e-12)
        counts = Counter(labels)
        # Pure d in c2v: z2 and x2-y2 -> a1, xy -> a2, xz/yz -> b1/b2.
        # a1 total: O s + O pz + H+ + two d = 5. b's: px/py (2) + H- (1)
        # + xz/yz (2) split 2/3 depending on the in-plane convention.
        self.assertEqual(counts["a1"], 5)
        self.assertEqual(counts["a2"], 1)
        self.assertEqual(sorted([counts["b1"], counts["b2"]]), [2, 3])

    def test_mo_labels_invariant_under_input_rotation(self):
        rotation = random_rotation(42)
        rotated = WATER_COORDS @ rotation.T
        det_rot = self.detect.detect_point_group(
            WATER_CHARGES, rotated, tolerance=1.0e-6
        )

        u_ref, labels_ref = self.symmetry.build_symmetry_adapted_transform(
            WATER_SPH_SHELLS, self.detection["operations"],
            self.detection["character_table"],
        )
        # Express the standard-frame SALCs in the rotated input frame.
        r_ref = np.asarray(self.detection["orientation"])
        w3 = rotation @ r_ref.T
        w = np.eye(WATER_SPH_NAO)
        w[1:4, 1:4] = self.symmetry._shell_block_any(1, False, w3)
        w[4:9, 4:9] = self.symmetry._shell_block_any(2, True, w3)
        mos = w @ u_ref

        result = self.symmetry.assign_mo_irreps(
            mos, np.eye(WATER_SPH_NAO), WATER_SPH_SHELLS,
            det_rot["operations"], det_rot["character_table"],
            matrix_key="matrix_input_frame",
        )
        self.assertLess(result["max_deviation"], 1.0e-8)
        self.assertEqual(Counter(result["labels"]), Counter(labels_ref))

    def test_reduction_maps_with_pure_shells(self):
        maps = self.symmetry.build_reduction_maps(
            WATER_SPH_SHELLS, self.detection["operations"]
        )
        self.assertEqual(maps["n_ao"], WATER_SPH_NAO)
        for target in maps["ao_target"]:
            self.assertEqual(sorted(target), list(range(WATER_SPH_NAO)))
        for sign_row in maps["ao_sign"]:
            self.assertTrue(all(s in (-1, 1) for s in sign_row))


if __name__ == "__main__":
    unittest.main()
