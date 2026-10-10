import importlib.util
import sys
from pathlib import Path

import numpy as np
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load_detect_module():
    spec = importlib.util.spec_from_file_location(
        "openqp_symmetry_detect_under_test",
        str(ROOT / "pyoqp/oqp/library/symmetry_detect.py"),
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def water():
    return [8, 1, 1], [
        [0.0, 0.0, 0.1173],
        [0.0, 0.7572, -0.4692],
        [0.0, -0.7572, -0.4692],
    ]


def ammonia():
    coords = [[0.0, 0.0, 0.38]]
    for angle_deg in (90.0, 210.0, 330.0):
        angle = np.deg2rad(angle_deg)
        coords.append([0.95 * np.cos(angle), 0.95 * np.sin(angle), -0.13])
    return [7, 1, 1, 1], coords


def methane():
    s = 1.0 / np.sqrt(3.0)
    return [6, 1, 1, 1, 1], [
        [0.0, 0.0, 0.0],
        [s, s, s],
        [s, -s, -s],
        [-s, s, -s],
        [-s, -s, s],
    ]


def benzene():
    charges = []
    coords = []
    for k in range(6):
        angle = np.deg2rad(60.0 * k)
        charges.append(6)
        coords.append([1.39 * np.cos(angle), 1.39 * np.sin(angle), 0.0])
        charges.append(1)
        coords.append([2.49 * np.cos(angle), 2.49 * np.sin(angle), 0.0])
    return charges, coords


class TestPointGroupDetection(unittest.TestCase):
    def setUp(self):
        self.detect = load_detect_module()

    def check(self, charges, coords, point_group, subgroup):
        result = self.detect.detect_point_group(charges, coords, tolerance=1.0e-6)
        self.assertEqual(result["point_group"], point_group)
        self.assertEqual(result["abelian_subgroup"], subgroup)
        return result

    def test_single_atom_is_kh_with_d2h_subgroup(self):
        self.check([8], [[1.0, 2.0, 3.0]], "kh", "d2h")

    def test_water_is_c2v(self):
        charges, coords = water()
        result = self.check(charges, coords, "c2v", "c2v")
        self.assertEqual(sorted(result["irreps"]), ["a1", "a2", "b1", "b2"])

    def test_ammonia_is_c3v_with_cs_subgroup(self):
        charges, coords = ammonia()
        self.check(charges, coords, "c3v", "cs")

    def test_methane_is_td_with_d2_subgroup(self):
        charges, coords = methane()
        self.check(charges, coords, "td", "d2")

    def test_benzene_is_d6h_with_d2h_subgroup(self):
        charges, coords = benzene()
        self.check(charges, coords, "d6h", "d2h")

    def test_co2_is_dooh_with_d2h_subgroup(self):
        self.check(
            [8, 6, 8],
            [[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.16]],
            "dooh",
            "d2h",
        )

    def test_hcn_is_coov_with_c2v_subgroup(self):
        self.check(
            [1, 6, 7],
            [[0.0, 0.0, -1.064], [0.0, 0.0, 0.0], [0.0, 0.0, 1.156]],
            "coov",
            "c2v",
        )

    def test_hydrogen_peroxide_skew_is_c2(self):
        self.check(
            [8, 8, 1, 1],
            [
                [0.7, 0.0, 0.0],
                [-0.7, 0.0, 0.0],
                [1.2, 0.8, 0.5],
                [-1.2, -0.8, 0.5],
            ],
            "c2",
            "c2",
        )

    def test_trans_difluoroethylene_is_c2h(self):
        self.check(
            [6, 6, 9, 9, 1, 1],
            [
                [0.66, 0.0, 0.0],
                [-0.66, 0.0, 0.0],
                [1.4, 1.1, 0.0],
                [-1.4, -1.1, 0.0],
                [1.2, -1.0, 0.0],
                [-1.2, 1.0, 0.0],
            ],
            "c2h",
            "c2h",
        )

    def test_ethylene_is_d2h(self):
        self.check(
            [6, 6, 1, 1, 1, 1],
            [
                [0.665, 0.0, 0.0],
                [-0.665, 0.0, 0.0],
                [1.23, 0.92, 0.0],
                [1.23, -0.92, 0.0],
                [-1.23, 0.92, 0.0],
                [-1.23, -0.92, 0.0],
            ],
            "d2h",
            "d2h",
        )

    def test_allene_is_d2d_with_d2_subgroup(self):
        self.check(
            [6, 6, 6, 1, 1, 1, 1],
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.31],
                [0.0, 0.0, -1.31],
                [0.93, 0.0, 1.87],
                [-0.93, 0.0, 1.87],
                [0.0, 0.93, -1.87],
                [0.0, -0.93, -1.87],
            ],
            "d2d",
            "d2",
        )

    def test_inversion_only_molecule_is_ci(self):
        charges = [6, 6, 9, 9, 17, 17, 1, 1]
        plus = [
            [0.7, 0.1, 0.2],
            [1.3, 1.0, 0.4],
            [1.5, -0.9, -0.8],
            [0.6, 0.5, -1.1],
        ]
        coords = []
        for vec in plus:
            coords.append(vec)
            coords.append([-x for x in vec])
        self.check(charges, coords, "ci", "ci")

    def test_generic_molecule_is_c1(self):
        self.check(
            [6, 1, 9, 17, 35],
            [
                [0.0, 0.0, 0.0],
                [0.6, 0.6, 0.6],
                [-1.0, 0.2, 0.1],
                [0.1, -1.2, 0.3],
                [0.2, 0.3, 1.5],
            ],
            "c1",
            "c1",
        )

    def test_operations_carry_valid_permutations(self):
        charges, coords = water()
        result = self.detect.detect_point_group(charges, coords, tolerance=1.0e-6)
        charges_arr = np.asarray(charges, dtype=float)
        centered = np.asarray(coords, dtype=float)
        centered -= np.einsum(
            "i,ij->j", charges_arr, centered
        ) / float(np.sum(charges_arr))
        rotation = np.asarray(result["orientation"])
        rotated = centered @ rotation.T
        for op in result["operations"]:
            matrix = np.asarray(op["matrix"])
            perm = op["permutation"]
            self.assertEqual(sorted(perm), [0, 1, 2])
            transformed = rotated @ matrix.T
            for i, j in enumerate(perm):
                self.assertLess(
                    float(np.linalg.norm(transformed[i] - rotated[j])), 1.0e-6
                )
                self.assertEqual(charges[i], charges[j])

    def test_detection_payload_is_json_serializable(self):
        import json

        charges, coords = water()
        result = self.detect.detect_point_group(charges, coords, tolerance=1.0e-6)
        json.dumps(result)

    def test_attach_detection_resolves_auto_requests(self):
        charges, coords = water()
        metadata = {
            "requested_point_group": "auto",
            "requested_subgroup": "auto",
            "point_group": "c1",
            "subgroup": "c1",
            "tolerance": 1.0e-6,
        }
        self.detect.attach_detection_metadata(metadata, charges, coords)
        self.assertEqual(metadata["point_group"], "c2v")
        self.assertEqual(metadata["subgroup"], "c2v")
        self.assertEqual(metadata["detected_point_group"], "c2v")
        self.assertTrue(metadata["requested_matches_detected"])
        self.assertFalse(metadata["use_integral_symmetry"])
        self.assertFalse(metadata["use_response_symmetry"])

    def test_attach_detection_reports_request_mismatch(self):
        charges, coords = water()
        metadata = {
            "requested_point_group": "d2h",
            "requested_subgroup": "auto",
            "point_group": "d2h",
            "subgroup": "c1",
            "tolerance": 1.0e-6,
        }
        self.detect.attach_detection_metadata(metadata, charges, coords)
        # Explicit request is honored, mismatch recorded.
        self.assertEqual(metadata["point_group"], "d2h")
        self.assertEqual(metadata["detected_point_group"], "c2v")
        self.assertFalse(metadata["requested_matches_detected"])

    def test_sf6_is_oh_with_d2h_subgroup(self):
        self.check(
            [16, 9, 9, 9, 9, 9, 9],
            [
                [0.0, 0.0, 0.0],
                [1.56, 0.0, 0.0],
                [-1.56, 0.0, 0.0],
                [0.0, 1.56, 0.0],
                [0.0, -1.56, 0.0],
                [0.0, 0.0, 1.56],
                [0.0, 0.0, -1.56],
            ],
            "oh",
            "d2h",
        )

    def test_staggered_ethane_is_d3d_with_c2h_subgroup(self):
        coords = [[0.0, 0.0, 0.77], [0.0, 0.0, -0.77]]
        for k in range(3):
            angle = np.deg2rad(120.0 * k)
            coords.append([1.02 * np.cos(angle), 1.02 * np.sin(angle), 1.16])
        for k in range(3):
            angle = np.deg2rad(120.0 * k + 60.0)
            coords.append([1.02 * np.cos(angle), 1.02 * np.sin(angle), -1.16])
        self.check([6, 6, 1, 1, 1, 1, 1, 1], coords, "d3d", "c2h")

    def test_bf3_is_d3h_with_c2v_subgroup(self):
        coords = [[0.0, 0.0, 0.0]]
        for k in range(3):
            angle = np.deg2rad(120.0 * k)
            coords.append([1.3 * np.cos(angle), 1.3 * np.sin(angle), 0.0])
        self.check([5, 9, 9, 9], coords, "d3h", "c2v")

    def test_xef4_is_d4h_with_d2h_subgroup(self):
        self.check(
            [54, 9, 9, 9, 9],
            [
                [0.0, 0.0, 0.0],
                [1.95, 0.0, 0.0],
                [-1.95, 0.0, 0.0],
                [0.0, 1.95, 0.0],
                [0.0, -1.95, 0.0],
            ],
            "d4h",
            "d2h",
        )

    def test_noisy_geometry_detected_within_loose_tolerance(self):
        rng = np.random.RandomState(0)
        charges, coords = water()
        noisy = np.asarray(coords) + rng.randn(3, 3) * 1.0e-5
        result = self.detect.detect_point_group(charges, noisy, tolerance=1.0e-3)
        self.assertEqual(result["point_group"], "c2v")
        self.assertEqual(result["abelian_subgroup"], "c2v")

    def test_tolerance_must_be_positive(self):
        charges, coords = water()
        with self.assertRaises(ValueError):
            self.detect.detect_point_group(charges, coords, tolerance=0.0)


if __name__ == "__main__":
    unittest.main()
