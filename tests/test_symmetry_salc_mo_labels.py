import importlib.util
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relpath):
    spec = importlib.util.spec_from_file_location(name, str(ROOT / relpath))
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def water_detection(detect):
    charges = [8, 1, 1]
    coords = [
        [0.0, 0.0, 0.1173],
        [0.0, 0.7572, -0.4692],
        [0.0, -0.7572, -0.4692],
    ]
    return detect.detect_point_group(charges, coords, tolerance=1.0e-6)


# Minimal water basis: O 1s, 2s, 2p; H 1s each.
WATER_SHELLS = [(0, 0), (0, 0), (0, 1), (1, 0), (2, 0)]
WATER_NAO = 7


class TestSymmetryAdaptedTransform(unittest.TestCase):
    def setUp(self):
        self.symmetry = load_module(
            "openqp_symmetry_salc_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )
        self.detect = load_module(
            "openqp_symmetry_detect_salc_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )
        self.detection = water_detection(self.detect)

    def build_water_transform(self):
        return self.symmetry.build_symmetry_adapted_transform(
            WATER_SHELLS,
            self.detection["operations"],
            self.detection["character_table"],
        )

    def test_water_transform_is_complete_and_orthogonal(self):
        u, labels = self.build_water_transform()
        self.assertEqual(u.shape, (WATER_NAO, WATER_NAO))
        self.assertEqual(len(labels), WATER_NAO)
        deviation = np.max(np.abs(u.T @ u - np.eye(WATER_NAO)))
        self.assertLess(float(deviation), 1.0e-12)

    def test_water_salc_irrep_counts(self):
        _, labels = self.build_water_transform()
        counts = Counter(labels)
        # O 1s, 2s, p_z and the symmetric H combination are a1.
        self.assertEqual(counts["a1"], 4)
        self.assertEqual(counts.get("a2", 0), 0)
        # In-plane p + antisymmetric H combination (2) and out-of-plane p (1)
        # split between b1/b2 depending on the standard orientation.
        b_counts = sorted([counts.get("b1", 0), counts.get("b2", 0)])
        self.assertEqual(b_counts, [1, 2])

    def test_invariant_matrix_is_block_diagonal_in_salc_basis(self):
        u, labels = self.build_water_transform()
        operations = self.detection["operations"]

        rng = np.random.RandomState(7)
        seed_matrix = rng.rand(WATER_NAO, WATER_NAO)
        seed_matrix = 0.5 * (seed_matrix + seed_matrix.T)

        # Group-average the seed so it commutes with every operation.
        n_ao, maps = self.symmetry._ao_operator_maps(
            self.symmetry._normalize_shells(WATER_SHELLS), operations
        )
        self.assertEqual(n_ao, WATER_NAO)
        averaged = np.zeros_like(seed_matrix)
        for target, sign in maps:
            t = np.zeros((WATER_NAO, WATER_NAO))
            t[target, np.arange(WATER_NAO)] = sign
            averaged += t @ seed_matrix @ t.T
        averaged /= len(operations)

        payload = self.symmetry.one_electron_block_diagnostics(
            averaged, u, labels, tolerance=1.0e-10
        )
        self.assertTrue(payload["within_tolerance"])
        self.assertTrue(payload["transform_orthogonality"]["ok"])

    def test_salc_columns_label_consistently_as_mos(self):
        u, labels = self.build_water_transform()
        result = self.symmetry.assign_mo_irreps(
            u,
            np.eye(WATER_NAO),
            WATER_SHELLS,
            self.detection["operations"],
            self.detection["character_table"],
        )
        self.assertEqual(result["labels"], labels)
        self.assertLess(result["max_deviation"], 1.0e-10)
        self.assertEqual(result["status"], "label_only_no_reductions")

    def test_mixed_orbital_is_flagged(self):
        u, labels = self.build_water_transform()
        ia1 = labels.index("a1")
        ib = next(i for i, lbl in enumerate(labels) if lbl.startswith("b"))
        mixed = (u[:, ia1] + u[:, ib]) / np.sqrt(2.0)
        result = self.symmetry.assign_mo_irreps(
            mixed[:, None],
            np.eye(WATER_NAO),
            WATER_SHELLS,
            self.detection["operations"],
            self.detection["character_table"],
        )
        self.assertEqual(result["labels"], ["mixed"])

    def test_equivalent_atoms_must_share_shell_lists(self):
        bad_shells = [(0, 0), (1, 0), (1, 0), (2, 0)]  # H1 has 2 s, H2 has 1
        with self.assertRaises(ValueError):
            self.symmetry.build_symmetry_adapted_transform(
                bad_shells,
                self.detection["operations"],
                self.detection["character_table"],
            )

    def test_non_sign_operation_matrices_are_rejected(self):
        operations = [
            {
                "name": "C4z",
                "matrix": [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                "permutation": [0, 1, 2],
            }
        ]
        with self.assertRaises(ValueError):
            self.symmetry.build_symmetry_adapted_transform(
                WATER_SHELLS, operations, {"a": [1]}
            )

    def test_d_shell_components_transform_with_product_signs(self):
        # Single atom at the origin with one Cartesian d shell under C2z.
        detection = self.detect.detect_point_group([6], [[0.0, 0.0, 0.0]])
        u, labels = self.symmetry.build_symmetry_adapted_transform(
            [(0, 2)],
            detection["operations"],
            detection["character_table"],
        )
        self.assertEqual(u.shape, (6, 6))
        counts = Counter(labels)
        # xx, yy, zz -> ag; xy -> b1g; xz -> b2g; yz -> b3g in D2h.
        self.assertEqual(counts["ag"], 3)
        self.assertEqual(counts["b1g"], 1)
        self.assertEqual(counts["b2g"], 1)
        self.assertEqual(counts["b3g"], 1)


class TestMoleculeDetectionWiring(unittest.TestCase):
    def test_initialize_symmetry_metadata_detects_water(self):
        detect = load_module(
            "oqp_symmetry_detect_wiring_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )

        metadata_tests = load_module(
            "openqp_symmetry_metadata_tests_for_wiring",
            "tests/test_symmetry_metadata.py",
        )
        molecule_mod = metadata_tests.load_molecule_module()

        # Register the detection module where molecule.py imports it from.
        import types

        library_pkg = types.ModuleType("oqp.library")
        library_pkg.__path__ = [str(ROOT / "pyoqp/oqp/library")]
        library_pkg.symmetry_detect = detect
        sys.modules["oqp.library"] = library_pkg
        sys.modules["oqp.library.symmetry_detect"] = detect

        mol = object.__new__(molecule_mod.Molecule)
        mol.config = {"symmetry": {"enabled": "true"}}
        mol.symmetry_metadata = {}
        mol.get_atoms = lambda: np.array([8.0, 1.0, 1.0])
        mol.get_system = lambda: np.array(
            [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692]
        )

        metadata = mol.initialize_symmetry_metadata()
        self.assertEqual(metadata["status"], "enabled")
        self.assertEqual(metadata["detected_point_group"], "c2v")
        self.assertEqual(metadata["detected_subgroup"], "c2v")
        self.assertEqual(metadata["point_group"], "c2v")
        self.assertNotIn("detection_error", metadata)
        self.assertFalse(metadata["use_integral_symmetry"])
        self.assertFalse(metadata["use_response_symmetry"])

    def test_strict_mismatch_raises(self):
        detect = load_module(
            "oqp_symmetry_detect_strict_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )
        metadata_tests = load_module(
            "openqp_symmetry_metadata_tests_for_strict",
            "tests/test_symmetry_metadata.py",
        )
        molecule_mod = metadata_tests.load_molecule_module()

        import types

        library_pkg = types.ModuleType("oqp.library")
        library_pkg.__path__ = [str(ROOT / "pyoqp/oqp/library")]
        library_pkg.symmetry_detect = detect
        sys.modules["oqp.library"] = library_pkg
        sys.modules["oqp.library.symmetry_detect"] = detect

        mol = object.__new__(molecule_mod.Molecule)
        mol.config = {
            "symmetry": {
                "enabled": "true",
                "point_group": "d6h",
                "strict": "true",
            }
        }
        mol.symmetry_metadata = {}
        mol.get_atoms = lambda: np.array([8.0, 1.0, 1.0])
        mol.get_system = lambda: np.array(
            [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692]
        )

        with self.assertRaises(ValueError):
            mol.initialize_symmetry_metadata()


if __name__ == "__main__":
    unittest.main()
