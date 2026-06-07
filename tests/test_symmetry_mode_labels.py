import importlib.util
import sys
import types
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


WATER_CHARGES = [8, 1, 1]
WATER_COORDS = np.array(
    [
        [0.0, 0.0, 0.1173],
        [0.0, 0.7572, -0.4692],
        [0.0, -0.7572, -0.4692],
    ]
)


def water_modes():
    """Symmetry-exact model modes: two a1-type and one in-plane b-type.

    Water sits in the x=0 plane with the C2 axis along z, so a1 modes are
    invariant and the asymmetric stretch flips under C2z.
    """
    sym_stretch = np.array(
        [
            [0.0, 0.0, 0.1],
            [0.0, 0.4, -0.3],
            [0.0, -0.4, -0.3],
        ]
    ).ravel()
    bend = np.array(
        [
            [0.0, 0.0, -0.2],
            [0.0, 0.3, 0.4],
            [0.0, -0.3, 0.4],
        ]
    ).ravel()
    asym_stretch = np.array(
        [
            [0.0, 0.2, 0.0],
            [0.0, 0.4, 0.3],
            [0.0, 0.4, -0.3],
        ]
    ).ravel()
    return np.vstack([sym_stretch, bend, asym_stretch])


class TestModeIrrepAssignment(unittest.TestCase):
    def setUp(self):
        self.symmetry = load_module(
            "openqp_symmetry_modes_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )
        self.detect = load_module(
            "openqp_symmetry_detect_modes_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )
        self.detection = self.detect.detect_point_group(
            WATER_CHARGES, WATER_COORDS, tolerance=1.0e-6
        )

    def test_water_modes_get_expected_labels(self):
        result = self.symmetry.assign_mode_irreps(
            water_modes(),
            self.detection["operations"],
            self.detection["character_table"],
            matrix_key="matrix_input_frame",
        )
        self.assertEqual(result["labels"][0], "a1")
        self.assertEqual(result["labels"][1], "a1")
        self.assertIn(result["labels"][2], ("b1", "b2"))
        self.assertLess(result["max_deviation"], 1.0e-10)

    def test_labels_invariant_under_input_rotation(self):
        rng = np.random.RandomState(11)
        q, _ = np.linalg.qr(rng.randn(3, 3))
        if np.linalg.det(q) < 0:
            q[:, 0] = -q[:, 0]

        rotated_coords = WATER_COORDS @ q.T
        rotated_detection = self.detect.detect_point_group(
            WATER_CHARGES, rotated_coords, tolerance=1.0e-6
        )
        modes = water_modes().reshape(3, 3, 3)
        rotated_modes = np.einsum("ij,maj->mai", q, modes).reshape(3, 9)

        result = self.symmetry.assign_mode_irreps(
            rotated_modes,
            rotated_detection["operations"],
            rotated_detection["character_table"],
            matrix_key="matrix_input_frame",
        )
        self.assertEqual(result["labels"][:2], ["a1", "a1"])
        self.assertIn(result["labels"][2], ("b1", "b2"))
        self.assertLess(result["max_deviation"], 1.0e-8)

    def test_symmetry_broken_mode_is_mixed(self):
        modes = water_modes().copy()
        modes[0] = modes[0] + 0.5 * modes[2]  # contaminate a1 with b-type
        result = self.symmetry.assign_mode_irreps(
            modes,
            self.detection["operations"],
            self.detection["character_table"],
            matrix_key="matrix_input_frame",
        )
        self.assertEqual(result["labels"][0], "mixed")

    def test_molecule_wiring_labels_modes(self):
        detect = load_module(
            "oqp_symmetry_detect_modewiring_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )
        symmetry = load_module(
            "oqp_symmetry_modewiring_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )
        metadata_tests = load_module(
            "openqp_symmetry_metadata_tests_for_modes",
            "tests/test_symmetry_metadata.py",
        )
        molecule_mod = metadata_tests.load_molecule_module()

        library_pkg = types.ModuleType("oqp.library")
        library_pkg.__path__ = [str(ROOT / "pyoqp/oqp/library")]
        library_pkg.symmetry_detect = detect
        library_pkg.symmetry = symmetry
        sys.modules["oqp.library"] = library_pkg
        sys.modules["oqp.library.symmetry_detect"] = detect
        sys.modules["oqp.library.symmetry"] = symmetry

        mol = object.__new__(molecule_mod.Molecule)
        mol.config = {"symmetry": {"enabled": "true"}}
        mol.symmetry_metadata = {}
        mol.get_atoms = lambda: np.array(WATER_CHARGES, dtype=float)
        mol.get_system = lambda: WATER_COORDS.ravel().copy()
        mol.initialize_symmetry_metadata()
        mol.modes = water_modes()

        result = mol.label_normal_modes()
        self.assertIsNotNone(result)
        self.assertEqual(result["status"], "ok")
        self.assertEqual(result["labels"][:2], ["a1", "a1"])
        self.assertEqual(
            mol.symmetry_metadata["mode_labels"]["labels"], result["labels"]
        )

    def test_label_modes_false_skips(self):
        detect = load_module(
            "oqp_symmetry_detect_modeskip_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )
        metadata_tests = load_module(
            "openqp_symmetry_metadata_tests_for_modeskip",
            "tests/test_symmetry_metadata.py",
        )
        molecule_mod = metadata_tests.load_molecule_module()
        sys.modules["oqp.library.symmetry_detect"] = detect

        mol = object.__new__(molecule_mod.Molecule)
        mol.config = {"symmetry": {"enabled": "true", "label_modes": "false"}}
        mol.symmetry_metadata = {}
        mol.get_atoms = lambda: np.array(WATER_CHARGES, dtype=float)
        mol.get_system = lambda: WATER_COORDS.ravel().copy()
        mol.initialize_symmetry_metadata()
        mol.modes = water_modes()

        self.assertIsNone(mol.label_normal_modes())
        self.assertNotIn("mode_labels", mol.symmetry_metadata)


if __name__ == "__main__":
    unittest.main()
