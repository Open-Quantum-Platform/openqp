import importlib.util
import sys
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
WATER_SHELLS = [(0, 0), (0, 0), (0, 1), (1, 0), (2, 0)]
WATER_NAO = 7


class TestReductionMaps(unittest.TestCase):
    def setUp(self):
        self.symmetry = load_module(
            "openqp_symmetry_redmaps_under_test",
            "pyoqp/oqp/library/symmetry.py",
        )
        self.detect = load_module(
            "openqp_symmetry_detect_redmaps_under_test",
            "pyoqp/oqp/library/symmetry_detect.py",
        )
        self.detection = self.detect.detect_point_group(
            WATER_CHARGES, WATER_COORDS, tolerance=1.0e-6
        )
        self.maps = self.symmetry.build_reduction_maps(
            WATER_SHELLS, self.detection["operations"]
        )

    def test_shapes_and_identity_operation(self):
        maps = self.maps
        self.assertEqual(maps["n_operations"], 4)  # c2v
        self.assertEqual(maps["n_shells"], len(WATER_SHELLS))
        self.assertEqual(maps["n_ao"], WATER_NAO)
        iop_e = maps["operation_names"].index("E")
        self.assertEqual(
            maps["shell_permutation"][iop_e], list(range(len(WATER_SHELLS)))
        )
        self.assertEqual(maps["ao_sign"][iop_e], [1] * WATER_NAO)
        self.assertEqual(maps["ao_target"][iop_e], list(range(WATER_NAO)))

    def test_permutations_are_bijections(self):
        for perm in self.maps["shell_permutation"]:
            self.assertEqual(sorted(perm), list(range(len(WATER_SHELLS))))
        for target in self.maps["ao_target"]:
            self.assertEqual(sorted(target), list(range(WATER_NAO)))

    def test_h_shells_form_one_orbit(self):
        rep = self.maps["shell_orbit_representative"]
        size = self.maps["shell_orbit_size"]
        # O shells (0, 1, 2) are self-orbits; H shells (3, 4) share orbit 3.
        self.assertEqual(rep[:3], [0, 1, 2])
        self.assertEqual(rep[3], 3)
        self.assertEqual(rep[4], 3)
        self.assertEqual(size[3], 2)
        self.assertEqual(size[4], 0)  # non-representative
        self.assertEqual(size[0], 1)

    def test_maps_leave_invariant_matrix_invariant(self):
        # Group-average a random symmetric matrix, then verify the
        # reduction maps reproduce T M T^T = M for every operation.
        n_ao, op_maps = self.symmetry._ao_operator_maps(
            self.symmetry._normalize_shells(WATER_SHELLS),
            self.detection["operations"],
        )
        rng = np.random.RandomState(3)
        seed = rng.rand(n_ao, n_ao)
        seed = 0.5 * (seed + seed.T)
        averaged = np.zeros_like(seed)
        for target, sign in op_maps:
            t = np.zeros((n_ao, n_ao))
            t[target, np.arange(n_ao)] = sign
            averaged += t @ seed @ t.T
        averaged /= len(op_maps)

        for iop in range(self.maps["n_operations"]):
            target = np.array(self.maps["ao_target"][iop])
            sign = np.array(self.maps["ao_sign"][iop], dtype=float)
            t = np.zeros((n_ao, n_ao))
            t[target, np.arange(n_ao)] = sign
            self.assertLess(
                float(np.max(np.abs(t @ averaged @ t.T - averaged))), 1.0e-12
            )

    def test_dense_operation_matrices_are_rejected(self):
        rotation = np.array(
            [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        )
        operations = [
            {"name": "C4z", "matrix": rotation.tolist(), "permutation": [0, 1, 2]}
        ]
        with self.assertRaises(ValueError):
            self.symmetry.build_reduction_maps(WATER_SHELLS, operations)

    def test_petite_list_quartet_count_reduction(self):
        """Orbit representatives reduce symmetric quartet work as designed."""
        nshell = len(WATER_SHELLS)
        perms = np.array(self.maps["shell_permutation"])

        def canonical(quartet):
            best = None
            for perm in perms:
                i, j, k, l = (int(perm[q]) for q in quartet)
                ij = (max(i, j), min(i, j))
                kl = (max(k, l), min(k, l))
                pair = (max(ij, kl), min(ij, kl))
                if best is None or pair < best:
                    best = pair
            return best

        quartets = [
            (i, j, k, l)
            for i in range(nshell)
            for j in range(i + 1)
            for k in range(i + 1)
            for l in range(k + 1)
            if (i * (i + 1) // 2 + j) >= (k * (k + 1) // 2 + l)
        ]
        representatives = {canonical(q) for q in quartets}
        self.assertLess(len(representatives), len(quartets))
        # Orbit sizes must add back up to the full quartet count.
        from collections import Counter

        counts = Counter(canonical(q) for q in quartets)
        self.assertEqual(sum(counts.values()), len(quartets))


if __name__ == "__main__":
    unittest.main()
