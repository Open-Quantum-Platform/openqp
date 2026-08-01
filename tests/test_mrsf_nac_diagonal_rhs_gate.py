import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
GATE_PATH = ROOT / "tools" / "nac_lagrangian" / "diagonal_rhs_gate.py"
SPEC = importlib.util.spec_from_file_location("diagonal_rhs_gate", GATE_PATH)
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)


class TestMrsfNacDiagonalRhsGate(unittest.TestCase):
    def test_rohf_dual_projection_has_three_native_blocks(self):
        matrix = np.arange(25.0).reshape(5, 5)
        packed = GATE.pack_rohf_dual(matrix, nocb=1, noca=3)
        expected = np.array(
            [
                matrix[1, 0] - matrix[0, 1],
                matrix[2, 0] - matrix[0, 2],
                matrix[3, 0] - matrix[0, 3],
                matrix[4, 0] - matrix[0, 4],
                matrix[3, 1] - matrix[1, 3],
                matrix[3, 2] - matrix[2, 3],
                matrix[4, 1] - matrix[1, 4],
                matrix[4, 2] - matrix[2, 4],
            ]
        )
        np.testing.assert_array_equal(packed, expected)

    def test_duplicate_uses_distinct_labels_with_equal_amplitudes(self):
        raw = np.arange(12.0).reshape(4, 3)
        duplicated, amplitude = GATE.duplicate_state_slots(
            raw, nstate=3, nij=4, state=2, auxiliary=0
        )
        flat = duplicated.reshape(-1)
        np.testing.assert_array_equal(amplitude, np.arange(8.0, 12.0))
        np.testing.assert_array_equal(flat[0:4], amplitude)
        np.testing.assert_array_equal(flat[8:12], amplitude)
        np.testing.assert_array_equal(flat[4:8], np.arange(4.0, 8.0))
        np.testing.assert_array_equal(raw, np.arange(12.0).reshape(4, 3))

    def test_duplicate_rejects_the_physical_diagonal_label(self):
        with self.assertRaisesRegex(ValueError, "must differ"):
            GATE.duplicate_state_slots(
                np.zeros((2, 3)), nstate=2, nij=3, state=1, auxiliary=1
            )

    def test_symmetric_pack_roundtrip(self):
        matrix = np.array(
            [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]
        )
        packed = GATE.pack_symmetric(matrix)
        np.testing.assert_array_equal(
            GATE.unpack_symmetric(packed, matrix.shape[0]), matrix
        )

    def test_diagonal_source_closes_with_lee_opposite_sign(self):
        ell_pair = np.array([1.5, -2.0, 0.25])
        rhs_lee = -ell_pair
        np.testing.assert_array_equal(
            GATE.diagonal_source_closure(ell_pair, rhs_lee),
            np.zeros_like(ell_pair),
        )
        with self.assertRaisesRegex(ValueError, "same shape"):
            GATE.diagonal_source_closure(ell_pair, rhs_lee[:2])


if __name__ == "__main__":
    unittest.main()
