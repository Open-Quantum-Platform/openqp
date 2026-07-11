import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "openqp_json_utils", ROOT / "pyoqp/oqp/utils/json_utils.py"
)
JSON_UTILS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(JSON_UTILS)


class TestJsonTdLayout(unittest.TestCase):
    def test_td_response_vectors_are_written_as_drf_by_state(self):
        # Fortran buffer for three DRFs and two states: state 1, then state 2.
        bridged = np.array([[11.0, 12.0], [13.0, 21.0], [22.0, 23.0]])

        result = JSON_UTILS.json_array("OQP::td_bvec_mo", bridged)

        self.assertEqual(result, [[11.0, 21.0], [12.0, 22.0], [13.0, 23.0]])

    def test_unrelated_matrices_keep_their_layout(self):
        matrix = np.array([[1.0, 2.0], [3.0, 4.0]])

        self.assertEqual(JSON_UTILS.json_array("OQP::SM", matrix), matrix.tolist())


if __name__ == "__main__":
    unittest.main()
