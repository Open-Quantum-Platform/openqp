import importlib.util
import ast
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

    def test_td_response_vectors_round_trip_back_to_tag_bridge_layout(self):
        bridged = np.array([[11.0, 12.0], [13.0, 21.0], [22.0, 23.0]])

        documented = JSON_UTILS.json_array("OQP::td_bvec_mo", bridged)
        restored = JSON_UTILS.tag_array_from_json(
            "OQP::td_bvec_mo", documented
        )

        np.testing.assert_array_equal(restored, bridged)

    def test_unrelated_matrices_keep_their_layout(self):
        matrix = np.array([[1.0, 2.0], [3.0, 4.0]])

        self.assertEqual(JSON_UTILS.json_array("OQP::SM", matrix), matrix.tolist())

    def test_conversion_is_confined_to_json_save_path(self):
        source = (ROOT / "pyoqp/oqp/molecule/molecule.py").read_text()
        tree = ast.parse(source)
        methods = {
            node.name: node
            for cls in tree.body
            if isinstance(cls, ast.ClassDef) and cls.name == "Molecule"
            for node in cls.body
            if isinstance(node, ast.FunctionDef)
        }

        get_data_calls = [
            node for node in ast.walk(methods["get_data"])
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id == "json_array"
        ]
        save_data_calls = [
            node for node in ast.walk(methods["save_data"])
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id == "json_array"
        ]

        self.assertEqual(get_data_calls, [])
        self.assertEqual(len(save_data_calls), 1)

    def test_full_json_preserves_mrsf_state_interaction_data(self):
        source = (ROOT / "pyoqp/oqp/molecule/molecule.py").read_text()
        tree = ast.parse(source)
        molecule = next(
            node for node in tree.body
            if isinstance(node, ast.ClassDef) and node.name == "Molecule"
        )
        init = next(
            node for node in molecule.body
            if isinstance(node, ast.FunctionDef) and node.name == "__init__"
        )
        tag_assignment = next(
            node for node in ast.walk(init)
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Attribute) and target.attr == "tag"
                for target in node.targets
            )
        )
        tags = {item.value for item in tag_assignment.value.elts}

        self.assertTrue({
            "OQP::td_trans_density_mo",
            "OQP::td_trans_dipole",
            "OQP::td_dip_ao",
        }.issubset(tags))


if __name__ == "__main__":
    unittest.main()
