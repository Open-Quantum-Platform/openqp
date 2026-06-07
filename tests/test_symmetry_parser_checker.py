import importlib.util
import re
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_input_checker():
    oqp_mod = sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    utils_mod = sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    setattr(oqp_mod, "utils", utils_mod)
    setattr(utils_mod, "mpi_utils", mpi_utils)

    path = ROOT / "pyoqp/oqp/utils/input_checker.py"
    spec = importlib.util.spec_from_file_location("openqp_symmetry_checker_under_test", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class TestSymmetryParserAndMetadataGates(unittest.TestCase):
    def setUp(self):
        self.input_checker = load_input_checker()

    def test_schema_file_declares_symmetry_section_with_safe_defaults(self):
        text = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()

        self.assertIn("'symmetry': {", text)
        self.assertIn("'enabled': {'type': string, 'default': 'false'}", text)
        self.assertIn("'point_group': {'type': string, 'default': 'auto'}", text)
        self.assertIn("'subgroup': {'type': string, 'default': 'auto'}", text)
        self.assertIn("'label_mo': {'type': bool, 'default': 'True'}", text)
        self.assertIn("'label_states': {'type': bool, 'default': 'True'}", text)
        self.assertIn("'label_modes': {'type': bool, 'default': 'True'}", text)
        self.assertIn("'use_integral_symmetry': {'type': bool, 'default': 'False'}", text)
        self.assertIn("'use_response_symmetry': {'type': bool, 'default': 'False'}", text)
        self.assertIn("'tolerance': {'type': float, 'default': '1.0e-5'}", text)
        self.assertIn("'strict': {'type': bool, 'default': 'False'}", text)

    def test_schema_does_not_wire_symmetry_handlers_by_default(self):
        data = (ROOT / "pyoqp/oqp/molecule/oqpdata.py").read_text()

        self.assertIn('"symmetry": {', data)
        self.assertNotIn('"symmetry": {\n        "enabled"', data[data.index('"symmetry"'):])

    def test_symmetry_check_rejects_invalid_enabled_value(self):
        report = self.input_checker.CheckReport()
        config = {"symmetry": {"enabled": "maybe", "label_mo": True, "label_states": True,
                                "label_modes": True, "tolerance": 1.0e-5,
                                "use_integral_symmetry": False,
                                "use_response_symmetry": False,
                                "strict": False}}

        self.input_checker._check_symmetry(config, report)

        messages = "\n".join(item.path for item in report.errors)
        self.assertIn("symmetry.enabled", messages)

    def test_symmetry_tolerance_must_be_positive(self):
        report = self.input_checker.CheckReport()
        config = {"symmetry": {"enabled": "auto", "label_mo": True, "label_states": True,
                                "label_modes": True, "tolerance": 0.0,
                                "use_integral_symmetry": False,
                                "use_response_symmetry": False,
                                "strict": False}}

        self.input_checker._check_symmetry(config, report)

        self.assertFalse(report.ok)
        self.assertIn("symmetry.tolerance", "\n".join(item.path for item in report.errors))

    def test_symmetry_reduction_flags_stay_off_in_metadata_gate(self):
        report = self.input_checker.CheckReport()
        config = {
            "symmetry": {
                "enabled": False,
                "point_group": "auto",
                "subgroup": "auto",
                "label_mo": True,
                "label_states": True,
                "label_modes": True,
                "tolerance": 1e-5,
                "use_integral_symmetry": "False",
                "use_response_symmetry": "False",
                "strict": False,
            }
        }

        self.input_checker._check_symmetry(config, report)

        self.assertTrue(report.ok, report.to_text())
        self.assertNotIn("symmetry.use_integral_symmetry", "\n".join(item.path for item in report.errors))
        self.assertNotIn("symmetry.use_response_symmetry", "\n".join(item.path for item in report.errors))

    def test_reduction_flag_policy(self):
        """Integral reduction is experimental (warn); response stays rejected."""
        report = self.input_checker.CheckReport()
        config = {
            "symmetry": {
                "enabled": False,
                "label_mo": True,
                "label_states": True,
                "label_modes": True,
                "tolerance": 1e-5,
                "use_integral_symmetry": "True",
                "use_response_symmetry": "yes",
                "strict": False,
            }
        }

        self.input_checker._check_symmetry(config, report)

        errors = "\n".join(item.path for item in report.errors)
        self.assertNotIn("symmetry.use_integral_symmetry", errors)
        self.assertIn("symmetry.use_response_symmetry", errors)
        warnings = "\n".join(item.path for item in report.warnings)
        self.assertIn("symmetry.use_integral_symmetry", warnings)

    def test_symmetry_metadata_files_stay_backend_free(self):
        root = ROOT
        targets = [
            root / "pyoqp/oqp/molecule/oqpdata.py",
            root / "pyoqp/oqp/molecule/molecule.py",
            root / "pyoqp/oqp/utils/input_checker.py",
            root / "pyoqp/oqp/library/symmetry.py",
        ]
        forbidden = ["libintx", "metc", "cuda", "opencl", "gpu"]

        for path in targets:
            source = path.read_text().lower()
            for term in forbidden:
                pattern = rf"\\b{re.escape(term)}\\b"
                self.assertIsNone(
                    re.search(pattern, source),
                    msg=f"Found {term} in backend-sensitive symmetry file {path}",
                )
