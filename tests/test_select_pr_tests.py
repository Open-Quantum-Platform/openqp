"""CI pull-request test selection (.github/scripts/select_pr_tests.py)."""
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_selector():
    spec = importlib.util.spec_from_file_location(
        "select_pr_tests_under_test", ROOT / ".github/scripts/select_pr_tests.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SelectPrTestsTests(unittest.TestCase):
    def setUp(self):
        self.sel = load_selector()

    def test_harness_build_and_workflow_changes_run_the_full_suite(self):
        for path in ("tests/conftest.py", "tests/oqp_module_isolation.py",
                     "tests/known_failures.txt", ".github/workflows/CI.yml",
                     "CMakeLists.txt", "pyproject.toml", "cmake/x.cmake",
                     "external/CMakeLists.txt"):
            with self.subTest(path=path):
                self.assertIsNone(self.sel.select([path]))

    def test_changed_test_file_and_smoke_set_are_selected(self):
        chosen = self.sel.select(["tests/test_tdhf_input_checker.py"])
        self.assertIsNotNone(chosen)
        self.assertIn("tests/test_tdhf_input_checker.py", chosen)
        for smoke in self.sel.SMOKE:
            if (ROOT / smoke).exists():
                self.assertIn(smoke, chosen)

    def test_python_module_change_selects_tests_that_name_it(self):
        chosen = self.sel.select(["pyoqp/oqp/utils/input_checker.py"])
        self.assertIsNotNone(chosen)
        self.assertIn("tests/test_tdhf_input_checker.py", chosen)

    def test_native_source_change_selects_tests_naming_the_stem(self):
        chosen = self.sel.select(["source/ecp.F90"])
        self.assertIsNotNone(chosen)
        self.assertIn("tests/test_ecp_integrals.py", chosen)

    def test_deleted_test_file_is_not_selected(self):
        chosen = self.sel.select(["tests/test_file_that_does_not_exist.py"])
        self.assertIsNotNone(chosen)
        self.assertNotIn("tests/test_file_that_does_not_exist.py", chosen)

    def test_selection_above_the_cap_falls_back_to_the_full_suite(self):
        self.sel.MAX_SELECTED = 1
        self.assertIsNone(self.sel.select(["tests/test_tdhf_input_checker.py"]))


if __name__ == "__main__":
    unittest.main()
