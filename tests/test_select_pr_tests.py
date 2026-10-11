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

    def test_changed_test_module_selects_the_tests_that_import_it(self):
        chosen = self.sel.select(["tests/test_symmetry_metadata.py"])
        self.assertIsNotNone(chosen)
        for importer in ("tests/test_full_group_tier_declines_for_dft.py",
                         "tests/test_response_symmetry_staging_gates.py"):
            self.assertIn(importer, chosen)

    def test_python_module_change_runs_the_full_suite(self):
        # runfunc is reached through Runner (oqp.pyoqp) in test_casci.py,
        # which never names it.
        for path in ("pyoqp/oqp/library/runfunc.py", "pyoqp/oqp/utils/input_checker.py",
                     "pyoqp/requirements.txt"):
            with self.subTest(path=path):
                self.assertIsNone(self.sel.select([path]))

    def test_native_source_change_runs_the_full_suite(self):
        for path in ("source/ecp.F90", "source/tdhf_hessian.F90", "include/oqp.h"):
            with self.subTest(path=path):
                self.assertIsNone(self.sel.select([path]))

    def test_deleted_test_file_is_not_selected(self):
        chosen = self.sel.select(["tests/test_file_that_does_not_exist.py"])
        self.assertIsNotNone(chosen)
        self.assertNotIn("tests/test_file_that_does_not_exist.py", chosen)

    def test_selection_above_the_cap_falls_back_to_the_full_suite(self):
        self.sel.MAX_SELECTED = 1
        self.assertIsNone(self.sel.select(["tests/test_tdhf_input_checker.py"]))

    def test_unclassified_path_runs_the_full_suite(self):
        for path in ("Dockerfile", "docker/base-images.lock.json", "examples/TDDFT/x.inp"):
            with self.subTest(path=path):
                self.assertIsNone(self.sel.select([path]))

    def test_changed_document_selects_the_tests_that_read_it(self):
        for path, reader in (("README.md", "tests/test_readme_log_metadata.py"),
                             ("docs/logging.md", "tests/test_log_format.py"),
                             (".gitignore", "tests/test_dftd4_shared_interface.py")):
            with self.subTest(path=path):
                chosen = self.sel.select([path])
                self.assertIsNotNone(chosen)
                self.assertIn(reader, chosen)

    def test_document_no_test_reads_runs_just_the_smoke_set(self):
        # Built at run time so no test source (this one included) names it.
        unread = "docs/" + "-".join(("zz", "note", "no", "test", "reads")) + ".md"
        chosen = self.sel.select([unread])
        self.assertIsNotNone(chosen)
        self.assertEqual(sorted(chosen), sorted(s for s in self.sel.SMOKE if (ROOT / s).exists()))

    def test_tests_subdirectory_file_selects_its_driver(self):
        chosen = self.sel.select(["tests/fortran/test_guess_orbital_phase.F90"])
        if chosen is not None:
            self.assertIn("tests/test_guess_orbital_phase.py", chosen)

    def test_every_grd2_digest_source_selects_the_spherical_basis_detector(self):
        # AGENTS.md rule 6: tests/test_mrsf_nac_spherical_basis.py is the
        # detector for every grd2 digest, whether or not a test names the file.
        import re
        digest = re.compile(r"extends\s*\(\s*grd2_\w+_t\s*\)", re.IGNORECASE)
        sources = sorted(p for p in (ROOT / "source").rglob("*.F90")
                         if digest.search(p.read_text(errors="ignore")))
        self.assertTrue(sources)
        for path in [ROOT / "source/integrals/grd2.F90", *sources]:
            rel = str(path.relative_to(ROOT))
            with self.subTest(path=rel):
                chosen = self.sel.select([rel])
                if chosen is not None:
                    self.assertIn("tests/test_mrsf_nac_spherical_basis.py", chosen)

    def test_unknown_base_falls_back_to_the_full_suite(self):
        self.assertIsNone(self.sel.changed_files("0" * 40))


if __name__ == "__main__":
    unittest.main()
