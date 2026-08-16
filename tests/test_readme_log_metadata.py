import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ReadmeLogMetadataTests(unittest.TestCase):
    def test_readme_web_copy_omits_phase_one(self):
        readme = (ROOT / "README.md").read_text()

        self.assertNotIn("Phase " + "1", readme)

    def test_requested_contributors_are_listed(self):
        readme = (ROOT / "README.md").read_text()

        self.assertIn("[VladimirMakhnev](https://github.com/VladimirMakhnev)", readme)
        self.assertIn("[Alireza Lashkaripour](https://github.com/Alireza-Lashkaripour)", readme)

    def test_log_banner_lists_requested_contributors(self):
        banner = (ROOT / "source" / "modules" / "oqp_banner.F90").read_text()

        self.assertIn("Vladimir Makhnev", banner)
        self.assertIn("Alireza Lashkaripour", banner)

    def test_package_metadata_lists_alireza_without_invented_contact_data(self):
        metadata = (ROOT / "pyproject.toml").read_text()

        self.assertIn('{ name = "Alireza Lashkaripour" }', metadata)
        self.assertNotRegex(
            metadata,
            r'Alireza Lashkaripour"\s*,\s*email\s*=',
        )

    def test_mrsf_methods_are_presented_first(self):
        readme = (ROOT / "README.md").read_text()
        methods = readme.split("#### Electronic-Structure Methods", 1)[1]
        methods = methods.split("#### Properties & Spectroscopy", 1)[0]
        method_rows = [
            line.split("|", 2)[1].strip()
            for line in methods.splitlines()
            if line.startswith("|") and not line.startswith("| ---")
        ]

        self.assertEqual(
            method_rows[1:4],
            ["**MRSF-TDDFT**", "UMRSF-TDDFT", "MRSF-EKT"],
        )
        tutorials = next(
            line for line in methods.splitlines()
            if line.startswith("**Tutorials:**")
        )
        self.assertTrue(tutorials.startswith("**Tutorials:** [MRSF-TDDFT]"))

    def test_first_log_section_prints_git_head(self):
        runner = (ROOT / "pyoqp" / "oqp" / "pyoqp.py").read_text()
        file_utils = (ROOT / "pyoqp" / "oqp" / "utils" / "file_utils.py").read_text()

        self.assertIn("def _openqp_build_label()", runner)
        self.assertIn("git HEAD", runner)
        self.assertIn("section='start'", runner)
        self.assertIn('info={"build": _openqp_build_label()}', runner)
        self.assertIn("PyOQP build: %s", file_utils)


if __name__ == "__main__":
    unittest.main()
