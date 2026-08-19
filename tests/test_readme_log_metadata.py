import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def markdown_table_rows(text, start_heading, end_heading):
    block = text.split(start_heading, 1)[1].split(end_heading, 1)[0]
    rows = []
    for line in block.splitlines():
        if not line.startswith("|") or line.startswith("| ---"):
            continue
        rows.append([cell.strip() for cell in line.strip("|").split("|")])
    return rows


class ReadmeLogMetadataTests(unittest.TestCase):
    def test_readme_web_copy_omits_phase_one(self):
        readme = (ROOT / "README.md").read_text()

        self.assertNotIn("Phase " + "1", readme)

    def test_opening_sentence_retains_the_core_program_scope(self):
        readme = (ROOT / "README.md").read_text()
        opening = readme.split("\n\n", 2)[1]

        for scope in (
            "MRSF-TDDFT",
            "correlated wavefunction methods",
            "linear-response excited states",
            "nuclear derivatives",
            "reaction paths",
            "nonadiabatic dynamics",
            "QM/MM calculations",
        ):
            self.assertIn(scope, opening)

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
        method_rows = markdown_table_rows(
            readme,
            "#### Electronic-Structure Methods",
            "#### Capabilities",
        )

        self.assertEqual(
            [row[0] for row in method_rows[1:4]],
            ["**MRSF-TDDFT**", "**UMRSF-TDDFT**", "**MRSF-EKT**"],
        )

    def test_method_rows_have_consistent_content_and_learning_resources(self):
        readme = (ROOT / "README.md").read_text()
        method_rows = markdown_table_rows(
            readme,
            "#### Electronic-Structure Methods",
            "#### Capabilities",
        )

        self.assertEqual(
            method_rows[0],
            ["Method", "References / variants", "Available calculations", "Learn"],
        )
        for row in method_rows[1:]:
            self.assertEqual(len(row), 4, row)
            self.assertGreaterEqual(len(row[2]), 50, row)
            self.assertIn("](", row[3], row)
        methods = {row[0].replace("**", "") for row in method_rows[1:]}
        self.assertEqual(
            methods,
            {
                "MRSF-TDDFT",
                "UMRSF-TDDFT",
                "MRSF-EKT",
                "SF-TDDFT",
                "TDHF / TDDFT",
                "Hartree–Fock",
                "Density functional theory",
                "MP2",
                "Coupled cluster",
                "Configuration interaction",
                "CASSCF",
                "Multireference PT2",
            },
        )

    def test_readme_reports_current_multireference_gradients(self):
        readme = (ROOT / "README.md").read_text()

        self.assertIn("analytic state-specific CASSCF gradients", readme)
        self.assertIn("central-difference SA-CASSCF gradients", readme)
        self.assertIn("SA-CASSCF and CASPT2/NEVPT2/QDPT2 families", readme)

    def test_capabilities_include_learning_resources_and_geometry_workflows(self):
        readme = (ROOT / "README.md").read_text()
        capability_rows = markdown_table_rows(
            readme,
            "#### Capabilities",
            "#### Ecosystem & Integrations",
        )

        for row in capability_rows[1:]:
            self.assertEqual(len(row), 4, row)
            self.assertIn("](", row[3], row)
        capabilities = [row[0] for row in capability_rows[1:]]
        self.assertIn("Geometry optimization and transition states", capabilities)
        self.assertIn("Conical intersections", capabilities)
        self.assertIn("Reaction paths", capabilities)
        self.assertNotIn("#### Geometry & Reaction Paths", readme)

    def test_detailed_controls_are_delegated_to_the_manual(self):
        readme = (ROOT / "README.md").read_text()

        self.assertNotIn("#### SCF, Initial Guesses & Performance", readme)
        self.assertNotIn("### Upcoming Features", readme)
        self.assertNotIn("Omitting the selector", readme)
        self.assertIn("openqp-docs/keywords/tests/", readme)
        self.assertIn("openqp-docs/keywords/scf/", readme)
        self.assertIn("openqp-docs/keywords/guess/", readme)
        self.assertIn("openqp-docs/performance/", readme)

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
