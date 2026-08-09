"""Dependency-light guards for the complete removal of NLopt."""

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class NLoptRemovalTests(unittest.TestCase):
    def test_source_module_and_scf_calls_are_absent(self):
        self.assertFalse((ROOT / "source" / "nlopt.F90").exists())
        scf = (ROOT / "source" / "scf_converger.F90").read_text()
        scf_lower = scf.lower()
        self.assertNotIn("use nlopt", scf_lower)
        self.assertNotIn("nlo_", scf_lower)
        self.assertNotIn("nlopt_", scf_lower)
        self.assertIn("use simplex_qp, only: solve_simplex_qp", scf_lower)
        self.assertIn("hqp = -0.5_dp*", scf_lower)
        self.assertIn("gqp = 2.0_dp*self%b(1:na)", scf_lower)
        self.assertIn("self%nlog", scf_lower)
        self.assertIn("self%solution_current", scf_lower)

    def test_cmake_does_not_fetch_include_or_link_nlopt(self):
        source_cmake = (ROOT / "source" / "CMakeLists.txt").read_text()
        self.assertNotIn("NLOPT_BUILD_DIR", source_cmake)
        self.assertNotIn("NLOPT_LIBRARY", source_cmake)

        external = (ROOT / "external" / "CMakeLists.txt").read_text()
        self.assertNotRegex(
            external, r"ExternalProject_Add\s*\(\s*nlopt", msg=external
        )
        self.assertNotRegex(
            external, r"oqp_reuse_or_build\s*\(\s*nlopt", msg=external
        )
        self.assertNotIn("_OQP_NLOPT_VERSION", external)
        self.assertNotRegex(external, r"\bNLOPT_(?!LEGACY)")

    def test_legacy_cache_token_preserves_the_existing_namespace(self):
        external = (ROOT / "external" / "CMakeLists.txt").read_text()
        self.assertIn(
            'set(_OQP_LEGACY_NLOPT_CACHE_TOKEN "nlopt2.9.1")', external
        )
        self.assertIn(
            "-libint${_OQP_LIBINT2_VERSION}-"
            "${_OQP_LEGACY_NLOPT_CACHE_TOKEN}-libxc",
            external,
        )

    def test_distribution_notice_for_removed_component_is_absent(self):
        self.assertFalse(
            (ROOT / "licenses" / "third_party" / "nlopt-notices.txt").exists()
        )
        notices = (ROOT / "THIRD_PARTY_NOTICES.md").read_text().lower()
        self.assertNotIn("nlopt-notices.txt", notices)
        self.assertIn("nlopt", notices)
        self.assertIn("removed", notices)

    def test_binary_artifact_guards_are_present(self):
        wheel = (ROOT / ".github" / "scripts" / "wheel_smoke_test.py").read_text()
        self.assertIn('"nlopt" not in oqp_deps.lower()', wheel)
        self.assertRegex(wheel, r"nlopt\|nlo_")

        docker = (ROOT / "Dockerfile").read_text()
        container_smoke = (
            ROOT / ".github" / "scripts" / "container_runtime_smoke.py"
        ).read_text()
        self.assertIn("container_runtime_smoke.py", docker)
        self.assertIn("ldd_dependencies", container_smoke)
        self.assertIn("NLopt dependency leaked", container_smoke)
        self.assertIn("NLopt symbols/strings leaked", container_smoke)
        self.assertIn('b"nlopt_"', container_smoke)
        self.assertIn('b"libnlopt"', container_smoke)


if __name__ == "__main__":
    unittest.main()
