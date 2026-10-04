"""Dependency-light guards for the complete removal of libecpint.

Effective-core-potential integrals and their first and second nuclear
derivatives come from source/ecp.F90 (validated by
tests/test_ecp_integrals.py).
"""

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class LibecpintRemovalTests(unittest.TestCase):
    def test_wrapper_sources_are_absent(self):
        self.assertFalse((ROOT / "source" / "ecpint.F90").exists())
        self.assertFalse((ROOT / "source" / "wrapper" / "libecpint_wrapper.cpp").exists())
        self.assertFalse((ROOT / "external" / "fix_libecpint_accuracy.py").exists())
        ecp = (ROOT / "source" / "ecp.F90").read_text().lower()
        self.assertNotIn("libecpint_wrapper", ecp)
        self.assertNotIn("libecp_result", ecp)
        self.assertNotIn("set_integrator", ecp)
        self.assertIn("subroutine ecp_raw_ints(", ecp)
        self.assertFalse((ROOT / "source" / "ecp_native.F90").exists())
        basis_api = (ROOT / "source" / "basis_api.F90").read_text().lower()
        self.assertNotIn("libecpint", basis_api)

    def test_cmake_does_not_fetch_include_or_link_libecpint(self):
        source_cmake = (ROOT / "source" / "CMakeLists.txt").read_text()
        for token in ("LIBECPINT", "FADDEEVA", "ecpint", "Faddeeva", "libecpint_wrapper"):
            self.assertNotIn(token, source_cmake)
        external = (ROOT / "external" / "CMakeLists.txt").read_text()
        self.assertNotRegex(external, r"ExternalProject_Add\s*\(\s*libecpint")
        self.assertNotRegex(external, r"oqp_reuse_or_build\s*\(\s*libecpint")
        self.assertNotIn("_OQP_LIBECPINT_VERSION", external)
        self.assertNotRegex(external, r"\bLIBECPINT_(?!CACHE)")

    def test_legacy_cache_token_preserves_the_existing_namespace(self):
        external = (ROOT / "external" / "CMakeLists.txt").read_text()
        self.assertIn('set(_OQP_LEGACY_LIBECPINT_CACHE_TOKEN "ecp1.0.7")', external)
        self.assertIn(
            "-tag${_OQP_TAGARRAY_VERSION}-${_OQP_LEGACY_LIBECPINT_CACHE_TOKEN}-lapack",
            external,
        )

    def test_distribution_notices_for_removed_component_are_absent(self):
        for name in ("libecpint-mit.txt", "faddeeva-mit.txt"):
            self.assertFalse((ROOT / "licenses" / "third_party" / name).exists())
        notices = (ROOT / "THIRD_PARTY_NOTICES.md").read_text().lower()
        self.assertNotIn("libecpint-mit.txt", notices)
        self.assertNotIn("faddeeva-mit.txt", notices)
        self.assertIn("libecpint", notices)
        self.assertIn("removed", notices)


if __name__ == "__main__":
    unittest.main()
