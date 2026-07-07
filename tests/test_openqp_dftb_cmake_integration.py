"""Lightweight checks for the optional openqp-dftb external hook."""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


class OpenQPDFTBCMakeIntegrationTests(unittest.TestCase):
    def test_top_level_cmake_defines_openqp_dftb_option(self):
        text = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("option(ENABLE_OPENQP_DFTB", text)
        self.assertIn("OPENQP_DFTB_GIT_REPOSITORY", text)
        self.assertIn("Open-Quantum-Platform/openqp-dftb.git", text)
        self.assertIn("OPENQP_DFTB_GIT_TAG", text)
        self.assertIn("OPENQP_DFTB_SOURCE_DIR", text)

    def test_external_project_can_use_git_or_local_source(self):
        text = (ROOT / "external" / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("if(ENABLE_OPENQP_DFTB)", text)
        self.assertIn("ExternalProject_Add(openqp_dftb", text)
        self.assertIn("OPENQP_DFTB_SOURCE_DIR", text)
        self.assertIn("OPENQP_DFTB_GIT_REPOSITORY", text)
        self.assertIn("OPENQP_DFTB_GIT_TAG", text)
        self.assertIn("-DOPENQP_DFTB_BUILD_TESTS=OFF", text)
        self.assertIn("OQP_EXTERNAL_GENERATOR_ARGS", text)
        self.assertIn("OPENQP_DFTB_INSTALL_LIBDIR", text)
        self.assertIn("-DCMAKE_INSTALL_LIBDIR=${OPENQP_DFTB_INSTALL_LIBDIR}", text)
        self.assertIn("OPENQP_DFTB_INCLUDE_DIR", text)
        self.assertIn("OPENQP_DFTB_LIBRARY", text)

    def test_oqp_links_openqp_dftb_only_when_enabled(self):
        text = (ROOT / "source" / "CMakeLists.txt").read_text(encoding="utf-8")
        match = re.search(r"if\(ENABLE_OPENQP_DFTB\)(.*?)endif\(\)", text, re.S)
        self.assertIsNotNone(
            match, "source/CMakeLists.txt must guard openqp-dftb behind ENABLE_OPENQP_DFTB"
        )
        block = match.group(1)
        self.assertIn("OQP_ENABLE_OPENQP_DFTB", block)
        self.assertIn("add_dependencies(oqp openqp_dftb)", block)
        self.assertIn("target_include_directories(oqp PUBLIC ${OPENQP_DFTB_INCLUDE_DIR})", block)
        self.assertIn("target_link_libraries(oqp ${OPENQP_DFTB_LIBRARY})", block)

    def test_openqp_repo_contains_only_thin_dftb_bridge(self):
        source_files = [path.name for path in (ROOT / "source").glob("openqp_dftb*.F90")]
        self.assertEqual(["openqp_dftb_bridge.F90"], source_files)
        bridge = (ROOT / "source" / "openqp_dftb_bridge.F90").read_text(encoding="utf-8")
        self.assertIn("#ifdef OQP_ENABLE_OPENQP_DFTB", bridge)
        self.assertIn("openqp_dftb_api", bridge)
        self.assertIn('bind(C, name="oqp_dftb_state_gradient")', bridge)

    def test_dftb_disabled_stub_uses_unconditional_status_writer(self):
        bridge = (ROOT / "source" / "openqp_dftb_bridge.F90").read_text(encoding="utf-8")
        writer_pos = bridge.index("subroutine write_c_string")
        guard_pos = bridge.index("#ifdef OQP_ENABLE_OPENQP_DFTB", bridge.index("end subroutine oqp_dftb_state_gradient_C"))
        self.assertLess(writer_pos, guard_pos)

    def test_cffi_header_declares_dftb_bridge(self):
        header = (ROOT / "include" / "oqp.h").read_text(encoding="utf-8")
        self.assertIn("void oqp_dftb_state_gradient", header)
        self.assertIn("status=-9001", header)


if __name__ == "__main__":
    unittest.main()
