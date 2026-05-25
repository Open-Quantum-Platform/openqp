"""Lightweight checks for optional ddX CMake scaffolding."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class DDXCMakeScaffoldTests(unittest.TestCase):
    def test_top_level_cmake_defines_optional_ddx_backend(self):
        text = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("option(ENABLE_DDX", text)
        self.assertIn("find_package(DDX REQUIRED)", text)
        self.assertIn("oqp_ddx_link_smoke", text)
        self.assertIn("oqp_ddx_adapter_smoke", text)

    def test_find_ddx_module_creates_imported_target(self):
        text = (ROOT / "cmake" / "FindDDX.cmake").read_text(encoding="utf-8")
        self.assertIn("find_path(DDX_INCLUDE_DIR", text)
        self.assertIn("find_library(DDX_LIBRARY", text)
        self.assertIn("DDX::ddx", text)

    def test_oqp_links_ddx_only_when_enabled(self):
        text = (ROOT / "source" / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("if(ENABLE_DDX)", text)
        self.assertIn("OQP_ENABLE_DDX", text)
        self.assertIn("target_link_libraries(oqp DDX::ddx)", text)


if __name__ == "__main__":
    unittest.main()
