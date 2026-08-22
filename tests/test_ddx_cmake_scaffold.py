"""Lightweight checks for optional ddX CMake scaffolding."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class DDXCMakeScaffoldTests(unittest.TestCase):
    def test_top_level_cmake_defines_optional_ddx_backend(self):
        text = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("option(ENABLE_DDX", text)
        self.assertNotIn("ENABLE_DDX is not supported on Windows", text)
        self.assertIn("find_package(DDX REQUIRED)", text)
        self.assertIn("oqp_ddx_link_smoke", text)
        self.assertIn("oqp_ddx_adapter_smoke", text)
        self.assertIn("source/solvent_ddx_adapter.c", text)

    def test_find_ddx_module_creates_imported_target(self):
        text = (ROOT / "cmake" / "FindDDX.cmake").read_text(encoding="utf-8")
        self.assertIn("find_path(DDX_INCLUDE_DIR", text)
        self.assertIn("find_library(DDX_LIBRARY", text)
        self.assertIn("DDX_RUNTIME_LIBRARY", text)
        self.assertIn("IMPORTED_IMPLIB", text)
        self.assertIn("add_library(DDX::ddx SHARED IMPORTED)", text)
        self.assertIn("DDX::ddx", text)

    def test_oqp_links_ddx_only_when_enabled(self):
        text = (ROOT / "source" / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("if(ENABLE_DDX)", text)
        self.assertIn("OQP_ENABLE_DDX", text)
        self.assertIn("target_link_libraries(oqp DDX::ddx)", text)
        self.assertIn(
            "install(FILES ${DDX_RUNTIME_LIBRARY} DESTINATION ${CMAKE_INSTALL_LIBDIR})",
            text,
        )

    def test_autobuild_registers_distinct_windows_runtime_once(self):
        text = (ROOT / "external" / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn('set(_OQP_DDX_CACHE_NAME "ddx-runtime6")', text)
        self.assertIn("oqp_set_external_paths(DDX ${_OQP_DDX_CACHE_NAME}", text)
        self.assertIn("-DCMAKE_INSTALL_LIBDIR=lib", text)
        self.assertIn("-DDDX_DRIVER=OFF", text)
        patch = (ROOT / "cmake" / "patches" / "ddx-v0.8.0-ilp64.patch").read_text(
            encoding="utf-8"
        )
        self.assertIn('option(DDX_DRIVER "Build the ddX command-line driver" ON)', patch)
        self.assertIn("ARCHIVE DESTINATION lib", patch)
        self.assertIn("if(NOT DDX_RUNTIME_LIBRARY STREQUAL DDX_LIBRARY)", text)
        self.assertIn("BUILD_BYPRODUCTS ${_ddx_build_byproducts}", text)
        self.assertIn("INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR>", text)
        self.assertIn("add_library(DDX::ddx SHARED IMPORTED GLOBAL)", text)

    def test_oqp_owned_ddx_adapter_api_exists(self):
        header = (ROOT / "source" / "solvent_ddx_adapter.h").read_text(encoding="utf-8")
        source = (ROOT / "source" / "solvent_ddx_adapter.c").read_text(encoding="utf-8")
        self.assertIn("oqp_ddx_smoke_result_t", header)
        self.assertIn("oqp_ddx_run_point_charge_smoke", header)
        self.assertIn("oqp_ddx_run_explicit_pcm_smoke", header)
        self.assertIn("oqp_ddx_run_explicit_pcm_reaction_field_smoke", header)
        self.assertIn("cavity_xyz", header)
        self.assertIn("q_cav", header)
        self.assertIn("q_cav_norm", header)
        self.assertIn("q_cav_fd_derivative", header)
        self.assertIn("q_cav_fd_direct_abs_error", header)
        self.assertIn("q_cav_fd_abs_error", header)
        self.assertIn("oqp_ddx_run_point_charge_smoke", source)
        self.assertIn("oqp_ddx_run_explicit_pcm_smoke", source)
        self.assertIn("oqp_ddx_run_explicit_pcm_reaction_field_smoke", source)
        self.assertIn("ddx_pcm_setup", source)
        self.assertIn("ddx_pcm_solve", source)
        self.assertIn("ddx_pcm_solve_adjoint", source)
        self.assertIn("q_cav_out", source)
        self.assertIn("cavity_xyz_out", source)
        self.assertIn("finite_difference_delta", source)
        self.assertIn("q_cav_fd_derivative", source)
        self.assertIn("q_cav_norm", source)
        self.assertIn("ddx_get_xi", source)
        self.assertIn("ddx_get_cavity", source)
        self.assertIn("OQP_PI", source)
        self.assertNotIn("M_PI", source)


if __name__ == "__main__":
    unittest.main()
