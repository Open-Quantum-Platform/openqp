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
        self.assertIn("source/solvent_ddx_adapter.c", text)

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

    def test_mapping_doc_records_ddx_fock_uncertainty(self):
        text = (ROOT / "docs" / "solvent_ddx_scf_integration_seam.md").read_text(encoding="utf-8")
        self.assertIn("state%q", text)
        self.assertIn("ddx_get_xi` projects ddPCM `state%q", text)
        self.assertIn("Fock/Kohn-Sham operator", text)
        self.assertIn("unweighted total solute electrostatic potential", text)
        self.assertIn("cavity-projected `state%q`", text)

    def test_mapping_doc_records_reviewed_payload_handoff_chain(self):
        text = (ROOT / "docs" / "solvent_ddx_scf_integration_seam.md").read_text(encoding="utf-8")
        self.assertIn("reference_scf_pcm_runtime_payload", text)
        self.assertIn("reference_scf_pcm_calc_fock_handoff_from_molecule", text)
        self.assertIn("reference_scf_pcm_calc_fock_request", text)
        self.assertIn("non-incremental Fock path", text)
        self.assertIn("reference PCM incremental Fock is not validated", text)
        self.assertIn("pcm_reaction_potential_in", text)
        self.assertIn("runtime PCM remains disabled", text)
        self.assertIn("no state-specific or nonequilibrium MRSF solvent response", text)


if __name__ == "__main__":
    unittest.main()
