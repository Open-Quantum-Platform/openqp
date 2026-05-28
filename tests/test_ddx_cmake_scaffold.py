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
        self.assertIn("reference_scf_pcm_calc_fock_request_from_scf_state", text)
        self.assertIn("dens_old", text)
        self.assertIn("f_old", text)
        self.assertIn("non-incremental Fock path", text)
        self.assertIn("reference PCM incremental Fock is not validated", text)
        self.assertIn("pcm_reaction_potential_in", text)
        self.assertIn("runtime PCM remains disabled", text)
        self.assertIn("no state-specific or nonequilibrium MRSF solvent response", text)

    def test_validation_matrix_records_payload_shape_metadata_contract(self):
        text = (ROOT / "docs" / "solvent_pcm_validation_matrix.md").read_text(encoding="utf-8")
        self.assertIn("reference_scf_pcm_runtime_payload()", text)
        self.assertIn("reference_scf_pcm_reaction_potential_from_payload()", text)
        self.assertIn("packed_ao_length", text)
        self.assertIn("expected_packed_ao_length", text)
        self.assertIn("packed_ao_shape_formula", text)
        self.assertIn("consumer must recompute and validate", text)
        self.assertIn("before exposing `pcm_reaction_potential_in`", text)

    def test_validation_matrix_records_backend_validation_status_payload_gate(self):
        text = (ROOT / "docs" / "solvent_pcm_validation_matrix.md").read_text(encoding="utf-8")
        self.assertIn("backend_validation_status", text)
        self.assertIn("pending PySCF/ddX/reference cross-check", text)
        self.assertIn("required before exposing `pcm_reaction_potential_in`", text)

    def test_validation_matrix_records_old_buffer_provenance_for_incremental_fock_guard(self):
        text = (ROOT / "docs" / "solvent_pcm_validation_matrix.md").read_text(encoding="utf-8")
        self.assertIn("dens_old_present", text)
        self.assertIn("f_old_present", text)
        self.assertIn("incremental_trigger_fields=dens_old,f_old", text)
        self.assertIn("incremental_trigger_fields=dens_old", text)
        self.assertIn("incremental_trigger_fields=f_old", text)

    def test_mapping_doc_records_payload_consumer_required_fields_and_nonmapping_guard(self):
        text = (ROOT / "docs" / "solvent_ddx_scf_integration_seam.md").read_text(encoding="utf-8")
        self.assertIn("PCM runtime payload must be a mapping", text)
        self.assertIn("required, not optional defaults", text)
        self.assertIn("nbf", text)
        self.assertIn("packed_ao_length", text)
        self.assertIn("expected_packed_ao_length", text)
        self.assertIn("packed_ao_shape_formula", text)
        self.assertIn("boolean numeric", text)
        self.assertIn("backend_validation_status", text)

    def test_validation_matrix_records_opt_in_calc_jk_xc_prototype_boundary(self):
        text = (ROOT / "docs" / "solvent_pcm_validation_matrix.md").read_text(encoding="utf-8")
        self.assertIn("calc_jk_xc", text)
        self.assertIn("opt-in/prototype-only", text)
        self.assertIn("optional `pcm_reaction_potential(:)`", text)
        self.assertIn("after `hcore`", text)
        self.assertIn("before SCF energy accumulation", text)
        self.assertIn("keep `pcm%enabled` out", text)

    def test_validation_matrix_records_molecule_payload_nonmapping_guard(self):
        text = (ROOT / "docs" / "solvent_pcm_validation_matrix.md").read_text(encoding="utf-8")
        self.assertIn("Molecule.set_pcm_runtime_payload()", text)
        self.assertIn("Molecule-level setter", text)
        self.assertIn("restored/prototype payloads", text)
        self.assertIn("PCM runtime payload must be a mapping", text)
        self.assertIn("before applying the allowlist", text)

    def test_ddx_seam_doc_records_final_calc_fock_call_site_bridge(self):
        text = (ROOT / "docs" / "solvent_ddx_scf_integration_seam.md").read_text(encoding="utf-8")
        self.assertIn("reference_scf_pcm_calc_fock_call_site_bridge()", text)
        self.assertIn("forward_pcm_reaction_potential", text)
        self.assertIn("disabled/no-payload", text)
        self.assertIn("pcm_reaction_potential_in", text)
        self.assertIn("non-incremental path", text)

    def test_ddx_seam_doc_records_molecule_payload_roundtrip_allowlist(self):
        text = (ROOT / "docs" / "solvent_ddx_scf_integration_seam.md").read_text(encoding="utf-8")
        self.assertIn("Molecule.set_pcm_runtime_payload()", text)
        self.assertIn("JSON round-trip persistence", text)
        self.assertIn("explicit allowlist", text)
        self.assertIn("backend_validation_status", text)
        self.assertIn("must not be stripped", text)
        self.assertIn("before exposing `pcm_reaction_potential_in`", text)

    def test_ddx_seam_doc_records_call_site_shape_audit_metadata(self):
        text = (ROOT / "docs" / "solvent_ddx_scf_integration_seam.md").read_text(encoding="utf-8")
        bridge_section = text.split("reference_scf_pcm_calc_fock_call_site_bridge()", 1)[1].split(
            "This is intentionally not a production solvent switch", 1
        )[0]
        self.assertIn("expected_packed_ao_length", bridge_section)
        self.assertIn("packed_ao_shape_formula", bridge_section)
        self.assertIn("call_site_shape_validated", bridge_section)
        self.assertIn("len(pcm_reaction_potential_in) == expected_packed_ao_length", bridge_section)

    def test_validation_matrix_records_call_site_shape_audit_boolean_boundary(self):
        text = (ROOT / "docs" / "solvent_pcm_validation_matrix.md").read_text(encoding="utf-8")
        bridge_section = text.split("reference_scf_pcm_calc_fock_call_site_bridge()", 1)[1].split(
            "The native `calc_fock", 1
        )[0]
        self.assertIn("call_site_shape_validated=True", bridge_section)
        self.assertIn("call_site_shape_validated=False", bridge_section)
        self.assertIn("forwarded payloads", bridge_section)
        self.assertIn("disabled/no-payload", bridge_section)
        self.assertIn("expected_packed_ao_length", bridge_section)
        self.assertIn("pcm_reaction_potential_in", bridge_section)


if __name__ == "__main__":
    unittest.main()
