"""Checks for the ddX-to-SCF integration seam and the single canonical PCM path.

After the dual-path reconciliation there is exactly one runtime PCM path:

    infos%control%pcm_enabled  ->  add_pcm_reaction_field  ->  E%e_pcm

The retired reference-supplied-potential prototype (the optional
``pcm_reaction_potential``/``pcm_reaction_potential_in`` arguments, the separate
``epcm`` field, the ``add_reference_pcm_reaction_field`` hook, and the Python
reviewed-payload/``calc_fock`` handoff chain) has been removed. The helpers left
in ``pyoqp/oqp/library/solvent.py`` are validation/cross-check helpers only and
are not imported by the runtime SCF.
"""

from pathlib import Path
import importlib.util
import unittest


ROOT = Path(__file__).resolve().parents[1]


def _load_solvent(tag):
    module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
    spec = importlib.util.spec_from_file_location(f"solvent_under_test_{tag}", module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load {module_path}")
    solvent = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(solvent)
    return solvent


class DDXSCFIntegrationSeamTests(unittest.TestCase):
    # ----- retained validation/cross-check helpers in solvent.py --------------

    def test_provisional_ddx_q_cav_to_external_charges_is_opt_in(self):
        solvent = _load_solvent("q_cav")
        with self.assertRaisesRegex(ValueError, "provisional"):
            solvent.provisional_ddx_external_charges([0.2, -0.4])
        charges = solvent.provisional_ddx_external_charges([0.2, -0.4], allow_provisional=True)
        self.assertEqual(charges, [-0.1, 0.2])

    def test_provisional_ddx_reaction_field_inputs_validate_cavity_shape(self):
        solvent = _load_solvent("shapes")
        with self.assertRaisesRegex(ValueError, r"3 \* len\(q_cav\)"):
            solvent.provisional_ddx_reaction_field_inputs(
                [0.2, -0.4],
                [0.0, 0.0, 0.0, 1.0, 0.0],
                allow_provisional=True,
            )
        with self.assertRaisesRegex(ValueError, "q_cav must contain at least one"):
            solvent.provisional_ddx_reaction_field_inputs([], [], allow_provisional=True)
        inputs = solvent.provisional_ddx_reaction_field_inputs(
            [0.2, -0.4],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            allow_provisional=True,
        )
        self.assertEqual(inputs["charges"], [-0.1, 0.2])
        self.assertEqual(inputs["cavity_xyz"], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        self.assertEqual(inputs["ncav"], 2)

    def test_provisional_ddx_reaction_field_inputs_split_external_charge_arrays(self):
        solvent = _load_solvent("external_arrays")
        inputs = solvent.provisional_ddx_reaction_field_inputs(
            [0.2, -0.4],
            [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
            allow_provisional=True,
        )
        self.assertEqual(inputs["x"], [0.0, 1.0])
        self.assertEqual(inputs["y"], [0.1, 1.1])
        self.assertEqual(inputs["z"], [0.2, 1.2])
        self.assertEqual(inputs["chg"], [-0.1, 0.2])

    def test_provisional_ddx_reaction_field_inputs_label_sign_scale_and_runtime_scope(self):
        solvent = _load_solvent("provisional_scope")
        inputs = solvent.provisional_ddx_reaction_field_inputs(
            [0.2], [0.0, 0.1, 0.2], allow_provisional=True
        )
        self.assertTrue(inputs["provisional_sign_scale"])
        self.assertEqual(inputs["sign_scale_convention"], "chg = -0.5 * q_cav")
        self.assertEqual(
            inputs["validation_status"],
            "requires independent-reference/ddX cross-check before runtime use",
        )
        self.assertFalse(inputs["runtime_pcm_enabled"])

    def test_provisional_ddx_reaction_field_inputs_preserve_first_scope_metadata(self):
        solvent = _load_solvent("provisional_first_scope")
        inputs = solvent.provisional_ddx_reaction_field_inputs(
            [0.2], [0.0, 0.1, 0.2], allow_provisional=True
        )
        self.assertEqual(inputs["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(inputs["reference_target"], "RHF/ROHF reference density")
        self.assertEqual(inputs["response_solvent_coupling"], "not enabled")
        self.assertEqual(inputs["gradient_support"], "not enabled")
        self.assertEqual(inputs["backend_validation_status"], "pending independent-reference/ddX cross-check")

    def test_solvent_helpers_reject_nonfinite_numeric_inputs(self):
        solvent = _load_solvent("finite_values")
        with self.assertRaisesRegex(ValueError, "finite"):
            solvent.reference_scf_phi_cav_inputs([[1.0, float("nan"), 0.25]], [0.0, 0.1, 0.2])
        with self.assertRaisesRegex(ValueError, "finite"):
            solvent.provisional_ddx_reaction_field_inputs(
                [0.2, float("inf")],
                [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
                allow_provisional=True,
            )

    def test_reference_scf_total_density_sums_spin_blocks_for_ddx_phi_cav(self):
        solvent = _load_solvent("total_density")
        self.assertEqual(solvent.reference_scf_total_density([[1.0, 0.25, -0.5]]), [1.0, 0.25, -0.5])
        self.assertEqual(
            solvent.reference_scf_total_density([[1.0, 0.25, -0.5], [0.75, -0.25, 0.5]]),
            [1.75, 0.0, 0.0],
        )
        with self.assertRaisesRegex(ValueError, "same packed length"):
            solvent.reference_scf_total_density([[1.0, 0.25], [0.75]])
        with self.assertRaisesRegex(ValueError, "one or two density blocks"):
            solvent.reference_scf_total_density([[1.0], [0.5], [0.25]])
        with self.assertRaisesRegex(ValueError, "at least one density block"):
            solvent.reference_scf_total_density([])
        with self.assertRaisesRegex(ValueError, "packed density blocks must not be empty"):
            solvent.reference_scf_total_density([[]])

    def test_reference_scf_reaction_field_contract_validates_packed_matrix(self):
        solvent = _load_solvent("rf_contract")
        with self.assertRaisesRegex(ValueError, "same packed length"):
            solvent.reference_scf_reaction_field_contract([[1.0, 0.5, 0.25]], [0.1, 0.2])
        with self.assertRaisesRegex(ValueError, "triangular packed AO"):
            solvent.reference_scf_reaction_field_contract([[1.0, 0.5]], [0.1, 0.2])
        contract = solvent.reference_scf_reaction_field_contract(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.1, 0.2, 0.3],
        )
        self.assertEqual(contract["nbf"], 2)
        self.assertEqual(contract["nfocks"], 2)
        self.assertEqual(contract["total_density"], [1.75, 0.0, 0.25])
        self.assertEqual(contract["reaction_potential"], [0.1, 0.2, 0.3])

    def test_reference_scf_phi_cav_inputs_validate_density_and_cavity_points(self):
        solvent = _load_solvent("phi_cav")
        with self.assertRaisesRegex(ValueError, r"3 \* ncav"):
            solvent.reference_scf_phi_cav_inputs([[1.0, 0.5, 0.25]], [0.0, 0.1])
        with self.assertRaisesRegex(ValueError, "at least one cavity point"):
            solvent.reference_scf_phi_cav_inputs([[1.0, 0.5, 0.25]], [])
        inputs = solvent.reference_scf_phi_cav_inputs(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
        )
        self.assertEqual(inputs["nbf"], 2)
        self.assertEqual(inputs["ncav"], 2)
        self.assertEqual(inputs["total_density"], [1.75, 0.0, 0.25])
        self.assertEqual(inputs["density_packed"], [1.75, 0.0, 0.25])
        self.assertNotIn("density_blocks", inputs)
        self.assertEqual(inputs["x"], [0.0, 1.0])
        self.assertEqual(inputs["y"], [0.1, 1.1])
        self.assertEqual(inputs["z"], [0.2, 1.2])

    def test_reference_scf_pcm_energy_terms_use_validated_packed_density_and_potential(self):
        solvent = _load_solvent("energy_terms")
        with self.assertRaisesRegex(ValueError, "same packed length"):
            solvent.reference_scf_pcm_energy_terms([[1.0, 0.5, 0.25]], [0.1, 0.2])
        terms = solvent.reference_scf_pcm_energy_terms(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.1, 0.2, 0.3],
        )
        self.assertEqual(terms["nbf"], 2)
        self.assertEqual(terms["density_packed"], [1.75, 0.0, 0.25])
        self.assertEqual(terms["reaction_potential"], [0.1, 0.2, 0.3])
        self.assertAlmostEqual(terms["density_reaction_dot"], 0.25)
        self.assertAlmostEqual(terms["candidate_polarization_energy"], 0.125)
        self.assertEqual(terms["energy_convention"], "0.5 * dot(reference_density_packed, reaction_potential)")

    def test_reference_scf_pcm_energy_terms_preserve_energy_only_scope_metadata(self):
        solvent = _load_solvent("energy_scope")
        terms = solvent.reference_scf_pcm_energy_terms([[1.0, 0.5, 0.25]], [0.1, 0.2, 0.2])
        self.assertEqual(terms["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(terms["reference_target"], "RHF/ROHF reference density")
        self.assertEqual(terms["response_solvent_coupling"], "not enabled")
        self.assertEqual(terms["gradient_support"], "not enabled")
        self.assertFalse(terms["runtime_pcm_enabled"])
        self.assertEqual(terms["backend_validation_status"], "pending independent-reference/ddX cross-check")

    def test_provisional_ddx_reaction_field_inputs_requires_literal_boolean_opt_in(self):
        solvent = _load_solvent("provisional_opt_in")
        with self.assertRaisesRegex(ValueError, "allow_provisional must be the boolean True"):
            solvent.provisional_ddx_external_charges([0.2], allow_provisional="yes")
        with self.assertRaisesRegex(ValueError, "allow_provisional must be the boolean True"):
            solvent.provisional_ddx_reaction_field_inputs([0.2], [0.0, 0.0, 0.0], allow_provisional=1)
        self.assertEqual(solvent.provisional_ddx_external_charges([0.2], allow_provisional=True), [-0.1])

    # ----- unchanged int1.F90 / C adapter seam pieces -------------------------

    def test_unweighted_electrostatic_potential_is_public(self):
        text = (ROOT / "source" / "integrals" / "int1.F90").read_text(encoding="utf-8")
        self.assertIn("public electrostatic_potential_unweighted", text)
        self.assertIn("subroutine electrostatic_potential_unweighted", text.lower())
        wrapper = text.split("subroutine electrostatic_potential_unweighted", 1)[1].split(
            "end subroutine electrostatic_potential_unweighted", 1
        )[0]
        self.assertIn("call int1_el_pot", wrapper)
        self.assertNotIn("pot = pot*wt", wrapper)
        self.assertIn("invnrm = 1.0_real64 / basis%bfnrm", wrapper)
        self.assertIn("call bas_norm_matrix(d, invnrm", wrapper)

    def test_external_charge_potential_is_public(self):
        text = (ROOT / "source" / "integrals" / "int1.F90").read_text(encoding="utf-8")
        self.assertIn("public external_charge_potential", text)
        self.assertIn("subroutine external_charge_potential", text.lower())
        self.assertIn("call int1_coul_ext_chg", text)
        self.assertIn("call bas_norm_matrix(v, basis%bfnrm, basis%nbf)", text)

    def test_ddx_adapter_exposes_cavity_q_for_external_charge_potential(self):
        header = (ROOT / "source" / "solvent_ddx_adapter.h").read_text(encoding="utf-8")
        source = (ROOT / "source" / "solvent_ddx_adapter.c").read_text(encoding="utf-8")
        self.assertIn("oqp_ddx_run_explicit_pcm_reaction_field_smoke", header)
        self.assertIn("double* cavity_xyz", header)
        self.assertIn("double* q_cav", header)
        self.assertIn("oqp_ddx_run_explicit_pcm_reaction_field_smoke", source)
        self.assertIn("cavity_xyz_out", source)
        self.assertIn("q_cav_out", source)
        self.assertIn("q_cav_out[icav] = xi[icav]", source)
        self.assertIn("q_cav_fd_derivative", header)
        self.assertIn("q_cav_fd_direct_abs_error", header)
        self.assertIn("q_cav_fd_abs_error", header)
        self.assertIn("finite_difference_delta", source)
        self.assertIn("fabs(q_cav_fd_derivative", source)

    def test_scf_fock_builder_is_identified_for_solvent_hook(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        self.assertIn("subroutine calc_jk_xc", text.lower())
        self.assertIn("f(:,ii) =  f(:,ii) + hcore", text)
        self.assertIn("E%ehf1 = E%ehf1 + traceprod_sym_packed", text)

    # ----- single canonical runtime PCM path (dual-path reconciliation) -------

    def test_runtime_pcm_has_single_gate_hook_and_energy_field(self):
        """The canonical live path is pcm_enabled -> add_pcm_reaction_field -> e_pcm."""
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        calc_jk_xc = text.split("subroutine calc_jk_xc", 1)[1].split("end subroutine calc_jk_xc", 1)[0]
        # Exactly one gate, one hook, one energy field in the live Fock builder.
        self.assertIn("if (infos%control%pcm_enabled) then", calc_jk_xc)
        self.assertIn("call add_pcm_reaction_field(basis, infos, d, nfocks, f, E%e_pcm)", calc_jk_xc)
        self.assertIn("E%etot = E%etot + E%e_pcm", calc_jk_xc)
        self.assertEqual(calc_jk_xc.count("call add_pcm_reaction_field"), 1)
        self.assertEqual(calc_jk_xc.count("infos%control%pcm_enabled"), 1)

    def test_no_second_runtime_pcm_path_remains(self):
        """Prove the duplicate/double-count path cannot silently execute."""
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        # The retired reference-supplied-potential hook, argument, and energy
        # field must be gone everywhere in the SCF Fock builder module.
        self.assertNotIn("add_reference_pcm_reaction_field", text)
        self.assertNotIn("pcm_reaction_potential", text)
        self.assertNotIn("%epcm", text)
        self.assertNotIn(":: epcm", text)
        # calc_fock takes no PCM argument; PCM is driven only inside calc_jk_xc.
        calc_fock = text.split("subroutine calc_fock", 1)[1].split("end subroutine calc_fock", 1)[0]
        self.assertNotIn("pcm_reaction_potential_in", calc_fock)

    def test_solvent_py_keeps_only_validation_helpers(self):
        """solvent.py must not carry runtime handoff plumbing toward calc_fock."""
        solvent = _load_solvent("public_surface")
        # Preserved validation/cross-check helpers.
        for kept in (
            "reference_scf_total_density",
            "reference_scf_reaction_field_contract",
            "reference_scf_phi_cav_inputs",
            "reference_scf_pcm_energy_terms",
            "provisional_ddx_external_charges",
            "provisional_ddx_reaction_field_inputs",
        ):
            self.assertTrue(hasattr(solvent, kept), f"missing retained helper {kept}")
        # Removed runtime-handoff plumbing.
        for removed in (
            "reference_scf_pcm_runtime_payload",
            "reference_scf_pcm_reaction_potential_from_payload",
            "reference_scf_pcm_calc_fock_handoff",
            "reference_scf_pcm_calc_fock_handoff_from_molecule",
            "reference_scf_pcm_calc_fock_request",
            "reference_scf_pcm_calc_fock_request_from_scf_state",
            "reference_scf_pcm_calc_fock_call_site_bridge",
            "reference_scf_pcm_incremental_fock_audit",
            "reference_scf_reaction_fock_updates",
            "reference_scf_pcm_energy_handoff",
            "reference_scf_pcm_coupling_contract",
        ):
            self.assertFalse(hasattr(solvent, removed), f"runtime plumbing not removed: {removed}")

    def test_molecule_has_no_pcm_runtime_payload_persistence(self):
        text = (ROOT / "pyoqp" / "oqp" / "molecule" / "molecule.py").read_text(encoding="utf-8")
        self.assertNotIn("pcm_runtime_payload", text)
        self.assertNotIn("set_pcm_runtime_payload", text)
        self.assertNotIn("OQP::pcm_reaction_potential", text)


if __name__ == "__main__":
    unittest.main()
