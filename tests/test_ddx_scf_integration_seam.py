"""Checks for the planned ddX-to-SCF integration seam."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class DDXSCFIntegrationSeamTests(unittest.TestCase):
    def test_provisional_ddx_q_cav_to_external_charges_is_opt_in(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        with self.assertRaisesRegex(ValueError, "provisional"):
            solvent.provisional_ddx_external_charges([0.2, -0.4])

        charges = solvent.provisional_ddx_external_charges([0.2, -0.4], allow_provisional=True)
        self.assertEqual(charges, [-0.1, 0.2])

    def test_provisional_ddx_reaction_field_inputs_validate_cavity_shape(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_shapes", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        with self.assertRaisesRegex(ValueError, "3 \* len\(q_cav\)"):
            solvent.provisional_ddx_reaction_field_inputs(
                [0.2, -0.4],
                [0.0, 0.0, 0.0, 1.0, 0.0],
                allow_provisional=True,
            )

        with self.assertRaisesRegex(ValueError, "q_cav must contain at least one"):
            solvent.provisional_ddx_reaction_field_inputs(
                [],
                [],
                allow_provisional=True,
            )

        inputs = solvent.provisional_ddx_reaction_field_inputs(
            [0.2, -0.4],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            allow_provisional=True,
        )
        self.assertEqual(inputs["charges"], [-0.1, 0.2])
        self.assertEqual(inputs["cavity_xyz"], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        self.assertEqual(inputs["ncav"], 2)

    def test_provisional_ddx_reaction_field_inputs_split_external_charge_arrays(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_external_arrays", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        inputs = solvent.provisional_ddx_reaction_field_inputs(
            [0.2, -0.4],
            [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
            allow_provisional=True,
        )

        self.assertEqual(inputs["x"], [0.0, 1.0])
        self.assertEqual(inputs["y"], [0.1, 1.1])
        self.assertEqual(inputs["z"], [0.2, 1.2])
        self.assertEqual(inputs["chg"], [-0.1, 0.2])

    def test_solvent_helpers_reject_nonfinite_numeric_inputs(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_finite_values", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        with self.assertRaisesRegex(ValueError, "finite"):
            solvent.reference_scf_phi_cav_inputs([[1.0, float("nan"), 0.25]], [0.0, 0.1, 0.2])
        with self.assertRaisesRegex(ValueError, "finite"):
            solvent.provisional_ddx_reaction_field_inputs(
                [0.2, float("inf")],
                [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
                allow_provisional=True,
            )

    def test_reference_scf_total_density_sums_spin_blocks_for_ddx_phi_cav(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_total_density", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

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
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_rf_contract", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

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
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_phi_cav", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        with self.assertRaisesRegex(ValueError, "3 \\* ncav"):
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

    def test_reference_scf_pcm_coupling_contract_combines_phi_and_reaction_inputs(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_pcm_contract", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        with self.assertRaisesRegex(ValueError, "same packed length"):
            solvent.reference_scf_pcm_coupling_contract(
                [[1.0, 0.5, 0.25]],
                [0.0, 0.1, 0.2],
                [0.1, 0.2],
            )

        contract = solvent.reference_scf_pcm_coupling_contract(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
            [0.1, 0.2, 0.3],
        )
        self.assertEqual(contract["nbf"], 2)
        self.assertEqual(contract["nfocks"], 2)
        self.assertEqual(contract["ncav"], 2)
        self.assertEqual(contract["density_packed"], [1.75, 0.0, 0.25])
        self.assertEqual(contract["reaction_potential"], [0.1, 0.2, 0.3])
        self.assertEqual(contract["x"], [0.0, 1.0])
        self.assertEqual(contract["y"], [0.1, 1.1])
        self.assertEqual(contract["z"], [0.2, 1.2])
        self.assertEqual(contract["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(contract["reference_target"], "RHF/ROHF reference density")
        self.assertEqual(contract["response_solvent_coupling"], "not enabled")
        self.assertEqual(contract["gradient_support"], "not enabled")
        self.assertNotIn("density_blocks", contract)
        self.assertNotIn("state_density", contract)

    def test_reference_scf_pcm_energy_terms_use_validated_packed_density_and_potential(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_energy_terms", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

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

    def test_reference_scf_reaction_fock_updates_replicate_reaction_potential_per_fock(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_fock_updates", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        updates = solvent.reference_scf_reaction_fock_updates(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.1, 0.2, 0.3],
        )

        self.assertEqual(updates["nbf"], 2)
        self.assertEqual(updates["nfocks"], 2)
        self.assertEqual(updates["density_packed"], [1.75, 0.0, 0.25])
        self.assertEqual(updates["reaction_potential"], [0.1, 0.2, 0.3])
        self.assertEqual(updates["fock_updates"], [[0.1, 0.2, 0.3], [0.1, 0.2, 0.3]])
        self.assertAlmostEqual(updates["candidate_polarization_energy"], 0.125)
        self.assertEqual(updates["application_scope"], "add reaction_potential to each reference SCF Fock block")

    def test_unweighted_electrostatic_potential_is_public(self):
        text = (ROOT / "source" / "integrals" / "int1.F90").read_text(encoding="utf-8")
        self.assertIn("public electrostatic_potential_unweighted", text)
        self.assertIn("subroutine electrostatic_potential_unweighted", text.lower())
        wrapper = text.split("subroutine electrostatic_potential_unweighted", 1)[1].split(
            "end subroutine electrostatic_potential_unweighted", 1
        )[0]
        self.assertIn("call int1_el_pot", wrapper)
        self.assertNotIn("pot = pot*wt", wrapper)
        self.assertIn("call bas_denorm_matrix", wrapper)

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

    def test_calc_jk_xc_has_opt_in_reference_pcm_reaction_field_handoff(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        self.assertIn("call add_reference_pcm_reaction_field(f, pcm_reaction_potential, nfocks)", text)
        helper_pos = text.index("call add_reference_pcm_reaction_field(f, pcm_reaction_potential, nfocks)")
        energy_pos = text.index("E%ehf1 = E%ehf1 + traceprod_sym_packed")
        self.assertLess(helper_pos, energy_pos)
        self.assertIn("pcm_reaction_potential", text)
        self.assertIn("real(dp),          intent(in),   optional :: pcm_reaction_potential(:)", text)
        self.assertIn("if (present(pcm_reaction_potential)) then", text)
        self.assertNotIn("pcm%enabled", text[text.index("subroutine calc_jk_xc") : text.index("end subroutine calc_jk_xc")].lower())

    def test_reference_pcm_reaction_field_helper_replicates_per_fock_block(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        self.assertIn("public :: add_reference_pcm_reaction_field", text)
        self.assertIn("subroutine add_reference_pcm_reaction_field", text.lower())
        helper = text.split("subroutine add_reference_pcm_reaction_field", 1)[1].split(
            "end subroutine add_reference_pcm_reaction_field", 1
        )[0]
        self.assertIn("size(f,1) /= size(reaction_potential)", helper)
        self.assertIn("size(f,2) < nfocks", helper)
        self.assertIn("do ii = 1, nfocks", helper)
        self.assertIn("f(:,ii) = f(:,ii) + reaction_potential", helper)
        self.assertNotIn("pcm%enabled", helper.lower())


if __name__ == "__main__":
    unittest.main()
