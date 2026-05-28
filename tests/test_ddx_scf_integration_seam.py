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

    def test_provisional_ddx_reaction_field_inputs_label_sign_scale_and_runtime_scope(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_provisional_scope", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        inputs = solvent.provisional_ddx_reaction_field_inputs(
            [0.2],
            [0.0, 0.1, 0.2],
            allow_provisional=True,
        )

        self.assertTrue(inputs["provisional_sign_scale"])
        self.assertEqual(inputs["sign_scale_convention"], "chg = -0.5 * q_cav")
        self.assertEqual(inputs["validation_status"], "requires PySCF/ddX/reference cross-check before runtime use")
        self.assertFalse(inputs["runtime_pcm_enabled"])

    def test_provisional_ddx_reaction_field_inputs_preserve_first_scope_metadata(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_provisional_first_scope", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        inputs = solvent.provisional_ddx_reaction_field_inputs(
            [0.2],
            [0.0, 0.1, 0.2],
            allow_provisional=True,
        )

        self.assertEqual(inputs["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(inputs["reference_target"], "RHF/ROHF reference density")
        self.assertEqual(inputs["response_solvent_coupling"], "not enabled")
        self.assertEqual(inputs["gradient_support"], "not enabled")
        self.assertEqual(inputs["backend_validation_status"], "pending PySCF/ddX/reference cross-check")

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
        self.assertFalse(contract["runtime_pcm_enabled"])
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

    def test_reference_scf_pcm_energy_terms_preserve_energy_only_scope_metadata(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_energy_scope", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        terms = solvent.reference_scf_pcm_energy_terms(
            [[1.0, 0.5, 0.25]],
            [0.1, 0.2, 0.2],
        )

        self.assertEqual(terms["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(terms["reference_target"], "RHF/ROHF reference density")
        self.assertEqual(terms["response_solvent_coupling"], "not enabled")
        self.assertEqual(terms["gradient_support"], "not enabled")
        self.assertFalse(terms["runtime_pcm_enabled"])
        self.assertEqual(terms["backend_validation_status"], "pending PySCF/ddX/reference cross-check")

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
        self.assertEqual(updates["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(updates["reference_target"], "RHF/ROHF reference density")
        self.assertEqual(updates["response_solvent_coupling"], "not enabled")
        self.assertEqual(updates["gradient_support"], "not enabled")
        self.assertFalse(updates["runtime_pcm_enabled"])
        self.assertEqual(updates["backend_validation_status"], "pending PySCF/ddX/reference cross-check")

    def test_reference_scf_pcm_energy_handoff_packages_phi_fock_and_energy_terms(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_energy_handoff", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        handoff = solvent.reference_scf_pcm_energy_handoff(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
            [0.1, 0.2, 0.3],
        )

        self.assertEqual(handoff["nbf"], 2)
        self.assertEqual(handoff["nfocks"], 2)
        self.assertEqual(handoff["ncav"], 2)
        self.assertEqual(handoff["phi_cav_inputs"]["density_packed"], [1.75, 0.0, 0.25])
        self.assertEqual(handoff["fock_updates"], [[0.1, 0.2, 0.3], [0.1, 0.2, 0.3]])
        self.assertAlmostEqual(handoff["energy_terms"]["candidate_polarization_energy"], 0.125)
        self.assertEqual(handoff["handoff_stage"], "reference_scf_pcm_energy_prototype")
        self.assertEqual(handoff["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(handoff["response_solvent_coupling"], "not enabled")
        self.assertEqual(handoff["gradient_support"], "not enabled")
        self.assertFalse(handoff["runtime_pcm_enabled"])
        self.assertEqual(handoff["backend_validation_status"], "pending PySCF/ddX/reference cross-check")
        self.assertNotIn("density_blocks", handoff)
        self.assertNotIn("state_density", handoff)

    def test_reference_scf_pcm_runtime_payload_records_validated_potential_and_epcm(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_runtime_payload", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        payload = solvent.reference_scf_pcm_runtime_payload(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.1, 0.2, 0.3],
        )

        self.assertEqual(payload["OQP::pcm_reaction_potential"], [0.1, 0.2, 0.3])
        self.assertAlmostEqual(payload["OQP::pcm_epcm"], 0.125)
        self.assertEqual(payload["pcm_runtime_payload_version"], 1)
        self.assertEqual(payload["nbf"], 2)
        self.assertEqual(payload["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(payload["reference_target"], "RHF/ROHF reference density")
        self.assertEqual(payload["response_solvent_coupling"], "not enabled")
        self.assertEqual(payload["gradient_support"], "not enabled")
        self.assertFalse(payload["runtime_pcm_enabled"])
        self.assertEqual(payload["backend_validation_status"], "pending PySCF/ddX/reference cross-check")
        self.assertNotIn("density_blocks", payload)
        self.assertNotIn("state_density", payload)

    def test_molecule_round_trips_pcm_runtime_payload_outside_fortran_tags(self):
        text = (ROOT / "pyoqp" / "oqp" / "molecule" / "molecule.py").read_text(encoding="utf-8")
        self.assertIn("self.pcm_runtime_payload = {}", text)
        self.assertIn("def set_pcm_runtime_payload", text)
        self.assertIn("def get_pcm_runtime_payload", text)
        self.assertIn("def _restore_pcm_runtime_payload", text)
        self.assertIn("OQP::pcm_reaction_potential", text)
        self.assertIn("OQP::pcm_epcm", text)
        self.assertIn("pcm_runtime_payload_version", text)
        self.assertIn("data.update(self.get_pcm_runtime_payload())", text)
        self.assertIn("self._restore_pcm_runtime_payload(data)", text)
        self.assertNotIn("'OQP::pcm_reaction_potential',\n            'OQP::DM_A'", text)

    def test_reference_scf_pcm_reaction_potential_from_payload_requires_reviewed_scope(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_payload_consumer", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        payload = solvent.reference_scf_pcm_runtime_payload(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.1, 0.2, 0.3],
        )

        reviewed = solvent.reference_scf_pcm_reaction_potential_from_payload(payload)

        self.assertEqual(reviewed["reaction_potential"], [0.1, 0.2, 0.3])
        self.assertEqual(reviewed["nbf"], 2)
        self.assertAlmostEqual(reviewed["candidate_polarization_energy"], 0.125)
        self.assertEqual(reviewed["pcm_runtime_payload_version"], 1)
        self.assertEqual(reviewed["pcm_scope"], "reference_scf_energy_only")
        self.assertFalse(reviewed["runtime_pcm_enabled"])
        self.assertEqual(reviewed["handoff_target"], "calc_fock pcm_reaction_potential_in")

        with self.assertRaisesRegex(ValueError, "reference_scf_energy_only"):
            solvent.reference_scf_pcm_reaction_potential_from_payload({**payload, "pcm_scope": "state_specific"})
        with self.assertRaisesRegex(ValueError, "runtime_pcm_enabled must be False"):
            solvent.reference_scf_pcm_reaction_potential_from_payload({**payload, "runtime_pcm_enabled": True})
        with self.assertRaisesRegex(ValueError, "RHF/ROHF reference density"):
            solvent.reference_scf_pcm_reaction_potential_from_payload(
                {**payload, "reference_target": "MRSF state density"}
            )
        missing_potential = {
            "pcm_scope": "reference_scf_energy_only",
            "reference_target": "RHF/ROHF reference density",
            "runtime_pcm_enabled": False,
            "response_solvent_coupling": "not enabled",
            "gradient_support": "not enabled",
            "pcm_runtime_payload_version": 1,
            "OQP::pcm_epcm": 0.125,
        }
        with self.assertRaisesRegex(ValueError, "OQP::pcm_reaction_potential"):
            solvent.reference_scf_pcm_reaction_potential_from_payload(missing_potential)
        with self.assertRaisesRegex(ValueError, "nbf"):
            solvent.reference_scf_pcm_reaction_potential_from_payload({**payload, "nbf": 3})
        with self.assertRaisesRegex(ValueError, "nbf must be an integer"):
            solvent.reference_scf_pcm_reaction_potential_from_payload({**payload, "nbf": 2.5})

    def test_reference_scf_pcm_reaction_potential_from_payload_rejects_density_leakage(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_payload_leakage", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        payload = solvent.reference_scf_pcm_runtime_payload(
            [[1.0, 0.25, 0.75]],
            [0.1, 0.2, 0.3],
        )
        with self.assertRaisesRegex(ValueError, "density_blocks"):
            solvent.reference_scf_pcm_reaction_potential_from_payload({**payload, "density_blocks": [[1.0]]})
        with self.assertRaisesRegex(ValueError, "state_density"):
            solvent.reference_scf_pcm_reaction_potential_from_payload({**payload, "state_density": [1.0]})

    def test_reference_scf_pcm_calc_fock_handoff_exposes_only_guarded_keyword_argument(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_calc_fock_handoff", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        payload = solvent.reference_scf_pcm_runtime_payload(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.1, 0.2, 0.3],
        )

        handoff = solvent.reference_scf_pcm_calc_fock_handoff(payload)

        self.assertEqual(handoff["calc_fock_kwargs"], {"pcm_reaction_potential_in": [0.1, 0.2, 0.3]})
        self.assertEqual(handoff["nbf"], 2)
        self.assertEqual(handoff["packed_ao_length"], 3)
        self.assertAlmostEqual(handoff["candidate_polarization_energy"], 0.125)
        self.assertEqual(handoff["handoff_target"], "calc_fock pcm_reaction_potential_in")
        self.assertEqual(handoff["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(handoff["response_solvent_coupling"], "not enabled")
        self.assertEqual(handoff["gradient_support"], "not enabled")
        self.assertFalse(handoff["runtime_pcm_enabled"])
        self.assertNotIn("density_blocks", handoff)
        self.assertNotIn("state_density", handoff)
        with self.assertRaisesRegex(ValueError, "runtime_pcm_enabled must be False"):
            solvent.reference_scf_pcm_calc_fock_handoff({**payload, "runtime_pcm_enabled": True})

    def test_reference_scf_pcm_calc_fock_handoff_from_molecule_requires_explicit_payload(self):
        import importlib.util
        import types

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_molecule_calc_fock_handoff", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        empty_mol = types.SimpleNamespace(get_pcm_runtime_payload=lambda: {})
        disabled = solvent.reference_scf_pcm_calc_fock_handoff_from_molecule(empty_mol)
        self.assertEqual(disabled["calc_fock_kwargs"], {})
        self.assertFalse(disabled["payload_present"])
        self.assertEqual(disabled["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(disabled["reference_target"], "RHF/ROHF reference density")
        self.assertFalse(disabled["runtime_pcm_enabled"])
        self.assertEqual(
            disabled["backend_validation_status"],
            "pending PySCF/ddX/reference cross-check",
        )
        self.assertEqual(disabled["handoff_target"], "calc_fock pcm_reaction_potential_in")

        payload = solvent.reference_scf_pcm_runtime_payload(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.1, 0.2, 0.3],
        )
        payload_mol = types.SimpleNamespace(get_pcm_runtime_payload=lambda: payload)
        handoff = solvent.reference_scf_pcm_calc_fock_handoff_from_molecule(payload_mol)
        self.assertEqual(handoff["calc_fock_kwargs"], {"pcm_reaction_potential_in": [0.1, 0.2, 0.3]})
        self.assertTrue(handoff["payload_present"])
        self.assertEqual(handoff["nbf"], 2)
        self.assertEqual(handoff["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(handoff["response_solvent_coupling"], "not enabled")
        self.assertFalse(handoff["runtime_pcm_enabled"])
        with self.assertRaisesRegex(ValueError, "runtime_pcm_enabled must be False"):
            solvent.reference_scf_pcm_calc_fock_handoff_from_molecule(
                types.SimpleNamespace(get_pcm_runtime_payload=lambda: {**payload, "runtime_pcm_enabled": True})
            )

    def test_reference_scf_pcm_calc_fock_request_blocks_incremental_fock(self):
        import importlib.util
        import types

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_calc_fock_request", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        payload = solvent.reference_scf_pcm_runtime_payload(
            [[1.0, 0.5, 0.25], [0.75, -0.5, 0.0]],
            [0.1, 0.2, 0.3],
        )
        request = solvent.reference_scf_pcm_calc_fock_request(
            types.SimpleNamespace(get_pcm_runtime_payload=lambda: payload),
            incremental_fock=False,
        )
        self.assertEqual(request["calc_fock_kwargs"], {"pcm_reaction_potential_in": [0.1, 0.2, 0.3]})
        self.assertEqual(request["call_mode"], "non_incremental_only")
        self.assertTrue(request["requires_non_incremental_fock"])
        self.assertFalse(request["incremental_fock_allowed"])
        self.assertTrue(request["payload_present"])
        self.assertFalse(request["runtime_pcm_enabled"])
        self.assertEqual(request["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(request["handoff_target"], "calc_fock pcm_reaction_potential_in")

        disabled = solvent.reference_scf_pcm_calc_fock_request(
            types.SimpleNamespace(get_pcm_runtime_payload=lambda: {}),
            incremental_fock=True,
        )
        self.assertEqual(disabled["calc_fock_kwargs"], {})
        self.assertFalse(disabled["payload_present"])
        self.assertEqual(disabled["call_mode"], "disabled_no_payload")

        with self.assertRaisesRegex(ValueError, "reference PCM incremental Fock is not validated"):
            solvent.reference_scf_pcm_calc_fock_request(
                types.SimpleNamespace(get_pcm_runtime_payload=lambda: payload),
                incremental_fock=True,
            )

    def test_reference_scf_pcm_calc_fock_request_requires_boolean_incremental_flag(self):
        import importlib.util
        import types

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_calc_fock_request_bool", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        payload = solvent.reference_scf_pcm_runtime_payload(
            [[1.0, 0.5, 0.25]],
            [0.1, 0.2, 0.3],
        )

        with self.assertRaisesRegex(ValueError, "incremental_fock must be boolean"):
            solvent.reference_scf_pcm_calc_fock_request(
                types.SimpleNamespace(get_pcm_runtime_payload=lambda: payload),
                incremental_fock="unknown",
            )

    def test_reference_scf_pcm_incremental_fock_audit_preserves_old_buffer_metadata(self):
        import importlib.util

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_incremental_audit", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        audit = solvent.reference_scf_pcm_incremental_fock_audit(dens_old=[0.0], f_old=None)
        self.assertTrue(audit["scf_state_incremental_fock"])
        self.assertTrue(audit["scf_state_dens_old_present"])
        self.assertFalse(audit["scf_state_f_old_present"])
        self.assertEqual(audit["incremental_trigger_fields"], ["dens_old"])
        self.assertEqual(
            audit["scf_state_guard"],
            "dens_old/f_old presence blocks reviewed reference PCM payloads",
        )

        empty = solvent.reference_scf_pcm_incremental_fock_audit(dens_old=None, f_old=None)
        self.assertFalse(empty["scf_state_incremental_fock"])
        self.assertEqual(empty["incremental_trigger_fields"], [])

    def test_reference_scf_pcm_calc_fock_request_from_scf_state_blocks_old_buffers(self):
        import importlib.util
        import types

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_calc_fock_request_state", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        payload = solvent.reference_scf_pcm_runtime_payload(
            [[1.0, 0.5, 0.25]],
            [0.1, 0.2, 0.3],
        )
        mol = types.SimpleNamespace(get_pcm_runtime_payload=lambda: payload)
        request = solvent.reference_scf_pcm_calc_fock_request_from_scf_state(
            mol,
            dens_old=None,
            f_old=None,
        )
        self.assertEqual(request["call_mode"], "non_incremental_only")
        self.assertEqual(request["scf_state_incremental_fock"], False)
        self.assertFalse(request["scf_state_dens_old_present"])
        self.assertFalse(request["scf_state_f_old_present"])
        self.assertEqual(request["incremental_trigger_fields"], [])
        self.assertEqual(request["calc_fock_kwargs"], {"pcm_reaction_potential_in": [0.1, 0.2, 0.3]})

        with self.assertRaisesRegex(
            ValueError,
            "reference PCM incremental Fock is not validated.*dens_old_present=true.*f_old_present=false",
        ):
            solvent.reference_scf_pcm_calc_fock_request_from_scf_state(mol, dens_old=[0.0], f_old=None)
        with self.assertRaisesRegex(
            ValueError,
            "reference PCM incremental Fock is not validated.*dens_old_present=false.*f_old_present=true",
        ):
            solvent.reference_scf_pcm_calc_fock_request_from_scf_state(mol, dens_old=None, f_old=[0.0])
        with self.assertRaisesRegex(
            ValueError,
            "reference PCM incremental Fock is not validated.*dens_old_present=true.*f_old_present=true.*incremental_trigger_fields=dens_old,f_old",
        ):
            solvent.reference_scf_pcm_calc_fock_request_from_scf_state(mol, dens_old=[0.0], f_old=[0.0])

        disabled = solvent.reference_scf_pcm_calc_fock_request_from_scf_state(
            types.SimpleNamespace(get_pcm_runtime_payload=lambda: {}),
            dens_old=[0.0],
            f_old=[0.0],
        )
        self.assertEqual(disabled["call_mode"], "disabled_no_payload")
        self.assertFalse(disabled["payload_present"])
        self.assertTrue(disabled["scf_state_incremental_fock"])
        self.assertTrue(disabled["scf_state_dens_old_present"])
        self.assertTrue(disabled["scf_state_f_old_present"])

    def test_reference_scf_pcm_calc_fock_call_site_bridge_forwards_only_reviewed_nonincremental_payload(self):
        import importlib.util
        import types

        module_path = ROOT / "pyoqp" / "oqp" / "library" / "solvent.py"
        spec = importlib.util.spec_from_file_location("solvent_under_test_calc_fock_bridge", module_path)
        if spec is None or spec.loader is None:
            self.fail(f"Unable to load {module_path}")
        solvent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(solvent)

        payload = solvent.reference_scf_pcm_runtime_payload(
            [[1.0, 0.5, 0.25]],
            [0.1, 0.2, 0.3],
        )
        bridge = solvent.reference_scf_pcm_calc_fock_call_site_bridge(
            types.SimpleNamespace(get_pcm_runtime_payload=lambda: payload),
            dens_old=None,
            f_old=None,
        )
        self.assertTrue(bridge["forward_pcm_reaction_potential"])
        self.assertEqual(bridge["calc_fock_kwargs"], {"pcm_reaction_potential_in": [0.1, 0.2, 0.3]})
        self.assertEqual(bridge["call_site_bridge"], "reference_scf_calc_fock")
        self.assertEqual(bridge["call_mode"], "non_incremental_only")
        self.assertFalse(bridge["runtime_pcm_enabled"])
        self.assertEqual(bridge["pcm_scope"], "reference_scf_energy_only")
        self.assertEqual(bridge["response_solvent_coupling"], "not enabled")
        self.assertEqual(bridge["gradient_support"], "not enabled")

        disabled = solvent.reference_scf_pcm_calc_fock_call_site_bridge(
            types.SimpleNamespace(get_pcm_runtime_payload=lambda: {}),
            dens_old=[0.0],
            f_old=[0.0],
        )
        self.assertFalse(disabled["forward_pcm_reaction_potential"])
        self.assertEqual(disabled["calc_fock_kwargs"], {})
        self.assertEqual(disabled["call_mode"], "disabled_no_payload")
        self.assertTrue(disabled["scf_state_incremental_fock"])

        with self.assertRaisesRegex(
            ValueError,
            "reference PCM incremental Fock is not validated.*incremental_trigger_fields=dens_old",
        ):
            solvent.reference_scf_pcm_calc_fock_call_site_bridge(
                types.SimpleNamespace(get_pcm_runtime_payload=lambda: payload),
                dens_old=[0.0],
                f_old=None,
            )

    def test_calc_fock_pcm_incremental_guard_checks_both_old_density_and_fock_state(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        body = text.split("subroutine calc_fock", 1)[1].split("end subroutine calc_fock", 1)[0]
        guard = body.split("call calc_jk_xc", 1)[0]
        self.assertIn("present(pcm_reaction_potential_in)", guard)
        self.assertIn("present(dens_old)", guard)
        self.assertIn("present(f_old)", guard)
        self.assertIn("reference PCM incremental Fock is not validated", guard)
        self.assertIn("dens_old/f_old old-buffer state", guard)
        self.assertIn("dens_old_present=true", guard)
        self.assertIn("f_old_present=true", guard)
        self.assertIn("dens_old_present=false", guard)
        self.assertIn("f_old_present=false", guard)
        self.assertIn("incremental_trigger_fields=dens_old,f_old", guard)
        self.assertIn("incremental_trigger_fields=dens_old", guard)
        self.assertIn("incremental_trigger_fields=f_old", guard)

    def test_calc_fock_validates_pcm_reaction_potential_shape_before_fock_build(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        body = text.split("subroutine calc_fock", 1)[1].split("end subroutine calc_fock", 1)[0]
        guard = body.split("call calc_jk_xc", 1)[0]
        self.assertIn("size(pcm_reaction_potential_in) /= nbf_tri", guard)
        self.assertIn("reference PCM reaction potential length must match packed AO dimension", guard)

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
        self.assertIn("nfocks < 1", helper)
        self.assertIn("size(f,2) < nfocks", helper)
        self.assertIn("do ii = 1, nfocks", helper)
        self.assertIn("f(:,ii) = f(:,ii) + reaction_potential", helper)
        self.assertNotIn("pcm%enabled", helper.lower())

    def test_calc_fock_has_guarded_reference_pcm_reaction_potential_argument(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        calc_fock = text.split("subroutine calc_fock", 1)[1].split("end subroutine calc_fock", 1)[0]
        self.assertIn("pcm_reaction_potential_in", calc_fock)
        normalized = " ".join(calc_fock.split())
        self.assertIn("real(dp), intent(in), optional", normalized)
        self.assertIn("if (present(pcm_reaction_potential_in)) then", calc_fock)
        self.assertIn("pcm_reaction_potential=pcm_reaction_potential_in", calc_fock)
        self.assertNotIn("pcm%enabled", calc_fock.lower())

    def test_calc_fock_rejects_reference_pcm_incremental_fock_shortcut(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        calc_fock = text.split("subroutine calc_fock", 1)[1].split("end subroutine calc_fock", 1)[0]
        self.assertIn("reference PCM incremental Fock is not validated", calc_fock)
        self.assertIn("present(pcm_reaction_potential_in)", calc_fock)
        self.assertIn("present(dens_old)", calc_fock)
        incremental_branch = calc_fock.split("if (present(dens_old)) then", 1)[1].split("else", 1)[0]
        self.assertNotIn("pcm_reaction_potential=pcm_reaction_potential_in", incremental_branch)

    def test_scf_energy_type_has_dedicated_pcm_bookkeeping_field(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        energy_type = text.split("type :: scf_energy_t", 1)[1].split("end type scf_energy_t", 1)[0]
        self.assertIn("epcm", energy_type)
        self.assertIn("PCM reaction-field energy", energy_type)
        printer = text.split("subroutine print_scf_energy", 1)[1].split("end subroutine print_scf_energy", 1)[0]
        self.assertIn("PCM reaction-field energy", printer)
        calc_jk_xc = text.split("subroutine calc_jk_xc", 1)[1].split("end subroutine calc_jk_xc", 1)[0]
        self.assertIn("E%epcm = 0.0_dp", calc_jk_xc)

    def test_calc_jk_xc_records_opt_in_pcm_energy_without_double_counting_total(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        calc_jk_xc = text.split("subroutine calc_jk_xc", 1)[1].split("end subroutine calc_jk_xc", 1)[0]
        self.assertIn("E%epcm = E%epcm + 0.5_dp * traceprod_sym_packed(d(:,ii), pcm_reaction_potential, nbf)", calc_jk_xc)
        self.assertIn("E%etot = E%ehf + E%nenergy", calc_jk_xc)
        self.assertNotIn("E%etot = E%ehf + E%nenergy + E%epcm", calc_jk_xc)


if __name__ == "__main__":
    unittest.main()
