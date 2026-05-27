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


if __name__ == "__main__":
    unittest.main()
