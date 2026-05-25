"""Checks for the planned ddX-to-SCF integration seam."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class DDXSCFIntegrationSeamTests(unittest.TestCase):
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

    def test_scf_fock_builder_is_identified_for_solvent_hook(self):
        text = (ROOT / "source" / "scf_addons.F90").read_text(encoding="utf-8")
        self.assertIn("subroutine calc_jk_xc", text.lower())
        self.assertIn("f(:,ii) =  f(:,ii) + hcore", text)
        self.assertIn("E%ehf1 = E%ehf1 + traceprod_sym_packed", text)


if __name__ == "__main__":
    unittest.main()
