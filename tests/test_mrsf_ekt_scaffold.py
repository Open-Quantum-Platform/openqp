import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


class TestMRSFEKTScaffold(unittest.TestCase):
    def test_c_header_declares_mrsf_ekt_entry_points(self):
        header = read("include/oqp.h")

        self.assertIn("void tdhf_mrsf_ekt_ip(struct oqp_handle_t *inf);", header)
        self.assertIn("void tdhf_mrsf_ekt_ea(struct oqp_handle_t *inf);", header)

    def test_python_dispatch_accepts_mrsf_ekt_ip_and_ea(self):
        source = read("pyoqp/oqp/library/single_point.py")

        self.assertIn("'mrsf_ekt_ip': oqp.tdhf_mrsf_ekt_ip", source)
        self.assertIn("'mrsf_ekt_ea': oqp.tdhf_mrsf_ekt_ea", source)
        self.assertIn("'mrsf_ekt_ip'", source)
        self.assertIn("'mrsf_ekt_ea'", source)

    def test_input_validation_knows_mrsf_ekt_types(self):
        checker = read("pyoqp/oqp/utils/input_checker.py")
        data = read("pyoqp/oqp/molecule/oqpdata.py")

        self.assertIn('"mrsf_ekt_ip"', checker)
        self.assertIn('"mrsf_ekt_ea"', checker)
        self.assertIn("'mrsf_ekt_ip'", data)
        self.assertIn("'mrsf_ekt_ea'", data)

    def test_fortran_module_ports_gamess_ekt_equations(self):
        source = read("source/modules/tdhf_mrsf_ekt.F90")

        self.assertIn('bind(C, name="tdhf_mrsf_ekt_ip")', source)
        self.assertIn('bind(C, name="tdhf_mrsf_ekt_ea")', source)
        self.assertIn("EKT: W * x = P * x * lambda", source)
        self.assertIn("EKT-EA: (F - W) * x = (I - P) * x * lambda", source)
        self.assertIn("call solve_symmetric_generalized", source)
        self.assertIn("ekt_metric = density_mo", source)
        self.assertIn("ekt_operator = fock_mo - lagrangian_mo", source)

    def test_fortran_module_uses_relaxed_target_state_density_and_lagrangian(self):
        source = read("source/modules/tdhf_mrsf_ekt.F90")

        self.assertIn("use tdhf_mrsf_z_vector_mod, only: tdhf_mrsf_z_vector", source)
        self.assertIn("call tdhf_mrsf_z_vector(infos)", source)
        self.assertIn("OQP_td_p", source)
        self.assertIn("OQP_WAO", source)
        self.assertIn("density_alpha_mo", source)
        self.assertIn("density_beta_mo", source)
        self.assertIn("density_mo = density_alpha_mo + density_beta_mo", source)
        self.assertIn("lagrangian_mo", source)
        self.assertNotIn("Fock/Lagrangian proxy", source)
        self.assertNotIn("call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, lag_pack)", source)

    def test_examples_cover_mrsf_ekt_ip_and_ea(self):
        cases = {
            "examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp.inp": "type=mrsf_ekt_ip",
            "examples/other/h2o_rohf_mrsf_ekt_ea_6-31g_bhhlyp.inp": "type=mrsf_ekt_ea",
        }

        for relpath, td_type in cases.items():
            with self.subTest(relpath=relpath):
                path = ROOT / relpath
                ref = path.with_suffix(".json")
                self.assertTrue(path.exists(), f"missing example input: {relpath}")
                self.assertTrue(ref.exists(), f"missing example reference: {ref.relative_to(ROOT)}")
                text = path.read_text()
                self.assertIn("method=tdhf", text)
                self.assertIn("runtype=energy", text)
                self.assertIn("type=rohf", text)
                self.assertIn(td_type, text)


if __name__ == "__main__":
    unittest.main()
