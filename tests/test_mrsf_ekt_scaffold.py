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


if __name__ == "__main__":
    unittest.main()
