import importlib.util
import json
import sys
import types
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(relpath):
    return (ROOT / relpath).read_text()


def minimal_config(**overrides):
    config = {
        "input": {
            "method": "tdhf",
            "runtype": "ekt",
            "basis": "6-31g*",
            "system": "\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 1.0 0.0 0.0",
        },
        "scf": {"type": "rohf", "multiplicity": 3},
        "tdhf": {"type": "mrsf", "multiplicity": 1, "nstate": 1, "nvdav": 10},
        "ekt": {"ip": True, "ea": False},
        "guess": {},
        "properties": {"grad": "1"},
    }
    for section, values in overrides.items():
        config.setdefault(section, {}).update(values)
    return config


def load_input_checker():
    """Load the backend-free input checker with a narrow MPI stub."""
    oqp_mod = sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    utils_mod = sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_mod = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    setattr(mpi_mod, "MPIManager", MPIManager)
    sys.modules["oqp.utils.mpi_utils"] = mpi_mod
    setattr(oqp_mod, "utils", utils_mod)
    setattr(utils_mod, "mpi_utils", mpi_mod)

    path = ROOT / "pyoqp/oqp/utils/input_checker.py"
    spec = importlib.util.spec_from_file_location("openqp_pr161_input_checker", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


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

    def test_ekt_has_dedicated_runtype_and_input_group(self):
        checker = read("pyoqp/oqp/utils/input_checker.py")
        data = read("pyoqp/oqp/molecule/oqpdata.py")
        runner = read("pyoqp/oqp/pyoqp.py")
        single_point = read("pyoqp/oqp/library/single_point.py")

        self.assertIn('"ekt"', checker)
        self.assertIn("'ekt':", data)
        self.assertIn("'ip': {'type': bool", data)
        self.assertIn("'ea': {'type': bool", data)
        self.assertIn("'ekt': compute_energy", runner)
        self.assertIn("self.runtype == 'ekt'", single_point)
        self.assertIn("tdhf.type=mrsf", checker)
        self.assertIn("MRSF-TDDFT", checker)
        self.assertIn("IP, EA, or both", checker)

    def test_legacy_mrsf_ekt_tdhf_types_are_energy_only(self):
        checker = read("pyoqp/oqp/utils/input_checker.py")

        self.assertIn("EKT analysis must use [input] runtype=ekt", checker)
        self.assertIn("Legacy tdhf.type=mrsf_ekt_ip/mrsf_ekt_ea", checker)
        self.assertIn("runtype != \"energy\"", checker)

    def test_legacy_mrsf_ekt_tdhf_types_error_for_grad_runtype(self):
        checker = load_input_checker()
        config = minimal_config(
            input={"runtype": "grad"},
            tdhf={"type": "mrsf_ekt_ip", "multiplicity": 1, "nstate": 1, "nvdav": 10},
        )

        report = checker.check_input_values(config, raise_error=False, emit=False)
        messages = "\n".join(item.message for item in report.errors)

        self.assertIn("Legacy tdhf.type=mrsf_ekt_ip/mrsf_ekt_ea is energy-only", messages)
        self.assertIn("EKT analysis must use [input] runtype=ekt", messages)

    def test_dedicated_ekt_runtype_requires_mrsf_tddft_and_ip_or_ea(self):
        checker = load_input_checker()
        config = minimal_config(
            tdhf={"type": "rpa", "multiplicity": 1, "nstate": 1, "nvdav": 10},
            ekt={"ip": False, "ea": False},
        )

        report = checker.check_input_values(config, raise_error=False, emit=False)
        messages = "\n".join(item.message for item in report.errors)

        self.assertIn("EKT runtype only supports MRSF-TDDFT", messages)
        self.assertIn("EKT runtype must request IP, EA, or both", messages)

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
        self.assertIn("density_ip_mo = density_alpha_mo", source)
        self.assertIn("density_ea_mo = density_beta_mo", source)
        self.assertNotIn("density_mo = density_alpha_mo + density_beta_mo", source)
        self.assertIn("lagrangian_mo", source)
        self.assertNotIn("Fock/Lagrangian proxy", source)
        self.assertNotIn("call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, lag_pack)", source)

    def test_examples_cover_mrsf_ekt_ip_and_ea(self):
        cases = {
            "examples/other/h2o_rohf_mrsf_ekt_ip_6-31g_bhhlyp.inp": "ip=True",
            "examples/other/h2o_rohf_mrsf_ekt_ea_6-31g_bhhlyp.inp": "ea=True",
        }

        for relpath, td_type in cases.items():
            with self.subTest(relpath=relpath):
                path = ROOT / relpath
                ref = path.with_suffix(".json")
                self.assertTrue(path.exists(), f"missing example input: {relpath}")
                self.assertTrue(ref.exists(), f"missing example reference: {ref.relative_to(ROOT)}")
                text = path.read_text()
                ref_data = json.loads(ref.read_text())
                self.assertIn("method=tdhf", text)
                self.assertIn("runtype=ekt", text)
                self.assertIn("type=rohf", text)
                self.assertIn("type=mrsf", text)
                self.assertIn("[ekt]", text)
                self.assertIn(td_type, text)
                self.assertIn("mrsf_ekt", ref_data)
                self.assertIn("eigenvalues_hartree", ref_data["mrsf_ekt"])
                self.assertIn("pole_strengths", ref_data["mrsf_ekt"])

    def test_final_json_persists_mrsf_ekt_root_results(self):
        fortran = read("source/modules/tdhf_mrsf_ekt.F90")
        tags = read("source/tagarray_driver.F90")
        molecule = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("OQP_mrsf_ekt_eigenvalues", tags)
        self.assertIn("OQP::mrsf_ekt_eigenvalues", molecule)
        self.assertIn("OQP::mrsf_ekt_strengths", molecule)
        self.assertIn("OQP::mrsf_ekt_orbitals_mo", molecule)
        self.assertIn("mrsf_ekt", molecule)
        self.assertIn("ebe_ev", molecule)
        self.assertIn("call infos%dat%reserve_data(OQP_mrsf_ekt_eigenvalues", fortran)
        self.assertIn("call tagarray_get_data(infos%dat, OQP_mrsf_ekt_eigenvalues", fortran)
        self.assertIn("print_ekt_dyson_orbitals", fortran)
        self.assertIn("EKT DYSON ORBITALS, ENERGIES AND NORMS", fortran)
        self.assertIn("basis%bf_label", fortran)

    def test_run_tests_checks_structured_mrsf_ekt_values(self):
        molecule = read("pyoqp/oqp/molecule/molecule.py")

        self.assertIn("'mrsf_ekt_ip'", molecule)
        self.assertIn("'mrsf_ekt_ea'", molecule)
        self.assertIn("required_ref_keys", molecule)
        self.assertIn("'mrsf_ekt'", molecule)
        self.assertIn("missing reference {key", molecule)
        self.assertIn("OQP::mrsf_ekt_density_mo", molecule)
        self.assertIn("OQP::mrsf_ekt_eigenvalues", molecule)
        self.assertIn("def compare_data", molecule)
        self.assertIn("isinstance(data_1, dict)", molecule)
        self.assertIn("arr_1.shape != arr_2.shape", molecule)
        self.assertIn("np.max(np.abs(arr_1 - arr_2))", molecule)
        self.assertIn("key == 'orbitals_mo'", molecule)


if __name__ == "__main__":
    unittest.main()
