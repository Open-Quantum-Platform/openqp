import importlib.util
import json
import sys
import types
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_external_module():
    spec = importlib.util.spec_from_file_location(
        "external_native_exporter_under_test",
        ROOT / "pyoqp/oqp/library/external.py",
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def install_external_stubs():
    oqp = types.ModuleType("oqp")
    oqp.library = types.SimpleNamespace(
        set_basis=lambda mol: None,
        project_basis=lambda mol: None,
    )
    sys.modules["oqp"] = oqp

    utils = types.ModuleType("oqp.utils")
    utils.__path__ = []
    sys.modules["oqp.utils"] = utils

    constants = types.ModuleType("oqp.utils.constants")
    constants.ANGSTROM_TO_BOHR = 1.8897259886
    sys.modules["oqp.utils.constants"] = constants

    periodic_table_spec = importlib.util.spec_from_file_location(
        "oqp.periodic_table",
        ROOT / "pyoqp/oqp/periodic_table/__init__.py",
    )
    assert periodic_table_spec is not None
    assert periodic_table_spec.loader is not None
    periodic_table = importlib.util.module_from_spec(periodic_table_spec)
    sys.modules["oqp.periodic_table"] = periodic_table
    periodic_table_spec.loader.exec_module(periodic_table)

    # The native exporter must not depend on mokit being importable.
    sys.modules.pop("mokit", None)
    sys.modules.pop("mokit.lib", None)


class FakePySCFMol:
    def energy_tot(self):
        return -40.5


class FakeMeanField:
    def __init__(self):
        self.mol = FakePySCFMol()
        self.mo_coeff = np.array([[1.0, 0.2], [0.0, 0.8]])
        self.mo_energy = np.array([-0.7, 0.1])
        self.mo_occ = np.array([2.0, 0.0])

    def make_rdm1(self):
        return np.array([[2.0, 0.4], [0.4, 0.1]])


class FakeOQPMol:
    def __init__(self):
        self.project_name = "fake"
        self.config = {
            "input": {"basis": "sto-3g", "library": ""},
            "scf": {"type": "rhf"},
        }

    def get_atoms(self):
        return np.array([6, 1])

    def get_system(self):
        return np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0])


class TestNativePySCFExporter(unittest.TestCase):
    def setUp(self):
        install_external_stubs()
        self.external = load_external_module()

    def test_exporter_writes_openqp_restart_json_without_mokit(self):
        output = ROOT / "tests" / "tmp_native_pyscf_guess.json"
        try:
            result = self.external.export_pyscf_guess_to_openqp_json(
                FakeMeanField(),
                FakeOQPMol(),
                output,
            )
            data = json.loads(output.read_text())
        finally:
            output.unlink(missing_ok=True)

        self.assertEqual(result, output)
        self.assertEqual(data["atoms"], [6, 1])
        self.assertEqual(data["coord"], [0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
        self.assertEqual(data["json"], {"scf_type": "rhf", "basis": "sto-3g", "library": ""})
        self.assertEqual(data["OQP::VEC_MO_A"], [[1.0, 0.2], [0.0, 0.8]])
        self.assertEqual(data["OQP::VEC_MO_B"], [[1.0, 0.2], [0.0, 0.8]])
        self.assertEqual(data["OQP::E_MO_A"], [-0.7, 0.1])
        self.assertEqual(data["OQP::E_MO_B"], [-0.7, 0.1])
        self.assertEqual(data["OQP::DM_A"], [1.0, 0.0, 0.0])
        self.assertEqual(data["OQP::DM_B"], [1.0, 0.0, 0.0])
        self.assertEqual(data["energy"], -40.5)

    def test_exporter_handles_cr2_basis_with_g_functions(self):
        try:
            from pyscf import gto, scf
        except ModuleNotFoundError:
            self.skipTest("PySCF is not installed")

        pyscf_mol = gto.Mole()
        pyscf_mol.atom = "Cr 0.000000000 0.000000000 -0.840000000; Cr 0.000000000 0.000000000 0.840000000"
        pyscf_mol.unit = "Angstrom"
        pyscf_mol.basis = "cc-pvtz"
        pyscf_mol.charge = 0
        pyscf_mol.spin = 0
        pyscf_mol.verbose = 0
        pyscf_mol.build(cart=True)
        self.assertEqual(max(pyscf_mol.bas_angular(i) for i in range(pyscf_mol.nbas)), 4)

        mf = scf.RHF(pyscf_mol)
        mf.max_cycle = 0
        mf.kernel(dm0=mf.get_init_guess(key="atom"))

        oqp_mol = FakeOQPMol()
        oqp_mol.config["input"]["basis"] = "cc-pvtz"
        oqp_mol.get_atoms = lambda: np.array([24, 24])
        oqp_mol.get_system = lambda: np.array([
            0.000000000, 0.000000000, -0.840000000,
            0.000000000, 0.000000000, 0.840000000,
        ])

        output = ROOT / "tests" / "tmp_native_pyscf_cr2_g_guess.json"
        try:
            self.external.export_pyscf_guess_to_openqp_json(mf, oqp_mol, output)
            data = json.loads(output.read_text())
        finally:
            output.unlink(missing_ok=True)

        nbf = pyscf_mol.nao_nr()
        self.assertEqual(nbf, 168)
        self.assertEqual(len(data["OQP::VEC_MO_A"]), nbf)
        self.assertEqual(len(data["OQP::VEC_MO_A"][0]), nbf)
        self.assertEqual(len(data["OQP::E_MO_A"]), nbf)
        self.assertEqual(len(data["OQP::DM_A"]), nbf * (nbf + 1) // 2)
        self.assertEqual(data["atoms"], [24, 24])
        self.assertEqual(data["json"]["basis"], "cc-pvtz")


if __name__ == "__main__":
    unittest.main()
