import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def install_single_point_stubs():
    oqp = types.ModuleType("oqp")
    oqp.__path__ = []
    setattr(oqp, "hf_energy", lambda mol: None)
    setattr(oqp, "tdhf_energy", lambda mol: None)
    setattr(oqp, "tdhf_sf_energy", lambda mol: None)
    setattr(oqp, "tdhf_mrsf_energy", lambda mol: None)
    setattr(oqp, "tdhf_umrsf_energy", lambda mol: None)
    setattr(oqp, "library", types.SimpleNamespace(
        set_basis=lambda mol: None,
        ints_1e=lambda mol: None,
        guess=lambda mol: None,
        project_basis=lambda mol: None,
    ))
    sys.modules["oqp"] = oqp

    molecule = types.ModuleType("oqp.molecule")
    setattr(molecule, "Molecule", type("Molecule", (), {}))
    sys.modules["oqp.molecule"] = molecule

    utils = types.ModuleType("oqp.utils")
    utils.__path__ = []
    sys.modules["oqp.utils"] = utils

    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")
    setattr(mpi_utils, "MPIManager", type("MPIManager", (), {}))
    setattr(mpi_utils, "MPIPool", type("MPIPool", (), {}))
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils

    matrix = types.ModuleType("oqp.utils.matrix")
    setattr(matrix, "DampingParam", type("DampingParam", (), {}))
    setattr(matrix, "DispersionModel", type("DispersionModel", (), {}))
    sys.modules["oqp.utils.matrix"] = matrix

    file_utils = types.ModuleType("oqp.utils.file_utils")
    setattr(file_utils, "dump_log", lambda *args, **kwargs: None)
    setattr(file_utils, "dump_data", lambda *args, **kwargs: None)
    setattr(file_utils, "write_config", lambda *args, **kwargs: None)
    setattr(file_utils, "write_xyz", lambda *args, **kwargs: None)
    sys.modules["oqp.utils.file_utils"] = file_utils

    library = types.ModuleType("oqp.library")
    library.__path__ = []
    sys.modules["oqp.library"] = library

    frequency = types.ModuleType("oqp.library.frequency")
    setattr(frequency, "normal_mode", lambda *args, **kwargs: None)
    setattr(frequency, "thermal_analysis", lambda *args, **kwargs: None)
    sys.modules["oqp.library.frequency"] = frequency


class FakeData(dict):
    def __init__(self):
        super().__init__()
        self.convergers = []
        self.sd_scf_flags = []

    def set_scf_converger_type(self, converger):
        self.convergers.append(converger)

    def set_sd_scf(self, flag):
        self.sd_scf_flags.append(flag)

    def set_trah_stability(self, flag):
        self.trah_stability = flag


class FakeMol:
    def __init__(self):
        self.config = {
            "scf": {"converger_type": "diis", "trh_stab": False},
        }
        self.data = FakeData()
        self.mol_energy = types.SimpleNamespace(energy=-1.0, SCF_converged=False)
        self.energies = []

    def write_molden(self, filename):
        raise AssertionError("save_molden is disabled in this test")


class TestSinglePointScfFallback(unittest.TestCase):
    def setUp(self):
        install_single_point_stubs()
        self.single_point = load_module(
            "single_point_under_test",
            "pyoqp/oqp/library/single_point.py",
        )

    def make_calculator(self):
        calc = self.single_point.SinglePoint.__new__(self.single_point.SinglePoint)
        calc.mol = FakeMol()
        calc.method = "hf"
        calc.init_scf = "no"
        calc.forced_attempt = 2
        calc.converger_type = "diis"
        calc.alternative_scf = "trah"
        calc.stability = False
        calc.save_molden = False
        calc.exception = False
        calc._prep_guess = lambda: None
        calc.swapmo = lambda: None
        calc.ixcore_shift = lambda: None
        calc.scf_calls = 0

        def scf():
            calc.scf_calls += 1
            calc.mol.mol_energy.energy = -1.0 - calc.scf_calls
            calc.mol.mol_energy.SCF_converged = calc.scf_calls == 2

        calc.scf = scf
        return calc

    def test_energy_restores_fallback_converger_after_successful_recovery(self):
        calc = self.make_calculator()

        energy = calc.energy(do_init_scf=False)

        self.assertEqual(energy, [-3.0])
        self.assertEqual(calc.mol.data.convergers, ["diis", "trah", "diis", "diis"])
        self.assertEqual(calc.mol.data.sd_scf_flags, [False])

    def test_energy_restores_primary_converger_after_failed_recovery(self):
        calc = self.make_calculator()

        def scf_never_converges():
            calc.scf_calls += 1
            calc.mol.mol_energy.energy = -1.0 - calc.scf_calls
            calc.mol.mol_energy.SCF_converged = False

        calc.scf = scf_never_converges

        with self.assertRaises(RuntimeError):
            calc.energy(do_init_scf=False)

        self.assertEqual(calc.mol.data.convergers, ["diis", "trah", "diis", "diis"])

    def test_reference_restores_primary_converger_for_matching_gradient(self):
        calc = self.make_calculator()

        energy = calc.reference(do_init_scf=False)

        self.assertEqual(energy, [-3.0])
        self.assertEqual(calc.mol.data.convergers, ["diis", "trah", "diis"])

    def test_stability_noop_restores_pre_trah_energy_metadata(self):
        calc = self.make_calculator()
        calc.stability = True

        def scf_stable_after_trah_check():
            calc.scf_calls += 1
            if calc.scf_calls == 1:
                calc.mol.mol_energy.energy = -2.0
                calc.mol.mol_energy.SCF_converged = True
            else:
                calc.mol.mol_energy.energy = -1.99999999
                calc.mol.mol_energy.SCF_converged = True

        calc.scf = scf_stable_after_trah_check

        energy = calc.reference(do_init_scf=False)

        self.assertEqual(energy, [-2.0])
        self.assertEqual(calc.mol.mol_energy.energy, -2.0)


if __name__ == "__main__":
    unittest.main()
