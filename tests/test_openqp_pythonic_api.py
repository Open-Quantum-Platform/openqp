import importlib.util
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _string(value):
    return str(value).lower()


SCHEMA = {
    "input": {
        "charge": {"type": int, "default": "0"},
        "basis": {"type": _string, "default": "6-31g*"},
        "functional": {"type": _string, "default": ""},
        "method": {"type": _string, "default": "hf"},
        "runtype": {"type": _string, "default": "energy"},
        "system": {"type": str, "default": ""},
        "system2": {"type": str, "default": ""},
        "ispher": {"type": _string, "default": "auto"},
        "omp_threads": {"type": int, "default": "0"},
    },
    "scf": {
        "type": {"type": _string, "default": "rhf"},
        "multiplicity": {"type": int, "default": "1"},
    },
    "tdhf": {
        "type": {"type": _string, "default": "rpa"},
        "nstate": {"type": int, "default": "1"},
    },
}


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_openqp_module():
    stub_names = (
        "oqp",
        "oqp.molecule",
        "oqp.molecule.oqpdata",
        "oqp.pyoqp",
        "oqp.utils",
        "oqp.utils.constants",
        "oqp.utils.input_parser",
        "oqp.utils.kword_map",
        "openqp_under_test",
    )
    saved_modules = {name: sys.modules.get(name) for name in stub_names}

    try:
        oqp = types.ModuleType("oqp")
        oqp.__path__ = []
        sys.modules["oqp"] = oqp

        molecule = types.ModuleType("oqp.molecule")
        molecule.__path__ = []
        sys.modules["oqp.molecule"] = molecule

        oqpdata = types.ModuleType("oqp.molecule.oqpdata")
        setattr(oqpdata, "OQP_CONFIG_SCHEMA", SCHEMA)
        sys.modules["oqp.molecule.oqpdata"] = oqpdata

        utils = types.ModuleType("oqp.utils")
        utils.__path__ = []
        sys.modules["oqp.utils"] = utils

        constants = types.ModuleType("oqp.utils.constants")
        setattr(constants, "ANGSTROM_TO_BOHR", 0.529177210903)
        sys.modules["oqp.utils.constants"] = constants

        _load_module("oqp.utils.input_parser", ROOT / "pyoqp/oqp/utils/input_parser.py")
        _load_module("oqp.utils.kword_map", ROOT / "pyoqp/oqp/utils/kword_map.py")

        class FakeMol:
            def __init__(self):
                self.loaded_configs = []

            def load_config(self, config):
                self.loaded_configs.append(config)

        class FakeRunner:
            instances = []

            def __init__(
                self,
                project=None,
                input_file=None,
                log=None,
                input_dict=None,
                silent=0,
                usempi=True,
            ):
                self.project = project
                self.input_file = input_file
                self.log = log
                self.input_dict = input_dict
                self.silent = silent
                self.usempi = usempi
                self.ran = False
                self.mol = FakeMol()
                self.__class__.instances.append(self)

            def run(self):
                self.ran = True

        pyoqp = types.ModuleType("oqp.pyoqp")
        setattr(pyoqp, "Runner", FakeRunner)
        sys.modules["oqp.pyoqp"] = pyoqp

        module = _load_module("openqp_under_test", ROOT / "pyoqp/oqp/openqp.py")
        module.Runner.instances.clear()
        return module
    finally:
        for name, module in saved_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


class TestOpenQPPythonicAPI(unittest.TestCase):
    def test_constructor_accepts_pyscf_style_molecule_arguments(self):
        openqp = load_openqp_module()

        job = openqp.OpenQP(
            atom=[("H", (0, 0, 0)), ("H", (0, 0, 1.4))],
            basis="6-31g*",
            charge=0,
            spin=2,
            project="h2_triplet",
            usempi=False,
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["system"], "\nH 0 0 0\nH 0 0 1.4")
        self.assertEqual(config["input"]["basis"], "6-31g*")
        self.assertEqual(config["input"]["charge"], "0")
        self.assertEqual(config["scf"]["multiplicity"], "3")
        self.assertEqual(job.scf.multiplicity, 3)

    def test_section_proxy_updates_openqp_keywords(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(atom="H 0 0 0; H 0 0 0.74")

        job.input(method="tdhf", runtype="energy")
        job.scf(type="rohf", multiplicity=3)
        job.tdhf.type = "mrsf"
        job.tdhf.nstate = 4

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["scf"]["type"], "rohf")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["nstate"], "4")
        self.assertEqual(job.tdhf.nstate, 4)

    def test_run_builds_runner_lazily_and_returns_molecule(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(atom="H 0 0 0", basis="sto-3g", project="h_atom", usempi=False)

        mol = job.run(run_type="grad")

        runner = openqp.Runner.instances[-1]
        self.assertIs(mol, runner.mol)
        self.assertIs(job.runner, runner)
        self.assertTrue(runner.ran)
        self.assertEqual(runner.project, "h_atom")
        self.assertEqual(runner.log, "h_atom.log")
        self.assertEqual(runner.input_dict["input"]["runtype"], "grad")
        self.assertFalse(runner.usempi)

    def test_from_pyscf_maps_spin_and_bohr_coordinates(self):
        openqp = load_openqp_module()

        class PySCFMol:
            atom = [["H", (0, 0, 0)], ["H", (0, 0, 1.0)]]
            basis = "sto-3g"
            charge = 1
            spin = 1
            unit = "Bohr"

        job = openqp.OpenQP.from_pyscf(PySCFMol())
        config = job.to_input_dict()

        self.assertEqual(config["input"]["system"], "\nH 0 0 0\nH 0 0 0.529177210903")
        self.assertEqual(config["input"]["basis"], "sto-3g")
        self.assertEqual(config["input"]["charge"], "1")
        self.assertEqual(config["scf"]["multiplicity"], "2")

    def test_legacy_openqp_wrapper_still_constructs_runner_immediately(self):
        openqp = load_openqp_module()

        wrapper = openqp.OPENQP(
            {
                "input.system": "H 0 0 0; H 0 0 0.74",
                "input.basis": "6-31g*",
                "input.method": "hf",
                "input.runtype": "energy",
                "scf.type": "rhf",
            }
        )
        runner = openqp.Runner.instances[-1]

        self.assertEqual(runner.input_dict["input"]["system"], "\nH 0 0 0\nH 0 0 0.74")
        self.assertFalse(runner.ran)
        mol = wrapper.run()
        self.assertIs(mol, runner.mol)
        self.assertTrue(runner.ran)


if __name__ == "__main__":
    unittest.main()
