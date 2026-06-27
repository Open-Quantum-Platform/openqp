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
        "conv": {"type": float, "default": "1.0e-6"},
    },
    "tdhf": {
        "type": {"type": _string, "default": "rpa"},
        "nstate": {"type": int, "default": "1"},
    },
    "optimize": {
        "lib": {"type": _string, "default": "oqp"},
        "istate": {"type": int, "default": "1"},
        "maxit": {"type": int, "default": "30"},
    },
    "oqp": {
        "coordsys": {"type": _string, "default": "tric"},
        "trust": {"type": float, "default": "0.2"},
    },
    "geometric": {
        "coordsys": {"type": _string, "default": "tric"},
        "trust": {"type": float, "default": "0.1"},
        "constraints_file": {"type": str, "default": ""},
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
        "oqp.utils.geometry",
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

        _load_module("oqp.utils.geometry", ROOT / "pyoqp/oqp/utils/geometry.py")
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


class TestOpenQPNativeAPI(unittest.TestCase):
    def test_builtin_geometry_resolves_common_names(self):
        openqp = load_openqp_module()

        geometry = openqp.get_geometry("water", source="builtin")

        self.assertIn("\nO ", geometry)
        self.assertEqual(len(geometry.strip().splitlines()), 3)

    def test_molecule_accepts_named_geometry(self):
        openqp = load_openqp_module()

        job = openqp.OpenQP(project="methane").molecule(geometry="ch4", basis="6-31g*")
        config = job.to_input_dict()

        self.assertEqual(len(config["input"]["system"].strip().splitlines()), 5)
        self.assertTrue(config["input"]["system"].startswith("\nC "))
        self.assertEqual(config["input"]["basis"], "6-31g*")

    def test_pubchem_sdf_parser_returns_openqp_geometry(self):
        openqp = load_openqp_module()
        sdf = """water
  OpenQP

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0
    0.7586    0.0000    0.5043 H   0  0  0  0  0  0
   -0.7586    0.0000    0.5043 H   0  0  0  0  0  0
M  END
$$$$
"""

        geometry = openqp.geometry_from_sdf(sdf, "water")

        self.assertEqual(
            geometry,
            "\nO 0 0 0\nH 0.7586 0 0.5043\nH -0.7586 0 0.5043",
        )

    def test_molecule_and_hf_helpers_build_openqp_input(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2", usempi=False)
            .molecule([("H", (0, 0, 0)), ("H", (0, 0, 1.4))], basis="6-31g*", charge=0)
            .hf()
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["system"], "\nH 0 0 0\nH 0 0 1.4")
        self.assertEqual(config["input"]["basis"], "6-31g*")
        self.assertEqual(config["input"]["charge"], "0")
        self.assertEqual(config["input"]["method"], "hf")
        self.assertEqual(config["input"]["runtype"], "energy")
        self.assertEqual(config["scf"]["type"], "rhf")

    def test_dft_helper_sets_functional_separately_from_hf(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2o_pbe")
            .molecule(geometry="water", basis="6-31g*")
            .dft("pbe", reference="rhf", runtype="grad", conv=1.0e-7)
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "hf")
        self.assertEqual(config["input"]["functional"], "pbe")
        self.assertEqual(config["input"]["runtype"], "grad")
        self.assertEqual(config["scf"]["type"], "rhf")
        self.assertEqual(config["scf"]["conv"], "1e-07")

    def test_mrsf_helper_uses_openqp_defaults(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(project="h2_mrsf").molecule("H 0 0 0; H 0 0 0.74").mrsf(nstate=4)

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["input"]["runtype"], "energy")
        self.assertEqual(config["scf"]["type"], "rohf")
        self.assertEqual(config["scf"]["multiplicity"], "3")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["nstate"], "4")

    def test_mrsf_helper_accepts_inline_functional(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2o_mrsf")
            .molecule(geometry="water", basis="6-31g*")
            .mrsf(nstate=5, functional="bhhlyp")
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["input"]["functional"], "bhhlyp")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["nstate"], "5")

    def test_mrsf_helper_preserves_existing_functional(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2o_mrsf")
            .molecule(geometry="water", basis="6-31g*")
            .input(functional="pbe0")
            .mrsf(nstate=4)
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["functional"], "pbe0")
        self.assertEqual(config["tdhf"]["nstate"], "4")

    def test_section_proxy_updates_openqp_keywords(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP().molecule("H 0 0 0; H 0 0 0.74")

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

    def test_optimize_helper_routes_native_backend_options(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(project="h2o_opt").molecule(geometry="water", basis="6-31g*")

        job.optimize(lib="oqp", istate=0, maxit=10, coordsys="dlc", trust=0.25)

        config = job.to_input_dict()
        self.assertEqual(config["optimize"]["lib"], "oqp")
        self.assertEqual(config["optimize"]["istate"], "0")
        self.assertEqual(config["optimize"]["maxit"], "10")
        self.assertEqual(config["oqp"]["coordsys"], "dlc")
        self.assertEqual(config["oqp"]["trust"], "0.25")
        self.assertEqual(job.optimize.coordsys, "dlc")

    def test_optimize_helper_routes_geometric_backend_options(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(project="h2o_geometric").molecule(geometry="water", basis="6-31g*")

        job.optimize(
            lib="geometric",
            maxit=8,
            coordsys="tric",
            trust=0.12,
            constraints_file="bond.constraints",
        )

        config = job.to_input_dict()
        self.assertEqual(config["optimize"]["lib"], "geometric")
        self.assertEqual(config["optimize"]["maxit"], "8")
        self.assertEqual(config["geometric"]["coordsys"], "tric")
        self.assertEqual(config["geometric"]["trust"], "0.12")
        self.assertEqual(config["geometric"]["constraints_file"], "bond.constraints")
        self.assertEqual(job.optimize.constraints_file, "bond.constraints")

    def test_run_builds_runner_lazily_and_returns_molecule(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(project="h_atom", usempi=False).molecule("H 0 0 0", basis="sto-3g")

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
