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
        "soc_2e": {"type": int, "default": "1"},
        "system": {"type": str, "default": ""},
        "system2": {"type": str, "default": ""},
        "library": {"type": str, "default": ""},
        "ispher": {"type": _string, "default": "auto"},
        "omp_threads": {"type": int, "default": "0"},
        "qmmm_flag": {"type": bool, "default": "False"},
    },
    "qmmm": {
        "forcefield_files": {"type": str, "default": ""},
        "pdb_file": {"type": str, "default": ""},
        "qm_atoms": {"type": str, "default": ""},
        "cutoff": {"type": _string, "default": "NoCutoff"},
        "embedding": {"type": _string, "default": "electrostatic"},
        "rigidwater": {"type": bool, "default": "False"},
        "frontier_scheme": {"type": _string, "default": "none"},
    },
    "md": {
        "nstep": {"type": int, "default": "100"},
        "dt": {"type": float, "default": "0.5"},
        "active": {"type": int, "default": "1"},
        "soc": {"type": bool, "default": "False"},
        "soc_basis": {"type": _string, "default": "adiabatic"},
        "init_state": {"type": _string, "default": ""},
        "thrshe": {"type": float, "default": "1.0e9"},
        "init_temp": {"type": float, "default": "300.0"},
        "seed": {"type": int, "default": "1"},
        "rng_stream": {"type": int, "default": "0"},
        "first_hop_step": {"type": int, "default": "2"},
        "nacme_check": {"type": _string, "default": "off"},
        "ba_gap_max": {"type": float, "default": "0.0734986443513"},
        "nacme_gate": {"type": _string, "default": "warn"},
        "nacme_gate_invariant_tol": {"type": float, "default": "1.0e-10"},
        "nacme_gate_abs_tol": {"type": float, "default": "1.0e-4"},
        "nacme_gate_rel_tol": {"type": float, "default": "1.0"},
        "nacme_gate_consecutive": {"type": int, "default": "3"},
        "nve_gate": {"type": _string, "default": "off"},
        "nve_gate_abs_tol": {"type": float, "default": "5.0e-3"},
        "nve_gate_step_tol": {"type": float, "default": "1.0e-3"},
        "nve_gate_transition_tol": {"type": float, "default": "1.0e-6"},
        "nve_gate_consecutive": {"type": int, "default": "3"},
        "trajectory_interval": {"type": int, "default": "1"},
        "restart_interval": {"type": int, "default": "1"},
        "trajectory_file": {"type": _string, "default": ""},
        "nacme_audit_file": {"type": _string, "default": ""},
        "restart_file": {"type": _string, "default": ""},
    },
    "odp": {
        "enabled": {"type": bool, "default": "False"},
        "cv": {"type": str, "default": ""},
        "scale": {"type": str, "default": ""},
        "reference_r": {"type": str, "default": ""},
        "reference_p": {"type": str, "default": ""},
        "center": {"type": float, "default": "0.0"},
        "k_parallel": {"type": float, "default": "0.0"},
        "k_perpendicular": {"type": float, "default": "0.0"},
        "window": {"type": int, "default": "0"},
    },
    "scf": {
        "type": {"type": _string, "default": "rhf"},
        "multiplicity": {"type": int, "default": "1"},
        "conv": {"type": float, "default": "1.0e-6"},
        "scal_rel": {"type": int, "default": "0"},
    },
    "tdhf": {
        "type": {"type": _string, "default": "rpa"},
        "nstate": {"type": int, "default": "1"},
        "multiplicity": {"type": int, "default": "1"},
    },
    "dftb": {
        "backend": {"type": _string, "default": "native"},
        "type": {"type": _string, "default": "auto"},
        "parameter_path": {"type": str, "default": ""},
        "nstate": {"type": int, "default": "3"},
    },
    "xtb": {
        "backend": {"type": _string, "default": "native"},
        "type": {"type": _string, "default": "auto"},
        "parameter_path": {"type": str, "default": ""},
        "model": {"type": _string, "default": "gfn1"},
        "lc_ground_state": {"type": bool, "default": "False"},
        "nstate": {"type": int, "default": "3"},
    },
    "mp2": {
        "variant": {"type": _string, "default": "mp2"},
        "same_spin_scale": {"type": float, "default": "1.0"},
        "opposite_spin_scale": {"type": float, "default": "1.0"},
    },
    "properties": {
        "grad": {"type": int, "default": "0"},
        "scf_prop": {"type": _string, "default": ""},
        "nmr_gauge": {"type": _string, "default": "cgo"},
    },
    "hess": {
        "type": {"type": _string, "default": "numerical"},
        "state": {"type": int, "default": "0"},
    },
    "nac": {
        "type": {"type": _string, "default": "nacme"},
        "states": {"type": str, "default": ""},
    },
    "ekt": {
        "ip": {"type": bool, "default": "False"},
        "ea": {"type": bool, "default": "False"},
    },
    "pcm": {
        "enabled": {"type": bool, "default": "False"},
        "backend": {"type": _string, "default": "ddx"},
        "mode": {"type": _string, "default": "reference_scf"},
        "model": {"type": _string, "default": "ddpcm"},
        "epsilon": {"type": float, "default": "78.3553"},
    },
    "optimize": {
        "lib": {"type": _string, "default": "oqp"},
        "istate": {"type": int, "default": "1"},
        "jstate": {"type": int, "default": "2"},
        "states": {"type": str, "default": ""},
        "meci_search": {"type": _string, "default": "penalty"},
        "pen_sigma": {"type": float, "default": "1.0"},
        "pen_alpha": {"type": float, "default": "0.0"},
        "pen_delta": {"type": float, "default": "0.025"},
        "pen_jump": {"type": str, "default": "10,25"},
        "energy_gap": {"type": float, "default": "1e-5"},
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
        "oqp.utils.tb_backends",
        "oqp.utils.state_labels",
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
        _load_module("oqp.utils.tb_backends", ROOT / "pyoqp/oqp/utils/tb_backends.py")
        _load_module("oqp.utils.state_labels", ROOT / "pyoqp/oqp/utils/state_labels.py")

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
    def test_odp_section_is_exposed_through_settings_api(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(project="odp_window")
        job.settings.odp(
            enabled=True,
            cv="distance(1,2);angle(1,2,3)",
            scale="0.5,1.0",
            reference_r="2.0,1.5",
            reference_p="3.0,2.1",
            center=0.4,
            k_parallel=0.08,
            k_perpendicular=0.01,
            window=7,
        )
        config = job.to_input_dict()
        self.assertEqual(config["odp"]["enabled"], "True")
        self.assertEqual(config["odp"]["cv"], "distance(1,2);angle(1,2,3)")
        self.assertEqual(config["odp"]["scale"], "0.5,1.0")
        self.assertEqual(config["odp"]["window"], "7")

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

    def test_molecule_accepts_second_geometry_and_multiplicity(self):
        openqp = load_openqp_module()

        system = "H 0 0 0; H 0 0 0.74"
        system2 = "H 0 0 0; H 0 0 0.80"
        job = openqp.OpenQP(project="h2_nac").molecule(
            system,
            system2,
            charge=0,
            multiplicity=3,
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["system"], "\nH 0 0 0\nH 0 0 0.74")
        self.assertEqual(config["input"]["system2"], "\nH 0 0 0\nH 0 0 0.80")
        self.assertEqual(config["input"]["charge"], "0")
        self.assertEqual(config["scf"]["multiplicity"], "3")

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
            openqp.OpenQP(project="h2")
            .control(usempi=False)
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

    def test_hf_helper_clears_prior_dft_functional(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="reuse_as_hf")
            .molecule(geometry="water", basis="6-31g*")
            .dft("pbe")
            .hf()
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "hf")
        self.assertEqual(config["input"]["functional"], "")
        self.assertEqual(config["scf"]["type"], "rhf")

    def test_mp2_helper_sets_energy_reference_and_variant(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2o_mp2")
            .molecule(geometry="water", basis="6-31g", charge=0, multiplicity=1)
            .mp2(reference="uhf", variant="scs-mp2", conv=1.0e-10)
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "mp2")
        self.assertEqual(config["input"]["functional"], "")
        self.assertEqual(config["input"]["runtype"], "energy")
        self.assertEqual(config["input"]["basis"], "6-31g")
        self.assertEqual(config["scf"]["type"], "uhf")
        self.assertEqual(config["scf"]["conv"], "1e-10")
        self.assertEqual(config["mp2"]["variant"], "scs-mp2")

    def test_mp2_helper_clears_prior_dft_and_sets_custom_scales(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="reuse_as_mp2")
            .molecule(geometry="water", basis="6-31g*")
            .dft("pbe", runtype="grad")
            .mp2(same_spin_scale=0.5, opposite_spin_scale=1.1)
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "mp2")
        self.assertEqual(config["input"]["functional"], "")
        self.assertEqual(config["input"]["runtype"], "energy")
        self.assertEqual(config["mp2"]["variant"], "custom")
        self.assertEqual(config["mp2"]["same_spin_scale"], "0.5")
        self.assertEqual(config["mp2"]["opposite_spin_scale"], "1.1")

    def test_mp2_helper_rejects_non_energy_and_functional(self):
        openqp = load_openqp_module()

        with self.assertRaisesRegex(ValueError, "runtype='energy'"):
            openqp.OpenQP(project="bad_mp2_runtype").mp2(runtype="grad")

        with self.assertRaisesRegex(ValueError, "do not pass functional"):
            openqp.OpenQP(project="bad_mp2_functional").theory("mp2", functional="pbe")

        with self.assertRaisesRegex(ValueError, "variant='custom'"):
            openqp.OpenQP(project="bad_mp2_scales").mp2(
                variant="sos-mp2",
                opposite_spin_scale=1.1,
            )

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

    def test_theory_helper_sets_basis_and_method(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2o_theory")
            .molecule(geometry="water", charge=0)
            .theory("mrsf-tddft", functional="bhhlyp", basis="6-31g*", nstate=6)
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["basis"], "6-31g*")
        self.assertEqual(config["input"]["functional"], "bhhlyp")
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["scf"]["type"], "rohf")
        self.assertEqual(config["scf"]["multiplicity"], "3")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["nstate"], "6")

    def test_theory_namespace_helpers_set_models(self):
        openqp = load_openqp_module()

        dft = (
            openqp.OpenQP(project="h2o_dft_namespace")
            .molecule(geometry="water", charge=0, multiplicity=1)
            .theory.dft(functional="pbe0", basis="6-31g*")
        )
        config = dft.to_input_dict()
        self.assertEqual(config["input"]["method"], "hf")
        self.assertEqual(config["input"]["functional"], "pbe0")
        self.assertEqual(config["input"]["basis"], "6-31g*")
        self.assertEqual(config["scf"]["type"], "rhf")

        mp2 = (
            openqp.OpenQP(project="h2o_mp2_namespace")
            .molecule(geometry="water", charge=0, multiplicity=1)
            .theory.mp2(basis="6-31g", reference="uhf", variant="sos-mp2")
        )
        config = mp2.to_input_dict()
        self.assertEqual(config["input"]["method"], "mp2")
        self.assertEqual(config["input"]["functional"], "")
        self.assertEqual(config["input"]["basis"], "6-31g")
        self.assertEqual(config["scf"]["type"], "uhf")
        self.assertEqual(config["mp2"]["variant"], "sos-mp2")

        mrsf = (
            openqp.OpenQP(project="h2o_mrsf_namespace")
            .molecule(geometry="water", charge=0)
            .theory.mrsf(functional="bhhlyp", basis="6-31g*", nstate=4)
        )
        config = mrsf.to_input_dict()
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["input"]["functional"], "bhhlyp")
        self.assertEqual(config["scf"]["type"], "rohf")
        self.assertEqual(config["scf"]["multiplicity"], "3")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["nstate"], "4")

        tddft = (
            openqp.OpenQP(project="h2o_tddft_namespace")
            .molecule(geometry="water", charge=0, multiplicity=1)
            .theory.tddft(functional="b3lyp5", basis="6-31g*", nstate=2)
        )
        config = tddft.to_input_dict()
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["input"]["functional"], "b3lyp5")
        self.assertEqual(config["tdhf"]["nstate"], "2")

        with self.assertRaisesRegex(ValueError, "DFT theory requires"):
            openqp.OpenQP(project="bad_dft_namespace").theory.dft()

    def test_theory_helper_sets_response_theories(self):
        openqp = load_openqp_module()

        tdhf = (
            openqp.OpenQP(project="h2o_tdhf")
            .molecule(geometry="water", charge=0)
            .theory("tdhf", basis="6-31g*", nstate=4)
        )
        config = tdhf.to_input_dict()
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["input"]["functional"], "")
        self.assertEqual(config["input"]["basis"], "6-31g*")
        self.assertEqual(config["scf"]["type"], "rhf")
        self.assertEqual(config["scf"]["multiplicity"], "1")
        self.assertEqual(config["tdhf"]["nstate"], "4")

        tddft = (
            openqp.OpenQP(project="h2o_tddft")
            .molecule(geometry="water", charge=0)
            .theory("tddft", functional="b3lyp5", basis="6-31g*", nstate=5)
        )
        config = tddft.to_input_dict()
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["input"]["functional"], "b3lyp5")
        self.assertEqual(config["tdhf"]["nstate"], "5")

        sf = (
            openqp.OpenQP(project="h2o_sf")
            .molecule(geometry="water", charge=0)
            .theory("sf-tddft", functional="bhhlyp", basis="6-31g*", nstate=3)
        )
        config = sf.to_input_dict()
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["input"]["functional"], "bhhlyp")
        self.assertEqual(config["scf"]["type"], "rohf")
        self.assertEqual(config["scf"]["multiplicity"], "3")
        self.assertEqual(config["tdhf"]["type"], "sf")
        self.assertEqual(config["tdhf"]["nstate"], "3")

        mp2 = (
            openqp.OpenQP(project="h2o_mp2")
            .molecule(geometry="water", charge=0)
            .theory("mp2", basis="6-31g", reference="rohf", multiplicity=1)
        )
        config = mp2.to_input_dict()
        self.assertEqual(config["input"]["method"], "mp2")
        self.assertEqual(config["input"]["basis"], "6-31g")
        self.assertEqual(config["scf"]["type"], "rohf")
        self.assertEqual(config["scf"]["multiplicity"], "1")

        with self.assertRaisesRegex(ValueError, "TDDFT theory requires"):
            openqp.OpenQP(project="bad_tddft").theory("tddft")
        with self.assertRaisesRegex(ValueError, "SF-TDDFT theory requires"):
            openqp.OpenQP(project="bad_sf").theory("sf-tddft")

    def test_control_sets_runtype_threads_and_optimizer_options(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2o_opt")
            .molecule(geometry="water", charge=0, multiplicity=1)
        )
        job.control(omp_threads=8, usempi=False)
        job.workflow.optimize(
            lib="oqp",
            maxit=12,
            coordsys="dlc",
            trust=0.25,
        ).theory("dft", functional="bhhlyp", basis="6-31g*")

        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "optimize")
        self.assertEqual(config["input"]["omp_threads"], "8")
        self.assertFalse(job.usempi)
        self.assertEqual(config["input"]["basis"], "6-31g*")
        self.assertEqual(config["input"]["functional"], "bhhlyp")
        self.assertEqual(config["optimize"]["lib"], "oqp")
        self.assertEqual(config["optimize"]["maxit"], "12")
        self.assertEqual(config["oqp"]["coordsys"], "dlc")
        self.assertEqual(config["oqp"]["trust"], "0.25")

    def test_control_meci_sets_crossing_runtype_and_options(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2o_meci")
            .molecule(geometry="water", charge=0)
            .theory("mrsf-tddft", functional="bhhlyp", basis="6-31g*", nstate=5)
        )

        job.workflow.meci(lib="oqp", istate=1, jstate=2)

        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "meci")
        self.assertEqual(config["optimize"]["lib"], "oqp")
        self.assertEqual(config["optimize"]["istate"], "1")
        self.assertEqual(config["optimize"]["jstate"], "2")

    def test_meci_public_baeka_aliases_map_to_optimizer_schema(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(project="baeka").molecule(geometry="water")

        job.workflow.meci(
            states=[1, 2, 3],
            algorithm="baeka",
            sigma=2.0,
            alpha=0.02,
            delta_beta=0.05,
            beta_schedule=[10, 25],
            gap=1.0e-4,
        )

        optimize = job.to_input_dict()["optimize"]
        self.assertEqual(optimize["meci_search"], "baeka")
        self.assertEqual(optimize["states"], "1,2,3")
        self.assertEqual(optimize["pen_sigma"], "2.0")
        self.assertEqual(optimize["pen_alpha"], "0.02")
        self.assertEqual(optimize["pen_delta"], "0.05")
        self.assertEqual(optimize["pen_jump"], "10,25")
        self.assertEqual(optimize["energy_gap"], "0.0001")

    def test_workflow_sublevels_set_runtype_and_sections(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2o_workflows")
            .molecule(geometry="water", charge=0)
            .theory("mrsf-tddft", functional="bhhlyp", basis="6-31g*", nstate=6)
        )

        job.workflow.energy()
        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "energy")

        job.workflow.gradient(state=3)
        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "grad")
        self.assertEqual(config["properties"]["grad"], "3")

        with self.assertRaisesRegex(ValueError, "either state"):
            job.workflow.gradient(state=1, grad=1)

        job.workflow.hessian(type="analytical", state=0)
        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "hess")
        self.assertEqual(config["hess"]["type"], "analytical")
        self.assertEqual(config["hess"]["state"], "0")

        job.workflow.nacme(states="1,2")
        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "nacme")
        self.assertEqual(config["nac"]["states"], "1,2")

        job.workflow.ekt(ip=True, ea=False)
        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "ekt")
        self.assertEqual(config["ekt"]["ip"], "True")
        self.assertEqual(config["ekt"]["ea"], "False")

        pcm_job = (
            openqp.OpenQP(project="h2o_pcm")
            .molecule(geometry="water", charge=0, multiplicity=1)
            .theory("hf", basis="6-31g*")
        )
        pcm_job.workflow.pcm(
            enabled=True,
            backend="ddx",
            mode="reference_scf",
            model="ddpcm",
            epsilon=78.3553,
        )
        config = pcm_job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "energy")
        self.assertEqual(config["pcm"]["enabled"], "True")
        self.assertEqual(config["pcm"]["epsilon"], "78.3553")

    def test_workflow_pcm_requires_reference_scf_theory(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="bad_pcm")
            .molecule(geometry="water", charge=0)
            .theory("mrsf-tddft", functional="bhhlyp", basis="6-31g*", nstate=3)
        )

        with self.assertRaisesRegex(ValueError, "HF/DFT reference-SCF"):
            job.workflow.pcm(enabled=True, backend="ddx")

    def test_workflow_pcm_blocks_unsupported_scope(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="bad_pcm_scope")
            .molecule(geometry="water", charge=0, multiplicity=2)
            .theory("hf", reference="uhf", basis="6-31g*")
        )

        with self.assertRaisesRegex(ValueError, "RHF/ROHF"):
            job.workflow.pcm(enabled=True, backend="ddx")

        job = (
            openqp.OpenQP(project="bad_pcm_backend")
            .molecule(geometry="water", charge=0, multiplicity=1)
            .theory("hf", basis="6-31g*")
        )
        with self.assertRaisesRegex(ValueError, "backend='ddx'"):
            job.workflow.pcm(enabled=True, backend="pcmsolver")
        with self.assertRaisesRegex(ValueError, "mode='reference_scf'"):
            job.workflow.pcm(enabled=True, backend="ddx", mode="post_state_correction")

    def test_workflow_nmr_requires_reference_scf_theory(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="bad_nmr")
            .molecule(geometry="water", charge=0)
            .theory("mrsf-tddft", functional="bhhlyp", basis="6-31g*", nstate=3)
        )

        with self.assertRaisesRegex(ValueError, "HF/DFT reference-SCF"):
            job.workflow.nmr(gauge="cgo")

    def test_workflow_nmr_sets_properties_and_blocks_cgo_open_shell(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2o_nmr")
            .molecule(geometry="water", charge=0, multiplicity=1)
            .theory("dft", functional="bhhlyp", basis="6-31g*")
        )
        job.workflow.nmr(gauge="cgo")

        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "energy")
        self.assertEqual(config["properties"]["scf_prop"], "nmr")
        self.assertEqual(config["properties"]["nmr_gauge"], "cgo")

        open_shell = (
            openqp.OpenQP(project="bad_open_shell_nmr")
            .molecule(geometry="water", charge=0, multiplicity=3)
            .theory("hf", reference="rohf", basis="6-31g*")
        )
        with self.assertRaisesRegex(ValueError, "CGO NMR"):
            open_shell.workflow.nmr(gauge="cgo")
        open_shell.workflow.nmr(gauge="giao")

    def test_workflow_ekt_requires_mrsf_and_channel(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="bad_ekt")
            .molecule(geometry="water", charge=0, multiplicity=1)
            .theory("dft", functional="bhhlyp", basis="6-31g*")
        )

        with self.assertRaisesRegex(ValueError, "MRSF-TDDFT"):
            job.workflow.ekt(ip=True)

        mrsf = (
            openqp.OpenQP(project="bad_ekt_channel")
            .molecule(geometry="water", charge=0)
            .theory("mrsf-tddft", functional="bhhlyp", basis="6-31g*", nstate=5)
        )
        with self.assertRaisesRegex(ValueError, "requires ip=True"):
            mrsf.workflow.ekt()

    def test_control_call_remains_compatible_for_explicit_runtype(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2o_opt")
            .molecule(geometry="water", charge=0, multiplicity=1)
            .control(
                runtype="optimize",
                omp_threads=8,
                lib="oqp",
                maxit=12,
                coordsys="dlc",
                trust=0.25,
            )
            .theory("dft", functional="bhhlyp", basis="6-31g*")
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "optimize")
        self.assertEqual(config["input"]["omp_threads"], "8")
        self.assertEqual(config["input"]["basis"], "6-31g*")
        self.assertEqual(config["input"]["functional"], "bhhlyp")
        self.assertEqual(config["optimize"]["lib"], "oqp")
        self.assertEqual(config["optimize"]["maxit"], "12")
        self.assertEqual(config["oqp"]["coordsys"], "dlc")
        self.assertEqual(config["oqp"]["trust"], "0.25")

    def test_control_rejects_optimizer_options_for_nonoptimizer_runtype(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(project="bad_control").molecule(geometry="water")

        with self.assertRaisesRegex(KeyError, "known workflow"):
            job.control(runtype="energy", maxit=10)

    def test_soc_helper_sets_soc_without_response_multiplicity(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2o_soc")
            .molecule(geometry="water", charge=0)
            .theory(
                "mrsf-tddft",
                functional="bhhlyp",
                basis="6-31G(2df,p)",
                nstate=12,
            )
        )
        job.workflow.soc(soc_2e=1)

        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "soc")
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["input"]["basis"], "6-31G(2df,p)")
        self.assertEqual(config["input"]["functional"], "bhhlyp")
        self.assertEqual(config["input"]["soc_2e"], "1")
        self.assertEqual(config["scf"]["type"], "rohf")
        self.assertEqual(config["scf"]["multiplicity"], "3")
        self.assertEqual(config["scf"]["scal_rel"], "2")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["nstate"], "12")
        self.assertEqual(config["tdhf"]["multiplicity"], "1")

        job.workflow.soc(soc_2e=0, scal_rel=1)
        config = job.to_input_dict()
        self.assertEqual(config["input"]["soc_2e"], "0")
        self.assertEqual(config["scf"]["scal_rel"], "1")

    def test_workflow_soc_requires_mrsf_theory(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="bad_soc")
            .molecule(geometry="water", charge=0, multiplicity=1)
            .theory("dft", functional="bhhlyp", basis="6-31g*")
        )

        with self.assertRaisesRegex(ValueError, "only with MRSF-TDDFT"):
            job.workflow.soc(soc_2e=1)

    def test_workflow_soc_rejects_theory_options(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="bad_soc_options")
            .molecule(geometry="water", charge=0)
            .theory("mrsf-tddft", functional="bhhlyp", basis="6-31g*", nstate=12)
        )

        with self.assertRaisesRegex(ValueError, "Move these options"):
            job.workflow.soc(functional="bhhlyp")

    def test_qmmm_enables_flag_and_section(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="qmmm_energy")
            .molecule("ala.pdb 9 10 17 18 19", basis="6-31g*")
            .qmmm(embedding="electrostatic")
        )
        job.workflow.energy()
        config = job.to_input_dict()
        self.assertEqual(config["input"]["qmmm_flag"], "True")
        self.assertEqual(config["input"]["runtype"], "energy")
        self.assertEqual(config["qmmm"]["embedding"], "electrostatic")
        self.assertEqual(config["qmmm"]["pdb_file"], "ala.pdb")
        self.assertEqual(config["qmmm"]["qm_atoms"], "9 10 17 18 19")

    def test_qmmm_frontier_scheme_sets_section_key(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="qmmm_rcd")
            .molecule("ala.pdb 8 9 16 17 18", basis="6-31g*")
            .qmmm(pdb_file="ala.pdb", forcefield="amber14-all.xml",
                  qm_atoms=[8, 9, 16, 17, 18], embedding="electrostatic",
                  frontier_scheme="rcd")
        )
        config = job.to_input_dict()
        self.assertEqual(config["qmmm"]["frontier_scheme"], "rcd")

    def test_qmmm_normalizes_lists_and_rejects_duplicate_forcefield(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(project="qmmm_md").molecule("m.pdb 0 1 2", basis="6-31g")
        job.qmmm(
            forcefield=["amber14-all.xml", "amber14/tip3p.xml"],
            qm_atoms=[0, 1, 2],
            cutoff="PME",
            rigidwater=True,
        )
        config = job.to_input_dict()
        self.assertEqual(
            config["qmmm"]["forcefield_files"], "amber14-all.xml,amber14/tip3p.xml"
        )
        self.assertEqual(config["qmmm"]["pdb_file"], "m.pdb")
        self.assertEqual(config["qmmm"]["qm_atoms"], "0 1 2")
        self.assertEqual(config["qmmm"]["cutoff"], "PME")
        with self.assertRaisesRegex(ValueError, "either forcefield or forcefield_files"):
            job.qmmm(forcefield="a.xml", forcefield_files="b.xml")

    def test_workflow_namd_builds_soc_qmmm_deck(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="socnamd_qmmm")
            .molecule("chromo.pdb 0-4", basis="6-31g*")
            .theory("mrsf-tddft", functional="bhhlyp", nstate=3)
            .qmmm(cutoff="PME")
        )
        job.workflow.namd(
            soc=True,
            soc_basis="mch",
            nstep=200,
            dt=0.5,
            init_state="S1",
            seed=20260803,
            rng_stream=9,
            first_hop_step=2,
            nacme_check="baeck_an",
            ba_gap_max=0.05,
            nacme_gate="error",
            nacme_gate_abs_tol=2.0e-4,
            nacme_gate_rel_tol=0.5,
            nacme_gate_consecutive=4,
        )
        config = job.to_input_dict()
        self.assertEqual(config["input"]["qmmm_flag"], "True")
        self.assertEqual(config["input"]["runtype"], "namd")
        self.assertEqual(config["qmmm"]["pdb_file"], "chromo.pdb")
        self.assertEqual(config["qmmm"]["qm_atoms"], "0-4")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["md"]["soc"], "True")
        self.assertEqual(config["md"]["soc_basis"], "mch")
        self.assertEqual(config["md"]["nstep"], "200")
        self.assertEqual(config["md"]["init_state"], "S1")
        self.assertEqual(config["md"]["seed"], "20260803")
        self.assertEqual(config["md"]["rng_stream"], "9")
        self.assertEqual(config["md"]["first_hop_step"], "2")
        self.assertEqual(config["md"]["nacme_check"], "baeck_an")
        self.assertEqual(config["md"]["ba_gap_max"], "0.05")
        self.assertEqual(config["md"]["nacme_gate"], "error")
        self.assertEqual(config["md"]["nacme_gate_abs_tol"], "0.0002")
        self.assertEqual(config["md"]["nacme_gate_rel_tol"], "0.5")
        self.assertEqual(config["md"]["nacme_gate_consecutive"], "4")

    def test_workflow_namd_requires_mrsf_theory(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="bad_namd")
            .molecule(geometry="water")
            .theory("dft", functional="bhhlyp", basis="6-31g*")
        )
        with self.assertRaisesRegex(ValueError, "only with MRSF-TDDFT"):
            job.workflow.namd(nstep=10)

    def test_soc_helper_rejects_response_multiplicity(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP(project="h2o_soc").molecule(geometry="water")

        with self.assertRaisesRegex(ValueError, "do not set tdhf.multiplicity"):
            job.soc(nstate=12, functional="bhhlyp", **{"multiplicity": 3})

        job.soc(nstate=12, functional="bhhlyp", basis="6-31G(2df,p)")
        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "soc")
        self.assertEqual(config["input"]["soc_2e"], "1")
        self.assertEqual(config["scf"]["scal_rel"], "2")

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

    def test_settings_proxy_updates_openqp_keywords(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP().molecule("H 0 0 0; H 0 0 0.74")

        job.settings.input(method="tdhf", basis="6-31g*")
        job.settings.scf(type="rohf", multiplicity=3)
        job.settings.tdhf(type="sf", nstate=4)
        job.settings.tdhf.nstate = 5

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "tdhf")
        self.assertEqual(config["input"]["basis"], "6-31g*")
        self.assertEqual(config["scf"]["type"], "rohf")
        self.assertEqual(config["scf"]["multiplicity"], "3")
        self.assertEqual(config["tdhf"]["type"], "sf")
        self.assertEqual(config["tdhf"]["nstate"], "5")
        self.assertEqual(job.settings.tdhf.nstate, 5)

    def test_settings_basis_sets_atom_wise_basis_assignments(self):
        openqp = load_openqp_module()
        job = openqp.OpenQP().molecule("Br 0 0 0; H 0 0 1.4")

        job.settings.basis(["LANL2DZ", "6-31g*"])
        config = job.to_input_dict()
        self.assertEqual(config["input"]["basis"], "LANL2DZ;6-31g*")

        tagged = openqp.OpenQP().molecule(
            "C 0 0 0 c1; H 0 0 1 h1; H 1 0 0 h1"
        )
        tagged.settings.basis(c1="cc-pvdz", h1="6-31g*")
        config = tagged.to_input_dict()
        self.assertEqual(config["input"]["basis"], "library")
        self.assertEqual(config["input"]["library"], "c1 cc-pvdz\nh1 6-31g*")

        with self.assertRaisesRegex(ValueError, "single global basis"):
            openqp.OpenQP().settings.basis("6-31g*")

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
        job = (
            openqp.OpenQP(project="h_atom")
            .control(usempi=False)
            .molecule("H 0 0 0", basis="sto-3g")
        )

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
            unit = "B"

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


    def test_dftb_helper_builds_mrsf_tddftb_input(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2_dftb")
            .molecule([("H", (0, 0, 0)), ("H", (0, 0, 1.4))], basis="sto-3g", charge=0)
            .dftb(runtype="grad", response_type="mrsf", nstate=3,
                  parameter_path="/tmp/minimal_hh.opdftb")
        )
        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "dftb")
        self.assertEqual(config["input"]["runtype"], "grad")
        self.assertEqual(config["dftb"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["dftb"]["parameter_path"], "/tmp/minimal_hh.opdftb")

    def test_dftb_tda_alias_canonicalizes_backend_type(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2_dftb_tda")
            .molecule([("H", (0, 0, 0)), ("H", (0, 0, 1.4))], basis="sto-3g")
            .dftb(response_type="tda", parameter_path="/tmp/minimal_hh.opdftb")
        )
        config = job.to_input_dict()
        # backend method name is canonicalized to tddftb; tdhf.type stays tda.
        self.assertEqual(config["dftb"]["type"], "tddftb")
        self.assertEqual(config["tdhf"]["type"], "tda")

    def test_dftb_mrsf_permits_soc_and_namd_workflows(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2_dftb_soc")
            .molecule([("H", (0, 0, 0)), ("H", (0, 0, 1.4))], basis="sto-3g")
            .dftb(response_type="mrsf", parameter_path="/tmp/minimal_hh.opdftb")
        )
        # The MRSF-TDDFTB workflow guard must accept SOC/NAMD (they are wired for
        # method=dftb), matching the input-file validation path.
        job._require_mrsf_theory_for("SOC")
        job._require_mrsf_theory_for("NAMD")

    def test_xtb_helper_builds_mrsf_tddftb_input(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2_xtb")
            .molecule([("H", (0, 0, 0)), ("H", (0, 0, 1.4))], basis="sto-3g", charge=0)
            .xtb(runtype="grad", response_type="mrsf", nstate=3,
                 parameter_path="/tmp/gfn1.opxtb")
        )
        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "xtb")
        self.assertEqual(config["input"]["runtype"], "grad")
        self.assertEqual(config["xtb"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["xtb"]["parameter_path"], "/tmp/gfn1.opxtb")

    def test_xtb_tda_alias_canonicalizes_backend_type(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2_xtb_tda")
            .molecule([("H", (0, 0, 0)), ("H", (0, 0, 1.4))], basis="sto-3g")
            .xtb(response_type="tda", parameter_path="/tmp/gfn1.opxtb")
        )
        config = job.to_input_dict()
        # backend method name is canonicalized to tddftb; tdhf.type stays tda.
        self.assertEqual(config["xtb"]["type"], "tddftb")
        self.assertEqual(config["tdhf"]["type"], "tda")

    def test_xtb_helper_drains_gfn1_model_keyword_into_section(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2_xtb_lc")
            .molecule([("H", (0, 0, 0)), ("H", (0, 0, 1.4))], basis="sto-3g")
            .xtb(response_type="mrsf", parameter_path="/tmp/gfn1.opxtb",
                 model="gfn1", lc_ground_state=True)
        )
        config = job.to_input_dict()
        # [xtb]-only schema keywords are routed into the [xtb] section rather
        # than leaking into the [tdhf] response block.
        self.assertEqual(config["xtb"]["model"], "gfn1")
        self.assertEqual(config["xtb"]["lc_ground_state"], "True")
        self.assertNotIn("model", config["tdhf"])
        self.assertNotIn("lc_ground_state", config["tdhf"])

    def test_xtb_mrsf_permits_soc_and_namd_workflows(self):
        openqp = load_openqp_module()

        job = (
            openqp.OpenQP(project="h2_xtb_soc")
            .molecule([("H", (0, 0, 0)), ("H", (0, 0, 1.4))], basis="sto-3g")
            .xtb(response_type="mrsf", parameter_path="/tmp/gfn1.opxtb")
        )
        # The MRSF-xTB workflow guard must accept SOC/NAMD (they are wired for
        # method=xtb through the shared TB dispatch), matching the input-file
        # validation path.
        job._require_mrsf_theory_for("SOC")
        job._require_mrsf_theory_for("NAMD")

if __name__ == "__main__":
    unittest.main()
