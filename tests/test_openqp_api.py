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


if __name__ == "__main__":
    unittest.main()
