"""Python-API coverage for the optional OpenQP-xTB backend helpers.

Reuses the stubbed-schema module loader from ``tests/test_openqp_api.py`` so the
compact ``job.xtb(...)`` / ``job.theory('xtb', ...)`` entry points and the
xTB-MRSF workflow guards are exercised without the compiled library. Mirrors
``tests/test_openqp_dftb_api.py`` for method=xtb / the [xtb] section.
"""

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load_api_helpers():
    """Load tests/test_openqp_api.py under a unique name and return it.

    A distinct module name avoids clashing with pytest's own import of
    ``test_openqp_api`` during collection.
    """
    path = ROOT / "tests/test_openqp_api.py"
    spec = importlib.util.spec_from_file_location("openqp_api_helpers", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["openqp_api_helpers"] = module
    spec.loader.exec_module(module)
    return module


load_openqp_module = _load_api_helpers().load_openqp_module


class TestOpenQPXtbAPI(unittest.TestCase):
    def test_xtb_helper_defaults_to_mrsf_response(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2_xtb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .xtb(nstate=4, parameter_path="gfn1.opxtb")
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "xtb")
        self.assertEqual(config["input"]["functional"], "")
        self.assertEqual(config["xtb"]["type"], "mrsf")
        self.assertEqual(config["xtb"]["parameter_path"], "gfn1.opxtb")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["nstate"], "4")

    def test_xtb_helper_maps_response_and_noscc_aliases(self):
        openqp = load_openqp_module()

        ground = (
            openqp.OpenQP(project="h2_ground")
            .molecule("H 0 0 0; H 0 0 0.74")
            .xtb(response_type="ground", parameter_path="gfn1.opxtb")
        )
        config = ground.to_input_dict()
        self.assertEqual(config["xtb"]["type"], "ground")
        self.assertEqual(config["tdhf"]["type"], "tda")

        noscc = (
            openqp.OpenQP(project="h2_noscc")
            .molecule("H 0 0 0; H 0 0 0.74")
            .xtb(response_type="noscc", parameter_path="gfn1.opxtb")
        )
        config = noscc.to_input_dict()
        self.assertEqual(config["xtb"]["type"], "noscc")
        self.assertEqual(config["tdhf"]["type"], "tda")

        ground_noscc = (
            openqp.OpenQP(project="h2_ground_noscc")
            .molecule("H 0 0 0; H 0 0 0.74")
            .xtb(type="ground_noscc", parameter_path="gfn1.opxtb")
        )
        config = ground_noscc.to_input_dict()
        self.assertEqual(config["xtb"]["type"], "ground_noscc")
        self.assertEqual(config["tdhf"]["type"], "tda")

        with self.assertRaisesRegex(ValueError, "response_type"):
            openqp.OpenQP(project="bad_xtb").xtb(response_type="bogus")

    def test_theory_xtb_defaults_to_ground_state(self):
        openqp = load_openqp_module()

        ground = (
            openqp.OpenQP(project="h2_theory_xtb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory("xtb", parameter_path="gfn1.opxtb")
        )
        config = ground.to_input_dict()
        self.assertEqual(config["input"]["method"], "xtb")
        self.assertEqual(config["xtb"]["type"], "ground")
        self.assertEqual(config["tdhf"]["type"], "tda")

        tdxtb = (
            openqp.OpenQP(project="h2_theory_tdxtb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory("tdxtb", nstate=5, parameter_path="gfn1.opxtb")
        )
        config = tdxtb.to_input_dict()
        self.assertEqual(config["xtb"]["type"], "tddftb")
        self.assertEqual(config["tdhf"]["type"], "tda")
        self.assertEqual(config["tdhf"]["nstate"], "5")

    def test_theory_namespace_xtb_response_families(self):
        openqp = load_openqp_module()

        mrsf = (
            openqp.OpenQP(project="h2_mrsf_xtb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory.xtb(response_type="mrsf", parameter_path="gfn1.opxtb")
        )
        config = mrsf.to_input_dict()
        self.assertEqual(config["xtb"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["type"], "mrsf")

        sf = (
            openqp.OpenQP(project="h2_sf_xtb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory("sf-xtb", parameter_path="gfn1.opxtb")
        )
        config = sf.to_input_dict()
        self.assertEqual(config["xtb"]["type"], "sf")
        self.assertEqual(config["tdhf"]["type"], "sf")

    def test_workflow_soc_accepts_xtb_mrsf(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="xtb_soc")
            .molecule("H 0 0 0; H 0 0 0.74")
            .xtb(response_type="mrsf", parameter_path="gfn1.opxtb")
        )
        job.workflow.soc(soc_2e=1)

        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "soc")
        self.assertEqual(config["input"]["soc_2e"], "1")
        self.assertEqual(config["xtb"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["type"], "mrsf")

    def test_workflow_namd_accepts_xtb_mrsf(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="xtb_namd")
            .molecule("H 0 0 0; H 0 0 0.74")
            .xtb(response_type="mrsf", parameter_path="gfn1.opxtb")
        )
        job.workflow.namd(nstep=10, dt=0.5, init_state="S1")

        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "namd")
        self.assertEqual(config["md"]["nstep"], "10")
        self.assertEqual(config["tdhf"]["type"], "mrsf")

    def test_workflow_soc_rejects_ground_state_xtb(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="xtb_soc_bad")
            .molecule("H 0 0 0; H 0 0 0.74")
            .xtb(response_type="ground", parameter_path="gfn1.opxtb")
        )
        with self.assertRaisesRegex(ValueError, "only with MRSF"):
            job.workflow.soc(soc_2e=1)


if __name__ == "__main__":
    unittest.main()
