"""Python-API coverage for the optional OpenQP-DFTB backend helpers.

Reuses the stubbed-schema module loader from ``tests/test_openqp_api.py`` so the
compact ``job.dftb(...)`` / ``job.theory('dftb', ...)`` entry points and the
DFTB-MRSF workflow guards are exercised without the compiled library.
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


class TestOpenQPDftbAPI(unittest.TestCase):
    def test_dftb_helper_preserves_legacy_mrsf_default(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="h2_dftb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .dftb(nstate=4, parameter_path="mio-1-1")
        )

        config = job.to_input_dict()
        self.assertEqual(config["input"]["method"], "dftb")
        self.assertEqual(config["input"]["functional"], "")
        self.assertEqual(config["dftb"]["type"], "mrsf")
        self.assertEqual(config["dftb"]["parameter_path"], "mio-1-1")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["nstate"], "4")
        self.assertEqual(config["scf"]["type"], "rohf")

        ground = (
            openqp.OpenQP(project="h2_explicit_ground_dftb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .ground_dftb(parameter_path="mio-1-1")
        )
        self.assertEqual(ground.to_input_dict()["dftb"]["type"], "ground")

    def test_dftb_helper_maps_response_and_noscc_aliases(self):
        openqp = load_openqp_module()

        ground = (
            openqp.OpenQP(project="h2_ground")
            .molecule("H 0 0 0; H 0 0 0.74")
            .dftb(response_type="ground", parameter_path="mio-1-1")
        )
        config = ground.to_input_dict()
        self.assertEqual(config["dftb"]["type"], "ground")
        self.assertEqual(config["tdhf"]["type"], "tda")

        noscc = (
            openqp.OpenQP(project="h2_noscc")
            .molecule("H 0 0 0; H 0 0 0.74")
            .dftb(response_type="noscc", parameter_path="mio-1-1")
        )
        config = noscc.to_input_dict()
        self.assertEqual(config["dftb"]["type"], "ground_noscc")
        self.assertEqual(config["tdhf"]["type"], "tda")

        ground_noscc = (
            openqp.OpenQP(project="h2_ground_noscc")
            .molecule("H 0 0 0; H 0 0 0.74")
            .dftb(type="ground_noscc", parameter_path="mio-1-1")
        )
        config = ground_noscc.to_input_dict()
        self.assertEqual(config["dftb"]["type"], "ground_noscc")
        self.assertEqual(config["tdhf"]["type"], "tda")

        with self.assertRaisesRegex(ValueError, "response_type"):
            openqp.OpenQP(project="bad_dftb").dftb(response_type="bogus")

    def test_theory_dftb_defaults_to_ground_state(self):
        openqp = load_openqp_module()

        ground = (
            openqp.OpenQP(project="h2_theory_dftb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory("dftb", parameter_path="mio-1-1")
        )
        config = ground.to_input_dict()
        self.assertEqual(config["input"]["method"], "dftb")
        self.assertEqual(config["dftb"]["type"], "ground")
        self.assertEqual(config["tdhf"]["type"], "tda")

        tddftb = (
            openqp.OpenQP(project="h2_theory_tddftb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory("tddftb", nstate=5, parameter_path="mio-1-1")
        )
        config = tddftb.to_input_dict()
        self.assertEqual(config["dftb"]["type"], "tddftb")
        self.assertEqual(config["tdhf"]["type"], "tda")
        self.assertEqual(config["tdhf"]["nstate"], "5")

    def test_theory_namespace_dftb_response_families(self):
        openqp = load_openqp_module()

        legacy_mrsf = (
            openqp.OpenQP(project="h2_dftb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory.dftb(parameter_path="mio-1-1")
        )
        self.assertEqual(legacy_mrsf.to_input_dict()["dftb"]["type"], "mrsf")

        ground = (
            openqp.OpenQP(project="h2_ground_dftb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory.ground_dftb(parameter_path="mio-1-1")
        )
        self.assertEqual(ground.to_input_dict()["dftb"]["type"], "ground")

        mrsf = (
            openqp.OpenQP(project="h2_mrsf_dftb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory.mrsf_tddftb(parameter_path="mio-1-1")
        )
        config = mrsf.to_input_dict()
        self.assertEqual(config["dftb"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["type"], "mrsf")
        self.assertEqual(config["scf"]["type"], "rohf")
        self.assertEqual(config["scf"]["multiplicity"], "3")

        sf = (
            openqp.OpenQP(project="h2_sf_dftb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory.sf_tddftb(parameter_path="mio-1-1")
        )
        config = sf.to_input_dict()
        self.assertEqual(config["dftb"]["type"], "sf")
        self.assertEqual(config["tdhf"]["type"], "sf")

        conventional = (
            openqp.OpenQP(project="h2_tddftb")
            .molecule("H 0 0 0; H 0 0 0.74")
            .theory.tddftb(nstate=2, parameter_path="mio-1-1")
        )
        config = conventional.to_input_dict()
        self.assertEqual(config["dftb"]["type"], "tddftb")
        self.assertEqual(config["tdhf"]["type"], "tda")

    def test_top_level_dftb_response_helpers_are_symmetric(self):
        openqp = load_openqp_module()

        for helper, expected in (
            ("ground_dftb", "ground"),
            ("tddftb", "tddftb"),
            ("sf_tddftb", "sf"),
            ("mrsf_tddftb", "mrsf"),
        ):
            job = openqp.OpenQP(project=helper)
            getattr(job, helper)(nstate=2, parameter_path="mio-1-1")
            self.assertEqual(job.to_input_dict()["dftb"]["type"], expected)

    def test_workflow_soc_accepts_dftb_mrsf(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="dftb_soc")
            .molecule("H 0 0 0; H 0 0 0.74")
            .dftb(response_type="mrsf", parameter_path="mio-1-1")
        )
        job.workflow.soc(soc_2e=1)

        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "soc")
        self.assertEqual(config["input"]["soc_2e"], "1")
        self.assertEqual(config["dftb"]["type"], "mrsf")
        self.assertEqual(config["tdhf"]["type"], "mrsf")

    def test_workflow_namd_accepts_dftb_mrsf(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="dftb_namd")
            .molecule("H 0 0 0; H 0 0 0.74")
            .dftb(response_type="mrsf", parameter_path="mio-1-1")
        )
        job.workflow.namd(nstep=10, dt=0.5, init_state="S1")

        config = job.to_input_dict()
        self.assertEqual(config["input"]["runtype"], "namd")
        self.assertEqual(config["md"]["nstep"], "10")
        self.assertEqual(config["tdhf"]["type"], "mrsf")

    def test_workflow_soc_rejects_ground_state_dftb(self):
        openqp = load_openqp_module()
        job = (
            openqp.OpenQP(project="dftb_soc_bad")
            .molecule("H 0 0 0; H 0 0 0.74")
            .dftb(response_type="ground", parameter_path="mio-1-1")
        )
        with self.assertRaisesRegex(ValueError, "only with MRSF"):
            job.workflow.soc(soc_2e=1)


if __name__ == "__main__":
    unittest.main()
