"""Schema and dispatch checks for the OpenQP-DFTB runtime hook."""

import importlib
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]


def _import_or_skip(module_name):
    try:
        return importlib.import_module(module_name)
    except RuntimeError as exc:
        text = str(exc)
        if "Cannot locate OpenQP runtime files" in text or "OPENQP_ROOT does not contain" in text:
            for name in list(sys.modules):
                if name == "oqp" or name.startswith("oqp."):
                    sys.modules.pop(name, None)
            raise unittest.SkipTest(f"OpenQP native runtime is not available: {exc}") from exc
        raise


def _base_dftb_config(*, runtype="meci", grad=None):
    if grad is None:
        grad = [1, 2]
    return {
        "input": {
            "method": "dftb",
            "runtype": runtype,
            "system": "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7",
            "basis": "6-31g*",
            "functional": "",
        },
        "guess": {"type": "huckel"},
        "pcm": {"enabled": False, "backend": "ddx", "mode": "reference_scf", "model": "ddpcm", "epsilon": 78.3553},
        "symmetry": {"enabled": "false", "tolerance": 1.0e-5, "point_group": "auto", "subgroup": "auto", "use_integral_symmetry": "false"},
        "scf": {"type": "rohf", "multiplicity": 3, "converger_type": "diis", "diis_type": "cdiis", "alternative_scf": "trah", "init_scf": "no"},
        "tdhf": {"type": "mrsf", "nstate": 2, "multiplicity": 1},
        "properties": {"grad": grad, "scf_prop": [], "td_prop": False, "nmr_gauge": "cgo", "nac": ""},
        "optimize": {"lib": "oqp", "istate": 1, "jstate": 2, "meci_search": "penalty"},
        "nac": {"states": [[1, 2]], "nproc": 1, "align": "reorder"},
        "dftb": {
            "backend": "native",
            "type": "mrsf",
            "parameter_path": "/tmp/minimal.opdftb",
            "scc_tolerance": 1.0e-8,
            "scc_mixer": "auto",
            "scc_mixing": 0.35,
            "scc_history": 12,
            "scc_max_step": 0.5,
            "max_scc_iterations": 1200,
            "response_tolerance": 1.0e-6,
            "response_max_iterations": 50,
            "response_max_subspace": 50,
            "spc": 0.5,
            "omega": 0.3,
            "cam_alpha": 0.0,
            "cam_beta": 1.0,
            "timeout": 300,
        },
    }


class OpenQPDFTBSchemaHookTests(unittest.TestCase):
    def test_input_schema_defines_dftb_section(self):
        OQP_CONFIG_SCHEMA = _import_or_skip("oqp.molecule.oqpdata").OQP_CONFIG_SCHEMA

        self.assertIn("dftb", OQP_CONFIG_SCHEMA)
        dftb = OQP_CONFIG_SCHEMA["dftb"]
        for key in (
            "backend",
            "type",
            "parameter_path",
            "scc_mixer",
            "scc_tolerance",
            "spc",
            "spc_coco",
            "spc_ovov",
            "spc_coov",
            "omega",
            "lc_gamma",
            "lc_ground_state",
            "zvector",
            "spin_complete",
            "target_multiplicity",
            "response_solver",
            "mrsf_shift_oo",
            "timeout",
        ):
            self.assertIn(key, dftb)
        self.assertEqual("native", dftb["backend"]["default"])
        self.assertEqual(0.5, dftb["spc"]["type"](dftb["spc"]["default"]))
        self.assertEqual("yukawa", dftb["lc_gamma"]["default"])
        self.assertEqual("auto", dftb["response_solver"]["default"])

    def test_parser_accepts_dftb_input_section(self):
        OQP_CONFIG_SCHEMA = _import_or_skip("oqp.molecule.oqpdata").OQP_CONFIG_SCHEMA
        OQPConfigParser = _import_or_skip("oqp.utils.input_parser").OQPConfigParser

        parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA, allow_no_value=True)
        parser.read_string(
            """
[input]
method = dftb
runtype = meci
system =
  H 0.0 0.0 0.0
  H 0.0 0.0 0.7

[tdhf]
type = mrsf
nstate = 2

[optimize]
istate = 1
jstate = 2

[dftb]
backend = native
type = mrsf
parameter_path = /tmp/minimal.opdftb
scc_mixer = auto
spc = 0.5
"""
        )
        config = parser.validate()

        self.assertEqual("dftb", config["input"]["method"])
        self.assertEqual("meci", config["input"]["runtype"])
        self.assertEqual("native", config["dftb"]["backend"])
        self.assertEqual("mrsf", config["dftb"]["type"])
        self.assertEqual(0.5, config["dftb"]["spc"])

    def test_input_checker_allows_dftb_meci_response(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        config = _base_dftb_config()

        report = check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

    def test_input_checker_allows_dftb_data_for_namd_surface_data(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        config = _base_dftb_config(runtype="data", grad=[1, 2])

        report = check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

    def test_input_checker_accepts_dftb_nacme_with_mrsf_overlap_backend(self):
        # The MRSF-TDDFTB state-overlap (TLF) backend provides nacme/nac/namd/soc;
        # nacme still needs previous-step geometry, and non-MRSF response types
        # remain rejected.
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        # nacme with previous geometry passes.
        config = _base_dftb_config(runtype="nacme")
        config["input"]["system2"] = "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.71"
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

        # nacme without previous geometry is rejected (BasisOverlap needs it).
        config = _base_dftb_config(runtype="nacme")
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("previous-step information", report.to_text())

        # non-MRSF response type is still rejected for nacme.
        config = _base_dftb_config(runtype="nacme")
        config["input"]["system2"] = "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.71"
        config["tdhf"]["type"] = "sf"
        config["dftb"]["type"] = "sf"
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("state-overlap backend", report.to_text())

    def test_single_point_routes_dftb_to_adapter(self):
        text = (ROOT / "pyoqp" / "oqp" / "library" / "single_point.py").read_text(encoding="utf-8")
        self.assertIn("from oqp.library.openqp_dftb import OpenQPDFTBAdapter", text)
        self.assertIn("if self.method == 'dftb':", text)
        self.assertIn("OpenQPDFTBAdapter(self.mol).energy()", text)
        self.assertIn("OpenQPDFTBAdapter(self.mol).gradient(self.grads)", text)
        self.assertIn("OpenQPDFTBAdapter(self.mol).reference()", text)
        self.assertIn("OpenQPDFTBAdapter(self.mol).excitation(ref_energy)", text)

    def test_adapter_exposes_reference_excitation_for_namd_driver_shape(self):
        text = (ROOT / "pyoqp" / "oqp" / "library" / "openqp_dftb.py").read_text(encoding="utf-8")
        self.assertIn("def reference(self):", text)
        self.assertIn("def excitation(self, _ref_energy=None):", text)
        self.assertIn('if runtype in {"grad", "data"}:', text)
        self.assertIn("need_grad=False", text)
        self.assertIn("need_grad=True", text)
        self.assertIn('elif backend == "probe":', text)
        self.assertNotIn("except _NativeUnavailable", text)


if __name__ == "__main__":
    unittest.main()
