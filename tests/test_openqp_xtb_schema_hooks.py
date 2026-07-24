"""Schema, input-checker, and dispatch checks for the OpenQP-xTB runtime hook.

Mirrors tests/test_openqp_dftb_schema_hooks.py for method=xtb / the [xtb]
section, plus the GFN1-specific keys of the C ABI v3 (model, dispersion,
halogen_bond, third_order, spin_scale, and the 'ok' lc_gamma kind).
"""

import importlib
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
PYOQP = ROOT / "pyoqp"
if PYOQP.is_dir() and str(PYOQP) not in sys.path:
    sys.path.insert(0, str(PYOQP))


def _import_or_skip(module_name):
    try:
        return importlib.import_module(module_name)
    except RuntimeError as exc:
        text = str(exc)
        # Skip (rather than fail) when the native runtime is unavailable: either
        # liboqp cannot be located, or a dependency-light environment has mpi4py
        # importable but no libmpi to load.
        if (
            "Cannot locate OpenQP runtime files" in text
            or "OPENQP_ROOT does not contain" in text
            or "cannot load MPI library" in text
        ):
            for name in list(sys.modules):
                if name == "oqp" or name.startswith("oqp."):
                    sys.modules.pop(name, None)
            raise unittest.SkipTest(f"OpenQP native runtime is not available: {exc}") from exc
        raise


def _base_xtb_config(*, runtype="meci", grad=None):
    if grad is None:
        grad = [1, 2]
    return {
        "input": {
            "method": "xtb",
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
        "xtb": {
            "backend": "native",
            "type": "mrsf",
            "parameter_path": "/tmp/gfn1.opxtb",
            "model": "gfn1",
            "dispersion": True,
            "halogen_bond": True,
            "third_order": True,
            "spin_scale": 1.0,
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
            "lc_gamma": "ok",
            "timeout": 300,
        },
    }


class _FakeMol:
    """Minimal molecule stand-in for adapter construction (no runtime calls)."""

    def __init__(self, method="xtb", section=None):
        self.config = {
            "input": {"method": method, "runtype": "energy"},
            "tdhf": {"type": "tda", "nstate": 1},
            "properties": {"grad": [0]},
            method: dict(section or {}),
        }
        self.data = {"natom": 2}

    def get_atoms(self):
        return [1, 1]

    def get_system(self):
        return [0.0, 0.0, 0.0, 0.0, 0.0, 1.4]


class OpenQPXTBSchemaHookTests(unittest.TestCase):
    def test_input_schema_defines_xtb_section(self):
        OQP_CONFIG_SCHEMA = _import_or_skip("oqp.molecule.oqpdata").OQP_CONFIG_SCHEMA

        self.assertIn("xtb", OQP_CONFIG_SCHEMA)
        xtb = OQP_CONFIG_SCHEMA["xtb"]
        for key in (
            "backend",
            "type",
            "parameter_path",
            "library_path",
            "model",
            "dispersion",
            "halogen_bond",
            "third_order",
            "spin_scale",
            "scc_mixer",
            "scc_tolerance",
            "spc",
            "spc_coco",
            "spc_ovov",
            "spc_coov",
            "omega",
            "cam_alpha",
            "cam_beta",
            "lc_gamma",
            "lc_ground_state",
            "zvector",
            "spin_complete",
            "target_multiplicity",
            "response_solver",
            "mrsf_shift_oo",
            "timeout",
        ):
            self.assertIn(key, xtb)
        # DFTB-only probe executable must NOT leak into [xtb].
        self.assertNotIn("executable", xtb)
        self.assertEqual("native", xtb["backend"]["default"])
        self.assertEqual("gfn1", xtb["model"]["default"])
        self.assertEqual("ok", xtb["lc_gamma"]["default"])
        self.assertEqual(0.3, xtb["omega"]["type"](xtb["omega"]["default"]))
        self.assertEqual(0.0, xtb["cam_alpha"]["type"](xtb["cam_alpha"]["default"]))
        self.assertEqual(1.0, xtb["cam_beta"]["type"](xtb["cam_beta"]["default"]))
        self.assertEqual(1.0, xtb["spin_scale"]["type"](xtb["spin_scale"]["default"]))
        self.assertEqual("auto", xtb["response_solver"]["default"])

    def test_parser_accepts_xtb_input_section(self):
        OQP_CONFIG_SCHEMA = _import_or_skip("oqp.molecule.oqpdata").OQP_CONFIG_SCHEMA
        OQPConfigParser = _import_or_skip("oqp.utils.input_parser").OQPConfigParser

        parser = OQPConfigParser(schema=OQP_CONFIG_SCHEMA, allow_no_value=True)
        parser.read_string(
            """
[input]
method = xtb
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

[xtb]
backend = native
type = mrsf
parameter_path = /tmp/gfn1.opxtb
model = gfn1
dispersion = false
spin_scale = 1.25
lc_gamma = ok
"""
        )
        config = parser.validate()

        self.assertEqual("xtb", config["input"]["method"])
        self.assertEqual("meci", config["input"]["runtype"])
        self.assertEqual("native", config["xtb"]["backend"])
        self.assertEqual("mrsf", config["xtb"]["type"])
        self.assertEqual("gfn1", config["xtb"]["model"])
        self.assertFalse(config["xtb"]["dispersion"])
        self.assertEqual(1.25, config["xtb"]["spin_scale"])
        self.assertEqual("ok", config["xtb"]["lc_gamma"])

    def test_input_checker_allows_xtb_meci_response(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        config = _base_xtb_config()

        report = check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

    def test_input_checker_allows_xtb_data_for_namd_surface_data(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        config = _base_xtb_config(runtype="data", grad=[1, 2])

        report = check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

    def test_input_checker_accepts_xtb_nacme_with_mrsf_overlap_backend(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        # nacme with previous geometry passes.
        config = _base_xtb_config(runtype="nacme")
        config["input"]["system2"] = "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.71"
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

        # nacme without previous geometry is rejected (BasisOverlap needs it).
        config = _base_xtb_config(runtype="nacme")
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("previous-step information", report.to_text())

        # non-MRSF response type is still rejected for nacme.
        config = _base_xtb_config(runtype="nacme")
        config["input"]["system2"] = "\nH 0.0 0.0 0.0\nH 0.0 0.0 0.71"
        config["tdhf"]["type"] = "sf"
        config["xtb"]["type"] = "sf"
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("state-overlap backend", report.to_text())

    def test_input_checker_allows_xtb_soc_and_namd_for_mrsf(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        config = _base_xtb_config(runtype="soc")
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

        config = _base_xtb_config(runtype="namd")
        config["md"] = {"active": 1, "nstep": 5}
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

    def test_input_checker_blocks_xtb_md_like_dftb(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        for method in ("xtb", "dftb"):
            config = _base_xtb_config(runtype="md")
            config["input"]["method"] = method
            if method == "dftb":
                config["dftb"] = {"backend": "native", "type": "mrsf",
                                  "parameter_path": "/tmp/minimal.opdftb"}
            report = check_input_values(config, raise_error=False, emit=False)
            self.assertFalse(report.ok)
            self.assertIn("recognized but not implemented", report.to_text())

    def test_input_checker_rejects_xtb_probe_backend(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        config = _base_xtb_config()
        config["xtb"]["backend"] = "probe"
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("probe backend", report.to_text())
        self.assertIn("DFTB-only", report.to_text())

    def test_input_checker_rejects_bad_xtb_model_and_spin_scale(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        config = _base_xtb_config()
        config["xtb"]["model"] = "gfn2"
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("xtb.model", report.to_text())

        config = _base_xtb_config()
        config["xtb"]["spin_scale"] = 0.0
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("xtb.spin_scale", report.to_text())

    def test_input_checker_lc_gamma_ok_is_xtb_only(self):
        check_input_values = _import_or_skip("oqp.utils.input_checker").check_input_values

        # 'ok' accepted for [xtb] (it is the default) ...
        config = _base_xtb_config()
        config["xtb"]["lc_gamma"] = "ok"
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertTrue(report.ok, report.to_text())

        # ... a bogus kernel is rejected for [xtb] ...
        config = _base_xtb_config()
        config["xtb"]["lc_gamma"] = "gauss"
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("xtb.lc_gamma", report.to_text())

        # ... and 'ok' stays rejected for [dftb].
        config = _base_xtb_config()
        config["input"]["method"] = "dftb"
        config["dftb"] = {"backend": "native", "type": "mrsf",
                          "parameter_path": "/tmp/minimal.opdftb", "lc_gamma": "ok"}
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertFalse(report.ok)
        self.assertIn("dftb.lc_gamma", report.to_text())

    def test_method_resolution_table(self):
        tb = _import_or_skip("oqp.utils.tb_backends")

        self.assertEqual({"dftb", "xtb"}, set(tb.TB_METHODS))
        for method, expected in (
            ("dftb", True), ("xtb", True), ("XTB", True), ("DFTB", True),
            ("hf", False), ("tdhf", False), ("mp2", False), ("", False),
        ):
            self.assertEqual(expected, tb.is_tb_method(method), method)

        cfg_x = {"input": {"method": "xtb"}, "xtb": {"omega": 0.25}}
        cfg_d = {"input": {"method": "dftb"}, "dftb": {"omega": 0.30}}
        self.assertEqual("xtb", tb.tb_section_name(cfg_x))
        self.assertEqual("dftb", tb.tb_section_name(cfg_d))
        self.assertEqual({"omega": 0.25}, tb.tb_config(cfg_x))
        self.assertEqual({"omega": 0.30}, tb.tb_config(cfg_d))
        self.assertEqual({}, tb.tb_config({"input": {"method": "xtb"}}))
        with self.assertRaises(ValueError):
            tb.tb_section_name({"input": {"method": "hf"}})

    def test_make_tb_adapter_resolves_backend_class(self):
        tb = _import_or_skip("oqp.utils.tb_backends")
        openqp_dftb = _import_or_skip("oqp.library.openqp_dftb")
        openqp_xtb = _import_or_skip("oqp.library.openqp_xtb")

        mol_x = _FakeMol("xtb")
        adapter_x = tb.make_tb_adapter(mol_x)
        self.assertIsInstance(adapter_x, openqp_xtb.OpenQPXTBAdapter)
        # is_tb_method also accepts a molecule-like object.
        self.assertTrue(tb.is_tb_method(mol_x))

        mol_d = _FakeMol("dftb")
        adapter_d = tb.make_tb_adapter(mol_d)
        self.assertIsInstance(adapter_d, openqp_dftb.OpenQPDFTBAdapter)
        self.assertNotIsInstance(adapter_d, openqp_xtb.OpenQPXTBAdapter)

    def test_xtb_adapter_backend_identity_attributes(self):
        openqp_xtb = _import_or_skip("oqp.library.openqp_xtb")
        cls = openqp_xtb.OpenQPXTBAdapter

        self.assertEqual("xtb", cls.SECTION)
        self.assertEqual("openqp_xtb", cls.SYMBOL_PREFIX)
        self.assertEqual(("libopenqp_xtb_c.dylib", "libopenqp_xtb_c.so"), cls.LIB_BASENAMES)
        self.assertEqual("OPENQP_XTB_LIBRARY", cls.ENV_LIBRARY)
        self.assertEqual("OPENQP_XTB_PARAMETER_PATH", cls.ENV_PARAMETER)
        self.assertEqual("openqp_xtb", cls.PIP_LOCATOR)

    def test_xtb_adapter_lc_gamma_codes_and_probe_rejection(self):
        openqp_xtb = _import_or_skip("oqp.library.openqp_xtb")

        # default 'ok' -> 2; explicit kinds map per spec section 8a.
        self.assertEqual(2, openqp_xtb.OpenQPXTBAdapter(_FakeMol("xtb"))._lc_gamma_code())
        for kind, code in (("yukawa", 0), ("erf", 1), ("ok", 2)):
            adapter = openqp_xtb.OpenQPXTBAdapter(_FakeMol("xtb", {"lc_gamma": kind}))
            self.assertEqual(code, adapter._lc_gamma_code())
        adapter = openqp_xtb.OpenQPXTBAdapter(_FakeMol("xtb", {"lc_gamma": "gauss"}))
        with self.assertRaisesRegex(ValueError, "lc_gamma"):
            adapter._lc_gamma_code()

        adapter = openqp_xtb.OpenQPXTBAdapter(_FakeMol("xtb"))
        with self.assertRaisesRegex(ValueError, "backend=probe is not supported"):
            adapter._run_probe("mrsf", 1)

    def test_runfunc_routes_soc_through_tb_helpers(self):
        text = (ROOT / "pyoqp" / "oqp" / "library" / "runfunc.py").read_text(encoding="utf-8")
        self.assertIn("if is_tb_method(mol.config['input']['method']):", text)
        self.assertIn("from oqp.library.openqp_xtb import xtb_soc as tb_soc", text)
        self.assertIn("from oqp.library.openqp_dftb import dftb_soc as tb_soc", text)


if __name__ == "__main__":
    unittest.main()
