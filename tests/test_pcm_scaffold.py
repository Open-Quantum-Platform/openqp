"""Input-scaffold tests for planned PCM solvent support."""

from pathlib import Path
import importlib.util
import sys
import types
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load_module(module_name, relative_path):
    spec = importlib.util.spec_from_file_location(module_name, ROOT / relative_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load {relative_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def install_minimal_oqp_stubs():
    numpy_stub = sys.modules.setdefault("numpy", types.ModuleType("numpy"))
    for name in (
        "void",
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
        "float32",
        "float64",
        "complex64",
        "complex128",
    ):
        setattr(numpy_stub, name, type(name, (), {}))
    setattr(numpy_stub, "dtype", lambda value: value)

    oqp = sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    setattr(oqp, "ffi", object())
    setattr(oqp, "lib", object())

    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    constants = types.ModuleType("oqp.utils.constants")
    setattr(constants, "ANGSTROM_TO_BOHR", 1.8897261254576558)
    sys.modules["oqp.utils.constants"] = constants

    periodic_table = types.ModuleType("oqp.periodic_table")
    setattr(periodic_table, "MASSES", {})
    setattr(periodic_table, "SYMBOL_MAP", {})
    sys.modules["oqp.periodic_table"] = periodic_table

    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        use_mpi = False
        size = 1

    setattr(mpi_utils, "MPIManager", MPIManager)
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils


class PCMScaffoldTests(unittest.TestCase):
    def setUp(self):
        install_minimal_oqp_stubs()

    def test_schema_accepts_pcm_section_defaults(self):
        oqpdata = load_module(
            "oqpdata_pcm_schema_under_test",
            "pyoqp/oqp/molecule/oqpdata.py",
        )
        parser_mod = load_module(
            "input_parser_pcm_schema_under_test",
            "pyoqp/oqp/utils/input_parser.py",
        )

        parser = parser_mod.OQPConfigParser(schema=oqpdata.OQP_CONFIG_SCHEMA)
        parser.read_string(
            """
[input]
system=
 O 0.000000 0.000000 0.000000
 H 0.000000 0.757000 0.587000
 H 0.000000 -0.757000 0.587000
basis=6-31g*

[pcm]
enabled=false
backend=ddx
mode=reference_scf
model=ddpcm
epsilon=78.3553
"""
        )

        config = parser.validate()
        self.assertFalse(config["pcm"]["enabled"])
        self.assertEqual(config["pcm"]["backend"], "ddx")
        self.assertEqual(config["pcm"]["mode"], "reference_scf")
        self.assertAlmostEqual(config["pcm"]["epsilon"], 78.3553)

    def test_checker_blocks_enabled_pcm_until_runtime_exists(self):
        input_checker = load_module(
            "input_checker_pcm_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {
                "system": "\n O 0.0 0.0 0.0\n H 0.0 0.757 0.587\n H 0.0 -0.757 0.587",
                "basis": "6-31g*",
                "method": "hf",
                "runtype": "energy",
            },
            "guess": {"type": "huckel"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 1},
            "pcm": {
                "enabled": True,
                "backend": "ddx",
                "mode": "reference_scf",
                "model": "ddpcm",
                "epsilon": 78.3553,
            },
        }

        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        errors = {item.path: item.message for item in report.errors}
        self.assertIn("pcm.enabled", errors)
        self.assertIn("runtime solvent coupling is not implemented", errors["pcm.enabled"])

    def test_checker_rejects_pcm_gradients_for_first_scope(self):
        input_checker = load_module(
            "input_checker_pcm_grad_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {
                "system": "\n O 0.0 0.0 0.0\n H 0.0 0.757 0.587\n H 0.0 -0.757 0.587",
                "basis": "6-31g*",
                "method": "hf",
                "runtype": "grad",
            },
            "guess": {"type": "huckel"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 1},
            "pcm": {
                "enabled": True,
                "backend": "ddx",
                "mode": "reference_scf",
                "model": "ddpcm",
                "epsilon": 78.3553,
            },
        }

        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        errors = {item.path: item.message for item in report.errors}
        self.assertIn("input.runtype", errors)
        self.assertIn("single-point energies only", errors["input.runtype"])

    def test_checker_rejects_backend_model_mismatch(self):
        input_checker = load_module(
            "input_checker_pcm_model_backend_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {
                "system": "\n O 0.0 0.0 0.0\n H 0.0 0.757 0.587\n H 0.0 -0.757 0.587",
                "basis": "6-31g*",
                "method": "hf",
                "runtype": "energy",
            },
            "guess": {"type": "huckel"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 1},
            "pcm": {
                "enabled": False,
                "backend": "ddx",
                "mode": "reference_scf",
                "model": "iefpcm",
                "epsilon": 78.3553,
            },
        }

        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        errors = {item.path: item.message for item in report.errors}
        self.assertIn("pcm.model", errors)
        self.assertIn("is not supported by backend ddx", errors["pcm.model"])

    def test_checker_rejects_non_reference_scf_pcm_modes_until_post_state_is_implemented(self):
        input_checker = load_module(
            "input_checker_pcm_mode_scope_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {
                "system": "\n O 0.0 0.0 0.0\n H 0.0 0.757 0.587\n H 0.0 -0.757 0.587",
                "basis": "6-31g*",
                "method": "tdhf",
                "runtype": "energy",
            },
            "guess": {"type": "huckel"},
            "scf": {"type": "rohf", "multiplicity": 3},
            "tdhf": {"type": "mrsf", "nstate": 3},
            "pcm": {
                "enabled": False,
                "backend": "ddx",
                "mode": "post_state_correction",
                "model": "ddpcm",
                "epsilon": 78.3553,
            },
        }

        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        errors = {item.path: item.message for item in report.errors}
        self.assertIn("pcm.mode", errors)
        self.assertIn("Only reference_scf", errors["pcm.mode"])

    def test_checker_rejects_uhf_reference_for_first_pcm_scope(self):
        input_checker = load_module(
            "input_checker_pcm_scf_scope_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {
                "system": "\n O 0.0 0.0 0.0\n H 0.0 0.757 0.587\n H 0.0 -0.757 0.587",
                "basis": "6-31g*",
                "method": "hf",
                "runtype": "energy",
            },
            "guess": {"type": "huckel"},
            "scf": {"type": "uhf", "multiplicity": 3},
            "tdhf": {"type": "rpa", "nstate": 1},
            "pcm": {
                "enabled": False,
                "backend": "ddx",
                "mode": "reference_scf",
                "model": "ddpcm",
                "epsilon": 78.3553,
            },
        }

        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        errors = {item.path: item.message for item in report.errors}
        self.assertIn("scf.type", errors)
        self.assertIn("PCM first scope supports RHF/ROHF", errors["scf.type"])

    def test_parser_preserves_malformed_pcm_dielectric_for_checker_diagnostic(self):
        oqpdata = load_module(
            "oqpdata_pcm_bad_epsilon_schema_under_test",
            "pyoqp/oqp/molecule/oqpdata.py",
        )
        parser_mod = load_module(
            "input_parser_pcm_bad_epsilon_under_test",
            "pyoqp/oqp/utils/input_parser.py",
        )
        input_checker = load_module(
            "input_checker_pcm_bad_epsilon_parser_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )

        parser = parser_mod.OQPConfigParser(schema=oqpdata.OQP_CONFIG_SCHEMA)
        parser.read_string(
            """
[input]
system=
 O 0.000000 0.000000 0.000000
 H 0.000000 0.757000 0.587000
 H 0.000000 -0.757000 0.587000
basis=6-31g*
method=hf
runtype=energy

[scf]
type=rhf
multiplicity=1

[pcm]
enabled=false
backend=ddx
mode=reference_scf
model=ddpcm
epsilon=water
"""
        )

        config = parser.validate()
        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        errors = {item.path: item.message for item in report.errors}
        self.assertIn("pcm.epsilon", errors)
        self.assertIn("must be numeric", errors["pcm.epsilon"])

    def test_checker_rejects_non_numeric_pcm_dielectric_without_crashing(self):
        input_checker = load_module(
            "input_checker_pcm_epsilon_under_test",
            "pyoqp/oqp/utils/input_checker.py",
        )
        config = {
            "input": {
                "system": "\n O 0.0 0.0 0.0\n H 0.0 0.757 0.587\n H 0.0 -0.757 0.587",
                "basis": "6-31g*",
                "method": "hf",
                "runtype": "energy",
            },
            "guess": {"type": "huckel"},
            "scf": {"type": "rhf", "multiplicity": 1},
            "tdhf": {"type": "rpa", "nstate": 1},
            "pcm": {
                "enabled": False,
                "backend": "ddx",
                "mode": "reference_scf",
                "model": "ddpcm",
                "epsilon": "water",
            },
        }

        report = input_checker.check_input_values(config, raise_error=False, emit=False)
        errors = {item.path: item.message for item in report.errors}
        self.assertIn("pcm.epsilon", errors)
        self.assertIn("must be numeric", errors["pcm.epsilon"])

    def test_validation_matrix_documents_energy_only_mrsf_reference_target(self):
        matrix_path = ROOT / "docs" / "solvent_pcm_validation_matrix.md"
        self.assertTrue(matrix_path.exists(), "PCM validation matrix document is required")
        text = matrix_path.read_text(encoding="utf-8")
        self.assertIn("MRSF-TDDFT with PCM-solvated ROHF reference", text)
        self.assertIn("runtype=energy", text)
        self.assertIn("PySCF", text)
        self.assertIn("analytic PCM gradients", text)
        self.assertIn("ddX q_cav sign/scale is provisional", text)

    def test_validation_matrix_documents_incremental_fock_blocker(self):
        matrix_path = ROOT / "docs" / "solvent_pcm_validation_matrix.md"
        text = matrix_path.read_text(encoding="utf-8")
        self.assertIn("incremental Fock", text)
        self.assertIn("dens_old", text)
        self.assertIn("f_old", text)
        self.assertIn("reference PCM incremental Fock is not validated", text)

    def test_validation_matrix_documents_call_site_bridge_shape_contract(self):
        matrix_path = ROOT / "docs" / "solvent_pcm_validation_matrix.md"
        text = matrix_path.read_text(encoding="utf-8")
        self.assertIn("reference_scf_pcm_calc_fock_call_site_bridge", text)
        self.assertIn("disabled/no-payload", text)
        self.assertIn("pcm_reaction_potential_in", text)
        self.assertIn("packed AO length", text)
        self.assertIn("size(pcm_reaction_potential_in) == nbf * (nbf + 1) / 2", text)


if __name__ == "__main__":
    unittest.main()
