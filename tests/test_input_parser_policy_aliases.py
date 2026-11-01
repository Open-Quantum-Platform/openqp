import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load_input_parser():
    path = ROOT / "pyoqp/oqp/utils/input_parser.py"
    spec = importlib.util.spec_from_file_location(
        "openqp_policy_alias_parser_under_test", path
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SCHEMA = {
    "input": {
        "system": {"type": str, "default": "H 0.0 0.0 0.0"},
    },
    "md": {
        "nacme_policy": {"type": str, "default": "off"},
        "nacme_policy_abs_tol": {"type": float, "default": "1.0e-4"},
        "nve_policy": {"type": str, "default": "warn"},
        "nve_policy_step_tol": {"type": float, "default": "1.0e-3"},
    },
}


def test_sectioned_input_accepts_public_policy_names():
    parser_module = _load_input_parser()
    parser = parser_module.OQPConfigParser(schema=SCHEMA)
    parser.read_string(
        "[md]\n"
        "nacme_policy = error\n"
        "nacme_policy_abs_tol = 2.0e-4\n"
        "nve_policy = off\n"
        "nve_policy_step_tol = 8.0e-4\n"
    )

    config = parser.validate()

    assert config["md"]["nacme_policy"] == "error"
    assert config["md"]["nacme_policy_abs_tol"] == pytest.approx(2.0e-4)
    assert config["md"]["nve_policy"] == "off"
    assert config["md"]["nve_policy_step_tol"] == pytest.approx(8.0e-4)
    assert not {key for key in parser["md"] if "gate" in key}


def test_sectioned_input_rejects_policy_and_nondefault_legacy_name():
    parser_module = _load_input_parser()
    parser = parser_module.OQPConfigParser(schema=SCHEMA)
    parser.read_string(
        "[md]\n"
        "nve_policy = error\n"
        "nve_gate = off\n"
    )

    with pytest.raises(ValueError, match="specify the same policy"):
        parser.validate()


def test_sectioned_input_lowers_legacy_gate_name_to_policy():
    parser_module = _load_input_parser()
    parser = parser_module.OQPConfigParser(schema=SCHEMA)
    parser.read_string("[md]\nnve_gate = off\n")

    config = parser.validate()

    assert config["md"]["nve_policy"] == "off"
    assert "nve_gate" not in parser["md"]
