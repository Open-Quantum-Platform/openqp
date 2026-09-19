"""Static NAC input acceptance through the real schema and checker, without liboqp."""

import ast
import importlib.util
from pathlib import Path
import sys
import types

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load(name, relative_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def inputs(monkeypatch):
    # Only MPI discovery is replaced. Parser defaults, type conversion and all
    # scientific input checks below execute the production Python source.
    mpi = types.ModuleType("oqp.utils.mpi_utils")
    mpi.MPIManager = lambda: types.SimpleNamespace(use_mpi=False, size=1)
    monkeypatch.setitem(sys.modules, "oqp.utils.mpi_utils", mpi)
    checker = _load("_static_nac_checker", "pyoqp/oqp/utils/input_checker.py", monkeypatch)
    semantic = _load("_static_nac_semantic", "pyoqp/oqp/utils/oqp_input.py", monkeypatch)
    parser = _load("_static_nac_parser", "pyoqp/oqp/utils/input_parser.py", monkeypatch)

    # The schema and its converters precede the native OQPData class. Load
    # that exact source without importing the CFFI library or molecular data.
    path = ROOT / "pyoqp/oqp/molecule/oqpdata.py"
    tree = ast.parse(path.read_text(), filename=str(path))
    nodes = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            nodes.append(node)
        elif isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == "OQP_CONFIG_SCHEMA"
                for target in node.targets):
            nodes.append(node)
            break
    namespace = {"Path": Path}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), namespace)
    schema = namespace["OQP_CONFIG_SCHEMA"]

    def convert(legacy):
        config = parser.OQPConfigParser(schema=schema)
        config.load_dict(legacy)
        return config.validate()

    return types.SimpleNamespace(checker=checker, semantic=semantic, convert=convert)


def _legacy():
    return {
        "input": {"method": "tdhf", "runtype": "nac", "basis": "6-31g",
                  "functional": "bhhlyp",
                  "system": str(ROOT / "examples/geometries/H2O-0381125c86f2.xyz")},
        "scf": {"type": "rohf", "multiplicity": "3", "conv": "1e-10"},
        "tdhf": {"type": "mrsf", "multiplicity": "1", "nstate": "3", "conv": "1e-10"},
        "nac": {"type": "analytical", "states": "1 2"},
    }


def _report(inputs, legacy):
    config = inputs.convert(legacy)
    return inputs.checker.check_input_values(config, raise_error=False, emit=False)


@pytest.mark.parametrize("driver", ["nac", "bp"])
@pytest.mark.parametrize("nac_type", ["analytical", "ANALYTICAL"])
def test_semantic_static_nac_passes_full_checker(inputs, driver, nac_type):
    text = (
        'mrsf(nstate=3)/bhhlyp/6-31g geom="h2o.xyz" '
        f'{driver}(S0,S1,type={nac_type}) scf(conv=1e-8) tdhf(conv=1e-10)'
    )
    spec = inputs.semantic.parse_canonical_oqp(text)
    legacy = inputs.semantic.lower_to_legacy(spec)
    legacy["input"]["system"] = _legacy()["input"]["system"]
    config = inputs.convert(legacy)
    assert config["nac"]["type"] == "analytical"
    assert config["nac"]["states"] == [[1, 2]]
    report = inputs.checker.check_input_values(config, raise_error=False, emit=False)
    assert report.ok, report.to_text()


def test_example_passes_full_checker(inputs):
    path = ROOT / "examples/other/h2o_nac_analytical_mrsf.oqp"
    spec = inputs.semantic.parse_canonical_oqp(path.read_text())
    legacy = inputs.semantic.lower_to_legacy(spec, source_dir=path.parent)
    report = _report(inputs, legacy)
    assert report.ok, report.to_text()


@pytest.mark.parametrize("runtype", ["nac", "bp"])
@pytest.mark.parametrize("conv", ["1e-8", "1e-10"])
def test_legacy_static_analytic_nac_accepts_threshold_boundary(inputs, runtype, conv):
    legacy = _legacy()
    legacy["input"]["runtype"] = runtype
    legacy["scf"]["conv"] = legacy["tdhf"]["conv"] = conv
    report = _report(inputs, legacy)
    assert report.ok, report.to_text()


@pytest.mark.parametrize("section, key, value, message", [
    ("input", "method", "hf", "NAC workflows require method=tdhf"),
    ("tdhf", "type", "rpa", "require tdhf.type=mrsf"),
    ("tdhf", "type", "umrsf", "UMRSF-TDDFT only supports runtype=energy"),
    ("scf", "type", "uhf", "requires an ROHF/ROKS reference"),
    ("scf", "multiplicity", "5", "requires a two-SOMO triplet reference"),
    ("tdhf", "multiplicity", "3", "implements singlet states only"),
    ("tdhf", "nstate", "1", "requires at least two states"),
    ("nac", "type", "analytic", "Unknown NAC vector type"),
    ("nac", "states", "1 1", "requires two distinct states"),
    ("nac", "states", "0 1", "indices >= 1"),
    ("nac", "states", "1 4", "exceeds the number"),
])
def test_legacy_analytic_nac_rejects_unsupported_cases(inputs, section, key, value, message):
    legacy = _legacy()
    legacy[section][key] = value
    report = _report(inputs, legacy)
    assert not report.ok
    assert message in report.to_text()


@pytest.mark.parametrize("section", ["scf", "tdhf"])
@pytest.mark.parametrize("conv", [None, "1e-6", "0", "-1e-10", "nan", "inf"])
def test_legacy_analytic_nac_rejects_bad_convergence(inputs, section, conv):
    legacy = _legacy()
    if conv is None:
        del legacy[section]["conv"]  # Exercise the real schema default.
    else:
        legacy[section]["conv"] = conv
    report = _report(inputs, legacy)
    assert not report.ok
    assert f"{section}.conv" in report.to_text()
    assert "finite, positive convergence threshold <= 1e-8" in report.to_text()


@pytest.mark.parametrize("runtype", ["nac", "bp", "nacme"])
def test_numerical_nac_and_overlap_nacme_remain_accepted(inputs, runtype):
    legacy = _legacy()
    legacy["input"]["runtype"] = runtype
    legacy["input"]["system2"] = legacy["input"]["system"]
    legacy["nac"]["type"] = "numerical"
    del legacy["scf"]["conv"]
    del legacy["tdhf"]["conv"]
    report = _report(inputs, legacy)
    assert report.ok, report.to_text()


def test_analytical_type_is_not_a_scalar_nacme_option(inputs):
    legacy = _legacy()
    legacy["input"]["runtype"] = "nacme"
    legacy["input"]["system2"] = legacy["input"]["system"]
    report = _report(inputs, legacy)
    assert not report.ok
    assert "NACME uses overlap time-derivative couplings" in report.to_text()
