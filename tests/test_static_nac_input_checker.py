"""Static NAC input checks through the real schema and checker, without liboqp."""

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
        "nac": {"type": "numerical", "states": "1 2"},
    }


def _report(inputs, legacy):
    config = inputs.convert(legacy)
    return inputs.checker.check_input_values(config, raise_error=False, emit=False)


@pytest.mark.parametrize("runtype", ["nac", "bp"])
def test_analytical_nac_is_refused(inputs, runtype):
    legacy = _legacy()
    legacy["input"]["runtype"] = runtype
    legacy["nac"]["type"] = "analytical"
    report = _report(inputs, legacy)
    assert not report.ok
    assert "analytical NAC vectors are not available" in report.to_text()


@pytest.mark.parametrize("driver", ["nac", "bp"])
def test_semantic_analytical_nac_is_refused(inputs, driver):
    with pytest.raises(inputs.semantic.OQPInputError, match="type must be numerical"):
        inputs.semantic.parse_canonical_oqp(
            f'mrsf(nstate=3)/bhhlyp/6-31g geom="h2o.xyz" {driver}(S0,S1,type=analytical)'
        )


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
    legacy["nac"]["type"] = "analytical"
    report = _report(inputs, legacy)
    assert not report.ok
    assert "analytical NAC vectors are not available" in report.to_text()
