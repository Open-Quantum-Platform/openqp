"""ACID input acceptance through the real schema and checker, without liboqp.

The .oqp and Python surfaces have their own tests; this covers the legacy
[properties] path, where the ordering rule, the GIAO requirement and the
grid values are enforced before any electronic structure is computed.
"""

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
    checker = _load("_acid_checker", "pyoqp/oqp/utils/input_checker.py", monkeypatch)
    semantic = _load("_acid_semantic", "pyoqp/oqp/utils/oqp_input.py", monkeypatch)
    parser = _load("_acid_parser", "pyoqp/oqp/utils/input_parser.py", monkeypatch)

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


def _legacy(**properties):
    props = {"scf_prop": "nmr,acid", "nmr_gauge": "giao"}
    props.update(properties)
    return {
        "input": {"method": "hf", "runtype": "energy", "basis": "sto-3g",
                  "system": str(ROOT / "examples/geometries/H2O-0381125c86f2.xyz")},
        "scf": {"type": "rhf", "multiplicity": "1"},
        "properties": props,
    }


def _report(inputs, legacy):
    config = inputs.convert(legacy)
    return inputs.checker.check_input_values(config, raise_error=False, emit=False)


def _errors(inputs, **properties):
    report = _report(inputs, _legacy(**properties))
    return [e for e in report.to_text().splitlines() if "properties." in e]


def test_acid_with_giao_and_sane_grid_passes(inputs):
    report = _report(inputs, _legacy())
    assert report.ok, report.to_text()
    report = _report(inputs, _legacy(acid_spacing="0.5", acid_padding="3.0"))
    assert report.ok, report.to_text()


def test_acid_needs_nmr_first(inputs):
    assert any("scf_prop" in e for e in _errors(inputs, scf_prop="acid"))
    assert any("scf_prop" in e for e in _errors(inputs, scf_prop="acid,nmr"))


def test_acid_requires_giao(inputs):
    assert any("nmr_gauge" in e for e in _errors(inputs, nmr_gauge="cgo"))


@pytest.mark.parametrize("spacing", ["0.0", "-0.2", "nan", "inf"])
def test_bad_spacing_is_rejected_before_any_scf(inputs, spacing):
    """The run would otherwise pay for the SCF and the GIAO response first."""
    assert any("acid_spacing" in e for e in _errors(inputs, acid_spacing=spacing))


@pytest.mark.parametrize("padding", ["-1.0", "nan", "-inf"])
def test_bad_padding_is_rejected_before_any_scf(inputs, padding):
    assert any("acid_padding" in e for e in _errors(inputs, acid_padding=padding))


def test_zero_padding_is_allowed(inputs):
    """A zero-padding box is tight but well defined; only negative is wrong."""
    assert not any("acid_padding" in e for e in _errors(inputs, acid_padding="0.0"))
