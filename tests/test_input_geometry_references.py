"""Regression tests for optional QM/MM input and geometry references."""

import importlib.util
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load_input_checker():
    oqp = sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    utils = sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils
    oqp.utils = utils
    utils.mpi_utils = mpi_utils

    path = ROOT / "pyoqp/oqp/utils/input_checker.py"
    spec = importlib.util.spec_from_file_location(
        "openqp_geometry_reference_checker_under_test", path
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_oqpdata():
    """Load the schema without requiring the compiled OpenQP extension."""
    oqp = types.ModuleType("oqp")
    oqp.ffi = types.SimpleNamespace()
    oqp.lib = types.SimpleNamespace()
    periodic_table = types.ModuleType("oqp.periodic_table")
    periodic_table.MASSES = [0.0, 1.0]
    periodic_table.SYMBOL_MAP = {"H": 1, "h": 1, "1": 1}
    constants = types.ModuleType("oqp.utils.constants")
    constants.ANGSTROM_TO_BOHR = 0.529177210903
    qmmm = types.ModuleType("oqp.utils.qmmm")

    saved = {
        name: sys.modules.get(name)
        for name in (
            "oqp", "oqp.periodic_table", "oqp.utils",
            "oqp.utils.constants", "oqp.utils.qmmm",
        )
    }
    try:
        utils = types.ModuleType("oqp.utils")
        sys.modules.update({
            "oqp": oqp,
            "oqp.periodic_table": periodic_table,
            "oqp.utils": utils,
            "oqp.utils.constants": constants,
            "oqp.utils.qmmm": qmmm,
        })
        path = ROOT / "pyoqp/oqp/molecule/oqpdata.py"
        spec = importlib.util.spec_from_file_location(
            "openqp_optional_qmmm_schema_under_test", path
        )
        assert spec and spec.loader
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module
    finally:
        for name, previous in saved.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous


def _load_input_parser():
    path = ROOT / "pyoqp/oqp/utils/input_parser.py"
    spec = importlib.util.spec_from_file_location(
        "openqp_optional_qmmm_parser_under_test", path
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_optional_qmmm_xyz_and_index_defaults_are_none():
    oqpdata = _load_oqpdata()
    input_parser = _load_input_parser()
    parser = input_parser.OQPConfigParser(
        schema=oqpdata.OQP_CONFIG_SCHEMA, allow_no_value=True
    )

    config = parser.validate()
    assert config["qmmm"]["qm_atoms_xyz"] is None
    assert config["qmmm"]["qm_list"] is None

    parser["qmmm"]["qm_atoms_xyz"] = " qm-region.xyz "
    parser["qmmm"]["qm_list"] = "0,2,4"
    config = parser.validate()
    assert config["qmmm"]["qm_atoms_xyz"] == "qm-region.xyz"
    assert config["qmmm"]["qm_list"] == "0,2,4"


def test_pdb_reference_with_qm_indices_is_not_inline_geometry(tmp_path):
    checker = _load_input_checker()
    pdb_path = tmp_path / "mol.pdb"
    pdb_path.write_text("END\n", encoding="utf-8")
    system = f"{pdb_path} 0 1 2"

    inline, reference = checker._iter_coordinate_lines(system)
    assert inline == []
    assert reference == system

    report = checker.CheckReport()
    checker._check_system({"input": {"system": system}}, report)
    assert report.errors == []


def test_xyz_reference_still_checks_file_existence(tmp_path):
    checker = _load_input_checker()
    missing = tmp_path / "missing.xyz"
    report = checker.CheckReport()

    checker._check_system({"input": {"system": f"{missing} 0 1 2"}}, report)

    assert any("does not exist" in item.message for item in report.errors)


def test_pdb_reference_still_validates_qm_indices(tmp_path):
    checker = _load_input_checker()
    pdb_path = tmp_path / "mol.pdb"
    pdb_path.write_text("END\n", encoding="utf-8")
    report = checker.CheckReport()

    checker._check_system(
        {"input": {"system": f"{pdb_path} 0 1 1"}}, report
    )

    assert any("Repeated entries" in item.message for item in report.errors)


def test_symbol_and_atomic_number_inline_geometries_are_preserved():
    checker = _load_input_checker()

    for system in (
        "H 0 0 0\nO 0 0 1",
        "8 0 0 0\n1 0 0 1",
    ):
        inline, reference = checker._iter_coordinate_lines(system)
        assert reference is None
        assert inline == system.splitlines()
