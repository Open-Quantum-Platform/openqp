"""Verify that devkit tools resolve resources from the private superset."""

import importlib.util
from pathlib import Path


DEVKIT_ROOT = Path(__file__).resolve().parents[1]
PATHS_MODULE = DEVKIT_ROOT / "tools" / "_repo_paths.py"


def _load_paths_module():
    spec = importlib.util.spec_from_file_location("devkit_repo_paths", PATHS_MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_integrated_paths_point_to_engine_resources():
    paths = _load_paths_module()
    engine_root = DEVKIT_ROOT.parent

    assert paths.DEVKIT_ROOT == DEVKIT_ROOT
    assert paths.ENGINE_ROOT == engine_root
    assert paths.EXAMPLES_ROOT == engine_root / "examples"
    assert paths.EXAMPLES_ROOT.is_dir()
    assert paths.SCHEMA_FILE == (
        engine_root / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
    )
    assert paths.SCHEMA_FILE.is_file()
