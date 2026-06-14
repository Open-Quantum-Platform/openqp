from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_python_runtime_prefers_package_local_root_over_openqp_root():
    source = (ROOT / "pyoqp" / "oqp" / "__init__.py").read_text()

    assert "def _resolve_oqp_root" in source
    assert "if _is_oqp_root(package_root, suffix):" in source
    assert "os.environ[\"OPENQP_ROOT\"] = package_root" in source
    assert "if _is_oqp_root(env_root, suffix):" in source


def test_cmake_forces_ilp64_when_stale_cache_sets_lp64():
    source = (ROOT / "CMakeLists.txt").read_text()

    assert "Forcing LINALG_LIB_INT64=ON" in source
    assert "set(LINALG_LIB_INT64 ON CACHE BOOL" in source
    assert "FORCE)" in source
