"""The import-state guard in tests/conftest.py must fail every kind of leak.

Each case runs in a child pytest session that loads only the guard, so the
leaks it provokes never reach this session.  A leak that the guard misses
would again make results depend on test-file order (a test importing the
unbuilt source tree, or a later test running against a stub).
"""
import os
import subprocess
import sys
import textwrap
from pathlib import Path

CONFTEST = Path(__file__).resolve().with_name("conftest.py")
PYOQP = CONFTEST.parents[1] / "pyoqp"

CHILD_CONFTEST = f"""
import importlib.util

_spec = importlib.util.spec_from_file_location("oqp_tests_conftest", {str(CONFTEST)!r})
_guard = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_guard)
pytest_make_collect_report = _guard.pytest_make_collect_report
_guard_import_state = _guard._guard_import_state
"""

MODULE_LEVEL_INSERT = f"""
import sys
sys.path.insert(0, {str(PYOQP)!r})

def test_never_reached():
    pass
"""

PER_TEST = f"""
import os
import sys
import types


def test_leaves_source_tree_on_sys_path():
    sys.path.insert(0, {str(PYOQP)!r})


def test_leaves_openqp_root_set():
    os.environ["OPENQP_ROOT"] = "/nonexistent/openqp-root"


def test_leaves_stub_module():
    sys.modules["oqp.zz_guard_stub"] = types.ModuleType("oqp.zz_guard_stub")


def test_scoped_changes_are_allowed(monkeypatch):
    monkeypatch.syspath_prepend({str(PYOQP)!r})
    monkeypatch.setenv("OPENQP_ROOT", "/nonexistent/openqp-root")
    monkeypatch.setitem(sys.modules, "oqp.zz_guard_stub",
                        types.ModuleType("oqp.zz_guard_stub"))
"""


def _run_child(tmp_path, files):
    (tmp_path / "conftest.py").write_text(CHILD_CONFTEST)
    for name, body in files.items():
        (tmp_path / name).write_text(textwrap.dedent(body))
    env = {k: v for k, v in os.environ.items()
           if k not in ("OPENQP_ROOT", "PYTEST_ADDOPTS", "PYTHONPATH")}
    # importlib mode: the child's own test directory is outside tests/, so the
    # default prepend mode would put it on sys.path and trip the guard itself.
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "-p", "no:xdist",
         "--import-mode=importlib", "--rootdir", str(tmp_path),
         "--confcutdir", str(tmp_path), "-rA", "-q", str(tmp_path)],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=300,
    )
    return result.returncode, result.stdout + result.stderr


def test_module_that_edits_sys_path_at_import_fails_collection(tmp_path):
    code, out = _run_child(tmp_path, {"test_child_insert.py": MODULE_LEVEL_INSERT})
    assert code != 0, out
    assert "ERROR test_child_insert.py" in out, out
    assert "importing this module changed process-global import state" in out, out
    assert "sys.path changed" in out, out


def test_tests_that_leak_fail_and_scoped_changes_pass(tmp_path):
    code, out = _run_child(tmp_path, {"test_child_leaks.py": PER_TEST})
    assert code != 0, out
    for name in ("test_leaves_source_tree_on_sys_path", "test_leaves_openqp_root_set",
                 "test_leaves_stub_module"):
        assert f"ERROR test_child_leaks.py::{name}" in out, out
    for reason in ("sys.path changed", "OPENQP_ROOT changed",
                   "oqp.zz_guard_stub (stub: no __file__)"):
        assert reason in out, out
    assert "PASSED test_child_leaks.py::test_scoped_changes_are_allowed" in out, out
    assert "ERROR test_child_leaks.py::test_scoped_changes_are_allowed" not in out, out
