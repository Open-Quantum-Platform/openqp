"""isolate_modules() and stub_module() must leave real packages as they found them."""
import importlib.machinery
import importlib.util
import os
import re
import subprocess
import sys
import types
from pathlib import Path

from oqp_module_isolation import child_env, isolate_modules, stub_module


def _package(root, name):
    (root / name).mkdir()
    (root / name / "__init__.py").write_text("")
    (root / name / "sub.py").write_text("VALUE = 1\n")


def test_submodule_first_imported_while_isolated_is_unbound_from_its_parent(
        tmp_path, monkeypatch):
    _package(tmp_path, "zz_isolation_pkg")
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.delitem(sys.modules, "zz_isolation_pkg", raising=False)
    import zz_isolation_pkg

    restore = isolate_modules("zz_isolation_pkg")
    import zz_isolation_pkg.sub  # noqa: F401
    restore()

    assert "zz_isolation_pkg.sub" not in sys.modules
    assert not hasattr(zz_isolation_pkg, "sub")
    # A later import gives one consistent copy again.
    import zz_isolation_pkg.sub as again
    assert zz_isolation_pkg.sub is again is sys.modules["zz_isolation_pkg.sub"]
    for name in ("zz_isolation_pkg.sub", "zz_isolation_pkg"):
        sys.modules.pop(name, None)


def test_stub_module_never_returns_a_real_module(monkeypatch):
    real = types.ModuleType("zz_real")
    real.__file__ = "/somewhere/zz_real.py"
    monkeypatch.setitem(sys.modules, "zz_real", real)

    stub = stub_module("zz_real")
    assert stub is not real
    assert sys.modules["zz_real"] is stub
    assert stub_module("zz_real") is stub


# A child interpreter must import the oqp this session uses (child_env), never
# the unbuilt source tree, and a test must not point the session itself at it.
_SOURCE_TREE_PATTERNS = tuple(re.compile(p, re.IGNORECASE) for p in (
    r"""PYTHONPATH.*ROOT\s*/\s*['"]pyoqp['"]""",
    r"""pythonpath\s*=\s*str\(\s*ROOT\s*/\s*['"]pyoqp['"]\s*\)""",
    # The bare source package directory as a string has no use but a search path.
    r"""str\(\s*ROOT\s*/\s*['"]pyoqp['"]\s*\)""",
    r"""sys\.path\.(insert|append)\(.*pyoqp""",
    r"""environ\.setdefault\(\s*['"]OPENQP_ROOT""",
))

_SOURCE_TREE_EXCLUDED = {
    # Holds the patterns themselves.
    "test_oqp_module_isolation.py",
    # Deliberately starts a child with the source tree on sys.path to trip the conftest guard.
    "test_import_state_guard.py",
}


def source_tree_path_offences(text):
    """Return (line number, line) for every line that puts the source oqp on a search path."""
    return [(number, line.strip())
            for number, line in enumerate(text.splitlines(), 1)
            if any(pattern.search(line) for pattern in _SOURCE_TREE_PATTERNS)]


def test_no_test_puts_the_source_tree_on_a_search_path():
    offences = []
    for path in sorted(Path(__file__).resolve().parent.glob("test_*.py")):
        if path.name in _SOURCE_TREE_EXCLUDED:
            continue
        for number, line in source_tree_path_offences(path.read_text(encoding="utf-8")):
            offences.append(f"{path.name}:{number}: {line}")
    assert not offences, (
        "use oqp_module_isolation.child_env() for a child's environment instead of "
        "the source pyoqp directory:\n" + "\n".join(offences))


def test_no_example_rewrites_a_childs_search_path():
    # Examples are installed with the package and run against it; a child they
    # start is the same interpreter and finds the same oqp without help.
    examples = Path(__file__).resolve().parents[1] / "examples"
    pattern = re.compile(r"""['"]PYTHONPATH['"]\s*\]\s*=|sys\.path\.(insert|append)\(""")
    offences = [f"{path.relative_to(examples)}:{number}: {line.strip()}"
                for path in sorted(examples.rglob("*.py"))
                for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
                if pattern.search(line)]
    assert not offences, "examples must not edit a search path:\n" + "\n".join(offences)



def _oqp_checkout(root, label):
    package = root / label / "oqp"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(f"WHERE = {label!r}\n")
    return package


def test_child_env_follows_the_loaded_oqp_not_another_checkout_on_sys_path(
        tmp_path, monkeypatch):
    loaded = _oqp_checkout(tmp_path, "loaded")
    other = _oqp_checkout(tmp_path, "other")
    spec = importlib.util.spec_from_file_location(
        "oqp", loaded / "__init__.py", submodule_search_locations=[str(loaded)])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setitem(sys.modules, "oqp", module)
    monkeypatch.syspath_prepend(str(other.parent))  # another checkout, searched first

    env = child_env({"PATH": os.environ.get("PATH", "")})
    child = subprocess.run([sys.executable, "-c", "import oqp; print(oqp.WHERE)"],
                           env=env, cwd=tmp_path, capture_output=True, text=True,
                           check=True)
    assert child.stdout.strip() == "loaded"


class _EditableFinder:
    """Stands in for the meta-path finder an editable install registers."""

    def __init__(self, location):
        self.location = location

    def find_spec(self, name, path=None, target=None):
        if name != "oqp":
            return None
        spec = importlib.machinery.ModuleSpec(name, loader=None, is_package=True)
        spec.submodule_search_locations = [str(self.location)]
        return spec


def test_child_env_leaves_an_editable_install_to_its_finder(tmp_path, monkeypatch):
    other = _oqp_checkout(tmp_path, "other")
    editable = tmp_path / "editable" / "oqp"
    monkeypatch.delitem(sys.modules, "oqp", raising=False)
    monkeypatch.syspath_prepend(str(other.parent))
    monkeypatch.setattr(sys, "meta_path", [_EditableFinder(editable), *sys.meta_path])

    env = child_env({"PYTHONPATH": "/kept"})
    assert env["PYTHONPATH"] == "/kept"
