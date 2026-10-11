"""isolate_modules() and stub_module() must leave real packages as they found them."""
import re
import sys
import types
from pathlib import Path

from oqp_module_isolation import isolate_modules, stub_module


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
