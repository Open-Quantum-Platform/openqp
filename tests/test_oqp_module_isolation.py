"""isolate_modules() and stub_module() must leave real packages as they found them."""
import sys
import types

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
