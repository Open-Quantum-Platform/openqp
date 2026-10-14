"""Scope ``sys.modules`` stubs to the test that installs them.

Many tests load one source file from ``pyoqp/oqp`` with its ``oqp.*``
dependencies replaced by ``types.ModuleType`` stubs.  ``sys.modules`` is
process-wide, so a stub that is not removed again is what every later test
imports instead of the real package (see the guard in ``conftest.py``).

Call :func:`isolate_modules` before installing stubs and register what it
returns as cleanup, so the stubs stay active for the whole test -- including
imports the loaded module makes lazily at call time -- and are gone after it::

    def setUp(self):
        self.addCleanup(isolate_modules())
        sys.modules["oqp.utils.mpi_utils"] = mpi_utils_stub
"""
import importlib.machinery
import os
import sys
import types


def _session_oqp_spec():
    """Spec of the ``oqp`` package this session uses.

    The loaded package wins.  If none is loaded (or a test's stub stands in
    for it), resolve ``oqp`` through every meta-path finder, as ``import``
    would: an editable install is found by its own finder, not by a
    ``sys.path`` directory, so a ``sys.path``-only search can pick another
    checkout instead.
    """
    module = sys.modules.get("oqp")
    if module is not None and getattr(module, "__path__", None) \
            and getattr(module, "__spec__", None) is not None:
        return module.__spec__
    for finder in sys.meta_path:
        find_spec = getattr(finder, "find_spec", None)
        spec = find_spec("oqp", None) if find_spec is not None else None
        if spec is not None:
            return spec
    return None


def child_env(env=None):
    """Environment for a child interpreter that imports the same ``oqp`` as this one.

    Prepending ``ROOT/pyoqp`` to the child's ``PYTHONPATH`` made it import the
    unbuilt source tree (no liboqp) instead of the installed package; CI hid
    that only because it exported ``OPENQP_ROOT``.  The child is this
    interpreter, so it already finds an installed or editable ``oqp`` through
    the same finders; only a plain package directory found on a ``sys.path``
    entry (conftest's source-tree fallback) needs that entry put first.
    """
    env = dict(os.environ if env is None else env)
    spec = _session_oqp_spec()
    if spec is not None and spec.submodule_search_locations \
            and isinstance(spec.loader, importlib.machinery.SourceFileLoader):
        location = os.path.dirname(list(spec.submodule_search_locations)[0])
        env["PYTHONPATH"] = os.pathsep.join(
            p for p in (location, env.get("PYTHONPATH", "")) if p)
    return env


def stub_module(name):
    """Return the stub installed under ``name``, installing a fresh one if needed.

    Use this instead of ``sys.modules.setdefault(name, ModuleType(name))``:
    when the real package is already imported, ``setdefault`` returns it, and
    the attributes a test then sets (``oqp.ffi``, ``utils.mpi_utils``,
    ``__path__``) rewrite the real package for the rest of the session --
    restoring ``sys.modules`` cannot undo that.
    """
    module = sys.modules.get(name)
    if module is None or getattr(module, "__file__", None) is not None:
        module = types.ModuleType(name)
        sys.modules[name] = module
    return module


def _in_namespace(name, prefixes):
    return any(name == prefix or name.startswith(prefix + ".") for prefix in prefixes)


def isolate_modules(*extra_prefixes):
    """Snapshot the ``oqp`` namespace (and ``extra_prefixes``); return its restorer.

    The snapshot covers whole namespaces rather than a list of stubbed names,
    because such a list drifts from the stubs it is meant to undo.
    """
    prefixes = ("oqp",) + extra_prefixes
    snapshot = {name: module for name, module in sys.modules.items()
                if _in_namespace(name, prefixes)}

    def restore():
        for name in [name for name in sys.modules if _in_namespace(name, prefixes)]:
            if name not in snapshot:
                module = sys.modules.pop(name)
                # Importing a submodule also binds it on its parent package.
                # Leave that binding and the next import creates a second
                # copy while the parent still hands out this one.
                parent, _, child = name.rpartition(".")
                if parent and getattr(sys.modules.get(parent), child, None) is module:
                    delattr(sys.modules[parent], child)
        sys.modules.update(snapshot)

    return restore
