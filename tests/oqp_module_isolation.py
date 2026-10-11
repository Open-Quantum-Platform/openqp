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
import sys
import types


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
                del sys.modules[name]
        sys.modules.update(snapshot)

    return restore
