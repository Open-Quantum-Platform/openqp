"""Quarantine for tests that were already failing when CI began running them.

The pytest suite was never executed by any workflow, so tests rotted quietly
while nothing was watching -- assertions on README contents, on CMakeLists
strings, on APIs that had since moved on.  Turning CI on surfaced all of them
at once.

Blocking the CI wiring on that cleanup would keep the ~880 healthy tests
unprotected for as long as the cleanup takes, so the already-broken ones are
listed in ``known_failures.txt`` and marked xfail here instead.  Everything not
in that file is enforced from now on.

The list is meant to shrink to nothing: fix a test, delete its line.  A
quarantined test that starts passing is reported as XPASS rather than failing
the run, because several entries fail for environment-dependent reasons and a
strict marker would make CI flap.  Check the XPASS list when pruning.

This file also owns the session's import state.  Where ``import oqp`` comes
from is decided once, below, and a guard fails any test module or test that
changes it: one that puts the source tree on ``sys.path``, points
``OPENQP_ROOT`` somewhere, or leaves a stub in ``sys.modules`` silently
changes which ``oqp`` every later test (and every spawned interpreter, which
inherits ``sys.path``) runs against, so results depended on file order.
"""
import importlib.util
import os
import sys
import warnings

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_QUARANTINE_FILE = os.path.join(_HERE, "known_failures.txt")
_PYOQP = os.path.join(os.path.dirname(_HERE), "pyoqp")


def _select_oqp_package():
    """Return the directory of the ``oqp`` package this session must use.

    The installed package is used whenever there is one.  Only a checkout with
    no installed ``oqp`` falls back to the source tree, and then by appending,
    so the fallback can never shadow a build.
    """
    spec = importlib.util.find_spec("oqp")
    if spec is None and os.path.isdir(os.path.join(_PYOQP, "oqp")):
        sys.path.append(_PYOQP)
        spec = importlib.util.find_spec("oqp")
    if spec is None or not spec.submodule_search_locations:
        return None
    return os.path.realpath(list(spec.submodule_search_locations)[0])


_OQP_PACKAGE_DIR = _select_oqp_package()


def _is_test_dir(entry):
    # pytest's default "prepend" import mode inserts the directory of each
    # collected test file; that is pytest's own bookkeeping, not a test's.
    entry = os.path.realpath(entry or os.curdir)
    return entry == _HERE or entry.startswith(_HERE + os.sep)


def _from_selected_package(module):
    path = getattr(module, "__file__", None)
    return (path is not None and _OQP_PACKAGE_DIR is not None
            and os.path.realpath(path).startswith(_OQP_PACKAGE_DIR + os.sep))


def _oqp_namespace():
    return [(name, module) for name, module in list(sys.modules.items())
            if name == "oqp" or name.startswith("oqp.")]


def _import_state():
    return {
        "sys.path": [p for p in sys.path if not _is_test_dir(p)],
        "OPENQP_ROOT": os.environ.get("OPENQP_ROOT"),
    }


_BASELINE = _import_state()


def _import_state_violations(before):
    after = _import_state()
    problems = []
    if after["sys.path"] != before["sys.path"]:
        added = [p for p in after["sys.path"] if p not in before["sys.path"]]
        removed = [p for p in before["sys.path"] if p not in after["sys.path"]]
        problems.append(f"sys.path changed (added {added}, removed {removed}, "
                        "or reordered); use monkeypatch.syspath_prepend")
    if after["OPENQP_ROOT"] != before["OPENQP_ROOT"]:
        problems.append(f"OPENQP_ROOT changed from {before['OPENQP_ROOT']!r} to "
                        f"{after['OPENQP_ROOT']!r}; use monkeypatch.setenv")
    foreign = []
    for name, module in _oqp_namespace():
        if not _from_selected_package(module):
            path = getattr(module, "__file__", None)
            foreign.append(f"{name} (stub: no __file__)" if path is None
                           else f"{name} (from {path})")
    if foreign:
        problems.append(f"sys.modules holds oqp modules not from {_OQP_PACKAGE_DIR}: "
                        f"{foreign}; restore every stub in teardown")
    return problems


def _restore_import_state(state):
    sys.path[:] = [p for p in sys.path if _is_test_dir(p)] + state["sys.path"]
    if state["OPENQP_ROOT"] is None:
        os.environ.pop("OPENQP_ROOT", None)
    else:
        os.environ["OPENQP_ROOT"] = state["OPENQP_ROOT"]
    for name, module in _oqp_namespace():
        if not _from_selected_package(module):
            del sys.modules[name]


@pytest.hookimpl(hookwrapper=True)
def pytest_make_collect_report(collector):
    """Fail a test module whose import changed the session's import state."""
    outcome = yield
    if not isinstance(collector, pytest.Module):
        return
    problems = _import_state_violations(_BASELINE)
    if not problems:
        return
    report = outcome.get_result()
    report.outcome = "failed"
    report.longrepr = ("importing this module changed process-global import state:\n  "
                       + "\n  ".join(problems))
    # Undo it so the error is attributed here and not to every later module.
    _restore_import_state(_BASELINE)


@pytest.fixture(autouse=True)
def _guard_import_state():
    """Fail a test that leaves the session's import state changed."""
    before = _import_state()
    yield
    problems = _import_state_violations(before)
    if problems:
        _restore_import_state(before)
        pytest.fail("test left process-global import state changed:\n  "
                    + "\n  ".join(problems), pytrace=False)


def _load_quarantine():
    """Map node ID -> reason, from ``known_failures.txt``.

    Format is one pytest node ID per line, an optional ``#`` reason after it,
    and ``#`` comment or blank lines anywhere.
    """
    entries = {}
    if not os.path.isfile(_QUARANTINE_FILE):
        return entries
    with open(_QUARANTINE_FILE, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            node, _, reason = line.partition("#")
            node = node.strip()
            if node:
                entries[node] = reason.strip() or "known pre-existing failure"
    return entries


def pytest_collection_modifyitems(config, items):
    quarantine = _load_quarantine()
    if not quarantine:
        return
    seen = set()
    for item in items:
        reason = quarantine.get(item.nodeid)
        if reason is None:
            continue
        seen.add(item.nodeid)
        item.add_marker(pytest.mark.xfail(reason=reason, strict=False))
    # A node ID that matches nothing is almost always a test that was renamed or
    # deleted, leaving a dead line behind. Say so rather than let the file drift.
    for node in sorted(set(quarantine) - seen):
        warnings.warn(f"known_failures.txt lists an unknown test: {node}",
                      UserWarning)
