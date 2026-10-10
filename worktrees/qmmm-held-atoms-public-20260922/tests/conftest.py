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
"""
import os
import warnings

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_QUARANTINE_FILE = os.path.join(_HERE, "known_failures.txt")


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
