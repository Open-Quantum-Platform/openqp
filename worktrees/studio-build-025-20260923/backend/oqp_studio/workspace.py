"""Where the app keeps the results it produces.

Every run writes a log, a JSON export and often a molden file into a
directory of its own. Those are the user's data, not the app's, so where they
land is the user's decision: this module stores that choice and resolves it,
falling back to a location the app can be sure of when it has not been made.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

from . import network


def settings_path() -> Path:
    """Beside the app's other settings, not beside the results themselves."""
    return network.settings_path().parent / "workspace.json"


def configured() -> str:
    """The directory the user chose, or an empty string."""
    try:
        stored = json.loads(settings_path().read_text())
    except (OSError, ValueError):
        return ""
    return str(stored.get("jobs_dir") or "") if isinstance(stored, dict) else ""


def save(jobs_dir: str) -> str:
    """Record the user's choice. Returns the path as stored."""
    path = settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"jobs_dir": jobs_dir}, indent=2))
    return jobs_dir


def usable(directory: Path) -> Path:
    """The directory, created and proven writable.

    Creating it is not enough to know results can be written there: a network
    share can be mounted read-only, and macOS can hand back a directory it
    will then refuse writes to. So a file is actually written.
    """
    directory.mkdir(parents=True, exist_ok=True)
    probe = directory / ".oqp-studio-write-test"
    probe.write_text("")
    probe.unlink()
    return directory


def default_candidates() -> list[Path]:
    """Where results go when the user has not said, best first.

    Documents, because the outputs are the user's own data: they will want to
    open them, copy them to a cluster, back them up. Application Support is
    where the app's settings and the engine belong, but Finder hides
    ~/Library, so results put there are results the user cannot find.

    It is a list because the first choice can fail -- macOS asks permission
    before an app may write to Documents, and that can be refused.
    """
    home = Path.home()
    roots = []
    documents = home / "Documents"
    if documents.is_dir():
        roots.append(documents / "OQP Studio" / "jobs")
    roots.append(home / "OQP Studio" / "jobs")
    roots.append(network.settings_path().parent / "jobs")
    roots.append(Path(tempfile.gettempdir()) / "oqp-studio-jobs")
    return roots


def preferred() -> Path:
    """Where results should go, worked out without touching the disk.

    Nothing here may create a directory or write a file. macOS asks the user
    before an app may write to Documents and BLOCKS until they answer, and
    this runs while the module is being imported -- ahead of the server
    binding its port, which the shell gives up waiting for after thirty
    seconds. Startup must not depend on a filesystem the app may not have
    permission for yet; `ensure()` does that part, on first use.
    """
    override = os.environ.get("OQP_STUDIO_JOBS")
    if override:
        return Path(override).expanduser()
    chosen = configured()
    if chosen:
        return Path(chosen).expanduser()
    return Path.home() / "Documents" / "OQP Studio" / "jobs"


def ensure(start: Path | None = None) -> Path:
    """The first directory this process can actually write results into.

    Called when a run is submitted, not at import: on macOS the consent
    prompt then arrives in response to something the user did, with the
    window already up.

    Never the working directory: an app launched from Finder or the Dock
    inherits launchd's cwd, which is "/", so a relative default resolved to
    /jobs_data -- a path macOS refuses an unprivileged process, and every run
    failed at submission with nothing but a 500.
    """
    wanted = start or preferred()
    for candidate in [wanted, *default_candidates()]:
        try:
            return usable(candidate).resolve()
        except OSError:
            continue
    return usable(Path(tempfile.gettempdir()) / "oqp-studio-jobs").resolve()


def resolve() -> Path:
    """Backwards-compatible name for `ensure()`."""
    return ensure()


def status(active: Path) -> dict:
    """What the settings dialog shows."""
    return {
        "jobs_dir": configured(),
        "active": str(active),
        "default": str(default_candidates()[0]),
        "overridden": bool(os.environ.get("OQP_STUDIO_JOBS")),
        "settings_path": str(settings_path()),
    }
