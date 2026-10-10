"""Install and locate the standalone OpenQP compute engine.

The release carries an archive of the engine for each platform — unzip and
run, with no Python, compiler or BLAS on the user's machine. This module
downloads the one that matches, unpacks it beside the app's settings, and
tells the runners where it went, so a fresh install can compute without the
user assembling a toolchain first.

The same archive is a standalone command-line program: anyone who prefers a
terminal can download it from the release page and run ./openqp directly,
which is why it is published as its own asset rather than hidden inside the
installer.
"""

from __future__ import annotations

import os
import platform
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
import zipfile
from pathlib import Path

from . import network
from .archive import extract_tar, extract_zip, validate_links

EXECUTABLE = "openqp.exe" if os.name == "nt" else "openqp"
_VERSION_LINE = re.compile(r"^OpenQP version\s*:\s*(.+)$", re.MULTILINE)
_versions: dict[Path, str | None] = {}


def data_dir() -> Path:
    """Where the engine is unpacked: beside the app's own settings."""
    return network.settings_path().parent / "engine"


def archive_suffix() -> str | None:
    """The release asset that matches this machine."""
    machine = platform.machine().lower()
    if sys.platform == "darwin":
        return "macos-arm64.zip" if machine in ("arm64", "aarch64") else "macos-x86_64.zip"
    if os.name == "nt":
        return "windows-x86_64.zip"
    if sys.platform.startswith("linux"):
        return "linux-aarch64.tar.gz" if machine in ("aarch64", "arm64") else "linux-x86_64.tar.gz"
    return None


def pick_asset(assets: list[dict]) -> dict | None:
    suffix = archive_suffix()
    if not suffix:
        return None
    for asset in assets:
        name = str(asset.get("name", ""))
        if name.startswith("openqp-") and name.endswith(suffix):
            return asset
    return None


_bundled: Path | None = None
_bundled_resolved = False


def bundled_dir() -> Path | None:
    """The engine an all-in-one installer shipped inside the app, if any.

    The shell passes its resource directory in OQP_STUDIO_RESOURCES, because
    that directory is somewhere different on every platform -- beside the
    executable on Windows, Contents/Resources on macOS, /usr/lib elsewhere.
    Falling back to the sidecar's own location covers a build that predates
    the variable, and running the backend straight out of a source tree.
    """
    global _bundled, _bundled_resolved

    if _bundled_resolved:
        return _bundled
    roots: list[Path] = []
    declared = os.environ.get("OQP_STUDIO_RESOURCES")
    if declared:
        roots.append(Path(declared))
    executable = Path(sys.executable).resolve().parent
    roots += [executable, executable.parent / "Resources"]
    for root in roots:
        candidate = root / "engine"
        if (candidate / EXECUTABLE).is_file():
            _bundled = candidate
            break
    _bundled_resolved = True
    if _bundled is not None:
        _make_runnable(_bundled)
    return _bundled


def _make_runnable(directory: Path) -> None:
    """Clear what would stop the bundled engine from being executed.

    An unsigned app downloaded through a browser arrives quarantined, and
    Gatekeeper applies that to every executable inside it -- so the app can be
    approved and opened while the engine it carries is still refused when the
    backend tries to run it. Clearing the flag on our own directory is the
    same thing the downloaded engine already does, done once per process.
    """
    if sys.platform != "darwin":
        return
    import subprocess

    try:
        subprocess.run(["xattr", "-dr", "com.apple.quarantine", str(directory)],
                       check=False, capture_output=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        # Best effort: a quarantined engine still runs once the user approves
        # it, and there is nothing useful to do here if xattr is unavailable.
        pass


# Where an engine may already be, in the order worth trying: one the user
# downloaded through the app (deliberately chosen, so it wins), the copy an
# all-in-one installer shipped, then the places an archive gets unpacked by
# hand.
def search_paths() -> list[Path]:
    candidates = [data_dir()]
    bundled = bundled_dir()
    if bundled is not None:
        candidates.append(bundled)
    home = Path.home()
    candidates += [
        home / "openqp",
        home / "Applications" / "openqp",
        home / ".local" / "share" / "openqp",
        Path("/opt/openqp"),
        Path("/usr/local/openqp"),
    ]
    return candidates


def locate() -> str | None:
    """An engine executable this machine can run, or None.

    PATH wins: a user who installed OpenQP themselves means that one.
    """
    from . import environment

    on_path = environment.locate(os.environ.get("OQP_COMMAND", "openqp"))
    if on_path:
        return on_path
    for directory in search_paths():
        candidate = directory / EXECUTABLE
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


def bundled_or_downloaded() -> str | None:
    """The engine managed by Studio, deliberately excluding PATH.

    This keeps the bundled-runner choice meaningful when a different OpenQP
    command is also installed on the computer.
    """
    for directory in (bundled_dir(), data_dir()):
        if directory is None:
            continue
        candidate = directory / EXECUTABLE
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


def version(command: str | None) -> str | None:
    """Best-effort OpenQP version for a specific executable.

    The standalone engine records its version in README.txt. A pip-installed
    command is a Python script, so its interpreter can report the OpenQP
    distribution metadata. The executable has no --version flag.
    """
    if not command:
        return None
    path = Path(command).resolve()
    if path in _versions:
        return _versions[path]
    found: str | None = None
    readme = path.parent / "README.txt"
    try:
        match = _VERSION_LINE.search(readme.read_text(errors="replace"))
        if match:
            found = match.group(1).strip()
    except OSError:
        pass
    if found is None:
        try:
            first_line = path.read_text(errors="replace").splitlines()[0]
            if first_line.startswith("#!"):
                interpreter = first_line[2:].strip().split()[0]
                result = subprocess.run(
                    [interpreter, "-c", "from importlib.metadata import version; print(version('OpenQP'))"],
                    capture_output=True,
                    text=True,
                    timeout=3,
                    check=False,
                )
                if result.returncode == 0 and result.stdout.strip():
                    found = result.stdout.strip().splitlines()[0]
        except (OSError, IndexError, subprocess.SubprocessError):
            pass
    _versions[path] = found
    return found


def _source(found: str | None) -> str:
    """Which of the three ways this engine got here, for the UI to show."""
    if not found:
        return ""
    path = Path(found).resolve()
    if path.parent == data_dir().resolve():
        return "downloaded by this app"
    bundled = bundled_dir()
    if bundled is not None and path.parent == bundled.resolve():
        return "included with the installer"
    return "found on this machine"


def status() -> dict:
    found = locate()
    return {
        "installed": bool(found),
        "path": found,
        "source": _source(found),
        "install_dir": str(data_dir()),
        "archive": archive_suffix(),
        "supported": archive_suffix() is not None,
    }


def _strip_single_root(target: Path) -> None:
    """Flatten an archive that unpacked into one wrapper directory.

    The wrapper is called openqp and so is the executable inside it, so the
    contents cannot be moved up in place — the first move would land on the
    wrapper itself. The wrapper is moved aside first.
    """
    entries = [path for path in target.iterdir() if path.name != "__MACOSX"]
    if len(entries) != 1 or entries[0].is_symlink() or not entries[0].is_dir():
        return
    staging = target.parent / f"{target.name}.unpack"
    if staging.exists():
        shutil.rmtree(staging)
    entries[0].rename(staging)
    for path in list(staging.iterdir()):
        shutil.move(str(path), str(target / path.name))
    staging.rmdir()


def install(url: str, progress=None) -> str:
    """Download and atomically install an engine archive.

    The current engine remains available until the replacement has downloaded,
    unpacked safely, and passed its executable check.
    """
    target = data_dir()
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{target.name}-", dir=target.parent))
    backup = target.with_name(f".{target.name}-previous")

    suffix = ".zip" if url.endswith(".zip") else ".tar.gz"
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as handle:
        archive = Path(handle.name)
    try:
        with urllib.request.urlopen(url, timeout=300, context=network.context()) as response:
            total = int(response.headers.get("Content-Length") or 0)
            done = 0
            with archive.open("wb") as out:
                while True:
                    chunk = response.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
                    done += len(chunk)
                    if progress:
                        progress(done, total)

        if suffix == ".zip":
            with zipfile.ZipFile(archive) as zipped:
                extract_zip(zipped, staging)
        else:
            with tarfile.open(archive) as tarred:
                extract_tar(tarred, staging)
        _strip_single_root(staging)
        validate_links(staging)

        executable = staging / EXECUTABLE
        if not executable.is_file():
            raise RuntimeError(f"the archive contained no {EXECUTABLE}")
        # Zip archives do not carry the executable bit on every platform.
        executable.chmod(executable.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
        if sys.platform == "darwin":
            # Downloaded archives are quarantined; the engine is unsigned.
            subprocess.run(
                ["xattr", "-dr", "com.apple.quarantine", str(staging)],
                check=False,
                capture_output=True,
            )

        if backup.exists():
            shutil.rmtree(backup)
        if target.exists():
            target.rename(backup)
        try:
            staging.rename(target)
        except Exception:
            if backup.exists() and not target.exists():
                backup.rename(target)
            raise
        if backup.exists():
            shutil.rmtree(backup)
        return str(target / EXECUTABLE)
    finally:
        archive.unlink(missing_ok=True)
        if staging.exists():
            shutil.rmtree(staging)
