"""Helpers for locating OpenQP runtime files."""

import os
import platform
from pathlib import Path


_LIBDIR_METADATA = Path("share") / "openqp" / "openqp-runtime-libdir.txt"


def library_suffix():
    if platform.uname()[0] == "Windows":
        return "dll"
    if platform.uname()[0] == "Darwin":
        return "dylib"
    return "so"


def library_directories(root):
    """Return native-library directories for an OpenQP runtime root."""
    root = Path(root)
    candidates = []

    metadata = root / _LIBDIR_METADATA
    if metadata.is_file():
        configured = metadata.read_text(encoding="utf-8").strip()
        if configured:
            configured_path = Path(configured)
            candidates.append(
                configured_path if configured_path.is_absolute()
                else root / configured_path
            )

    # Source-tree builds and older installs have no metadata.  Keep their
    # canonical layout, then recognize common GNUInstallDirs alternatives.
    candidates.extend((root / "lib", root / "lib64"))
    multiarch_root = root / "lib"
    if multiarch_root.is_dir():
        candidates.extend(path for path in sorted(multiarch_root.iterdir()) if path.is_dir())

    seen = set()
    unique = []
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved not in seen:
            unique.append(resolved)
            seen.add(resolved)
    return unique


def library_path(root, suffix=None):
    """Locate liboqp in the configured or a compatible legacy libdir."""
    suffix = suffix or library_suffix()
    for directory in library_directories(root):
        candidate = directory / f"liboqp.{suffix}"
        if candidate.is_file():
            return candidate
    return None


def _is_oqp_root(path, suffix):
    if not path:
        return False

    root = Path(path)
    return (
        (root / "include" / "oqp.h").exists()
        and library_path(root, suffix) is not None
    )


def _source_root_from_package(package_root):
    package_root = Path(package_root)
    if package_root.parent.name == "pyoqp":
        return package_root.parent.parent
    return None


def _candidate_roots(package_root=None):
    package_root = Path(package_root or Path(__file__).resolve().parent)
    roots = []

    # Prefer the runtime files shipped with the current Python package (and the
    # inferred source tree in a dev checkout) so a stale OPENQP_ROOT pointing at
    # another OpenQP build cannot shadow the matching liboqp/data of THIS
    # package. OPENQP_ROOT is only a compatibility fallback (see resolve_oqp_root).
    roots.append(package_root)

    source_root = _source_root_from_package(package_root)
    if source_root is not None:
        roots.append(source_root)

    env_root = os.environ.get("OPENQP_ROOT")
    if env_root:
        roots.append(Path(env_root).expanduser())

    seen = set()
    unique_roots = []
    for root in roots:
        resolved = root.resolve()
        if resolved not in seen:
            unique_roots.append(resolved)
            seen.add(resolved)
    return unique_roots


def resolve_oqp_root(package_root=None):
    """Choose one root containing the matching header and native library.

    Installed wheels keep the Python package, header, native library, and data
    files together under the package directory. Source-tree development usually
    keeps the Python package under pyoqp/oqp and runtime files under the repo
    root. OPENQP_ROOT is retained only as a compatibility fallback for layouts
    that cannot be inferred from the package location.
    """
    suffix = library_suffix()
    for root in _candidate_roots(package_root):
        if _is_oqp_root(root, suffix):
            return str(root), suffix

    env_root = os.environ.get("OPENQP_ROOT")
    if env_root:
        raise RuntimeError(
            "OPENQP_ROOT does not contain matching include/oqp.h and "
            f"liboqp.{suffix} in its configured library directory: {env_root}"
        )

    searched = ", ".join(str(root) for root in _candidate_roots(package_root))
    raise RuntimeError(
        "Cannot locate OpenQP runtime files. Install OpenQP as a package, run "
        "from a built source tree, or set OPENQP_ROOT to a tree containing "
        f"include/oqp.h and liboqp.{suffix}. Searched: {searched}"
    )


def basis_search_paths(root=None, package_root=None):
    if root is None:
        root, _ = resolve_oqp_root(package_root=package_root)
    root = Path(root)
    return [
        str(root / "share" / "basis_sets"),
        str(root / "basis_sets"),
    ]
