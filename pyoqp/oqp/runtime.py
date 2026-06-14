"""Helpers for locating OpenQP runtime files."""

import os
import platform
from pathlib import Path


def library_suffix():
    if platform.uname()[0] == "Windows":
        return "dll"
    if platform.uname()[0] == "Darwin":
        return "dylib"
    return "so"


def _is_oqp_root(path, suffix):
    if not path:
        return False

    root = Path(path)
    return (
        (root / "include" / "oqp.h").exists()
        and (root / "lib" / f"liboqp.{suffix}").exists()
    )


def _source_root_from_package(package_root):
    package_root = Path(package_root)
    if package_root.parent.name == "pyoqp":
        return package_root.parent.parent
    return None


def _candidate_roots(package_root=None):
    package_root = Path(package_root or Path(__file__).resolve().parent)
    roots = [package_root]

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
            f"lib/liboqp.{suffix}: {env_root}"
        )

    searched = ", ".join(str(root) for root in _candidate_roots(package_root))
    raise RuntimeError(
        "Cannot locate OpenQP runtime files. Install OpenQP as a package, run "
        "from a built source tree, or set OPENQP_ROOT to a tree containing "
        f"include/oqp.h and lib/liboqp.{suffix}. Searched: {searched}"
    )


def basis_search_paths(root=None, package_root=None):
    if root is None:
        root, _ = resolve_oqp_root(package_root=package_root)
    root = Path(root)
    return [
        str(root / "share" / "basis_sets"),
        str(root / "basis_sets"),
    ]
