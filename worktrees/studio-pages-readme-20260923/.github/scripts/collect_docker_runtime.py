#!/usr/bin/env python3
"""Stage the non-glibc shared-library closure for the Docker runtime.

The OpenQP wheel contains its own native library and canonical DFT-D4 stack,
but the Docker builder also links to an ILP64 OpenBLAS and GNU runtime
libraries.  This script walks ``ldd`` recursively, copies only that external
closure, and records exact hashes/package licenses.  Any unresolved library,
unowned non-OpenBLAS library, or colliding SONAME aborts the image build.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from collections import deque
from pathlib import Path


# These are supplied by the digest-pinned glibc runtime base and must not be
# copied from Ubuntu into the Debian runtime.  GCC/OpenMP/Fortran libraries are
# deliberately absent from this set and are staged from the builder instead.
GLIBC_RUNTIME_LIBRARIES = {
    "ld-linux-x86-64.so.2",
    "libanl.so.1",
    "libc.so.6",
    "libdl.so.2",
    "libm.so.6",
    "libnss_compat.so.2",
    "libnss_dns.so.2",
    "libnss_files.so.2",
    "libnss_hesiod.so.2",
    "libpthread.so.0",
    "libresolv.so.2",
    "librt.so.1",
    "libthread_db.so.1",
    "libutil.so.1",
}

_ARROW_LINE = re.compile(r"^\s*(\S+)\s+=>\s+(\S+)\s+\(0x[0-9a-fA-F]+\)\s*$")
_DIRECT_LINE = re.compile(r"^\s*(/\S+)\s+\(0x[0-9a-fA-F]+\)\s*$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_ldd(output: str) -> list[tuple[str, Path]]:
    """Parse resolved ``ldd`` entries and reject missing dependencies."""

    resolved: list[tuple[str, Path]] = []
    for raw_line in output.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("linux-vdso.so"):
            continue
        if "=> not found" in line:
            raise RuntimeError(f"unresolved shared-library dependency: {line}")
        arrow = _ARROW_LINE.match(raw_line)
        if arrow:
            soname, value = arrow.groups()
            path = Path(value)
            if not path.is_absolute():
                raise RuntimeError(f"ldd returned a non-absolute dependency: {line}")
            resolved.append((soname, path))
            continue
        direct = _DIRECT_LINE.match(raw_line)
        if direct:
            path = Path(direct.group(1))
            resolved.append((path.name, path))
            continue
        if line in {"statically linked", "not a dynamic executable"}:
            raise RuntimeError(f"expected a dynamic shared object, got: {line}")
        raise RuntimeError(f"unrecognized ldd output: {line}")
    return resolved


def ldd(
    path: Path, library_paths: tuple[Path, ...] = ()
) -> list[tuple[str, Path]]:
    environment = os.environ.copy()
    environment.pop("LD_LIBRARY_PATH", None)
    if library_paths:
        environment["LD_LIBRARY_PATH"] = os.pathsep.join(
            str(directory) for directory in library_paths
        )
    result = subprocess.run(
        ["ldd", str(path)], capture_output=True, text=True, check=False,
        env=environment,
    )
    output = result.stdout + result.stderr
    if result.returncode != 0:
        raise RuntimeError(f"ldd failed for {path}:\n{output}")
    return parse_ldd(output)


def dpkg_owner(path: Path) -> tuple[str, str, str, str, Path]:
    """Return binary/source package metadata and copyright for ``path``."""

    candidates = [path, path.resolve()]
    search_output = ""
    for candidate in candidates:
        result = subprocess.run(
            ["dpkg-query", "--search", str(candidate)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            search_output = result.stdout
            break
    if not search_output:
        raise RuntimeError(f"runtime dependency is not owned by a package: {path}")

    packages = sorted(
        {
            line.split(": ", 1)[0]
            for line in search_output.splitlines()
            if ": " in line
        }
    )
    if not packages:
        raise RuntimeError(f"could not parse dpkg owner for {path}: {search_output}")
    package = packages[0]
    result = subprocess.run(
        [
            "dpkg-query",
            "--show",
            "--showformat=${binary:Package}\t${Version}\t"
            "${source:Package}\t${source:Version}\n",
            package,
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    fields = result.stdout.rstrip("\n").split("\t")
    if len(fields) != 4:
        raise RuntimeError(f"could not parse dpkg metadata for {package}: {fields}")
    binary_package, version, source_package, source_version = fields
    source_package = source_package or binary_package.split(":", 1)[0]
    source_version = source_version or version
    doc_name = binary_package.split(":", 1)[0]
    copyright_file = Path("/usr/share/doc") / doc_name / "copyright"
    if not copyright_file.is_file():
        raise RuntimeError(
            f"package {binary_package} has no installed copyright file: "
            f"{copyright_file}"
        )
    return binary_package, version, source_package, source_version, copyright_file


def openblas_runtime_config(path: Path) -> tuple[str, dict[str, object]]:
    """Read and validate version/ABI features from the actual OpenBLAS binary."""

    library = ctypes.CDLL(str(path.resolve()))
    get_config = library.openblas_get_config
    get_config.argtypes = []
    get_config.restype = ctypes.c_char_p
    raw_config = get_config()
    if not raw_config:
        raise RuntimeError(f"openblas_get_config returned no data for {path}")
    config = raw_config.decode("ascii", errors="strict")
    match = re.search(r"\bOpenBLAS\s+([0-9]+(?:\.[0-9]+)+)\b", config)
    if not match:
        raise RuntimeError(f"could not parse OpenBLAS version from {config!r}")
    required_tokens = {"DYNAMIC_ARCH", "USE64BITINT", "USE_OPENMP"}
    tokens = set(config.split())
    missing_tokens = sorted(required_tokens - tokens)
    if missing_tokens:
        raise RuntimeError(
            f"OpenBLAS runtime config lacks required ILP64/OpenMP features "
            f"{missing_tokens}: {config!r}"
        )
    try:
        getattr(library, "dgemm_")
    except AttributeError as error:
        raise RuntimeError("OpenBLAS lacks the required unsuffixed dgemm_ ABI") from error
    try:
        getattr(library, "dgemm_64_")
    except AttributeError:
        pass
    else:
        raise RuntimeError("OpenBLAS unexpectedly exports suffixed dgemm_64_ ABI")
    return match.group(1), {
        "dynamic_arch": True,
        "interface64": True,
        "openmp": True,
        "symbol_suffix": "none",
    }


def _copy_unique(source: Path, destination: Path, seen: dict[Path, str]) -> str:
    source_hash = sha256(source.resolve())
    previous_hash = seen.get(destination)
    if previous_hash is not None and previous_hash != source_hash:
        raise RuntimeError(
            f"different libraries collide at runtime path {destination}: "
            f"{previous_hash} != {source_hash}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    if previous_hash is None:
        shutil.copy2(source.resolve(), destination)
        seen[destination] = source_hash
    return source_hash


def collect(
    output: Path,
    expected_openblas_version: str,
    expected_openblas_abi: dict[str, object],
    builder_lock: dict[str, object],
) -> dict[str, object]:
    from oqp.runtime import resolve_oqp_root

    oqp_root, _ = resolve_oqp_root()
    package_root = Path(oqp_root).resolve()
    package_lib = package_root / "lib"
    if not package_lib.is_dir():
        raise RuntimeError(f"installed OpenQP library directory missing: {package_lib}")
    openblas_lib = Path("/opt/openblas/lib")
    if not openblas_lib.is_dir():
        raise RuntimeError(f"pinned OpenBLAS library directory missing: {openblas_lib}")
    runtime_search_paths = (package_lib, openblas_lib)

    seeds = sorted(
        {
            path.resolve()
            for path in package_lib.rglob("*")
            if path.is_file() and (".so" in path.name or path.name.endswith(".dylib"))
        }
    )
    if not seeds:
        raise RuntimeError(f"no installed OpenQP shared objects found in {package_lib}")

    expected_python = str(builder_lock["python_version"])
    actual_python = f"{sys.version_info.major}.{sys.version_info.minor}"
    if actual_python != expected_python:
        raise RuntimeError(
            f"pinned builder Python {actual_python} != lock {expected_python}"
        )
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "amd64"}:
        raise RuntimeError(
            f"pinned builder platform is not Linux/amd64: "
            f"{platform.system()}/{platform.machine()}"
        )
    toolchain_lock = dict(builder_lock["toolchain"])
    toolchain_versions = {}
    for role in ("cc", "cxx", "fc"):
        executable = str(toolchain_lock[role])
        version = subprocess.check_output(
            [executable, "-dumpfullversion", "-dumpversion"], text=True
        ).strip()
        if version.split(".", 1)[0] != str(toolchain_lock["major_version"]):
            raise RuntimeError(
                f"pinned builder {executable} version {version} does not match "
                f"lock major {toolchain_lock['major_version']}"
            )
        toolchain_versions[role] = {"executable": executable, "version": version}

    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"runtime staging directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)

    queue: deque[Path] = deque(seeds)
    inspected: set[Path] = set()
    copied: dict[Path, str] = {}
    packages: dict[str, dict[str, str]] = {}
    entries: dict[str, dict[str, object]] = {}
    source_packages: dict[tuple[str, str], None] = {}
    verified_openblas_version: str | None = None
    verified_openblas_abi: dict[str, object] | None = None

    while queue:
        owner = queue.popleft().resolve()
        if owner in inspected:
            continue
        inspected.add(owner)
        for soname, dependency in ldd(owner, runtime_search_paths):
            dependency = dependency.resolve()
            if not dependency.is_file():
                raise RuntimeError(f"ldd dependency does not exist: {dependency}")
            basename = soname.rsplit("/", 1)[-1]
            if basename in GLIBC_RUNTIME_LIBRARIES:
                continue
            if dependency == package_root or package_root in dependency.parents:
                queue.append(dependency)
                continue

            is_openblas = (
                dependency == Path("/opt/openblas")
                or Path("/opt/openblas") in dependency.parents
                or basename.startswith("libopenblas")
            )
            if is_openblas:
                actual_version, actual_abi = openblas_runtime_config(dependency)
                if actual_version != expected_openblas_version:
                    raise RuntimeError(
                        "pinned builder OpenBLAS version mismatch: "
                        f"binary reports {actual_version}, lock expects "
                        f"{expected_openblas_version}"
                    )
                if actual_abi != expected_openblas_abi:
                    raise RuntimeError(
                        "pinned builder OpenBLAS ABI mismatch: "
                        f"binary reports {actual_abi}, lock expects "
                        f"{expected_openblas_abi}"
                    )
                verified_openblas_version = actual_version
                verified_openblas_abi = actual_abi
                destination = output / "openblas" / basename
                package_name = "OpenBLAS"
                package_version = actual_version
                license_name = f"openblas-{actual_version}-bsd-3-clause.txt"
                license_source = (
                    Path.cwd()
                    / "licenses"
                    / "third_party"
                    / "openblas-bsd-3-clause.txt"
                )
                if not license_source.is_file():
                    raise RuntimeError(f"OpenBLAS license missing: {license_source}")
                license_destination = output / "licenses" / license_name
                _copy_unique(license_source, license_destination, copied)
            else:
                destination = output / "lib" / basename
                (
                    package_name,
                    package_version,
                    source_package,
                    source_version,
                    license_source,
                ) = dpkg_owner(dependency)
                source_packages[(source_package, source_version)] = None
                safe_package = package_name.replace(":", "-")
                license_name = f"{safe_package}-{package_version}-copyright"
                license_destination = output / "licenses" / license_name
                _copy_unique(license_source, license_destination, copied)

            library_hash = _copy_unique(dependency, destination, copied)
            queue.append(dependency)
            packages[package_name] = {
                "version": package_version,
                "license_file": f"runtime-packages/{license_name}",
            }
            entry = entries.setdefault(
                str(destination.relative_to(output)),
                {
                    "sha256": library_hash,
                    "builder_soname": basename,
                    "package": package_name,
                    "package_version": package_version,
                    "required_by": [],
                },
            )
            entry["required_by"].append(owner.name)

    if not any(name.startswith("openblas/") for name in entries):
        raise RuntimeError("OpenQP runtime closure did not include OpenBLAS")
    if verified_openblas_version is None:
        raise RuntimeError("OpenBLAS binary version was not verified")
    if verified_openblas_abi is None:
        raise RuntimeError("OpenBLAS binary ABI was not verified")

    for entry in entries.values():
        entry["required_by"] = sorted(set(entry["required_by"]))
    manifest: dict[str, object] = {
        "schema_version": 1,
        "collector": ".github/scripts/collect_docker_runtime.py",
        "builder_verified": {
            "platform": "linux/amd64",
            "python_version": actual_python,
            "toolchain": toolchain_versions,
        },
        "openblas_version_verified_from_binary": verified_openblas_version,
        "openblas_abi_verified_from_binary": verified_openblas_abi,
        "libraries": dict(sorted(entries.items())),
        "packages": dict(sorted(packages.items())),
        "publication_gate": {
            "ready": False,
            "reason": (
                "Matching source packages for copied builder runtime libraries "
                "are recorded but not bundled in this candidate image."
            ),
            "required_source_packages": [
                {"name": name, "version": version}
                for name, version in sorted(source_packages)
            ],
        },
    }
    manifest_path = output / "runtime-library-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--base-lock", type=Path, required=True)
    args = parser.parse_args()
    lock = json.loads(args.base_lock.read_text())
    builder = lock["images"]["builder"]
    manifest = collect(
        args.output,
        str(builder["openblas_version"]),
        dict(builder["openblas_abi"]),
        builder,
    )
    print(
        "Staged Docker runtime closure: "
        f"{len(manifest['libraries'])} libraries from "
        f"{len(manifest['packages'])} packages"
    )


if __name__ == "__main__":
    main()
