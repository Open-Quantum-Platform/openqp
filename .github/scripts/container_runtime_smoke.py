#!/usr/bin/env python3
"""Fail-closed verification of the final OpenQP Docker runtime stage."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path


ENERGY_REFERENCE = -1.1167593075
D4_ENERGY_REFERENCE = -0.00019542038860006095
D4_CHILD_MARKER = "OQP_D4_CONTAINER_RESULT="
MRSF_GRADIENT_TOLERANCE = 1.0e-5
PT_LOAD = 1
PT_DYNAMIC = 2
DT_NULL = 0
DT_NEEDED = 1
DT_STRTAB = 5
DT_STRSZ = 10


def require_mrsf_minres_result(output: str, reference: list[float]) -> None:
    """Require a converged, finite MRSF MINRES reference gradient."""
    if not re.search(r"\bZ-Vector converged\b", output):
        raise AssertionError("MRSF MINRES z-vector did not converge")
    marker = "PyOQP dispersion corrected gradients"
    if marker not in output:
        raise AssertionError("MRSF MINRES smoke produced no final gradient")
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
    rows = re.findall(
        rf"^\s*[A-Za-z]{{1,3}}\s+({number})\s+({number})\s+({number})\s*$",
        output.rsplit(marker, 1)[1],
        flags=re.MULTILINE,
    )
    observed = [float(value) for row in rows[:3] for value in row]
    if len(reference) != 9 or len(observed) != len(reference):
        raise AssertionError("MRSF MINRES smoke produced an incomplete gradient")
    if not all(math.isfinite(value) for value in observed):
        raise AssertionError("MRSF MINRES smoke produced a non-finite gradient")
    difference = max(abs(actual - expected) for actual, expected in zip(observed, reference))
    if difference >= MRSF_GRADIENT_TOLERANCE:
        raise AssertionError(
            "MRSF MINRES gradient does not match the packaged reference: "
            f"max difference {difference:.3e}"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ldd_dependencies(path: Path) -> dict[str, Path]:
    result = subprocess.run(
        ["ldd", str(path)], capture_output=True, text=True, check=False
    )
    output = result.stdout + result.stderr
    if result.returncode != 0 or "not found" in output:
        raise AssertionError(f"unresolved runtime dependency for {path}:\n{output}")
    dependencies: dict[str, Path] = {}
    for line in output.splitlines():
        match = re.match(r"\s*(\S+)\s+=>\s+(/\S+)\s+\(0x[0-9a-fA-F]+\)", line)
        if match:
            name, resolved = match.groups()
            dependencies[name] = Path(resolved).resolve()
    return dependencies


def elf_needed(path: Path) -> list[str]:
    """Read direct DT_NEEDED entries without requiring binutils at runtime."""
    data = path.read_bytes()
    if len(data) < 16 or data[:4] != b"\x7fELF":
        raise AssertionError(f"not an ELF runtime library: {path}")
    if data[5] == 1:
        endian = "<"
    elif data[5] == 2:
        endian = ">"
    else:
        raise AssertionError(f"unsupported ELF byte order for {path}: {data[5]}")

    elf_class = data[4]
    if elf_class == 1:
        header_size, phoff_offset, phentsize_offset, phnum_offset = 52, 28, 42, 44
        address_format, dynamic_format = "I", "II"
        program_header_size = 32
    elif elf_class == 2:
        header_size, phoff_offset, phentsize_offset, phnum_offset = 64, 32, 54, 56
        address_format, dynamic_format = "Q", "QQ"
        program_header_size = 56
    else:
        raise AssertionError(f"unsupported ELF class for {path}: {elf_class}")
    if len(data) < header_size:
        raise AssertionError(f"truncated ELF header: {path}")

    phoff = struct.unpack_from(endian + address_format, data, phoff_offset)[0]
    phentsize = struct.unpack_from(endian + "H", data, phentsize_offset)[0]
    phnum = struct.unpack_from(endian + "H", data, phnum_offset)[0]
    if phnum == 0xFFFF or phentsize < program_header_size:
        raise AssertionError(f"unsupported ELF program-header table: {path}")
    if phoff + phentsize * phnum > len(data):
        raise AssertionError(f"truncated ELF program-header table: {path}")

    loads: list[tuple[int, int, int]] = []
    dynamic: tuple[int, int] | None = None
    for index in range(phnum):
        entry = phoff + index * phentsize
        segment_type = struct.unpack_from(endian + "I", data, entry)[0]
        if elf_class == 1:
            offset, virtual_address, file_size = struct.unpack_from(
                endian + "III", data, entry + 4
            )
            file_size = struct.unpack_from(endian + "I", data, entry + 16)[0]
        else:
            offset, virtual_address = struct.unpack_from(endian + "QQ", data, entry + 8)
            file_size = struct.unpack_from(endian + "Q", data, entry + 32)[0]
        if offset + file_size > len(data):
            raise AssertionError(f"ELF segment extends past end of file: {path}")
        if segment_type == PT_LOAD:
            loads.append((virtual_address, offset, file_size))
        elif segment_type == PT_DYNAMIC:
            dynamic = (offset, file_size)
    if dynamic is None:
        raise AssertionError(f"ELF library has no PT_DYNAMIC segment: {path}")

    dynamic_size = struct.calcsize(endian + dynamic_format)
    string_table_address = None
    string_table_size = None
    needed_offsets: list[int] = []
    dynamic_offset, dynamic_bytes = dynamic
    for entry in range(dynamic_offset, dynamic_offset + dynamic_bytes, dynamic_size):
        if entry + dynamic_size > len(data):
            raise AssertionError(f"truncated ELF dynamic table: {path}")
        tag, value = struct.unpack_from(endian + dynamic_format, data, entry)
        if tag == DT_NULL:
            break
        if tag == DT_NEEDED:
            needed_offsets.append(value)
        elif tag == DT_STRTAB:
            string_table_address = value
        elif tag == DT_STRSZ:
            string_table_size = value
    if string_table_address is None or string_table_size is None:
        raise AssertionError(f"ELF dynamic string table is incomplete: {path}")

    string_table_offset = None
    for virtual_address, offset, file_size in loads:
        if virtual_address <= string_table_address < virtual_address + file_size:
            string_table_offset = offset + string_table_address - virtual_address
            break
    if string_table_offset is None or string_table_offset + string_table_size > len(data):
        raise AssertionError(f"ELF dynamic string table is outside load segments: {path}")

    names = []
    for needed_offset in needed_offsets:
        if needed_offset >= string_table_size:
            raise AssertionError(f"ELF DT_NEEDED offset is out of range: {path}")
        start = string_table_offset + needed_offset
        end = data.find(b"\0", start, string_table_offset + string_table_size)
        if end < 0:
            raise AssertionError(f"unterminated ELF DT_NEEDED name: {path}")
        names.append(data[start:end].decode("utf-8"))
    return names


def assert_stack_edges(
    owner: Path, needed: list[str], expected: dict[str, Path]
) -> None:
    prefixes = ("libdftd4", "libmulticharge", "libmctc-lib")
    actual = {name for name in needed if name.startswith(prefixes)}
    expected_names = set(expected)
    if actual != expected_names:
        raise AssertionError(
            f"noncanonical DFT-D4 dependency graph for {owner}: "
            f"{actual} != {expected_names}"
        )


D4_CHILD_SCRIPT = r'''
import json
import os
import numpy as np
from oqp.library.single_point import dftd4_native_disp

z_values = [8, 1, 1]
xyz_values = np.asarray([
    0.0, 0.0, 0.2225904,
    0.0, 1.427594, -0.890371,
    0.0, -1.427594, -0.890371,
], dtype=np.float64).reshape(3, 3)
energy, _ = dftd4_native_disp(
    z_values, xyz_values, "pbe", False, total_charge=0.0
)
reference = -0.00019542038860006095
if abs(energy - reference) >= 1.0e-12:
    raise AssertionError(f"DFT-D4 energy {energy} != {reference}")

loaded = set()
with open("/proc/self/maps", encoding="utf-8") as maps:
    for line in maps:
        fields = line.rstrip("\n").split(None, 5)
        if len(fields) == 6 and fields[5].startswith("/"):
            path = fields[5].removesuffix(" (deleted)")
            if os.path.basename(path).startswith(
                ("libdftd4", "libmulticharge", "libmctc-lib")
            ):
                loaded.add(os.path.realpath(path))
print("OQP_D4_CONTAINER_RESULT=" + json.dumps({
    "energy": energy,
    "loaded": sorted(loaded),
}, sort_keys=True))
'''


def run_d4_child() -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    for variable in (
        "PYTHONPATH",
        "OPENQP_ROOT",
        "DYLD_LIBRARY_PATH",
        "DYLD_FALLBACK_LIBRARY_PATH",
    ):
        environment.pop(variable, None)
    return subprocess.run(
        [sys.executable, "-c", D4_CHILD_SCRIPT],
        capture_output=True,
        text=True,
        env=environment,
        timeout=120,
    )


def require_d4_child(expected_paths: set[Path]) -> None:
    result = run_d4_child()
    if result.returncode != 0:
        raise AssertionError(
            f"isolated DFT-D4 runtime failed:\n{result.stdout}\n{result.stderr}"
        )
    payloads = [
        json.loads(line[len(D4_CHILD_MARKER) :])
        for line in result.stdout.splitlines()
        if line.startswith(D4_CHILD_MARKER)
    ]
    if len(payloads) != 1:
        raise AssertionError(f"DFT-D4 child emitted no unique result: {result.stdout}")
    payload = payloads[0]
    if abs(float(payload["energy"]) - D4_ENERGY_REFERENCE) >= 1.0e-12:
        raise AssertionError(f"unexpected DFT-D4 energy: {payload}")
    loaded = {Path(path).resolve() for path in payload["loaded"]}
    if loaded != expected_paths:
        raise AssertionError(
            f"DFT-D4 loaded alternate libraries: {loaded} != {expected_paths}"
        )


def main() -> None:
    from oqp.runtime import resolve_oqp_root

    oqp_root, native_suffix = resolve_oqp_root()
    root = Path(oqp_root).resolve()
    libdir = root / "lib"
    d4_paths = {
        "libdftd4.so.3": libdir / "libdftd4.so.3",
        "libmulticharge.so.0": libdir / "libmulticharge.so.0",
        "libmctc-lib.so.0": libdir / "libmctc-lib.so.0",
    }
    liboqp = libdir / f"liboqp.{native_suffix}"
    for path in (liboqp, *d4_paths.values()):
        if not path.is_file() or path.is_symlink():
            raise AssertionError(f"missing regular replaceable runtime library: {path}")

    graph = {
        liboqp: {
            "libdftd4.so.3": d4_paths["libdftd4.so.3"],
            "libmctc-lib.so.0": d4_paths["libmctc-lib.so.0"],
        },
        d4_paths["libdftd4.so.3"]: {
            "libmulticharge.so.0": d4_paths["libmulticharge.so.0"],
            "libmctc-lib.so.0": d4_paths["libmctc-lib.so.0"],
        },
        d4_paths["libmulticharge.so.0"]: {
            "libmctc-lib.so.0": d4_paths["libmctc-lib.so.0"],
        },
        d4_paths["libmctc-lib.so.0"]: {},
    }
    all_dependencies: dict[Path, dict[str, Path]] = {}
    for owner, expected in graph.items():
        dependencies = ldd_dependencies(owner)
        if any("nlopt" in name.lower() for name in dependencies):
            raise AssertionError(f"NLopt dependency leaked into runtime: {owner}")
        all_dependencies[owner] = dependencies
        assert_stack_edges(owner, elf_needed(owner), expected)

    openblas_paths = {
        resolved
        for dependencies in all_dependencies.values()
        for name, resolved in dependencies.items()
        if name.startswith("libopenblas")
    }
    expected_openblas = Path("/opt/openblas/lib/libopenblas.so.0").resolve()
    if openblas_paths != {expected_openblas}:
        raise AssertionError(
            f"OpenQP did not load the staged ILP64 OpenBLAS: "
            f"{openblas_paths} != {{{expected_openblas}}}"
        )

    # Verify every staged external library and its package notice by hash.
    license_root = Path("/usr/share/licenses/openqp")
    manifest_path = license_root / "runtime-library-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("publication_gate", {}).get("ready") is not False:
        raise AssertionError("runtime corresponding-source publication gate is not closed")
    if not manifest["publication_gate"].get("required_source_packages"):
        raise AssertionError("copied system runtimes lack exact source-package records")
    if manifest.get("openblas_version_verified_from_binary") != "0.3.28":
        raise AssertionError("container OpenBLAS binary was not verified as 0.3.28")
    verified_builder = manifest.get("builder_verified", {})
    if verified_builder.get("platform") != "linux/amd64":
        raise AssertionError(f"builder platform was not verified: {verified_builder}")
    if verified_builder.get("python_version") != "3.12":
        raise AssertionError(f"builder Python was not verified: {verified_builder}")
    for compiler in verified_builder.get("toolchain", {}).values():
        if not str(compiler.get("version", "")).startswith("14."):
            raise AssertionError(f"builder GCC 14 toolchain was not verified: {compiler}")
    expected_openblas_abi = {
        "dynamic_arch": True,
        "interface64": True,
        "openmp": True,
        "symbol_suffix": "none",
    }
    if manifest.get("openblas_abi_verified_from_binary") != expected_openblas_abi:
        raise AssertionError(
            "container OpenBLAS ILP64/OpenMP ABI was not verified: "
            f"{manifest.get('openblas_abi_verified_from_binary')}"
        )
    serialized_manifest = json.dumps(manifest, sort_keys=True)
    if "/usr/local/" in serialized_manifest or "/opt/openqp/" in serialized_manifest:
        raise AssertionError("transient builder paths leaked into runtime manifest")
    for relative, record in manifest["libraries"].items():
        if relative.startswith("openblas/"):
            installed = Path("/opt/openblas/lib") / Path(relative).name
        elif relative.startswith("lib/"):
            installed = Path("/opt/openqp-runtime/lib") / Path(relative).name
        else:
            raise AssertionError(f"unexpected runtime manifest path: {relative}")
        if not installed.is_file() or sha256(installed) != record["sha256"]:
            raise AssertionError(f"runtime dependency hash mismatch: {installed}")
        dependencies = ldd_dependencies(installed)
        if any("nlopt" in name.lower() for name in dependencies):
            raise AssertionError(f"NLopt dependency leaked into runtime: {installed}")
    for package in manifest["packages"].values():
        notice = license_root / package["license_file"]
        if not notice.is_file() or not notice.read_bytes():
            raise AssertionError(f"runtime package notice missing: {notice}")

    wheelhouse_manifest = json.loads(
        (license_root / "wheelhouse-manifest.json").read_text()
    )
    if wheelhouse_manifest.get("publication_gate", {}).get("ready") is not False:
        raise AssertionError("unlocked wheelhouse publication gate is not closed")
    build_wheels = wheelhouse_manifest.get("build_wheels", [])
    build_names = {
        re.sub(r"[-_.]+", "-", str(wheel.get("name", ""))).lower()
        for wheel in build_wheels
    }
    required_build_names = {
        "scikit-build-core", "cmake", "ninja", "numpy", "cffi", "setuptools"
    }
    if not required_build_names.issubset(build_names):
        raise AssertionError(
            "build-system wheel inventory is incomplete: "
            f"{sorted(required_build_names - build_names)}"
        )
    openqp_wheels = [
        wheel
        for wheel in wheelhouse_manifest.get("wheels", [])
        if wheel.get("name", "").lower() == "openqp"
    ]
    if len(openqp_wheels) != 1:
        raise AssertionError(f"invalid recorded OpenQP wheel set: {openqp_wheels}")
    for wheel in [*wheelhouse_manifest.get("wheels", []), *build_wheels]:
        if not re.fullmatch(r"[0-9a-f]{64}", str(wheel.get("sha256", ""))):
            raise AssertionError(f"wheel input lacks exact SHA-256: {wheel}")

    for filename in (
        "LICENSE",
        "LICENSING.md",
        "SUSTAINABILITY.md",
        "THIRD_PARTY_NOTICES.md",
        "THIRD_PARTY_RUNTIME.md",
    ):
        path = license_root / filename
        if not path.is_file() or not path.read_bytes():
            raise AssertionError(f"OpenQP legal file missing: {path}")
    third_party_notice = (license_root / "THIRD_PARTY_NOTICES.md").read_text().lower()
    if "| nlopt |" in third_party_notice:
        raise AssertionError("THIRD_PARTY_NOTICES.md still describes NLopt as shipped")
    if "statically linked into `liboqp` as part of the dft-d4 stack" in third_party_notice:
        raise AssertionError("THIRD_PARTY_NOTICES.md still describes static DFT-D4")
    if "dft-d4" not in third_party_notice or "dynamically linked" not in third_party_notice:
        raise AssertionError("THIRD_PARTY_NOTICES.md lacks dynamic DFT-D4 disclosure")

    corresponding_source = root / "share/corresponding-source/dftd4-stack"
    required_source = (
        corresponding_source / "README.md",
        corresponding_source / "BUILD-INFO.json",
        corresponding_source / "apply-patch.cmake",
        corresponding_source / "generate-build-info.cmake",
        corresponding_source / "openqp-external-build.cmake",
        corresponding_source / "patches/mctc-lib-0.4.2-disable-tests.patch",
        corresponding_source / "patches/dftd4-3.7.0-disable-tests-and-mstore.patch",
        corresponding_source / "mctc-lib-0.4.2/LICENSE",
        corresponding_source / "multicharge-0.3.0/LICENSE",
        corresponding_source / "dftd4-3.7.0/COPYING",
        corresponding_source / "dftd4-3.7.0/COPYING.LESSER",
    )
    missing_source = [str(path) for path in required_source if not path.is_file()]
    if missing_source:
        raise AssertionError(f"DFT-D4 corresponding source missing: {missing_source}")

    build_info_path = corresponding_source / "BUILD-INFO.json"
    build_info_text = build_info_path.read_text(encoding="utf-8")
    build_info = json.loads(build_info_text)
    expected_runtime_names = {
        "dftd4": "libdftd4.so.3",
        "multicharge": "libmulticharge.so.0",
        "mctc-lib": "libmctc-lib.so.0",
    }
    if (
        build_info.get("schema") != "org.open-quantum-platform.dftd4-build-info"
        or build_info.get("schema_version") != 1
        or build_info.get("canonical_runtime_names") != expected_runtime_names
    ):
        raise AssertionError("DFT-D4 BUILD-INFO schema/runtime names are invalid")
    expected_revision = os.environ.get("OPENQP_EXPECTED_REVISION", "")
    if expected_revision and not re.fullmatch(r"[0-9a-f]{40}", expected_revision):
        raise AssertionError(
            f"container expected revision is not a full Git SHA: {expected_revision!r}"
        )
    openqp_build = build_info.get("openqp", {})
    if expected_revision:
        revision_matches = (
            openqp_build.get("source_revision") == expected_revision
            and openqp_build.get("source_tree_dirty") is False
        )
    else:
        # A direct local ``docker build .`` has no Git metadata in its context.
        # It may omit the revision, but it must not invent a clean source ID.
        revision_matches = (
            openqp_build.get("source_revision") is None
            and openqp_build.get("source_tree_dirty") is None
        )
    if not revision_matches:
        raise AssertionError(
            "DFT-D4 BUILD-INFO does not identify the exact clean container source: "
            f"{openqp_build}"
        )
    expected_components = {
        "mctc-lib": (
            "0.4.2",
            "https://github.com/grimme-lab/mctc-lib/archive/refs/tags/v0.4.2.tar.gz",
            "c7aa45c0a3e6f96e3316e15fc6cdbe48b15234940d3773927a57bb7bfe9744ac",
            "Apache-2.0",
        ),
        "multicharge": (
            "0.3.0",
            "https://github.com/grimme-lab/multicharge/archive/refs/tags/v0.3.0.tar.gz",
            "2fcc1f80871f404f005e9db458ffaec95bb28a19516a0245278cd3175b63a6b2",
            "Apache-2.0",
        ),
        "dftd4": (
            "3.7.0",
            "https://github.com/dftd4/dftd4/archive/refs/tags/v3.7.0.tar.gz",
            "f00b244759eff2c4f54b80a40673440ce951b6ddfa5eee1f46124297e056f69c",
            "LGPL-3.0-or-later",
        ),
    }
    components = {item["name"]: item for item in build_info.get("components", [])}
    if len(build_info.get("components", [])) != len(expected_components):
        raise AssertionError("DFT-D4 BUILD-INFO has duplicate/extra components")
    for name, expected in expected_components.items():
        component = components.get(name, {})
        actual = tuple(
            component.get(field)
            for field in ("version", "archive_url", "sha256", "license")
        )
        if actual != expected:
            raise AssertionError(f"DFT-D4 BUILD-INFO component is invalid: {name}")
    build = build_info.get("build", {})
    blas = build.get("blas", {})
    if (
        build.get("build_type") != "Release"
        or build.get("build_shared_libs") is not True
        or not isinstance(build.get("openmp"), bool)
        or blas.get("integer_bytes") != 8
        or not blas.get("resolved_blas_libraries")
        or not blas.get("resolved_lapack_libraries")
    ):
        raise AssertionError("DFT-D4 BUILD-INFO resolved build data are incomplete")
    expected_patches = {
        "mctc-lib-0.4.2-disable-tests.patch": "mctc-lib",
        "dftd4-3.7.0-disable-tests-and-mstore.patch": "dftd4",
    }
    patches = {item["file"]: item for item in build_info.get("patches", [])}
    if len(build_info.get("patches", [])) != len(expected_patches):
        raise AssertionError("DFT-D4 BUILD-INFO has duplicate/extra patches")
    for patch_name, component in expected_patches.items():
        patch_path = corresponding_source / "patches" / patch_name
        if patches.get(patch_name) != {
            "component": component,
            "file": patch_name,
            "sha256": sha256(patch_path),
            "strip": 1,
        }:
            raise AssertionError(f"DFT-D4 patch record is invalid: {patch_name}")
    for forbidden in (
        "/private/tmp/",
        "/tmp/",
        "/private/var/folders/",
        "/root/.cache/",
        "/.cache/openqp/",
        "/host-cache/",
        "/Library/Caches/openqp/",
        "/opt/openqp",
    ):
        if forbidden in build_info_text:
            raise AssertionError(f"transient path leaked into DFT-D4 BUILD-INFO: {forbidden}")

    static_archives = []
    for base in (root, Path("/opt/openblas"), Path("/opt/openqp-runtime")):
        static_archives.extend(base.rglob("*.a"))
    if static_archives:
        raise AssertionError(f"static archives leaked into runtime: {static_archives}")
    nlopt_artifacts = []
    for base in (root, Path("/opt/openblas"), Path("/opt/openqp-runtime")):
        nlopt_artifacts.extend(
            path for path in base.rglob("*") if "nlopt" in path.name.lower()
        )
    if nlopt_artifacts:
        raise AssertionError(f"NLopt files leaked into runtime: {nlopt_artifacts}")
    for native in root.rglob("*"):
        if native.is_file() and ".so" in native.name:
            lowered = native.read_bytes().lower()
            if b"nlopt_" in lowered or b"libnlopt" in lowered:
                raise AssertionError(f"NLopt symbols/strings leaked into {native}")
    for forbidden in ("gcc", "g++", "gfortran", "cc", "cmake", "ninja", "make", "git"):
        if shutil.which(forbidden):
            raise AssertionError(f"build tool leaked into runtime image: {forbidden}")

    # End-to-end installed example (no source checkout in the final image).
    source_input = root / "share/examples/QUANTUM/h2.inp"
    if not source_input.is_file():
        raise AssertionError(f"packaged smoke input missing: {source_input}")
    with tempfile.TemporaryDirectory(prefix="openqp-container-smoke-") as temporary:
        temporary_path = Path(temporary)
        input_path = temporary_path / "h2.inp"
        shutil.copy2(source_input, input_path)
        result = subprocess.run(
            ["openqp", input_path.name],
            cwd=temporary_path,
            capture_output=True,
            text=True,
            timeout=900,
        )
        log_path = temporary_path / "h2.log"
        log = log_path.read_text() if log_path.is_file() else ""
        if result.returncode != 0:
            raise AssertionError(
                f"container OpenQP smoke failed:\n{result.stdout}\n"
                f"{result.stderr}\n{log[-4000:]}"
            )
        match = re.search(r"Final RHF energy is\s+(-?\d+\.\d+)", log)
        if not match or abs(float(match.group(1)) - ENERGY_REFERENCE) >= 1.0e-6:
            raise AssertionError(f"unexpected container RHF result: {log[-4000:]}")

    # Exercise the module-procedure callbacks used by the MRSF MINRES z-vector
    # path. A basic RHF or DFT-D4 smoke calculation cannot reach this code.
    mrsf_source = (
        root
        / "share/examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.inp"
    )
    if not mrsf_source.is_file():
        raise AssertionError(f"packaged MRSF smoke input missing: {mrsf_source}")
    mrsf_reference_path = mrsf_source.with_suffix(".json")
    if not mrsf_reference_path.is_file():
        raise AssertionError(
            f"packaged MRSF smoke reference missing: {mrsf_reference_path}"
        )
    mrsf_reference = json.loads(mrsf_reference_path.read_text(encoding="utf-8"))
    reference_gradient = mrsf_reference.get("grad")
    if not isinstance(reference_gradient, list):
        raise AssertionError("packaged MRSF smoke reference has no gradient")
    with tempfile.TemporaryDirectory(prefix="openqp-mrsf-minres-smoke-") as temporary:
        temporary_path = Path(temporary)
        input_path = temporary_path / "mrsf-minres.inp"
        input_text = mrsf_source.read_text(encoding="utf-8")
        section = "[tdhf]\ntype=mrsf\n"
        if input_text.count(section) != 1:
            raise AssertionError("packaged MRSF input has an unexpected [tdhf] section")
        input_path.write_text(
            input_text.replace(section, section + "z_solver=2\n"),
            encoding="utf-8",
        )
        result = subprocess.run(
            ["openqp", input_path.name],
            cwd=temporary_path,
            capture_output=True,
            text=True,
            timeout=900,
        )
        log_path = temporary_path / "mrsf-minres.log"
        log = log_path.read_text() if log_path.is_file() else ""
        if result.returncode != 0:
            raise AssertionError(
                f"container MRSF MINRES smoke failed:\n{result.stdout}\n"
                f"{result.stderr}\n{log[-4000:]}"
            )
        if not re.search(r"Solver method\s+is\s+MINRES", log):
            raise AssertionError(f"MRSF smoke did not select MINRES: {log[-4000:]}")
        if "MINRES total iterations" not in log:
            raise AssertionError(f"MRSF MINRES callback did not run: {log[-4000:]}")
        require_mrsf_minres_result(
            "\n".join((result.stdout, result.stderr, log)),
            reference_gradient,
        )

    expected_d4_paths = {path.resolve() for path in d4_paths.values()}
    require_d4_child(expected_d4_paths)
    for name, path in d4_paths.items():
        hidden = path.with_name(f".{path.name}.container-smoke-hidden")
        if hidden.exists():
            raise AssertionError(f"stale DFT-D4 smoke file: {hidden}")
        path.rename(hidden)
        try:
            result = run_d4_child()
            if result.returncode == 0:
                raise AssertionError(
                    f"DFT-D4 still ran after canonical {name} was removed; "
                    "an alternate library copy is present"
                )
        finally:
            hidden.rename(path)
    require_d4_child(expected_d4_paths)

    print(
        "OpenQP minimal container verified: offline wheel install, numerical "
        "runtime, canonical replaceable DFT-D4, complete shared closure/licenses"
    )


if __name__ == "__main__":
    main()
