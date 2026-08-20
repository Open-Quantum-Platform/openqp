#!/usr/bin/env python3
"""Numerical smoke test run against every built wheel (CIBW_TEST_COMMAND).

Goes beyond the old import/file-exists check: runs a real RHF/STO-3G
calculation, asserts the raw basis-metadata integer ABI, and exercises the
replaceable DFT-D4 shared-library stack.  This covers failure classes that an
``import oqp`` check alone cannot catch:
  * a wrong/slow BLAS selection (silent NetLib fallback) -- caught here by
    requiring the calculation to complete with the reference energy;
  * the LP64 integer-ABI misread (oqp_get_basis exporting build-width
    integers read as int64 by Python), which corrupted basis metadata
    (centers [0,1] -> [8589934592, -1]) and crashed the molden writer.
  * missing/relocated DFT-D4 runtime libraries, wrong loader paths, lost
    corresponding source, or a changed legacy D4 numerical ABI.

Usage: wheel_smoke_test.py <project-dir>
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

EREF = -1.1167593075  # RHF/STO-3G h2 at 0.74 A (examples/QUANTUM/h2.inp)

project = sys.argv[1]
src = os.path.join(project, "examples", "QUANTUM", "h2.inp")
tmp = tempfile.mkdtemp(prefix="oqp_wheel_test_")
inp = os.path.join(tmp, "h2.inp")
shutil.copy(src, inp)

# 1. Full end-to-end run through the console entry point.
r = subprocess.run(["openqp", "h2.inp"], cwd=tmp, capture_output=True,
                   text=True, timeout=900)
log_path = os.path.join(tmp, "h2.log")
log = open(log_path).read() if os.path.exists(log_path) else ""
assert r.returncode == 0, (
    f"openqp h2.inp failed rc={r.returncode}\n--- stdout ---\n{r.stdout[-2000:]}"
    f"\n--- stderr ---\n{r.stderr[-2000:]}\n--- log tail ---\n{log[-2000:]}")

m = re.search(r"Final RHF energy is\s+(-?\d+\.\d+)", log)
assert m, f"no 'Final RHF energy is' line in h2.log\n{log[-2000:]}"
energy = float(m.group(1))
assert abs(energy - EREF) < 1e-6, f"energy {energy} != reference {EREF}"

# 2. The molden writer must have produced output (end-to-end guard for the
#    basis-export path; scf.save_molden defaults to True).
assert any(f.endswith(".molden") for f in os.listdir(tmp)), (
    f"no .molden file written; dir: {os.listdir(tmp)}")

# 3. Raw basis-metadata ABI check: the integer arrays crossing the
#    Fortran->C->Python boundary must be exact small integers.
os.chdir(tmp)
from oqp.pyoqp import Runner  # noqa: E402  (import after wheel install)
run = Runner(project="abi_probe", input_file=inp,
             log=os.path.join(tmp, "probe.log"), silent=1)
run.run(test_mod=True)
basis = run.mol.data.get_basis()
got = {k: list(map(int, basis[k])) for k in ("centers", "ncontr", "angs")}
expect = {"centers": [0, 1], "ncontr": [3, 3], "angs": [0, 0]}
assert got == expect, f"basis metadata ABI corrupt: {got} != {expect}"

# 4. DFT-D4 must remain a package-local, replaceable shared-library stack after
#    auditwheel/delocate repair.  Check the files, dependency graph, relative
#    loader paths, and absence of persistent external-cache paths.
from oqp.runtime import resolve_oqp_root  # noqa: E402

oqp_root, native_suffix = resolve_oqp_root()
libdir = Path(oqp_root) / "lib"
if sys.platform == "darwin":
    d4_names = {
        "dftd4": "libdftd4.3.dylib",
        "multicharge": "libmulticharge.0.dylib",
        "mctc": "libmctc-lib.0.dylib",
    }
    inspect_command = lambda path: ["otool", "-L", str(path)]
    rpath_command = lambda path: ["otool", "-l", str(path)]
    local_rpath = "@loader_path"
else:
    d4_names = {
        "dftd4": "libdftd4.so.3",
        "multicharge": "libmulticharge.so.0",
        "mctc": "libmctc-lib.so.0",
    }
    inspect_command = lambda path: ["readelf", "-d", str(path)]
    rpath_command = inspect_command
    local_rpath = "$ORIGIN"

d4_paths = {name: libdir / filename for name, filename in d4_names.items()}
for name, path in d4_paths.items():
    assert path.is_file(), f"missing package-local {name} shared library: {path}"
    assert not path.is_symlink(), f"wheel must contain a regular replaceable file: {path}"
assert not list(libdir.glob("libdftd4.a")), "static libdftd4 archive leaked into wheel"
assert not list(libdir.glob("libmulticharge.a")), "static multicharge archive leaked into wheel"
assert not list(libdir.glob("libmctc-lib.a")), "static mctc-lib archive leaked into wheel"

liboqp_path = libdir / f"liboqp.{native_suffix}"


def native_metadata(path, command):
    result = subprocess.run(command(path), capture_output=True, text=True, check=True)
    metadata = result.stdout + result.stderr
    assert "/root/.cache/openqp" not in metadata, (
        f"build cache path leaked into {path}:\n{metadata}"
    )
    assert "Library/Caches/openqp" not in metadata, (
        f"macOS build cache path leaked into {path}:\n{metadata}"
    )
    return metadata


def parse_dynamic_dependencies(metadata, platform_name):
    """Return exact DT_NEEDED/install-name entries from native-tool output."""
    if platform_name == "darwin":
        dependencies = []
        for line in metadata.splitlines()[1:]:
            match = re.match(
                r"\s*(.+?)\s+\(compatibility version [^)]+\)$", line
            )
            if match:
                dependencies.append(match.group(1))
        return dependencies
    return re.findall(
        r"\(NEEDED\)\s+Shared library:\s*\[([^]]+)\]", metadata
    )


def parse_runtime_search_paths(metadata, platform_name):
    """Return individual LC_RPATH/RPATH/RUNPATH entries."""
    if platform_name == "darwin":
        return re.findall(
            r"^\s*path\s+(\S+)\s+\(offset\s+\d+\)\s*$",
            metadata,
            re.MULTILINE,
        )
    encoded_paths = re.findall(
        r"\((?:RPATH|RUNPATH)\).*?\[([^]]*)\]", metadata
    )
    return [
        entry
        for encoded_path in encoded_paths
        for entry in encoded_path.split(":")
        if entry
    ]


def assert_package_local_runtime_search_paths(
    path, metadata, platform_name, local_rpath, dependencies=()
):
    runtime_search_paths = parse_runtime_search_paths(metadata, platform_name)
    if not runtime_search_paths:
        # ELF has no per-dependency location: DT_NEEDED carries bare sonames,
        # and RUNPATH/$ORIGIN is the only thing keeping resolution inside the
        # package.  Its absence stays fatal, with no exception to inspect.
        assert platform_name == "darwin", (
            f"{path} lacks package-local RPATH {local_rpath}:\n{metadata}"
        )
        # Mach-O does: delocate >= 0.13 rewrites every ``@rpath/x`` load
        # command to an explicit ``@loader_path/x`` and then drops the LC_RPATH
        # it no longer needs.  That is stricter than an RPATH, which is a
        # search list the loader could satisfy from elsewhere -- but only if
        # EVERY load command is genuinely loader-relative, so verify that
        # rather than assuming it.  Anything else -- a surviving ``@rpath/``
        # with nothing left to resolve it, an ``@executable_path`` reference,
        # a bare install name, or an absolute build-machine path that happens
        # to exist on this runner and will not on a user's Mac -- fails here.
        assert dependencies, (
            f"{path} has no runtime search path and no dependency list was "
            f"supplied, so nothing was actually verified:\n{metadata}"
        )
        # ``otool -L`` prints the library's own LC_ID_DYLIB first.  That ID is
        # the name consumers link against, not an edge this file resolves, so
        # it needs no search path of its own.
        own_basename = str(path).rsplit("/", 1)[-1]
        system_prefixes = ("/usr/lib/", "/System/Library/")
        nonlocal_dependencies = [
            d for d in dependencies
            if d.rsplit("/", 1)[-1] != own_basename
            and not d.startswith("@loader_path/")
            and not d.startswith(system_prefixes)
        ]
        assert not nonlocal_dependencies, (
            f"{path} has no runtime search path, so every dependency must be "
            f"loader-relative or a macOS system library; these are neither: "
            f"{nonlocal_dependencies}\n{metadata}"
        )
        return
    nonlocal_search_paths = [
        entry for entry in runtime_search_paths
        if entry != local_rpath and not entry.startswith(f"{local_rpath}/")
    ]
    assert not nonlocal_search_paths, (
        f"{path} contains non-package-local runtime search paths "
        f"{nonlocal_search_paths}:\n{metadata}"
    )


def assert_canonical_dependency_graph(
    owner, dependencies, required, canonical_names, platform_name
):
    """Require exactly the canonical SOVERSION name at every D4-stack edge."""
    prefixes = {
        "dftd4": "libdftd4",
        "multicharge": "libmulticharge",
        "mctc": "libmctc-lib",
    }
    stack_entries = {}
    owner_basename = str(owner).rsplit("/", 1)[-1]
    for dependency in dependencies:
        basename = dependency.rsplit("/", 1)[-1]
        logical_name = next(
            (name for name, prefix in prefixes.items()
             if basename.startswith(prefix)),
            None,
        )
        if logical_name is None:
            continue
        canonical = canonical_names[logical_name]
        allowed = {canonical}
        if platform_name == "darwin":
            allowed = {f"@rpath/{canonical}", f"@loader_path/{canonical}"}
        assert dependency in allowed, (
            f"{owner} has noncanonical {logical_name} dependency "
            f"{dependency!r}; expected one of {sorted(allowed)}"
        )
        # ``otool -L`` reports a dylib's LC_ID_DYLIB before its actual load
        # commands.  Validate that install ID above, but do not count it as a
        # dependency edge from the library to itself.
        if platform_name == "darwin" and basename == owner_basename:
            continue
        stack_entries.setdefault(logical_name, []).append(dependency)

    assert set(stack_entries) == set(required), (
        f"{owner} DFT-D4 dependency edges are {stack_entries}; "
        f"expected exactly {sorted(required)}"
    )
    duplicates = {
        name: entries for name, entries in stack_entries.items()
        if len(entries) != 1
    }
    assert not duplicates, f"{owner} has duplicate DFT-D4 edges: {duplicates}"


assert liboqp_path.is_file(), f"missing package-local OpenQP library: {liboqp_path}"
dependency_graph = {
    liboqp_path: {"dftd4", "mctc"},
    d4_paths["dftd4"]: {"multicharge", "mctc"},
    d4_paths["multicharge"]: {"mctc"},
    d4_paths["mctc"]: set(),
}
oqp_deps = ""
dependencies_by_owner = {}
for owner, required_edges in dependency_graph.items():
    metadata = native_metadata(owner, inspect_command)
    if owner == liboqp_path:
        oqp_deps = metadata
    dependencies = parse_dynamic_dependencies(metadata, sys.platform)
    dependencies_by_owner[owner] = dependencies
    assert_canonical_dependency_graph(
        owner, dependencies, required_edges, d4_names, sys.platform
    )

# NLopt was replaced by OpenQP's deterministic simplex-QP solver.  Inspect both
# the native dependency table and defined/undefined dynamic symbols so a stale
# link or an accidentally retained call cannot pass the wheel gate.
assert "nlopt" not in oqp_deps.lower(), (
    f"NLopt dependency leaked into {liboqp_path}:\n{oqp_deps}"
)
if sys.platform == "darwin":
    symbol_commands = [
        ["nm", "-gU", str(liboqp_path)],
        ["nm", "-u", str(liboqp_path)],
    ]
else:
    symbol_commands = [["nm", "-D", str(liboqp_path)]]
symbol_text = "\n".join(
    subprocess.run(command, check=True, capture_output=True, text=True).stdout
    for command in symbol_commands
)
assert not re.search(r"(?i)(?:nlopt|nlo_[a-z0-9_]+)", symbol_text), symbol_text

for path in (liboqp_path, *d4_paths.values()):
    rpath_metadata = native_metadata(path, rpath_command)
    assert_package_local_runtime_search_paths(
        path, rpath_metadata, sys.platform, local_rpath,
        dependencies_by_owner.get(path, ()),
    )

# A repaired wheel must not retain a canonical file while secretly relinking to
# an auditwheel/delocate hash copy in ``*.libs`` or ``.dylibs``.  Inspect the
# complete installed artifact, not just oqp/lib.
artifact_root = Path(oqp_root).parent
stack_prefixes = {
    "dftd4": "libdftd4",
    "multicharge": "libmulticharge",
    "mctc": "libmctc-lib",
}
for logical_name, prefix in stack_prefixes.items():
    candidates = set()
    for candidate in artifact_root.rglob(f"{prefix}*"):
        filename = candidate.name
        is_native = (
            filename.endswith(".dylib") if sys.platform == "darwin"
            else ".so" in filename
        )
        if candidate.is_file() and is_native:
            candidates.add(candidate.absolute())
    expected = {d4_paths[logical_name].absolute()}
    assert candidates == expected, (
        f"installed wheel contains alternate {logical_name} binaries: "
        f"{sorted(map(str, candidates))}; expected only {sorted(map(str, expected))}"
    )


D4_CHILD_MARKER = "OQP_D4_CHILD_RESULT="
D4_CHILD_SCRIPT = r'''
import ctypes
import json
import os
import sys

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
    raise AssertionError(f"DFT-D4 child energy {energy} != {reference}")

if sys.platform == "darwin":
    dyld = ctypes.CDLL(None)
    dyld._dyld_image_count.restype = ctypes.c_uint32
    dyld._dyld_get_image_name.argtypes = [ctypes.c_uint32]
    dyld._dyld_get_image_name.restype = ctypes.c_char_p
    image_paths = []
    for index in range(dyld._dyld_image_count()):
        value = dyld._dyld_get_image_name(index)
        if value:
            image_paths.append(value.decode("utf-8", errors="surrogateescape"))
else:
    image_paths = []
    with open("/proc/self/maps", encoding="utf-8") as maps:
        for line in maps:
            fields = line.rstrip("\n").split(None, 5)
            if len(fields) == 6 and fields[5].startswith("/"):
                image_paths.append(fields[5].removesuffix(" (deleted)"))

prefixes = ("libdftd4", "libmulticharge", "libmctc-lib")
loaded_stack = sorted({
    os.path.realpath(path) for path in image_paths
    if os.path.basename(path).startswith(prefixes)
})
print("OQP_D4_CHILD_RESULT=" + json.dumps({
    "energy": energy,
    "loaded_stack": loaded_stack,
}, sort_keys=True))
'''


def run_d4_child():
    child_env = os.environ.copy()
    external_runtime_path = child_env.pop(
        "OQP_WHEEL_SMOKE_EXTERNAL_RUNTIME_PATH", ""
    )
    for variable in (
        "PYTHONPATH", "OPENQP_ROOT", "LD_LIBRARY_PATH",
        "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH",
    ):
        child_env.pop(variable, None)
    if external_runtime_path:
        runtime_directories = [
            Path(entry) for entry in external_runtime_path.split(os.pathsep)
            if entry
        ]
        assert runtime_directories, "external runtime path is empty"
        for directory in runtime_directories:
            assert directory.is_absolute() and directory.is_dir(), (
                f"invalid external runtime directory: {directory}"
            )
            assert not any(
                candidate.name.startswith(
                    ("libdftd4", "libmulticharge", "libmctc-lib")
                )
                for candidate in directory.iterdir()
            ), f"external runtime directory contains a DFT-D4 library: {directory}"
        child_env["LD_LIBRARY_PATH"] = os.pathsep.join(
            str(directory) for directory in runtime_directories
        )
    return subprocess.run(
        [sys.executable, "-c", D4_CHILD_SCRIPT],
        cwd=tmp,
        env=child_env,
        capture_output=True,
        text=True,
        timeout=120,
    )


def require_successful_d4_child():
    result = run_d4_child()
    assert result.returncode == 0, (
        f"isolated DFT-D4 child failed rc={result.returncode}"
        f"\n--- stdout ---\n{result.stdout[-4000:]}"
        f"\n--- stderr ---\n{result.stderr[-4000:]}"
    )
    payloads = [
        line[len(D4_CHILD_MARKER):]
        for line in result.stdout.splitlines()
        if line.startswith(D4_CHILD_MARKER)
    ]
    assert len(payloads) == 1, (
        f"isolated DFT-D4 child emitted no unique result: {result.stdout[-4000:]}"
    )
    payload = json.loads(payloads[0])
    loaded = {Path(path).resolve() for path in payload["loaded_stack"]}
    expected = {path.resolve() for path in d4_paths.values()}
    assert loaded == expected, (
        "DFT-D4 child loaded noncanonical/alternate libraries: "
        f"{sorted(map(str, loaded))}; expected {sorted(map(str, expected))}"
    )
    return payload


# Prove a fresh interpreter can run D4 and that every canonical file is
# necessary.  Rename one file at a time and always restore it, even when the
# negative child check itself fails.
require_successful_d4_child()
for logical_name, path in d4_paths.items():
    hidden_path = path.with_name(f".{path.name}.wheel-smoke-hidden")
    assert not hidden_path.exists(), f"stale wheel-smoke file: {hidden_path}"
    path.rename(hidden_path)
    try:
        missing_result = run_d4_child()
        assert missing_result.returncode != 0, (
            f"DFT-D4 child still worked with canonical {logical_name} removed; "
            "a hidden alternate copy is satisfying the load\n"
            f"--- stdout ---\n{missing_result.stdout[-4000:]}\n"
            f"--- stderr ---\n{missing_result.stderr[-4000:]}"
        )
    finally:
        hidden_path.rename(path)
require_successful_d4_child()

# 5. The legacy neutral ABI must preserve the validated numerical result.  The
#    high-level charge-aware v2 path must reproduce it for charge=0 and respond
#    to a nonzero molecular charge.
import oqp  # noqa: E402
from oqp.library.single_point import dftd4_native_disp  # noqa: E402

d4_z_values = [8, 1, 1]
d4_xyz_values = [
    0.0, 0.0, 0.2225904,
    0.0, 1.427594, -0.890371,
    0.0, -1.427594, -0.890371,
]
d4_z = oqp.ffi.new("int[]", d4_z_values)
d4_xyz = oqp.ffi.new("double[]", d4_xyz_values)
d4_energy = oqp.ffi.new("double*")
d4_gradient = oqp.ffi.new("double[]", 9)
d4_status = oqp.ffi.new("int*")
oqp.lib.oqp_dftd4_disp(
    3, d4_z, d4_xyz, b"pbe", 3, 1,
    d4_energy, d4_gradient, d4_status,
)
assert d4_status[0] == 0, f"legacy DFT-D4 failed with {d4_status[0]}"
d4_reference_energy = -0.00019542038860006095
d4_reference_gradient = np.array([
    0.0, 1.8882671126329643e-20, -5.9042243157216118e-05,
    0.0, -2.3607262650612779e-05, 2.9521121578608069e-05,
    0.0, 2.3607262650612759e-05, 2.9521121578608056e-05,
])
legacy_gradient = np.frombuffer(
    oqp.ffi.buffer(d4_gradient, 9 * 8), dtype=np.float64
).copy()
assert abs(d4_energy[0] - d4_reference_energy) < 1.0e-12
np.testing.assert_allclose(legacy_gradient, d4_reference_gradient,
                           rtol=0.0, atol=1.0e-12)

xyz_by_atom = np.asarray(d4_xyz_values, dtype=np.float64).reshape(3, 3)
neutral_energy, neutral_gradient = dftd4_native_disp(
    d4_z_values, xyz_by_atom, "pbe", True, total_charge=0.0
)
assert abs(neutral_energy - d4_reference_energy) < 1.0e-12
np.testing.assert_allclose(neutral_gradient.reshape(-1), d4_reference_gradient,
                           rtol=0.0, atol=1.0e-12)
charged_energy, _ = dftd4_native_disp(
    d4_z_values, xyz_by_atom, "pbe", False, total_charge=1.0
)
assert abs(charged_energy - neutral_energy) > 1.0e-12, (
    "DFT-D4 v2 ignored the nonzero molecular charge"
)

# 6. Complete, patched corresponding source is part of the wheel and therefore
#    also survives a Docker pip install when BuildKit cache mounts disappear.
source_root = Path(oqp_root) / "share" / "corresponding-source" / "dftd4-stack"
required_sources = [
    source_root / "README.md",
    source_root / "BUILD-INFO.json",
    source_root / "apply-patch.cmake",
    source_root / "generate-build-info.cmake",
    source_root / "openqp-external-build.cmake",
    source_root / "patches" / "mctc-lib-0.4.2-disable-tests.patch",
    source_root / "patches" / "dftd4-3.7.0-disable-tests-and-mstore.patch",
    source_root / "mctc-lib-0.4.2" / "LICENSE",
    source_root / "multicharge-0.3.0" / "LICENSE",
    source_root / "dftd4-3.7.0" / "COPYING",
    source_root / "dftd4-3.7.0" / "COPYING.LESSER",
]
missing_sources = [str(path) for path in required_sources if not path.is_file()]
assert not missing_sources, f"DFT-D4 corresponding source missing: {missing_sources}"

build_info_path = source_root / "BUILD-INFO.json"
build_info_text = build_info_path.read_text(encoding="utf-8")
build_info = json.loads(build_info_text)
assert build_info["schema"] == "org.open-quantum-platform.dftd4-build-info"
assert build_info["schema_version"] == 1
expected_components = {
    "mctc-lib": {
        "version": "0.4.2",
        "archive_url": (
            "https://github.com/grimme-lab/mctc-lib/archive/refs/tags/"
            "v0.4.2.tar.gz"
        ),
        "sha256": "c7aa45c0a3e6f96e3316e15fc6cdbe48b15234940d3773927a57bb7bfe9744ac",
        "license": "Apache-2.0",
    },
    "multicharge": {
        "version": "0.3.0",
        "archive_url": (
            "https://github.com/grimme-lab/multicharge/archive/refs/tags/"
            "v0.3.0.tar.gz"
        ),
        "sha256": "2fcc1f80871f404f005e9db458ffaec95bb28a19516a0245278cd3175b63a6b2",
        "license": "Apache-2.0",
    },
    "dftd4": {
        "version": "3.7.0",
        "archive_url": (
            "https://github.com/dftd4/dftd4/archive/refs/tags/v3.7.0.tar.gz"
        ),
        "sha256": "f00b244759eff2c4f54b80a40673440ce951b6ddfa5eee1f46124297e056f69c",
        "license": "LGPL-3.0-or-later",
    },
}
components = {
    entry["name"]: {key: value for key, value in entry.items() if key != "name"}
    for entry in build_info["components"]
}
assert components == expected_components, components
assert len(build_info["components"]) == len(expected_components)

build = build_info["build"]
assert build["cmake_version"] and build["generator"]
assert build["system"]["name"] and build["system"]["processor"]
assert build["build_type"] == "Release"
assert build["build_shared_libs"] is True
assert isinstance(build["openmp"], bool)
assert build["blas"]["requested_provider"]
assert build["blas"]["resolved_provider"]
assert build["blas"]["integer_bytes"] == 8
for library_kind in ("resolved_blas_libraries", "resolved_lapack_libraries"):
    library_names = build["blas"][library_kind]
    assert library_names and all(
        "/" not in name and "\\" not in name for name in library_names
    ), (library_kind, library_names)
for compiler in build["compilers"].values():
    assert compiler["id"] and compiler["version"]
    assert compiler["executable"] == Path(compiler["executable"]).name
for flag_name in ("c", "c_release", "fortran", "fortran_release"):
    assert isinstance(build["forwarded_flags"][flag_name], str)

manifest_runtime_names = build_info["canonical_runtime_names"]
assert manifest_runtime_names == {
    "dftd4": d4_names["dftd4"],
    "multicharge": d4_names["multicharge"],
    "mctc-lib": d4_names["mctc"],
}
revision = build_info["openqp"]["source_revision"]
assert revision is None or re.fullmatch(r"[0-9a-fA-F]{40}", revision)
dirty = build_info["openqp"]["source_tree_dirty"]
assert dirty is None or isinstance(dirty, bool)
if revision is not None:
    assert dirty is False

patch_records = {entry["file"]: entry for entry in build_info["patches"]}
expected_patch_names = {
    "mctc-lib-0.4.2-disable-tests.patch": "mctc-lib",
    "dftd4-3.7.0-disable-tests-and-mstore.patch": "dftd4",
}
assert set(patch_records) == set(expected_patch_names)
assert len(build_info["patches"]) == len(expected_patch_names)
for patch_name, component in expected_patch_names.items():
    patch_path = source_root / "patches" / patch_name
    digest = hashlib.sha256(patch_path.read_bytes()).hexdigest()
    assert patch_records[patch_name] == {
        "component": component,
        "file": patch_name,
        "sha256": digest,
        "strip": 1,
    }

for forbidden_path in (
    "/private/tmp/", "/tmp/", "/private/var/folders/", "/root/.cache/",
    "/.cache/openqp/", "/host-cache/", "/Library/Caches/openqp/", "/opt/openqp",
):
    assert forbidden_path not in build_info_text, (
        f"transient build/cache/source path leaked into BUILD-INFO.json: "
        f"{forbidden_path}"
    )

print(
    f"wheel smoke test OK: energy {energy} | basis ABI clean | "
    "replaceable DFT-D4 shared stack and corresponding source verified"
)
