"""Dependency-free checks for the replaceable, charge-aware DFT-D4 seam."""

from __future__ import annotations

import array
import ast
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SINGLE_POINT = ROOT / "pyoqp" / "oqp" / "library" / "single_point.py"
WHEEL_SMOKE = ROOT / ".github" / "scripts" / "wheel_smoke_test.py"
COMPLIANCE = ROOT / "external" / "dftd4-corresponding-source"
BUILD_INFO_GENERATOR = COMPLIANCE / "generate-build-info.cmake"
APPLY_PATCH_HELPER = COMPLIANCE / "apply-patch.cmake"
MCTC_PATCH = COMPLIANCE / "patches" / "mctc-lib-0.4.2-disable-tests.patch"
DFTD4_PATCH = (
    COMPLIANCE / "patches" / "dftd4-3.7.0-disable-tests-and-mstore.patch"
)


class _FakeFFI:
    def new(self, cdecl, initializer=None):
        if cdecl == "int[]":
            return array.array("i", initializer)
        if cdecl == "double[]":
            values = [0.0] * initializer if isinstance(initializer, int) else initializer
            return array.array("d", values)
        if cdecl == "double*":
            return array.array("d", [0.0])
        if cdecl == "int*":
            return array.array("i", [0])
        raise AssertionError(f"unexpected cdecl: {cdecl}")

    @staticmethod
    def buffer(value, size):
        assert size <= value.buffer_info()[1] * value.itemsize
        return memoryview(value)


class _V2Library:
    def __init__(self):
        self.calls = []

    def oqp_dftd4_disp_v2(self, *args):
        self.calls.append(args)
        energy, gradient, status = args[-3:]
        energy[0] = -0.25
        for index in range(len(gradient)):
            gradient[index] = float(index + 1)
        status[0] = 0


class _LegacyLibrary:
    def __init__(self):
        self.calls = []

    def oqp_dftd4_disp(self, *args):
        self.calls.append(args)
        energy, gradient, status = args[-3:]
        energy[0] = -0.125
        status[0] = 0


def _load_dftd4_functions(library):
    """Execute only the two small D4 helpers, avoiding a native oqp import."""
    tree = ast.parse(SINGLE_POINT.read_text(encoding="utf-8"), filename=str(SINGLE_POINT))
    selected = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "_D4_DAMPING_KEYS"
            for target in node.targets
        ):
            selected.append(node)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in {
            "_dftd4_damping_values", "dftd4_native_disp"
        }:
            selected.append(node)
    namespace = {
        "np": np,
        "oqp": SimpleNamespace(ffi=_FakeFFI(), lib=library),
    }
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(SINGLE_POINT), "exec"), namespace)
    return namespace


def _load_wheel_artifact_helpers():
    """Load pure metadata helpers without executing the wheel smoke script."""
    tree = ast.parse(WHEEL_SMOKE.read_text(encoding="utf-8"), filename=str(WHEEL_SMOKE))
    selected = [
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in {
            "parse_dynamic_dependencies", "assert_canonical_dependency_graph",
            "parse_runtime_search_paths",
            "assert_package_local_runtime_search_paths",
        }
    ]
    namespace = {"re": re}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(WHEEL_SMOKE), "exec"), namespace)
    return namespace


def test_v2_receives_charge_and_explicit_damping():
    library = _V2Library()
    functions = _load_dftd4_functions(library)
    energy, gradient = functions["dftd4_native_disp"](
        [8, 1], np.zeros((2, 3)), "pbe", True,
        total_charge=-1.0,
        damping_params={
            "s6": 1.0, "s8": 2.0, "s9": 3.0,
            "a1": 4.0, "a2": 5.0, "alp": 16.0,
        },
    )

    assert energy == -0.25
    np.testing.assert_array_equal(gradient, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    call = library.calls[0]
    assert call[3] == -1.0
    assert call[6] == 1
    assert list(call[7]) == [1.0, 2.0, 3.0, 4.0, 5.0, 16.0]


def test_legacy_fallback_never_discards_charge_or_parameters():
    library = _LegacyLibrary()
    functions = _load_dftd4_functions(library)
    neutral_energy, _ = functions["dftd4_native_disp"](
        [1], np.zeros((1, 3)), "pbe", False
    )
    assert neutral_energy == -0.125
    assert len(library.calls) == 1

    try:
        functions["dftd4_native_disp"](
            [1], np.zeros((1, 3)), "pbe", False, total_charge=1.0
        )
    except RuntimeError as exc:
        assert "charge-aware DFT-D4" in str(exc)
    else:
        raise AssertionError("legacy DFT-D4 silently discarded molecular charge")


def test_shared_stack_packaging_contract_is_declared():
    top_level = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    external = (ROOT / "external" / "CMakeLists.txt").read_text(encoding="utf-8")
    source = (ROOT / "source" / "CMakeLists.txt").read_text(encoding="utf-8")
    header = (ROOT / "include" / "oqp.h").read_text(encoding="utf-8")
    manifest = (
        ROOT / "external" / "dftd4-corresponding-source" / "README.md"
    ).read_text(encoding="utf-8")

    assert '-DBUILD_SHARED_LIBS=ON' in external
    assert 'd4stack-${_OQP_DFTD4_LINKAGE}1' in external
    assert 'd4stack-${_OQP_DFTD4_LINKAGE}1-rpathclean1' in external
    assert '-omp${ENABLE_OPENMP}' in external
    assert '-patch${_OQP_DFTD4_PATCHSET_KEY}' in external
    assert "find_program(OQP_DFTD4_PATCH_EXECUTABLE NAMES patch)" in external
    assert external.count('"-DPATCH_EXECUTABLE=${OQP_DFTD4_PATCH_EXECUTABLE}"') == 2
    assert external.count("apply-patch.cmake") == 2
    assert "DFTD4_DROP_TESTS_PERL" not in external
    assert "DFTD4_DROP_MSTORE_PERL" not in external
    assert "CMAKE_C_FLAGS_RELEASE" in external
    assert "CMAKE_Fortran_FLAGS_RELEASE" in external
    assert external.count("CMAKE_ARGS ${DFTD4_COMMON_ARGS}") == 3
    assert "flags${_oqp_external_flags_hash}" in external
    assert 'libdftd4.a' not in source
    assert 'libmulticharge.a' not in source
    assert 'libmctc-lib.a' not in source
    assert 'INSTALL_RPATH_USE_LINK_PATH FALSE' in source
    assert 'INSTALL_REMOVE_ENVIRONMENT_RPATH TRUE' in source
    assert '-DCMAKE_INSTALL_REMOVE_ENVIRONMENT_RPATH=TRUE' in external
    assert 'set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)' in top_level
    assert 'set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)' not in top_level
    d4_link = re.search(
        r'target_link_libraries\(oqp\s+'
        r'"\$<BUILD_INTERFACE:\$\{DFTD4_DFTD4_LIB\}>"'
        r'.*?\)',
        source,
        re.DOTALL,
    )
    assert d4_link is not None
    assert "DFTD4_MCTC_LIB" in d4_link.group(0)
    assert "DFTD4_MULTICHARGE_LIB" not in d4_link.group(0)
    assert 'oqp_dftd4_disp_v2' in header
    assert 'mctc-lib-0.4.2/' in manifest
    assert 'multicharge-0.3.0/' in manifest
    assert 'dftd4-3.7.0/' in manifest
    assert "BUILD-INFO.json" in source
    assert "DFTD4_MCTC_PATCH" in source
    assert "DFTD4_DFTD4_PATCH" in source


def test_checked_in_build_file_patches_are_stable_and_declared():
    mctc_patch = MCTC_PATCH.read_text(encoding="utf-8")
    dftd4_patch = DFTD4_PATCH.read_text(encoding="utf-8")
    for patch in (mctc_patch, dftd4_patch):
        assert patch.startswith("diff --git a/CMakeLists.txt b/CMakeLists.txt\n")
        assert "--- a/CMakeLists.txt\n+++ b/CMakeLists.txt\n" in patch
        assert 'add_subdirectory("test")' in patch
    assert '-  find_package("mstore")' in dftd4_patch
    assert '+  # find_package("mstore")' in dftd4_patch

    external = (ROOT / "external" / "CMakeLists.txt").read_text(encoding="utf-8")
    for patch_path in (MCTC_PATCH, DFTD4_PATCH):
        digest = hashlib.sha256(patch_path.read_bytes()).hexdigest()
        assert re.fullmatch(r"[0-9a-f]{64}", digest)
        assert patch_path.name in external


def test_patch_helper_is_idempotent_and_rejects_unknown_source_state():
    patch_executable = shutil.which("patch")
    assert patch_executable, "the DFT-D4 build requires the patch executable"
    cmake = "/opt/homebrew/bin/cmake"
    if not Path(cmake).is_file():
        cmake = "cmake"

    with tempfile.TemporaryDirectory(prefix="oqp_d4_patch_helper_") as tmp:
        tmp_path = Path(tmp)
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        target = source_dir / "CMakeLists.txt"
        target.write_text("before\n", encoding="utf-8")
        patch_file = tmp_path / "change.patch"
        patch_file.write_text(
            "diff --git a/CMakeLists.txt b/CMakeLists.txt\n"
            "--- a/CMakeLists.txt\n"
            "+++ b/CMakeLists.txt\n"
            "@@ -1 +1 @@\n"
            "-before\n"
            "+after\n",
            encoding="utf-8",
        )
        command = [
            cmake,
            f"-DPATCH_EXECUTABLE={patch_executable}",
            f"-DSOURCE_DIR={source_dir}",
            f"-DPATCH_FILE={patch_file}",
            "-DPATCH_STRIP=1",
            "-P", str(APPLY_PATCH_HELPER),
        ]
        first = subprocess.run(command, capture_output=True, text=True, timeout=30)
        assert first.returncode == 0, first.stdout + first.stderr
        assert target.read_text(encoding="utf-8") == "after\n"
        second = subprocess.run(command, capture_output=True, text=True, timeout=30)
        assert second.returncode == 0, second.stdout + second.stderr
        assert "already applied" in second.stdout + second.stderr
        assert target.read_text(encoding="utf-8") == "after\n"

        target.write_text("neither\n", encoding="utf-8")
        invalid = subprocess.run(command, capture_output=True, text=True, timeout=30)
        assert invalid.returncode != 0
        assert "matches neither side" in invalid.stdout + invalid.stderr


def test_build_info_generation_is_valid_json_and_scrubs_transient_paths():
    driver_template = r'''
set(CMAKE_SOURCE_DIR "/private/tmp/secret-source")
set(CMAKE_BINARY_DIR "/private/tmp/secret-build")
set(CMAKE_GENERATOR "Ninja")
set(CMAKE_SYSTEM_NAME "Darwin")
set(CMAKE_SYSTEM_PROCESSOR "arm64")
set(CMAKE_C_COMPILER "/toolchain/bin/gcc-15")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "15.2.0")
set(CMAKE_CXX_COMPILER "/toolchain/bin/g++-15")
set(CMAKE_CXX_COMPILER_ID "GNU")
set(CMAKE_CXX_COMPILER_VERSION "15.2.0")
set(CMAKE_Fortran_COMPILER "/toolchain/bin/gfortran-15")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "15.2.0")
set(CMAKE_C_FLAGS [=[-O2 -I/private/tmp/secret-source/include -DNAME="quoted" -DWIN=C:\cache]=])
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_Fortran_FLAGS [=[-O2
-fopenmp]=])
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(PROJECT_VERSION "1.3.0")
set(ENABLE_OPENMP ON)
set(LINALG_LIB "OpenBLAS")
set(_LINALG_LIB_TYPE "other")
set(BLA_SIZEOF_INTEGER 8)
set(BLAS_LIBRARIES "/root/.cache/openqp/lib/libopenblas64.dylib")
set(LAPACK_LIBRARIES "/host-cache/lib/libopenblas64.dylib")
set(OQP_EXTERNALS_ROOT "/root/.cache/openqp/externals")
set(DFTD4_INSTALL "/host-cache/openqp/dftd4")
set(_OQP_MCTC_LIB_VERSION "0.4.2")
set(_OQP_MCTC_LIB_URL "https://example.test/mctc-lib-v0.4.2.tar.gz")
set(_OQP_MCTC_LIB_SHA256 "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
set(_OQP_MULTICHARGE_VERSION "0.3.0")
set(_OQP_MULTICHARGE_URL "https://example.test/multicharge-v0.3.0.tar.gz")
set(_OQP_MULTICHARGE_SHA256 "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb")
set(_OQP_DFTD4_VERSION "3.7.0")
set(_OQP_DFTD4_URL "https://example.test/dftd4-v3.7.0.tar.gz")
set(_OQP_DFTD4_SHA256 "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc")
set(DFTD4_MCTC_PATCH "@MCTC_PATCH@")
set(DFTD4_MCTC_PATCH_SHA256 "@MCTC_PATCH_SHA256@")
set(DFTD4_DFTD4_PATCH "@DFTD4_PATCH@")
set(DFTD4_DFTD4_PATCH_SHA256 "@DFTD4_PATCH_SHA256@")
set(DFTD4_MCTC_RUNTIME_NAME "libmctc-lib.0.dylib")
set(DFTD4_MULTICHARGE_RUNTIME_NAME "libmulticharge.0.dylib")
set(DFTD4_DFTD4_RUNTIME_NAME "libdftd4.3.dylib")
include("@GENERATOR@")
oqp_generate_dftd4_build_info("@OUTPUT@")
'''
    with tempfile.TemporaryDirectory(prefix="oqp_d4_build_info_") as tmp:
        tmp_path = Path(tmp)
        output = tmp_path / "BUILD-INFO.json"
        replacements = {
            "@MCTC_PATCH@": str(MCTC_PATCH),
            "@MCTC_PATCH_SHA256@": hashlib.sha256(MCTC_PATCH.read_bytes()).hexdigest(),
            "@DFTD4_PATCH@": str(DFTD4_PATCH),
            "@DFTD4_PATCH_SHA256@": hashlib.sha256(DFTD4_PATCH.read_bytes()).hexdigest(),
            "@GENERATOR@": str(BUILD_INFO_GENERATOR),
            "@OUTPUT@": str(output),
        }
        driver = driver_template
        for marker, value in replacements.items():
            driver = driver.replace(marker, value)
        driver_path = tmp_path / "generate.cmake"
        driver_path.write_text(driver, encoding="utf-8")
        cmake = "/opt/homebrew/bin/cmake"
        if not Path(cmake).is_file():
            cmake = "cmake"
        result = subprocess.run(
            [cmake, "-P", str(driver_path)],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        raw = output.read_text(encoding="utf-8")
        info = json.loads(raw)

    assert info["schema"] == "org.open-quantum-platform.dftd4-build-info"
    assert info["schema_version"] == 1
    assert info["openqp"] == {
        "version": "1.3.0", "source_revision": None, "source_tree_dirty": None
    }
    assert info["build"]["build_type"] == "Release"
    assert info["build"]["build_shared_libs"] is True
    assert info["build"]["openmp"] is True
    assert info["build"]["blas"]["integer_bytes"] == 8
    assert info["build"]["blas"]["resolved_blas_libraries"] == [
        "libopenblas64.dylib"
    ]
    assert info["build"]["compilers"]["fortran"]["executable"] == "gfortran-15"
    assert '<CMAKE_SOURCE_DIR>/include' in info["build"]["forwarded_flags"]["c"]
    assert '-DNAME="quoted"' in info["build"]["forwarded_flags"]["c"]
    assert "C:\\cache" in info["build"]["forwarded_flags"]["c"]
    assert "\n" in info["build"]["forwarded_flags"]["fortran"]
    for leaked in (
        "/private/tmp/", "/root/.cache/", "/host-cache/",
        str(ROOT), str(output.parent),
    ):
        assert leaked not in raw


def test_high_level_callers_forward_actual_input_charge():
    source = SINGLE_POINT.read_text(encoding="utf-8")
    assert "total_charge = float(mol.config.get('input', {}).get('charge', 0))" in source
    assert "total_charge = float(self.mol.config.get('input', {}).get('charge', 0))" in source
    assert "total_charge=total_charge, damping_params=self.d4_param" in source


def test_wheel_metadata_requires_exact_canonical_soversion_edges():
    helpers = _load_wheel_artifact_helpers()
    parse = helpers["parse_dynamic_dependencies"]
    validate = helpers["assert_canonical_dependency_graph"]

    linux_names = {
        "dftd4": "libdftd4.so.3",
        "multicharge": "libmulticharge.so.0",
        "mctc": "libmctc-lib.so.0",
    }
    linux_metadata = """
 0x0000000000000001 (NEEDED)             Shared library: [libdftd4.so.3]
 0x0000000000000001 (NEEDED)             Shared library: [libmctc-lib.so.0]
 0x0000000000000001 (NEEDED)             Shared library: [libgfortran.so.5]
"""
    linux_dependencies = parse(linux_metadata, "linux")
    assert linux_dependencies == [
        "libdftd4.so.3", "libmctc-lib.so.0", "libgfortran.so.5"
    ]
    validate(
        "liboqp.so", linux_dependencies, {"dftd4", "mctc"},
        linux_names, "linux"
    )
    try:
        validate(
            "liboqp.so",
            ["libdftd4-deadbeef.so.3", "libmctc-lib.so.0"],
            {"dftd4", "mctc"}, linux_names, "linux"
        )
    except AssertionError as exc:
        assert "noncanonical dftd4 dependency" in str(exc)
    else:
        raise AssertionError("hashed auditwheel dependency passed the exact-name gate")

    mac_names = {
        "dftd4": "libdftd4.3.dylib",
        "multicharge": "libmulticharge.0.dylib",
        "mctc": "libmctc-lib.0.dylib",
    }
    mac_metadata = """/wheel/oqp/lib/libdftd4.3.dylib:
\t@rpath/libdftd4.3.dylib (compatibility version 3.0.0, current version 3.7.0)
\t@rpath/libmulticharge.0.dylib (compatibility version 0.0.0, current version 0.3.0)
\t@loader_path/libmctc-lib.0.dylib (compatibility version 0.0.0, current version 0.4.2)
\t/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1345.120.2)
"""
    mac_dependencies = parse(mac_metadata, "darwin")
    assert mac_dependencies[:3] == [
        "@rpath/libdftd4.3.dylib", "@rpath/libmulticharge.0.dylib",
        "@loader_path/libmctc-lib.0.dylib"
    ]
    validate(
        "libdftd4.3.dylib", mac_dependencies, {"multicharge", "mctc"},
        mac_names, "darwin"
    )
    try:
        validate(
            "libdftd4.3.dylib",
            [
                "@loader_path/../.dylibs/libmulticharge-deadbeef.0.dylib",
                "@rpath/libmctc-lib.0.dylib",
            ],
            {"multicharge", "mctc"}, mac_names, "darwin"
        )
    except AssertionError as exc:
        assert "noncanonical multicharge dependency" in str(exc)
    else:
        raise AssertionError("delocate .dylibs dependency passed the exact-name gate")


def test_wheel_metadata_rejects_every_nonlocal_runtime_search_path():
    helpers = _load_wheel_artifact_helpers()
    parse = helpers["parse_runtime_search_paths"]
    validate = helpers["assert_package_local_runtime_search_paths"]

    mac_metadata = """
Load command 20
          cmd LC_RPATH
      cmdsize 32
         path @loader_path (offset 12)
Load command 21
          cmd LC_RPATH
      cmdsize 64
         path @loader_path/vendor (offset 12)
"""
    assert parse(mac_metadata, "darwin") == [
        "@loader_path", "@loader_path/vendor"
    ]
    validate("liboqp.dylib", mac_metadata, "darwin", "@loader_path")

    linux_metadata = """
 0x000000000000001d (RUNPATH)            Library runpath: [$ORIGIN:$ORIGIN/vendor]
"""
    assert parse(linux_metadata, "linux") == ["$ORIGIN", "$ORIGIN/vendor"]
    validate("liboqp.so", linux_metadata, "linux", "$ORIGIN")

    for platform_name, local_rpath, metadata in (
        ("darwin", "@loader_path", mac_metadata + "\n         path /tmp/cache (offset 12)\n"),
        ("linux", "$ORIGIN", linux_metadata.replace("]", ":/tmp/cache]")),
    ):
        try:
            validate("liboqp", metadata, platform_name, local_rpath)
        except AssertionError as exc:
            assert "non-package-local runtime search paths" in str(exc)
        else:
            raise AssertionError("absolute runtime search path passed wheel gate")


def test_wheel_smoke_probes_loaded_paths_and_removal_failures():
    source = WHEEL_SMOKE.read_text(encoding="utf-8")
    assert 'D4_CHILD_MARKER = "OQP_D4_CHILD_RESULT="' in source
    assert 'with open("/proc/self/maps"' in source
    assert "_dyld_image_count" in source
    assert 'artifact_root.rglob(f"{prefix}*")' in source
    assert "path.rename(hidden_path)" in source
    assert "finally:" in source
    assert "hidden_path.rename(path)" in source
    assert "require_successful_d4_child()" in source


def test_distribution_gates_require_build_info_and_exact_patches():
    top_cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    wheel = WHEEL_SMOKE.read_text(encoding="utf-8")
    docker = (ROOT / "Dockerfile").read_text(encoding="utf-8")
    source = (ROOT / "source" / "CMakeLists.txt").read_text(encoding="utf-8")
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    generator = BUILD_INFO_GENERATOR.read_text(encoding="utf-8")
    for gate in (wheel, docker, source):
        assert "BUILD-INFO.json" in gate
        assert "apply-patch.cmake" in gate
        assert "mctc-lib-0.4.2-disable-tests.patch" in gate or (
            "DFTD4_MCTC_PATCH" in gate
        )
        assert "dftd4-3.7.0-disable-tests-and-mstore.patch" in gate or (
            "DFTD4_DFTD4_PATCH" in gate
        )
    assert "json.loads(build_info_text)" in wheel
    assert "hashlib.sha256(patch_path.read_bytes()).hexdigest()" in wheel
    assert "json.loads(build_info_text)" in docker
    assert "hashlib.sha256(patch_path.read_bytes()).hexdigest()" in docker
    assert "recursive-include external/dftd4-corresponding-source *" in manifest
    assert "include(GNUInstallDirs)" in top_cmake
    assert top_cmake.index("include(GNUInstallDirs)") < top_cmake.index(
        "add_subdirectory(source)"
    )
    assert "string(JSON _schema" in generator
    assert "source_revision" in generator
    assert "--untracked-files=normal" in generator
    helper = APPLY_PATCH_HELPER.read_text(encoding="utf-8")
    assert "--dry-run -N" in helper
    assert "--dry-run -R" in helper


if __name__ == "__main__":
    test_v2_receives_charge_and_explicit_damping()
    test_legacy_fallback_never_discards_charge_or_parameters()
    test_shared_stack_packaging_contract_is_declared()
    test_checked_in_build_file_patches_are_stable_and_declared()
    test_patch_helper_is_idempotent_and_rejects_unknown_source_state()
    test_build_info_generation_is_valid_json_and_scrubs_transient_paths()
    test_high_level_callers_forward_actual_input_charge()
    test_wheel_metadata_requires_exact_canonical_soversion_edges()
    test_wheel_metadata_rejects_every_nonlocal_runtime_search_path()
    test_wheel_smoke_probes_loaded_paths_and_removal_failures()
    test_distribution_gates_require_build_info_and_exact_patches()
    print("DFT-D4 shared-interface unit tests passed")
