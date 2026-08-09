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


oqp_deps = native_metadata(liboqp_path, inspect_command)
assert "libdftd4" in oqp_deps and "libmctc-lib" in oqp_deps, oqp_deps
d4_deps = native_metadata(d4_paths["dftd4"], inspect_command)
assert "libmulticharge" in d4_deps and "libmctc-lib" in d4_deps, d4_deps
multicharge_deps = native_metadata(d4_paths["multicharge"], inspect_command)
assert "libmctc-lib" in multicharge_deps, multicharge_deps

for path in (liboqp_path, *d4_paths.values()):
    rpath_metadata = native_metadata(path, rpath_command)
    assert local_rpath in rpath_metadata, (
        f"{path} lacks package-local RPATH {local_rpath}:\n{rpath_metadata}"
    )

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
    source_root / "openqp-external-build.cmake",
    source_root / "mctc-lib-0.4.2" / "LICENSE",
    source_root / "multicharge-0.3.0" / "LICENSE",
    source_root / "dftd4-3.7.0" / "COPYING",
    source_root / "dftd4-3.7.0" / "COPYING.LESSER",
]
missing_sources = [str(path) for path in required_sources if not path.is_file()]
assert not missing_sources, f"DFT-D4 corresponding source missing: {missing_sources}"

print(
    f"wheel smoke test OK: energy {energy} | basis ABI clean | "
    "replaceable DFT-D4 shared stack and corresponding source verified"
)
