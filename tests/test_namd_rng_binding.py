"""ABI and frozen-vector guards for the resident Fortran NAMD counter RNG."""

import json
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_namd_rng_is_fortran_resident_and_c_interoperable():
    source = (ROOT / "source" / "modules" / "namd.F90").read_text()
    header = (ROOT / "include" / "oqp.h").read_text()
    wrapper = (ROOT / "pyoqp" / "oqp" / "__init__.py").read_text()

    assert "pure function namd_counter_random(seed, stream, step)" in source
    assert 'bind(C, name="oqp_namd_counter_random")' in source
    assert 'bind(C, name="oqp_namd_counter_normal_fill")' in source
    assert "double oqp_namd_counter_random(int64_t seed, int64_t stream, int64_t step);" in header
    assert "void oqp_namd_counter_normal_fill(int64_t seed, int64_t stream, int64_t count," in header
    assert "'oqp_namd_counter_random'" in wrapper
    assert "'oqp_namd_counter_normal_fill'" in wrapper
    assert "pure function namd_add64_mod(a, b)" in source
    assert "pure function namd_mul64_mod(a, b)" in source


def test_built_namd_rng_matches_frozen_splitmix64_vectors():
    """Exercise the actual C ABI when a local OpenQP build is available."""
    script = """
import json
import numpy as np
import oqp
cases = [(0, 0, 0), (1, 0, 1), (20260803, 0, 1),
         (20260803, 0, 2), (20260803, 1, 2), (-1, 0, 2)]
normal0 = np.empty(6, dtype=np.float64)
normal1 = np.empty(6, dtype=np.float64)
normal0_repeat = np.empty(6, dtype=np.float64)
for stream, values in ((0, normal0), (1, normal1), (0, normal0_repeat)):
    oqp.oqp_namd_counter_normal_fill(
        20260803, stream, values.size,
        oqp.ffi.cast('double *', values.ctypes.data))
print('NAMD_RNG=' + json.dumps({
    'uniform': [oqp.oqp_namd_counter_random(*case) for case in cases],
    'normal0': normal0.tolist(),
    'normal1': normal1.tolist(),
    'normal0_repeat': normal0_repeat.tolist(),
}))
"""
    env = os.environ.copy()
    pythonpath = str(ROOT / "pyoqp")
    if env.get("PYTHONPATH"):
        pythonpath += os.pathsep + env["PYTHONPATH"]
    env["PYTHONPATH"] = pythonpath
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        pytest.skip("compiled OpenQP runtime is not available")

    marker = next(
        (line for line in result.stdout.splitlines() if line.startswith("NAMD_RNG=")),
        None,
    )
    assert marker is not None, result.stdout + result.stderr
    values = json.loads(marker.removeprefix("NAMD_RNG="))
    assert values["uniform"] == [
        0.42377815504956573,
        0.712959830204546,
        0.27463452015027734,
        0.8410083610050165,
        0.6348865656316214,
        0.9164863378832538,
    ]
    expected_normal0 = [
        -0.05626842588421754,
        -1.3031353299334301,
        0.33707926150078443,
        0.3411069899205209,
        -0.9355022113288782,
        0.25568063822139614,
    ]
    assert values["normal0"] == pytest.approx(expected_normal0, abs=2e-15)
    assert values["normal0_repeat"] == values["normal0"]
    assert values["normal1"] != values["normal0"]
    assert all(math.isfinite(value) for value in values["normal0"])
