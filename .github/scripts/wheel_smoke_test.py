#!/usr/bin/env python3
"""Numerical smoke test run against every built wheel (CIBW_TEST_COMMAND).

Goes beyond the old import/file-exists check: runs a real RHF/STO-3G
calculation and asserts the energy AND the raw basis-metadata integer ABI.
These are exactly the two failure classes that previously shipped in wheels
while ``import oqp`` still succeeded:
  * a wrong/slow BLAS selection (silent NetLib fallback) -- caught here by
    requiring the calculation to complete with the reference energy;
  * the LP64 integer-ABI misread (oqp_get_basis exporting build-width
    integers read as int64 by Python), which corrupted basis metadata
    (centers [0,1] -> [8589934592, -1]) and crashed the molden writer.

Usage: wheel_smoke_test.py <project-dir>
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

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

print(f"wheel smoke test OK: energy {energy} | basis metadata ABI clean")
