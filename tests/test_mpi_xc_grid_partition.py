"""Ranks with different OpenMP thread counts must integrate the XC grid once.

run_xc deals the grid slices to the MPI ranks in chunks, by chunk index, so
every rank must cut the grid at the same places.  The chunk size used to be
computed inside the OpenMP region from the rank's own thread count; a rank
with 1 thread and a rank with 3 threads then cut differently, and some slices
were integrated twice and others not at all.  CH2O BLYP/6-31G* ROKS with one
rank at OMP_NUM_THREADS=1 and one at 3 integrated 16.76 instead of 16
electrons and converged to -113.63 instead of -114.34 Eh.

The test runs RKS CH2O on one rank and then as an MPMD job whose two ranks
have 1 and 3 OpenMP threads, and requires the XC grid to hold exactly the 16
electrons and the one-rank SCF energy.  It needs mpiexec and an MPI build
(skipped otherwise; CI runs pytest only without MPI).
"""

import os
import re
import shutil
import subprocess
import sys

import pytest

DECK = """[input]
system=
   C 0.000000 0.000000 -0.528500
   O 0.000000 0.000000  0.677500
   H 0.000000 0.935100 -1.112000
   H 0.000000 -0.935100 -1.112000
charge=0
runtype=energy
basis=6-31g*
functional=blyp
method=hf

[guess]
type=huckel

[scf]
type=rhf
multiplicity=1
conv=1e-8
"""

FINAL = re.compile(r"Final RHF energy is\s+(-?\d+\.\d+)")
NELEC = re.compile(r"DFT: total electron density =\s+(-?\d+\.\d+)")
RANKS = re.compile(r"MPI Processors:\s+(\d+)")

pytestmark = pytest.mark.skipif(shutil.which("mpiexec") is None,
                                reason="mpiexec not available")


def _oversubscribe():
    version = subprocess.run(["mpiexec", "--version"], capture_output=True, text=True)
    text = version.stdout + version.stderr
    return ["--oversubscribe"] if ("Open MPI" in text or "OpenRTE" in text) else []


def _run(tmp_path, name, threads):
    work = tmp_path / name
    work.mkdir()
    (work / "ch2o.inp").write_text(DECK)
    cmd = ["mpiexec"] + _oversubscribe()
    for i, nthreads in enumerate(threads):
        if i:
            cmd.append(":")
        cmd += ["-n", "1", "env", f"OMP_NUM_THREADS={nthreads}", "MKL_NUM_THREADS=1",
                sys.executable, "-m", "oqp.pyoqp", "ch2o.inp"]
    proc = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=300)
    log_path = work / "ch2o.log"
    log = log_path.read_text() if log_path.exists() else ""
    assert proc.returncode == 0, log[-3000:] + proc.stdout[-2000:] + proc.stderr[-2000:]
    return log


def test_xc_grid_is_integrated_once_with_unequal_thread_counts(tmp_path):
    mixed = _run(tmp_path, "threads_1_3", [1, 3])
    ranks = RANKS.findall(mixed)
    if not ranks or int(ranks[-1]) != 2:
        pytest.skip("OpenQP was built without MPI (the run reports no 2-rank job)")
    single = _run(tmp_path, "one_rank", [2])
    assert float(NELEC.findall(single)[-1]) == pytest.approx(16.0, abs=1e-3)
    assert float(NELEC.findall(mixed)[-1]) == pytest.approx(16.0, abs=1e-3)
    assert float(FINAL.findall(mixed)[-1]) == pytest.approx(
        float(FINAL.findall(single)[-1]), abs=1e-8)
