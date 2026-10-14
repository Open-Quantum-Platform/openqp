"""The open-shell XC response must be split and reduced over the MPI ranks.

utddft_fxc (the UHF/ROHF XC response kernel) left its consumer's parallel
environment at the default -- one rank, no MPI -- so under mpiexec every rank
integrated the whole grid on its own and nothing reduced the results.  With
OpenMP > 1 the ranks' responses then differed in the last bits.  Native TRAH
runs its Davidson on every rank, so near convergence the ranks took different
Newton steps and the reduced Fock of their inconsistent orbitals was garbage:
ROKS BLYP/6-31G* CH2O, mpiexec -n 2, OMP_NUM_THREADS=2 went -114.34 -> -153.03
Eh and hung, while the serial run converged to -114.3412610598 Eh.

The test runs exactly that case on one rank and then three times on two
ranks, and requires the converged one-rank energy every time.  Whether the
ranks' last bits split a run is timing-dependent; on chc r630 (OpenMPI 4.1.5,
OMP 2) the unfixed build failed 9 of 9 counted two-rank runs (14 of 15 in
all) and the fixed build passed 9 of 9, so three runs make a silent pass of
the defect unlikely.  Both runs go through mpiexec: an MPI build cannot
initialise MPI without a launcher, so the test does not import oqp in the
pytest process either.  Without mpiexec, or with a build whose log does not
report two ranks (built without MPI), it is skipped; CI runs pytest only on
its non-MPI legs, so this test is executed by hand on an MPI build.
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
type=rohf
multiplicity=3
maxit=100
conv=1e-11
converger_type=trah
trh_impl=native
escalation=trah
"""

FINAL = re.compile(r"Final ROHF energy is\s+(-?\d+\.\d+)")
RANKS = re.compile(r"MPI Processors:\s+(\d+)")


pytestmark = pytest.mark.skipif(shutil.which("mpiexec") is None,
                                reason="mpiexec not available")


def _launcher(nranks):
    version = subprocess.run(["mpiexec", "--version"], capture_output=True, text=True)
    cmd = ["mpiexec", "-n", str(nranks)]
    text = version.stdout + version.stderr
    if "Open MPI" in text or "OpenRTE" in text:
        cmd.insert(1, "--oversubscribe")
    return cmd


def _run(tmp_path, name, launcher):
    work = tmp_path / name
    work.mkdir()
    (work / "ch2o.inp").write_text(DECK)
    env = dict(os.environ, OMP_NUM_THREADS="2", MKL_NUM_THREADS="1")
    proc = subprocess.run(launcher + [sys.executable, "-m", "oqp.pyoqp", "ch2o.inp"],
                          cwd=work, env=env, capture_output=True, text=True,
                          timeout=300)
    log_path = work / "ch2o.log"
    log = log_path.read_text() if log_path.exists() else ""
    return proc, log


def test_rohf_trah_dft_converges_identically_on_two_ranks(tmp_path):
    single, slog = _run(tmp_path, "one_rank", _launcher(1))
    assert single.returncode == 0, slog[-3000:] + single.stdout[-2000:] + single.stderr[-2000:]
    e_serial = float(FINAL.findall(slog)[-1])
    for attempt in range(3):
        parallel, plog = _run(tmp_path, f"two_ranks_{attempt}", _launcher(2))
        # A launcher or run failure (divergence, timeout) is a failure, not a skip.
        assert parallel.returncode == 0, plog[-3000:] + parallel.stdout[-2000:] + parallel.stderr[-2000:]
        ranks = RANKS.findall(plog)
        if not ranks or int(ranks[-1]) != 2:
            pytest.skip("OpenQP was built without MPI (the run reports no 2-rank job)")
        assert "CONVERGED" in plog, plog[-3000:]
        assert float(FINAL.findall(plog)[-1]) == pytest.approx(e_serial, abs=1e-9)
