"""A level-shifted DIIS SCF may only report convergence on a converged iteration.

With a level shift the DIIS loop first converges with the shift, prints
"Performing a last SCF with zero VSHIFT." and builds one more Fock without it.
The loop used to leave on that unshifted build whatever its DIIS error was
(``elseif (vshift_last_iter) ... exit`` in scf_driver): an error above conv
skipped the convergence test, yet ``iter <= maxit`` after the loop printed
"SCF convergence achieved", set SCF_converged and passed the orbitals on as a
supplied guess.  The density was not rebuilt from the last diagonalisation
either.

The error after the shift is removed is the error of one more DIIS step, so it
exceeds conv whenever DIIS is oscillating.  NO (doublet) ROHF/PBE/6-31G with a
large shift oscillates by one to two decades near convergence.  On GitHub main
5cdbddc68 (chc r630, GCC 12.3, MKL ILP64) the zero-shift iteration had, e.g.

    vshift=2.8 conv=1e-6  1 thread   error 5.83e-06  "SCF convergence achieved"
    vshift=2.8 conv=1e-6  4 threads  error 1.23e-05  "SCF convergence achieved"
    vshift=2.6 conv=1e-6  4 threads  error 1.49e-04  "SCF convergence achieved"

Which shift lands on an oscillation peak depends on BLAS and thread count, so
several shifts are run.  The asserted property does not depend on that: every
"SCF convergence achieved" must close a DIIS table whose last iteration is
below conv, and the final energy must be the converged ROHF energy.
"""

import os
import re
import subprocess
import sys

import pytest

NO_ROHF_PBE = """\
[input]
system=
   7   0.000000000   0.000000000   0.000000000
   8   0.000000000   0.000000000   1.150000000
charge=0
runtype=energy
basis=6-31g
method=hf
functional=pbe

[guess]
type=huckel

[scf]
type=rohf
multiplicity=2
maxit=200
conv=1.0e-6
vshift={vshift}
"""

# Converged ROHF/PBE/6-31G energy of this geometry (unshifted runs, chc r630).
E_REF = -129.7042413182

LAST_SCF = "Performing a last SCF with zero VSHIFT."
CONV = 1.0e-6
ITER = re.compile(r"^\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+\.\d+)\s+(\S+)\s*$")
# SOSCF/BFGS rows (rms gradient and rms step instead of the DIIS error).
ITER_OTHER = re.compile(r"^\s+\d+\s+-?\d\.\d+E[+-]\d+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\.\d+\s+\S+\s*$")
TABLE_BREAKS = ("iterations begin", "Macro", "stalled", "did not converge")
FINAL = re.compile(r"Final \S+ energy is\s+(-?\d+\.\d+) after")


def _have_native_oqp():
    try:
        import oqp  # noqa: F401
    except Exception:
        return False
    return True


pytestmark = pytest.mark.skipif(
    not _have_native_oqp(),
    reason="native OpenQP library not importable",
)


def _run(tmp_path, name, deck, threads):
    inp = tmp_path / f"{name}.inp"
    inp.write_text(deck)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    proc = subprocess.run(
        [sys.executable, "-m", "oqp.pyoqp", str(inp)],
        capture_output=True, text=True, cwd=str(tmp_path), env=env, timeout=900,
    )
    log_path = inp.with_suffix(".log")
    log = log_path.read_text() if log_path.exists() else ""
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    return log


def _convergence_claims(log):
    """DIIS error of the last iteration before each 'SCF convergence achieved'.

    None when the claim does not close a DIIS table (SOSCF prints gradient
    columns instead of the DIIS error, TRAH prints its own macro table).
    """
    claims = []
    last_err = None
    for line in log.splitlines():
        m = ITER.match(line)
        if m:
            last_err = float(m.group(5))
        elif "SCF convergence achieved" in line:
            claims.append(last_err)
            last_err = None
        elif any(k in line for k in TABLE_BREAKS) or ITER_OTHER.match(line):
            last_err = None
    return claims


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("vshift", [2.4, 2.6, 2.8, 3.0])
def test_no_rohf_large_vshift_reports_convergence_only_below_conv(
        tmp_path, vshift, threads):
    log = _run(tmp_path, f"no_rohf_vs{vshift}_t{threads}",
               NO_ROHF_PBE.format(vshift=vshift), threads)
    # The zero-shift step must really have been taken.
    assert LAST_SCF in log
    claims = _convergence_claims(log)
    assert claims, log[-3000:]
    for err in claims:
        # A claim that does not follow an iteration line (TRAH prints none) is
        # judged by TRAH's own criterion; only DIIS tables are checked here.
        if err is not None:
            assert err < CONV, (
                f"'SCF convergence achieved' after an iteration with DIIS "
                f"error {err:.3e} >= conv {CONV:.0e}\n" + log[-4000:])
    energies = [float(e) for e in FINAL.findall(log)]
    assert energies, log[-3000:]
    # Sanity only (DIIS conv=1e-6 leaves ~1e-6 Eh); the error check above is the gate.
    assert energies[-1] == pytest.approx(E_REF, abs=1.0e-5)
