"""TDDFT/TDHF RPA and TDA excitation energies under MPI equal the serial ones.

The XC kernel contraction distributes grid slices over the MPI ranks
(dft_gridint: ``mod(iChunk/chunkSize, pe%size) == pe%rank``) and reduces the
response Fock matrices afterwards.  A buffer missing from that reduction, or a
rank-local quantity used before it, changes the excitation energies with the
rank count.  Covered: the closed-shell singlet and triplet branches, the
CONTIGUOUS MO-pointer consumer (GGA and the meta-GGA ``moG1`` branch), the
spin-antisymmetric triplet contraction and the CAM two-pass integrals.

Runs only with mpi4py and an MPI launcher.  ``OQP_REQUIRE_MPI_TESTS=1`` (set
in the CI MPI leg) turns the skip into a failure, so a broken MPI environment
cannot pass silently.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

H2O_C1 = ("   O 0.010000 -0.020000 0.117300\n   H 0.080000 0.770200 -0.449200\n"
          "   H -0.050000 -0.731200 -0.509200\n")

CASES = [  # functional, [tdhf] type, multiplicity
    ("b3lypv5", "rpa", 1),
    ("b3lypv5", "rpa", 3),
    ("b3lypv5", "tda", 3),
    ("m06-2x", "rpa", 3),
    ("cam-b3lyp", "rpa", 3),
]
RANKS = (2, 4)
REQUIRED = os.environ.get("OQP_REQUIRE_MPI_TESTS", "") == "1"


def _launcher():
    try:
        import mpi4py  # noqa: F401
        import oqp  # noqa: F401
    except (ImportError, OSError):  # not installed, or liboqp/libmpi not loadable
        return None
    return shutil.which("mpiexec") or shutil.which("mpirun")


def _deck(functional, td_type, mult):
    return (f"[input]\nsystem=\n{H2O_C1}charge=0\nruntype=energy\nbasis=6-31g*\n"
            f"method=tdhf\nfunctional={functional}\n\n[guess]\ntype=huckel\n\n"
            "[scf]\ntype=rhf\nmultiplicity=1\nmaxit=100\nconv=1e-10\n\n"
            f"[tdhf]\ntype={td_type}\nnstate=4\nmaxit=100\nconv=1e-10\nmultiplicity={mult}\n\n"
            "[dftgrid]\nrad_npts=96\nang_npts=302\n")


def _run(text, nrank, launcher):
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "job.inp").write_text(text)
        env = dict(os.environ, OMP_NUM_THREADS="1")
        cmd = [sys.executable, "-m", "oqp.pyoqp", "job.inp"]
        if nrank > 1:
            cmd = [launcher, "-n", str(nrank), "--oversubscribe", "--bind-to", "none"] + cmd
        proc = subprocess.run(cmd, cwd=tmp, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=900)
        log_path = Path(tmp) / "job.log"
        log = log_path.read_text(errors="ignore") if log_path.exists() else ""
        summary = log.split("Summary of the TD-DFT calculation")[-1]
        ev = [float(x) for x in re.findall(r"^\s*0 ->\s*\d+\s+(-?[\d.]+)", summary, re.M)]
        return proc.returncode, ev, proc.stdout.decode(errors="ignore")[-1500:] + log[-2000:]


class TestTddftUnderMpi(unittest.TestCase):
    def setUp(self):
        self.launcher = _launcher()
        if self.launcher is None:
            if REQUIRED:
                self.fail("OQP_REQUIRE_MPI_TESTS=1 but mpi4py, oqp or an MPI launcher is missing")
            self.skipTest("mpi4py, the compiled OpenQP runtime and an MPI launcher are needed")

    def test_excitation_energies_independent_of_rank_count(self):
        for functional, td_type, mult in CASES:
            text = _deck(functional, td_type, mult)
            rc0, e0, d0 = _run(text, 1, self.launcher)
            self.assertEqual(rc0, 0, d0)
            self.assertEqual(len(e0), 4, d0)
            for nrank in RANKS:
                with self.subTest(functional=functional, td=td_type, mult=mult, ranks=nrank):
                    rc, e, d = _run(text, nrank, self.launcher)
                    self.assertEqual(rc, 0, d)
                    self.assertEqual(len(e), 4, d)
                    for a, b in zip(e0, e):
                        self.assertAlmostEqual(a, b, delta=1e-6)


if __name__ == "__main__":
    unittest.main()
