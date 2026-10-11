"""Finite-difference displacement workers must not re-import oqp.

The numerical Hessian and numerical NAC drivers fan their displaced-geometry
jobs out over a pool.  Without MPI each job is an external ``python -m
oqp.pyoqp`` subprocess, so the pool workers only write an input file and wait.
They used to be ``multiprocessing.Pool`` processes.  Under the spawn start
method (macOS default; Linux from Python 3.14) every worker is a fresh
interpreter that re-imports ``oqp`` with the parent's ``sys.path`` just to
unpickle ``grad_wrapper``.  When that import failed the worker died inside the
pool's task ``get()``, the pool replaced it without resubmitting the task, and
``imap_unordered`` waited forever.  Observed 2026-10-11 on macOS: the full
pytest suite hung in ``test_cam_hessian`` because earlier-collected test
modules had put the unbuilt source tree ``pyoqp/`` first on ``sys.path``.

The end-to-end test below recreates that condition deliberately -- spawn start
method, and an ``oqp`` package that raises on import placed first on
``sys.path`` after the real one is loaded -- so it fails (times out) on any
platform if the workers ever again depend on re-importing ``oqp``.
"""
import os
import subprocess
import sys
import textwrap
import unittest


def _runtime_available():
    # Probe in a clean interpreter: this process's sys.path may already carry
    # the source tree ahead of the installed package.
    try:
        probe = subprocess.run([sys.executable, "-c", "import oqp.pyoqp"],
                               capture_output=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return probe.returncode == 0


RUNTIME = _runtime_available()

DRIVER = textwrap.dedent('''
    import multiprocessing
    import pathlib
    import sys

    import numpy as np

    multiprocessing.set_start_method("spawn", force=True)
    from oqp.pyoqp import Runner  # the real package, imported first

    work = pathlib.Path(sys.argv[1])
    poison = work / "poison" / "oqp"
    poison.mkdir(parents=True)
    (poison / "__init__.py").write_text(
        "raise ImportError('a displacement worker re-imported oqp')\\n")
    sys.path.insert(0, str(poison.parent))

    inp = work / "h2.inp"
    inp.write_text("""[input]
    system=
       1   0.000000000   0.000000000   0.000000000
       1   0.000000000   0.000000000   0.740000000
    charge=0
    runtype=hess
    basis=sto-3g
    method=hf
    [guess]
    type=huckel
    [scf]
    multiplicity=1
    type=rhf
    [hess]
    type=numerical
    state=0
    dx=0.005
    nproc=2
    """)
    runner = Runner(project="h2", input_file=str(inp), log=str(work / "h2.log"))
    runner.run()
    hessian = np.asarray(runner.mol.hessian, dtype=float)
    assert hessian.shape == (6, 6), hessian.shape
    assert np.all(np.isfinite(hessian))
    assert np.max(np.abs(hessian)) > 1.0e-3
    print("DISPLACEMENT_POOL_OK")
''')


@unittest.skipUnless(RUNTIME, "compiled OpenQP runtime not available")
class DisplacementPoolTests(unittest.TestCase):
    def test_numerical_hessian_completes_when_workers_cannot_import_oqp(self):
        import tempfile
        with tempfile.TemporaryDirectory(prefix="oqp_disp_pool_") as work:
            env = dict(os.environ)
            env.pop("OPENQP_ROOT", None)
            env.setdefault("OMP_NUM_THREADS", "1")
            try:
                proc = subprocess.run(
                    [sys.executable, "-c", DRIVER, work], cwd=work, env=env,
                    capture_output=True, text=True, timeout=600)
            except subprocess.TimeoutExpired as exc:
                self.fail("numerical Hessian hung: its displacement workers "
                          "depend on re-importing oqp in a fresh interpreter\n"
                          f"{(exc.stderr or b'')[-4000:]!r}")
        self.assertEqual(proc.returncode, 0, proc.stderr[-4000:])
        self.assertIn("DISPLACEMENT_POOL_OK", proc.stdout)

    def test_pool_runs_unpicklable_work_and_reports_worker_errors(self):
        from oqp.library.single_point import _displacement_pool

        class Serial:
            use_mpi = 0

        pool, ncpu = _displacement_pool(Serial(), 12, 4)
        self.assertEqual(ncpu, 4)
        offset = 7  # a closure: cannot be pickled to a process worker
        try:
            self.assertEqual(sorted(pool.imap_unordered(lambda x: x + offset, range(12))),
                             [x + offset for x in range(12)])

            def boom(x):
                if x == 5:
                    raise RuntimeError("worker failure")
                return x

            with self.assertRaisesRegex(RuntimeError, "worker failure"):
                list(pool.imap_unordered(boom, range(12)))
        finally:
            pool.close()
            pool.join()

        pool, ncpu = _displacement_pool(Serial(), 3, 8)
        pool.close()
        pool.join()
        self.assertEqual(ncpu, 3)

    def test_mpi_ranks_remain_the_workers(self):
        from oqp.library.single_point import _displacement_pool
        from oqp.utils.mpi_utils import MPIPool

        class Comm:
            size = 4

        class Parallel:
            use_mpi = 1
            comm = Comm()

        pool, ncpu = _displacement_pool(Parallel(), 12, 1)
        self.assertIsInstance(pool, MPIPool)
        self.assertEqual(ncpu, 4)


if __name__ == "__main__":
    unittest.main()
