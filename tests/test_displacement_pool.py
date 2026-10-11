"""Finite-difference displacement workers must not re-import oqp.

The numerical Hessian and numerical NAC drivers fan their displaced-geometry
jobs out over a pool.  Without MPI each job is an external ``python -m
oqp.pyoqp`` subprocess, so the pool workers only write an input file and wait.
They used to be ``multiprocessing.Pool`` processes.  Under the spawn start
method (macOS default) or forkserver (Linux default from Python 3.14) every
worker unpickles its task in a process that has not imported ``oqp``, so it
re-imports ``oqp`` with the parent's ``sys.path``.  When that import failed the
worker died inside the pool's task ``get()``, the pool replaced it without
resubmitting the task, and ``imap_unordered`` waited forever.  Observed
2026-10-11 on macOS: the full pytest suite hung in ``test_cam_hessian`` because
earlier-collected test modules had put the unbuilt source tree ``pyoqp/`` first
on ``sys.path``.

Every check runs in a fresh interpreter, so it sees the installed package and
is independent of what this pytest process has done to ``sys.path`` or
``sys.modules``.  The end-to-end checks recreate the failure condition
deliberately -- spawn start method, and an ``oqp`` package that raises on
import placed first on ``sys.path`` after the real one is loaded -- so they
fail (time out) on any platform if either driver's workers ever again depend
on re-importing ``oqp``.
"""
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest


def _runtime_available():
    try:
        probe = subprocess.run([sys.executable, "-c", "import oqp.pyoqp"],
                               capture_output=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return probe.returncode == 0


RUNTIME = _runtime_available()

POISON = textwrap.dedent('''
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
''')

HESSIAN_DRIVER = POISON + textwrap.dedent('''
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

NAC_DRIVER = POISON + textwrap.dedent('''
    from oqp.library.single_point import NAC

    inp = work / "h2o_nac.inp"
    inp.write_text("""[input]
    system=
       8   0.000000000   0.000000000  -0.041061554
       1  -0.533194329   0.533194329  -0.614469223
       1   0.533194329  -0.533194329  -0.614469223
    charge=0
    runtype=energy
    basis=6-31g*
    functional=bhhlyp
    method=tdhf
    [guess]
    type=huckel
    [scf]
    multiplicity=3
    type=rohf
    [tdhf]
    type=mrsf
    nstate=3
    [nac]
    type=numerical
    states=1 2
    dx=0.001
    nproc=2
    """)
    runner = Runner(input_file=str(inp), log=str(work / "h2o_nac.log"))
    runner.run()
    nacv, dcv, flags = NAC(runner.mol).numerical_nac()
    assert set(flags) == {"computed"}, flags
    assert np.all(np.isfinite(dcv))
    print("DISPLACEMENT_POOL_OK")
''')

UNIT_DRIVER = textwrap.dedent('''
    from oqp.library.single_point import _displacement_pool
    from oqp.utils.mpi_utils import MPIPool

    class Serial:
        use_mpi = 0

    pool, ncpu = _displacement_pool(Serial(), 12, 4)
    assert ncpu == 4, ncpu
    offset = 7  # a closure: cannot be pickled to a process worker
    try:
        got = sorted(pool.imap_unordered(lambda x: x + offset, range(12)))
        assert got == [x + offset for x in range(12)], got

        def boom(x):
            if x == 5:
                raise RuntimeError("worker failure")
            return x

        try:
            list(pool.imap_unordered(boom, range(12)))
        except RuntimeError as exc:
            assert "worker failure" in str(exc)
        else:
            raise AssertionError("worker exception did not reach the caller")
    finally:
        pool.close()
        pool.join()

    pool, ncpu = _displacement_pool(Serial(), 3, 8)
    pool.close()
    pool.join()
    assert ncpu == 3, ncpu

    class Comm:
        size = 4

    class Parallel:
        use_mpi = 1
        comm = Comm()

    pool, ncpu = _displacement_pool(Parallel(), 12, 1)
    assert isinstance(pool, MPIPool), type(pool)
    assert ncpu == 4, ncpu
    print("DISPLACEMENT_POOL_OK")
''')


TERMINATE_DRIVER = textwrap.dedent('''
    import pathlib
    import sys
    import time

    from oqp.pyoqp import Runner
    from oqp.library import single_point

    work = pathlib.Path(sys.argv[1])
    inp = work / "h2.inp"
    inp.write_text("""[input]
    system=
       1   0.000000000   0.000000000   0.000000000
       1   0.000000000   0.000000000   0.740000000
    charge=0
    runtype=energy
    basis=sto-3g
    method=hf
    [guess]
    type=huckel
    [scf]
    multiplicity=1
    type=rhf
    [hess]
    type=numerical
    nproc=1
    """)
    runner = Runner(project="h2", input_file=str(inp), log=str(work / "h2.log"))
    runner.run()

    # Starts are recorded as files, so they are counted whether the pool runs
    # its workers as threads or as processes.
    marks = work / "started"
    marks.mkdir()

    def failing_worker(task):
        (marks / str(task["idx"])).touch()
        if task["idx"] == 0:
            raise RuntimeError("displacement worker failure")
        time.sleep(0.3)
        return task["idx"], None, "computed", None

    single_point.grad_wrapper = failing_worker
    try:
        single_point.Hessian(runner.mol).numerical_hess()
    except RuntimeError as exc:
        assert "displacement worker failure" in str(exc)
    else:
        raise AssertionError("worker exception did not reach numerical_hess")
    # 12 queued displacements x 0.3 s: an abandoned pool would have started
    # all of them within this wait; a terminated one starts at most the one
    # the worker thread had already taken.
    time.sleep(4.0)
    started = sorted(int(m.name) for m in marks.iterdir())
    assert len(started) <= 2, started
    print("DISPLACEMENT_POOL_OK")
''')


STOP_RUNNING_DRIVER = textwrap.dedent('''
    import os
    import pathlib
    import sys
    import time

    from oqp.pyoqp import Runner
    from oqp.library import single_point

    work = pathlib.Path(sys.argv[1])
    inp = work / "h2.inp"
    inp.write_text("""[input]
    system=
       1   0.000000000   0.000000000   0.000000000
       1   0.000000000   0.000000000   0.740000000
    charge=0
    runtype=energy
    basis=sto-3g
    method=hf
    [guess]
    type=huckel
    [scf]
    multiplicity=1
    type=rhf
    [hess]
    type=numerical
    nproc=2
    """)
    runner = Runner(project="h2", input_file=str(inp), log=str(work / "h2.log"))
    runner.run()

    pid_file = work / "child.pid"
    grandchild_pid_file = work / "grandchild.pid"
    returned = work / "worker1.returned"
    # The child starts a grandchild of its own, as [dftb] backend=probe does,
    # and both sleep: stopping only the direct child would orphan the other.
    child = ("import os, pathlib, subprocess, sys, time; "
             "g = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(120)']); "
             "pathlib.Path(sys.argv[2]).write_text(str(g.pid)); "
             "pathlib.Path(sys.argv[1]).write_text(str(os.getpid())); time.sleep(120)")

    def worker(task):
        if task["idx"] == 1:
            # Stands in for a displaced OpenQP job: a long-running child
            # started through the driver's job registry.
            task["jobs"].run([sys.executable, "-c", child, str(pid_file),
                              str(grandchild_pid_file)], dict(os.environ))
            returned.touch()
            return task["idx"], None, "failed", None
        if task["idx"] == 0:
            deadline = time.time() + 60
            while not pid_file.exists() and time.time() < deadline:
                time.sleep(0.05)
            raise RuntimeError("displacement worker failure")
        return task["idx"], None, "computed", None

    single_point.grad_wrapper = worker
    start = time.time()
    try:
        single_point.Hessian(runner.mol).numerical_hess()
    except RuntimeError as exc:
        assert "displacement worker failure" in str(exc)
    else:
        raise AssertionError("worker exception did not reach numerical_hess")
    elapsed = time.time() - start
    def alive(pid):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        return True

    # Both processes sleep 120 s: an un-stopped one is still alive here, and a
    # driver that does not wait for its workers raises before they return.
    # The grandchild is not ours to reap, so give init a moment to collect it.
    pid = int(pid_file.read_text())
    grandchild = int(grandchild_pid_file.read_text())
    deadline = time.time() + 10
    while alive(grandchild) and time.time() < deadline:
        time.sleep(0.1)
    assert not alive(pid), f"displacement child {pid} outlived the failed driver"
    assert not alive(grandchild), f"grandchild {grandchild} outlived the failed driver"
    assert returned.exists(), "driver raised before its active worker returned"
    assert elapsed < 60, elapsed
    print("DISPLACEMENT_POOL_OK")
''')


CTRL_C_DRIVER = textwrap.dedent('''
    import os
    import pathlib
    import sys
    import time

    from oqp.pyoqp import Runner  # Runner leaves SIGINT at SIG_DFL
    from oqp.library import single_point

    work = pathlib.Path(sys.argv[1])
    inp = work / "h2.inp"
    inp.write_text("""[input]
    system=
       1   0.000000000   0.000000000   0.000000000
       1   0.000000000   0.000000000   0.740000000
    charge=0
    runtype=energy
    basis=sto-3g
    method=hf
    [guess]
    type=huckel
    [scf]
    multiplicity=1
    type=rhf
    [hess]
    type=numerical
    nproc=2
    """)
    runner = Runner(project="h2", input_file=str(inp), log=str(work / "h2.log"))
    runner.run()

    child = ("import os, pathlib, subprocess, sys, time; "
             "g = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(120)']); "
             "pathlib.Path(sys.argv[2]).write_text(str(g.pid)); "
             "pathlib.Path(sys.argv[1]).write_text(str(os.getpid())); time.sleep(120)")

    def worker(task):
        if task["idx"] == 1:
            task["jobs"].run([sys.executable, "-c", child, str(work / "child.pid"),
                              str(work / "grandchild.pid")], dict(os.environ))
        else:
            time.sleep(120)
        return task["idx"], None, "computed", None

    single_point.grad_wrapper = worker
    (work / "driver.pid").write_text(str(os.getpid()))
    single_point.Hessian(runner.mol).numerical_hess()
''')


STOP_ORDER_DRIVER = textwrap.dedent('''
    import os
    import pathlib
    import signal
    import subprocess
    import sys
    import threading
    import time

    from oqp.library import single_point

    work = pathlib.Path(sys.argv[1])
    # The child ignores SIGTERM, so stop() must wait out the grace period and
    # fall back to SIGKILL; below it are a grandchild and a great-grandchild,
    # so freezing has to proceed over more than one generation.
    leaf = "import time; time.sleep(120)"
    middle = ("import pathlib, subprocess, sys, time; "
              f"h = subprocess.Popen([sys.executable, '-c', {leaf!r}]); "
              "pathlib.Path(sys.argv[1]).write_text(str(h.pid)); time.sleep(120)")
    child = ("import os, pathlib, signal, subprocess, sys, time; "
             "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
             f"g = subprocess.Popen([sys.executable, '-c', {middle!r}, sys.argv[3]]); "
             "pathlib.Path(sys.argv[2]).write_text(str(g.pid)); "
             "pathlib.Path(sys.argv[1]).write_text(str(os.getpid())); time.sleep(120)")

    events = []
    real_signal_pids = single_point._signal_pids
    real_wait = subprocess.Popen.wait

    stopped = set()

    def parent_of(pid):
        out = subprocess.run(["ps", "-o", "ppid=", "-p", str(pid)],
                             capture_output=True, text=True).stdout.split()
        return int(out[0]) if out else None

    def recording_signal_pids(pids, sig):
        events.append(("raw-signal", sig))
        if sig == signal.SIGSTOP:
            # Each PID must be a direct child of an already-stopped process,
            # which cannot reap it, so the PID cannot have been recycled.
            for pid in pids:
                assert parent_of(pid) in stopped, ("parent not frozen", pid, stopped)
            stopped.update(pids)
        return real_signal_pids(pids, sig)

    def recording_wait(self, timeout=None):
        if timeout is not None:
            events.append(("grace-wait", None))
        return real_wait(self, timeout=timeout)

    single_point._signal_pids = recording_signal_pids
    subprocess.Popen.wait = recording_wait

    jobs = single_point._ExternalJobs()
    worker = threading.Thread(target=jobs.run, args=(
        [sys.executable, "-c", child, str(work / "child.pid"), str(work / "grandchild.pid"),
         str(work / "great.pid")],
        dict(os.environ)))
    worker.start()
    deadline = time.time() + 60
    while not ((work / "child.pid").exists() and (work / "great.pid").exists()) \
            and time.time() < deadline:
        time.sleep(0.05)
    stopped.add(int((work / "child.pid").read_text()))  # frozen via its Popen handle
    jobs.stop(grace=1.0)
    worker.join(30)

    kinds = [kind for kind, _ in events]
    assert "grace-wait" in kinds, events
    first_wait = kinds.index("grace-wait")
    # A raw PID signalled after a wait could belong to a reaped and recycled
    # process; every raw signal must precede the first wait.
    assert "raw-signal" not in kinds[first_wait:], events

    def alive(pid):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        return True

    tree = [int((work / f).read_text()) for f in ("child.pid", "grandchild.pid", "great.pid")]
    deadline = time.time() + 10
    while any(alive(p) for p in tree) and time.time() < deadline:
        time.sleep(0.1)
    assert not worker.is_alive()
    assert not any(alive(p) for p in tree), tree
    print("DISPLACEMENT_POOL_OK")
''')


@unittest.skipUnless(RUNTIME, "compiled OpenQP runtime not available")
class DisplacementPoolTests(unittest.TestCase):
    def _run(self, driver, timeout):
        with tempfile.TemporaryDirectory(prefix="oqp_disp_pool_") as work:
            env = dict(os.environ)
            env.pop("OPENQP_ROOT", None)
            env.setdefault("OMP_NUM_THREADS", "1")
            try:
                proc = subprocess.run(
                    [sys.executable, "-c", driver, work], cwd=work, env=env,
                    capture_output=True, text=True, timeout=timeout)
            except subprocess.TimeoutExpired as exc:
                self.fail(f"hung for {timeout} s: the displacement workers "
                          "depend on re-importing oqp in a fresh interpreter\n"
                          f"{(exc.stderr or '')[-4000:]!r}")
        self.assertEqual(proc.returncode, 0, proc.stderr[-4000:])
        self.assertIn("DISPLACEMENT_POOL_OK", proc.stdout)

    def test_numerical_hessian_completes_when_workers_cannot_import_oqp(self):
        self._run(HESSIAN_DRIVER, timeout=180)

    def test_numerical_nac_completes_when_workers_cannot_import_oqp(self):
        self._run(NAC_DRIVER, timeout=300)

    def test_pool_runs_unpicklable_work_and_reports_worker_errors(self):
        self._run(UNIT_DRIVER, timeout=120)

    def test_failed_worker_stops_queued_displacements(self):
        self._run(TERMINATE_DRIVER, timeout=120)

    def test_failed_worker_stops_running_displacement_jobs(self):
        self._run(STOP_RUNNING_DRIVER, timeout=180)

    @unittest.skipUnless(os.name == "posix", "process-tree stop is POSIX here")
    def test_stop_signals_no_raw_pid_after_waiting(self):
        self._run(STOP_ORDER_DRIVER, timeout=120)

    @unittest.skipUnless(os.name == "posix", "terminal process groups are POSIX")
    def test_ctrl_c_reaches_running_displacement_jobs(self):
        # A terminal Ctrl-C is SIGINT to the foreground process group.  Runner
        # leaves SIGINT at SIG_DFL, so the driver dies without any cleanup of
        # its own; its displacement jobs must therefore share its group.
        import signal
        import time

        def alive(pid):
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                return False
            return True

        with tempfile.TemporaryDirectory(prefix="oqp_disp_pool_") as work:
            env = dict(os.environ)
            env.pop("OPENQP_ROOT", None)
            env.setdefault("OMP_NUM_THREADS", "1")
            driver = subprocess.Popen([sys.executable, "-c", CTRL_C_DRIVER, work],
                                      cwd=work, env=env, start_new_session=True,
                                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            pids = {}
            try:
                deadline = time.time() + 120
                while time.time() < deadline and len(pids) < 2:
                    for name in ("child", "grandchild"):
                        f = os.path.join(work, f"{name}.pid")
                        if name not in pids and os.path.exists(f) and open(f).read():
                            pids[name] = int(open(f).read())
                    time.sleep(0.1)
                self.assertEqual(len(pids), 2, "displacement job never started")
                os.killpg(driver.pid, signal.SIGINT)
                driver.wait(timeout=30)
                deadline = time.time() + 15
                while any(alive(p) for p in pids.values()) and time.time() < deadline:
                    time.sleep(0.1)
                survivors = {k: v for k, v in pids.items() if alive(v)}
                self.assertEqual(survivors, {},
                                 "Ctrl-C left displacement processes running")
            finally:
                for pid in list(pids.values()) + [driver.pid]:
                    try:
                        os.kill(pid, signal.SIGKILL)
                    except (ProcessLookupError, PermissionError):
                        pass


if __name__ == "__main__":
    unittest.main()
