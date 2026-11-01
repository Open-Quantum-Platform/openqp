"""conftest.py refuses ``--forked`` in a session that ``import oqp`` puts under MPI.

``oqp.utils.mpi_utils`` calls MPI_Init at import, so with mpi4py installed the
pytest process that pytest-forked forks once per test is already an MPI
process.  On macOS that MPI_Init also starts libdispatch, and Accelerate runs
its threaded BLAS through dispatch_apply, which a forked child cannot use:
75 tests in 28 files died with SIGSEGV inside libBLAS, all of them passing
without --forked.  The refusal turns that into one usage error naming the
cause.  These tests pin both sides of the decision and the import it mirrors.
"""
import importlib.util
import os
import re
import sys
import types

import pytest

import conftest

_MPI_UTILS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "pyoqp", "oqp", "utils", "mpi_utils.py")


class _Config:
    def __init__(self, forked):
        self._forked = forked

    def getoption(self, name, default=None):
        assert name == "forked"
        return self._forked


def test_forked_session_is_refused_when_import_oqp_initialises_mpi(monkeypatch):
    monkeypatch.setattr(conftest, "_import_oqp_initialises_mpi", lambda: True)
    with pytest.raises(pytest.UsageError) as refusal:
        conftest.pytest_configure(_Config(forked=True))
    message = str(refusal.value)
    assert "MPI_Init" in message and "fork()" in message
    assert "dispatch_apply" in message


def test_forked_session_without_mpi_is_allowed(monkeypatch):
    monkeypatch.setattr(conftest, "_import_oqp_initialises_mpi", lambda: False)
    conftest.pytest_configure(_Config(forked=True))


def test_unforked_session_never_probes_mpi(monkeypatch):
    # The probe imports mpi4py.MPI, i.e. calls MPI_Init; an unforked session
    # must not pay for that or have its MPI state touched by the conftest.
    def probe():
        raise AssertionError("MPI probed in an unforked session")
    monkeypatch.setattr(conftest, "_import_oqp_initialises_mpi", probe)
    conftest.pytest_configure(_Config(forked=False))


class _MissingLibmpiFinder:
    """Meta-path finder that fails like mpi4py 4.x without an MPI library.

    mpi4py >= 4 imports ``mpi4py.MPI`` through ``mpi4py._mpiabi._Finder``,
    whose ``find_spec`` dlopens libmpi and raises
    ``RuntimeError("cannot load MPI library ...")`` when it is absent
    (``_mpiabi._dlopen_libmpi``; reproduced with mpi4py 4.1.2 on zeus via
    ``MPI4PY_LIBMPI=/nonexistent/libmpi.dylib``).
    """

    def __init__(self, error, message):
        self.error = error
        self.message = message

    def find_spec(self, fullname, path, target=None):
        if fullname == "mpi4py.MPI":
            raise self.error(self.message)
        return None


_NO_LIBMPI = "cannot load MPI library\ndlopen(libmpi.dylib): no such file"


def _install_failing_mpi4py(monkeypatch, error, message=_NO_LIBMPI):
    """Make ``from mpi4py import MPI`` raise ``error``; monkeypatch restores it."""
    package = types.ModuleType("mpi4py")
    package.__path__ = []
    monkeypatch.setitem(sys.modules, "mpi4py", package)
    monkeypatch.delitem(sys.modules, "mpi4py.MPI", raising=False)
    monkeypatch.setattr(sys, "meta_path", [_MissingLibmpiFinder(error, message)] + sys.meta_path)


_NO_MPI_ERRORS = [ImportError, OSError, RuntimeError]


@pytest.mark.parametrize("error", _NO_MPI_ERRORS)
def test_probe_treats_a_missing_mpi_like_mpi_utils_does(monkeypatch, error):
    _install_failing_mpi4py(monkeypatch, error)
    assert conftest._import_oqp_initialises_mpi() is False


def test_forked_session_without_libmpi_is_allowed(monkeypatch):
    # mpi4py installed, MPI library missing: a --forked run must start, not
    # end in INTERNALERROR from the probe.
    _install_failing_mpi4py(monkeypatch, RuntimeError)
    conftest.pytest_configure(_Config(forked=True))


@pytest.mark.parametrize("error", _NO_MPI_ERRORS)
def test_mpi_utils_runs_serial_when_mpi_is_missing(monkeypatch, error):
    # The probe is only right if ``import oqp`` survives the same failure.
    _install_failing_mpi4py(monkeypatch, error)
    spec = importlib.util.spec_from_file_location("_mpi_utils_without_mpi", _MPI_UTILS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.mpi4py_available is False


# RuntimeErrors that mpi4py raises only after libmpi has loaded: they mean a
# broken MPI, not an absent one, and must not turn an mpiexec rank serial.
_MPI_FAILURES_AFTER_LOAD = [
    "MPI_Abi_get_version [ierr=1]",
    "MPI_Init_thread failed: error code 15",
]


@pytest.mark.parametrize("message", _MPI_FAILURES_AFTER_LOAD)
def test_probe_propagates_other_mpi_runtime_errors(monkeypatch, message):
    _install_failing_mpi4py(monkeypatch, RuntimeError, message)
    with pytest.raises(RuntimeError, match=re.escape(message)):
        conftest._import_oqp_initialises_mpi()


@pytest.mark.parametrize("message", _MPI_FAILURES_AFTER_LOAD)
def test_mpi_utils_propagates_other_mpi_runtime_errors(monkeypatch, message):
    _install_failing_mpi4py(monkeypatch, RuntimeError, message)
    spec = importlib.util.spec_from_file_location("_mpi_utils_broken_mpi", _MPI_UTILS)
    module = importlib.util.module_from_spec(spec)
    with pytest.raises(RuntimeError, match=re.escape(message)):
        spec.loader.exec_module(module)


def test_probe_mirrors_the_import_in_mpi_utils():
    # If mpi_utils stops initialising MPI at import (or catches other errors),
    # the probe above no longer describes what collection does.
    with open(_MPI_UTILS, encoding="utf-8") as handle:
        source = handle.read()
    assert re.search(r"^try:\n\s+from mpi4py import MPI\n", source, re.M)
    assert re.search(r"^except \(ImportError, OSError\)", source, re.M)
    assert re.search(r'^_NO_LIBMPI_MESSAGE = "cannot load MPI library"$', source, re.M)
    assert re.search(r"^except RuntimeError as \w+:\n\s+if not str\(\w+\)\."
                     r"startswith\(_NO_LIBMPI_MESSAGE\):\n\s+raise\n", source, re.M)


_GITLAB_CI = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          ".gitlab-ci.yml")


def _gitlab_jobs(text):
    """Split .gitlab-ci.yml into {top-level key: block text}."""
    jobs, name, lines = {}, None, []
    for line in text.splitlines():
        match = re.match(r"^([.\w-]+):\s*$", line)
        if match:
            if name is not None:
                jobs[name] = "\n".join(lines)
            name, lines = match.group(1), []
        elif name is not None:
            lines.append(line)
    if name is not None:
        jobs[name] = "\n".join(lines)
    return jobs


def _script_items(block):
    """Each ``- ...`` item of a job, folded continuation lines joined."""
    items, current = [], None
    for line in block.splitlines():
        if re.match(r"^\s{4}- ", line):
            if current is not None:
                items.append(current)
            current = line.strip()[2:]
        elif current is not None and line.startswith("      "):
            current += " " + line.strip()
        else:
            if current is not None:
                items.append(current)
            current = None
    if current is not None:
        items.append(current)
    return items


def test_gitlab_forked_pytest_jobs_install_no_mpi4py():
    # The conftest refusal above would stop such a job before any test runs.
    # A job that runs pytest --forked must keep mpi4py out of its venv, as
    # the GitHub pytest job does; MPI runs belong in their own job.
    with open(_GITLAB_CI, encoding="utf-8") as handle:
        jobs = _gitlab_jobs(handle.read())
    forked = [name for name, block in jobs.items()
              if re.search(r"-m pytest\b[^\n]*(\n[^\n]*){0,2}--forked", block)]
    assert forked, "no GitLab job runs pytest --forked; update this test"
    for name in forked:
        for item in _script_items(jobs[name]):
            if "pip install" not in item:
                continue
            assert "mpi4py" not in item, (name, item)
            assert "pyoqp/requirements.txt" not in item, (name, item)
