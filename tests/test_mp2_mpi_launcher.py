"""The MPI regression launcher must work without initializing parent MPI."""
import builtins
from pathlib import Path
import runpy
import shutil
import subprocess

import pytest


def test_parent_reaches_worker_without_importing_mpi(monkeypatch, tmp_path):
    parent = runpy.run_path(str(Path(__file__).with_name('mp2_mpi_regression.py')))
    original_import = builtins.__import__

    def reject_mpi(name, *args, **kwargs):
        if name == 'mpi4py' or name.startswith('mpi4py.'):
            raise AssertionError('launcher imported MPI before starting a worker')
        return original_import(name, *args, **kwargs)

    class WorkerLaunched(Exception):
        pass

    commands = []

    def run(command, **kwargs):
        commands.append(command)
        if command[0] == 'cc':
            return subprocess.CompletedProcess(command, 0)
        assert command[0] == '/test/mpiexec'
        assert command[1:4] == ['--oversubscribe', '-n', '1']
        assert command[6] == '--memory-worker'
        raise WorkerLaunched

    monkeypatch.setattr(builtins, '__import__', reject_mpi)
    monkeypatch.setattr(shutil, 'which', lambda name: '/test/mpiexec')
    monkeypatch.setattr(subprocess, 'run', run)
    with pytest.raises(WorkerLaunched):
        parent['main'](tmp_path / 'parent')
    assert len(commands) == 2
