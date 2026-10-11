"""Real two-rank MP2 failure consensus; run by CI's MPI-enabled leg.

Run as ``python tests/mp2_mpi_regression.py OUTPUT_DIRECTORY``. The parent
enforces both OpenMP widths and a timeout; only rank 1 arms the native shim.
"""
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def worker(directory, stage, library):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    assert comm.Get_size() == 2
    tests = runpy.run_path(str(ROOT / 'test_mp2_dispatch.py'))
    tests['_allocation_child'](Path(directory), stage, library, comm=comm, fault_rank=1)
    comm.Barrier()
    if comm.Get_rank() == 0:
        print('MULTIRANK_ONE_RANK_FAILURE_PASS', flush=True)


def main(directory):
    from mpi4py import MPI
    assert MPI.COMM_WORLD.Get_size() == 1, 'launch the parent outside mpiexec'
    launcher = shutil.which('mpiexec') or shutil.which('mpirun')
    assert launcher, 'MPI launcher required'
    assert sys.platform in ('linux', 'darwin'), 'native interposer platform unsupported'
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    library = directory / ('mp2_alloc_fail.dylib' if sys.platform == 'darwin' else 'mp2_alloc_fail.so')
    flags = ['-dynamiclib'] if sys.platform == 'darwin' else ['-shared', '-fPIC']
    subprocess.run(['cc', *flags, '-O2', str(ROOT / 'native/mp2_alloc_fail.c'),
                    '-o', str(library)], check=True)
    for width in (1, 4):
        for stage in ('full_scratch', 'alpha_batched_scratch',
                      'after_alpha_same_spin', 'beta_store', 'beta_scratch'):
            case = directory / f'{stage}-{width}'
            case.mkdir()
            env = os.environ.copy()
            env['OMP_NUM_THREADS'] = str(width)
            env['MKL_NUM_THREADS'] = env['OPENBLAS_NUM_THREADS'] = '1'
            env['DYLD_INSERT_LIBRARIES' if sys.platform == 'darwin' else 'LD_PRELOAD'] = str(library)
            result = subprocess.run([launcher, '--oversubscribe', '-n', '2', sys.executable,
                                     str(Path(__file__).resolve()), '--worker', str(case),
                                     stage, str(library)], env=env, text=True,
                                    capture_output=True, timeout=900)
            (case / 'launcher.log').write_text(result.stdout + result.stderr)
            assert result.returncode == 0, result.stdout + result.stderr
            assert 'MULTIRANK_ONE_RANK_FAILURE_PASS' in result.stdout
            print(f'{stage}: 2 ranks x {width} threads PASS', flush=True)


if __name__ == '__main__':
    if sys.argv[1] == '--worker':
        worker(*sys.argv[2:])
    else:
        main(sys.argv[1])
