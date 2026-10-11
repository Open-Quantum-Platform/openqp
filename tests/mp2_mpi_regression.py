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
    tests['_assert_bounds_runtime']()
    tests['_allocation_child'](Path(directory), stage, library, comm=comm, fault_rank=1)
    comm.Barrier()
    if comm.Get_rank() == 0:
        print('MULTIRANK_ONE_RANK_FAILURE_PASS', flush=True)


def memory_worker(directory, case):
    """Pin a node budget on opposite sides of a real replicated allocation."""
    from mpi4py import MPI
    import numpy as np
    import oqp
    comm = MPI.COMM_WORLD
    # Shared-memory splits can return singleton groups on one physical node.
    # Establish this fixture's co-location independently; production budgets
    # count matching processor identities rather than trusting that split.
    processors = comm.allgather(MPI.Get_processor_name())
    assert len(set(processors)) == 1, 'test requires one processor identity'
    tests = runpy.run_path(str(ROOT / 'test_mp2_dispatch.py'))
    tests['_assert_bounds_runtime']()
    runner, log = tests['_run'](Path(directory), 'memory', method='hf', usempi=True)
    reference_energy = runner.mol.data._data.mol_energy.energy
    os.environ['OQP_MP2_MAX_JBUILDS'] = '4000'
    oqp.mp2_energy(runner.mol)
    reference = tests['_components'](log) if comm.Get_rank() == 0 else None
    reference = comm.bcast(reference, root=0)
    nbf = int(runner.mol.data.get_basis()['nbf'])
    assert nbf < 120
    npair = nbf * (nbf + 1) // 2
    nt = int(os.environ['OMP_NUM_THREADS'])
    packed_bytes = 8 * npair * (npair + 1) // 2
    full_bytes = 8 * (nbf**4 + npair**2 + 3 * nbf**2 * nt) + 16*npair + packed_bytes
    # At this fixed node budget one full transform fits; two copies do not.
    budget = 1.5 * full_bytes / (0.25 * 0.9 * 2**30)
    if case == 'refuse':
        nocc, nvir = 5, nbf - 5
        first = 8 * nt * (nbf**2 + nbf*nvir + nocc*nvir)
        second = 8 * nt * (nbf**2 + nbf*nocc + nocc*nvir) + 8*(nocc*nvir)**2
        batch_bytes = packed_bytes + 8*npair*nocc*nvir + max(first, second) + 16*npair
        # One rank fits this node budget, but the two-rank replicated store
        # must be refused. Cap=0 prevents an unaffordable direct-J fallback.
        assert comm.Get_size() == 2
        budget = 1.5 * batch_bytes / (0.9 * 2**30)
    runner.mol.data._data.mol_energy.energy = reference_energy
    runner.mol.data._data.mol_energy.etot = reference_energy
    if comm.Get_rank() == 0:
        log.write_text('')
    comm.Barrier()
    os.environ['OQP_MEMORY_LIMIT_GB'] = str(budget)
    os.environ['OQP_MP2_MAX_JBUILDS'] = '0'
    oqp.mp2_energy(runner.mol)
    assert case != 'refuse', 'two-rank allocation exceeded the node budget'
    text = comm.bcast(log.read_text() if comm.Get_rank() == 0 else None, root=0)
    expected = 'full' if comm.Get_size() == 1 else 'occupied-virtual'
    assert 'MP2 path: n5 ' + expected in text, text
    if expected == 'occupied-virtual':
        assert 'MP2 path: n5 full' not in text, text
    values = comm.bcast(tests['_components'](log) if comm.Get_rank() == 0 else None, root=0)
    np.testing.assert_allclose(values, reference, rtol=0, atol=1e-9)
    np.testing.assert_allclose(comm.allgather(float(runner.mol.data._data.mol_energy.energy)),
                               reference[-1], rtol=0, atol=1e-9)
    if comm.Get_rank() == 0:
        print('NODE_MEMORY_DISPATCH_PASS', flush=True)


def main(directory):
    # Only workers initialize MPI. Singleton MPI_Init in this launcher can
    # stall before mpiexec, and initialized parents cannot safely spawn on
    # every supported MPI runtime.
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
        for ranks, case_name in ((1, 'select'), (2, 'select'), (2, 'refuse')):
            case = directory / f'memory-{case_name}-{ranks}-{width}'
            case.mkdir()
            env = os.environ.copy()
            env['OMP_NUM_THREADS'] = str(width)
            env['MKL_NUM_THREADS'] = env['OPENBLAS_NUM_THREADS'] = '1'
            result = subprocess.run([launcher, '--oversubscribe', '-n', str(ranks), sys.executable,
                                     str(Path(__file__).resolve()), '--memory-worker', str(case),
                                     case_name], env=env, text=True, capture_output=True, timeout=180)
            (case / 'launcher.log').write_text(result.stdout + result.stderr)
            if case_name == 'refuse':
                assert result.returncode != 0, result.stdout + result.stderr
                text = (case / 'memory.log').read_text()
                assert 'discard spin trials before direct-J fallback' in text, text
                assert 'MP2 not computed' in text, text
                assert 'MP2 path: n5 ' not in text, text
            else:
                assert result.returncode == 0, result.stdout + result.stderr
                assert 'NODE_MEMORY_DISPATCH_PASS' in result.stdout
            print(f'node memory {case_name}: {ranks} ranks x {width} threads PASS', flush=True)
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
    elif sys.argv[1] == '--memory-worker':
        memory_worker(*sys.argv[2:])
    else:
        main(sys.argv[1])
