"""Execute MP2 transforms, their spin components, and resource fallbacks.

Run this file in separate processes with OMP_NUM_THREADS=1 and 4. No backend
skip is allowed: a missing MP2 engine is a failed test, not coverage.
"""
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import oqp
from oqp.pyoqp import Runner

WATER = 'O 0 0 0; H 0 .757 .587; H 0 -.757 .587'
OH = 'O 0 0 0; H 0 0 .97'
COMPONENTS = r'E\(MP2, (?:same-spin aa|same-spin bb|opp-spin  ab|correlation|total)\)\s*=\s*([-+0-9.eE]+)'


def _run(tmp_path, name, scf='rhf', mult=1, atoms=WATER, variant='mp2',
         basis='cc-pvdz', method='mp2', usempi=False):
    assert oqp.BACKEND_AVAILABLE and hasattr(oqp, 'mp2_energy')
    log = tmp_path / (name + '.log')
    config = {
        'input': {'system': atoms.replace(';', '\n'), 'method': method,
                  'basis': basis, 'runtype': 'energy'},
        'guess': {'type': 'hcore', 'save_mol': 'False'},
        'scf': {'type': scf, 'multiplicity': str(mult), 'conv': '1.0e-11',
                'maxit': '200', 'save_molden': 'False'},
        'mp2': {'variant': variant},
        'tests': {'exception': 'True'},
    }
    runner = Runner(project=name, input_file=None, log=str(log),
                    input_dict=config, silent=1, usempi=usempi)
    runner.run()
    return runner, log


def _components(log):
    text = log.read_text()
    values = np.asarray([float(x) for x in re.findall(COMPONENTS, text)])
    assert len(values) == 5, text
    assert np.all(np.isfinite(values))
    return values


@pytest.mark.parametrize('scf,mult,atoms', [
    ('rhf', 1, WATER), ('uhf', 2, OH), ('rohf', 2, OH),
])
@pytest.mark.parametrize('variant', ['mp2', 'os-mp2', 'ss-mp2'])
def test_n5_matches_direct_j(tmp_path, monkeypatch, scf, mult, atoms, variant):
    results = []
    for name, cap, budget in [('direct', '4000', '16'),
                              ('n5', '0', '16'), ('batched', '0', '.005')]:
        monkeypatch.setenv('OQP_MP2_MAX_JBUILDS', cap)
        # Apply the MP2-only budget after SCF, so this tests the transform
        # selection rather than an unrelated SCF memory check.
        runner, log = _run(tmp_path, name, scf, mult, atoms, variant, method='hf')
        runner.mol.data.set_mp2_same_spin_scale(0 if variant == 'os-mp2' else 1)
        runner.mol.data.set_mp2_opposite_spin_scale(0 if variant == 'ss-mp2' else 1)
        monkeypatch.setenv('OQP_MEMORY_LIMIT_GB', budget)
        oqp.mp2_energy(runner.mol)
        text = log.read_text()
        if cap == '4000':
            assert 'MP2 path: direct-J' in text
        elif name == 'batched' or (scf != 'rhf' and variant != 'ss-mp2'):
            assert 'MP2 path: n5 occupied-virtual' in text
        else:
            assert 'MP2 path: n5 full' in text
        results.append(_components(log))
        monkeypatch.delenv('OQP_MEMORY_LIMIT_GB')
    np.testing.assert_allclose(results[0], results[1], rtol=0, atol=1e-9)
    np.testing.assert_allclose(results[0], results[2], rtol=0, atol=1e-9)


def test_input_dispatch_and_default_crossover(tmp_path, monkeypatch):
    """Use the actual default threshold, with 124 spherical AOs (>120).

    The small direct-J result pins the input route too. The large full and
    occupied-virtual calculations share the same converged SCF reference.
    """
    monkeypatch.delenv('OQP_MP2_MAX_JBUILDS', raising=False)
    monkeypatch.delenv('OQP_MEMORY_LIMIT_GB', raising=False)
    _, small_log = _run(tmp_path, 'small_input')
    assert 'MP2 path: direct-J' in small_log.read_text()
    runner, log = _run(tmp_path, 'large_input', atoms='H 0 0 0; H 0 0 .74', basis='d-aug-cc-pvqz')
    assert int(runner.mol.data.get_basis()['nbf']) > 120
    assert 'MP2 path: n5 ' in log.read_text()
    assert 'MP2 path: direct-J' not in log.read_text()
    reference = _components(log)
    # PySCF RHF/MP2 with the identical Basis Set Exchange d-aug-cc-pVQZ
    # basis, 0.74 Angstrom, SCF conv_tol=1e-12, conv_tol_grad=1e-9.
    np.testing.assert_allclose(reference[-2:], [-0.03325695678608737,
                                               -1.1667414879132405], rtol=0, atol=1e-9)
    # Restore the SCF energy before evaluating MP2 on the same native handle.
    scf_energy = float(re.findall(r'E\(reference, SCF\)\s*=\s*([-+0-9.eE]+)', log.read_text())[-1])
    runner.mol.data._data.mol_energy.energy = scf_energy
    runner.mol.data._data.mol_energy.etot = scf_energy
    monkeypatch.setenv('OQP_MEMORY_LIMIT_GB', '2')
    log.write_text('')
    oqp.mp2_energy(runner.mol)
    assert 'MP2 path: n5 occupied-virtual' in log.read_text()
    np.testing.assert_allclose(reference, _components(log), rtol=0, atol=1e-9)


def test_empty_beta_space_and_dimension_reuse(tmp_path, monkeypatch):
    monkeypatch.setenv('OQP_MP2_MAX_JBUILDS', '0')
    # Increasing then decreasing dimensions in one Python process exercises
    # repeated native creation/use/destruction and the empty beta space.
    for name, scf, mult, atoms, basis in [
        ('hydrogen1', 'uhf', 2, 'H 0 0 0', 'sto-3g'),
        ('water', 'rhf', 1, WATER, 'cc-pvdz'),
        ('hydrogen2', 'uhf', 2, 'H 0 0 0', 'cc-pvdz'),
    ]:
        _, log = _run(tmp_path, name, scf, mult, atoms, basis=basis)
        values = _components(log)
        if name.startswith('hydrogen'):
            np.testing.assert_allclose(values[:4], 0, rtol=0, atol=1e-12)


@pytest.mark.parametrize('stage', ['full_scratch', 'alpha_batched_scratch',
                                   'after_alpha_same_spin', 'beta_store', 'beta_scratch'])
def test_allocator_failure_runs_real_fallback(tmp_path, stage):
    """Fail malloc inside a native ALLOCATE, then execute the real fallback.

    The interposer is armed only after SCF. It returns NULL once, rather than
    mocking a Fortran result or asserting source text. The child must report
    the actual allocation failure and the subsequent production path.
    """
    fixture = Path(__file__).with_name('native') / 'mp2_alloc_fail.c'
    assert sys.platform in ('darwin', 'linux'), 'allocator fixture supports macOS and Linux/glibc'
    library = tmp_path / ('mp2_alloc_fail.dylib' if sys.platform == 'darwin' else 'mp2_alloc_fail.so')
    flags = ['-dynamiclib'] if sys.platform == 'darwin' else ['-shared', '-fPIC']
    subprocess.run(['cc', *flags, '-O2', str(fixture), '-o', str(library)], check=True)
    env = os.environ.copy()
    # The native scratch shape uses the OpenMP team width. Keep an explicit
    # caller's width, and make the child and the failure size agree by default.
    env.setdefault('OMP_NUM_THREADS', '1')
    env['DYLD_INSERT_LIBRARIES' if sys.platform == 'darwin' else 'LD_PRELOAD'] = str(library)
    child = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                            str(tmp_path), stage, str(library)], env=env,
                           text=True, capture_output=True, timeout=900)
    assert child.returncode == 0, child.stdout + child.stderr
    assert 'REAL_ALLOCATE_FAILURE_AND_FALLBACK_PASS' in child.stdout


def _native_thread_controls():
    import ctypes
    from oqp.runtime import library_path, resolve_oqp_root
    root, suffix = resolve_oqp_root()
    native = ctypes.CDLL(str(library_path(root, suffix)))
    native.oqp_blas_thread_count.restype = ctypes.c_int64
    native.oqp_blas_thread_set.argtypes = [ctypes.c_int64]
    native.omp_get_max_threads.restype = ctypes.c_int
    native.oqp_omp_set_num_threads.argtypes = [ctypes.c_int]
    return native


def _counts(native):
    return native.oqp_blas_thread_count(), native.omp_get_max_threads()


@pytest.mark.parametrize('width', [1, 4])
def test_thread_counts_are_restored_after_repeated_transforms(tmp_path, width):
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(width)
    env['MKL_NUM_THREADS'] = env['OPENBLAS_NUM_THREADS'] = '2'
    child = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                            str(tmp_path), 'threads', str(width)], env=env,
                           text=True, capture_output=True, timeout=300)
    assert child.returncode == 0, child.stdout + child.stderr
    assert 'THREAD_RESTORATION_PASS' in child.stdout


def _thread_child(tmp_path, width):
    native = _native_thread_controls()
    for scf, mult, atoms in [('rhf', 1, WATER), ('uhf', 2, OH)]:
        runner, log = _run(tmp_path, 'threads_' + scf, scf, mult, atoms, method='hf')
        energy = runner.mol.data._data.mol_energy.energy
        for budget in ['16', '.005', '16']:
            os.environ['OQP_MP2_MAX_JBUILDS'] = '0'
            os.environ['OQP_MEMORY_LIMIT_GB'] = budget
            native.oqp_blas_thread_set(2)
            native.oqp_omp_set_num_threads(width)
            before = _counts(native)
            assert before[1] == width
            runner.mol.data._data.mol_energy.energy = energy
            runner.mol.data._data.mol_energy.etot = energy
            oqp.mp2_energy(runner.mol)
            assert _counts(native) == before
    print('THREAD_RESTORATION_PASS')


def _allocation_child(tmp_path, stage, library, comm=None, fault_rank=1):
    import ctypes
    interposer = ctypes.CDLL(library)
    interposer.mp2_test_arm.argtypes = [ctypes.c_size_t, ctypes.c_int, ctypes.c_size_t]
    interposer.mp2_test_failures.restype = ctypes.c_int
    rank = 0 if comm is None else comm.Get_rank()
    native = _native_thread_controls()
    def components():
        value = _components(log) if rank == 0 else None
        return value if comm is None else comm.bcast(value, root=0)
    scf, mult, atoms = ('rhf', 1, WATER) if stage == 'full_scratch' else ('uhf', 2, 'Li 0 0 0')
    os.environ['OQP_MP2_MAX_JBUILDS'] = '4000'
    runner, log = _run(tmp_path, 'alloc_' + stage, scf, mult, atoms, method='hf', usempi=comm is not None)
    scf_energy = runner.mol.data._data.mol_energy.energy
    runner.mol.data.set_mp2_same_spin_scale(1)
    runner.mol.data.set_mp2_opposite_spin_scale(1)
    oqp.mp2_energy(runner.mol)
    reference = components()
    nbf = int(runner.mol.data.get_basis()['nbf'])
    npair = nbf * (nbf + 1) // 2
    nocca, noccb = 2, 1
    nt = int(os.environ.get('OMP_NUM_THREADS', '1'))
    # Full scratch failure falls back to occupied-virtual. Other failures
    # are induced at n5 even with an affordable direct-J fallback by using
    # the real >120 default crossover on a 127-AO Li basis.
    if stage != 'full_scratch':
        # The 127-AO default crossover leaves the direct-J guard affordable.
        # Li has nonzero alpha SS and OS terms, with an actual beta trial.
        runner, log = _run(tmp_path, 'alloc_large_' + stage, scf, mult, atoms,
                           basis='aug-cc-pv5z', method='hf', usempi=comm is not None)
        nbf = int(runner.mol.data.get_basis()['nbf'])
        assert nbf > 120
        npair = nbf * (nbf + 1) // 2
        scf_energy = runner.mol.data._data.mol_energy.energy
        # First run supplies the no-failure n5 reference. The independent
        # production-input pin and baseline gate check its numerical value.
        runner.mol.data.set_mp2_same_spin_scale(1)
        runner.mol.data.set_mp2_opposite_spin_scale(1)
        before = _counts(native)
        oqp.mp2_energy(runner.mol)
        assert _counts(native) == before
        reference = components()
    runner.mol.data._data.mol_energy.energy = scf_energy
    runner.mol.data._data.mol_energy.etot = scf_energy
    if rank == 0:
        log.write_text('')
    if comm is not None:
        comm.Barrier()
    if stage == 'full_scratch':
        # The packed full-half array is the first allocation of this size.
        os.environ['OQP_MP2_MAX_JBUILDS'] = '0'
        size, occurrence = 8*npair*npair, 1
    elif stage == 'alpha_batched_scratch':
        size, occurrence = 8*npair*nocca*(nbf-nocca), 1
    elif stage == 'after_alpha_same_spin':
        # Opposite-spin second-half scratch follows same-spin accumulation.
        size, occurrence = 8*noccb*nbf*nt, 1
    elif stage == 'beta_store':
        size, occurrence = 8*npair*(npair+1)//2, 2
    else:
        size, occurrence = 8*npair*noccb*(nbf-noccb), 1
        os.environ['OQP_MEMORY_LIMIT_GB'] = '2'
    native.oqp_blas_thread_set(2)
    native.oqp_omp_set_num_threads(nt)
    before = _counts(native)
    if comm is None or rank == fault_rank:
        interposer.mp2_test_arm(size, occurrence, 8*npair*(npair+1)//2)
    oqp.mp2_energy(runner.mol)
    assert _counts(native) == before
    failures = interposer.mp2_test_failures()
    if comm is None:
        assert failures == 1
    else:
        assert comm.allgather(failures) == [int(r == fault_rank) for r in range(comm.Get_size())]
        totals = comm.allgather(float(runner.mol.data._data.mol_energy.energy))
        np.testing.assert_allclose(totals, reference[-1], rtol=0, atol=1e-9)
    text = log.read_text() if rank == 0 else None
    if comm is not None:
        text = comm.bcast(text, root=0)
    expected = 'MP2 path: n5 occupied-virtual' if stage == 'full_scratch' else 'MP2 path: direct-J'
    assert expected in text, text
    if stage != 'full_scratch':
        assert 'discard spin trials before direct-J fallback' in text
    np.testing.assert_allclose(reference, components(), rtol=0, atol=1e-9)
    print('REAL_ALLOCATE_FAILURE_AND_FALLBACK_PASS')



@pytest.mark.parametrize('scf,mult,atoms', [('rhf', 1, WATER), ('uhf', 2, OH), ('rohf', 2, OH)])
def test_n5_pin_on_production_integrals(tmp_path, monkeypatch, scf, mult, atoms):
    """Independent NumPy contraction of the engine's own AO/MO inputs.

    The UHF alpha/beta occupied-virtual tensor has distinct pair dimensions;
    exchanging its free indices cannot hide behind spatial tensor symmetry.
    """
    monkeypatch.setenv('OQP_MP2_MAX_JBUILDS', '0')
    runner, log = _run(tmp_path, 'kernel_pin', scf, mult, atoms)
    result = _components(log)
    mol = runner.mol
    nbf = int(mol.data.get_basis()['nbf'])
    na, nb = (5, 5) if scf == 'rhf' else (5, 4)

    def semicanonical(spin, nocc):
        coeff = np.array(mol.data['OQP::VEC_MO_' + spin], copy=True).reshape(nbf, nbf).T
        packed = np.array(mol.data['OQP::FOCK_' + spin], copy=True)
        fock = np.zeros((nbf, nbf))
        rows, cols = np.tril_indices(nbf)
        fock[rows, cols] = packed
        fock[cols, rows] = packed
        fmo = coeff.T @ fock @ coeff
        eo, uo = np.linalg.eigh(fmo[:nocc, :nocc])
        ev, uv = np.linalg.eigh(fmo[nocc:, nocc:])
        semimo = np.column_stack((coeff[:, :nocc] @ uo, coeff[:, nocc:] @ uv))
        energies = np.concatenate((eo, ev))
        fov = (semimo.T @ fock @ semimo)[:nocc, nocc:]
        singles = np.sum(fov*fov / (eo[:, None]-ev[None, :]))
        return semimo, energies, singles

    ca, ea, sa = semicanonical('A', na)
    cb, eb, sb = semicanonical('A' if scf == 'rhf' else 'B', nb)
    # Own copies above must precede a new tag allocation: no borrowed NumPy
    # view is retained across native backing-record replacement.
    oqp.int2e(mol)
    eri = np.array(mol.data['OQP::ERI_AO'], copy=True).reshape((nbf,)*4)

    def ovov(cl, nl, cr, nr):
        return np.einsum('ui,va,wj,xb,uvwx->iajb', cl[:, :nl], cl[:, nl:],
                         cr[:, :nr], cr[:, nr:], eri, optimize=True)

    def denominator(el, nl, er, nr):
        return (el[:nl, None, None, None] - el[None, nl:, None, None] +
                er[None, None, :nr, None] - er[None, None, None, nr:])

    aa, bb, ab = ovov(ca, na, ca, na), ovov(cb, nb, cb, nb), ovov(ca, na, cb, nb)
    eaa = .25*np.sum((aa-aa.transpose(0, 3, 2, 1))**2 / denominator(ea, na, ea, na))
    ebb = .25*np.sum((bb-bb.transpose(0, 3, 2, 1))**2 / denominator(eb, nb, eb, nb))
    eab = np.sum(ab*ab / denominator(ea, na, eb, nb))
    if scf != 'rhf':
        assert ab.shape[0] != ab.shape[2]
        # Confirm that a one-spin OS transform gives a different energy.
        wrong = np.sum(ovov(ca, na, ca, nb)**2 / denominator(ea, na, eb, nb))
        assert abs(wrong-eab) > 1e-5
    np.testing.assert_allclose(result[:4], [eaa, ebb, eab, eaa+ebb+eab+sa+sb], rtol=0, atol=1e-9)

if __name__ == '__main__':
    if sys.argv[2] == 'threads':
        _thread_child(Path(sys.argv[1]), int(sys.argv[3]))
    else:
        _allocation_child(Path(sys.argv[1]), sys.argv[2], sys.argv[3])
