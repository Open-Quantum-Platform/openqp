"""Pin every SC-NEVPT2 subspace at its native boundary.

Unequal core/active/virtual dimensions and nonsymmetric operands expose axis
errors. Compare norm and Hamiltonian arrays before division by denominators,
including single-orbital diagonals and repeated calls with changing dimensions.
"""
import numpy as np
import pytest

from oqp.library import nevpt2_sc as sc
from oqp.library.fci import _determinants


SUBSPACES = ('Sr', 'Si', 'Sijrs', 'Sijr', 'Srsi', 'Srs', 'Sij', 'Sir')


def _python(monkeypatch, fn):
    with monkeypatch.context() as patch:
        patch.setattr(sc, '_koopmans_lib', lambda: None)
        return fn()


def _subspace_inputs(ncore, nact, nvirt, seed):
    rng = np.random.default_rng(seed)
    nbf = ncore + nact + nvirt
    h = rng.normal(size=(nbf, nbf)) * 0.1
    g = rng.normal(size=(nbf,) * 4) * 0.1
    eps = np.linspace(-3.0, 3.0, nbf)
    B = sc._blocks(h, g, ncore, nact, eps)
    dms = tuple(rng.normal(size=(nact,) * (2 * rank)) * 0.1
                for rank in range(1, 5))
    return B, dms


def _call_subspace(name, B, dms):
    d1, d2, d3, d4 = dms
    h, g, ec, ev = (B[key] for key in ('h1e', 'h2e', 'e_core', 'e_virt'))
    if name == 'Sr':
        v, hv = B[name]
        return sc._Sr(d1, d2, d3, d4, h, g, hv, v, ev)
    if name == 'Si':
        v, hv = B[name]
        return sc._Si(d1, d2, d3, d4, h, g, hv, v, ec)
    if name == 'Sijrs':
        return sc._Sijrs(ec, ev, B[name])
    if name in ('Sijr', 'Srsi'):
        return getattr(sc, '_' + name)(d1, d2, h, g, B[name], ec, ev)
    if name == 'Srs':
        return sc._Srs(d1, d2, d3, h, g, B[name], ev)
    if name == 'Sij':
        return sc._Sij(d1, d2, d3, h, g, B[name], ec)
    v1, v2, hv = B[name]
    return sc._Sir(d1, d2, d3, h, g, hv, v1, v2, ec, ev)


@pytest.mark.parametrize('name', SUBSPACES)
def test_final_contraction_matches_numpy(monkeypatch, name):
    backend = sc._koopmans_lib()
    assert backend is not None, 'Native NEVPT2 backend is required'
    assert hasattr(backend[0], 'nevpt2_' + name.lower())
    # Inspect the actual arrays contracted by the public subspace, before the
    # common denominator function reduces them to a scalar energy.
    monkeypatch.setattr(sc, '_norm_to_energy',
                        lambda norm, h, diff: (norm.copy(), h.copy()))
    for nc, na, nv in ((2, 3, 4), (1, 2, 1), (3, 2, 2), (0, 2, 3), (2, 2, 0)):
        B, dms = _subspace_inputs(nc, na, nv, 470 + nc + na + nv)
        run = lambda: _call_subspace(name, B, dms)
        native = run()
        reference = _python(monkeypatch, run)
        for actual, expected in zip(native, reference):
            assert np.shape(actual) == np.shape(expected), (name, nc, na, nv)
            difference = np.asarray(actual) - np.asarray(expected)
            assert np.max(np.abs(difference), initial=0.0) < 1e-11, (
                name, nc, na, nv, difference)


@pytest.mark.parametrize('nact,nelec', [(2, (1, 1)), (3, (2, 1)), (4, (2, 2))])
def test_fused_rdms_and_f3_match_numpy(monkeypatch, nact, nelec):
    backend = sc._koopmans_lib()
    assert backend is not None and hasattr(backend[0], 'nevpt2_make_f3')
    rng = np.random.default_rng(470 + nact)
    dets = _determinants(nact, nelec)
    ci = rng.normal(size=len(dets))
    ci /= np.linalg.norm(ci)
    g = rng.normal(size=(nact,) * 4) * 0.1
    native = sc._lib_make_f3(ci, nact, dets, g)
    assert native is not None
    # Disable both construction and contraction to get the independent
    # determinant-operator/NumPy implementation, rather than a second native call.
    with monkeypatch.context() as patch:
        patch.setattr(sc, '_lib_make_rdms', lambda *args: None)
        patch.setattr(sc, '_koopmans_lib', lambda: None)
        d1, d2, d3, d4 = sc.make_rdms(ci, nact, nelec)
        f3 = sc._f3ca_f3ac(g, d4)
    for actual, expected in zip((*native[:3], *native[3]), (d1, d2, d3, *f3)):
        assert actual.shape == expected.shape
        assert np.max(np.abs(actual - expected)) < 1e-11


def test_available_native_backend_error_is_not_silenced(monkeypatch):
    def failed(*args):
        raise ValueError('native construction defect')
    monkeypatch.setattr(sc, '_lib_make_f3', failed)
    h = np.eye(4)
    g = np.zeros((4,) * 4)
    with pytest.raises(ValueError, match='native construction defect'):
        sc.sc_nevpt2_energy(h, g, np.arange(4.0), 0, 2, (1, 1),
                           np.ones(4) / 2)


@pytest.mark.parametrize('norb,ndet', [(-1, 1), (0, 1), (32, 1), (1, 0), (1, -1)])
def test_fused_native_rejects_invalid_dimensions_before_allocation(norb, ndet):
    backend = sc._koopmans_lib()
    assert backend is not None
    lib, ffi = backend
    assert hasattr(lib, 'nevpt2_make_f3')
    # Invalid dimensions must be rejected before dereferencing any operand.
    status = lib.nevpt2_make_f3(norb, ndet, *(ffi.NULL for _ in range(8)))
    assert int(status) == -2


def test_fused_error_status_propagates(monkeypatch):
    lib, ffi = sc._koopmans_lib()
    class Failed:
        def nevpt2_make_f3(self, *args):
            return -1
    from oqp.library import fci
    monkeypatch.setattr(fci, '_lib_backend', lambda: (Failed(), ffi))
    with pytest.raises(RuntimeError, match='nevpt2_make_f3 failed with status -1'):
        sc._lib_make_f3(np.ones(1), 1, [3], np.zeros((1,)*4))


def test_memory_guard_precedes_native_allocation(monkeypatch):
    def should_not_run(*args):
        pytest.fail('allocation was reached before memory admission')
    monkeypatch.setattr(sc, '_lib_make_f3', should_not_run)
    # CAS(8,8) needs E2 plus outputs and excitation maps concurrently.
    # The previous estimate admitted this under a 300 MiB ceiling.
    with pytest.raises(ValueError, match='above the'):
        sc.sc_nevpt2_energy(None, None, None, 0, 8, (4, 4), None,
                           max_memory=300)


def test_unsupported_encoding_uses_explicit_fallback():
    # The wrapper must decline >31 without allocating a 32^8 tensor.
    assert sc._lib_make_f3([], 32, [], np.empty(0)) is None


@pytest.mark.parametrize('name', SUBSPACES)
@pytest.mark.parametrize('strided', [False, True])
def test_float32_inputs_are_promoted_at_native_boundary(monkeypatch, name, strided):
    assert sc._koopmans_lib() is not None
    monkeypatch.setattr(sc, '_norm_to_energy',
                        lambda norm, h, diff: (norm.copy(), h.copy()))
    B, dms = _subspace_inputs(2, 3, 4, 90470)
    def convert(x):
        if isinstance(x, tuple):
            return tuple(convert(v) for v in x)
        arr = np.asarray(x, dtype=np.float32)
        if strided:
            storage = np.empty((*arr.shape[:-1], 2*arr.shape[-1]), dtype=np.float32)
            storage[..., ::2] = arr
            arr = storage[..., ::2]
        return arr
    B = {key: convert(x) for key, x in B.items()}
    dms = convert(dms)
    native = _call_subspace(name, B, dms)
    def promote(x):
        if isinstance(x, tuple):
            return tuple(promote(v) for v in x)
        return np.asarray(x, dtype=np.float64)
    reference = _python(monkeypatch, lambda: _call_subspace(
        name, {key: promote(x) for key, x in B.items()}, promote(dms)))
    for actual, expected in zip(native, reference):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-11)


def test_external_workspace_is_admitted_before_packing(monkeypatch):
    # Use a zero-stride view to describe a large external space without
    # allocating the full integral tensor; admission must happen first.
    from numpy.lib.stride_tricks import as_strided
    h = as_strided(np.zeros(1), shape=(202, 202), strides=(0, 0))
    def should_not_run(*args):
        pytest.fail('integral packing reached before memory admission')
    monkeypatch.setattr(sc, '_blocks', should_not_run)
    with pytest.raises(ValueError, match='above the'):
        sc.sc_nevpt2_energy(h, None, None, 100, 2, (1, 1), None, max_memory=16)


@pytest.mark.parametrize('nc,na,nv,nelec', [(1, 2, 2, (1, 1)), (0, 3, 2, (2, 1))])
def test_whole_energy_matches_backend_unavailable(monkeypatch, nc, na, nv, nelec):
    rng = np.random.default_rng(20470 + na)
    nbf = nc + na + nv
    h = rng.normal(size=(nbf, nbf))*.1
    h = (h+h.T)/2
    g = rng.normal(size=(nbf,)*4)*.01
    dets = _determinants(na, nelec)
    ci = rng.normal(size=len(dets)); ci /= np.linalg.norm(ci)
    args = (h, g, np.linspace(-3, 3, nbf), nc, na, nelec, ci)
    native = sc.sc_nevpt2_energy(*args)
    with monkeypatch.context() as patch:
        patch.setattr(sc, '_koopmans_lib', lambda: None)
        patch.setattr(sc, '_lib_make_f3', lambda *args: None)
        patch.setattr(sc, '_lib_make_rdms', lambda *args: None)
        reference = sc.sc_nevpt2_energy(*args)
    assert native[1].keys() == reference[1].keys()
    np.testing.assert_allclose(native[0], reference[0], rtol=0, atol=1e-11)
    for name in SUBSPACES:
        np.testing.assert_allclose(native[1][name], reference[1][name], rtol=0, atol=1e-11)


def test_zero_norm_policy_is_shared():
    norm = np.array([0.0, sc.NUMERICAL_ZERO*.5, 2.0])
    h = np.array([0.0, 1.0, 4.0])
    diff = np.array([0.0, -1.0, 3.0])
    total_norm, energy = sc._norm_to_energy(norm, h, diff)
    assert total_norm == norm.sum()
    assert energy == -2.0/5.0


@pytest.mark.parametrize('nact', [-1, 0, 32, 2147483647])
def test_checked_f3_rejects_invalid_dimensions_before_dereference(nact):
    lib, ffi = sc._koopmans_lib()
    assert hasattr(lib, 'nevpt2_f3ca_f3ac_checked')
    assert int(lib.nevpt2_f3ca_f3ac_checked(
        nact, *(ffi.NULL for _ in range(4)))) == -2


@pytest.mark.parametrize('status', [-1, -2])
def test_checked_f3_failure_reaches_python(monkeypatch, status):
    _, ffi = sc._koopmans_lib()
    class Failed:
        def nevpt2_f3ca_f3ac_checked(self, *args):
            return status
    monkeypatch.setattr(sc, '_koopmans_lib', lambda: (Failed(), ffi))
    with pytest.raises(RuntimeError, match=f'f3ca_f3ac failed with status {status}'):
        sc._f3ca_f3ac(np.zeros((1,)*4), np.zeros((1,)*8))


def test_actual_native_allocation_failures_and_recovery(tmp_path):
    """Fail each f3 scratch allocation and the fused dm4, then reuse the library."""
    import os
    import platform
    import subprocess
    import sys
    from pathlib import Path
    if sys.platform != 'linux' or platform.libc_ver()[0] != 'glibc':
        pytest.skip('actual malloc fault injection requires Linux glibc')
    shim = tmp_path / 'malloc_fault.so'
    source = Path(__file__).with_name('native_malloc_fault.c')
    subprocess.run(['cc', '-shared', '-fPIC', '-O0', str(source),
                    '-o', str(shim)], check=True, capture_output=True)
    program = r"""
import ctypes, sys
import numpy as np
from oqp.library import nevpt2_sc as sc
from oqp.library.fci import _determinants
fault = ctypes.CDLL(sys.argv[1])
fault.oqp_test_arm_malloc.argtypes = [ctypes.c_size_t, ctypes.c_int]
fault.oqp_test_malloc_failures.restype = ctypes.c_int
lib, ffi = sc._koopmans_lib()
n = 3
rng = np.random.default_rng(470)
g = rng.normal(size=(n,)*4)
d4 = rng.normal(size=(n,)*8)
fca = np.empty((n,)*6); fac = np.empty_like(fca)
ptr = lambda a: ffi.cast('double *', a.ctypes.data)
args = (n, ptr(g), ptr(d4), ptr(fca), ptr(fac))
assert int(lib.nevpt2_f3ca_f3ac_checked(*args)) == 0
expected = fca.copy(), fac.copy()
for ordinal in range(3):
    fault.oqp_test_arm_malloc(8*n**4, ordinal)
    assert int(lib.nevpt2_f3ca_f3ac_checked(*args)) == -1, ordinal
    assert fault.oqp_test_malloc_failures() == 1, ordinal
    fault.oqp_test_arm_malloc(0, -1)
    assert int(lib.nevpt2_f3ca_f3ac_checked(*args)) == 0
    np.testing.assert_array_equal(fca, expected[0])
    np.testing.assert_array_equal(fac, expected[1])
dets = np.asarray(_determinants(n, (2, 1)), dtype=np.int64)
ci = np.ones(len(dets)); ci /= np.linalg.norm(ci)
d1 = np.empty((n,)*2); d2 = np.empty((n,)*4); d3 = np.empty((n,)*6)
fargs = (n, len(dets), ffi.cast('int64_t *', dets.ctypes.data),
         ptr(ci), ptr(g), ptr(d1), ptr(d2), ptr(d3), ptr(fca), ptr(fac))
fault.oqp_test_arm_malloc(8*n**8, 0)
assert int(lib.nevpt2_make_f3(*fargs)) == -1
assert fault.oqp_test_malloc_failures() == 1
fault.oqp_test_arm_malloc(0, -1)
assert int(lib.nevpt2_make_f3(*fargs)) == 0
"""
    env = os.environ.copy()
    env['LD_PRELOAD'] = str(shim)
    run = subprocess.run([sys.executable, '-c', program, str(shim)],
                         env=env, text=True, capture_output=True, timeout=60)
    assert run.returncode == 0, run.stdout + run.stderr


def test_legacy_f3_abi_compatibility(monkeypatch):
    lib, ffi = sc._koopmans_lib()
    class Legacy:
        def nevpt2_f3ca_f3ac(self, *args):
            lib.nevpt2_f3ca_f3ac(*args)
    rng = np.random.default_rng(47031)
    g = rng.normal(size=(3,)*4)
    d4 = rng.normal(size=(3,)*8)
    checked = sc._f3ca_f3ac(g, d4)
    monkeypatch.setattr(sc, '_koopmans_lib', lambda: (Legacy(), ffi))
    legacy = sc._f3ca_f3ac(g, d4)
    for actual, expected in zip(legacy, checked):
        np.testing.assert_array_equal(actual, expected)
