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
