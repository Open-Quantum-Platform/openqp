"""Analytic-Hessian failures retain their cause across the native/Python boundary."""
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace

import numpy as np
import pytest

from oqp import ffi, lib
from oqp.library import casscf as cs


@pytest.fixture
def status_bridge():
    assert hasattr(lib, "cas_hessian_driver_status"), "Missing Hessian status export"
    return lib.cas_hessian_driver_status


@pytest.mark.parametrize("workers", [1, 4])
def test_native_hessian_statuses_keep_distinct_causes(status_bridge, workers):
    # Unknown failures retain the input refusal, never claim allocation failure.
    cases = [(0, 0), (-1, -1), (-2, -2), (-5, -5), (-8, -8), (-9, -9), (-77, -1)]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        result = list(pool.map(status_bridge, [source for source, _ in cases]))
    assert result == [expected for _, expected in cases]


@pytest.mark.parametrize("hessian_status,driver_status", [(0, 0), (-1, -1), (-2, -2),
                                                         (-5, -5), (-8, -8), (-9, -9),
                                                         (-77, -1)])
def test_python_requests_fallback_for_nonallocation_failures(status_bridge, monkeypatch,
                                                   hessian_status, driver_status):
    # Use production conversion, then exercise the actual Python dispatcher.
    # Dummy packing isolates status handling from molecule/optimizer numerics.
    calls = []
    class FailingHessian:
        def casscf_energy(self, *args):
            calls.append(hessian_status)
            assert len(args) == 9
            return status_bridge(hessian_status)

    monkeypatch.setattr(cs, "_casscf_energy_backend", lambda: (FailingHessian(), ffi))
    monkeypatch.setattr(cs, "_cas_wavefunction_arrays", lambda *args, **kwargs: (
        np.zeros(len(cs._CAS_IOPT), dtype=np.int32),
        np.zeros(len(cs._CAS_DOPT)), SimpleNamespace(nroot=1)))
    mol = SimpleNamespace(data=SimpleNamespace(_data=ffi.NULL))
    args = (mol, None, cs.CASSCFOptions(), 2, 0, 2, (1, 1), 0.0, [1.0], [0])
    if driver_status == -2:
        with pytest.raises(MemoryError, match="Native CASSCF scratch allocation failed"):
            cs._lib_casscf_energy(*args, hessian=1)
    elif driver_status == 0:
        result = cs._lib_casscf_energy(*args, hessian=1)
        assert result is not None
        assert len(result[0]) == 1
    else:
        assert cs._lib_casscf_energy(*args, hessian=1) is None
    assert calls == [hessian_status]
