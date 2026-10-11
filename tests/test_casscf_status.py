"""Analytic-Hessian failures retain their cause across the native/Python boundary."""
from concurrent.futures import ThreadPoolExecutor
import os
import subprocess
import sys
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


@pytest.fixture
def fault_controller(monkeypatch):
    for name in ("oqp_casscf_test_fault_begin", "oqp_casscf_test_fault_end",
                 "oqp_casscf_test_fault_request"):
        assert hasattr(lib, name), "Missing native Hessian fault controller: " + name
    monkeypatch.setenv("OQP_CASSCF_TEST_FAULTS", "1")
    try:
        yield lib
    finally:
        assert lib.oqp_casscf_test_fault_end(ffi.new("int32_t *")) == 0


def test_hessian_fault_controller_is_gated_and_consumed(fault_controller, monkeypatch):
    ctl = fault_controller
    count = ffi.new("int32_t *")
    monkeypatch.delenv("OQP_CASSCF_TEST_FAULTS")
    assert ctl.oqp_casscf_test_fault_begin(1) == 1
    assert ctl.oqp_casscf_test_fault_request(1) == 0
    monkeypatch.setenv("OQP_CASSCF_TEST_FAULTS", "1")
    with monkeypatch.context() as guard:
        guard.delenv("PYTEST_CURRENT_TEST")
        assert ctl.oqp_casscf_test_fault_begin(1) == 1
        guard.setenv("PYTEST_CURRENT_TEST", "foreign.py::test_other (call)")
        assert ctl.oqp_casscf_test_fault_begin(1) == 1
        assert ctl.oqp_casscf_test_fault_request(1) == 0
    for site in (0, -1, 3):
        assert ctl.oqp_casscf_test_fault_begin(site) == 3
    assert ctl.oqp_casscf_test_fault_begin(1) == 0
    assert ctl.oqp_casscf_test_fault_begin(2) == 2
    assert ctl.oqp_casscf_test_fault_request(2) == 0
    assert ctl.oqp_casscf_test_fault_end(ffi.NULL) == 1
    assert ctl.oqp_casscf_test_fault_request(1) == 1
    assert ctl.oqp_casscf_test_fault_request(1) == 0
    assert ctl.oqp_casscf_test_fault_end(count) == 0 and count[0] == 1
    assert ctl.oqp_casscf_test_fault_end(count) == 0 and count[0] == 0
    assert ctl.oqp_casscf_test_fault_request(1) == 0

    def isolated(site):
        local = ffi.new("int32_t *")
        assert ctl.oqp_casscf_test_fault_begin(site) == 0
        try:
            assert ctl.oqp_casscf_test_fault_request(3-site) == 0
            assert ctl.oqp_casscf_test_fault_request(site) == 1
        finally:
            assert ctl.oqp_casscf_test_fault_end(local) == 0
        return local[0]

    assert ctl.oqp_casscf_test_fault_begin(1) == 0
    with ThreadPoolExecutor(max_workers=2) as pool:
        assert list(pool.map(isolated, [1, 2])) == [1, 1]
    assert ctl.oqp_casscf_test_fault_request(1) == 1
    assert ctl.oqp_casscf_test_fault_end(count) == 0 and count[0] == 1


@pytest.mark.parametrize("threads", [1, 4])
def test_real_hessian_producers_retry_and_recover_same_molecule(fault_controller, tmp_path, threads):
    """Zero-CI INPUT and injected solver EIGEN traverse the real energy driver."""
    child = r'''
import pathlib, sys
import numpy as np
from oqp.library import casscf as cs
from oqp.pyoqp import Runner
requested_threads = int(sys.argv[2])
real_backend = cs._casscf_energy_backend
real_wrapper = cs._lib_casscf_energy
real_optimize = cs._optimize
real_provider = cs._hessian_provider
lib, ffi = real_backend()
for name in ("casscf_energy", "cas_hessian_driver_status", "oqp_casscf_test_fault_begin",
             "oqp_casscf_test_fault_end", "oqp_casscf_test_fault_request"):
    assert hasattr(lib, name), "Missing native entry point: " + name
initial_mo = None
controls = None
phase, fault = "baseline", 0
statuses, retries, python_hessians = [], [], []
energies = {}
class NativeCall:
    def casscf_energy(self, *args):
        global controls
        current = tuple(bytes(ffi.buffer(args[i], n))
                        for i, n in ((1, len(cs._CAS_IOPT)*4),
                                     (2, len(cs._CAS_DOPT)*8), (3, 8), (4, 4)))
        iopt = np.frombuffer(current[0], dtype=np.int32)
        assert iopt[cs._CAS_IOPT_INDEX["nthreads"]] == requested_threads
        assert iopt[cs._CAS_IOPT_INDEX["hessian"]] == 1
        if controls is None:
            controls = current
        assert current == controls, "Native controls changed between calls"
        count = ffi.new("int32_t *")
        if fault:
            assert lib.oqp_casscf_test_fault_begin(fault) == 0
        try:
            # This invokes the actual producer, driver conversion and C ABI.
            status = int(lib.casscf_energy(*args))
        finally:
            assert lib.oqp_casscf_test_fault_end(count) == 0
        if fault:
            assert count[0] == 1, (phase, status, count[0])
            assert status == {1: -1, 2: -5}[fault], (phase, status)
            assert runner.mol.data["OQP::VEC_MO_A"].tobytes() == initial_mo.tobytes()
        else:
            assert count[0] == 0 and status == 0, (phase, status, count[0])
        statuses.append(status)
        return status
def pinned_wrapper(*args, **kwargs):
    global initial_mo
    mol = args[0]
    if initial_mo is None:
        initial_mo = mol.data["OQP::VEC_MO_A"].copy()
    mol.data["OQP::VEC_MO_A"][...] = initial_mo
    return real_wrapper(*args, **kwargs)
def retry(*args, **kwargs):
    retries.append(phase)
    return real_optimize(*args, **kwargs)
def provider(*args, **kwargs):
    actual = real_provider(*args, **kwargs)
    assert actual is not None
    def compute(*a, **kw):
        python_hessians.append(phase)
        return actual(*a, **kw)
    return compute
config = {
    "input": {"system": "\nLi 0 0 0\nH 0 0 1.6", "basis": "sto-3g",
              "method": "casscf", "runtype": "energy", "omp_threads": str(requested_threads)},
    "guess": {"type": "hcore"},
    "scf": {"type": "rhf", "multiplicity": "1", "save_molden": "False"},
    "cas": {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1"},
    "ci": {"nroot": "1", "solver": "dense", "target_spin": "any",
           "integral_backend": "native", "root_tracking": "energy"},
    "casscf": {"max_macro_iterations": "50", "canonicalize": "true",
               "converger": "ah", "hessian": "analytic", "gradient_norm_tol": "1.0e-7",
               "energy_decrease_tol": "1.0e-10", "step_norm_tol": "1.0e-8"},
    "tests": {"exception": "True"},
}
runner = Runner(project="hessian-status-recovery", input_file=None, silent=1, usempi=False,
                input_dict=config, log=str(pathlib.Path(sys.argv[1])/"hessian.log"))
molecule = runner.mol
cs._lib_casscf_energy = pinned_wrapper
cs._optimize = retry
cs._hessian_provider = provider
try:
    for phase, fault in (("baseline", 0), ("python-control", 0), ("input", 1),
                         ("input-recovery1", 0), ("input-recovery2", 0), ("eigen", 2),
                         ("eigen-recovery1", 0), ("eigen-recovery2", 0)):
        cs._casscf_energy_backend = (lambda: None) if phase == "python-control" else (lambda: (NativeCall(), ffi))
        runner.run(test_mod=True)
        assert runner.mol is molecule
        energies[phase] = np.array(runner.mol.data["OQP::CASSCF_ENERGIES"], copy=True)
finally:
    cs._casscf_energy_backend = real_backend
    cs._lib_casscf_energy = real_wrapper
    cs._optimize = real_optimize
    cs._hessian_provider = real_provider
assert statuses == [0, -1, 0, 0, -5, 0, 0], statuses
assert retries == ["python-control", "input", "eigen"], retries
for name in retries:
    assert name in python_hessians, (name, python_hessians)
for name in ("input", "eigen"):
    assert energies[name].tobytes() == energies["python-control"].tobytes(), (name, energies)
for name in ("input-recovery1", "input-recovery2", "eigen-recovery1", "eigen-recovery2"):
    assert energies[name].tobytes() == energies["baseline"].tobytes(), (name, energies)
# Existing native/Python orbital-optimizer agreement bound; no tolerance changes.
assert np.max(np.abs(energies["baseline"]-energies["python-control"])) < 1.0e-7, energies
print("INPUT/EIGEN producers reached; actual Python Hessian retries and two native recoveries passed")
'''
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OQP_FCI_THREADS"] = str(threads)
    result = subprocess.run([sys.executable, "-c", child, str(tmp_path), str(threads)],
                            env=env, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "actual Python Hessian retries and two native recoveries passed" in result.stdout
