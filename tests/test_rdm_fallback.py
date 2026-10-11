"""Native string RDM acceptance, numerical pin, and caller-state lifetime."""

import ctypes
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from threading import Barrier

import numpy as np
import pytest

from oqp.library import rdm
from oqp import runtime


@pytest.fixture
def native():
    backend = rdm._rdm_backend()
    assert backend is not None, "RDM tests require the installed native library"
    lib, ffi = backend
    assert hasattr(lib, "rdm12_strings_c"), "Missing string RDM entry point"
    assert hasattr(lib, "rdm2_spatial"), "Missing Gram RDM entry point"
    path = runtime.library_path(*runtime.resolve_oqp_root())
    ctl = ctypes.CDLL(str(path))
    ctl.oqp_blas_thread_count.restype = ctypes.c_int64
    ctl.oqp_blas_thread_set.argtypes = [ctypes.c_int64]
    return lib, ffi, ctl


@contextmanager
def _caller_blas(native):
    ctl = native[2]
    saved = ctl.oqp_blas_thread_count()
    try:
        if saved > 0:
            ctl.oqp_blas_thread_set(3)
        requested = ctl.oqp_blas_thread_count()
        if os.environ.get("OQP_RDM_REQUIRE_BLAS3") == "1":
            assert requested == 3, "This validation lane requires BLAS threads >1"
        yield
        assert ctl.oqp_blas_thread_count() == requested
    finally:
        ctl.oqp_blas_thread_set(saved)


def _ci(norb, nelec, seed=742):
    dets = np.ascontiguousarray(rdm.determinant_basis(norb, nelec), dtype=np.int64)
    c = np.random.default_rng(seed).standard_normal(len(dets))
    c /= np.linalg.norm(c)
    return dets, np.ascontiguousarray(c)


def _strings(native, norb, dets, c, threads, d1=None, d2=None):
    lib, ffi, _ = native
    if d1 is None:
        d1 = np.full((norb, norb), 23.0)
    if d2 is None:
        d2 = np.full((norb,) * 4, 29.0)
    info = lib.rdm12_strings_c(
        norb, len(dets), ffi.cast("int64_t *", dets.ctypes.data),
        ffi.cast("double *", c.ctypes.data),
        ffi.cast("double *", d1.ctypes.data), ffi.cast("double *", d2.ctypes.data), threads)
    return int(info), d1, d2


def _spin_gram(c, dets, norb):
    spin = rdm._lib_rdm2(c, dets, 2 * norb)
    assert spin is not None, "Independent spin-orbital Gram pin must execute"
    ref = np.zeros((norb,) * 4)
    for a in (0, norb):
        for b in (0, norb):
            ref += spin[a:a+norb, b:b+norb, a:a+norb, b:b+norb].transpose(0, 2, 1, 3)
    return ref


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("norb,nelec", [(4, (2, 1)), (8, (4, 4))])
def test_canonical_strings_match_signed_gram(native, threads, norb, nelec):
    # 24*16 < 65536, 4900*64 > 65536: both sides of the OpenMP condition.
    dets, c = _ci(norb, nelec)
    with _caller_blas(native):
        info, d1, d2 = _strings(native, norb, dets, c, threads)
    assert info == 0
    ref1 = rdm._lib_rdm1_spatial(c, dets, norb)
    ref2 = _spin_gram(c, dets, norb)
    assert np.max(np.abs(ref2 - ref2.swapaxes(0, 1))) > 1e-4
    np.testing.assert_allclose(d1, ref1, rtol=0, atol=2e-12)
    np.testing.assert_allclose(d2, ref2, rtol=0, atol=2e-12)
    np.testing.assert_allclose(np.einsum("pqrr->pq", d2), (sum(nelec)-1)*d1,
                               rtol=0, atol=2e-12)
    # The default Gram route retains its stronger, pre-existing bit contract.
    assert np.array_equal(rdm.make_rdm2_spatial(c, dets, norb), ref2)


@pytest.mark.parametrize("kind", ["truncated", "permuted"])
@pytest.mark.parametrize("threads", [1, 4])
def test_noncanonical_declines_without_mutating_caller(native, kind, threads):
    dets, c = _ci(4, (2, 2))
    if kind == "truncated":
        dets, c = dets[:7].copy(), c[:7].copy()
        c /= np.linalg.norm(c)
    else:
        order = np.random.default_rng(63).permutation(len(dets))
        dets, c = dets[order].copy(), c[order].copy()
    with _caller_blas(native):
        info, d1, d2 = _strings(native, 4, dets, c, threads)
        assert rdm.make_rdm12_spatial_strings(c, dets, 4) is None
    assert info > 0
    assert np.all(d1 == 23) and np.all(d2 == 29)
    assert np.array_equal(rdm.make_rdm2_spatial(c, dets, 4), _spin_gram(c, dets, 4))


@pytest.mark.parametrize("norb,dets", [(0, [0]), (-1, [0]), (32, [0]),
                                      (64, [0]), (2, [1 << 10]), (2, [])])
def test_native_declining_boundaries_preserve_buffers_and_threads(native, norb, dets):
    dets = np.ascontiguousarray(dets, dtype=np.int64)
    c = np.ones(len(dets))
    d1, d2 = np.array([23.0]), np.array([29.0])
    with _caller_blas(native):
        info, _, _ = _strings(native, norb, dets, c, 4, d1, d2)
    assert info > 0
    assert d1[0] == 23 and d2[0] == 29


@pytest.mark.parametrize("nelec", [(0, 0), (1, 0), (31, 31)])
def test_supported_upper_boundary(native, nelec):
    dets, c = _ci(31, nelec)
    with _caller_blas(native):
        info, d1, d2 = _strings(native, 31, dets, c, 4)
    assert info == 0
    assert np.trace(d1) == pytest.approx(sum(nelec), abs=2e-12)
    if sum(nelec) <= 1:
        np.testing.assert_allclose(d2, 0, rtol=0, atol=2e-12)
    else:
        eye = np.eye(31)
        ref = 4*np.einsum("pq,rs->pqrs", eye, eye) - 2*np.einsum("ps,rq->pqrs", eye, eye)
        assert np.array_equal(d1, 2*eye)
        assert np.array_equal(d2, ref)


def test_python_propagates_native_allocation_failure(native, monkeypatch):
    class AllocationFailure:
        def rdm12_strings_c(self, *args):
            return -1

    monkeypatch.setattr(rdm, "_rdm_backend", lambda: (AllocationFailure(), native[1]))
    dets, c = _ci(4, (2, 1))
    with pytest.raises(MemoryError, match="RDM scratch"):
        rdm.make_rdm12_spatial_strings(c, dets, 4)


@pytest.mark.parametrize("threads", [1, 2])
def test_overlapping_rdm_calls_restore_blas(native, threads):
    inputs = [_ci(8, (4, 4), seed=742 + i) for i in range(4)]
    with _caller_blas(native):
        reference = [_strings(native, 8, dets, c, threads) for dets, c in inputs]
        with ThreadPoolExecutor(max_workers=2) as pool:
            runs = list(pool.map(lambda pair: _strings(native, 8, *pair, threads), inputs))
    assert all(result[0] == 0 for result in runs)
    assert all(result[0] == 0 for result in reference)
    for result, ref in zip(runs, reference):
        assert result[1].tobytes() == ref[1].tobytes()
        assert result[2].tobytes() == ref[2].tobytes()


def test_repeated_dimensions_and_short_final_block(native):
    lib = native[0]
    saved = lib.fci_get_work_bytes_cap()
    try:
        for norb, nelec in [(4, (2, 1)), (10, (5, 5)), (2, (1, 1)), (8, (4, 4)), (4, (2, 1))]:
            dets, c = _ci(norb, nelec)
            lib.fci_set_work_bytes_cap(0)
            with _caller_blas(native):
                info, ref1, ref2 = _strings(native, norb, dets, c, 1)
            assert info == 0
            # The capped CAS(10,10) call traverses several alpha-string
            # blocks, including a short final block; smaller cases fit whole.
            lib.fci_set_work_bytes_cap(16*1024**2)
            with _caller_blas(native):
                info, d1, d2 = _strings(native, norb, dets, c, 4)
            assert info == 0
            np.testing.assert_allclose(d1, ref1, rtol=0, atol=2e-12)
            np.testing.assert_allclose(d2, ref2, rtol=0, atol=2e-12)
    finally:
        lib.fci_set_work_bytes_cap(saved)
    assert lib.fci_get_work_bytes_cap() == saved


@pytest.fixture
def allocation_probe(native, monkeypatch):
    """Portable synthetic failures at owned ALLOCATE sites, before workers start."""
    lib, ffi, _ = native
    for name in ("oqp_rdm_test_alloc_begin", "oqp_rdm_test_alloc_end",
                 "oqp_rdm_test_alloc_active", "oqp_rdm_test_alloc_request"):
        assert hasattr(lib, name), f"Missing linked allocation controller: {name}"
    monkeypatch.setenv("OQP_RDM_TEST_ALLOCATIONS", "1")
    assert lib.oqp_rdm_test_alloc_active() == 0
    try:
        yield lib, ffi
    finally:
        stats = ffi.new("struct oqp_rdm_test_alloc_stats *")
        assert lib.oqp_rdm_test_alloc_end(stats) == 0


@contextmanager
def _allocation_request(probe, ordinal):
    lib, ffi = probe
    stats = ffi.new("struct oqp_rdm_test_alloc_stats *")
    assert lib.oqp_rdm_test_alloc_begin(ordinal) == 0
    try:
        yield stats
    finally:
        assert lib.oqp_rdm_test_alloc_end(stats) == 0


def test_allocation_controller_guards_and_cleanup(allocation_probe, monkeypatch):
    lib, ffi = allocation_probe
    stats = ffi.new("struct oqp_rdm_test_alloc_stats *")
    with monkeypatch.context() as patch:
        patch.delenv("OQP_RDM_TEST_ALLOCATIONS")
        assert lib.oqp_rdm_test_alloc_begin(1) == 1
    with monkeypatch.context() as patch:
        patch.setenv("PYTEST_CURRENT_TEST", "tests/test_other.py::test_other (call)")
        assert lib.oqp_rdm_test_alloc_begin(1) == 1
    for ordinal in (-1, 19):
        assert lib.oqp_rdm_test_alloc_begin(ordinal) == 3
        assert lib.oqp_rdm_test_alloc_active() == 0
    with _allocation_request(allocation_probe, 1) as consumed:
        assert lib.oqp_rdm_test_alloc_begin(2) == 2
        assert lib.oqp_rdm_test_alloc_end(ffi.NULL) == 1
        assert lib.oqp_rdm_test_alloc_active() == 1
        assert lib.oqp_rdm_test_alloc_request(2, 3, 4, 8) == 1
        assert lib.oqp_rdm_test_alloc_request(2, 3, 4, 8) == 0
    assert (consumed.count, consumed.total_bytes, consumed.largest_bytes) == (2, 384, 192)
    assert lib.oqp_rdm_test_alloc_active() == 0
    assert lib.oqp_rdm_test_alloc_end(stats) == 0
    assert (stats.count, stats.total_bytes, stats.largest_bytes) == (0, 0, 0)
    # Unarmed bridge calls cannot inject or accumulate.
    assert lib.oqp_rdm_test_alloc_request(2, 3, 4, 8) == 0
    assert lib.oqp_rdm_test_alloc_end(stats) == 0
    assert stats.count == 0


def test_allocation_controller_checked_sizes(allocation_probe):
    lib, _ = allocation_probe
    with _allocation_request(allocation_probe, 0) as stats:
        for dims in ((-1, 1, 1, 8), (1, -1, 1, 8), (1, 1, -1, 8),
                     (1, 1, 1, 0), (2**62, 4, 1, 8),
                     (2**62, 1, 4, 8), (2**61, 1, 1, 8)):
            assert lib.oqp_rdm_test_alloc_request(*dims) == 1
        assert lib.oqp_rdm_test_alloc_request(2, 3, 4, 8) == 0
        # Each is valid alone, but the third would overflow the aggregate.
        for _ in range(2):
            assert lib.oqp_rdm_test_alloc_request(2**61, 1, 1, 1) == 0
        assert lib.oqp_rdm_test_alloc_request(2**62, 1, 1, 1) == 1
    assert stats.count == 3
    assert stats.total_bytes == 192 + 2**62
    assert stats.largest_bytes == 2**61


def test_injected_native_failure_raises_without_fallback(native, allocation_probe):
    dets, c = _ci(4, (2, 1))
    with _caller_blas(native):
        with pytest.raises(MemoryError, match="RDM scratch"):
            with _allocation_request(allocation_probe, 1) as stats:
                rdm.make_rdm12_spatial_strings(c, dets, 4)
        assert stats.count == 1
        assert native[0].oqp_rdm_test_alloc_active() == 0
        recovered = rdm.make_rdm12_spatial_strings(c, dets, 4)
        assert recovered is not None
        np.testing.assert_allclose(recovered[1], _spin_gram(c, dets, 4), rtol=0, atol=2e-12)


@pytest.mark.parametrize("threads", [1, 4])
def test_each_native_allocation_failure_restores_caller(native, allocation_probe, threads):
    dets, c = _ci(4, (2, 1))
    with _caller_blas(native):
        requested_blas = native[2].oqp_blas_thread_count()
        info, ref1, ref2 = _strings(native, 4, dets, c, threads)
        assert info == 0
        for ordinal in range(1, 19):
            with _allocation_request(allocation_probe, ordinal) as stats:
                info, d1, d2 = _strings(native, 4, dets, c, threads)
            assert stats.count == ordinal, (ordinal, stats.count)
            assert info == -1, (ordinal, info)
            assert native[2].oqp_blas_thread_count() == requested_blas, ordinal
            assert np.all(d1 == 23) and np.all(d2 == 29), ordinal
            # Fresh output buffers and no armed controller on recovery.
            for _ in range(2):
                info, d1, d2 = _strings(native, 4, dets, c, threads)
                assert info == 0
                assert native[2].oqp_blas_thread_count() == requested_blas, ordinal
                assert d1.tobytes() == ref1.tobytes()
                assert d2.tobytes() == ref2.tobytes()


def test_concurrent_allocation_controllers_are_thread_local(native, allocation_probe):
    lib, _ = allocation_probe
    barrier = Barrier(2)
    dets, c = _ci(4, (2, 1))
    def run(ordinal):
        with _allocation_request(allocation_probe, ordinal) as stats:
            barrier.wait(timeout=10)
            failed, d1, d2 = _strings(native, 4, dets, c, 1)
            assert np.all(d1 == 23) and np.all(d2 == 29)
            # The failure is consumed once; recovery on this OS thread stays armed.
            recovered = _strings(native, 4, dets, c, 1)[0]
        return failed, recovered, stats.count
    with _caller_blas(native):
        with ThreadPoolExecutor(max_workers=2) as pool:
            result = list(pool.map(run, (1, 18)))
    assert result == [(-1, 0, 19), (-1, 0, 36)]
    assert lib.oqp_rdm_test_alloc_active() == 0


def test_requested_width_respects_scratch_budget(native, allocation_probe):
    """Count requested payload bytes; this is not allocator overhead or RSS."""
    lib, _ = allocation_probe
    dets, c = _ci(31, (1, 1))
    serial = os.environ.get("OQP_RDM_SERIAL_BUILD") == "1"
    saved = lib.fci_get_work_bytes_cap()
    observed = []
    try:
        lib.fci_set_work_bytes_cap(512*1024**2 if serial else 16*1024**2)
        with _caller_blas(native):
            for width in (1, 64):
                with _allocation_request(allocation_probe, 0) as stats:
                    info, d1, d2 = _strings(native, 31, dets, c, width)
                assert info == 0
                assert stats.count == 18
                observed.append((stats.total_bytes, stats.largest_bytes))
                assert stats.largest_bytes <= 31**4*8, observed
                if width == 1:
                    ref1, ref2 = d1.copy(), d2.copy()
                else:
                    np.testing.assert_allclose(d1, ref1, rtol=0, atol=2e-12)
                    np.testing.assert_allclose(d2, ref2, rtol=0, atol=2e-12)
        if serial:
            assert observed[0] == observed[1], observed
    finally:
        lib.fci_set_work_bytes_cap(saved)


def test_minimum_block_budget_failure_preserves_caller(native):
    lib = native[0]
    saved = lib.fci_get_work_bytes_cap()
    try:
        # One beta interval, Gram, and one partial exceed 16 MiB here.
        dets, c = _ci(31, (1, 2))
        lib.fci_set_work_bytes_cap(16*1024**2)
        with _caller_blas(native):
            info, d1, d2 = _strings(native, 31, dets, c, 4)
        assert info == -1
        assert np.all(d1 == 23) and np.all(d2 == 29)
        dets, c = _ci(4, (2, 1))
        assert _strings(native, 4, dets, c, 4)[0] == 0
    finally:
        lib.fci_set_work_bytes_cap(saved)


@pytest.mark.parametrize("threads", [1, 4])
def test_native_casscf_propagates_rdm_allocation_failure(native, allocation_probe, tmp_path, threads):
    """One Runner/Molecule: success, native failure, then two top-level recoveries."""
    child = r'''
import pathlib, sys
import numpy as np
from oqp.library import casscf as cs
from oqp.pyoqp import Runner
requested_threads = int(sys.argv[2])
real_backend = cs._casscf_energy_backend
lib, ffi = real_backend()
for name in ("casscf_energy", "oqp_rdm_test_alloc_begin", "oqp_rdm_test_alloc_end"):
    assert hasattr(lib, name), "Missing native entry point: " + name
initial_mo = None
controls = None
seen = []
phase = "success"
class FailRDM:
    def casscf_energy(self, *args):
        global initial_mo, controls
        # All four calls share a Molecule. Pin the native input orbitals after
        # the top-level SCF, so separate SCF reductions cannot mask recovery.
        if initial_mo is None:
            initial_mo = runner.mol.data["OQP::VEC_MO_A"].copy()
        runner.mol.data["OQP::VEC_MO_A"][...] = initial_mo
        current = tuple(bytes(ffi.buffer(args[i], n))
                        for i, n in ((1, len(cs._CAS_IOPT)*4),
                                     (2, len(cs._CAS_DOPT)*8), (3, 8), (4, 4)))
        iopt = np.frombuffer(current[0], dtype=np.int32)
        assert iopt[cs._CAS_IOPT_INDEX["nthreads"]] == requested_threads
        if controls is None:
            controls = current
        assert current == controls, "Recovery changed native controls"
        if phase == "failure":
            stats = ffi.new("struct oqp_rdm_test_alloc_stats *")
            assert lib.oqp_rdm_test_alloc_begin(1) == 0
            try:
                status = int(lib.casscf_energy(*args))
            finally:
                assert lib.oqp_rdm_test_alloc_end(stats) == 0
            assert stats.count == 1, (status, stats.count)
            assert status == -2, status
            # Reacquire the native record; never keep a borrowed view across calls.
            assert runner.mol.data["OQP::VEC_MO_A"].tobytes() == initial_mo.tobytes()
        else:
            status = int(lib.casscf_energy(*args))
            assert status == 0, status
            native_energies.append(bytes(ffi.buffer(args[5], 8)))
        seen.append(status)
        return status
config = {
    "input": {"system": "\nLi 0 0 0\nH 0 0 1.6", "basis": "sto-3g",
              "method": "casscf", "runtype": "energy", "omp_threads": requested_threads},
    "guess": {"type": "hcore"},
    "scf": {"type": "rhf", "multiplicity": "1", "save_molden": "False"},
    "cas": {"active_electrons": "2", "active_orbitals": "2", "frozen_core": "1"},
    "ci": {"nroot": "1", "solver": "dense", "target_spin": "any",
           "integral_backend": "native", "root_tracking": "energy"},
    "casscf": {"max_macro_iterations": "50", "canonicalize": "true"},
    "tests": {"exception": "True"},
}
runner = Runner(project="persistent-recovery", input_file=None, silent=1, usempi=False,
                input_dict=config, log=str(pathlib.Path(sys.argv[1])/"recovery.log"))
molecule = runner.mol
native_energies = []
real_optimize = cs._optimize
def forbidden_fallback(*args, **kwargs):
    raise AssertionError("Allocation failure entered the Python optimizer")
cs._optimize = forbidden_fallback
cs._casscf_energy_backend = lambda: (FailRDM(), ffi)
try:
    for phase in ("success", "failure", "recovery1", "recovery2"):
        if phase == "failure":
            try:
                runner.run(test_mod=True)
            except MemoryError as exc:
                assert "Native CASSCF scratch" in str(exc)
            else:
                raise AssertionError("Native allocation failure was swallowed")
        else:
            runner.run(test_mod=True)
        assert runner.mol is molecule
finally:
    cs._casscf_energy_backend = real_backend
    cs._optimize = real_optimize
assert seen == [0, -2, 0, 0], seen
assert native_energies[0] == native_energies[1] == native_energies[2], native_energies
print("native CASSCF status -2 propagated; persistent top-level recoveries compared")
'''
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OQP_FCI_THREADS"] = str(threads)
    result = subprocess.run([sys.executable, "-c", child, str(tmp_path), str(threads)],
                            env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "persistent top-level recoveries compared" in result.stdout
